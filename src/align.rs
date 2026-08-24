//! Competitive minimap2 alignment, ambiguity auditing, and bounded ordered parallelism.
//!
//! Host, target, decoy, and contaminant roles compete in one single-shard index. The ordinary
//! mapper handles clear hits; selected reads use ALL_CHAINS auditing so hidden alternatives cannot
//! become false empty evidence. Role margins, MAPQ, edit-distance, split, and discordant thresholds
//! are centralized here. Worker results are emitted in input order, with caller-runs backpressure
//! and explicit propagation of stream, worker, sink, and audit-overflow failures.

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::mpsc::{sync_channel, TrySendError};
use std::sync::Arc;
use std::thread;

use minimap2::{Aligner, Built, Mapping};

use crate::reference::Role;

pub const MIN_MAPQ: u32 = 20;
pub const MIN_AS_DIFF: i32 = 12;
pub const MAX_NM: i32 = 8;

pub const SPLIT_SOFTCLIP: i32 = 20;
pub const SPLIT_MAPQ: u32 = 10;

pub const BEST_N: i32 = 15;

pub const MAX_AUDIT_HITS: usize = 4096;

const TOTAL_KALLOC_BUDGET: i64 = 1_000_000_000;

pub struct CompetitiveAligner {
    aligner: Arc<Aligner<Built>>,

    audit_aligner: Arc<Aligner<Built>>,

    pub threads: usize,
}

impl CompetitiveAligner {
    pub fn open(mmi_path: &Path, workers: usize) -> Result<Self, &'static str> {
        let mut builder = Aligner::builder().sr().with_cigar();
        builder.mapopt.best_n = BEST_N;

        let mapper_concurrency = workers.saturating_add(1);
        let mapper_concurrency = i64::try_from(mapper_concurrency).unwrap_or(i64::MAX);
        builder.mapopt.cap_kalloc = (TOTAL_KALLOC_BUDGET / mapper_concurrency).max(1);
        let aligner = builder.with_index(mmi_path, None)?;
        crate::reference::ensure_single_part_index(aligner.idx_parts.len())?;

        let mut audit_aligner = aligner.clone();
        audit_aligner.mapopt.flag |= minimap2::ffi::MM_F_ALL_CHAINS as i64;
        Ok(Self {
            aligner: Arc::new(aligner),
            audit_aligner: Arc::new(audit_aligner),
            threads: workers,
        })
    }

    pub fn map_pair(
        &self,
        r1: &[u8],
        r2: &[u8],
        name: &str,
    ) -> Result<(Vec<Mapping>, Vec<Mapping>), &'static str> {
        self.aligner
            .map_pair(r1, r2, false, false, None, None, Some(name.as_bytes()))
    }

    pub fn map_seq(&self, seq: &[u8], name: &str) -> Result<Vec<Mapping>, &'static str> {
        self.aligner
            .map(seq, false, false, None, None, Some(name.as_bytes()))
    }

    pub fn audit_map_seq(&self, seq: &[u8], name: &str) -> Result<Vec<Mapping>, &'static str> {
        self.audit_aligner
            .map(seq, false, false, None, None, Some(name.as_bytes()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    pub role: Role,
    pub contig: String,
    pub as_score: i32,
    pub mapq: u32,
    pub is_primary: bool,
    pub is_supplementary: bool,
    pub nm: i32,
    pub qstart: i32,
    pub qend: i32,
    pub qlen: i32,
    pub tstart: i32,
    pub tend: i32,
    pub strand: char,
}

pub fn hits_of(mappings: &[Mapping], roles: &HashMap<String, Role>) -> Vec<Hit> {
    mappings
        .iter()
        .filter_map(|m| {
            let contig = m.target_name.as_ref()?.to_string();
            let role = *roles.get(&contig)?;
            let aln = m.alignment.as_ref()?;
            Some(Hit {
                role,
                contig,
                as_score: aln.alignment_score.unwrap_or(0),
                mapq: m.mapq,
                is_primary: m.is_primary,
                is_supplementary: m.is_supplementary,
                nm: aln.nm,
                qstart: m.query_start,
                qend: m.query_end,
                qlen: m.query_len.map(|l| l.get()).unwrap_or(m.query_end.max(1)),
                tstart: m.target_start,
                tend: m.target_end,
                strand: if matches!(m.strand, minimap2::Strand::Forward) {
                    '+'
                } else {
                    '-'
                },
            })
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AuditDecision {
    Resolved { role: Role, hits: Vec<Hit> },

    NoEvidence,

    Overflow { hit_count: usize, limit: usize },
}

fn audit_hit_is_better(candidate: &Hit, current: &Hit) -> bool {
    candidate
        .as_score
        .cmp(&current.as_score)
        .then_with(|| current.nm.cmp(&candidate.nm))
        .then_with(|| candidate.mapq.cmp(&current.mapq))
        .then_with(|| candidate.is_primary.cmp(&current.is_primary))
        .then_with(|| current.is_supplementary.cmp(&candidate.is_supplementary))
        .is_gt()
}

fn resolve_audit_role(hits: &[Hit], role: Role) -> Option<Vec<Hit>> {
    let role_best = hits
        .iter()
        .filter(|hit| hit.role == role)
        .map(|hit| hit.as_score)
        .max()?;
    let other_best = hits
        .iter()
        .filter(|hit| hit.role != role)
        .map(|hit| hit.as_score)
        .max()
        .unwrap_or(0);
    if i64::from(role_best) - i64::from(other_best) < i64::from(MIN_AS_DIFF) {
        return None;
    }

    let mut by_contig: BTreeMap<&str, Hit> = BTreeMap::new();
    for hit in hits
        .iter()
        .filter(|hit| hit.role == role && hit.nm <= MAX_NM && hit.as_score == role_best)
    {
        match by_contig.get_mut(hit.contig.as_str()) {
            Some(current) if audit_hit_is_better(hit, current) => *current = hit.clone(),
            Some(_) => {}
            None => {
                by_contig.insert(hit.contig.as_str(), hit.clone());
            }
        }
    }
    let mut resolved = by_contig.into_values().collect::<Vec<_>>();
    resolved.sort_by(|a, b| {
        b.as_score
            .cmp(&a.as_score)
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.nm.cmp(&b.nm))
            .then_with(|| b.mapq.cmp(&a.mapq))
    });
    (!resolved.is_empty()).then_some(resolved)
}

/// Conservatively adjudicate one read's ALL_CHAINS hits across reference roles.
pub fn adjudicate_audit(hits: &[Hit]) -> AuditDecision {
    if hits.len() > MAX_AUDIT_HITS {
        return AuditDecision::Overflow {
            hit_count: hits.len(),
            limit: MAX_AUDIT_HITS,
        };
    }
    if let Some(hits) = resolve_audit_role(hits, Role::Target) {
        return AuditDecision::Resolved {
            role: Role::Target,
            hits,
        };
    }
    if let Some(hits) = resolve_audit_role(hits, Role::Decoy) {
        return AuditDecision::Resolved {
            role: Role::Decoy,
            hits,
        };
    }
    AuditDecision::NoEvidence
}

fn best_hit(hits: &[Hit], role: Role) -> Option<&Hit> {
    hits.iter()
        .filter(|h| h.role == role)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
}

fn best_non_target(hits: &[Hit]) -> Option<&Hit> {
    hits.iter()
        .filter(|h| h.role != Role::Target)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
}

#[derive(Debug, Clone, Default)]
pub struct ReadDecision {
    pub confident_target: Option<Hit>,

    pub best_host: Option<Hit>,

    pub split_target: Option<Hit>,

    pub split_host: Option<Hit>,

    pub best_overall: Option<Hit>,
}

pub fn adjudicate_read(hits: &[Hit]) -> ReadDecision {
    let best_other = best_non_target(hits).cloned();
    let best_host = best_hit(hits, Role::Host).cloned();
    let best_overall = hits
        .iter()
        .filter(|h| h.mapq >= MIN_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
        .cloned();

    let confident_target = {
        let best_target_conf = hits
            .iter()
            .filter(|h| h.role == Role::Target && h.mapq >= MIN_MAPQ)
            .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
            .cloned();
        best_target_conf.as_ref().and_then(|t| {
            let other_as = best_other.as_ref().map(|o| o.as_score).unwrap_or(0);
            if t.nm <= MAX_NM && t.as_score - other_as >= MIN_AS_DIFF {
                Some(t.clone())
            } else {
                None
            }
        })
    };

    let split_target = hits
        .iter()
        .filter(|h| h.role == Role::Target && h.mapq >= SPLIT_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score))
        .cloned();
    let split_host = hits
        .iter()
        .filter(|h| h.role == Role::Host && h.mapq >= SPLIT_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score))
        .cloned();

    ReadDecision {
        confident_target,
        best_host,
        split_target,
        split_host,
        best_overall,
    }
}

/// Return whether an ordinary mapping requires escalation to an ALL_CHAINS ambiguity audit.
pub fn needs_ambiguity_audit(hits: &[Hit], ordinary: &ReadDecision) -> bool {
    let has_low_mapq_primary = hits.iter().any(|hit| hit.is_primary && hit.mapq < MIN_MAPQ);
    let secondary_saturated = hits
        .iter()
        .filter(|hit| !hit.is_primary && !hit.is_supplementary)
        .count()
        >= BEST_N as usize;
    let has_unresolved_target =
        ordinary.confident_target.is_none() && hits.iter().any(|hit| hit.role == Role::Target);
    has_low_mapq_primary || secondary_saturated || has_unresolved_target
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitEvent {
    pub target_contig: String,
    pub target_pos: i32,
    pub host_contig: String,
    pub host_pos: i32,
    pub direction: String,

    pub qname: String,
}

fn split_event(read: &ReadDecision, qname: &str) -> Option<SplitEvent> {
    let t = read.split_target.as_ref()?;
    let h = read.split_host.as_ref()?;
    let t_left = t.qstart >= SPLIT_SOFTCLIP;
    let t_right = t.qlen - t.qend >= SPLIT_SOFTCLIP;
    let h_left = h.qstart >= SPLIT_SOFTCLIP;
    let h_right = h.qlen - h.qend >= SPLIT_SOFTCLIP;
    if !(t_left || t_right) || !(h_left || h_right) {
        return None;
    }
    let target_pos = if t_left { t.tstart } else { t.tend };
    let host_pos = if h_left { h.tstart } else { h.tend };
    Some(SplitEvent {
        target_contig: t.contig.clone(),
        target_pos,
        host_contig: h.contig.clone(),
        host_pos,
        direction: format!("host{}_{}", h.strand, t.strand),
        qname: qname.to_string(),
    })
}

#[derive(Debug, Clone, Default)]
pub struct FragmentEvidence {
    pub r1: ReadDecision,
    pub r2: ReadDecision,
    pub split_events: Vec<SplitEvent>,

    pub discordant: Option<String>,
}

pub fn adjudicate_fragment(key: &str, hits1: &[Hit], hits2: &[Hit]) -> FragmentEvidence {
    let r1 = adjudicate_read(hits1);
    let r2 = adjudicate_read(hits2);
    let mut split_events = Vec::new();
    split_events.extend(split_event(&r1, key));
    split_events.extend(split_event(&r2, key));
    let discordant = if r1.confident_target.is_some() && r2.confident_target.is_some() {
        None
    } else if let (Some(t), Some(h)) = (&r1.confident_target, &r2.best_host) {
        if h.mapq >= MIN_MAPQ {
            Some(t.contig.clone())
        } else {
            None
        }
    } else if let (Some(h), Some(t)) = (&r1.best_host, &r2.confident_target) {
        if h.mapq >= MIN_MAPQ {
            Some(t.contig.clone())
        } else {
            None
        }
    } else {
        None
    };
    FragmentEvidence {
        r1,
        r2,
        split_events,
        discordant,
    }
}

fn accept_ordered<R, G>(
    pending: &mut BTreeMap<usize, Vec<R>>,
    next: &mut usize,
    sink_err: &mut Option<String>,
    sink: &mut G,
    allow_sink: bool,
    idx: usize,
    results: Vec<R>,
) where
    G: FnMut(Vec<R>) -> Result<(), String>,
{
    pending.insert(idx, results);
    while let Some(ready) = pending.remove(next) {
        if allow_sink && sink_err.is_none() {
            if let Err(err) = sink(ready) {
                *sink_err = Some(err);
            }
        }
        *next += 1;
    }
}

/// Process chunks with bounded workers and deliver completed results in input order.
pub(crate) fn par_map_chunks<T, R, F, G>(
    items: impl Iterator<Item = Result<T, String>>,
    workers: usize,
    chunk_size: usize,
    process: F,
    mut sink: G,
) -> Result<(), String>
where
    T: Send,
    R: Send,
    F: Fn(&[T]) -> Vec<R> + Sync + Send,
    G: FnMut(Vec<R>) -> Result<(), String>,
{
    let chunk_size = chunk_size.max(1);

    if workers == 0 {
        let mut first_err: Option<String> = None;
        let mut chunk: Vec<T> = Vec::with_capacity(chunk_size);
        for item in items {
            match item {
                Ok(t) => {
                    chunk.push(t);
                    if chunk.len() >= chunk_size {
                        let results = process(&std::mem::take(&mut chunk));
                        if let Err(e) = sink(results) {
                            first_err = Some(e);
                            break;
                        }
                    }
                }
                Err(e) => {
                    first_err = Some(e);
                    break;
                }
            }
        }
        if !chunk.is_empty() && first_err.is_none() {
            if let Err(e) = sink(process(&chunk)) {
                first_err = Some(e);
            }
        }
        return if let Some(e) = first_err {
            Err(e)
        } else {
            Ok(())
        };
    }
    let (tx_chunk, rx_chunk) = sync_channel::<(usize, Vec<T>)>(workers * 2);
    let (tx_done, rx_done) = sync_channel::<(usize, Vec<R>)>(workers * 2);
    let rx_chunk = std::sync::Mutex::new(rx_chunk);
    let mut first_err: Option<String> = None;

    thread::scope(|scope| {
        for _ in 0..workers {
            let rx = &rx_chunk;

            let tx = tx_done.clone();
            let process = &process;
            scope.spawn(move || loop {
                let item = {
                    let guard = match rx.lock() {
                        Ok(g) => g,
                        Err(_) => return,
                    };
                    guard.recv()
                };
                let (idx, chunk) = match item {
                    Ok(x) => x,
                    Err(_) => return,
                };
                let results: Vec<R> = process(&chunk);
                if tx.send((idx, results)).is_err() {
                    return;
                }
            });
        }

        let mut idx = 0usize;
        let mut chunk: Vec<T> = Vec::with_capacity(chunk_size);
        let mut pending: BTreeMap<usize, Vec<R>> = BTreeMap::new();
        let mut next = 0usize;
        let mut sink_err: Option<String> = None;

        let drain = |pending: &mut BTreeMap<usize, Vec<R>>,
                     next: &mut usize,
                     sink_err: &mut Option<String>,
                     sink: &mut G| {
            while let Ok((i, results)) = rx_done.try_recv() {
                accept_ordered(pending, next, sink_err, sink, true, i, results);
            }
        };
        let submit = |i: usize,
                      chunk: Vec<T>,
                      pending: &mut BTreeMap<usize, Vec<R>>,
                      next: &mut usize,
                      sink_err: &mut Option<String>,
                      sink: &mut G| {
            match tx_chunk.try_send((i, chunk)) {
                Ok(()) => Ok(()),
                Err(TrySendError::Full((i, chunk))) => {
                    let results = process(&chunk);
                    accept_ordered(pending, next, sink_err, sink, true, i, results);
                    Ok(())
                }
                Err(TrySendError::Disconnected(_)) => {
                    Err("processing worker exited unexpectedly".to_string())
                }
            }
        };
        for item in items {
            match item {
                Ok(t) => {
                    chunk.push(t);
                    if chunk.len() >= chunk_size {
                        drain(&mut pending, &mut next, &mut sink_err, &mut sink);
                        let full = std::mem::replace(&mut chunk, Vec::with_capacity(chunk_size));
                        if let Err(err) =
                            submit(idx, full, &mut pending, &mut next, &mut sink_err, &mut sink)
                        {
                            first_err = Some(err);
                            break;
                        }
                        idx += 1;
                        drain(&mut pending, &mut next, &mut sink_err, &mut sink);
                    }
                }
                Err(e) => {
                    first_err = Some(e);
                    break;
                }
            }
        }
        if !chunk.is_empty() && first_err.is_none() {
            drain(&mut pending, &mut next, &mut sink_err, &mut sink);
            if let Err(err) = submit(
                idx,
                chunk,
                &mut pending,
                &mut next,
                &mut sink_err,
                &mut sink,
            ) {
                first_err = Some(err);
            }
            drain(&mut pending, &mut next, &mut sink_err, &mut sink);
        }
        drop(tx_chunk);

        drop(tx_done);

        while let Ok((i, results)) = rx_done.recv() {
            accept_ordered(
                &mut pending,
                &mut next,
                &mut sink_err,
                &mut sink,
                first_err.is_none(),
                i,
                results,
            );
        }
        if let Some(e) = first_err.take() {
            Err(e)
        } else if let Some(e) = sink_err.take() {
            Err(e)
        } else {
            Ok(())
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn hit(
        role: Role,
        as_score: i32,
        mapq: u32,
        qstart: i32,
        qend: i32,
        tstart: i32,
        tend: i32,
    ) -> Hit {
        Hit {
            role,
            contig: format!("{}_{}", role.prefix(), 0),
            as_score,
            mapq,
            is_primary: true,
            is_supplementary: false,
            nm: 0,
            qstart,
            qend,
            qlen: 150,
            tstart,
            tend,
            strand: '+',
        }
    }

    fn named_hit(role: Role, contig: &str, as_score: i32, mapq: u32, nm: i32) -> Hit {
        let mut hit = hit(role, as_score, mapq, 0, 150, 0, 150);
        hit.contig = contig.to_string();
        hit.nm = nm;
        hit
    }

    #[test]
    fn adjudicates_confident_target_with_as_margin() {
        let hits = vec![
            hit(Role::Target, 300, 60, 0, 150, 0, 150),
            hit(Role::Host, 280, 60, 0, 150, 0, 150),
        ];
        let d = adjudicate_read(&hits);
        assert!(d.confident_target.is_some());
    }

    #[test]
    fn rejects_target_without_as_margin() {
        let hits = vec![
            hit(Role::Target, 300, 60, 0, 150, 0, 150),
            hit(Role::Host, 295, 60, 0, 150, 0, 150),
        ];
        let d = adjudicate_read(&hits);
        assert!(d.confident_target.is_none());
    }

    #[test]
    fn audit_recovers_compatible_mapq_zero_targets_over_weak_non_target() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 300, 0, 1),
            named_hit(Role::Host, "host_0", 280, 0, 0),
        ];

        match adjudicate_audit(&hits) {
            AuditDecision::Resolved { role, hits } => {
                assert_eq!(role, Role::Target);
                assert_eq!(
                    hits.iter()
                        .map(|hit| hit.contig.as_str())
                        .collect::<Vec<_>>(),
                    ["target_a", "target_b"]
                );
            }
            decision => panic!("target should be recovered from ALL_CHAINS; got {decision:?}"),
        }
    }

    #[test]
    fn audit_excludes_non_top_target_until_delta_as_is_calibrated() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 299, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&hits) else {
            panic!("best target should be recovered");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(
            hits.len(),
            1,
            "the delta-AS>0 window must remain disabled until calibration"
        );
        assert_eq!(hits[0].contig, "target_a");
    }

    #[test]
    fn audit_tied_host_or_near_contaminant_blocks_target() {
        for (competitor, competitor_score) in [(Role::Host, 300), (Role::Contaminant, 289)] {
            let hits = vec![
                named_hit(Role::Target, "target_a", 300, 0, 0),
                named_hit(competitor, "non_target", competitor_score, 0, 0),
            ];
            assert!(matches!(adjudicate_audit(&hits), AuditDecision::NoEvidence));
        }
    }

    #[test]
    fn audit_returns_only_targets_hit_directly_by_this_read() {
        let prior_read = vec![
            named_hit(Role::Target, "target_b", 300, 0, 0),
            named_hit(Role::Target, "target_c", 300, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        assert!(matches!(
            adjudicate_audit(&prior_read),
            AuditDecision::Resolved { .. }
        ));

        let this_read = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 300, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&this_read) else {
            panic!("direct A and B hits should be recovered");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(
            hits.iter()
                .map(|hit| hit.contig.as_str())
                .collect::<Vec<_>>(),
            ["target_a", "target_b"]
        );
    }

    #[test]
    fn audit_deduplicates_contigs_and_keeps_the_best_hit() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 294, 0, 0),
            named_hit(Role::Target, "target_a", 300, 0, 2),
            named_hit(Role::Target, "target_b", 300, 0, 1),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&hits) else {
            panic!("direct target hit should be recovered");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].contig, "target_a");
        assert_eq!(hits[0].as_score, 300);
        assert_eq!(hits[1].contig, "target_b");
    }

    #[test]
    fn audit_decoy_resolution_is_symmetric_and_conservative() {
        let winning = vec![
            named_hit(Role::Decoy, "decoy_a", 300, 0, 0),
            named_hit(Role::Decoy, "decoy_b", 300, 0, 1),
            named_hit(Role::Decoy, "decoy_c", 299, 0, 0),
            named_hit(Role::Host, "host_0", 280, 0, 0),
            named_hit(Role::Target, "target_a", 279, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&winning) else {
            panic!("a clearly winning decoy should produce symmetric evidence");
        };
        assert_eq!(role, Role::Decoy);
        assert_eq!(
            hits.iter()
                .map(|hit| hit.contig.as_str())
                .collect::<Vec<_>>(),
            ["decoy_a", "decoy_b"],
            "the decoy delta-AS>0 window must also remain disabled until calibration"
        );

        let near_tie = vec![
            named_hit(Role::Decoy, "decoy_a", 300, 0, 0),
            named_hit(Role::Contaminant, "contam_0", 289, 0, 0),
        ];
        assert!(matches!(
            adjudicate_audit(&near_tie),
            AuditDecision::NoEvidence
        ));
    }

    #[test]
    fn audit_overflow_is_explicit() {
        let hits = (0..=MAX_AUDIT_HITS)
            .map(|i| named_hit(Role::Target, &format!("target_{i}"), 300, 0, 0))
            .collect::<Vec<_>>();
        assert_eq!(
            adjudicate_audit(&hits),
            AuditDecision::Overflow {
                hit_count: MAX_AUDIT_HITS + 1,
                limit: MAX_AUDIT_HITS,
            }
        );
    }

    #[test]
    fn best_n_saturation_requests_audit_before_a_hidden_host_can_be_missed() {
        let mut ordinary_hits = vec![named_hit(Role::Target, "target_primary", 300, 60, 0)];
        ordinary_hits[0].is_primary = true;
        for i in 0..BEST_N {
            let mut secondary = named_hit(
                Role::Target,
                &format!("target_secondary_{i}"),
                299 - i,
                0,
                0,
            );
            secondary.is_primary = false;
            ordinary_hits.push(secondary);
        }
        let ordinary = adjudicate_read(&ordinary_hits);
        assert!(ordinary.confident_target.is_some());
        assert!(needs_ambiguity_audit(&ordinary_hits, &ordinary));

        let mut all_chains = ordinary_hits;
        all_chains.push(named_hit(Role::Target, "target_after_best_n", 280, 0, 0));
        all_chains.push(named_hit(Role::Host, "host_after_best_n", 295, 0, 0));
        assert!(matches!(
            adjudicate_audit(&all_chains),
            AuditDecision::NoEvidence
        ));
    }

    #[test]
    fn audit_trigger_covers_low_mapq_primary_and_unresolved_target_but_not_clear_host() {
        let mut low_mapq_primary = named_hit(Role::Host, "host_0", 300, MIN_MAPQ - 1, 0);
        low_mapq_primary.is_primary = true;
        let hits = vec![low_mapq_primary];
        assert!(needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));

        let mut unresolved_target = named_hit(Role::Target, "target_a", 300, 0, 0);
        unresolved_target.is_primary = false;
        let hits = vec![unresolved_target];
        assert!(needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));

        let mut clear_host = named_hit(Role::Host, "host_0", 300, 60, 0);
        clear_host.is_primary = true;
        let hits = vec![clear_host];
        assert!(!needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));
    }

    #[test]
    fn hits_preserve_primary_and_supplementary_flags() {
        let mapping = Mapping {
            target_name: Some(Arc::new("target_a".to_string())),
            is_primary: false,
            is_supplementary: true,
            alignment: Some(minimap2::Alignment {
                nm: 0,
                cigar: None,
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: Some(300),
            }),
            ..Mapping::default()
        };
        let roles = HashMap::from([("target_a".to_string(), Role::Target)]);
        let hits = hits_of(&[mapping], &roles);
        assert_eq!(hits.len(), 1);
        assert!(!hits[0].is_primary);
        assert!(hits[0].is_supplementary);
    }

    #[test]
    fn split_event_requires_both_side_softclips() {
        let hits1 = vec![
            hit(Role::Target, 200, 40, 0, 100, 0, 100),
            hit(Role::Host, 180, 40, 100, 150, 500, 550),
        ];
        let hits2 = vec![
            hit(Role::Target, 200, 40, 0, 150, 0, 150),
            hit(Role::Host, 180, 40, 100, 150, 500, 550),
        ];
        let r1 = adjudicate_read(&hits1);
        let r2 = adjudicate_read(&hits2);
        assert!(split_event(&r1, "q1").is_some());
        assert!(split_event(&r2, "q2").is_none());
    }

    #[test]
    fn parallel_map_is_deterministic() {
        let items: Vec<Result<u64, String>> = (0..1000).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            4,
            64,
            |chunk| chunk.iter().map(|&v| v * 2).collect(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(out.len(), 1000);
        for (i, v) in out.iter().enumerate() {
            assert_eq!(*v, (i * 2) as u64);
        }
    }

    #[test]
    fn saturated_queue_uses_caller_and_preserves_order() {
        let caller = std::thread::current().id();
        let caller_used = std::sync::atomic::AtomicBool::new(false);
        let items: Vec<Result<u64, String>> = (0..32).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            1,
            1,
            |chunk| {
                let value = chunk[0];
                if std::thread::current().id() == caller {
                    caller_used.store(true, std::sync::atomic::Ordering::Release);
                } else if value == 0 {
                    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(1);
                    while !caller_used.load(std::sync::atomic::Ordering::Acquire) {
                        if std::time::Instant::now() >= deadline {
                            break;
                        }
                        std::thread::yield_now();
                    }
                }
                vec![value]
            },
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert!(caller_used.load(std::sync::atomic::Ordering::Acquire));
        assert_eq!(out, (0..32).collect::<Vec<_>>());
    }

    #[test]
    fn parallel_map_propagates_stream_error() {
        let mut items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        items.push(Err("stream error".into()));
        let mut out = Vec::new();
        let result = par_map_chunks(
            items.into_iter(),
            2,
            4,
            |chunk| chunk.to_vec(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        );
        assert_eq!(result.unwrap_err(), "stream error");
    }

    #[test]
    fn parallel_map_propagates_sink_error() {
        let items: Vec<Result<u64, String>> = (0..100).map(Ok).collect();
        let mut calls = 0u64;
        let result = par_map_chunks(
            items.into_iter(),
            4,
            8,
            |chunk| chunk.to_vec(),
            |rs| {
                calls += rs.len() as u64;
                if calls > 30 {
                    Err("sink failure".into())
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result.unwrap_err(), "sink failure");
        assert_eq!(
            calls, 32,
            "no later chunks may be delivered after sink failure"
        );
    }

    #[test]
    fn inline_map_zero_workers_is_deterministic() {
        let items: Vec<Result<u64, String>> = (0..1000).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            0,
            64,
            |chunk| chunk.iter().map(|&v| v * 2).collect(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(out.len(), 1000);
        for (i, v) in out.iter().enumerate() {
            assert_eq!(*v, (i * 2) as u64);
        }
    }

    #[test]
    fn inline_map_propagates_errors() {
        let items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        assert!(
            par_map_chunks(items.into_iter(), 0, 4, |chunk| chunk.to_vec(), |_| Ok(()),).is_ok()
        );
        let mut items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        items.push(Err("stream error".into()));
        assert!(
            par_map_chunks(items.into_iter(), 0, 4, |chunk| chunk.to_vec(), |_| Ok(()),).is_err()
        );
        let items: Vec<Result<u64, String>> = (0..100).map(Ok).collect();
        let mut calls = 0u64;
        assert!(par_map_chunks(
            items.into_iter(),
            0,
            8,
            |chunk| chunk.to_vec(),
            |rs| {
                calls += rs.len() as u64;
                if calls > 30 {
                    Err("sink failure".into())
                } else {
                    Ok(())
                }
            },
        )
        .is_err());
    }

    #[test]
    fn open_rejects_multi_part_index() {
        let dir = std::env::temp_dir().join(format!("vf_multipart_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("multi.fa");
        let mmi = dir.join("multi.mmi");
        let mut file = std::fs::File::create(&fasta).unwrap();
        let seq = "ACGT".repeat(200);
        for i in 0..4 {
            writeln!(file, ">contig_{i}\n{seq}").unwrap();
        }
        drop(file);

        let mut builder = Aligner::builder().sr().with_index_threads(1);
        builder.idxopt.mini_batch_size = 1;
        builder.idxopt.batch_size = 1;
        let indexed = builder
            .with_index(&fasta, Some(mmi.to_str().unwrap()))
            .unwrap();
        assert!(indexed.idx_parts.len() > 1);
        drop(indexed);

        let err = match CompetitiveAligner::open(&mmi, 1) {
            Ok(_) => panic!("a multi-shard index must be rejected"),
            Err(err) => err,
        };
        assert!(err.contains("single-shard"), "err={err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn open_shares_one_index_with_the_all_chains_audit_aligner() {
        let dir = std::env::temp_dir().join(format!("vf_audit_aligner_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("single.fa");
        let mmi = dir.join("single.mmi");
        let mut file = std::fs::File::create(&fasta).unwrap();
        writeln!(file, ">target_0\n{}", "ACGT".repeat(200)).unwrap();
        drop(file);

        let indexed = Aligner::builder()
            .sr()
            .with_index_threads(1)
            .with_index(&fasta, Some(mmi.to_str().unwrap()))
            .unwrap();
        assert_eq!(indexed.idx_parts.len(), 1);
        drop(indexed);

        let aligners = CompetitiveAligner::open(&mmi, 1).unwrap();
        let all_chains = minimap2::ffi::MM_F_ALL_CHAINS as i64;
        assert_eq!(aligners.aligner.mapopt.best_n, BEST_N);
        assert_eq!(aligners.aligner.mapopt.flag & all_chains, 0);
        assert_ne!(aligners.audit_aligner.mapopt.flag & all_chains, 0);
        assert!(Arc::ptr_eq(
            &aligners.aligner.idx_parts[0],
            &aligners.audit_aligner.idx_parts[0]
        ));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
