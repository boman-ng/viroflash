use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;
use std::sync::{mpsc::sync_channel, Arc};

use minimap2::{Aligner, Built, Mapping, Strand};

use crate::fastq::Fragment;
use crate::index::{ReferenceContig, ReferenceRole};

const ALIGNMENT_QUEUE_FRAGMENTS_PER_THREAD: usize = 1;

#[cfg(test)]
#[derive(Debug, Clone, Copy, Default)]
struct AlignmentRetention {
    pending_fragments: usize,
    pending_sequence_bytes: usize,
    maximum_pending_fragments: usize,
    maximum_pending_sequence_bytes: usize,
    maximum_fragment_bytes: usize,
    total_sequence_bytes: usize,
}

#[cfg(test)]
impl AlignmentRetention {
    fn claim(&mut self, fragment_bytes: usize, threads: usize) {
        self.pending_fragments += 1;
        self.pending_sequence_bytes += fragment_bytes;
        self.maximum_pending_fragments = self.maximum_pending_fragments.max(self.pending_fragments);
        self.maximum_pending_sequence_bytes = self
            .maximum_pending_sequence_bytes
            .max(self.pending_sequence_bytes);
        self.maximum_fragment_bytes = self.maximum_fragment_bytes.max(fragment_bytes);
        self.total_sequence_bytes += fragment_bytes;
        assert!(self.pending_fragments <= threads);
        assert!(self.pending_sequence_bytes <= threads.saturating_mul(self.maximum_fragment_bytes));
    }

    fn release(&mut self, fragment_bytes: usize) {
        self.pending_fragments -= 1;
        self.pending_sequence_bytes -= fragment_bytes;
    }
}

struct CompletedAlignment {
    worker_index: usize,
    #[cfg(test)]
    fragment_bytes: usize,
    result: Result<FragmentAlignmentEvidence, String>,
}

pub struct CompetitiveAligner {
    aligner: Aligner<Built>,
}

impl CompetitiveAligner {
    pub fn open(path: &Path) -> Result<Self, String> {
        let mut aligner = Aligner::builder()
            .sr()
            .with_cigar()
            .with_index(path, None)
            .map_err(|error| format!("Cannot open minimap2 index: {error}"))?;
        if aligner.idx_parts.len() != 1 {
            return Err("minimap2 index must contain exactly one part".into());
        }
        aligner.mapopt.flag |= minimap2::ffi::MM_F_ALL_CHAINS as i64;
        Ok(Self { aligner })
    }

    pub fn align_fragment_competitively(
        &self,
        fragment: &Fragment,
        contigs: &std::collections::HashMap<String, ReferenceContig>,
    ) -> Result<FragmentAlignmentEvidence, String> {
        let r1 = self
            .aligner
            .map(
                &fragment.r1,
                false,
                false,
                None,
                None,
                Some(fragment.id.as_bytes()),
            )
            .map_err(|error| format!("Failed to align fragment {}: {error}", fragment.id))?;
        let r2 = fragment
            .r2
            .as_ref()
            .map(|sequence| {
                self.aligner.map(
                    sequence,
                    false,
                    false,
                    None,
                    None,
                    Some(fragment.id.as_bytes()),
                )
            })
            .transpose()
            .map_err(|error| format!("Failed to align fragment {}: {error}", fragment.id))?;
        Ok(adjudicate_fragment(
            &hits_of(&r1, contigs),
            r2.as_deref()
                .map(|mappings| hits_of(mappings, contigs))
                .as_deref()
                .unwrap_or(&[]),
        ))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AlignmentHit {
    pub role: ReferenceRole,
    pub target_group_ordinal: Option<usize>,
    pub alignment_score: i32,
    pub query_length: u32,
    pub query_start: u32,
    pub query_end: u32,
    pub target_start: u64,
    pub target_end: u64,
    pub forward: bool,
    pub supplementary: bool,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum FragmentAdjudication {
    Supporting(usize),
    Unresolved(BTreeSet<usize>),
    Confounded(BTreeSet<usize>),
    NoTargetEvidence,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FragmentAlignmentEvidence {
    pub adjudication: FragmentAdjudication,
    pub aligned: bool,
    pub target_intervals: BTreeMap<usize, Vec<(u64, u64)>>,
    pub split_groups: BTreeSet<usize>,
    pub discordant_groups: BTreeSet<usize>,
}

fn hits_of(
    mappings: &[Mapping],
    contigs: &std::collections::HashMap<String, ReferenceContig>,
) -> Vec<AlignmentHit> {
    mappings
        .iter()
        .filter_map(|mapping| {
            let name = mapping.target_name.as_ref()?;
            let contig = contigs.get(name.as_ref())?;
            let alignment = mapping.alignment.as_ref()?;
            let query_length = mapping.query_len?.get();
            if query_length == 0 || alignment.cigar.as_ref().is_none_or(Vec::is_empty) {
                return None;
            }
            Some(AlignmentHit {
                role: contig.role,
                target_group_ordinal: contig.target_group_ordinal,
                alignment_score: alignment.alignment_score.unwrap_or(0),
                query_length: query_length as u32,
                query_start: mapping.query_start.max(0) as u32,
                query_end: mapping.query_end.max(0) as u32,
                target_start: mapping.target_start.max(0) as u64,
                target_end: mapping.target_end.max(0) as u64,
                forward: matches!(mapping.strand, Strand::Forward),
                supplementary: mapping.is_supplementary,
            })
        })
        .collect()
}

fn normalized_score(left: &AlignmentHit, right: &AlignmentHit) -> Ordering {
    (i128::from(left.alignment_score) * i128::from(right.query_length))
        .cmp(&(i128::from(right.alignment_score) * i128::from(left.query_length)))
}

#[derive(Default)]
struct EndEvidence {
    targets: BTreeSet<usize>,
    target_hits: BTreeMap<usize, Vec<AlignmentHit>>,
    split_groups: BTreeSet<usize>,
    has_host: bool,
    host_confounded: bool,
}

fn adjudicate_end(hits: &[AlignmentHit]) -> EndEvidence {
    let best_host = hits
        .iter()
        .filter(|hit| hit.role == ReferenceRole::Host)
        .max_by(|left, right| normalized_score(left, right));
    let best_target = hits
        .iter()
        .filter(|hit| hit.role == ReferenceRole::Target)
        .max_by(|left, right| normalized_score(left, right));
    let Some(best_target) = best_target else {
        return EndEvidence {
            has_host: best_host.is_some(),
            ..EndEvidence::default()
        };
    };
    let host_confounded =
        best_host.is_some_and(|host| normalized_score(host, best_target) != Ordering::Less);
    let mut result = EndEvidence {
        has_host: best_host.is_some(),
        host_confounded,
        ..EndEvidence::default()
    };
    for hit in hits.iter().filter(|hit| {
        hit.role == ReferenceRole::Target && normalized_score(hit, best_target) == Ordering::Equal
    }) {
        let ordinal = hit.target_group_ordinal.expect("target contig has a group");
        result.targets.insert(ordinal);
        result
            .target_hits
            .entry(ordinal)
            .or_default()
            .push(hit.clone());
    }
    result.split_groups = host_target_split_groups(hits, &result.targets);
    result
}

pub fn adjudicate_fragment(r1: &[AlignmentHit], r2: &[AlignmentHit]) -> FragmentAlignmentEvidence {
    let left = adjudicate_end(r1);
    let right = adjudicate_end(r2);
    let mut target_intervals: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    let split_groups = left
        .split_groups
        .union(&right.split_groups)
        .copied()
        .collect();
    for end in [&left, &right] {
        for (group, hits) in &end.target_hits {
            target_intervals
                .entry(*group)
                .or_default()
                .extend(hits.iter().map(|hit| {
                    (
                        hit.target_start.min(hit.target_end),
                        hit.target_start.max(hit.target_end),
                    )
                }));
        }
    }
    let all_targets = left
        .targets
        .union(&right.targets)
        .copied()
        .collect::<BTreeSet<_>>();
    let adjudication = if all_targets.is_empty() {
        FragmentAdjudication::NoTargetEvidence
    } else if left.host_confounded || right.host_confounded {
        FragmentAdjudication::Confounded(all_targets.clone())
    } else {
        let surviving = if left.targets.is_empty() {
            right.targets.clone()
        } else if right.targets.is_empty() {
            left.targets.clone()
        } else {
            left.targets.intersection(&right.targets).copied().collect()
        };
        if surviving.len() == 1 {
            FragmentAdjudication::Supporting(*surviving.first().unwrap())
        } else {
            FragmentAdjudication::Unresolved(if surviving.is_empty() {
                all_targets.clone()
            } else {
                surviving
            })
        }
    };
    let discordant_groups = match &adjudication {
        FragmentAdjudication::Supporting(group)
            if (left.targets.contains(group) && right.targets.is_empty() && right.has_host)
                || (right.targets.contains(group) && left.targets.is_empty() && left.has_host) =>
        {
            BTreeSet::from([*group])
        }
        _ => BTreeSet::new(),
    };
    FragmentAlignmentEvidence {
        adjudication,
        aligned: !r1.is_empty() || !r2.is_empty(),
        target_intervals,
        split_groups,
        discordant_groups,
    }
}

fn host_target_split_groups(
    hits: &[AlignmentHit],
    target_groups: &BTreeSet<usize>,
) -> BTreeSet<usize> {
    let mut split_groups = BTreeSet::new();
    for target in hits.iter().filter(|hit| hit.role == ReferenceRole::Target) {
        let Some(group) = target.target_group_ordinal else {
            continue;
        };
        if !target_groups.contains(&group) {
            continue;
        }
        if hits.iter().any(|host| {
            host.role == ReferenceRole::Host
                && (target.supplementary || host.supplementary)
                && (target.query_end <= host.query_start || host.query_end <= target.query_start)
        }) {
            split_groups.insert(group);
        }
    }
    split_groups
}

pub fn align_fragments_bounded<N, S>(
    index_path: &Path,
    contigs: &std::collections::HashMap<String, ReferenceContig>,
    threads: usize,
    mut next_fragment: N,
    mut sink: S,
) -> Result<(), String>
where
    N: FnMut() -> Result<Option<Fragment>, String>,
    S: FnMut(FragmentAlignmentEvidence),
{
    align_fragments_bounded_inner(
        index_path,
        contigs,
        threads,
        &mut next_fragment,
        &mut sink,
        #[cfg(test)]
        |_| {},
    )
}

fn align_fragments_bounded_inner<N, S>(
    index_path: &Path,
    contigs: &std::collections::HashMap<String, ReferenceContig>,
    threads: usize,
    next_fragment: &mut N,
    sink: &mut S,
    #[cfg(test)] mut observe_retention: impl FnMut(AlignmentRetention),
) -> Result<(), String>
where
    N: FnMut() -> Result<Option<Fragment>, String>,
    S: FnMut(FragmentAlignmentEvidence),
{
    let aligner = Arc::new(CompetitiveAligner::open(index_path)?);
    let queue_capacity = threads * ALIGNMENT_QUEUE_FRAGMENTS_PER_THREAD;
    let (result_sender, result_receiver) = sync_channel(queue_capacity);
    std::thread::scope(|scope| -> Result<(), String> {
        let mut handles = Vec::new();
        let mut task_senders = Vec::with_capacity(threads);
        for worker_index in 0..threads {
            let (task_sender, task_receiver) =
                sync_channel::<Fragment>(ALIGNMENT_QUEUE_FRAGMENTS_PER_THREAD);
            task_senders.push(task_sender);
            let aligner = Arc::clone(&aligner);
            let result_sender = result_sender.clone();
            handles.push(scope.spawn(move || loop {
                let task = task_receiver.recv();
                let Ok(fragment) = task else { break };
                #[cfg(test)]
                let fragment_bytes =
                    fragment.r1.len() + fragment.r2.as_ref().map_or(0, |sequence| sequence.len());
                let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                    aligner.align_fragment_competitively(&fragment, contigs)
                }));
                drop(fragment);
                let panicked = result.is_err();
                let result = result.unwrap_or_else(|_| Err("Alignment worker panicked".into()));
                if result_sender
                    .send(CompletedAlignment {
                        worker_index,
                        #[cfg(test)]
                        fragment_bytes,
                        result,
                    })
                    .is_err()
                {
                    break;
                }
                if panicked {
                    break;
                }
            }));
        }
        drop(result_sender);
        let processing = (|| -> Result<(), String> {
            #[cfg(test)]
            let mut retention = AlignmentRetention::default();
            let mut active = 0;
            for sender in &task_senders {
                let Some(fragment) = next_fragment()? else {
                    break;
                };
                #[cfg(test)]
                {
                    let fragment_bytes = fragment.r1.len()
                        + fragment.r2.as_ref().map_or(0, |sequence| sequence.len());
                    retention.claim(fragment_bytes, threads);
                    observe_retention(retention);
                }
                sender
                    .send(fragment)
                    .map_err(|_| "Alignment worker queue closed".to_string())?;
                active += 1;
            }
            while active > 0 {
                let completed = result_receiver
                    .recv()
                    .map_err(|_| "Alignment worker result queue closed".to_string())?;
                active -= 1;
                #[cfg(test)]
                {
                    retention.release(completed.fragment_bytes);
                    observe_retention(retention);
                }
                sink(completed.result?);
                if let Some(fragment) = next_fragment()? {
                    #[cfg(test)]
                    {
                        let fragment_bytes = fragment.r1.len()
                            + fragment.r2.as_ref().map_or(0, |sequence| sequence.len());
                        retention.claim(fragment_bytes, threads);
                        observe_retention(retention);
                    }
                    task_senders[completed.worker_index]
                        .send(fragment)
                        .map_err(|_| "Alignment worker queue closed".to_string())?;
                    active += 1;
                }
            }
            #[cfg(test)]
            {
                assert_eq!(retention.pending_fragments, 0);
                assert_eq!(retention.pending_sequence_bytes, 0);
            }
            Ok(())
        })();
        drop(task_senders);
        let mut join_error = None;
        for handle in handles {
            if handle.join().is_err() {
                join_error = Some("Alignment worker panicked".to_string());
            }
        }
        processing.and(join_error.map_or(Ok(()), Err))
    })
}

#[cfg(test)]
mod tests {
    use std::cell::RefCell;
    use std::collections::VecDeque;
    use std::fs::File;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

    use super::*;
    use crate::index::{build_index, IndexOptions};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn retention_within_bound(stats: AlignmentRetention, threads: usize) -> bool {
        stats.maximum_pending_fragments <= threads
            && stats.maximum_pending_sequence_bytes
                <= threads.saturating_mul(stats.maximum_fragment_bytes)
    }

    fn measured_alignment(
        index_path: &Path,
        contigs: &std::collections::HashMap<String, ReferenceContig>,
        threads: usize,
        fragments: Vec<Fragment>,
    ) -> AlignmentRetention {
        let stats = RefCell::new(AlignmentRetention::default());
        let mut fragments = VecDeque::from(fragments);
        align_fragments_bounded_inner(
            index_path,
            contigs,
            threads,
            &mut || Ok(fragments.pop_front()),
            &mut |_| {},
            |observation| {
                *stats.borrow_mut() = observation;
            },
        )
        .unwrap();
        stats.into_inner()
    }

    fn fixture() -> (std::path::PathBuf, crate::index::ReferenceIndex) {
        let root = std::env::temp_dir().join(format!(
            "viroflash-phase5-retention-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, AtomicOrdering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        for (name, base) in [("host", b'A'), ("target", b'C')] {
            let mut file = File::create(root.join(format!("{name}.fa"))).unwrap();
            writeln!(
                file,
                ">{name}\n{}",
                char::from(base).to_string().repeat(8_192)
            )
            .unwrap();
        }
        let index = build_index(&IndexOptions {
            host_fa: root.join("host.fa"),
            target_fa: root.join("target.fa"),
            out_dir: root.join("index"),
            threads: 1,
        })
        .unwrap();
        (root, index)
    }

    #[test]
    fn completed_worker_is_refilled_before_batch_drains() {
        let (root, index) = fixture();
        let mut fragments = (0..6)
            .map(|ordinal| Fragment {
                ordinal,
                id: format!("fragment-{ordinal}"),
                r1: vec![b'G'; 120],
                r2: Some(vec![b'T'; 120]),
            })
            .collect::<VecDeque<_>>();
        let trace = RefCell::new(Vec::new());
        align_fragments_bounded_inner(
            &index.mmi_path,
            &index.contigs,
            2,
            &mut || Ok(fragments.pop_front()),
            &mut |_| {},
            |retention| trace.borrow_mut().push(retention.pending_fragments),
        )
        .unwrap();

        assert!(trace
            .into_inner()
            .windows(3)
            .any(|window| window == [2, 1, 2]));
        let _ = std::fs::remove_dir_all(root);
    }

    #[test]
    fn phase5_alignment_retention_is_threads_times_max_fragment_not_total_selected() {
        let (root, index) = fixture();
        for read_length in [40, 120, 4_096] {
            for threads in [1, 2, 4, 8] {
                let fragments = |count| {
                    (0..count)
                        .map(|ordinal| Fragment {
                            ordinal,
                            id: format!("fragment-{ordinal}"),
                            r1: vec![b'G'; read_length],
                            r2: Some(vec![b'T'; read_length]),
                        })
                        .collect::<Vec<_>>()
                };
                let small = measured_alignment(
                    &index.mmi_path,
                    &index.contigs,
                    threads,
                    fragments(threads as u64 * 2),
                );
                let large = measured_alignment(
                    &index.mmi_path,
                    &index.contigs,
                    threads,
                    fragments(threads as u64 * 20),
                );
                assert!(retention_within_bound(small, threads));
                assert!(retention_within_bound(large, threads));
                assert_eq!(
                    small.maximum_pending_sequence_bytes,
                    large.maximum_pending_sequence_bytes
                );
                assert!(large.total_sequence_bytes > small.total_sequence_bytes);

                let counterfactual = AlignmentRetention {
                    maximum_pending_sequence_bytes: large.total_sequence_bytes,
                    ..large
                };
                assert!(!retention_within_bound(counterfactual, threads));
            }
        }
        let _ = std::fs::remove_dir_all(root);
    }

    fn target(group: usize, score: i32) -> AlignmentHit {
        AlignmentHit {
            role: ReferenceRole::Target,
            target_group_ordinal: Some(group),
            alignment_score: score,
            query_length: 100,
            query_start: 0,
            query_end: 89,
            target_start: 1,
            target_end: 90,
            forward: true,
            supplementary: false,
        }
    }
    fn host(score: i32) -> AlignmentHit {
        AlignmentHit {
            role: ReferenceRole::Host,
            target_group_ordinal: None,
            alignment_score: score,
            query_length: 100,
            query_start: 0,
            query_end: 89,
            target_start: 1,
            target_end: 90,
            forward: true,
            supplementary: false,
        }
    }

    #[test]
    fn host_tie_or_advantage_confounds_fragment() {
        assert!(matches!(
            adjudicate_fragment(&[target(0, 90), host(90)], &[]).adjudication,
            FragmentAdjudication::Confounded(_)
        ));
        assert!(matches!(
            adjudicate_fragment(&[target(0, 90), host(91)], &[]).adjudication,
            FragmentAdjudication::Confounded(_)
        ));
    }

    #[test]
    fn within_group_resolves_but_cross_group_tie_does_not() {
        assert_eq!(
            adjudicate_fragment(&[target(0, 90), target(0, 90)], &[]).adjudication,
            FragmentAdjudication::Supporting(0)
        );
        assert!(
            matches!(adjudicate_fragment(&[target(0, 90), target(1, 90)], &[]).adjudication, FragmentAdjudication::Unresolved(groups) if groups == BTreeSet::from([0, 1]))
        );
    }

    #[test]
    fn paired_disjoint_groups_are_unresolved_once() {
        let evidence = adjudicate_fragment(&[target(0, 90)], &[target(1, 90)]);
        assert!(matches!(
            evidence.adjudication,
            FragmentAdjudication::Unresolved(_)
        ));
        assert!(evidence.discordant_groups.is_empty());
    }

    #[test]
    fn target_and_host_mates_support_group_with_discordance_diagnostic() {
        let evidence = adjudicate_fragment(&[target(0, 90)], &[host(90)]);
        assert_eq!(evidence.adjudication, FragmentAdjudication::Supporting(0));
        assert_eq!(evidence.discordant_groups, BTreeSet::from([0]));
    }

    #[test]
    fn split_requires_host_target_supplementary_disjoint_geometry() {
        let primary = AlignmentHit {
            query_start: 0,
            query_end: 40,
            ..target(0, 90)
        };
        let target_supplementary = AlignmentHit {
            alignment_score: 35,
            query_start: 60,
            query_end: 100,
            supplementary: true,
            ..target(0, 90)
        };
        let host_supplementary = AlignmentHit {
            query_start: 60,
            query_end: 100,
            supplementary: true,
            ..host(35)
        };
        assert!(
            adjudicate_fragment(&[primary.clone(), target_supplementary], &[])
                .split_groups
                .is_empty()
        );
        assert_eq!(
            adjudicate_fragment(&[primary, host_supplementary], &[]).split_groups,
            BTreeSet::from([0])
        );
    }

    #[test]
    fn host_target_alternatives_without_supplementary_geometry_are_not_split() {
        let target_hit = target(0, 90);
        let overlapping_host = host(35);
        let disjoint_secondary_host = AlignmentHit {
            query_start: 90,
            query_end: 100,
            ..host(35)
        };
        assert!(
            adjudicate_fragment(&[target_hit.clone(), overlapping_host], &[])
                .split_groups
                .is_empty()
        );
        assert!(
            adjudicate_fragment(&[target_hit, disjoint_secondary_host], &[])
                .split_groups
                .is_empty()
        );
    }
}
