use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use crate::fastq::Fragment;
use crate::index::{ReferenceContig, ReferenceRole};
use minimap2::{Aligner, Built, Mapping};

mod workers;
pub(crate) use workers::{align_fragments_bounded, AnalysisWorkerConfig};

pub struct CompetitiveAligner {
    aligner: Aligner<Built>,
}

impl CompetitiveAligner {
    pub(crate) fn short_read_chain_requirements(
        expected_kmer_length: usize,
    ) -> Result<(usize, usize), String> {
        let aligner = Aligner::builder().sr();
        let kmer_length = usize::try_from(aligner.idxopt.k)
            .map_err(|_| "minimap2 short-read k-mer length is invalid".to_string())?;
        if kmer_length != expected_kmer_length {
            return Err(format!(
                "minimap2 short-read k-mer length {kmer_length} does not match target Bloom k-mer length {expected_kmer_length}"
            ));
        }
        let minimum_hits = usize::try_from(aligner.mapopt.min_cnt)
            .map_err(|_| "minimap2 short-read minimum chain count is invalid".to_string())?;
        let minimum_covered_bases = usize::try_from(aligner.mapopt.min_chain_score)
            .map_err(|_| "minimap2 short-read minimum chain score is invalid".to_string())?;
        Ok((minimum_hits, minimum_covered_bases))
    }

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
        fragment: &Fragment<'_>,
        contigs: &std::collections::HashMap<String, ReferenceContig>,
    ) -> Result<FragmentAlignmentEvidence, String> {
        let r1 = self
            .aligner
            .map(
                fragment.r1,
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

#[cfg(test)]
mod tests;
