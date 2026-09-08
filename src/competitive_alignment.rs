use std::cmp::Ordering;
use std::collections::{BTreeMap, BTreeSet};
use std::path::Path;

use minimap2::{Aligner, Built, Mapping, Strand};

use crate::fastq_input::Fragment;
use crate::reference_index::{ReferenceContig, ReferenceRole};

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
    pub target_start: u64,
    pub target_end: u64,
    pub forward: bool,
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
                target_start: mapping.target_start.max(0) as u64,
                target_end: mapping.target_end.max(0) as u64,
                forward: matches!(mapping.strand, Strand::Forward),
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
    result
}

pub fn adjudicate_fragment(r1: &[AlignmentHit], r2: &[AlignmentHit]) -> FragmentAlignmentEvidence {
    let left = adjudicate_end(r1);
    let right = adjudicate_end(r2);
    let mut target_intervals: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    let mut split_groups = BTreeSet::new();
    for end in [&left, &right] {
        for (group, hits) in &end.target_hits {
            if hits.len() > 1 {
                split_groups.insert(*group);
            }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn target(group: usize, score: i32) -> AlignmentHit {
        AlignmentHit {
            role: ReferenceRole::Target,
            target_group_ordinal: Some(group),
            alignment_score: score,
            query_length: 100,
            target_start: 1,
            target_end: 90,
            forward: true,
        }
    }
    fn host(score: i32) -> AlignmentHit {
        AlignmentHit {
            role: ReferenceRole::Host,
            target_group_ordinal: None,
            alignment_score: score,
            query_length: 100,
            target_start: 1,
            target_end: 90,
            forward: true,
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
}
