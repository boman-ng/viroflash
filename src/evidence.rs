use std::collections::BTreeSet;

use statrs::distribution::{DiscreteCDF, Hypergeometric};

use crate::competitive_alignment::{FragmentAdjudication, FragmentAlignmentEvidence};
use crate::integration_evidence::IntegrationStatus;
use crate::reference_group::ReferenceGroup;

pub const INTERVAL_METHOD: &str = "EQUAL_TAILED_EXACT_HYPERGEOMETRIC_INVERSION";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceStatus {
    ReferenceSignalObserved,
    IndeterminateEvidence,
}
impl EvidenceStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ReferenceSignalObserved => "REFERENCE_SIGNAL_OBSERVED",
            Self::IndeterminateEvidence => "INDETERMINATE_EVIDENCE",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AttributionStatus {
    ResolvedToReferenceGroup,
    AmbiguousWithinGroup,
    UnresolvedAcrossGroups,
    ConfoundedWithHost,
}
impl AttributionStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ResolvedToReferenceGroup => "RESOLVED_TO_REFERENCE_GROUP",
            Self::AmbiguousWithinGroup => "AMBIGUOUS_WITHIN_GROUP",
            Self::UnresolvedAcrossGroups => "UNRESOLVED_ACROSS_GROUPS",
            Self::ConfoundedWithHost => "CONFOUNDED_WITH_HOST",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FinitePopulationInterval {
    pub lower: f64,
    pub upper: f64,
    pub level: f64,
}

#[derive(Debug, Clone)]
pub struct TargetGroupEvidence {
    pub group: ReferenceGroup,
    pub supporting_selected_fragments: u64,
    pub host_confounded_fragments: u64,
    pub cross_group_ambiguous_fragments: u64,
    intervals: Vec<(u64, u64)>,
    pub split_events: u64,
    pub discordant_fragments: u64,
}

impl TargetGroupEvidence {
    fn new(group: ReferenceGroup) -> Self {
        Self {
            group,
            supporting_selected_fragments: 0,
            host_confounded_fragments: 0,
            cross_group_ambiguous_fragments: 0,
            intervals: Vec::new(),
            split_events: 0,
            discordant_fragments: 0,
        }
    }
    pub fn observed_or_indeterminate(&self) -> bool {
        self.supporting_selected_fragments > 0
            || self.host_confounded_fragments > 0
            || self.cross_group_ambiguous_fragments > 0
    }
    pub fn evidence_status(&self) -> EvidenceStatus {
        if self.supporting_selected_fragments > 0 {
            EvidenceStatus::ReferenceSignalObserved
        } else {
            EvidenceStatus::IndeterminateEvidence
        }
    }
    pub fn attribution_status(&self) -> AttributionStatus {
        if self.supporting_selected_fragments > 0 {
            if self.group.member_ids.len() == 1 {
                AttributionStatus::ResolvedToReferenceGroup
            } else {
                AttributionStatus::AmbiguousWithinGroup
            }
        } else if self.cross_group_ambiguous_fragments > 0 {
            AttributionStatus::UnresolvedAcrossGroups
        } else {
            AttributionStatus::ConfoundedWithHost
        }
    }
    pub fn covered_bases(&self) -> u64 {
        merged_covered_bases(&self.intervals)
    }
    pub fn occupied_windows(&self, bins: usize) -> u64 {
        let mut occupied = BTreeSet::new();
        for &(start, end) in &self.intervals {
            for position in start..end.min(self.group.representative_length) {
                occupied.insert(
                    (position * bins as u64 / self.group.representative_length)
                        .min(bins as u64 - 1),
                );
            }
        }
        occupied.len() as u64
    }
    pub fn integration_status(&self) -> IntegrationStatus {
        if self.split_events > 0 || self.discordant_fragments > 0 {
            IntegrationStatus::DiagnosticEvidenceObserved
        } else {
            IntegrationStatus::NotObserved
        }
    }
}

pub struct EvidenceAccumulator {
    pub groups: Vec<TargetGroupEvidence>,
    pub aligned_fragments: u64,
    pub unassigned_fragments: u64,
}

impl EvidenceAccumulator {
    pub fn new(groups: &[ReferenceGroup]) -> Self {
        Self {
            groups: groups
                .iter()
                .cloned()
                .map(TargetGroupEvidence::new)
                .collect(),
            aligned_fragments: 0,
            unassigned_fragments: 0,
        }
    }
    pub fn accumulate_group_evidence(&mut self, evidence: FragmentAlignmentEvidence) {
        if evidence.aligned {
            self.aligned_fragments += 1;
        }
        match &evidence.adjudication {
            FragmentAdjudication::Supporting(group) => {
                let target = &mut self.groups[*group];
                target.supporting_selected_fragments += 1;
                target.intervals.extend(
                    evidence
                        .target_intervals
                        .get(group)
                        .into_iter()
                        .flatten()
                        .copied(),
                );
                target.split_events += u64::from(evidence.split_groups.contains(group));
                target.discordant_fragments +=
                    u64::from(evidence.discordant_groups.contains(group));
            }
            FragmentAdjudication::Unresolved(groups) => {
                self.unassigned_fragments += u64::from(evidence.aligned);
                for group in groups {
                    self.groups[*group].cross_group_ambiguous_fragments += 1;
                }
            }
            FragmentAdjudication::Confounded(groups) => {
                self.unassigned_fragments += u64::from(evidence.aligned);
                for group in groups {
                    self.groups[*group].host_confounded_fragments += 1;
                }
            }
            FragmentAdjudication::NoTargetEvidence => {
                self.unassigned_fragments += u64::from(evidence.aligned);
            }
        }
    }
}

pub fn finite_population_interval(
    population: u64,
    sample: u64,
    successes: u64,
    alpha: f64,
) -> Result<FinitePopulationInterval, String> {
    if sample > population || successes > sample || population == 0 || sample == 0 {
        return Err("Invalid finite-population interval inputs".into());
    }
    if sample == population {
        let exact = successes as f64 / population as f64;
        return Ok(FinitePopulationInterval {
            lower: exact,
            upper: exact,
            level: 1.0 - alpha,
        });
    }
    let feasible_low = successes;
    let feasible_high = population - (sample - successes);
    let tail = alpha / 2.0;
    let lower = if successes == 0 {
        0
    } else {
        first_true(feasible_low, feasible_high, |total_successes| {
            Ok(hypergeom_cdf(successes - 1, population, total_successes, sample)? < 1.0 - tail)
        })?
        .saturating_sub(1)
        .max(feasible_low)
    };
    let upper = if successes == sample {
        population
    } else {
        last_true(feasible_low, feasible_high, |total_successes| {
            Ok(hypergeom_cdf(successes, population, total_successes, sample)? >= tail)
        })?
        .saturating_add(1)
        .min(feasible_high)
    };
    Ok(FinitePopulationInterval {
        lower: lower as f64 / population as f64,
        upper: upper as f64 / population as f64,
        level: 1.0 - alpha,
    })
}

fn hypergeom_cdf(
    x: u64,
    population: u64,
    total_successes: u64,
    sample: u64,
) -> Result<f64, String> {
    Hypergeometric::new(population, total_successes, sample)
        .map(|distribution| distribution.cdf(x))
        .map_err(|error| format!("Invalid hypergeometric distribution: {error}"))
}

fn first_true(
    mut low: u64,
    mut high: u64,
    predicate: impl Fn(u64) -> Result<bool, String>,
) -> Result<u64, String> {
    while low < high {
        let middle = low + (high - low) / 2;
        if predicate(middle)? {
            high = middle;
        } else {
            low = middle + 1;
        }
    }
    Ok(low)
}

fn last_true(
    mut low: u64,
    mut high: u64,
    predicate: impl Fn(u64) -> Result<bool, String>,
) -> Result<u64, String> {
    while low < high {
        let middle = low + (high - low).div_ceil(2);
        if predicate(middle)? {
            low = middle;
        } else {
            high = middle - 1;
        }
    }
    Ok(low)
}

fn merged_covered_bases(intervals: &[(u64, u64)]) -> u64 {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut total = 0;
    let mut current: Option<(u64, u64)> = None;
    for (start, end) in sorted {
        match current {
            Some((left, right)) if start <= right => current = Some((left, right.max(end))),
            Some((left, right)) => {
                total += right - left;
                current = Some((start, end));
            }
            None => current = Some((start, end)),
        }
    }
    if let Some((left, right)) = current {
        total += right - left;
    }
    total
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn exact_interval_matches_small_population_enumeration() {
        let interval = finite_population_interval(50, 10, 3, 0.05).unwrap();
        assert_eq!((interval.lower * 50.0).round() as u64, 4);
        assert_eq!((interval.upper * 50.0).round() as u64, 32);
        for population in 2..=30 {
            for sample in 1..population {
                for successes in 0..=sample {
                    let actual =
                        finite_population_interval(population, sample, successes, 0.1).unwrap();
                    let lower = if successes == 0 {
                        0
                    } else {
                        let mut candidate = successes;
                        while candidate < population - (sample - successes)
                            && hypergeom_cdf(successes - 1, population, candidate, sample).unwrap()
                                >= 0.95
                        {
                            candidate += 1;
                        }
                        candidate.saturating_sub(1).max(successes)
                    };
                    let upper = if successes == sample {
                        population
                    } else {
                        let mut candidate = population - (sample - successes);
                        while candidate > successes
                            && hypergeom_cdf(successes, population, candidate, sample).unwrap()
                                < 0.05
                        {
                            candidate -= 1;
                        }
                        candidate
                            .saturating_add(1)
                            .min(population - (sample - successes))
                    };
                    assert_eq!(
                        (actual.lower * population as f64).round() as u64,
                        lower,
                        "N={population} n={sample} x={successes}"
                    );
                    assert_eq!(
                        (actual.upper * population as f64).round() as u64,
                        upper,
                        "N={population} n={sample} x={successes}"
                    );
                }
            }
        }
    }

    #[test]
    fn census_interval_is_exact() {
        assert_eq!(
            finite_population_interval(10, 10, 3, 0.05).unwrap(),
            FinitePopulationInterval {
                lower: 0.3,
                upper: 0.3,
                level: 0.95
            }
        );
    }
}
