use std::collections::{BTreeMap, BTreeSet};

use serde::Serialize;

use crate::alignment::{FragmentAdjudication, FragmentAlignmentEvidence};
use crate::index::reference::ReferenceGroup;
use crate::integration_evidence::IntegrationStatus;

pub const INTERVAL_METHOD: &str = "EQUAL_TAILED_EXACT_HYPERGEOMETRIC_INVERSION";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
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
    intervals: IntervalUnion,
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
            intervals: IntervalUnion::default(),
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
        self.intervals.covered_bases()
    }
    pub fn occupied_windows(&self, bins: usize) -> u64 {
        let mut occupied = BTreeSet::new();
        for (&start, &end) in &self.intervals.segments {
            let end = end.min(self.group.representative_length);
            if start >= end {
                continue;
            }
            let first = start * bins as u64 / self.group.representative_length;
            let last =
                ((end - 1) * bins as u64 / self.group.representative_length).min(bins as u64 - 1);
            occupied.extend(first..=last);
        }
        occupied.len() as u64
    }
    pub fn integration_status(&self) -> IntegrationStatus {
        if self.supporting_selected_fragments > 0
            && (self.split_events > 0 || self.discordant_fragments > 0)
        {
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
                target.split_events += u64::from(evidence.split_groups.contains(group));
                for &(start, end) in evidence.target_intervals.get(group).into_iter().flatten() {
                    target
                        .intervals
                        .insert(start, end.min(target.group.representative_length));
                }
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
    if sample > population
        || successes > sample
        || population == 0
        || sample == 0
        || !(0.0..1.0).contains(&alpha)
    {
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
            hypergeometric_tail_reaches(
                population,
                total_successes,
                sample,
                successes,
                Tail::Upper,
                tail,
                false,
            )
        })?
        .saturating_sub(1)
        .max(feasible_low)
    };
    let upper = if successes == sample {
        population
    } else {
        last_true(feasible_low, feasible_high, |total_successes| {
            hypergeometric_tail_reaches(
                population,
                total_successes,
                sample,
                successes,
                Tail::Lower,
                tail,
                true,
            )
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

#[derive(Clone, Copy)]
enum Tail {
    Lower,
    Upper,
}

#[derive(Default)]
struct CompensatedSum {
    sum: f64,
    correction: f64,
}

impl CompensatedSum {
    fn add(&mut self, value: f64) {
        let adjusted = value - self.correction;
        let next = self.sum + adjusted;
        self.correction = (next - self.sum) - adjusted;
        self.sum = next;
    }
}

fn ratio(numerator: u128, denominator: u128) -> Result<f64, String> {
    if numerator == 0 || denominator == 0 {
        return Err("Invalid hypergeometric recurrence".into());
    }
    Ok(numerator as f64 / denominator as f64)
}

fn negligible_remainder(weight: f64, next_ratio: f64, accumulated: f64) -> bool {
    next_ratio < 1.0
        && weight * next_ratio / (1.0 - next_ratio) <= f64::EPSILON * accumulated.max(1.0)
}

fn hypergeometric_tail_reaches(
    population: u64,
    total_successes: u64,
    sample: u64,
    observed: u64,
    tail: Tail,
    threshold: f64,
    inclusive: bool,
) -> Result<bool, String> {
    let support_low = sample.saturating_sub(population - total_successes);
    let support_high = sample.min(total_successes);
    match tail {
        Tail::Upper if observed <= support_low => return Ok(true),
        Tail::Upper if observed > support_high => return Ok(false),
        Tail::Lower if observed < support_low => return Ok(false),
        Tail::Lower if observed >= support_high => return Ok(true),
        _ => {}
    }

    let mode = (((u128::from(sample) + 1) * (u128::from(total_successes) + 1))
        / (u128::from(population) + 2)) as u64;
    let mode = mode.clamp(support_low, support_high);
    let mut total_weight = CompensatedSum::default();
    let mut tail_weight = CompensatedSum::default();
    total_weight.add(1.0);
    if matches!(tail, Tail::Upper) && mode >= observed
        || matches!(tail, Tail::Lower) && mode <= observed
    {
        tail_weight.add(1.0);
    }

    let mut weight = 1.0;
    let mut current = mode;
    while current > support_low {
        let step_ratio = ratio(
            u128::from(current) * u128::from(population - total_successes - (sample - current)),
            u128::from(total_successes - current + 1) * u128::from(sample - current + 1),
        )?;
        weight *= step_ratio;
        current -= 1;
        total_weight.add(weight);
        if matches!(tail, Tail::Lower) && current <= observed {
            tail_weight.add(weight);
        }
        if current == support_low {
            break;
        }
        let next_ratio = ratio(
            u128::from(current) * u128::from(population - total_successes - (sample - current)),
            u128::from(total_successes - current + 1) * u128::from(sample - current + 1),
        )?;
        if negligible_remainder(weight, next_ratio, total_weight.sum) {
            break;
        }
    }

    weight = 1.0;
    current = mode;
    while current < support_high {
        let step_ratio = ratio(
            u128::from(total_successes - current) * u128::from(sample - current),
            u128::from(current + 1)
                * u128::from(population - total_successes - (sample - current) + 1),
        )?;
        weight *= step_ratio;
        current += 1;
        total_weight.add(weight);
        if matches!(tail, Tail::Upper) && current >= observed {
            tail_weight.add(weight);
        }
        if current == support_high {
            break;
        }
        let next_ratio = ratio(
            u128::from(total_successes - current) * u128::from(sample - current),
            u128::from(current + 1)
                * u128::from(population - total_successes - (sample - current) + 1),
        )?;
        if negligible_remainder(weight, next_ratio, total_weight.sum) {
            break;
        }
    }

    let scaled_threshold = threshold * total_weight.sum;
    Ok(if inclusive {
        tail_weight.sum >= scaled_threshold
    } else {
        tail_weight.sum > scaled_threshold
    })
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

#[derive(Debug, Clone, Default)]
struct IntervalUnion {
    segments: BTreeMap<u64, u64>,
}

impl IntervalUnion {
    fn insert(&mut self, mut start: u64, mut end: u64) {
        if start >= end {
            return;
        }
        if let Some((&previous_start, &previous_end)) = self.segments.range(..=start).next_back() {
            if previous_end >= start {
                start = previous_start;
                end = end.max(previous_end);
                self.segments.remove(&previous_start);
            }
        }
        loop {
            let next = self
                .segments
                .range(start..)
                .next()
                .map(|(&left, &right)| (left, right));
            match next {
                Some((next_start, next_end)) if next_start <= end => {
                    end = end.max(next_end);
                    self.segments.remove(&next_start);
                }
                _ => break,
            }
        }
        self.segments.insert(start, end);
    }

    fn covered_bases(&self) -> u64 {
        self.segments.iter().map(|(start, end)| end - start).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use statrs::distribution::{DiscreteCDF, Hypergeometric};

    fn choose(total: u64, count: u64) -> u128 {
        let count = count.min(total - count);
        (1..=count).fold(1_u128, |value, step| {
            value * u128::from(total - count + step) / u128::from(step)
        })
    }

    fn hypergeometric_count(
        population: u64,
        total_successes: u64,
        sample: u64,
        observed: u64,
    ) -> u128 {
        if observed > total_successes
            || observed > sample
            || sample - observed > population - total_successes
        {
            return 0;
        }
        choose(total_successes, observed) * choose(population - total_successes, sample - observed)
    }

    fn exact_cdf_count(observed: u64, population: u64, total_successes: u64, sample: u64) -> u128 {
        (0..=observed)
            .map(|candidate| hypergeometric_count(population, total_successes, sample, candidate))
            .sum()
    }

    fn exact_interval_endpoints(
        population: u64,
        sample: u64,
        successes: u64,
        family_size: u64,
    ) -> (u64, u64) {
        if sample == population {
            return (successes, successes);
        }
        let denominator = choose(population, sample);
        let tail_denominator = 40 * u128::from(family_size);
        let feasible_low = successes;
        let feasible_high = population - (sample - successes);
        let lower = if successes == 0 {
            0
        } else {
            let mut candidate = feasible_low;
            while candidate < feasible_high
                && exact_cdf_count(successes - 1, population, candidate, sample) * tail_denominator
                    >= denominator * (tail_denominator - 1)
            {
                candidate += 1;
            }
            candidate.saturating_sub(1).max(feasible_low)
        };
        let upper = if successes == sample {
            population
        } else {
            let mut candidate = feasible_high;
            while candidate > feasible_low
                && exact_cdf_count(successes, population, candidate, sample) * tail_denominator
                    < denominator
            {
                candidate -= 1;
            }
            candidate.saturating_add(1).min(feasible_high)
        };
        (lower, upper)
    }

    fn group() -> ReferenceGroup {
        ReferenceGroup {
            ordinal: 0,
            target_group_id: "sha256:group".into(),
            representative_id: "target".into(),
            member_ids: vec!["target".into()],
            representative_length: 100,
            contig_name: "target_0".into(),
        }
    }

    #[test]
    fn exact_interval_matches_small_population_enumeration() {
        let interval = finite_population_interval(50, 10, 3, 0.05).unwrap();
        assert_eq!((interval.lower * 50.0).round() as u64, 4);
        assert_eq!((interval.upper * 50.0).round() as u64, 32);
        for population in 2..=30 {
            for sample in 1..population {
                for successes in 0..=sample {
                    let actual =
                        finite_population_interval(population, sample, successes, 0.05).unwrap();
                    let (lower, upper) = exact_interval_endpoints(population, sample, successes, 1);
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
    fn stable_tail_recurrence_fixes_real_lower_endpoint() {
        let population = 35_801_278;
        let sample = 1_268_602;
        let successes = 51;
        let alpha = 0.05 / 20_560.0;
        let interval = finite_population_interval(population, sample, successes, alpha).unwrap();
        assert_eq!(
            (
                (interval.lower * population as f64).round() as u64,
                (interval.upper * population as f64).round() as u64,
            ),
            (692, 2612)
        );

        let tail = alpha / 2.0;
        assert!(hypergeometric_tail_reaches(
            population,
            693,
            sample,
            successes,
            Tail::Upper,
            tail,
            false,
        )
        .unwrap());
        let old_distribution = Hypergeometric::new(population, 693, sample).unwrap();
        assert!(old_distribution.cdf(successes - 1) >= 1.0 - tail);
    }

    #[test]
    fn interval_endpoints_are_bounded_and_monotone() {
        for &(population, sample, alpha) in &[(101, 17, 0.05), (10_003, 997, 0.05 / 20_560.0)] {
            let mut previous = (0, 0);
            for successes in 0..=sample {
                let interval =
                    finite_population_interval(population, sample, successes, alpha).unwrap();
                let endpoints = (
                    (interval.lower * population as f64).round() as u64,
                    (interval.upper * population as f64).round() as u64,
                );
                let feasible_high = population - (sample - successes);
                assert!(successes <= endpoints.0);
                assert!(endpoints.0 <= endpoints.1);
                assert!(endpoints.1 <= feasible_high);
                assert!(previous.0 <= endpoints.0);
                assert!(previous.1 <= endpoints.1);
                previous = endpoints;
            }
        }
    }

    #[test]
    fn phase5_production_endpoints_bias_and_family_coverage_match_exact_enumeration() {
        let interval = finite_population_interval(50, 10, 3, 0.05).unwrap();
        assert_eq!(
            (
                (interval.lower * 50.0).round() as u64,
                (interval.upper * 50.0).round() as u64
            ),
            exact_interval_endpoints(50, 10, 3, 1)
        );
        assert_eq!(exact_interval_endpoints(50, 10, 3, 1), (4, 32));

        for family_size in [1_u64, 20_560] {
            for population in 2..=20 {
                for sample in 1..=population {
                    let denominator = choose(population, sample);
                    for total_successes in 0..=population {
                        let weighted_sum = (0..=sample)
                            .map(|successes| {
                                u128::from(successes)
                                    * hypergeometric_count(
                                        population,
                                        total_successes,
                                        sample,
                                        successes,
                                    )
                            })
                            .sum::<u128>();
                        assert_eq!(
                            weighted_sum * u128::from(population),
                            u128::from(total_successes) * u128::from(sample) * denominator
                        );

                        let noncoverage = (0..=sample)
                            .filter(|successes| {
                                let interval = finite_population_interval(
                                    population,
                                    sample,
                                    *successes,
                                    0.05 / family_size as f64,
                                )
                                .unwrap();
                                let lower = (interval.lower * population as f64).round() as u64;
                                let upper = (interval.upper * population as f64).round() as u64;
                                total_successes < lower || total_successes > upper
                            })
                            .map(|successes| {
                                hypergeometric_count(population, total_successes, sample, successes)
                            })
                            .sum::<u128>();
                        assert!(
                            noncoverage * 20 * u128::from(family_size) <= denominator,
                            "N={population} n={sample} M={total_successes} m={family_size}"
                        );
                    }
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

    #[test]
    fn coverage_state_is_bounded_by_merged_reference_intervals() {
        let mut intervals = IntervalUnion::default();
        for offset in 0..100_000 {
            intervals.insert(offset % 50, 100 + offset % 50);
        }
        assert_eq!(intervals.segments.len(), 1);
        assert_eq!(intervals.covered_bases(), 149);
    }

    #[test]
    fn host_confounded_fragment_cannot_create_integration_evidence() {
        let mut accumulator = EvidenceAccumulator::new(&[group()]);
        accumulator.accumulate_group_evidence(FragmentAlignmentEvidence {
            adjudication: FragmentAdjudication::Confounded(BTreeSet::from([0])),
            aligned: true,
            target_intervals: BTreeMap::new(),
            split_groups: BTreeSet::from([0]),
            discordant_groups: BTreeSet::new(),
        });
        assert_eq!(accumulator.groups[0].host_confounded_fragments, 1);
        assert_eq!(accumulator.groups[0].supporting_selected_fragments, 0);
        assert_eq!(accumulator.groups[0].split_events, 0);
        assert_eq!(
            accumulator.groups[0].integration_status(),
            IntegrationStatus::NotObserved
        );
    }

    #[test]
    fn supporting_host_target_split_creates_integration_evidence_once() {
        let mut accumulator = EvidenceAccumulator::new(&[group()]);
        accumulator.accumulate_group_evidence(FragmentAlignmentEvidence {
            adjudication: FragmentAdjudication::Supporting(0),
            aligned: true,
            target_intervals: BTreeMap::from([(0, vec![(10, 40)])]),
            split_groups: BTreeSet::from([0]),
            discordant_groups: BTreeSet::new(),
        });
        assert_eq!(accumulator.groups[0].supporting_selected_fragments, 1);
        assert_eq!(accumulator.groups[0].split_events, 1);
        assert_eq!(
            accumulator.groups[0].integration_status(),
            IntegrationStatus::DiagnosticEvidenceObserved
        );
    }
}
