//! Compound hypotheses derived from direct evidence equivalence classes.
//!
//! Each discovery evidence unit contributes a target-member hyperedge. Hyperedges sharing any
//! member belong to one connected component. A component's `members` union forms an OR hypothesis;
//! `explanation` is only a compact description and does not narrow attribution.

use std::collections::{BTreeMap, BTreeSet};

/// Compound hypothesis built from direct discovery evidence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositeHypothesis {
    /// All directly observed component members, sorted by ID; collectively an OR hypothesis.
    pub members: Vec<usize>,
    /// Compact weighted-greedy hitting-set approximation, sorted by ID.
    pub explanation: Vec<usize>,
    /// Number of non-empty discovery evidence units in the component. Callers currently use one
    /// read end per unit; every input unit counts once regardless of member count.
    pub discovery_observations: u64,
}

/// Build disjoint compound hypotheses from discovery target-member sets.
///
/// Members are sorted and deduplicated per observation; empty sets are ignored. Distinct
/// observations retain separate votes even when their normalized sets match. Results are sorted
/// lexicographically by `members`, independent of observation and member order.
pub fn build_composite_hypotheses(
    discovery_fragment_members: &[Vec<usize>],
) -> Vec<CompositeHypothesis> {
    let fragments = canonical_fragments(discovery_fragment_members);
    if fragments.is_empty() {
        return Vec::new();
    }

    // The inverted index encodes hypergraph adjacency directly. Expand member postings and then
    // their remaining members without materializing every pair in a large hyperedge.
    let mut member_fragments = BTreeMap::<usize, Vec<usize>>::new();
    for (fragment_index, members) in fragments.iter().enumerate() {
        for &member in members {
            member_fragments
                .entry(member)
                .or_default()
                .push(fragment_index);
        }
    }

    let mut visited_members = BTreeSet::new();
    let mut visited_fragments = vec![false; fragments.len()];
    let mut hypotheses = Vec::new();

    for &start_member in member_fragments.keys() {
        if visited_members.contains(&start_member) {
            continue;
        }

        let mut pending_members = vec![start_member];
        let mut members = Vec::new();
        let mut component_fragment_indices = Vec::new();

        while let Some(member) = pending_members.pop() {
            if !visited_members.insert(member) {
                continue;
            }
            members.push(member);

            if let Some(fragment_indices) = member_fragments.get(&member) {
                for &fragment_index in fragment_indices {
                    // One hyperedge occurs in multiple postings but belongs to and counts once.
                    if visited_fragments[fragment_index] {
                        continue;
                    }
                    visited_fragments[fragment_index] = true;
                    component_fragment_indices.push(fragment_index);
                    for &linked_member in &fragments[fragment_index] {
                        if !visited_members.contains(&linked_member) {
                            pending_members.push(linked_member);
                        }
                    }
                }
            }
        }

        members.sort_unstable();
        component_fragment_indices.sort_unstable();
        let component_fragments: Vec<&[usize]> = component_fragment_indices
            .iter()
            .map(|&index| fragments[index].as_slice())
            .collect();
        hypotheses.push(CompositeHypothesis {
            members,
            explanation: greedy_explanation(&component_fragments),
            discovery_observations: component_fragment_indices.len() as u64,
        });
    }

    hypotheses.sort_by(|left, right| left.members.cmp(&right.members));
    hypotheses
}

/// Return the discovered hypothesis index containing a validation evidence unit.
///
/// The validation set is sorted and deduplicated. It must be non-empty and every member must have
/// been discovered in one hypothesis. Unknown or cross-hypothesis sets return `None` and never
/// expand an existing hypothesis.
pub fn validation_hypothesis_index(
    hypotheses: &[CompositeHypothesis],
    validation_members: &[usize],
) -> Option<usize> {
    let mut canonical_members = validation_members.to_vec();
    canonical_members.sort_unstable();
    canonical_members.dedup();
    let first_member = *canonical_members.first()?;

    let hypothesis_index = hypotheses
        .iter()
        .position(|hypothesis| hypothesis.members.binary_search(&first_member).is_ok())?;
    canonical_members
        .iter()
        .all(|member| {
            hypotheses[hypothesis_index]
                .members
                .binary_search(member)
                .is_ok()
        })
        .then_some(hypothesis_index)
}

fn canonical_fragments(discovery_fragment_members: &[Vec<usize>]) -> Vec<Vec<usize>> {
    let mut fragments = Vec::with_capacity(discovery_fragment_members.len());
    for raw_members in discovery_fragment_members {
        let mut members = raw_members.clone();
        members.sort_unstable();
        members.dedup();
        if !members.is_empty() {
            fragments.push(members);
        }
    }
    // Do not deduplicate equal sets across fragments; every input fragment retains one vote.
    fragments.sort_unstable();
    fragments
}

/// Deterministic weighted-greedy minimum-hitting-set approximation for equal-weight fragments.
///
/// Each round chooses the member covering the most uncovered sets and breaks ties by lower ID.
/// The final values are ID-sorted because they form an explanation set, not a ranking.
fn greedy_explanation(fragment_sets: &[&[usize]]) -> Vec<usize> {
    let mut uncovered = vec![true; fragment_sets.len()];
    let mut remaining = fragment_sets.len();
    let mut explanation = Vec::new();

    while remaining > 0 {
        let mut votes = BTreeMap::<usize, usize>::new();
        for (is_uncovered, members) in uncovered.iter().zip(fragment_sets) {
            if !is_uncovered {
                continue;
            }
            // canonical_fragments already deduplicates members, so a fragment votes once per member.
            for &member in *members {
                *votes.entry(member).or_default() += 1;
            }
        }

        let mut best_member = usize::MAX;
        let mut best_votes = 0;
        for (member, member_votes) in votes {
            if member_votes > best_votes || (member_votes == best_votes && member < best_member) {
                best_member = member;
                best_votes = member_votes;
            }
        }
        // Inputs are non-empty; this guard prevents future private changes from looping forever.
        if best_votes == 0 {
            break;
        }

        explanation.push(best_member);
        for (is_uncovered, members) in uncovered.iter_mut().zip(fragment_sets) {
            if *is_uncovered && members.binary_search(&best_member).is_ok() {
                *is_uncovered = false;
                remaining -= 1;
            }
        }
    }

    explanation.sort_unstable();
    explanation
}

#[cfg(test)]
mod tests {
    use super::*;

    fn hypothesis(
        members: &[usize],
        explanation: &[usize],
        discovery_observations: u64,
    ) -> CompositeHypothesis {
        CompositeHypothesis {
            members: members.to_vec(),
            explanation: explanation.to_vec(),
            discovery_observations,
        }
    }

    #[test]
    fn overlapping_hyperedges_form_one_or_hypothesis() {
        let hypotheses = build_composite_hypotheses(&[vec![1, 2], vec![2, 3]]);

        assert_eq!(hypotheses, vec![hypothesis(&[1, 2, 3], &[2], 2)]);
    }

    #[test]
    fn disjoint_hyperedges_form_disjoint_hypotheses() {
        let hypotheses =
            build_composite_hypotheses(&[vec![8, 7], vec![2, 1], vec![3, 2], vec![9, 8]]);

        assert_eq!(
            hypotheses,
            vec![
                hypothesis(&[1, 2, 3], &[2], 2),
                hypothesis(&[7, 8, 9], &[8], 2),
            ]
        );
    }

    #[test]
    fn empty_duplicates_and_input_order_do_not_change_result() {
        let first = vec![
            vec![3, 2, 2],
            Vec::new(),
            vec![1, 2, 1],
            vec![3, 2, 3],
            vec![9, 9],
            Vec::new(),
        ];
        let reordered = vec![
            Vec::new(),
            vec![9],
            vec![2, 3, 2],
            vec![2, 1, 1],
            Vec::new(),
            vec![2, 2, 3],
        ];

        let expected = vec![hypothesis(&[1, 2, 3], &[2], 3), hypothesis(&[9], &[9], 1)];
        assert_eq!(build_composite_hypotheses(&first), expected);
        assert_eq!(build_composite_hypotheses(&reordered), expected);
        assert!(build_composite_hypotheses(&[Vec::new(), Vec::new()]).is_empty());
    }

    #[test]
    fn validation_accepts_subsets_and_rejects_unknown_or_cross_hypothesis_sets() {
        let hypotheses = build_composite_hypotheses(&[vec![1, 2], vec![2, 3], vec![10, 11]]);

        assert_eq!(validation_hypothesis_index(&hypotheses, &[1]), Some(0));
        assert_eq!(
            validation_hypothesis_index(&hypotheses, &[3, 1, 3]),
            Some(0)
        );
        assert_eq!(validation_hypothesis_index(&hypotheses, &[11, 10]), Some(1));
        assert_eq!(validation_hypothesis_index(&hypotheses, &[]), None);
        assert_eq!(validation_hypothesis_index(&hypotheses, &[99]), None);
        assert_eq!(validation_hypothesis_index(&hypotheses, &[1, 99]), None);
        assert_eq!(validation_hypothesis_index(&hypotheses, &[2, 10]), None);
    }

    #[test]
    fn one_observation_is_counted_and_weighted_once() {
        let hypotheses = build_composite_hypotheses(&[vec![4, 1, 4, 2, 3, 1]]);

        // One hyperedge is reachable through four postings but remains one evidence unit. Each
        // normalized member gets one vote and the tie resolves to the lowest ID, 1.
        assert_eq!(hypotheses, vec![hypothesis(&[1, 2, 3, 4], &[1], 1)]);
    }

    #[test]
    fn greedy_explanation_recomputes_votes_after_each_selection() {
        let hypotheses =
            build_composite_hypotheses(&[vec![1, 2], vec![1, 3], vec![2, 4], vec![2, 5]]);

        assert_eq!(hypotheses, vec![hypothesis(&[1, 2, 3, 4, 5], &[1, 2], 4)]);
    }
}
