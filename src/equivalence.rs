//! 直接证据等价类的复合假设。
//!
//! discovery 证据单元的目标成员集合是超边：共享任一成员的超边属于同一连通块。
//! 每个块的 `members` 是直接出现成员的并集，语义为 OR hypothesis，不把证据归因
//! 到任一单独成员；`explanation` 只是该块的简约解释，不改变假设边界。

use std::collections::{BTreeMap, BTreeSet};

/// 由 discovery 直接证据构成的复合假设。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CompositeHypothesis {
    /// 连通块内直接出现的全部成员，按 ID 升序排列；整体表示 OR hypothesis。
    pub members: Vec<usize>,
    /// weighted greedy hitting-set 近似得到的简约解释，按 ID 升序排列。
    pub explanation: Vec<usize>,
    /// 归入该连通块的非空 discovery 证据单元数。调用方当前以 read-end 为
    /// 单元；每个输入单元恰好计一次，与命中成员数无关。
    pub discovery_observations: u64,
}

/// 从 discovery 证据单元的目标成员集合构造互不重叠的复合假设。
///
/// 每个观测会先在函数内按成员 ID 排序、去重；空集合被忽略。不同观测即使
/// 规范化后集合相同，仍各占一票。返回结果按 `members` 字典序稳定排列，因此
/// 不依赖观测顺序或观测内成员顺序。
pub fn build_composite_hypotheses(
    discovery_fragment_members: &[Vec<usize>],
) -> Vec<CompositeHypothesis> {
    let fragments = canonical_fragments(discovery_fragment_members);
    if fragments.is_empty() {
        return Vec::new();
    }

    // 倒排表直接表达超图关联。遍历某个成员时展开包含它的 fragment，再从该
    // fragment 展开其余成员；无需把一条大超边物化为全部成员对。
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
                    // 同一超边会出现在多个成员的 posting 中，但只能归属、计数一次。
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

/// 返回 validation 证据单元所属的已发现假设下标。
///
/// validation 集合在函数内排序、去重。集合必须非空，且每个成员都已由 discovery
/// 发现并属于同一个假设；未知成员或跨假设集合返回 `None`，不会扩张已有假设。
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
    // 不对相同集合做跨 fragment 去重：每个输入 fragment 都保留一票。
    fragments.sort_unstable();
    fragments
}

/// 等权 fragment 的 deterministic weighted greedy minimum-hitting-set 近似。
///
/// 每轮选择覆盖最多未覆盖集合的成员；票数相同时选择较小 ID。返回值最终按 ID
/// 排序，因为它是解释集合而不是成员排名。
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
            // members 已在 canonical_fragments 中去重，所以一个 fragment 对同一
            // member 至多投一票。
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
        // 所有传入集合都非空；该分支只保护私有函数未来改动时不产生死循环。
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

        // 单条超边虽可从四个 posting 到达，证据单元数仍为 1；规范化后四个成员
        // 各得一票，tie 选择最小 ID 1。
        assert_eq!(hypotheses, vec![hypothesis(&[1, 2, 3, 4], &[1], 1)]);
    }

    #[test]
    fn greedy_explanation_recomputes_votes_after_each_selection() {
        let hypotheses =
            build_composite_hypotheses(&[vec![1, 2], vec![1, 3], vec![2, 4], vec![2, 5]]);

        assert_eq!(hypotheses, vec![hypothesis(&[1, 2, 3, 4, 5], &[1, 2], 4)]);
    }
}
