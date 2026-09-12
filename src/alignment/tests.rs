use super::*;

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
