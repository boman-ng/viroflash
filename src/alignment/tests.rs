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

mod chain_gate {
    use super::*;
    use crate::gate::{GateEvaluation, GateScratch};
    use crate::index::{build_index, IndexOptions};
    use std::fs::File;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};
    static NEXT: AtomicU64 = AtomicU64::new(0);
    fn pseudo_random_sequence(length: usize, mut state: u64) -> Vec<u8> {
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state as usize) & 3]
            })
            .collect()
    }

    #[test]
    fn emitted_target_mappings_satisfy_the_short_read_chain_gate() {
        let root = std::env::temp_dir().join(format!(
            "viroflash-chain-gate-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, AtomicOrdering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let host = pseudo_random_sequence(8_192, 11);
        let target = pseudo_random_sequence(8_192, 29);
        for (name, sequence) in [("host", &host), ("target", &target)] {
            let mut file = File::create(root.join(format!("{name}.fa"))).unwrap();
            writeln!(file, ">{name}\n{}", String::from_utf8_lossy(sequence)).unwrap();
        }
        let index = build_index(&IndexOptions {
            host_fa: root.join("host.fa"),
            target_fa: root.join("target.fa"),
            out_dir: root.join("index"),
            threads: 1,
        })
        .unwrap();
        let aligner = CompetitiveAligner::open(&index.mmi_path).unwrap();
        let (minimum_hits, minimum_covered_bases) =
            CompetitiveAligner::short_read_chain_requirements(21).unwrap();
        assert_eq!((minimum_hits, minimum_covered_bases), (2, 25));

        let mut mapped_reads = 0;
        let mut gate_scratch = GateScratch::default();
        for (number, start) in (0..target.len() - 150).step_by(257).enumerate() {
            let read = target[start..start + 150].to_vec();
            let mappings = aligner
                .aligner
                .map(
                    &read,
                    false,
                    false,
                    None,
                    None,
                    Some(format!("target-read-{number}").as_bytes()),
                )
                .unwrap();
            if hits_of(&mappings, &index.contigs)
                .iter()
                .any(|hit| hit.role == ReferenceRole::Target)
            {
                mapped_reads += 1;
                assert_eq!(
                    index.bloom.evaluate_fragment(
                        &read,
                        None,
                        minimum_hits,
                        minimum_covered_bases,
                        &mut gate_scratch,
                    ),
                    GateEvaluation::Pass
                );
            }
        }
        assert!(mapped_reads > 20);
        let _ = std::fs::remove_dir_all(root);
    }
}
