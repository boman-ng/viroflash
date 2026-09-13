use std::cell::RefCell;
use std::collections::VecDeque;
use std::fs::File;
use std::io::Write;
use std::sync::atomic::{AtomicU64, Ordering as AtomicOrdering};

use super::*;
use crate::alignment::hits_of;
use crate::gate::{GateEvaluation, GateScratch};
use crate::index::{build_index, IndexOptions, ReferenceRole};

static NEXT: AtomicU64 = AtomicU64::new(0);

fn retention_within_bound(stats: AlignmentRetention, threads: usize) -> bool {
    let maximum_pending_fragments = (threads + 1).saturating_mul(FASTQ_BATCH_RECORDS);
    stats.maximum_pending_fragments <= maximum_pending_fragments
        && stats.maximum_pending_sequence_bytes
            <= maximum_pending_fragments.saturating_mul(stats.maximum_fragment_bytes)
}

fn measured_alignment(
    index: &crate::index::ReferenceIndex,
    threads: usize,
    fragments: Vec<Fragment<'_>>,
) -> AlignmentRetention {
    let stats = RefCell::new(AlignmentRetention::default());
    let mut fragments = batches(&fragments);
    align_fragments_bounded_inner(
        &worker_config(index, threads),
        &mut || Ok(fragments.pop_front()),
        &mut |_| {},
        |observation| {
            *stats.borrow_mut() = observation;
        },
    )
    .unwrap();
    stats.into_inner()
}

fn batches(fragments: &[Fragment<'_>]) -> VecDeque<FragmentBatch> {
    fragments
        .chunks(FASTQ_BATCH_RECORDS)
        .map(FragmentBatch::from_fragments)
        .collect()
}

fn worker_config(index: &crate::index::ReferenceIndex, threads: usize) -> AnalysisWorkerConfig<'_> {
    AnalysisWorkerConfig {
        index_path: &index.mmi_path,
        contigs: &index.contigs,
        threads,
        bloom: None,
        minimum_hits: 2,
        minimum_covered_bases: 25,
    }
}

fn fixture() -> (std::path::PathBuf, crate::index::ReferenceIndex) {
    let root = std::env::temp_dir().join(format!(
        "viroflash-retention-{}-{}",
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

#[test]
fn workers_validate_pairs_before_bloom() {
    let (root, index) = fixture();
    let r1 = root.join("r1.fq");
    let r2 = root.join("r2.fq");
    let valid = "@pair\nACGT\n+\nIIII\n".repeat(FASTQ_BATCH_RECORDS * 2);
    std::fs::write(&r1, format!("{valid}@left\nACGT\n+\nIIII\n")).unwrap();
    std::fs::write(&r2, format!("{valid}@right\nACGT\n+\nIIII\n")).unwrap();
    let mut config = worker_config(&index, 4);
    config.bloom = Some(&index.bloom);
    let mut reader = crate::fastq::FragmentReader::open(&r1, Some(&r2)).unwrap();
    let error = align_fragments_bounded(
        config,
        || reader.next_batch(),
        |_| panic!("unselected alignment"),
    )
    .unwrap_err();
    assert!(error.contains("Paired IDs do not match"));
    std::fs::remove_dir_all(root).unwrap();
}

#[test]
fn completed_worker_is_refilled_before_batch_drains() {
    let (root, index) = fixture();
    let fragments = (0..(FASTQ_BATCH_RECORDS * 3) as u64)
        .map(|ordinal| Fragment {
            ordinal,
            id: "fragment",
            r1: &[b'G'; 20],
            r2: Some(&[b'T'; 20]),
        })
        .collect::<Vec<_>>();
    let mut fragments = batches(&fragments);
    let trace = RefCell::new(Vec::new());
    align_fragments_bounded_inner(
        &worker_config(&index, 2),
        &mut || Ok(fragments.pop_front()),
        &mut |_| {},
        |retention| trace.borrow_mut().push(retention.pending_fragments),
    )
    .unwrap();

    assert!(trace.into_inner().windows(3).any(|window| {
        window
            == [
                FASTQ_BATCH_RECORDS * 2,
                FASTQ_BATCH_RECORDS * 3,
                FASTQ_BATCH_RECORDS * 2,
            ]
    }));
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn alignment_retention_is_bounded_batches_not_total_selected() {
    let (root, index) = fixture();
    let threads = 4;
    let fragments = |count: usize| {
        (0..count as u64)
            .map(|ordinal| Fragment {
                ordinal,
                id: "fragment",
                r1: &[b'G'; 20],
                r2: Some(&[b'T'; 20]),
            })
            .collect::<Vec<_>>()
    };
    let small = measured_alignment(&index, threads, fragments(threads * FASTQ_BATCH_RECORDS));
    let large = measured_alignment(
        &index,
        threads,
        fragments((threads + 3) * FASTQ_BATCH_RECORDS),
    );
    assert!(retention_within_bound(small, threads));
    assert!(retention_within_bound(large, threads));
    assert!(large.total_sequence_bytes > small.total_sequence_bytes);

    let counterfactual = AlignmentRetention {
        maximum_pending_sequence_bytes: large.total_sequence_bytes,
        ..large
    };
    assert!(!retention_within_bound(counterfactual, threads));
    let _ = std::fs::remove_dir_all(root);
}
