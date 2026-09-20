use super::*;

const MINIMUM_HITS: usize = 2;
const MINIMUM_COVERED_BASES: usize = 25;

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|&symbol| {
            let (_, complement) = symbol_codes(symbol).unwrap();
            b"ACGTMRWSYKVHDBN"[complement as usize]
        })
        .collect()
}

fn encode_oracle(sequence: &[u8]) -> Option<u128> {
    let mut forward = 0_u128;
    let mut reverse = 0_u128;
    for &symbol in sequence {
        forward = (forward << 4) | u128::from(symbol_codes(symbol)?.0);
    }
    for &symbol in sequence.iter().rev() {
        reverse = (reverse << 4) | u128::from(symbol_codes(symbol)?.1);
    }
    Some(forward.min(reverse))
}

fn rolling(sequence: &[u8], k: usize) -> Vec<(usize, u128)> {
    let mut encoded = Vec::new();
    for_each_kmer(sequence, k, |start, code| encoded.push((start, code)));
    encoded
}

fn sdust_intervals(sequence: &[u8]) -> Vec<(usize, usize)> {
    let mut scratch = GateScratch::default();
    sdust_intervals_into(sequence, &mut scratch);
    scratch.intervals
}

#[test]
fn block_read_preserves_bloom_and_rejects_corruption_before_use() {
    let path = std::env::temp_dir().join(format!("vf-bloom-read-{}.bin", std::process::id()));
    let bloom = TargetKmerBloom {
        k: 21,
        words: (0..16384).map(|word| word as u64 * 0x10001).collect(),
        inserted_kmers: 1234,
    };
    bloom.write(&path).unwrap();
    let bytes = std::fs::read(&path).unwrap();
    let digest = format!("{:x}", Sha256::digest(&bytes));
    let loaded = TargetKmerBloom::read(&path, 21, &digest).unwrap();
    assert_eq!(loaded.words, bloom.words);
    assert_eq!(loaded.summary(), bloom.summary());
    let mut corrupt = bytes.clone();
    corrupt[32] ^= 1;
    std::fs::write(&path, &corrupt).unwrap();
    assert!(TargetKmerBloom::read(&path, 21, &digest)
        .unwrap_err()
        .contains("digest mismatch"));
    corrupt = bytes.clone();
    corrupt[24..32].copy_from_slice(&(1_u64 << 62).to_le_bytes());
    std::fs::write(&path, &corrupt).unwrap();
    assert!(TargetKmerBloom::read(&path, 21, &digest).is_err());
    std::fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
    assert!(TargetKmerBloom::read(&path, 21, &digest).is_err());
    corrupt = bytes;
    corrupt.push(0);
    std::fs::write(&path, &corrupt).unwrap();
    assert!(TargetKmerBloom::read(&path, 21, &digest).is_err());
    std::fs::remove_file(path).unwrap();
}

#[test]
fn fragment_gate_requires_a_current_hit_and_chain_compatible_end() {
    let records = [FastaRecord {
        description: String::new(),
        id: "target".into(),
        sequence: b"AGTCGATCCTAGGCTAACGTACGTA".to_vec(),
    }];
    let bloom = TargetKmerBloom::build(&records, 21);
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            b"AGTCGATCCTAGGCTAACGTA",
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Negative
    );
    assert_eq!(
        bloom.evaluate_fragment(
            b"XXXX",
            Some(b"TACGTACGTTAGCCTAGGATCGACT"),
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
    assert_eq!(
        bloom.evaluate_fragment(
            b"XXXX",
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::NotEvaluable
    );
}

#[test]
fn chain_compatible_gate_uses_mapper_count_and_coverage_requirements() {
    let single = b"AGTCGATCCTAGGCTAACGTA";
    let adjacent = b"AGTCGATCCTAGGCTAACGTAC";
    let covered = b"AGTCGATCCTAGGCTAACGTACGTA";
    let bloom = TargetKmerBloom::build(
        &[FastaRecord {
            description: String::new(),
            id: "target".into(),
            sequence: covered.to_vec(),
        }],
        21,
    );
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            single,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Negative
    );
    assert_eq!(
        bloom.evaluate_fragment(
            adjacent,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Negative
    );
    assert_eq!(
        bloom.evaluate_fragment(
            covered,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
}

#[test]
fn chain_condition_can_be_met_by_the_opposite_low_complexity_end() {
    let left = b"AGTCGATCCTAGGCTAACGTA";
    let right = [b'A'; 25];
    let bloom = TargetKmerBloom::build(
        &[
            FastaRecord {
                description: String::new(),
                id: "left".into(),
                sequence: left.to_vec(),
            },
            FastaRecord {
                description: String::new(),
                id: "right".into(),
                sequence: right.to_vec(),
            },
        ],
        21,
    );
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            left,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Negative
    );
    assert_eq!(
        bloom.evaluate_fragment(
            &right,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::NotEvaluable
    );
    assert_eq!(
        bloom.evaluate_fragment(
            left,
            Some(&right),
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
}

#[test]
fn exact_iupac_target_kmer_opens_gate() {
    let records = [FastaRecord {
        description: String::new(),
        id: "target".into(),
        sequence: b"ACGTMRWSYKVHDBNACGTMRWSYK".to_vec(),
    }];
    let bloom = TargetKmerBloom::build(&records, 21);
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            b"acgtmrwsykvhdbnacgtmrwsyk",
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
}

#[test]
fn short_invalid_and_fully_masked_fragments_are_not_evaluable() {
    let bloom = TargetKmerBloom::build(
        &[FastaRecord {
            description: String::new(),
            id: "target".into(),
            sequence: vec![b'A'; 100],
        }],
        21,
    );
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            b"ACGTACGTACGTACGTACGT",
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::NotEvaluable
    );
    assert_eq!(
        bloom.evaluate_fragment(
            b"XXXXXXXXXXXXXXXXXXXXX",
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::NotEvaluable
    );
    assert_eq!(
        bloom.evaluate_fragment(
            &[b'A'; 100],
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::NotEvaluable
    );
}

#[test]
fn rolling_iupac_encoder_matches_exact_slice_oracle() {
    let sequence = b"acgtmrwsykvhdbnXNVHDBKYWSRMACGT";
    for k in [1, 7, 21, 31] {
        let expected = sequence
            .windows(k)
            .enumerate()
            .filter_map(|(start, window)| encode_oracle(window).map(|code| (start, code)))
            .collect::<Vec<_>>();
        assert_eq!(rolling(sequence, k), expected, "k={k}");
    }
}

#[test]
fn rolling_iupac_encoder_is_reverse_complement_canonical() {
    let sequence = b"ACGTMRWSYKVHDBNACGTTGCATMRWSYK";
    let reverse = reverse_complement(sequence);
    let k = 21;
    let encoded = rolling(sequence, k);
    let reverse_encoded = rolling(&reverse, k);
    for (start, code) in encoded {
        let reverse_start = sequence.len() - k - start;
        assert_eq!(
            reverse_encoded
                .iter()
                .find(|(candidate, _)| *candidate == reverse_start)
                .map(|(_, candidate)| *candidate),
            Some(code),
            "start={start}"
        );
    }
}

#[test]
fn bloom_has_no_false_negatives_for_encoded_target_kmers() {
    let records = [FastaRecord {
        description: String::new(),
        id: "target".into(),
        sequence: b"AGTCGATCCTAGGCTAACGTATGCAGTACCGATGCTAGCATCGATCGTACGAT".to_vec(),
    }];
    let bloom = TargetKmerBloom::build(&records, 21);
    for (_, code) in rolling(&records[0].sequence, 21) {
        assert!(bloom.contains(code));
    }
    let mut scratch = GateScratch::default();
    assert_eq!(
        bloom.evaluate_fragment(
            &records[0].sequence,
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
    assert_eq!(
        bloom.evaluate_fragment(
            &reverse_complement(&records[0].sequence),
            None,
            MINIMUM_HITS,
            MINIMUM_COVERED_BASES,
            &mut scratch,
        ),
        GateEvaluation::Pass
    );
}

#[test]
fn sdust_matches_minimap2_reference_intervals() {
    assert_eq!(sdust_intervals(&[b'A'; 200]), vec![(0, 200)]);
    let split = [vec![b'A'; 80], vec![b'N'], vec![b'A'; 80]].concat();
    assert_eq!(sdust_intervals(&split), vec![(0, 80), (81, 161)]);
}

#[test]
fn sdust_is_symmetric_under_reverse_complement() {
    let sequence = b"AGTCGATCCTAGGCTAACGTA".repeat(4);
    let reverse = reverse_complement(&sequence);
    let mut left = vec![false; sequence.len()];
    let mut right = vec![false; sequence.len()];
    for (start, finish) in sdust_intervals(&sequence) {
        left[start..finish].fill(true);
    }
    for (start, finish) in sdust_intervals(&reverse) {
        right[start..finish].fill(true);
    }
    right.reverse();
    assert_eq!(left, right);
}

#[test]
fn early_exit_matches_exhaustive_fragment_predicates() {
    use std::collections::BTreeSet;
    let mut state = 91_u64;
    let reference: Vec<_> = (0..600)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect();
    let bloom = TargetKmerBloom::build(
        &[FastaRecord {
            id: "target".into(),
            description: String::new(),
            sequence: reference.clone(),
        }],
        21,
    );
    let exhaustive = |sequence: &[u8]| {
        let masks = sdust_intervals(sequence);
        let kmers = rolling(sequence, 21);
        let unmasked: Vec<_> = kmers
            .iter()
            .filter(|(start, _)| {
                !masks
                    .iter()
                    .any(|(left, right)| left <= start && start < right)
            })
            .collect();
        let hits: Vec<_> = kmers
            .iter()
            .filter(|(_, kmer)| bloom.contains(*kmer))
            .collect();
        let covered: BTreeSet<_> = hits
            .iter()
            .flat_map(|(start, _)| *start..start + 21)
            .collect();
        (
            unmasked.iter().any(|(_, kmer)| bloom.contains(*kmer)),
            !unmasked.is_empty(),
            hits.len() >= MINIMUM_HITS && covered.len() >= MINIMUM_COVERED_BASES,
        )
    };
    let mut reads = vec![
        vec![],
        vec![b'A'; 120],
        vec![b'N'; 120],
        reference[..20].to_vec(),
    ];
    for start in (0..450).step_by(17) {
        let exact = reference[start..start + 120].to_vec();
        let mut divergent = exact.clone();
        for offset in (7..120).step_by(19) {
            divergent[offset] = b'N';
        }
        reads.extend([exact, divergent]);
    }
    let mut scratch = GateScratch::default();
    for (index, left) in reads.iter().enumerate() {
        for right in [
            None,
            Some(reads[(index + 1) % reads.len()].as_slice()),
            Some(reads[1].as_slice()),
            Some(reads[4].as_slice()),
        ] {
            let (lh, le, lc) = exhaustive(left);
            let (rh, re, rc) = right.map_or((false, false, false), exhaustive);
            let expected = if (lh || rh) && (lc || rc) {
                GateEvaluation::Pass
            } else if le || re {
                GateEvaluation::Negative
            } else {
                GateEvaluation::NotEvaluable
            };
            assert_eq!(
                bloom.evaluate_fragment(
                    left,
                    right,
                    MINIMUM_HITS,
                    MINIMUM_COVERED_BASES,
                    &mut scratch
                ),
                expected
            );
        }
    }
}
