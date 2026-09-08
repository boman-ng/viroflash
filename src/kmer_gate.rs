use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::analysis_profile::{SDUST_THRESHOLD, SDUST_WINDOW};
use crate::reference_group::FastaRecord;

const HASH_COUNT: u32 = 4;
const BITS_PER_KMER: usize = 16;
const MAX_KMER_SYMBOLS: usize = 31;
const MAGIC: &[u8; 8] = b"VFBLOOM1";

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct BloomSummary {
    pub inserted_kmers: u64,
    pub bit_count: u64,
    pub hash_count: u32,
    pub fill_fraction: f64,
    pub theoretical_false_positive_rate: f64,
}

#[derive(Debug, Clone)]
pub struct TargetKmerBloom {
    k: usize,
    words: Vec<u64>,
    inserted_kmers: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GateEvaluation {
    Pass,
    Negative,
    NotEvaluable,
}

#[derive(Default)]
pub struct GateScratch {
    intervals: Vec<(usize, usize)>,
    perfect_intervals: Vec<PerfectInterval>,
}

impl TargetKmerBloom {
    pub fn build(records: &[FastaRecord], k: usize) -> Self {
        let estimated = records
            .iter()
            .map(|record| record.sequence.len().saturating_sub(k - 1))
            .sum::<usize>()
            .max(1);
        let bit_count = (estimated * BITS_PER_KMER).next_power_of_two().max(64);
        let mut bloom = Self {
            k,
            words: vec![0; bit_count / 64],
            inserted_kmers: 0,
        };
        for record in records {
            for_each_kmer(&record.sequence, k, |_, kmer| {
                bloom.insert(kmer);
                bloom.inserted_kmers += 1;
            });
        }
        bloom
    }

    pub fn evaluate_fragment(
        &self,
        r1: &[u8],
        r2: Option<&[u8]>,
        scratch: &mut GateScratch,
    ) -> GateEvaluation {
        let left = self.read_gate(r1, scratch);
        let right = r2.map(|sequence| self.read_gate(sequence, scratch));
        if left == GateEvaluation::Pass || right == Some(GateEvaluation::Pass) {
            GateEvaluation::Pass
        } else if left == GateEvaluation::NotEvaluable
            && right.is_none_or(|evaluation| evaluation == GateEvaluation::NotEvaluable)
        {
            GateEvaluation::NotEvaluable
        } else {
            GateEvaluation::Negative
        }
    }

    fn read_gate(&self, sequence: &[u8], scratch: &mut GateScratch) -> GateEvaluation {
        sdust_intervals_into(sequence, scratch);
        let mut present = false;
        let mut evaluable = false;
        let mut interval_index = 0;
        for_each_kmer(sequence, self.k, |start, kmer| {
            while interval_index < scratch.intervals.len()
                && scratch.intervals[interval_index].1 <= start
            {
                interval_index += 1;
            }
            if interval_index >= scratch.intervals.len()
                || scratch.intervals[interval_index].0 > start
            {
                evaluable = true;
                present |= self.contains(kmer);
            }
        });
        if present {
            GateEvaluation::Pass
        } else if evaluable {
            GateEvaluation::Negative
        } else {
            GateEvaluation::NotEvaluable
        }
    }

    fn insert(&mut self, kmer: u128) {
        for index in hash_indices(kmer, self.words.len() * 64) {
            self.words[index / 64] |= 1_u64 << (index % 64);
        }
    }

    fn contains(&self, kmer: u128) -> bool {
        hash_indices(kmer, self.words.len() * 64)
            .all(|index| self.words[index / 64] & (1_u64 << (index % 64)) != 0)
    }

    pub fn summary(&self) -> BloomSummary {
        let set_bits = self
            .words
            .iter()
            .map(|word| word.count_ones() as u64)
            .sum::<u64>();
        let bit_count = (self.words.len() * 64) as u64;
        let fill_fraction = set_bits as f64 / bit_count as f64;
        BloomSummary {
            inserted_kmers: self.inserted_kmers,
            bit_count,
            hash_count: HASH_COUNT,
            fill_fraction,
            theoretical_false_positive_rate: fill_fraction.powi(HASH_COUNT as i32),
        }
    }

    pub fn write(&self, path: &Path) -> Result<(), String> {
        let mut writer = BufWriter::new(
            File::create(path)
                .map_err(|error| format!("Cannot create {}: {error}", path.display()))?,
        );
        writer
            .write_all(MAGIC)
            .and_then(|_| writer.write_all(&(self.k as u64).to_le_bytes()))
            .and_then(|_| writer.write_all(&self.inserted_kmers.to_le_bytes()))
            .and_then(|_| writer.write_all(&(self.words.len() as u64).to_le_bytes()))
            .map_err(|error| format!("Failed to write {}: {error}", path.display()))?;
        for word in &self.words {
            writer
                .write_all(&word.to_le_bytes())
                .map_err(|error| format!("Failed to write {}: {error}", path.display()))?;
        }
        writer
            .flush()
            .map_err(|error| format!("Failed to write {}: {error}", path.display()))
    }

    pub fn read(path: &Path, expected_k: usize) -> Result<Self, String> {
        let mut reader = BufReader::new(
            File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?,
        );
        let mut magic = [0; 8];
        reader
            .read_exact(&mut magic)
            .map_err(|error| format!("Invalid Bloom file {}: {error}", path.display()))?;
        if &magic != MAGIC {
            return Err(format!(
                "{} is not a current Viroflash target Bloom",
                path.display()
            ));
        }
        let k = read_u64(&mut reader)? as usize;
        let inserted_kmers = read_u64(&mut reader)?;
        let word_count = usize::try_from(read_u64(&mut reader)?)
            .map_err(|_| "Bloom word count exceeds platform limits".to_string())?;
        if k != expected_k || word_count == 0 || !word_count.is_power_of_two() {
            return Err("Bloom metadata does not match the frozen profile".into());
        }
        let mut words = vec![0; word_count];
        for word in &mut words {
            *word = read_u64(&mut reader)?;
        }
        let mut trailing = [0];
        if reader
            .read(&mut trailing)
            .map_err(|error| error.to_string())?
            != 0
        {
            return Err("Bloom file has trailing bytes".into());
        }
        Ok(Self {
            k,
            words,
            inserted_kmers,
        })
    }
}

fn read_u64(reader: &mut impl Read) -> Result<u64, String> {
    let mut bytes = [0; 8];
    reader
        .read_exact(&mut bytes)
        .map_err(|error| format!("Invalid Bloom file: {error}"))?;
    Ok(u64::from_le_bytes(bytes))
}

fn for_each_kmer(sequence: &[u8], k: usize, mut visit: impl FnMut(usize, u128)) {
    if k == 0 || k > MAX_KMER_SYMBOLS || sequence.len() < k {
        return;
    }
    let mask = (1_u128 << (4 * k)) - 1;
    let reverse_shift = 4 * (k - 1);
    let mut forward = 0_u128;
    let mut reverse = 0_u128;
    let mut valid = 0;
    for (index, &symbol) in sequence.iter().enumerate() {
        if let Some((code, complement)) = symbol_codes(symbol) {
            forward = ((forward << 4) | u128::from(code)) & mask;
            reverse = (reverse >> 4) | (u128::from(complement) << reverse_shift);
            valid = (valid + 1).min(k);
            if valid == k {
                visit(index + 1 - k, forward.min(reverse));
            }
        } else {
            forward = 0;
            reverse = 0;
            valid = 0;
        }
    }
}

fn symbol_codes(symbol: u8) -> Option<(u8, u8)> {
    match symbol.to_ascii_uppercase() {
        b'A' => Some((0, 3)),
        b'C' => Some((1, 2)),
        b'G' => Some((2, 1)),
        b'T' => Some((3, 0)),
        b'M' => Some((4, 9)),
        b'R' => Some((5, 8)),
        b'W' => Some((6, 6)),
        b'S' => Some((7, 7)),
        b'Y' => Some((8, 5)),
        b'K' => Some((9, 4)),
        b'V' => Some((10, 13)),
        b'H' => Some((11, 12)),
        b'D' => Some((12, 11)),
        b'B' => Some((13, 10)),
        b'N' => Some((14, 14)),
        _ => None,
    }
}

fn hash_indices(value: u128, bit_count: usize) -> impl Iterator<Item = usize> {
    let low = value as u64;
    let high = (value >> 64) as u64;
    let folded = low ^ mix64(high);
    let first = mix64(folded);
    let second = mix64(folded ^ 0x9e37_79b9_7f4a_7c15) | 1;
    (0..HASH_COUNT).map(move |index| {
        first.wrapping_add(u64::from(index).wrapping_mul(second)) as usize & (bit_count - 1)
    })
}

#[derive(Debug)]
struct PerfectInterval {
    start: usize,
    finish: usize,
    score: i64,
    length: i64,
}

fn save_masked_regions(
    result: &mut Vec<(usize, usize)>,
    perfect: &mut Vec<PerfectInterval>,
    start: usize,
) {
    if perfect.is_empty() || perfect[perfect.len() - 1].start >= start {
        return;
    }
    let last = &perfect[perfect.len() - 1];
    if let Some(previous) = result.last_mut() {
        if last.start <= previous.1 {
            previous.1 = previous.1.max(last.finish);
        } else {
            result.push((last.start, last.finish));
        }
    } else {
        result.push((last.start, last.finish));
    }
    let mut keep = 0;
    for index in (0..perfect.len()).rev() {
        if perfect[index].start >= start {
            keep = index + 1;
            break;
        }
    }
    perfect.truncate(keep);
}

const SDUST_RING_CAPACITY: usize = 128;
const SDUST_RING_MASK: usize = SDUST_RING_CAPACITY - 1;

struct SdustRing {
    values: [u32; SDUST_RING_CAPACITY],
    head: usize,
    length: usize,
}

impl SdustRing {
    fn new() -> Self {
        Self {
            values: [0; SDUST_RING_CAPACITY],
            head: 0,
            length: 0,
        }
    }

    fn push(&mut self, value: u32) {
        self.values[(self.head + self.length) & SDUST_RING_MASK] = value;
        self.length += 1;
    }

    fn pop_front(&mut self) -> u32 {
        debug_assert!(self.length > 0);
        let value = self.values[self.head];
        self.head = (self.head + 1) & SDUST_RING_MASK;
        self.length -= 1;
        value
    }

    fn at(&self, index: usize) -> u32 {
        self.values[(self.head + index) & SDUST_RING_MASK]
    }
}

#[allow(clippy::too_many_arguments)]
fn shift_sdust_window(
    queue: &mut SdustRing,
    triplet: u32,
    active_length: &mut usize,
    window_score: &mut i64,
    suffix_score: &mut i64,
    window_counts: &mut [i64; 64],
    suffix_counts: &mut [i64; 64],
) {
    if queue.length >= SDUST_WINDOW - 2 {
        let symbol = (queue.pop_front() as usize) & 63;
        window_counts[symbol] -= 1;
        *window_score -= window_counts[symbol];
        if *active_length > queue.length {
            *active_length -= 1;
            suffix_counts[symbol] -= 1;
            *suffix_score -= suffix_counts[symbol];
        }
    }
    queue.push(triplet);
    *active_length += 1;
    let symbol = triplet as usize;
    *window_score += window_counts[symbol];
    window_counts[symbol] += 1;
    *suffix_score += suffix_counts[symbol];
    suffix_counts[symbol] += 1;
    if suffix_counts[symbol] * 10 > SDUST_THRESHOLD << 1 {
        loop {
            let removed = (queue.at(queue.length - *active_length) as usize) & 63;
            suffix_counts[removed] -= 1;
            *suffix_score -= suffix_counts[removed];
            *active_length -= 1;
            if removed == symbol {
                break;
            }
        }
    }
}

fn find_perfect_intervals(
    perfect: &mut Vec<PerfectInterval>,
    queue: &SdustRing,
    start: usize,
    active_length: usize,
    suffix_score: i64,
    suffix_counts: &[i64; 64],
) {
    let mut counts = *suffix_counts;
    let mut score = suffix_score;
    let mut maximum_score = 0;
    let mut maximum_length = 0;
    for index in (0..queue.length.saturating_sub(active_length)).rev() {
        let triplet = (queue.at(index) as usize) & 63;
        score += counts[triplet];
        counts[triplet] += 1;
        let new_score = score;
        let new_length = (queue.length - index - 1) as i64;
        if new_score * 10 > SDUST_THRESHOLD * new_length {
            let mut insertion = 0;
            while insertion < perfect.len() && perfect[insertion].start >= index + start {
                if maximum_score == 0
                    || perfect[insertion].score * maximum_length
                        > maximum_score * perfect[insertion].length
                {
                    maximum_score = perfect[insertion].score;
                    maximum_length = perfect[insertion].length;
                }
                insertion += 1;
            }
            if maximum_score == 0 || new_score * maximum_length >= maximum_score * new_length {
                maximum_score = new_score;
                maximum_length = new_length;
                perfect.insert(
                    insertion,
                    PerfectInterval {
                        start: index + start,
                        finish: queue.length + 2 + start,
                        score: new_score,
                        length: new_length,
                    },
                );
            }
        }
    }
}

fn sdust_base(symbol: u8) -> Option<u8> {
    match symbol.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn sdust_intervals_into(sequence: &[u8], scratch: &mut GateScratch) {
    scratch.intervals.clear();
    scratch.perfect_intervals.clear();
    let mut queue = SdustRing::new();
    let mut window_counts = [0_i64; 64];
    let mut suffix_counts = [0_i64; 64];
    let mut suffix_score = 0;
    let mut window_score = 0;
    let mut active_length = 0;
    let mut contiguous_length = 0_usize;
    let mut triplet = 0_u32;

    for index in 0..=sequence.len() {
        if let Some(base) = sequence.get(index).copied().and_then(sdust_base) {
            contiguous_length += 1;
            triplet = ((triplet << 2) | u32::from(base)) & 0x3f;
            if contiguous_length >= 3 {
                let start = contiguous_length.saturating_sub(SDUST_WINDOW)
                    + (index + 1 - contiguous_length);
                save_masked_regions(
                    &mut scratch.intervals,
                    &mut scratch.perfect_intervals,
                    start,
                );
                shift_sdust_window(
                    &mut queue,
                    triplet,
                    &mut active_length,
                    &mut window_score,
                    &mut suffix_score,
                    &mut window_counts,
                    &mut suffix_counts,
                );
                if window_score * 10 > active_length as i64 * SDUST_THRESHOLD {
                    find_perfect_intervals(
                        &mut scratch.perfect_intervals,
                        &queue,
                        start,
                        active_length,
                        suffix_score,
                        &suffix_counts,
                    );
                }
            }
        } else {
            let mut start = contiguous_length.saturating_sub(SDUST_WINDOW.saturating_sub(1))
                + (index + 1 - contiguous_length);
            while !scratch.perfect_intervals.is_empty() {
                save_masked_regions(
                    &mut scratch.intervals,
                    &mut scratch.perfect_intervals,
                    start,
                );
                start += 1;
            }
            contiguous_length = 0;
            triplet = 0;
        }
    }
}

fn mix64(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn any_target_kmer_opens_gate_for_whole_fragment() {
        let records = [FastaRecord {
            id: "target".into(),
            sequence: b"AGTCGATCCTAGGCTAACGTA".to_vec(),
        }];
        let bloom = TargetKmerBloom::build(&records, 21);
        let mut scratch = GateScratch::default();
        assert_eq!(
            bloom.evaluate_fragment(b"AGTCGATCCTAGGCTAACGTA", None, &mut scratch),
            GateEvaluation::Pass
        );
        assert_eq!(
            bloom.evaluate_fragment(b"XXXX", Some(b"TACGTTAGCCTAGGATCGACT"), &mut scratch,),
            GateEvaluation::Pass
        );
        assert_eq!(
            bloom.evaluate_fragment(b"XXXX", None, &mut scratch),
            GateEvaluation::NotEvaluable
        );
    }

    #[test]
    fn exact_iupac_target_kmer_opens_gate() {
        let records = [FastaRecord {
            id: "target".into(),
            sequence: b"ACGTMRWSYKVHDBNACGTMR".to_vec(),
        }];
        let bloom = TargetKmerBloom::build(&records, 21);
        let mut scratch = GateScratch::default();
        assert_eq!(
            bloom.evaluate_fragment(b"acgtmrwsykvhdbnacgtmr", None, &mut scratch),
            GateEvaluation::Pass
        );
    }

    #[test]
    fn short_invalid_and_fully_masked_fragments_are_not_evaluable() {
        let bloom = TargetKmerBloom::build(
            &[FastaRecord {
                id: "target".into(),
                sequence: vec![b'A'; 100],
            }],
            21,
        );
        let mut scratch = GateScratch::default();
        assert_eq!(
            bloom.evaluate_fragment(b"ACGTACGTACGTACGTACGT", None, &mut scratch),
            GateEvaluation::NotEvaluable
        );
        assert_eq!(
            bloom.evaluate_fragment(b"XXXXXXXXXXXXXXXXXXXXX", None, &mut scratch),
            GateEvaluation::NotEvaluable
        );
        assert_eq!(
            bloom.evaluate_fragment(&[b'A'; 100], None, &mut scratch),
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
            id: "target".into(),
            sequence: b"AGTCGATCCTAGGCTAACGTATGCAGTACCGATGCTAGCATCGATCGTACGAT".to_vec(),
        }];
        let bloom = TargetKmerBloom::build(&records, 21);
        for (_, code) in rolling(&records[0].sequence, 21) {
            assert!(bloom.contains(code));
        }
        let mut scratch = GateScratch::default();
        assert_eq!(
            bloom.evaluate_fragment(&records[0].sequence, None, &mut scratch),
            GateEvaluation::Pass
        );
        assert_eq!(
            bloom.evaluate_fragment(
                &reverse_complement(&records[0].sequence),
                None,
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
}
