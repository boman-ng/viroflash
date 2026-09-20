use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

mod sdust;
use crate::index::reference::FastaRecord;
use sdust::{sdust_intervals_into, PerfectInterval};

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

#[derive(Debug, Clone, Copy)]
struct ReadGateEvaluation {
    current: GateEvaluation,
    chain_compatible: bool,
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
        minimum_hits: usize,
        minimum_covered_bases: usize,
        scratch: &mut GateScratch,
    ) -> GateEvaluation {
        let left = self.read_gate(r1, minimum_hits, minimum_covered_bases, scratch);
        if left.current == GateEvaluation::Pass && left.chain_compatible {
            return GateEvaluation::Pass;
        }
        let right = r2
            .map(|sequence| self.read_gate(sequence, minimum_hits, minimum_covered_bases, scratch));
        let current = if left.current == GateEvaluation::Pass
            || right.is_some_and(|evaluation| evaluation.current == GateEvaluation::Pass)
        {
            GateEvaluation::Pass
        } else if left.current == GateEvaluation::NotEvaluable
            && right.is_none_or(|evaluation| evaluation.current == GateEvaluation::NotEvaluable)
        {
            GateEvaluation::NotEvaluable
        } else {
            GateEvaluation::Negative
        };
        if current == GateEvaluation::Pass
            && (left.chain_compatible
                || right.is_some_and(|evaluation| evaluation.chain_compatible))
        {
            GateEvaluation::Pass
        } else if current == GateEvaluation::Pass {
            GateEvaluation::Negative
        } else {
            current
        }
    }

    fn read_gate(
        &self,
        sequence: &[u8],
        minimum_hits: usize,
        minimum_covered_bases: usize,
        scratch: &mut GateScratch,
    ) -> ReadGateEvaluation {
        sdust_intervals_into(sequence, scratch);
        let mut present = false;
        let mut evaluable = false;
        let mut interval_index = 0;
        let mut target_hits = 0;
        let mut covered_bases = 0;
        let mut hit_interval_end = 0;
        for_each_kmer_while(sequence, self.k, |start, kmer| {
            let target_present = self.contains(kmer);
            if target_present {
                target_hits += 1;
                let end = start + self.k;
                if target_hits == 1 || start >= hit_interval_end {
                    covered_bases += self.k;
                } else if end > hit_interval_end {
                    covered_bases += end - hit_interval_end;
                }
                hit_interval_end = hit_interval_end.max(end);
            }
            while interval_index < scratch.intervals.len()
                && scratch.intervals[interval_index].1 <= start
            {
                interval_index += 1;
            }
            if interval_index >= scratch.intervals.len()
                || scratch.intervals[interval_index].0 > start
            {
                evaluable = true;
                present |= target_present;
            }
            !(present && target_hits >= minimum_hits && covered_bases >= minimum_covered_bases)
        });
        let current = if present {
            GateEvaluation::Pass
        } else if evaluable {
            GateEvaluation::Negative
        } else {
            GateEvaluation::NotEvaluable
        };
        ReadGateEvaluation {
            current,
            chain_compatible: target_hits >= minimum_hits && covered_bases >= minimum_covered_bases,
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

    pub fn read(path: &Path, expected_k: usize, expected_digest: &str) -> Result<Self, String> {
        let mut reader =
            File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?;
        let mut header = [0; 32];
        reader
            .read_exact(&mut header)
            .map_err(|error| format!("Invalid Bloom file {}: {error}", path.display()))?;
        if &header[..8] != MAGIC {
            return Err(format!(
                "{} is not a current Viroflash target Bloom",
                path.display()
            ));
        }
        let k = u64::from_le_bytes([
            header[8], header[9], header[10], header[11], header[12], header[13], header[14],
            header[15],
        ]) as usize;
        let inserted_kmers = u64::from_le_bytes([
            header[16], header[17], header[18], header[19], header[20], header[21], header[22],
            header[23],
        ]);
        let word_count = usize::try_from(u64::from_le_bytes([
            header[24], header[25], header[26], header[27], header[28], header[29], header[30],
            header[31],
        ]))
        .map_err(|_| "Bloom word count exceeds platform limits".to_string())?;
        if k != expected_k || word_count == 0 || !word_count.is_power_of_two() {
            return Err("Bloom metadata does not match the frozen profile".into());
        }
        let expected_bytes = u64::try_from(word_count)
            .ok()
            .and_then(|words| words.checked_mul(8))
            .and_then(|bytes| bytes.checked_add(32))
            .ok_or_else(|| "Bloom byte count exceeds platform limits".to_string())?;
        let actual_bytes = reader
            .metadata()
            .map_err(|error| format!("Cannot inspect {}: {error}", path.display()))?
            .len();
        if actual_bytes != expected_bytes {
            return Err("Bloom file length does not match its word count".into());
        }
        let mut hasher = Sha256::new();
        hasher.update(header);
        let mut words = vec![0; word_count];
        let mut buffer = [0; 64 * 1024];
        for chunk in words.chunks_mut(buffer.len() / 8) {
            let bytes = &mut buffer[..chunk.len() * 8];
            reader
                .read_exact(bytes)
                .map_err(|error| format!("Invalid Bloom file: {error}"))?;
            hasher.update(&*bytes);
            for (word, bytes) in chunk.iter_mut().zip(bytes.chunks_exact(8)) {
                *word = u64::from_le_bytes([
                    bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5], bytes[6], bytes[7],
                ]);
            }
        }
        let mut trailing = [0];
        if reader
            .read(&mut trailing)
            .map_err(|error| error.to_string())?
            != 0
        {
            return Err("Bloom file has trailing bytes".into());
        }
        if format!("{:x}", hasher.finalize()) != expected_digest {
            return Err(format!(
                "Index artifact digest mismatch: {}",
                path.display()
            ));
        }
        Ok(Self {
            k,
            words,
            inserted_kmers,
        })
    }
}

fn for_each_kmer(sequence: &[u8], k: usize, mut visit: impl FnMut(usize, u128)) {
    for_each_kmer_while(sequence, k, |start, kmer| {
        visit(start, kmer);
        true
    });
}

fn for_each_kmer_while(sequence: &[u8], k: usize, mut visit: impl FnMut(usize, u128) -> bool) {
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
            if valid == k && !visit(index + 1 - k, forward.min(reverse)) {
                break;
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

fn mix64(mut value: u64) -> u64 {
    value ^= value >> 30;
    value = value.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value ^= value >> 27;
    value = value.wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

#[cfg(test)]
mod tests;
