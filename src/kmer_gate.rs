use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::reference_group::FastaRecord;

const HASH_COUNT: u32 = 4;
const BITS_PER_KMER: usize = 16;
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
            for_each_kmer(&record.sequence, k, |kmer| {
                bloom.insert(kmer);
                bloom.inserted_kmers += 1;
            });
        }
        bloom
    }

    pub fn evaluate_fragment(&self, r1: &[u8], r2: Option<&[u8]>) -> GateEvaluation {
        let left = self.read_gate(r1);
        let right = r2.map(|sequence| self.read_gate(sequence));
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

    fn read_gate(&self, sequence: &[u8]) -> GateEvaluation {
        let mut present = false;
        let mut evaluable = false;
        for_each_kmer(sequence, self.k, |kmer| {
            evaluable = true;
            present |= self.contains(kmer);
        });
        if present {
            GateEvaluation::Pass
        } else if evaluable {
            GateEvaluation::Negative
        } else {
            GateEvaluation::NotEvaluable
        }
    }

    fn insert(&mut self, kmer: u64) {
        for index in hash_indices(kmer, self.words.len() * 64) {
            self.words[index / 64] |= 1_u64 << (index % 64);
        }
    }

    fn contains(&self, kmer: u64) -> bool {
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

fn for_each_kmer(mut sequence: &[u8], k: usize, mut visit: impl FnMut(u64)) {
    if sequence.len() < k {
        return;
    }
    while sequence.len() >= k {
        if let Some(encoded) = encode_canonical(&sequence[..k]) {
            visit(encoded);
        }
        sequence = &sequence[1..];
    }
}

fn encode_canonical(sequence: &[u8]) -> Option<u64> {
    let forward = sequence
        .iter()
        .map(u8::to_ascii_uppercase)
        .collect::<Vec<_>>();
    let reverse = forward
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => Some(b'T'),
            b'T' => Some(b'A'),
            b'C' => Some(b'G'),
            b'G' => Some(b'C'),
            b'M' => Some(b'K'),
            b'K' => Some(b'M'),
            b'R' => Some(b'Y'),
            b'Y' => Some(b'R'),
            b'W' => Some(b'W'),
            b'S' => Some(b'S'),
            b'V' => Some(b'B'),
            b'B' => Some(b'V'),
            b'H' => Some(b'D'),
            b'D' => Some(b'H'),
            b'N' => Some(b'N'),
            _ => None,
        })
        .collect::<Option<Vec<_>>>()?;
    let canonical = if forward <= reverse { forward } else { reverse };
    let digest = blake3::hash(&canonical);
    Some(u64::from_le_bytes(
        digest.as_bytes()[..8].try_into().unwrap(),
    ))
}

fn hash_indices(value: u64, bit_count: usize) -> impl Iterator<Item = usize> {
    let first = mix64(value);
    let second = mix64(value ^ 0x9e37_79b9_7f4a_7c15) | 1;
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
mod tests {
    use super::*;

    #[test]
    fn any_target_kmer_opens_gate_for_whole_fragment() {
        let records = [FastaRecord {
            id: "target".into(),
            sequence: b"ACGTACGTACGTACGTACGTA".to_vec(),
        }];
        let bloom = TargetKmerBloom::build(&records, 21);
        assert_eq!(
            bloom.evaluate_fragment(b"ACGTACGTACGTACGTACGTA", None),
            GateEvaluation::Pass
        );
        assert_eq!(
            bloom.evaluate_fragment(b"NNNN", Some(b"TACGTACGTACGTACGTACGT")),
            GateEvaluation::Pass
        );
        assert_eq!(
            bloom.evaluate_fragment(b"NNNN", None),
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
        assert_eq!(
            bloom.evaluate_fragment(b"acgtmrwsykvhdbnacgtmr", None),
            GateEvaluation::Pass
        );
    }

    #[test]
    fn fragments_without_an_encodable_kmer_are_not_evaluable() {
        let bloom = TargetKmerBloom::build(
            &[FastaRecord {
                id: "target".into(),
                sequence: b"ACGTACGTACGTACGTACGTACGT".to_vec(),
            }],
            21,
        );
        assert_eq!(
            bloom.evaluate_fragment(b"ACGTACGTACGTACGTACGT", None),
            GateEvaluation::NotEvaluable
        );
        assert_eq!(
            bloom.evaluate_fragment(b"XXXXXXXXXXXXXXXXXXXXX", None),
            GateEvaluation::NotEvaluable
        );
    }
}
