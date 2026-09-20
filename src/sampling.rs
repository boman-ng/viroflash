use std::collections::{BinaryHeap, VecDeque};

use crate::fastq::{Fragment, FASTQ_BATCH_RECORDS};
use crate::profile::AnalysisProfile;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RunMode {
    Full,
    #[default]
    Screen,
}
impl RunMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::Screen => "screen",
        }
    }
    pub fn population(self) -> &'static str {
        match self {
            Self::Full => "bloom_candidates",
            Self::Screen => "input_fragments",
        }
    }
}
impl std::str::FromStr for RunMode {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, String> {
        match value {
            "full" => Ok(Self::Full),
            "screen" => Ok(Self::Screen),
            _ => Err("--mode must be full or screen".into()),
        }
    }
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum Precision {
    Fast,
    #[default]
    Standard,
    Sensitive,
}
impl Precision {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fast => "fast",
            Self::Standard => "standard",
            Self::Sensitive => "sensitive",
        }
    }
    pub fn denominator(self) -> u64 {
        match self {
            Self::Fast => 10_000,
            Self::Standard => 100_000,
            Self::Sensitive => 1_000_000,
        }
    }
    pub fn minimum_fraction(self) -> f64 {
        1.0 / self.denominator() as f64
    }
}
impl std::str::FromStr for Precision {
    type Err = String;
    fn from_str(value: &str) -> Result<Self, String> {
        match value {
            "fast" => Ok(Self::Fast),
            "standard" => Ok(Self::Standard),
            "sensitive" => Ok(Self::Sensitive),
            _ => Err("--precision must be fast, standard, or sensitive".into()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SamplingDesign {
    pub precision: Precision,
    pub population_fragments: u64,
    pub selection_probability: f64,
    pub sample_capacity: u64,
}

/// Invert the binomial zero-hit probability, with a per-group miss budget.
pub(crate) fn sample_capacity(
    family_size: usize,
    profile: AnalysisProfile,
    precision: Precision,
) -> Result<usize, String> {
    if family_size == 0 {
        return Err("Reference index contains no target groups".into());
    }
    let miss = (profile.familywise_miss_probability / family_size as f64).next_down();
    let delta = precision.minimum_fraction().next_down();
    let mut n = (miss.ln() / (-delta).ln_1p()).ceil() as u64;
    while miss_probability_upper_bound(delta, n) > miss {
        n += 1;
    }
    usize::try_from(n).map_err(|_| "Sampling capacity exceeds this platform's address space".into())
}

fn miss_probability_upper_bound(probability: f64, mut exponent: u64) -> f64 {
    let mut factor = (1.0 - probability).next_up();
    let mut result = 1.0;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = (result * factor).next_up();
        }
        exponent >>= 1;
        if exponent > 0 {
            factor = (factor * factor).next_up();
        }
    }
    result
}

#[derive(Clone)]
pub(crate) struct FragmentKeys(blake3::Hasher);
impl FragmentKeys {
    pub fn new(profile: &str, index: &str) -> Self {
        let mut prefix = blake3::Hasher::new();
        prefix.update(b"viroflash-fragment-selection-v1\0");
        prefix.update(profile.as_bytes());
        prefix.update(index.as_bytes());
        Self(prefix)
    }
    pub fn key(&self, id: &str, ordinal: u64) -> u128 {
        let mut hasher = self.0.clone();
        hasher.update(&(id.len() as u64).to_be_bytes());
        hasher.update(id.as_bytes());
        hasher.update(&ordinal.to_be_bytes());
        let mut bytes = [0; 16];
        bytes.copy_from_slice(&hasher.finalize().as_bytes()[..16]);
        u128::from_be_bytes(bytes)
    }
}

struct StoredFragment {
    ordinal: u64,
    id_end: usize,
    r1_end: usize,
    data: Vec<u8>,
}
impl StoredFragment {
    fn new(f: Fragment<'_>) -> Self {
        let mut stored = Self {
            ordinal: 0,
            id_end: 0,
            r1_end: 0,
            data: Vec::with_capacity(f.id.len() + f.r1.len() + f.r2.map_or(0, <[u8]>::len)),
        };
        stored.replace(f);
        stored
    }
    fn replace(&mut self, f: Fragment<'_>) {
        self.ordinal = f.ordinal;
        self.id_end = f.id.len();
        self.r1_end = self.id_end + f.r1.len();
        self.data.clear();
        self.data.extend_from_slice(f.id.as_bytes());
        self.data.extend_from_slice(f.r1);
        if let Some(read) = f.r2 {
            self.data.extend_from_slice(read);
        }
    }
    fn fragment(&self, paired: bool) -> Fragment<'_> {
        Fragment {
            ordinal: self.ordinal,
            id: std::str::from_utf8(&self.data[..self.id_end]).expect("canonical FASTQ ID"),
            r1: &self.data[self.id_end..self.r1_end],
            r2: paired.then_some(&self.data[self.r1_end..]),
        }
    }
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
struct SampleKey {
    key: u128,
    ordinal: u64,
    slot: usize,
}

pub(crate) struct BottomKSampler {
    capacity: usize,
    heap: BinaryHeap<SampleKey>,
    fragments: Vec<StoredFragment>,
    paired: bool,
}
impl BottomKSampler {
    pub fn new(capacity: usize, paired: bool) -> Self {
        Self {
            capacity,
            heap: BinaryHeap::with_capacity(capacity),
            fragments: Vec::with_capacity(capacity),
            paired,
        }
    }
    pub fn consider(&mut self, key: u128, f: Fragment<'_>) {
        if self.heap.len() < self.capacity {
            self.heap.push(SampleKey {
                key,
                ordinal: f.ordinal,
                slot: self.fragments.len(),
            });
            self.fragments.push(StoredFragment::new(f));
        } else if let Some(mut largest) = self.heap.peek_mut() {
            if (key, f.ordinal) < (largest.key, largest.ordinal) {
                self.fragments[largest.slot].replace(f);
                largest.key = key;
                largest.ordinal = f.ordinal;
            }
        }
    }
    pub fn len(&self) -> usize {
        self.heap.len()
    }
    pub fn finish(self) -> VecDeque<SelectedBatch> {
        drop(self.heap);
        let mut retained = self.fragments;
        retained.sort_unstable_by_key(|f| f.ordinal);
        let mut fragments = retained.into_iter();
        std::iter::from_fn(|| {
            let batch = SelectedBatch {
                fragments: fragments.by_ref().take(FASTQ_BATCH_RECORDS).collect(),
                paired: self.paired,
            };
            (!batch.fragments.is_empty()).then_some(batch)
        })
        .collect()
    }
}

pub(crate) struct SelectedBatch {
    fragments: Vec<StoredFragment>,
    paired: bool,
}
impl SelectedBatch {
    pub fn len(&self) -> usize {
        self.fragments.len()
    }
    pub fn fragments(&self) -> impl Iterator<Item = Fragment<'_>> {
        self.fragments.iter().map(|f| f.fragment(self.paired))
    }
    pub fn retain(&mut self, mut keep: impl FnMut(Fragment<'_>) -> bool) {
        self.fragments.retain(|f| keep(f.fragment(self.paired)));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::BTreeSet;

    #[test]
    fn binomial_design_meets_budget_and_precision_changes_capacity() {
        assert_eq!(
            sample_capacity(20560, AnalysisProfile::FROZEN, Precision::Standard).unwrap(),
            1_292_678
        );
        let mut previous = 0;
        for precision in [Precision::Fast, Precision::Standard, Precision::Sensitive] {
            let n = sample_capacity(20560, AnalysisProfile::FROZEN, precision).unwrap();
            assert!(n > previous);
            assert!(
                miss_probability_upper_bound(precision.minimum_fraction().next_down(), n as u64)
                    <= (0.05 / 20560.0_f64).next_down()
            );
            previous = n;
        }
        assert!(sample_capacity(0, AnalysisProfile::FROZEN, Precision::Standard).is_err());
    }

    fn selected(
        order: impl Iterator<Item = u64>,
        capacity: usize,
        candidates_only: bool,
    ) -> BTreeSet<u64> {
        let keys = FragmentKeys::new("profile", "index");
        let mut sampler = BottomKSampler::new(capacity, true);
        for ordinal in order {
            if candidates_only && ordinal % 3 != 0 {
                continue;
            }
            let id = format!("pair-{}", ordinal % 7);
            sampler.consider(
                keys.key(&id, ordinal),
                Fragment {
                    ordinal,
                    id: &id,
                    r1: &b"ACGTN"[..ordinal as usize % 6],
                    r2: Some(&b"TGCA"[..ordinal as usize % 5]),
                },
            );
            assert!(sampler.len() <= capacity);
        }
        let mut result = BTreeSet::new();
        let mut sample = sampler.finish();
        while let Some(batch) = sample.pop_front() {
            for f in batch.fragments() {
                assert_eq!(f.id, format!("pair-{}", f.ordinal % 7));
                assert_eq!(f.r1, &b"ACGTN"[..f.ordinal as usize % 6]);
                assert_eq!(f.r2, Some(&b"TGCA"[..f.ordinal as usize % 5]));
                result.insert(f.ordinal);
            }
        }
        result
    }

    #[test]
    fn bottom_k_matches_sorted_keys_and_full_contains_screen_candidates() {
        let keys = FragmentKeys::new("profile", "index");
        let mut offline: Vec<_> = (0..2000)
            .map(|i| (keys.key(&format!("pair-{}", i % 7), i), i))
            .collect();
        offline.sort_unstable();
        let expected: BTreeSet<_> = offline.iter().take(100).map(|&(_, i)| i).collect();
        let screen = selected(0..2000, 100, false);
        assert_eq!(screen, expected);
        assert_eq!(screen, selected((0..2000).rev(), 100, false));
        let full = selected((0..2000).rev(), 100, true);
        assert!(screen
            .into_iter()
            .filter(|i| i % 3 == 0)
            .all(|i| full.contains(&i)));
        assert_eq!(selected(0..12, 100, false).len(), 12);
        assert!(selected(0..0, 100, false).is_empty());
    }
}
