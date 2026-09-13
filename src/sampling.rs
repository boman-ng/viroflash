use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::fastq::{Fragment, FragmentBatch, FASTQ_BATCH_RECORDS};
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
    key: u128,
    ordinal: u64,
    id: String,
    r1: Vec<u8>,
    r2: Option<Vec<u8>>,
}
impl StoredFragment {
    fn new(key: u128, f: Fragment<'_>) -> Self {
        Self {
            key,
            ordinal: f.ordinal,
            id: f.id.into(),
            r1: f.r1.into(),
            r2: f.r2.map(Vec::from),
        }
    }
    fn replace(&mut self, key: u128, f: Fragment<'_>) {
        self.key = key;
        self.ordinal = f.ordinal;
        self.id.clear();
        self.id.push_str(f.id);
        self.r1.clear();
        self.r1.extend_from_slice(f.r1);
        match (self.r2.as_mut(), f.r2) {
            (Some(buffer), Some(read)) => {
                buffer.clear();
                buffer.extend_from_slice(read);
            }
            (_, read) => self.r2 = read.map(Vec::from),
        }
    }
    fn fragment(&self) -> Fragment<'_> {
        Fragment {
            ordinal: self.ordinal,
            id: &self.id,
            r1: &self.r1,
            r2: self.r2.as_deref(),
        }
    }
}
impl PartialEq for StoredFragment {
    fn eq(&self, other: &Self) -> bool {
        (self.key, self.ordinal) == (other.key, other.ordinal)
    }
}
impl Eq for StoredFragment {}
impl PartialOrd for StoredFragment {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
impl Ord for StoredFragment {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.key, self.ordinal).cmp(&(other.key, other.ordinal))
    }
}

pub(crate) struct BottomKSampler {
    capacity: usize,
    heap: BinaryHeap<StoredFragment>,
    paired: bool,
}
impl BottomKSampler {
    pub fn new(capacity: usize, paired: bool) -> Self {
        Self {
            capacity,
            heap: BinaryHeap::with_capacity(capacity),
            paired,
        }
    }
    pub fn consider(&mut self, key: u128, f: Fragment<'_>) {
        if self.heap.len() < self.capacity {
            self.heap.push(StoredFragment::new(key, f));
        } else if let Some(mut largest) = self.heap.peek_mut() {
            if (key, f.ordinal) < (largest.key, largest.ordinal) {
                largest.replace(key, f);
            }
        }
    }
    pub fn len(&self) -> usize {
        self.heap.len()
    }
    pub fn finish(self) -> SelectedFragments {
        let mut retained = self.heap.into_vec();
        retained.sort_unstable_by_key(|f| f.ordinal);
        SelectedFragments {
            fragments: retained.into_iter(),
            paired: self.paired,
        }
    }
}

pub(crate) struct SelectedFragments {
    fragments: std::vec::IntoIter<StoredFragment>,
    paired: bool,
}
impl SelectedFragments {
    pub fn next_batch(&mut self) -> Option<FragmentBatch> {
        let mut batch = FragmentBatch::new(self.paired);
        for f in self.fragments.by_ref().take(FASTQ_BATCH_RECORDS) {
            batch.push(f.fragment());
        }
        (batch.len() > 0).then_some(batch)
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
                    r1: b"ACGT",
                    r2: Some(b"TGCA"),
                },
            );
            assert!(sampler.len() <= capacity);
        }
        let mut result = BTreeSet::new();
        let mut sample = sampler.finish();
        while let Some(batch) = sample.next_batch() {
            for f in batch.fragments() {
                let f = f.unwrap();
                assert_eq!(f.r2, Some(b"TGCA".as_slice()));
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
