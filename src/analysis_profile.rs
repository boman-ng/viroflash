use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub const PROFILE_BYTES: &[u8] = include_bytes!("../evaluation/phase0/analysis-profile.json");
pub const MINIMUM_RELEVANT_FRACTION_NUMERATOR: u64 = 1;
pub const MINIMUM_RELEVANT_FRACTION_DENOMINATOR: u64 = 100_000;
pub const SDUST_WINDOW: usize = 64;
pub const SDUST_THRESHOLD: i64 = 20;

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AnalysisProfile {
    pub minimum_relevant_fraction: f64,
    pub familywise_miss_probability: f64,
    pub familywise_interval_error: f64,
    pub kmer_length: usize,
    pub sdust_window: usize,
    pub sdust_threshold: i64,
    pub occupied_window_bins: usize,
}

impl AnalysisProfile {
    pub const FROZEN: Self = Self {
        minimum_relevant_fraction: 1e-5,
        familywise_miss_probability: 0.05,
        familywise_interval_error: 0.05,
        kmer_length: 21,
        sdust_window: SDUST_WINDOW,
        sdust_threshold: SDUST_THRESHOLD,
        occupied_window_bins: 10,
    };

    pub fn digest(self) -> String {
        hex_sha256(PROFILE_BYTES)
    }

    pub fn minimum_relevant_fragments(self, input_fragments: u64) -> u64 {
        input_fragments
            .saturating_mul(MINIMUM_RELEVANT_FRACTION_NUMERATOR)
            .div_ceil(MINIMUM_RELEVANT_FRACTION_DENOMINATOR)
    }
}

pub(crate) fn hex_sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_profile_and_digest_match_committed_bytes() {
        let profile = AnalysisProfile::FROZEN;
        assert_eq!(profile.minimum_relevant_fraction, 1e-5);
        assert_eq!(profile.familywise_miss_probability, 0.05);
        assert_eq!(profile.familywise_interval_error, 0.05);
        assert_eq!(profile.kmer_length, 21);
        assert_eq!(profile.sdust_window, 64);
        assert_eq!(profile.sdust_threshold, 20);
        assert_eq!(
            profile.digest(),
            "7974800cbdb062b3b7c331bcdd5222cbea968f6842c63c27161be26686f01821"
        );
    }
}
