use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub const PROFILE_BYTES: &[u8] = include_bytes!("../evaluation/phase0/analysis-profile.json");

#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct AnalysisProfile {
    pub minimum_relevant_fraction: f64,
    pub familywise_miss_probability: f64,
    pub familywise_interval_error: f64,
    pub kmer_length: usize,
    pub occupied_window_bins: usize,
}

impl AnalysisProfile {
    pub const FROZEN: Self = Self {
        minimum_relevant_fraction: 1e-5,
        familywise_miss_probability: 0.05,
        familywise_interval_error: 0.05,
        kmer_length: 21,
        occupied_window_bins: 10,
    };

    pub fn digest(self) -> String {
        hex_sha256(PROFILE_BYTES)
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
        assert_eq!(
            profile.digest(),
            "1eb5960774d118530f462d93f316e789cee836c51aedbbb973a11a6bac421088"
        );
    }
}
