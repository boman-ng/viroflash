use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

pub const PROFILE_BYTES: &[u8] = include_bytes!("analysis-profile.json");
pub const SDUST_WINDOW: usize = 64;
pub const SDUST_THRESHOLD: i64 = 20;

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
    fn constants_match_the_embedded_profile() {
        let document: serde_json::Value = serde_json::from_slice(PROFILE_BYTES).unwrap();
        let profile = AnalysisProfile::FROZEN;
        let parameters = &document["parameters"];
        assert_eq!(
            parameters["minimum_relevant_fraction"]["value"],
            profile.minimum_relevant_fraction
        );
        assert_eq!(
            parameters["familywise_miss_probability"]["value"],
            profile.familywise_miss_probability
        );
        assert_eq!(
            parameters["familywise_interval_error"]["value"],
            profile.familywise_interval_error
        );
        let alignment = &document["index_and_alignment"];
        assert_eq!(alignment["kmer_length"], profile.kmer_length);
        assert_eq!(
            alignment["occupied_window_bins"],
            profile.occupied_window_bins
        );
        assert_eq!(alignment["sdust_window"], SDUST_WINDOW);
        assert_eq!(alignment["sdust_threshold"], SDUST_THRESHOLD);
    }
}
