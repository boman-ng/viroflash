use crate::analysis_profile::AnalysisProfile;
use crate::fastq_input::InputCensus;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SamplingDesign {
    pub selection_probability: f64,
    pub minimum_relevant_fragments: u64,
}

pub fn derive_sampling_design(
    census: &InputCensus,
    target_family_size: usize,
    profile: AnalysisProfile,
) -> Result<SamplingDesign, String> {
    if target_family_size == 0 {
        return Err("Reference index contains no target groups".into());
    }
    let minimum_relevant_fragments = profile.minimum_relevant_fragments(census.fragments);
    let group_miss = profile.familywise_miss_probability / target_family_size as f64;
    let mut probability =
        -f64::exp_m1(group_miss.ln() / minimum_relevant_fragments as f64).min(1.0);
    while (1.0 - probability).powf(minimum_relevant_fragments as f64) > group_miss {
        probability = probability.next_up();
    }
    assert!(
        (1.0 - probability).powf(minimum_relevant_fragments as f64) <= group_miss,
        "rounded sampling probability exceeds the frozen miss budget"
    );
    Ok(SamplingDesign {
        selection_probability: probability,
        minimum_relevant_fragments,
    })
}

pub fn fragment_selection_key(
    profile_digest: &str,
    input_digest: &str,
    fragment_id: &str,
    ordinal: u64,
) -> [u8; 16] {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"viroflash-fragment-selection-v1\0");
    hasher.update(profile_digest.as_bytes());
    hasher.update(input_digest.as_bytes());
    hasher.update(&(fragment_id.len() as u64).to_be_bytes());
    hasher.update(fragment_id.as_bytes());
    hasher.update(&ordinal.to_be_bytes());
    let mut key = [0; 16];
    key.copy_from_slice(&hasher.finalize().as_bytes()[..16]);
    key
}

pub fn include_fragment(key: [u8; 16], probability: f64) -> bool {
    if probability >= 1.0 {
        return true;
    }
    let value = u128::from_be_bytes(key);
    value < (probability * u128::MAX as f64) as u128
}

#[cfg(test)]
mod tests {
    use super::*;

    fn census(n: u64) -> InputCensus {
        InputCensus {
            input_mode: "SE",
            fragments: n,
            input_digest: "i".repeat(64),
            read_ends_per_fragment: 1,
        }
    }

    #[test]
    fn probability_obeys_formula_and_boundaries() {
        let design =
            derive_sampling_design(&census(1_000_000), 20, AnalysisProfile::FROZEN).unwrap();
        let expected = 1.0 - (0.05_f64 / 20.0).powf(1.0 / 10.0);
        assert!((design.selection_probability - expected).abs() < 1e-15);
        let census_design = derive_sampling_design(&census(1), 1, AnalysisProfile::FROZEN).unwrap();
        assert!(census_design.selection_probability >= 0.95);
        assert!(1.0 - census_design.selection_probability <= 0.05);
    }

    #[test]
    fn decimal_boundary_uses_exact_rational_ceiling_and_meets_miss_budget() {
        let design =
            derive_sampling_design(&census(10_000_000), 20, AnalysisProfile::FROZEN).unwrap();
        assert_eq!(design.minimum_relevant_fragments, 100);
        let group_miss = AnalysisProfile::FROZEN.familywise_miss_probability / 20.0;
        let actual_miss = (1.0 - design.selection_probability).powi(100);
        assert!(actual_miss <= group_miss, "{actual_miss} > {group_miss}");
    }

    #[test]
    fn selection_is_deterministic_and_fragment_owned() {
        let key = fragment_selection_key("p", "i", "pair", 7);
        assert_eq!(key, fragment_selection_key("p", "i", "pair", 7));
        assert_eq!(include_fragment(key, 0.5), include_fragment(key, 0.5));
        assert!(include_fragment(key, 1.0));
    }
}
