use crate::fastq_input::InputCensus;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SamplingDesign {
    pub selection_probability: f64,
}

pub fn derive_sampling_design(
    census: &InputCensus,
    target_family_size: usize,
    delta: f64,
    beta: f64,
) -> Result<SamplingDesign, String> {
    if target_family_size == 0 {
        return Err("Reference index contains no target groups".into());
    }
    let minimum_fragments = (delta * census.fragments as f64).ceil().max(1.0);
    let group_miss = beta / target_family_size as f64;
    let probability = -f64::exp_m1(group_miss.ln() / minimum_fragments).min(1.0);
    Ok(SamplingDesign {
        selection_probability: probability,
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
        let design = derive_sampling_design(&census(1_000_000), 20, 1e-5, 0.05).unwrap();
        let expected = 1.0 - (0.05_f64 / 20.0).powf(1.0 / 10.0);
        assert!((design.selection_probability - expected).abs() < 1e-15);
        assert_eq!(
            derive_sampling_design(&census(1), 1, 1e-5, 0.05)
                .unwrap()
                .selection_probability,
            0.95
        );
    }

    #[test]
    fn selection_is_deterministic_and_fragment_owned() {
        let key = fragment_selection_key("p", "i", "pair", 7);
        assert_eq!(key, fragment_selection_key("p", "i", "pair", 7));
        assert_eq!(include_fragment(key, 0.5), include_fragment(key, 0.5));
        assert!(include_fragment(key, 1.0));
    }
}
