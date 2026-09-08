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
    let conservative_group_miss = group_miss.next_down();
    let mut probability =
        -f64::exp_m1(group_miss.ln() / minimum_relevant_fragments as f64).min(1.0);
    while miss_probability_upper_bound(probability, minimum_relevant_fragments)
        > conservative_group_miss
    {
        probability = probability.next_up();
    }
    assert!(
        miss_probability_upper_bound(probability, minimum_relevant_fragments)
            <= conservative_group_miss,
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
    match inclusion_threshold(probability) {
        None => true,
        Some(threshold) => u128::from_be_bytes(key) < threshold,
    }
}

fn miss_probability_upper_bound(probability: f64, relevant_fragments: u64) -> f64 {
    if probability >= 1.0 {
        return 0.0;
    }
    let mut exponent = relevant_fragments;
    let mut factor = (1.0 - probability).next_up();
    let mut result = 1.0;
    while exponent > 0 {
        if exponent & 1 == 1 {
            result = multiply_up(result, factor);
        }
        exponent >>= 1;
        if exponent > 0 {
            factor = multiply_up(factor, factor);
        }
    }
    result
}

fn multiply_up(left: f64, right: f64) -> f64 {
    let product = left * right;
    if product == 0.0 || product.is_infinite() {
        product
    } else {
        product.next_up()
    }
}

fn inclusion_threshold(probability: f64) -> Option<u128> {
    if probability >= 1.0 {
        return None;
    }
    if probability <= 0.0 {
        return Some(0);
    }
    let bits = probability.to_bits();
    let exponent_bits = ((bits >> 52) & 0x7ff) as i32;
    let fraction = bits & ((1_u64 << 52) - 1);
    let (significand, exponent) = if exponent_bits == 0 {
        (u128::from(fraction), -1022 - 52)
    } else {
        (
            u128::from(fraction | (1_u64 << 52)),
            exponent_bits - 1023 - 52,
        )
    };
    let shift = exponent + 128;
    if shift >= 0 {
        Some(significand << shift)
    } else {
        let right = (-shift) as u32;
        if right >= 128 {
            Some(1)
        } else {
            Some(significand.div_ceil(1_u128 << right))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone)]
    struct ExactNatural(Vec<u32>);

    impl ExactNatural {
        fn from_u128(value: u128) -> Self {
            let mut words = (0..4)
                .map(|shift| (value >> (shift * 32)) as u32)
                .collect::<Vec<_>>();
            while words.last() == Some(&0) && words.len() > 1 {
                words.pop();
            }
            Self(words)
        }

        fn multiply_u128(&mut self, factor: u128) {
            let factors = (0..4)
                .map(|shift| (factor >> (shift * 32)) as u32)
                .collect::<Vec<_>>();
            let mut products = vec![0_u128; self.0.len() + factors.len()];
            for (left_index, left) in self.0.iter().enumerate() {
                for (right_index, right) in factors.iter().enumerate() {
                    products[left_index + right_index] += u128::from(*left) * u128::from(*right);
                }
            }
            let mut carry = 0_u128;
            self.0 = products
                .into_iter()
                .map(|product| {
                    let value = product + carry;
                    carry = value >> 32;
                    value as u32
                })
                .collect();
            while carry > 0 {
                self.0.push(carry as u32);
                carry >>= 32;
            }
            while self.0.last() == Some(&0) && self.0.len() > 1 {
                self.0.pop();
            }
        }

        fn multiply_u64(&mut self, factor: u64) {
            let mut carry = 0_u128;
            for word in &mut self.0 {
                let value = u128::from(*word) * u128::from(factor) + carry;
                *word = value as u32;
                carry = value >> 32;
            }
            while carry > 0 {
                self.0.push(carry as u32);
                carry >>= 32;
            }
        }

        fn no_greater_than_power_of_two(&self, exponent: usize) -> bool {
            let highest = self.0.len() * 32 - self.0.last().unwrap().leading_zeros() as usize;
            if highest <= exponent {
                return true;
            }
            highest == exponent + 1
                && self.0.iter().enumerate().all(|(index, word)| {
                    let expected = if index == exponent / 32 {
                        1_u32 << (exponent % 32)
                    } else {
                        0
                    };
                    *word == expected
                })
        }
    }

    fn census(n: u64) -> InputCensus {
        InputCensus {
            input_mode: "SE",
            fragments: n,
            input_digest: "i".repeat(64),
            compressed_artifact_digest: "a".repeat(64),
            read_ends_per_fragment: 1,
        }
    }

    #[test]
    fn probability_obeys_formula_and_boundaries() {
        let design =
            derive_sampling_design(&census(1_000_000), 20, AnalysisProfile::FROZEN).unwrap();
        let expected = 1.0 - (0.05_f64 / 20.0).powf(1.0 / 10.0);
        assert!(design.selection_probability >= expected);
        let census_design = derive_sampling_design(&census(1), 1, AnalysisProfile::FROZEN).unwrap();
        assert!(census_design.selection_probability >= 0.95);
        assert!(1.0 - census_design.selection_probability <= 0.05);
    }

    #[test]
    fn phase5_production_probability_meets_exact_rational_miss_budget() {
        for (population, family_size) in [
            (1, 1),
            (100_000, 1),
            (100_001, 1),
            (1_000_000, 20_560),
            (66_500_000, 1),
        ] {
            let design =
                derive_sampling_design(&census(population), family_size, AnalysisProfile::FROZEN)
                    .unwrap();
            let threshold = inclusion_threshold(design.selection_probability).unwrap();
            let excluded = u128::MAX - threshold + 1;
            let mut miss_numerator = ExactNatural::from_u128(1);
            for _ in 0..design.minimum_relevant_fragments {
                miss_numerator.multiply_u128(excluded);
            }
            miss_numerator.multiply_u64(u64::try_from(20 * family_size).unwrap());
            assert!(
                miss_numerator.no_greater_than_power_of_two(
                    usize::try_from(128 * design.minimum_relevant_fragments).unwrap()
                ),
                "N={population} m={family_size} p={:.20}",
                design.selection_probability
            );
        }
    }

    #[test]
    fn phase5_bernoulli_selection_conditioned_on_realized_n_is_uniform() {
        let design =
            derive_sampling_design(&census(1_000_000), 20_560, AnalysisProfile::FROZEN).unwrap();
        let threshold = inclusion_threshold(design.selection_probability).unwrap();
        let excluded = u128::MAX - threshold + 1;
        for realized_sample in 0..=8 {
            let mut first_weight = None;
            for selected_positions in 0_u16..(1 << 8) {
                if selected_positions.count_ones() as usize != realized_sample {
                    continue;
                }
                let mut weight = ExactNatural::from_u128(1);
                for position in 0..8 {
                    weight.multiply_u128(if selected_positions & (1 << position) != 0 {
                        threshold
                    } else {
                        excluded
                    });
                }
                if let Some(expected) = &first_weight {
                    assert_eq!(weight.0, *expected, "realized n={realized_sample}");
                } else {
                    first_weight = Some(weight.0);
                }
            }
        }
    }

    #[test]
    fn selection_is_deterministic_and_fragment_owned() {
        let key = fragment_selection_key("p", "i", "pair", 7);
        assert_eq!(key, fragment_selection_key("p", "i", "pair", 7));
        assert_eq!(include_fragment(key, 0.5), include_fragment(key, 0.5));
        assert!(include_fragment(key, 1.0));
        let probability = derive_sampling_design(&census(66_500_000), 1, AnalysisProfile::FROZEN)
            .unwrap()
            .selection_probability;
        let threshold = inclusion_threshold(probability).unwrap();
        assert!(include_fragment((threshold - 1).to_be_bytes(), probability));
        assert!(!include_fragment(threshold.to_be_bytes(), probability));
    }
}
