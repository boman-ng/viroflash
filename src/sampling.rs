use crate::profile::AnalysisProfile;

/// The minimum relevant fraction is measured within the Bloom candidate pool.
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

    fn from_str(value: &str) -> Result<Self, Self::Err> {
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
    pub candidate_fragments: u64,
    pub selection_probability: f64,
    pub minimum_relevant_fragments: u64,
}

pub fn derive_sampling_design(
    candidate_fragments: u64,
    target_family_size: usize,
    profile: AnalysisProfile,
    precision: Precision,
) -> Result<SamplingDesign, String> {
    if target_family_size == 0 {
        return Err("Reference index contains no target groups".into());
    }
    let minimum_relevant_fragments = candidate_fragments.div_ceil(precision.denominator());
    if candidate_fragments == 0 {
        return Ok(SamplingDesign {
            precision,
            candidate_fragments,
            selection_probability: 0.0,
            minimum_relevant_fragments: 0,
        });
    }
    let group_miss = profile.familywise_miss_probability / target_family_size as f64;
    let conservative_group_miss = group_miss.next_down();
    let mut probability =
        -f64::exp_m1(group_miss.ln() / minimum_relevant_fragments as f64).min(1.0);
    while miss_probability_upper_bound(probability, minimum_relevant_fragments)
        > conservative_group_miss
    {
        probability = probability.next_up();
    }
    Ok(SamplingDesign {
        precision,
        candidate_fragments,
        selection_probability: probability,
        minimum_relevant_fragments,
    })
}

/// The probability and input identity are fixed once the census has completed.
#[derive(Clone)]
pub(crate) struct FragmentSelector {
    prefix: blake3::Hasher,
    threshold: Option<u128>,
}

impl FragmentSelector {
    pub fn new(profile_digest: &str, input_digest: &str, probability: f64) -> Self {
        let mut prefix = blake3::Hasher::new();
        prefix.update(b"viroflash-fragment-selection-v1\0");
        prefix.update(profile_digest.as_bytes());
        prefix.update(input_digest.as_bytes());
        Self {
            prefix,
            threshold: inclusion_threshold(probability),
        }
    }

    pub fn includes(&self, fragment_id: &str, ordinal: u64) -> bool {
        let Some(threshold) = self.threshold else {
            return true;
        };
        let mut hasher = self.prefix.clone();
        hasher.update(&(fragment_id.len() as u64).to_be_bytes());
        hasher.update(fragment_id.as_bytes());
        hasher.update(&ordinal.to_be_bytes());
        let mut key = [0; 16];
        key.copy_from_slice(&hasher.finalize().as_bytes()[..16]);
        u128::from_be_bytes(key) < threshold
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

    #[test]
    fn sampling_probability_meets_exact_rational_miss_budget() {
        for (population, family_size) in [
            (1, 1),
            (100_000, 1),
            (100_001, 1),
            (1_000_000, 20_560),
            (66_500_000, 1),
        ] {
            let design = derive_sampling_design(
                population,
                family_size,
                AnalysisProfile::FROZEN,
                Precision::Standard,
            )
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
    fn precision_uses_candidate_population_and_nested_selection() {
        assert_eq!(
            Precision::Standard.minimum_fraction(),
            AnalysisProfile::FROZEN.minimum_relevant_fraction
        );
        let mut previous = std::collections::BTreeSet::new();
        for precision in [Precision::Fast, Precision::Standard, Precision::Sensitive] {
            let design =
                derive_sampling_design(10_000_000, 20_560, AnalysisProfile::FROZEN, precision)
                    .unwrap();
            assert_eq!(
                design.minimum_relevant_fragments,
                10_000_000_u64.div_ceil(precision.denominator())
            );
            assert!(
                miss_probability_upper_bound(
                    design.selection_probability,
                    design.minimum_relevant_fragments
                ) <= (0.05 / 20_560.0_f64).next_down()
            );
            let selector = FragmentSelector::new("profile", "input", design.selection_probability);
            let selected: std::collections::BTreeSet<_> = (0..10_000)
                .filter(|&ordinal| selector.includes("repeated", ordinal))
                .collect();
            assert!(previous.is_subset(&selected));
            previous = selected;
        }
        let empty = derive_sampling_design(0, 20_560, AnalysisProfile::FROZEN, Precision::Standard)
            .unwrap();
        assert_eq!(empty.selection_probability, 0.0);
        assert_eq!(empty.minimum_relevant_fragments, 0);
    }

    #[test]
    fn selection_preserves_fragment_identity() {
        // Frozen outputs of the v1 key contract, including UTF-8 IDs and repeated names.
        let selector = FragmentSelector::new("profile-digest", "input-digest", 0.5);
        let selected = (0..64)
            .filter(|&ordinal| selector.includes(&format!("read-{}-α", ordinal % 37), ordinal))
            .collect::<Vec<_>>();
        assert_eq!(
            selected,
            [
                3, 6, 8, 9, 10, 11, 12, 13, 15, 16, 19, 20, 21, 24, 27, 28, 29, 30, 31, 32, 34, 35,
                37, 41, 48, 55, 56, 57, 58, 61, 63
            ]
        );
        assert!(!FragmentSelector::new("p", "i", 0.0).includes("read", 0));
        assert!(FragmentSelector::new("p", "i", 1.0).includes("read", 0));
        assert_eq!(inclusion_threshold(0.5), Some(1_u128 << 127));
        assert_eq!(inclusion_threshold(f64::from_bits(1)), Some(1));
    }
}
