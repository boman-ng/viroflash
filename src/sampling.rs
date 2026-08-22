//! PE fragment 的确定性 bottom-k 抽样。

use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// 默认抽样容量：在 500 RPM 门槛处期望约 262 个 read-side 证据
/// （相对 Poisson SE 约 6.2%）；不是 panel 经验阈值。
pub const DEFAULT_SAMPLE_PAIRS: usize = 524_288;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct OwnedPair {
    pub ordinal: u64,
    pub qname: String,
    pub r1_seq: Vec<u8>,
    pub r2_seq: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SamplingResult {
    pub seen_pairs: u64,
    pub selected_pairs: u64,
    pub selected_passed_pairs: u64,
    pub passed_pairs: Vec<OwnedPair>,
    pub inclusion_probability: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct InclusionKey {
    digest: [u8; 32],
    /// 仅在摘要相同（重复 ID 或哈希碰撞）时参与排序。
    ordinal: u64,
}

#[derive(Debug)]
struct SelectedPair {
    key: InclusionKey,
    pair: Option<OwnedPair>,
}

impl PartialEq for SelectedPair {
    fn eq(&self, other: &Self) -> bool {
        self.key == other.key
    }
}

impl Eq for SelectedPair {}

impl PartialOrd for SelectedPair {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for SelectedPair {
    fn cmp(&self, other: &Self) -> Ordering {
        self.key.cmp(&other.key)
    }
}

fn inclusion_key(r1_id: &[u8], ordinal: u64) -> InclusionKey {
    InclusionKey {
        digest: *blake3::hash(r1_id).as_bytes(),
        ordinal,
    }
}

#[derive(Debug)]
pub struct PairReservoir {
    capacity: usize,
    seen_pairs: u64,
    selected: BinaryHeap<SelectedPair>,
}

impl PairReservoir {
    /// 创建固定容量的 bottom-k reservoir。
    pub fn new(capacity: usize) -> Result<Self, String> {
        if capacity == 0 {
            return Err("sampling capacity 必须大于 0".to_string());
        }
        Ok(Self {
            capacity,
            seen_pairs: 0,
            selected: BinaryHeap::new(),
        })
    }

    /// 观察一个已完成配对校验的 PE fragment。
    ///
    /// `ordinal` 必须是 fragment 在原始输入中的稳定唯一序号。通常只按 R1 ID
    /// 的 BLAKE3 摘要决定是否入样；仅当摘要相同（重复 ID 或哈希碰撞）时，
    /// `ordinal` 才影响选择。未通过预筛的入选项只保留键，不复制 ID 或序列。
    pub fn observe(
        &mut self,
        ordinal: u64,
        r1_id: &[u8],
        r1_seq: &[u8],
        r2_seq: &[u8],
        prescreen_pass: bool,
    ) -> Result<(), String> {
        let qname = std::str::from_utf8(r1_id)
            .map_err(|error| format!("R1 fragment ID 不是有效 UTF-8: {error}"))?;
        let next_seen = self
            .seen_pairs
            .checked_add(1)
            .ok_or_else(|| "fragment 计数溢出".to_string())?;
        let key = inclusion_key(r1_id, ordinal);
        let enters_sample = self.selected.len() < self.capacity
            || self
                .selected
                .peek()
                .is_some_and(|largest| key < largest.key);

        self.seen_pairs = next_seen;
        if !enters_sample {
            return Ok(());
        }

        if self.selected.len() == self.capacity {
            drop(self.selected.pop());
        }
        let pair = prescreen_pass.then(|| OwnedPair {
            ordinal,
            qname: qname.to_owned(),
            r1_seq: r1_seq.to_vec(),
            r2_seq: r2_seq.to_vec(),
        });
        self.selected.push(SelectedPair { key, pair });
        Ok(())
    }

    /// 完成抽样；通过预筛的入选 pair 按原始 ordinal 升序返回。
    pub fn finish(self) -> SamplingResult {
        let selected_pairs = self.selected.len() as u64;
        let mut passed_pairs: Vec<OwnedPair> = self
            .selected
            .into_iter()
            .filter_map(|selected| selected.pair)
            .collect();
        passed_pairs.sort_unstable_by_key(|pair| pair.ordinal);
        let selected_passed_pairs = passed_pairs.len() as u64;
        let inclusion_probability = if u128::from(self.seen_pairs) <= self.capacity as u128 {
            1.0
        } else {
            self.capacity as f64 / self.seen_pairs as f64
        };

        SamplingResult {
            seen_pairs: self.seen_pairs,
            selected_pairs,
            selected_passed_pairs,
            passed_pairs,
            inclusion_probability,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Clone, Copy)]
    struct Input {
        ordinal: u64,
        id: &'static [u8],
        r1_seq: &'static [u8],
        r2_seq: &'static [u8],
    }

    fn inputs() -> Vec<Input> {
        vec![
            Input {
                ordinal: 8,
                id: b"frag-8/1",
                r1_seq: b"ACGT",
                r2_seq: b"TGCA",
            },
            Input {
                ordinal: 1,
                id: b"frag-1/1",
                r1_seq: b"AAAA",
                r2_seq: b"TTTT",
            },
            Input {
                ordinal: 13,
                id: b"frag-13/1",
                r1_seq: b"CCCC",
                r2_seq: b"GGGG",
            },
            Input {
                ordinal: 3,
                id: b"frag-3/1",
                r1_seq: b"AGAG",
                r2_seq: b"CTCT",
            },
            Input {
                ordinal: 21,
                id: b"frag-21/1",
                r1_seq: b"ATAT",
                r2_seq: b"TATA",
            },
            Input {
                ordinal: 5,
                id: b"frag-5/1",
                r1_seq: b"CGCG",
                r2_seq: b"GCGC",
            },
            Input {
                ordinal: 34,
                id: b"frag-34/1",
                r1_seq: b"AACC",
                r2_seq: b"GGTT",
            },
            Input {
                ordinal: 2,
                id: b"frag-2/1",
                r1_seq: b"CACA",
                r2_seq: b"TGTG",
            },
        ]
    }

    fn expected_ordinals(records: &[Input], capacity: usize) -> Vec<u64> {
        let mut keyed: Vec<(InclusionKey, u64)> = records
            .iter()
            .map(|record| (inclusion_key(record.id, record.ordinal), record.ordinal))
            .collect();
        keyed.sort_unstable_by_key(|(key, _)| *key);
        keyed.truncate(capacity.min(keyed.len()));
        let mut ordinals: Vec<u64> = keyed.into_iter().map(|(_, ordinal)| ordinal).collect();
        ordinals.sort_unstable();
        ordinals
    }

    fn observe_all(
        records: &[Input],
        capacity: usize,
        passes: impl Fn(&Input) -> bool,
    ) -> SamplingResult {
        let mut reservoir = PairReservoir::new(capacity).expect("有效容量");
        for record in records {
            reservoir
                .observe(
                    record.ordinal,
                    record.id,
                    record.r1_seq,
                    record.r2_seq,
                    passes(record),
                )
                .expect("有效 fragment");
        }
        reservoir.finish()
    }

    fn result_ordinals(result: &SamplingResult) -> Vec<u64> {
        result
            .passed_pairs
            .iter()
            .map(|pair| pair.ordinal)
            .collect()
    }

    fn permutations(mut records: Vec<Input>) -> Vec<Vec<Input>> {
        fn generate(start: usize, records: &mut [Input], output: &mut Vec<Vec<Input>>) {
            if start == records.len() {
                output.push(records.to_vec());
                return;
            }
            for index in start..records.len() {
                records.swap(start, index);
                generate(start + 1, records, output);
                records.swap(start, index);
            }
        }

        let mut output = Vec::new();
        generate(0, &mut records, &mut output);
        output
    }

    #[test]
    fn zero_capacity_is_rejected() {
        let error = PairReservoir::new(0).expect_err("零容量必须报错");
        assert!(error.contains("capacity"));
    }

    #[test]
    fn non_utf8_id_is_rejected_without_counting_a_pair() {
        let mut reservoir = PairReservoir::new(2).unwrap();
        let error = reservoir
            .observe(0, &[0xff], b"AC", b"GT", true)
            .expect_err("非 UTF-8 ID 必须报错");
        assert!(error.contains("UTF-8"));

        let result = reservoir.finish();
        assert_eq!(result.seen_pairs, 0);
        assert_eq!(result.selected_pairs, 0);
        assert_eq!(result.selected_passed_pairs, 0);
        assert!(result.passed_pairs.is_empty());
        assert_eq!(result.inclusion_probability, 1.0);
    }

    #[test]
    fn exhaustive_small_inputs_match_bottom_k_oracle() {
        let universe = inputs();
        for mask in 1usize..(1usize << universe.len()) {
            let records: Vec<Input> = universe
                .iter()
                .enumerate()
                .filter(|(index, _)| mask & (1usize << index) != 0)
                .map(|(_, record)| *record)
                .collect();
            for capacity in 1..=records.len() + 1 {
                let result = observe_all(&records, capacity, |_| true);
                assert_eq!(
                    result_ordinals(&result),
                    expected_ordinals(&records, capacity),
                    "mask={mask:#x}, capacity={capacity}"
                );
                assert_eq!(result.selected_pairs, records.len().min(capacity) as u64);
            }
        }
    }

    #[test]
    fn unique_id_selection_is_invariant_to_input_permutation() {
        let records = inputs()[..5].to_vec();
        let expected = expected_ordinals(&records, 3);
        for permutation in permutations(records) {
            let result = observe_all(&permutation, 3, |_| true);
            assert_eq!(result_ordinals(&result), expected);
        }
    }

    #[test]
    fn n_at_most_k_selects_every_pair() {
        let records = inputs()[..3].to_vec();
        let result = observe_all(&records, 5, |_| true);
        let mut expected: Vec<u64> = records.iter().map(|record| record.ordinal).collect();
        expected.sort_unstable();

        assert_eq!(result.seen_pairs, 3);
        assert_eq!(result.selected_pairs, 3);
        assert_eq!(result.selected_passed_pairs, 3);
        assert_eq!(result_ordinals(&result), expected);
        assert_eq!(result.inclusion_probability, 1.0);
    }

    #[test]
    fn failed_prescreen_pair_occupies_sample_without_storing_sequences() {
        let mut records = inputs()[..4].to_vec();
        records.sort_unstable_by_key(|record| inclusion_key(record.id, record.ordinal));
        let failed_ordinal = records[0].ordinal;
        let expected_passed = records[1];
        records.reverse();

        let mut reservoir = PairReservoir::new(2).unwrap();
        for record in &records {
            reservoir
                .observe(
                    record.ordinal,
                    record.id,
                    record.r1_seq,
                    record.r2_seq,
                    record.ordinal != failed_ordinal,
                )
                .unwrap();
        }
        assert_eq!(
            reservoir
                .selected
                .iter()
                .filter(|entry| entry.pair.is_none())
                .count(),
            1
        );

        let result = reservoir.finish();
        assert_eq!(result.seen_pairs, 4);
        assert_eq!(result.selected_pairs, 2);
        assert_eq!(result.selected_passed_pairs, 1);
        assert_eq!(result.inclusion_probability, 0.5);
        assert_eq!(result.passed_pairs.len(), 1);
        assert_eq!(result.passed_pairs[0].ordinal, expected_passed.ordinal);
        assert_eq!(result.passed_pairs[0].r1_seq, expected_passed.r1_seq);
        assert_eq!(result.passed_pairs[0].r2_seq, expected_passed.r2_seq);
    }

    #[test]
    fn duplicate_ids_use_ordinal_as_a_stable_tie_break() {
        let records = vec![
            Input {
                ordinal: 9,
                id: b"duplicate/1",
                r1_seq: b"AAAA",
                r2_seq: b"TTTT",
            },
            Input {
                ordinal: 2,
                id: b"duplicate/1",
                r1_seq: b"CCCC",
                r2_seq: b"GGGG",
            },
            Input {
                ordinal: 5,
                id: b"duplicate/1",
                r1_seq: b"ACAC",
                r2_seq: b"TGTG",
            },
        ];

        for permutation in permutations(records) {
            let result = observe_all(&permutation, 2, |_| true);
            assert_eq!(result_ordinals(&result), vec![2, 5]);
        }
    }

    #[test]
    fn finish_orders_passed_pairs_by_original_ordinal() {
        let mut records = inputs()[..6].to_vec();
        records.sort_unstable_by_key(|record| std::cmp::Reverse(record.ordinal));
        let result = observe_all(&records, 4, |_| true);
        let ordinals = result_ordinals(&result);

        assert_eq!(ordinals, expected_ordinals(&records, 4));
        assert!(ordinals.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn repeated_input_produces_identical_result() {
        let records = inputs();
        let first = observe_all(&records, 5, |record| record.ordinal % 2 == 0);
        let second = observe_all(&records, 5, |record| record.ordinal % 2 == 0);
        assert_eq!(first, second);
    }

    #[test]
    fn empty_finish_does_not_invent_pairs() {
        let result = PairReservoir::new(4).unwrap().finish();
        assert_eq!(result.seen_pairs, 0);
        assert_eq!(result.selected_pairs, 0);
        assert_eq!(result.selected_passed_pairs, 0);
        assert!(result.passed_pairs.is_empty());
        assert_eq!(result.inclusion_probability, 1.0);
    }
}
