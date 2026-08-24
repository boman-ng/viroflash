//! fragment 的确定性 bottom-k 抽样与发现/验证折分。

use std::cmp::Ordering;
use std::collections::BinaryHeap;

/// 默认抽样容量。它限制任意输入规模下保留的序列内存与后续比对工作量；
/// `N <= K` 时退化为全量处理，不是 panel 或阳性判定阈值。
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
    pub pairs: Vec<OwnedPair>,
    pub inclusion_probability: f64,
}

/// 同一 fragment 的两个 read-end 必须进入同一折，防止一端参与发现、另一端又
/// 参与验证。ordinal 让重复 qname 的独立输入记录不会永久绑定在同一折。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EvidenceFold {
    Discovery,
    Validation,
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
    pair: OwnedPair,
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
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"viroflash/bottom-k/v1\0");
    hasher.update(&(r1_id.len() as u64).to_le_bytes());
    hasher.update(r1_id);
    hasher.update(&ordinal.to_le_bytes());
    InclusionKey {
        digest: *hasher.finalize().as_bytes(),
        // BLAKE3 碰撞时仍以稳定唯一序号给出全序。
        ordinal,
    }
}

pub fn evidence_fold(qname: &str, ordinal: u64) -> EvidenceFold {
    let mut hasher = blake3::Hasher::new();
    hasher.update(b"viroflash/evidence-fold/v1\0");
    hasher.update(&(qname.len() as u64).to_le_bytes());
    hasher.update(qname.as_bytes());
    hasher.update(&ordinal.to_le_bytes());
    if hasher.finalize().as_bytes()[0] & 1 == 0 {
        EvidenceFold::Discovery
    } else {
        EvidenceFold::Validation
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
    /// `ordinal` 必须是 fragment 在原始输入中的稳定唯一序号。抽样层不承担
    /// k-mer 门控：先完成 bottom-k，再只门控最终 K 项，避免对全部输入随机访问
    /// 大 Bloom。R1 ID 与 ordinal 共同进入域分离 BLAKE3，重复 ID 仍独立抽样。
    pub fn observe(
        &mut self,
        ordinal: u64,
        r1_id: &[u8],
        r1_seq: &[u8],
        r2_seq: &[u8],
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
        let pair = OwnedPair {
            ordinal,
            qname: qname.to_owned(),
            r1_seq: r1_seq.to_vec(),
            r2_seq: r2_seq.to_vec(),
        };
        self.selected.push(SelectedPair { key, pair });
        Ok(())
    }

    /// 完成抽样；入选 pair 按原始 ordinal 升序返回。
    pub fn finish(self) -> SamplingResult {
        let selected_pairs = self.selected.len() as u64;
        let mut pairs: Vec<OwnedPair> = self
            .selected
            .into_iter()
            .map(|selected| selected.pair)
            .collect();
        pairs.sort_unstable_by_key(|pair| pair.ordinal);
        let inclusion_probability = if u128::from(self.seen_pairs) <= self.capacity as u128 {
            1.0
        } else {
            self.capacity as f64 / self.seen_pairs as f64
        };

        SamplingResult {
            seen_pairs: self.seen_pairs,
            selected_pairs,
            pairs,
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

    fn observe_all(records: &[Input], capacity: usize) -> SamplingResult {
        let mut reservoir = PairReservoir::new(capacity).expect("有效容量");
        for record in records {
            reservoir
                .observe(record.ordinal, record.id, record.r1_seq, record.r2_seq)
                .expect("有效 fragment");
        }
        reservoir.finish()
    }

    fn result_ordinals(result: &SamplingResult) -> Vec<u64> {
        result.pairs.iter().map(|pair| pair.ordinal).collect()
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
            .observe(0, &[0xff], b"AC", b"GT")
            .expect_err("非 UTF-8 ID 必须报错");
        assert!(error.contains("UTF-8"));

        let result = reservoir.finish();
        assert_eq!(result.seen_pairs, 0);
        assert_eq!(result.selected_pairs, 0);
        assert!(result.pairs.is_empty());
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
                let result = observe_all(&records, capacity);
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
            let result = observe_all(&permutation, 3);
            assert_eq!(result_ordinals(&result), expected);
        }
    }

    #[test]
    fn n_at_most_k_selects_every_pair() {
        let records = inputs()[..3].to_vec();
        let result = observe_all(&records, 5);
        let mut expected: Vec<u64> = records.iter().map(|record| record.ordinal).collect();
        expected.sort_unstable();

        assert_eq!(result.seen_pairs, 3);
        assert_eq!(result.selected_pairs, 3);
        assert_eq!(result.pairs.len(), 3);
        assert_eq!(result_ordinals(&result), expected);
        assert_eq!(result.inclusion_probability, 1.0);
    }

    #[test]
    fn selected_records_retain_both_sequences_for_post_sampling_gate() {
        let records = inputs()[..4].to_vec();
        let result = observe_all(&records, 2);
        assert_eq!(result.seen_pairs, 4);
        assert_eq!(result.selected_pairs, 2);
        assert_eq!(result.pairs.len(), 2);
        assert!(result
            .pairs
            .iter()
            .all(|pair| !pair.r1_seq.is_empty() && !pair.r2_seq.is_empty()));
    }

    #[test]
    fn duplicate_ids_are_hashed_as_distinct_stable_records() {
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

        let expected = expected_ordinals(&records, 2);
        for permutation in permutations(records) {
            let result = observe_all(&permutation, 2);
            assert_eq!(result_ordinals(&result), expected);
        }
    }

    #[test]
    fn finish_orders_passed_pairs_by_original_ordinal() {
        let mut records = inputs()[..6].to_vec();
        records.sort_unstable_by_key(|record| std::cmp::Reverse(record.ordinal));
        let result = observe_all(&records, 4);
        let ordinals = result_ordinals(&result);

        assert_eq!(ordinals, expected_ordinals(&records, 4));
        assert!(ordinals.windows(2).all(|window| window[0] < window[1]));
    }

    #[test]
    fn repeated_input_produces_identical_result() {
        let records = inputs();
        let first = observe_all(&records, 5);
        let second = observe_all(&records, 5);
        assert_eq!(first, second);
    }

    #[test]
    fn empty_finish_does_not_invent_pairs() {
        let result = PairReservoir::new(4).unwrap().finish();
        assert_eq!(result.seen_pairs, 0);
        assert_eq!(result.selected_pairs, 0);
        assert!(result.pairs.is_empty());
        assert_eq!(result.inclusion_probability, 1.0);
    }

    #[test]
    fn fold_is_deterministic_and_shared_by_both_read_ends() {
        let first = evidence_fold("fragment/1", 42);
        assert_eq!(first, evidence_fold("fragment/1", 42));
        // 调用方只按 fragment 调一次，因此两个 read-end 共享这一结果。
        assert!(matches!(
            first,
            EvidenceFold::Discovery | EvidenceFold::Validation
        ));
    }

    #[test]
    fn duplicate_qnames_are_independent_records() {
        let keys = (0..64)
            .map(|ordinal| inclusion_key(b"duplicate/1", ordinal).digest)
            .collect::<std::collections::HashSet<_>>();
        assert_eq!(keys.len(), 64);

        let folds = (0..64)
            .map(|ordinal| evidence_fold("duplicate/1", ordinal))
            .collect::<std::collections::HashSet<_>>();
        assert_eq!(folds.len(), 2);
    }
}
