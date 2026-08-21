//! k-mer 预筛：单列 Bloom 词典 + 比例门 + SDUST 掩蔽。
//!
//! - 词典：COBS 式单列 Bloom（f=0.1, h=1），使用 canonical 编码
//!   （min(正链, 反补链)）；位容量取 2 的幂，以掩码代替取模。按总碱基数
//!   上界预分配，重复插入只降低有效 FPR。
//! - 门控：`hits ≥ max(10, ⌈0.35·n_q_eff⌉)`；绝对命中下限与 KMCP
//!   `--min-kmers=10` 一致。
//! - 早拒：掩蔽只会减少 hits，raw hits < 10 的 read 可在 SDUST 前直接拒绝，
//!   其结果与完整路径一致。
//! - 掩蔽：按 minimap2 `sdust.c` 实现对称 DUST，使用相同的 W=64/T=20。
//!   掩蔽结果为区间列表；门控扫描使用双指针跳过起点落在区间内的 k-mer，
//!   read 本身保持不变并继续进入后续比对。热路径通过栈上环形队列和
//!   `GateScratch` 复用缓冲区。
//!
//! 算法出处：COBS (arXiv:1905.09624)；KMCP
//! (DOI 10.1093/bioinformatics/btac845)；DUST/SDUST (PMID 16796549)；
//! splitmix64 (DOI 10.1145/2660193.2660195)。

use crate::reference::Contig;

pub const DEFAULT_K: usize = 21;
/// 2-bit 编码上限：k>31 时 u64 装不下（2k>62 位）。入口校验见 lib.rs run_pipeline。
pub const K_MAX: usize = 31;

/// 单 k-mer 假阳性率（COBS/KMCP 风格 f=0.1 + 比例门）。
pub const GATE_FPR: f64 = 0.1;
/// 绝对命中下限（与 KMCP `--min-kmers=10` 一致）。
pub const GATE_MIN_HITS: u64 = 10;
/// 命中率门 0.35 = 7/20，为低一致度目标保留更多候选 reads。
pub const GATE_HIT_FRAC_NUM: u64 = 7;
pub const GATE_HIT_FRAC_DEN: u64 = 20;
/// SDUST 窗口/阈值（与 minimap2 sdust.c 相同）。
pub const SDUST_W: usize = 64;
pub const SDUST_T: i64 = 20;

pub fn dna_bits(c: u8) -> Option<u8> {
    match c.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub fn encode_kmer(seq: &[u8]) -> Option<u64> {
    if seq.is_empty() || seq.len() > K_MAX {
        return None;
    }
    let mut code = 0u64;
    for &c in seq {
        code = (code << 2) | dna_bits(c)? as u64;
    }
    Some(code)
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

/// splitmix64 终结器（Steele et al. 2014；纯函数，k-mer 码 → 均匀 u64）。
/// 诱饵生成也使用该函数派生确定性随机种子。
pub(crate) fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

// ---------------------------------------------------------------------------
// 单列 Bloom 词典
// ---------------------------------------------------------------------------

/// 单列 Bloom：`words` 打包 m 位（m 为 2 的幂，`mask = m−1` 代替取模）。
pub struct KmerBloom {
    words: Vec<u64>,
    mask: u64,
    pub k: usize,
    /// 插入次数（非去重计数；审计用）。
    pub n_inserted: u64,
}

impl KmerBloom {
    /// 从目标+诱饵 contigs 构建 canonical k-mer Bloom（诱饵 reads 须能通过比例门）。
    /// 无有效 k-mer（全部过短/全 N）返回 None。位容量 m = pow2(⌈Σ碱基数 × (−1/ln(1−fpr))⌉)，
    /// fpr 为上界（canonical 去重与含 N 窗口跳过后实际 FPR 更低）。
    pub fn build(contigs: &[&Contig], k: usize, fpr: f64) -> Option<Self> {
        if k == 0 || k > K_MAX {
            return None;
        }
        let bases: u64 = contigs.iter().map(|c| c.seq.len() as u64).sum();
        if bases == 0 {
            return None;
        }
        let m_min_bits = (bases as f64 * (-1.0 / (1.0 - fpr).ln())).ceil().max(1.0) as u64;
        let m = m_min_bits.max(64).next_power_of_two();
        let words = vec![0u64; (m as usize) / 64];
        let mut bloom = Self {
            words,
            mask: m - 1,
            k,
            n_inserted: 0,
        };
        for c in contigs {
            let rc = reverse_complement(&c.seq);
            let n = c.seq.len();
            // 窗口起点 i ∈ 0..n−k+1（n<k 时为空）；反补链镜像窗口 rc[n−k−i..n−i]。
            // canonical = min(正链码, 反补码)，正链扫一遍即可覆盖两条链
            // （词典减半，同 KMCP/COBS）。
            for i in 0..n.saturating_sub(k - 1) {
                let fwd = encode_kmer(&c.seq[i..i + k]);
                let rev = encode_kmer(&rc[n - k - i..n - i]);
                let (Some(f), Some(r)) = (fwd, rev) else {
                    continue;
                };
                let canon = f.min(r);
                let idx = (splitmix64(canon) & bloom.mask) as usize;
                bloom.words[idx >> 6] |= 1 << (idx & 63);
                bloom.n_inserted += 1;
            }
        }
        Some(bloom)
    }

    /// 探测一个（canonical）k-mer 码。无假阴性；假阳性率 ≤ 构建时的 fpr 上界。
    #[inline]
    pub fn probe(&self, code: u64) -> bool {
        let idx = (splitmix64(code) & self.mask) as usize;
        (self.words[idx >> 6] >> (idx & 63)) & 1 == 1
    }

    /// 实际填充位比例（审计/校准用；构建后应 ≈ 1−e^{−n/m}）。
    pub fn fill_frac(&self) -> f64 {
        let set: u64 = self.words.iter().map(|w| u64::from(w.count_ones())).sum();
        set as f64 / (self.words.len() * 64) as f64
    }
}

// ---------------------------------------------------------------------------
// SDUST 低复杂度掩蔽（minimap2 sdust.c 直译）
// ---------------------------------------------------------------------------

/// 单个 perfect interval（find_perfect 维护，按 start 降序 + 密度剪枝）。
struct PerfectIntv {
    start: usize,
    finish: usize,
    r: i64,
    l: i64,
}

/// `save_masked_regions` 直译：把已滑出窗口的区间定稿（与上一区间重叠/相邻
/// 则合并），并移除 start < 窗口起点的区间（P 按 start 降序 → 尾部连续段）。
fn save_masked_regions(res: &mut Vec<(usize, usize)>, p: &mut Vec<PerfectIntv>, start: usize) {
    if p.is_empty() || p[p.len() - 1].start >= start {
        return;
    }
    let last = &p[p.len() - 1];
    if let Some(prev) = res.last_mut() {
        if last.start <= prev.1 {
            prev.1 = prev.1.max(last.finish);
        } else {
            res.push((last.start, last.finish));
        }
    } else {
        res.push((last.start, last.finish));
    }
    let mut keep = 0usize;
    for i in (0..p.len()).rev() {
        if p[i].start >= start {
            keep = i + 1;
            break;
        }
    }
    p.truncate(keep);
}

/// 栈上环形队列（minimap2 kdq 语义）：容量 128（2 的幂，掩码代替取模），
/// 杜绝每 read 堆分配。队列长度上限 = W−2 = 62（pop 条件 len ≥ W−2 后 push）。
const RING_CAP: usize = 128;
const RING_MASK: usize = RING_CAP - 1;

struct Ring {
    buf: [u32; RING_CAP],
    head: usize,
    len: usize,
}

impl Ring {
    fn new() -> Self {
        Ring {
            buf: [0; RING_CAP],
            head: 0,
            len: 0,
        }
    }
    #[inline]
    fn len(&self) -> usize {
        self.len
    }
    #[inline]
    fn push(&mut self, t: u32) {
        self.buf[(self.head + self.len) & RING_MASK] = t;
        self.len += 1;
    }
    #[inline]
    fn pop_front(&mut self) -> u32 {
        debug_assert!(self.len > 0);
        let t = self.buf[self.head];
        self.head = (self.head + 1) & RING_MASK;
        self.len -= 1;
        t
    }
    /// 从队首数第 i 个元素（kdq_at 语义）。
    #[inline]
    fn at(&self, i: usize) -> u32 {
        self.buf[(self.head + i) & RING_MASK]
    }
}

/// `shift_window` 直译：滑入当前 triplet，维护 rw/rv 与对称修剪（SDUST 核心：
/// 某 triplet 计数超 2T/10 时从左裁到该 triplet 上一次出现，保证掩蔽与扫描方向无关）。
// 参数逐项对应 SDUST 原算法状态，封装为结构体会掩盖直译关系。
#[allow(clippy::too_many_arguments)]
fn shift_window(
    que: &mut Ring,
    t3: u32,
    t: i64,
    w: usize,
    lcap: &mut usize,
    rw: &mut i64,
    rv: &mut i64,
    cw: &mut [i64; 64],
    cv: &mut [i64; 64],
) {
    if que.len() >= w - 2 {
        // W - SD_WLEN + 1 = 62
        let s = (que.pop_front() as usize) & 63; // t3 为 6-bit，掩码消 bounds check
        cw[s] -= 1;
        *rw -= cw[s];
        if *lcap > que.len() {
            *lcap -= 1;
            cv[s] -= 1;
            *rv -= cv[s];
        }
    }
    que.push(t3);
    *lcap += 1;
    let s = t3 as usize; // t3 已 & 0x3F，编译器可证 <64
    *rw += cw[s];
    cw[s] += 1;
    *rv += cv[s];
    cv[s] += 1;
    if cv[s] * 10 > t << 1 {
        loop {
            let u = (que.at(que.len() - *lcap) as usize) & 63;
            cv[u] -= 1;
            *rv -= cv[u];
            *lcap -= 1;
            if u == s {
                break;
            }
        }
    }
}

/// `find_perfect` 直译：窗口左侧延伸扫描（size−L−1 .. 0），密度剪枝维护
/// perfect 区间集（P 按 start 降序，同 end = 窗口尾）。
fn find_perfect(
    p: &mut Vec<PerfectIntv>,
    que: &Ring,
    t: i64,
    start: usize,
    lcap: usize,
    rv: i64,
    cv: &[i64; 64],
) {
    let mut c = *cv;
    let mut r = rv;
    let mut max_r = 0i64;
    let mut max_l = 0i64;
    let size = que.len();
    for i in (0..size.saturating_sub(lcap)).rev() {
        let t3 = (que.at(i) as usize) & 63;
        r += c[t3];
        c[t3] += 1;
        let (new_r, new_l) = (r, (size - i - 1) as i64);
        if new_r * 10 > t * new_l {
            let mut j = 0usize;
            while j < p.len() && p[j].start >= i + start {
                if max_r == 0 || p[j].r * max_l > max_r * p[j].l {
                    max_r = p[j].r;
                    max_l = p[j].l;
                }
                j += 1;
            }
            if max_r == 0 || new_r * max_l >= max_r * new_l {
                max_r = new_r;
                max_l = new_l;
                p.insert(
                    j,
                    PerfectIntv {
                        start: i + start,
                        finish: size + 2 + start, // SD_WLEN−1 = 2
                        r: new_r,
                        l: new_l,
                    },
                );
            }
        }
    }
}

/// `sdust_core` 直译：把掩蔽区间 [start, finish) 追加写入 `out`（清空后复用，
/// 热路径零分配）。N 断开连续段（C 版在 N 处不清队列，此处照搬保持语义一致）。
pub fn sdust_intervals_into(seq: &[u8], w: usize, t: i64, out: &mut Vec<(usize, usize)>) {
    out.clear();
    let mut p: Vec<PerfectIntv> = Vec::new(); // 通常为空：无 perfect 区间时零分配
    let mut que = Ring::new();
    let mut cw = [0i64; 64];
    let mut cv = [0i64; 64];
    let (mut rv, mut rw) = (0i64, 0i64);
    let mut lcap = 0usize;
    let mut l = 0usize;
    let mut t3 = 0u32;

    for i in 0..=seq.len() {
        let b = if i < seq.len() {
            dna_bits(seq[i]).unwrap_or(4)
        } else {
            4
        };
        if b < 4 {
            l += 1;
            t3 = ((t3 << 2) | b as u32) & 0x3F;
            if l >= 3 {
                let start = l.saturating_sub(w) + (i + 1 - l);
                save_masked_regions(out, &mut p, start);
                shift_window(
                    &mut que, t3, t, w, &mut lcap, &mut rw, &mut rv, &mut cw, &mut cv,
                );
                if rw * 10 > lcap as i64 * t {
                    find_perfect(&mut p, &que, t, start, lcap, rv, &cv);
                }
            }
        } else {
            // N 或序列尾：清空未定稿区间（start 逐步推进，每次定稿一个）
            let mut start = l.saturating_sub(w.saturating_sub(1)) + (i + 1 - l);
            while !p.is_empty() {
                save_masked_regions(out, &mut p, start);
                start += 1;
            }
            l = 0;
            t3 = 0;
        }
    }
}

/// 掩蔽区间列表（分配版；测试/审计用，热路径用 `sdust_intervals_into`）。
pub fn sdust_intervals(seq: &[u8], w: usize, t: i64) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    sdust_intervals_into(seq, w, t, &mut out);
    out
}

/// 每碱基掩码：掩蔽区间内 base 置 true（finish 截断到序列长）。
pub fn sdust_mask(seq: &[u8], w: usize, t: i64) -> Vec<bool> {
    let mut mask = vec![false; seq.len()];
    for (s, f) in sdust_intervals(seq, w, t) {
        let f = f.min(seq.len());
        for b in mask.iter_mut().take(f).skip(s) {
            *b = true;
        }
    }
    mask
}

// ---------------------------------------------------------------------------
// 比例门
// ---------------------------------------------------------------------------

/// 比例门：`hits ≥ max(GATE_MIN_HITS, ⌈0.35·n_q_eff⌉)`（整数精确）。
/// n_q_eff = 0（无未掩蔽 k-mer 窗口）恒不过门。
pub fn gate_passes(hits: u64, n_eff: u64) -> bool {
    if n_eff == 0 {
        return false;
    }
    let frac = (GATE_HIT_FRAC_NUM * n_eff).div_ceil(GATE_HIT_FRAC_DEN);
    hits >= GATE_MIN_HITS.max(frac)
}

/// 滚动双码扫描（Bloom 查询 + SDUST 掩蔽集成）：正链码左移、反补码右移（低位=最近碱基），
/// canonical = min(正链, 反补)。返回 (命中 k-mer 数, 有效 k-mer 数 n_q_eff)。
/// `intervals` = sdust 掩蔽区间（按 start 递增、互不重叠）；k-mer **起点**落在
/// 区间内（或含 N）不探测、不计 n_eff。掩蔽区间内的 read 本身保留。
pub fn read_gate_hits(
    seq: &[u8],
    k: usize,
    bloom: &KmerBloom,
    intervals: &[(usize, usize)],
) -> (u64, u64) {
    let (mut hits, mut n_eff) = (0u64, 0u64);
    if k == 0 || k > K_MAX || seq.len() < k {
        return (0, 0);
    }
    let kmask = (1u64 << (2 * k)) - 1;
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut filled = 0usize;
    let mut it = 0usize; // 双指针：当前候选掩蔽区间
    for (i, &c) in seq.iter().enumerate() {
        match dna_bits(c) {
            Some(b) => {
                fwd = ((fwd << 2) | b as u64) & kmask;
                rev = (rev >> 2) | (((3 - b) as u64) << (2 * k - 2));
                if filled < k {
                    filled += 1;
                }
                if filled == k {
                    let start = i - (k - 1);
                    // 推进到可能覆盖 start 的区间（区间按 start 递增）
                    while it < intervals.len() && intervals[it].1 <= start {
                        it += 1;
                    }
                    let masked = it < intervals.len() && intervals[it].0 <= start;
                    if !masked {
                        n_eff += 1;
                        if bloom.probe(fwd.min(rev)) {
                            hits += 1;
                        }
                    }
                }
            }
            None => {
                fwd = 0;
                rev = 0;
                filled = 0;
            }
        }
    }
    (hits, n_eff)
}

/// 门控 scratch：worker 线程级复用，热路径零分配。
#[derive(Default)]
pub struct GateScratch {
    pub intervals: Vec<(usize, usize)>,
}

/// 单 read 完整门控（热路径版本）：
/// 早拒（raw hits < 10 恒拒，免 SDUST）→ SDUST 掩蔽（复用 scratch）→ 命中扫描 → 比例门。
pub fn read_passes_gate_scratched(
    seq: &[u8],
    k: usize,
    bloom: &KmerBloom,
    scratch: &mut GateScratch,
) -> bool {
    // 早拒依据：掩蔽只会把 k-mer 移出计数（hits、n_eff 同减），故 masked hits ≤ raw hits；
    // 门要求 hits ≥ GATE_MIN_HITS，raw < GATE_MIN_HITS 时无论掩蔽与否恒拒。
    let (raw_hits, _) = read_gate_hits(seq, k, bloom, &[]);
    if raw_hits < GATE_MIN_HITS {
        return false;
    }
    sdust_intervals_into(seq, SDUST_W, SDUST_T, &mut scratch.intervals);
    let (hits, n_eff) = read_gate_hits(seq, k, bloom, &scratch.intervals);
    gate_passes(hits, n_eff)
}

/// 单 read 完整门控（便捷版；测试/审计用，热路径用 `read_passes_gate_scratched`）。
pub fn read_passes_gate(seq: &[u8], k: usize, bloom: &KmerBloom) -> bool {
    let mut scratch = GateScratch::default();
    read_passes_gate_scratched(seq, k, bloom, &mut scratch)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn contig(seq: &[u8]) -> Contig {
        Contig {
            name: "target_0".into(),
            role: crate::reference::Role::Target,
            seq: seq.to_vec(),
            gc_frac: 0.5,
        }
    }

    /// xorshift64*（Vigna 2016）——确定性伪随机序列（测试用）。
    fn xorshift64_star(mut s: u64) -> u64 {
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        s.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn pseudo_random_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut s = seed | 1;
        let mut out = Vec::with_capacity(len);
        for _ in 0..len {
            s = xorshift64_star(s);
            out.push(b"ACGT"[(s >> 59) as usize % 4]);
        }
        out
    }

    /// 随机 contig（非低复杂度，SDUST 不掩蔽；确定性 seed）。
    fn random_contig(len: usize, seed: u64) -> Contig {
        contig(&pseudo_random_seq(len, seed))
    }

    #[test]
    fn kmer_encode_and_rc() {
        assert_eq!(dna_bits(b'a'), Some(0));
        assert_eq!(dna_bits(b'N'), None);
        let fwd = encode_kmer(b"ACGT").unwrap();
        let rc = encode_kmer(&reverse_complement(b"ACGT")).unwrap();
        assert_eq!(fwd, rc); // ACGT 的回文反向互补
        let a = encode_kmer(b"AAAA").unwrap();
        let t = encode_kmer(b"TTTT").unwrap();
        assert_ne!(a, t);
    }

    #[test]
    fn bloom_all_contig_kmers_probe_true() {
        // 无假阴性：正链/反补链全部 k-mer 的 canonical 码必命中。
        let c = random_contig(60, 1);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let rc = reverse_complement(&c.seq);
        for i in 0..60 - 6 + 1 {
            let f = encode_kmer(&c.seq[i..i + 6]).unwrap();
            let r = encode_kmer(&rc[60 - 6 - i..60 - i]).unwrap();
            assert!(bloom.probe(f.min(r)), "k-mer @{i} 应命中");
        }
    }

    #[test]
    fn bloom_fpr_within_bound() {
        // 随机 k-mer 假阳性率 ≤ 构建 fpr 上界（pow2 取整后 f_eff 更低，断言放宽到 0.2）。
        let c = contig(b"ACGTACGTACGTACGTACGTAC");
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let ff = bloom.fill_frac();
        assert!(ff > 0.0 && ff < 1.0, "填充率异常: {ff}");
        let n = 50_000;
        let mut fp = 0u64;
        let mut s = 99u64;
        for _ in 0..n {
            s = xorshift64_star(s);
            let code = (s >> 34) & ((1u64 << 12) - 1); // 随机 6-mer 2-bit 码
            if bloom.probe(code) {
                fp += 1;
            }
        }
        assert!((fp as f64) < 0.2 * n as f64, "FP 率超界: {fp}/{n}");
    }

    #[test]
    fn gate_passes_thresholds() {
        // ⌈0.35·100⌉=35，下限 10
        assert!(!gate_passes(34, 100));
        assert!(gate_passes(35, 100));
        // max(10, ⌈0.35·50⌉=18) = 18
        assert!(!gate_passes(9, 50));
        assert!(gate_passes(18, 50));
        // max(10, ⌈0.35·10⌉=4) = 10
        assert!(!gate_passes(9, 10));
        assert!(gate_passes(10, 10));
        // 退化
        assert!(!gate_passes(0, 0));
        assert!(!gate_passes(5, 3));
    }

    #[test]
    fn read_passes_gate_fwd_and_rc_reads() {
        // 随机 contig（非低复杂度，SDUST 不掩蔽）：正链 read n_eff=hits=55，
        // 阈值 max(10, ⌈0.45·55⌉=25) → 过门；反补链 read canonical 覆盖同样过门。
        let c = random_contig(60, 1);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        assert!(read_passes_gate(&c.seq, 6, &bloom));
        let rc = reverse_complement(&c.seq);
        assert!(read_passes_gate(&rc, 6, &bloom));
    }

    #[test]
    fn early_reject_identical_to_full_path() {
        // 早拒与全路径严格等价：无论 raw hits 是否 ≥ GATE_MIN_HITS，
        // 早拒分支结果必须与"始终跑 SDUST + 比例门"逐 read 一致。
        let c = random_contig(60, 5);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let mut scratch = GateScratch::default();
        for seed in [11u64, 12, 13, 14, 15] {
            let read = pseudo_random_seq(48, seed);
            let fast = read_passes_gate_scratched(&read, 6, &bloom, &mut scratch);
            let mut iv = Vec::new();
            sdust_intervals_into(&read, SDUST_W, SDUST_T, &mut iv);
            let (hits, n_eff) = read_gate_hits(&read, 6, &bloom, &iv);
            let full = gate_passes(hits, n_eff);
            assert_eq!(fast, full, "seed={seed} 早拒与全路径分歧");
        }
        // 低复杂度掺杂 read（SDUST 掩蔽生效）走同一断言
        for seed in [16u64, 17] {
            let mut read = pseudo_random_seq(60, seed);
            read[10..40].copy_from_slice(&[b'A'; 30]);
            let fast = read_passes_gate_scratched(&read, 6, &bloom, &mut scratch);
            let mut iv = Vec::new();
            sdust_intervals_into(&read, SDUST_W, SDUST_T, &mut iv);
            let (hits, n_eff) = read_gate_hits(&read, 6, &bloom, &iv);
            let full = gate_passes(hits, n_eff);
            assert_eq!(fast, full, "seed={seed} 低复杂度 read 早拒与全路径分歧");
        }
    }

    #[test]
    fn random_read_fails_gate() {
        let c = random_contig(60, 2);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let rnd = pseudo_random_seq(24, 42);
        assert!(!read_passes_gate(&rnd, 6, &bloom));
    }

    #[test]
    fn short_seq_never_hits() {
        let c = contig(b"ACGTACGTACGTACGTACGTAC");
        let bloom = KmerBloom::build(&[&c], 21, GATE_FPR).unwrap();
        assert!(!read_passes_gate(b"ACGT", 21, &bloom));
    }

    #[test]
    fn bloom_superset_of_hashset_on_reads() {
        // 对照 HashSet 等价验证：Bloom 无假阴性（hits_bloom ≥ hits_set），
        // 且纯假命中增量不超过 n_eff/3。
        let c = random_contig(60, 3);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let rc = reverse_complement(&c.seq);
        let mut set: HashSet<u64> = HashSet::new();
        for i in 0..60 - 6 + 1 {
            let f = encode_kmer(&c.seq[i..i + 6]).unwrap();
            let r = encode_kmer(&rc[60 - 6 - i..60 - i]).unwrap();
            set.insert(f.min(r));
        }
        for (seed, len) in [(7u64, 30usize), (8u64, 36usize), (9u64, 42usize)] {
            let read = pseudo_random_seq(len, seed);
            let (bh, bn) = read_gate_hits(&read, 6, &bloom, &[]);
            let mut sh = 0u64;
            for i in 0..len - 6 + 1 {
                let f = encode_kmer(&read[i..i + 6]).unwrap();
                let r = {
                    let mut win = read[i..i + 6].to_vec();
                    win = reverse_complement(&win);
                    encode_kmer(&win).unwrap()
                };
                if set.contains(&f.min(r)) {
                    sh += 1;
                }
            }
            assert!(bh >= sh, "Bloom 不得漏检: {bh} < {sh}");
            assert!(
                bh - sh <= bn / 3,
                "假命中增量超界: +{} (n_eff={bn})",
                bh - sh
            );
        }
    }

    #[test]
    fn sdust_masks_low_complexity() {
        // poly-A 与双碱基周期（AC 交替）全掩；均匀随机序列近零掩蔽。
        let poly_a = vec![b'A'; 200];
        let m = sdust_mask(&poly_a, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            masked as f64 >= 0.9 * 200.0,
            "poly-A 掩蔽不足: {masked}/200"
        );

        let ac: Vec<u8> = (0..200)
            .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
            .collect();
        let m = sdust_mask(&ac, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            masked as f64 >= 0.9 * 200.0,
            "AC 交替掩蔽不足: {masked}/200"
        );

        let rnd = pseudo_random_seq(500, 42);
        let m = sdust_mask(&rnd, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            (masked as f64) <= 0.05 * 500.0,
            "随机序列误掩蔽: {masked}/500"
        );
    }

    #[test]
    fn sdust_symmetric_under_reverse_complement() {
        // SDUST 对称性（Morgulis 2006 核心性质）：掩码在反补反转下逐位对应。
        // 这是直译正确性的强检验（打乱方向/平局裁决会破坏对称）。
        let seq = pseudo_random_seq(300, 7);
        let rc = reverse_complement(&seq);
        let m1 = sdust_mask(&seq, 64, 20);
        let m2 = sdust_mask(&rc, 64, 20);
        for i in 0..seq.len() {
            assert_eq!(m1[i], m2[seq.len() - 1 - i], "对称性破坏 @{i}");
        }
    }

    #[test]
    fn masked_kmers_excluded_from_n_eff() {
        // poly-A 前缀被掩：n_eff 收缩；随机 contig 段窗口保留（高复杂度不掩），
        // 且未掩蔽窗口全部命中（Bloom 无假阴性）。hits_u 可能略高于 hits_m
        // （poly-A 窗口对词典的随机假命中），故只做单向断言。
        let c = random_contig(60, 4);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let seq: Vec<u8> = [vec![b'A'; 60], c.seq.clone()].concat();
        let iv = sdust_intervals(&seq, 64, 20);
        let (hits_m, n_m) = read_gate_hits(&seq, 6, &bloom, &iv);
        let (_, n_u) = read_gate_hits(&seq, 6, &bloom, &[]);
        assert_eq!(n_u, (seq.len() - 6 + 1) as u64);
        assert!(n_m < n_u, "掩蔽应收缩 n_eff: {n_m} vs {n_u}");
        assert!(n_m >= 50, "contig 段 55 窗口应基本保留: n_eff={n_m}");
        assert!(hits_m >= 50, "未掩蔽窗口无假阴性: hits={hits_m}");
    }
}
