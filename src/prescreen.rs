//! Canonical k-mer Bloom prescreening with SDUST low-complexity masking.
//!
//! The Bloom filter is sized from observed occupancy to enforce the configured false-positive
//! bound without introducing false negatives. Rolling forward/reverse-complement encoders produce
//! canonical k-mers, while SDUST masking removes low-complexity windows from both the hit count and
//! effective denominator. Implementations preserve deterministic behavior across platforms.

use crate::reference::Contig;

pub const DEFAULT_K: usize = 21;

pub const K_MAX: usize = 31;

pub const GATE_FPR: f64 = 0.1;

pub const GATE_MIN_HITS: u64 = 10;

pub const GATE_HIT_FRAC_NUM: u64 = 7;
pub const GATE_HIT_FRAC_DEN: u64 = 20;

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

#[inline]
fn for_each_canonical_kmer<F>(seq: &[u8], k: usize, mut visit: F)
where
    F: FnMut(usize, u64),
{
    if k == 0 || k > K_MAX || seq.len() < k {
        return;
    }
    let kmask = (1u64 << (2 * k)) - 1;
    let mut fwd = 0u64;
    let mut rev = 0u64;
    let mut filled = 0usize;
    for (i, &c) in seq.iter().enumerate() {
        match dna_bits(c) {
            Some(b) => {
                fwd = ((fwd << 2) | u64::from(b)) & kmask;
                rev = (rev >> 2) | (u64::from(3 - b) << (2 * k - 2));
                filled = (filled + 1).min(k);
                if filled == k {
                    visit(i + 1 - k, fwd.min(rev));
                }
            }
            None => {
                fwd = 0;
                rev = 0;
                filled = 0;
            }
        }
    }
}

pub(crate) fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

#[derive(Debug)]
pub struct KmerBloom {
    words: Vec<u64>,
    mask: u64,
    pub k: usize,

    pub n_inserted: u64,
}

impl KmerBloom {
    pub fn build(contigs: &[&Contig], k: usize, fpr: f64) -> Option<Self> {
        let full_capacity = Self::capacity_bits(contigs, k, fpr)?;
        if full_capacity > 64 {
            let mut candidate = Self::build_with_capacity(contigs, k, full_capacity / 2);
            if candidate.fill_frac() <= fpr {
                candidate.fold_to_fpr(fpr);
                return Some(candidate);
            }
        }
        let mut bloom = Self::build_with_capacity(contigs, k, full_capacity);
        bloom.fold_to_fpr(fpr);
        Some(bloom)
    }

    #[cfg(test)]
    fn build_unfolded(contigs: &[&Contig], k: usize, fpr: f64) -> Option<Self> {
        let capacity = Self::capacity_bits(contigs, k, fpr)?;
        Some(Self::build_with_capacity(contigs, k, capacity))
    }

    fn capacity_bits(contigs: &[&Contig], k: usize, fpr: f64) -> Option<u64> {
        if k == 0 || k > K_MAX {
            return None;
        }
        let bases: u64 = contigs.iter().map(|c| c.seq.len() as u64).sum();
        if bases == 0 {
            return None;
        }
        let m_min_bits = (bases as f64 * (-1.0 / (1.0 - fpr).ln())).ceil().max(1.0) as u64;
        Some(m_min_bits.max(64).next_power_of_two())
    }

    fn build_with_capacity(contigs: &[&Contig], k: usize, capacity: u64) -> Self {
        let words = vec![0u64; (capacity as usize) / 64];
        let mut bloom = Self {
            words,
            mask: capacity - 1,
            k,
            n_inserted: 0,
        };
        for c in contigs {
            for_each_canonical_kmer(&c.seq, k, |_, canon| {
                let idx = (splitmix64(canon) & bloom.mask) as usize;
                bloom.words[idx >> 6] |= 1 << (idx & 63);
                bloom.n_inserted += 1;
            });
        }
        bloom
    }

    fn fold_to_fpr(&mut self, fpr: f64) {
        while self.words.len() > 1 {
            let half = self.words.len() / 2;
            let folded_set_bits: u64 = (0..half)
                .map(|i| u64::from((self.words[i] | self.words[i + half]).count_ones()))
                .sum();
            let folded_fill = folded_set_bits as f64 / (half * 64) as f64;
            if folded_fill > fpr {
                break;
            }
            for i in 0..half {
                self.words[i] |= self.words[i + half];
            }
            self.words.truncate(half);
            self.mask = (half * 64) as u64 - 1;
        }
        self.words.shrink_to_fit();
    }

    #[inline]
    pub fn probe(&self, code: u64) -> bool {
        let idx = (splitmix64(code) & self.mask) as usize;
        (self.words[idx >> 6] >> (idx & 63)) & 1 == 1
    }

    pub fn fill_frac(&self) -> f64 {
        let set: u64 = self.words.iter().map(|w| u64::from(w.count_ones())).sum();
        set as f64 / (self.words.len() * 64) as f64
    }

    pub(crate) fn to_raw(&self) -> (&[u64], u64) {
        (&self.words, self.mask)
    }

    pub(crate) fn from_raw(words: Vec<u64>, mask: u64, k: usize, n_inserted: u64) -> Self {
        Self {
            words,
            mask,
            k,
            n_inserted,
        }
    }
}

// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

struct PerfectIntv {
    start: usize,
    finish: usize,
    r: i64,
    l: i64,
}

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

    #[inline]
    fn at(&self, i: usize) -> u32 {
        self.buf[(self.head + i) & RING_MASK]
    }
}

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
        let s = (que.pop_front() as usize) & 63;
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
    let s = t3 as usize;
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

pub fn sdust_intervals_into(seq: &[u8], w: usize, t: i64, out: &mut Vec<(usize, usize)>) {
    out.clear();
    let mut p: Vec<PerfectIntv> = Vec::new();
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

pub fn sdust_intervals(seq: &[u8], w: usize, t: i64) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    sdust_intervals_into(seq, w, t, &mut out);
    out
}

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

// ---------------------------------------------------------------------------

pub fn gate_passes(hits: u64, n_eff: u64) -> bool {
    if n_eff == 0 {
        return false;
    }
    let frac = (GATE_HIT_FRAC_NUM * n_eff).div_ceil(GATE_HIT_FRAC_DEN);
    hits >= GATE_MIN_HITS.max(frac)
}

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
    let mut it = 0usize;
    for_each_canonical_kmer(seq, k, |start, canon| {
        while it < intervals.len() && intervals[it].1 <= start {
            it += 1;
        }
        let masked = it < intervals.len() && intervals[it].0 <= start;
        if !masked {
            n_eff += 1;
            if bloom.probe(canon) {
                hits += 1;
            }
        }
    });
    (hits, n_eff)
}

#[derive(Default)]
pub struct GateScratch {
    pub intervals: Vec<(usize, usize)>,
}

pub fn read_passes_gate_scratched(
    seq: &[u8],
    k: usize,
    bloom: &KmerBloom,
    scratch: &mut GateScratch,
) -> bool {
    let (raw_hits, _) = read_gate_hits(seq, k, bloom, &[]);
    if raw_hits < GATE_MIN_HITS {
        return false;
    }
    sdust_intervals_into(seq, SDUST_W, SDUST_T, &mut scratch.intervals);
    let (hits, n_eff) = read_gate_hits(seq, k, bloom, &scratch.intervals);
    gate_passes(hits, n_eff)
}

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

    fn random_contig(len: usize, seed: u64) -> Contig {
        contig(&pseudo_random_seq(len, seed))
    }

    fn build_slice_oracle(contigs: &[&Contig], k: usize, fpr: f64) -> Option<KmerBloom> {
        if k == 0 || k > K_MAX {
            return None;
        }
        let bases: u64 = contigs.iter().map(|c| c.seq.len() as u64).sum();
        if bases == 0 {
            return None;
        }
        let m_min_bits = (bases as f64 * (-1.0 / (1.0 - fpr).ln())).ceil().max(1.0) as u64;
        let m = m_min_bits.max(64).next_power_of_two();
        let mut bloom = KmerBloom {
            words: vec![0u64; (m as usize) / 64],
            mask: m - 1,
            k,
            n_inserted: 0,
        };
        for c in contigs {
            let rc = reverse_complement(&c.seq);
            let n = c.seq.len();
            for i in 0..n.saturating_sub(k - 1) {
                let (Some(fwd), Some(rev)) = (
                    encode_kmer(&c.seq[i..i + k]),
                    encode_kmer(&rc[n - k - i..n - i]),
                ) else {
                    continue;
                };
                let idx = (splitmix64(fwd.min(rev)) & bloom.mask) as usize;
                bloom.words[idx >> 6] |= 1 << (idx & 63);
                bloom.n_inserted += 1;
            }
        }
        Some(bloom)
    }

    #[test]
    fn kmer_encode_and_rc() {
        assert_eq!(dna_bits(b'a'), Some(0));
        assert_eq!(dna_bits(b'N'), None);
        let fwd = encode_kmer(b"ACGT").unwrap();
        let rc = encode_kmer(&reverse_complement(b"ACGT")).unwrap();
        assert_eq!(fwd, rc);
        let a = encode_kmer(b"AAAA").unwrap();
        let t = encode_kmer(b"TTTT").unwrap();
        assert_ne!(a, t);
    }

    #[test]
    fn rolling_bloom_matches_slice_oracle_bit_for_bit() {
        let contigs = [
            contig(b"ACGTNACGTTTTTACGCGTATACGATCGATCGATCGATC"),
            random_contig(96, 0xDEAD_BEEF),
            contig(b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCA"),
        ];
        let refs: Vec<&Contig> = contigs.iter().collect();
        for k in [1, 6, 21, 31] {
            let rolling = KmerBloom::build_unfolded(&refs, k, GATE_FPR).unwrap();
            let sliced = build_slice_oracle(&refs, k, GATE_FPR).unwrap();
            assert_eq!(rolling.k, sliced.k, "k={k}");
            assert_eq!(rolling.mask, sliced.mask, "k={k}");
            assert_eq!(rolling.n_inserted, sliced.n_inserted, "k={k}");
            assert_eq!(rolling.words, sliced.words, "k={k}");
        }
    }

    #[test]
    fn bloom_folding_reduces_capacity_without_false_negatives() {
        let seq = vec![b'A'; 4096];
        let c = contig(&seq);
        let mut unfolded = KmerBloom::build_unfolded(&[&c], 21, GATE_FPR).unwrap();
        let folded = KmerBloom::build(&[&c], 21, GATE_FPR).unwrap();
        assert!(folded.words.len() < unfolded.words.len());
        assert!(folded.fill_frac() <= GATE_FPR);

        for (word_index, &word) in unfolded.words.iter().enumerate() {
            let mut set = word;
            while set != 0 {
                let bit = set.trailing_zeros() as usize;
                let old_index = word_index * 64 + bit;
                let new_index = old_index & folded.mask as usize;
                assert_ne!(folded.words[new_index >> 6] & (1 << (new_index & 63)), 0);
                set &= set - 1;
            }
        }

        let code = encode_kmer(&c.seq[..21]).unwrap();
        assert!(folded.probe(code));

        unfolded.fold_to_fpr(GATE_FPR);
        assert_eq!(folded.mask, unfolded.mask);
        assert_eq!(folded.n_inserted, unfolded.n_inserted);
        assert_eq!(folded.words, unfolded.words);
    }

    #[test]
    fn bloom_half_capacity_falls_back_when_fill_would_exceed_bound() {
        let c = random_contig(4096, 0x0A11_CE55);
        let refs = [&c];
        let mut unfolded = KmerBloom::build_unfolded(&refs, 21, GATE_FPR).unwrap();
        let half = KmerBloom::build_with_capacity(&refs, 21, unfolded.mask.div_ceil(2));
        assert!(
            half.fill_frac() > GATE_FPR,
            "fixture must trigger conservative capacity fallback"
        );

        unfolded.fold_to_fpr(GATE_FPR);
        let built = KmerBloom::build(&refs, 21, GATE_FPR).unwrap();
        assert_eq!(built.mask, unfolded.mask);
        assert_eq!(built.n_inserted, unfolded.n_inserted);
        assert_eq!(built.words, unfolded.words);
    }

    #[test]
    fn bloom_all_contig_kmers_probe_true() {
        let c = random_contig(60, 1);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let rc = reverse_complement(&c.seq);
        for i in 0..60 - 6 + 1 {
            let f = encode_kmer(&c.seq[i..i + 6]).unwrap();
            let r = encode_kmer(&rc[60 - 6 - i..60 - i]).unwrap();
            assert!(bloom.probe(f.min(r)), "k-mer @{i} should be present");
        }
    }

    #[test]
    fn bloom_fpr_within_bound() {
        let c = contig(b"ACGTACGTACGTACGTACGTAC");
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let ff = bloom.fill_frac();
        assert!(ff > 0.0 && ff < 1.0, "invalid fill fraction: {ff}");
        let n = 50_000;
        let mut fp = 0u64;
        let mut s = 99u64;
        for _ in 0..n {
            s = xorshift64_star(s);
            let code = (s >> 34) & ((1u64 << 12) - 1);
            if bloom.probe(code) {
                fp += 1;
            }
        }
        assert!(
            (fp as f64) < 0.2 * n as f64,
            "false-positive rate exceeds limit: {fp}/{n}"
        );
    }

    #[test]
    fn gate_passes_thresholds() {
        assert!(!gate_passes(34, 100));
        assert!(gate_passes(35, 100));
        // max(10, ⌈0.35·50⌉=18) = 18
        assert!(!gate_passes(9, 50));
        assert!(gate_passes(18, 50));
        // max(10, ⌈0.35·10⌉=4) = 10
        assert!(!gate_passes(9, 10));
        assert!(gate_passes(10, 10));

        assert!(!gate_passes(0, 0));
        assert!(!gate_passes(5, 3));
    }

    #[test]
    fn read_passes_gate_fwd_and_rc_reads() {
        let c = random_contig(60, 1);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        assert!(read_passes_gate(&c.seq, 6, &bloom));
        let rc = reverse_complement(&c.seq);
        assert!(read_passes_gate(&rc, 6, &bloom));
    }

    #[test]
    fn early_reject_identical_to_full_path() {
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
            assert_eq!(
                fast, full,
                "seed={seed}: early rejection differs from full path"
            );
        }

        for seed in [16u64, 17] {
            let mut read = pseudo_random_seq(60, seed);
            read[10..40].copy_from_slice(&[b'A'; 30]);
            let fast = read_passes_gate_scratched(&read, 6, &bloom, &mut scratch);
            let mut iv = Vec::new();
            sdust_intervals_into(&read, SDUST_W, SDUST_T, &mut iv);
            let (hits, n_eff) = read_gate_hits(&read, 6, &bloom, &iv);
            let full = gate_passes(hits, n_eff);
            assert_eq!(
                fast, full,
                "seed={seed}: low-complexity early rejection differs from full path"
            );
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
            assert!(
                bh >= sh,
                "Bloom filter must not miss true hits: {bh} < {sh}"
            );
            assert!(
                bh - sh <= bn / 3,
                "false-hit increase exceeds limit: +{} (n_eff={bn})",
                bh - sh
            );
        }
    }

    #[test]
    fn sdust_masks_low_complexity() {
        let poly_a = vec![b'A'; 200];
        let m = sdust_mask(&poly_a, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            masked as f64 >= 0.9 * 200.0,
            "insufficient poly-A masking: {masked}/200"
        );

        let ac: Vec<u8> = (0..200)
            .map(|i| if i % 2 == 0 { b'A' } else { b'C' })
            .collect();
        let m = sdust_mask(&ac, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            masked as f64 >= 0.9 * 200.0,
            "insufficient alternating-AC masking: {masked}/200"
        );

        let rnd = pseudo_random_seq(500, 42);
        let m = sdust_mask(&rnd, 64, 20);
        let masked = m.iter().filter(|&&x| x).count();
        assert!(
            (masked as f64) <= 0.05 * 500.0,
            "random sequence was over-masked: {masked}/500"
        );
    }

    #[test]
    fn sdust_symmetric_under_reverse_complement() {
        let seq = pseudo_random_seq(300, 7);
        let rc = reverse_complement(&seq);
        let m1 = sdust_mask(&seq, 64, 20);
        let m2 = sdust_mask(&rc, 64, 20);
        for i in 0..seq.len() {
            assert_eq!(
                m1[i],
                m2[seq.len() - 1 - i],
                "reverse-complement symmetry broken at {i}"
            );
        }
    }

    #[test]
    fn masked_kmers_excluded_from_n_eff() {
        let c = random_contig(60, 4);
        let bloom = KmerBloom::build(&[&c], 6, GATE_FPR).unwrap();
        let seq: Vec<u8> = [vec![b'A'; 60], c.seq.clone()].concat();
        let iv = sdust_intervals(&seq, 64, 20);
        let (hits_m, n_m) = read_gate_hits(&seq, 6, &bloom, &iv);
        let (_, n_u) = read_gate_hits(&seq, 6, &bloom, &[]);
        assert_eq!(n_u, (seq.len() - 6 + 1) as u64);
        assert!(n_m < n_u, "masking should reduce n_eff: {n_m} vs {n_u}");
        assert!(
            n_m >= 50,
            "the 55-window contig segment should remain mostly available: n_eff={n_m}"
        );
        assert!(
            hits_m >= 50,
            "unmasked windows must have no false negatives: hits={hits_m}"
        );
    }
}
