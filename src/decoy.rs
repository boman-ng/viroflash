//! `viroflash decoy`：依据目标基因组离线、确定性地生成分层诱饵。
//!
//! - RNG：splitmix64(master_seed ⊕ fnv64(target_name) ⊕ layer ⊕ idx) 派生 →
//!   xorshift64* 采样流（Vigna 2016；同种子同输出、跨平台确定、零依赖）。
//! - 突变模型：SNP-only；每诱饵按率 r = 1−ani 逐位 iid；替换采样 = 目标 GC
//!   背景分布 + 抽中原碱基时同 GC 组交换（E[ΔGC]=0 精确、任意 gc）；**无 indel**
//!   （层带方差最小、覆盖统计稳定）；**无 RC**；**无同源 run 硬限**（iid 突变下
//!   150kb 最长同源 run 期望 ≈51/61/76bp 与真实近缘同分布，设上限反使诱饵偏离真实行为）。
//! - 门禁：GC 偏差 ≤1% 生成时软检查（超限记录警告，不静默）；层内互斥与
//!   ANI 精度可用 minimap2 验证（生成器不依赖比对）。
//! - k-mer 共享率：与目标精确 21-mer 共享率 = p²¹ = 1.5%/3.3%/6.8%
//!   （82/85/88 层）。

use std::path::{Path, PathBuf};

use crate::reference::{self, gc_fraction};

pub const DEFAULT_ANIS: [u8; 3] = [82, 85, 88];
pub const DEFAULT_PER_LAYER: usize = 4;

#[derive(Debug, Clone)]
pub struct DecoyOptions {
    pub target_fa: PathBuf,
    pub out: PathBuf,
    /// ANI 层（百分比整数，如 82/85/88）。
    pub anis: Vec<u8>,
    pub per_layer: usize,
    pub seed: u64,
    pub report: Option<PathBuf>,
}

impl Default for DecoyOptions {
    fn default() -> Self {
        Self {
            target_fa: PathBuf::new(),
            out: PathBuf::new(),
            anis: DEFAULT_ANIS.to_vec(),
            per_layer: DEFAULT_PER_LAYER,
            seed: 0,
            report: None,
        }
    }
}

/// FNV-1a 64（目标名 → 种子分量；零依赖）。
fn fnv64(s: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in s {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// xorshift64* 采样流（Vigna 2016；splitmix64 种子派生见 prescreen::splitmix64）。
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut s = self.0;
        s ^= s >> 12;
        s ^= s << 25;
        s ^= s >> 27;
        self.0 = s;
        s.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }
    /// [0,1) 均匀。
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

/// 替换采样：新碱基按目标 GC 含量的背景分布采样（A/T 各 (1−gc)/2，
/// C/G 各 gc/2）；若抽中原碱基则改为同 GC 组另一碱基（A↔T、C↔G）。
/// 该构造下 E[ΔGC]=0 **精确**（任意 gc，无系统性漂移），且每个突变事件
/// 必改变碱基（实际突变率 = r，ANI 层位无偏）。
fn sample_replacement(rng: &mut Rng, cur: u8, gc: f64) -> u8 {
    let pa = (1.0 - gc) / 2.0;
    let pc = gc / 2.0;
    let u = rng.unit();
    let b = if u < pa {
        b'A'
    } else if u < 0.5 {
        b'C'
    } else if u < 0.5 + pc {
        b'G'
    } else {
        b'T'
    };
    if b == cur.to_ascii_uppercase() {
        return match cur.to_ascii_uppercase() {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            _ => b'C',
        };
    }
    b
}

/// 单条诱饵：SNP-only iid 突变。N/非 ACGT 位保留不突变。
fn decoy_seq(src: &[u8], r: f64, gc: f64, rng: &mut Rng) -> Vec<u8> {
    src.iter()
        .map(|&b| {
            if crate::prescreen::dna_bits(b).is_none() {
                b
            } else if rng.unit() < r {
                sample_replacement(rng, b, gc)
            } else {
                b
            }
        })
        .collect()
}

fn write_fasta(path: &Path, entries: &[(String, Vec<u8>)]) -> Result<(), String> {
    let mut out = String::new();
    for (name, seq) in entries {
        out.push('>');
        out.push_str(name);
        out.push('\n');
        for chunk in seq.chunks(60) {
            out.push_str(std::str::from_utf8(chunk).map_err(|e| e.to_string())?);
            out.push('\n');
        }
    }
    std::fs::write(path, out).map_err(|e| format!("写 {} 失败: {e}", path.display()))
}

/// 生成诱饵面板。返回摘要（条数、GC 警告数、输出路径）。
pub fn generate(opt: &DecoyOptions) -> Result<String, String> {
    if opt.target_fa.as_os_str().is_empty() || opt.out.as_os_str().is_empty() {
        return Err("诱饵生成需要 --target-fa 与 --out".into());
    }
    if opt.anis.is_empty() {
        return Err("--ani 至少一层".into());
    }
    if opt.per_layer == 0 {
        return Err("--per-layer 必须大于 0".into());
    }
    for &a in &opt.anis {
        if a == 0 || a >= 100 {
            return Err(format!("--ani 层 {a} 非法（须在 1..99 之间）"));
        }
    }
    let fasta = reference::parse_fasta(&opt.target_fa)?;
    if fasta.is_empty() {
        return Err("目标 FASTA 中没有序列".into());
    }

    let mut entries: Vec<(String, Vec<u8>)> = Vec::new();
    let mut report_rows: Vec<String> =
        vec!["name\tsource\tani\tr\tlen\tgc_src\tgc_decoy\tseed".to_string()];
    let mut gc_warnings = 0usize;

    for (name, seq) in &fasta {
        let gc = gc_fraction(seq);
        let name_seed = fnv64(name.as_bytes());
        for &ani in &opt.anis {
            let r = 1.0 - f64::from(ani) / 100.0;
            for idx in 0..opt.per_layer {
                let seed = crate::prescreen::splitmix64(
                    opt.seed ^ name_seed ^ (u64::from(ani) << 48) ^ idx as u64,
                );
                let mut rng = Rng(seed);
                let decoy = decoy_seq(seq, r, gc, &mut rng);
                let decoy_gc = gc_fraction(&decoy);
                if (decoy_gc - gc).abs() > 0.01 {
                    gc_warnings += 1;
                }
                let decoy_name = format!("decoy:{name}:ani{ani}:i{idx}");
                report_rows.push(format!(
                    "{decoy_name}\t{name}\t{ani}\t{r:.6}\t{}\t{gc:.6}\t{decoy_gc:.6}\t{seed}",
                    decoy.len()
                ));
                entries.push((decoy_name, decoy));
            }
        }
    }

    write_fasta(&opt.out, &entries)?;
    let report_path = opt.report.clone().unwrap_or_else(|| {
        let mut p = opt.out.clone().into_os_string();
        p.push(".tsv");
        PathBuf::from(p)
    });
    std::fs::write(&report_path, report_rows.join("\n") + "\n")
        .map_err(|e| format!("写报告 {} 失败: {e}", report_path.display()))?;

    Ok(format!(
        "诱饵 {} 条（目标 {} 条 × 层 {} × 每层 {}）→ {}；报告 {}；GC 偏差超 1% 警告 {} 条",
        entries.len(),
        fasta.len(),
        opt.anis.len(),
        opt.per_layer,
        opt.out.display(),
        report_path.display(),
        gc_warnings
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    /// 确定性伪随机源序列（xorshift64* 派生）。
    fn pseudo_seq(len: usize, seed: u64) -> Vec<u8> {
        let mut rng = Rng(seed | 1);
        (0..len)
            .map(|_| b"ACGT"[(rng.next() >> 59) as usize & 3])
            .collect()
    }

    fn tmp_dir(tag: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("viroflash_decoy_test_{tag}_{}", std::process::id()));
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn parse_back(path: &Path) -> Vec<(String, Vec<u8>)> {
        reference::parse_fasta(path).unwrap()
    }

    #[test]
    fn same_seed_reproducible_byte_identical() {
        let dir = tmp_dir("repro");
        let fa = dir.join("t.fa");
        let src = pseudo_seq(50_000, 7);
        std::fs::write(&fa, format!(">t\n{}", String::from_utf8(src).unwrap())).unwrap();
        let out1 = dir.join("d1.fa");
        let out2 = dir.join("d2.fa");
        for (out, tag) in [(&out1, "r1"), (&out2, "r2")] {
            let opt = DecoyOptions {
                target_fa: fa.clone(),
                out: out.clone(),
                anis: vec![82, 88],
                per_layer: 2,
                seed: 12345,
                report: Some(dir.join(format!("{tag}.tsv"))),
            };
            generate(&opt).unwrap();
        }
        let a = std::fs::read(&out1).unwrap();
        let b = std::fs::read(&out2).unwrap();
        assert_eq!(a, b, "同 seed 两次生成必须逐字节一致");
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn mutation_rate_no_indel_always_differs() {
        // 100kb 周期源（GC 精确 0.5），r=0.15：diff 率 ∈ 0.15±0.5%（σ≈0.11%）。
        let src: Vec<u8> = (0..100_000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = Rng(crate::prescreen::splitmix64(999));
        let decoy = decoy_seq(&src, 0.15, 0.5, &mut rng);
        assert_eq!(decoy.len(), src.len(), "无 indel：长度必须相等");
        let mut diff = 0usize;
        for (a, b) in src.iter().zip(decoy.iter()) {
            if a != b {
                diff += 1;
                assert_ne!(a.to_ascii_uppercase(), b.to_ascii_uppercase());
            }
        }
        let rate = diff as f64 / src.len() as f64;
        assert!((rate - 0.15).abs() < 0.005, "突变率偏离: {rate:.4}");
    }

    #[test]
    fn gc_preserved_within_tolerance() {
        // gc=0.5 源替换无偏（E[ΔGC]=0），100kb 采样噪声 σ≈0.06% → 断言 ≤0.8%。
        let src: Vec<u8> = (0..100_000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = Rng(crate::prescreen::splitmix64(4242));
        let decoy = decoy_seq(&src, 0.18, 0.5, &mut rng);
        let d = (gc_fraction(&decoy) - 0.5).abs();
        assert!(d <= 0.008, "GC 漂移超界: {d:.4}");
    }

    #[test]
    fn gc_unbiased_for_extreme_gc() {
        // gc=0.3 的周期源（AAAAAAACCC × 10000）：E[ΔGC]=0 精确，
        // 验证 GC 偏差不超过 1% 软检查线。
        let src: Vec<u8> = (0..100_000).map(|i| b"AAAAAAACCC"[i % 10]).collect();
        let gc0 = gc_fraction(&src);
        assert!((gc0 - 0.3).abs() < 1e-9);
        let mut rng = Rng(crate::prescreen::splitmix64(777));
        let decoy = decoy_seq(&src, 0.18, gc0, &mut rng);
        let d = (gc_fraction(&decoy) - gc0).abs();
        assert!(d <= 0.008, "gc=0.3 漂移超界: {d:.4}");
    }

    #[test]
    fn naming_unique_and_layered() {
        let src = pseudo_seq(3_000, 3);
        let dir = tmp_dir("names");
        let fa = dir.join("t.fa");
        std::fs::write(&fa, format!(">t\n{}", String::from_utf8(src).unwrap())).unwrap();
        let opt = DecoyOptions {
            target_fa: fa,
            out: dir.join("d.fa"),
            anis: vec![82, 85, 88],
            per_layer: 4,
            seed: 0,
            report: None,
        };
        generate(&opt).unwrap();
        let entries = parse_back(&dir.join("d.fa"));
        assert_eq!(entries.len(), 12);
        let mut names: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for (n, _) in &entries {
            assert!(n.starts_with("decoy:t:ani"), "命名约定: {n}");
            assert!(names.insert(n), "名字重复: {n}");
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
