//! Offline deterministic generation of stratified decoys from target genomes.
//!
//! - RNG: derive splitmix64(master_seed XOR fnv64(target_name) XOR layer XOR index), then use a
//!   xorshift64* stream. Equal seeds produce cross-platform identical output without dependencies.
//! - Mutation model: SNP-only IID substitutions at `r = 1 - ANI`, sampled from target GC background.
//!   Selecting the original base swaps within its GC class, making expected GC drift exactly zero.
//!   There are no indels, reverse complements, or artificial homologous-run caps.
//! - Guard: GC deviation above 1% emits a warning. minimap2 may independently verify ANI and
//!   within-layer exclusivity, but generation itself does not depend on alignment.
//! - Default: one ANI-85 decoy per representative. Extra layers and replicas require explicit
//!   stress-test options so default index size does not multiply unnecessarily.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use crate::reference::{self, gc_fraction};

pub const DEFAULT_ANIS: [u8; 1] = [85];
pub const DEFAULT_PER_LAYER: usize = 1;

#[derive(Debug, Clone)]
pub struct DecoyOptions {
    pub target_fa: PathBuf,
    pub out: PathBuf,
    /// ANI layer as an integer percentage, such as 82, 85, or 88.
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

/// Dependency-free FNV-1a 64 mapping a target name to a seed component.
fn fnv64(s: &[u8]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325u64;
    for &b in s {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01B3);
    }
    h
}

/// xorshift64* sampling stream; see `prescreen::splitmix64` for seed derivation.
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
    /// Uniform value in [0, 1).
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

/// Sample replacement bases from target GC background: A/T each `(1-gc)/2`, C/G each `gc/2`.
/// If the original base is selected, swap within its GC class (A/T or C/G). This gives exactly zero
/// expected GC drift for any GC fraction while ensuring every mutation changes the base.
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

/// Generate one decoy with IID SNP-only mutations; retain N and non-ACGT positions unchanged.
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

fn write_fasta_record<W: Write>(out: &mut W, name: &str, seq: &[u8]) -> std::io::Result<()> {
    out.write_all(b">")?;
    out.write_all(name.as_bytes())?;
    out.write_all(b"\n")?;
    for chunk in seq.chunks(60) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

/// Generate a decoy panel and return counts, GC warnings, and output paths.
pub fn generate(opt: &DecoyOptions) -> Result<String, String> {
    if opt.target_fa.as_os_str().is_empty() || opt.out.as_os_str().is_empty() {
        return Err("Decoy generation requires --target-fa and --out".into());
    }
    if opt.anis.is_empty() {
        return Err("--ani requires at least one layer".into());
    }
    if opt.per_layer == 0 {
        return Err("--per-layer must be greater than 0".into());
    }
    for &a in &opt.anis {
        if a == 0 || a >= 100 {
            return Err(format!(
                "Invalid --ani layer {a} (must be between 1 and 99)"
            ));
        }
    }
    let fasta = reference::parse_fasta(&opt.target_fa)?;
    if fasta.is_empty() {
        return Err("Target FASTA contains no sequences".into());
    }

    let report_path = opt.report.clone().unwrap_or_else(|| {
        let mut p = opt.out.clone().into_os_string();
        p.push(".tsv");
        PathBuf::from(p)
    });
    let fasta_file =
        File::create(&opt.out).map_err(|e| format!("Cannot create {}: {e}", opt.out.display()))?;
    let report_file = File::create(&report_path)
        .map_err(|e| format!("Cannot create report {}: {e}", report_path.display()))?;
    let mut fasta_out = BufWriter::with_capacity(1 << 20, fasta_file);
    let mut report_out = BufWriter::with_capacity(1 << 20, report_file);
    report_out
        .write_all(b"name\tsource\tani\tr\tlen\tgc_src\tgc_decoy\tseed\n")
        .map_err(|e| format!("Failed to write report {}: {e}", report_path.display()))?;

    let mut entry_count = 0usize;
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
                write_fasta_record(&mut fasta_out, &decoy_name, &decoy)
                    .map_err(|e| format!("Failed to write {}: {e}", opt.out.display()))?;
                writeln!(
                    report_out,
                    "{decoy_name}\t{name}\t{ani}\t{r:.6}\t{}\t{gc:.6}\t{decoy_gc:.6}\t{seed}",
                    decoy.len()
                )
                .map_err(|e| format!("Failed to write report {}: {e}", report_path.display()))?;
                entry_count += 1;
            }
        }
    }

    fasta_out
        .flush()
        .map_err(|e| format!("Failed to write {}: {e}", opt.out.display()))?;
    report_out
        .flush()
        .map_err(|e| format!("Failed to write report {}: {e}", report_path.display()))?;

    Ok(format!(
        "Generated {} decoys ({} targets x {} layers x {} per layer) -> {}; report {}; {} warnings for GC deviation above 1%",
        entry_count,
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

    /// Deterministic pseudo-random source sequence derived from xorshift64*.
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

    fn parse_back(path: &std::path::Path) -> Vec<(String, Vec<u8>)> {
        reference::parse_fasta(path).unwrap()
    }

    #[test]
    fn fasta_record_writer_wraps_at_sixty_bases() {
        let seq = vec![b'A'; 121];
        let mut out = Vec::new();
        write_fasta_record(&mut out, "target", &seq).unwrap();
        let lines: Vec<&[u8]> = out.split(|&b| b == b'\n').collect();
        assert_eq!(lines[0], b">target");
        assert_eq!(lines[1].len(), 60);
        assert_eq!(lines[2].len(), 60);
        assert_eq!(lines[3].len(), 1);
        assert!(lines[4].is_empty());
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
        assert_eq!(
            a, b,
            "Repeated generation with the same seed must be byte-identical"
        );
        std::fs::remove_dir_all(&dir).unwrap();
    }

    #[test]
    fn mutation_rate_no_indel_always_differs() {
        // A 100 kb periodic source with exact GC=0.5; r=0.15 should remain within 0.5%.
        let src: Vec<u8> = (0..100_000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = Rng(crate::prescreen::splitmix64(999));
        let decoy = decoy_seq(&src, 0.15, 0.5, &mut rng);
        assert_eq!(decoy.len(), src.len(), "No indels: lengths must match");
        let mut diff = 0usize;
        for (a, b) in src.iter().zip(decoy.iter()) {
            if a != b {
                diff += 1;
                assert_ne!(a.to_ascii_uppercase(), b.to_ascii_uppercase());
            }
        }
        let rate = diff as f64 / src.len() as f64;
        assert!(
            (rate - 0.15).abs() < 0.005,
            "Mutation rate deviation: {rate:.4}"
        );
    }

    #[test]
    fn gc_preserved_within_tolerance() {
        // GC=0.5 substitutions are unbiased; 100 kb sampling noise should stay below 0.8%.
        let src: Vec<u8> = (0..100_000).map(|i| b"ACGT"[i % 4]).collect();
        let mut rng = Rng(crate::prescreen::splitmix64(4242));
        let decoy = decoy_seq(&src, 0.18, 0.5, &mut rng);
        let d = (gc_fraction(&decoy) - 0.5).abs();
        assert!(d <= 0.008, "GC drift exceeds limit: {d:.4}");
    }

    #[test]
    fn gc_unbiased_for_extreme_gc() {
        // A GC=0.3 periodic source has exactly zero expected drift; verify the 1% soft limit.
        let src: Vec<u8> = (0..100_000).map(|i| b"AAAAAAACCC"[i % 10]).collect();
        let gc0 = gc_fraction(&src);
        assert!((gc0 - 0.3).abs() < 1e-9);
        let mut rng = Rng(crate::prescreen::splitmix64(777));
        let decoy = decoy_seq(&src, 0.18, gc0, &mut rng);
        let d = (gc_fraction(&decoy) - gc0).abs();
        assert!(d <= 0.008, "gc=0.3 drift exceeds limit: {d:.4}");
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
            assert!(n.starts_with("decoy:t:ani"), "Naming convention: {n}");
            assert!(names.insert(n), "Duplicate name: {n}");
        }
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
