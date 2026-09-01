//! Command-line entry point for `viroflash index`, `viroflash run`, and `viroflash version`.
//! Arguments are parsed directly; the library's `sysinfo` module provides telemetry.

use std::path::PathBuf;
use std::process::ExitCode;

use viroflash::index::{self, IndexOptions};
use viroflash::reference::Role;
use viroflash::{run_pipeline, Options};

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("Error: {e}");
            ExitCode::FAILURE
        }
    }
}

#[derive(Debug)]
enum Command {
    Help,
    Version,
    Index(IndexOptions),
    Run(Options),
}

fn run() -> Result<(), String> {
    match parse_args()? {
        Command::Help => {
            print_usage();
            Ok(())
        }
        Command::Version => {
            println!("viroflash {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Command::Index(opt) => {
            let built = index::build_index(&opt)?;
            let n = |r| built.contigs.iter().filter(|c| c.role == r).count();
            eprintln!(
                "Index built: {} (contigs: host={} target={} decoy={} contam={}; k={})",
                opt.out_dir.display(),
                n(Role::Host),
                n(Role::Target),
                n(Role::Decoy),
                n(Role::Contaminant),
                opt.k
            );
            eprintln!(
                "  ref.mmi / bloom.bin / manifest.json / targets.fa (manifest_blake3={})",
                built.manifest_blake3
            );
            let perf = viroflash::perf::report_paths(&opt.out_dir);
            eprintln!(
                "Performance reports: {} / {}",
                perf.json.display(),
                perf.tsv.display()
            );
            Ok(())
        }
        Command::Run(opt) => {
            let summary = run_pipeline(&opt)?;
            eprintln!("Input pairs: {}", summary.input_pairs);
            eprintln!(
                "Sampled pairs: {} (p={:.6}); prescreened pairs in sample: {}",
                summary.sampling.selected_pairs,
                summary.sampling.inclusion_probability,
                summary.prescreen_pairs
            );
            if summary.map_errors > 0 {
                eprintln!(
                    "Warning: {} fragments failed alignment (see run.map_errors)",
                    summary.map_errors
                );
            }
            if summary.sampling.audit_overflows > 0 {
                eprintln!(
                    "Warning: {} read ends exceeded the whole-chain hit limit and were not treated as empty evidence (see sampling.audit_overflows)",
                    summary.sampling.audit_overflows
                );
            }
            eprintln!(
                "Fixed testing family: {}; reported candidates: {}",
                summary.test_family_size,
                summary.candidates.len()
            );
            eprintln!("Interpretation scope: candidate-level research gate; no sample-level conclusion; QC=NOT_EVALUATED");
            for c in &summary.candidates {
                eprintln!(
                    "  hypothesis={} OR({})  representative={}  adjusted-p={}  status={}  candidate_gate={}  breadth={:.2}  reads={}  integration={}  split={}  discordant={}",
                    c.contig,
                    c.hypothesis_members.len(),
                    c.representative,
                    if c.q_underflow {
                        format!(
                            "{:.3e} (underflow-range, log10={:.3})",
                            c.q_value,
                            c.ln_q_value / std::f64::consts::LN_10
                        )
                    } else {
                        format!("{:.3e}", c.q_value)
                    },
                    c.confidence(),
                    c.decision,
                    c.covered_frac,
                    c.reads,
                    c.integration_evidence,
                    c.split_events,
                    c.discordant
                );
            }
            eprintln!(
                "Reports: {} / {} / {} / {}",
                summary.result_json.display(),
                summary.result_tsv.display(),
                summary.result_html.display(),
                summary.result_csv.display()
            );
            let perf = viroflash::perf::report_paths(&opt.out);
            eprintln!(
                "Performance reports: {} / {}",
                perf.json.display(),
                perf.tsv.display()
            );
            Ok(())
        }
    }
}

fn parse_args() -> Result<Command, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    parse_args_from(&args)
}

fn parse_args_from(args: &[String]) -> Result<Command, String> {
    if args.is_empty() || args[0] == "-h" || args[0] == "--help" {
        return Ok(Command::Help);
    }
    match args[0].as_str() {
        "version" | "--version" => Ok(Command::Version),
        "index" => parse_index(&args[1..]),
        "run" => parse_run(&args[1..]),
        other => Err(format!("Unknown subcommand: {other}")),
    }
}

fn parse_num<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("Invalid value for {flag}: {value}"))
}

fn parse_anis(value: &str) -> Result<Vec<u8>, String> {
    value
        .split(',')
        .map(|s| {
            s.trim()
                .parse()
                .map_err(|_| format!("Invalid value for --decoy-ani: {s}"))
        })
        .collect()
}

/// Build a reusable index; generate decoys when `--decoy-fa` is omitted.
fn parse_index(args: &[String]) -> Result<Command, String> {
    let mut opt = IndexOptions::default();
    let mut decoy_tuning = false;
    let mut i = 0;
    while i < args.len() {
        let name = args[i].as_str();
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("Missing value for argument {name}"))?;
        match name {
            "--host-fa" => opt.host_fa = PathBuf::from(value),
            "--target-fa" => opt.target_fa = PathBuf::from(value),
            "--contam-fa" => opt.contam_fa = Some(PathBuf::from(value)),
            "--decoy-fa" => opt.decoy_fa = Some(PathBuf::from(value)),
            "--decoy-ani" => {
                opt.decoy_anis = parse_anis(value)?;
                decoy_tuning = true;
            }
            "--decoy-per-layer" => {
                opt.decoy_per_layer = parse_num("--decoy-per-layer", value)?;
                decoy_tuning = true;
            }
            "--decoy-seed" => {
                opt.decoy_seed = parse_num("--decoy-seed", value)?;
                decoy_tuning = true;
            }
            "--out" => opt.out_dir = PathBuf::from(value),
            "--k" => opt.k = parse_num("--k", value)?,
            "--threads" => opt.threads = parse_num("--threads", value)?,
            other => return Err(format!("Unknown argument: {other}")),
        }
        i += 2;
    }
    if opt.decoy_fa.is_some() && decoy_tuning {
        return Err(
            "--decoy-fa cannot be combined with --decoy-ani/--decoy-per-layer/--decoy-seed (generation options apply only to generated decoys)"
                .into(),
        );
    }
    if opt.host_fa.as_os_str().is_empty() || opt.target_fa.as_os_str().is_empty() {
        return Err("The index subcommand requires --host-fa and --target-fa".into());
    }
    if opt.out_dir.as_os_str().is_empty() {
        return Err("The index subcommand requires --out".into());
    }
    Ok(Command::Index(opt))
}

/// Run with a reusable `--index` or build one from FASTA inputs (mutually exclusive).
fn parse_run(args: &[String]) -> Result<Command, String> {
    let mut opt = Options::default();
    let mut i = 0;
    while i < args.len() {
        let name = args[i].as_str();
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("Missing value for argument {name}"))?;
        match name {
            "--r1" => opt.r1 = PathBuf::from(value),
            "--r2" => opt.r2 = Some(PathBuf::from(value)),
            "--index" => opt.index = Some(PathBuf::from(value)),
            "--host-fa" => opt.host_fa = Some(PathBuf::from(value)),
            "--target-fa" => opt.target_fa = Some(PathBuf::from(value)),
            "--decoy-fa" => opt.decoy_fa = Some(PathBuf::from(value)),
            "--contam-fa" => opt.contam_fa = Some(PathBuf::from(value)),
            "--threads" => opt.threads = parse_num("--threads", value)?,
            "--out" => opt.out = PathBuf::from(value),
            "--k" => opt.k = parse_num("--k", value)?,
            other => return Err(format!("Unknown argument: {other}")),
        }
        i += 2;
    }
    if opt.r1.as_os_str().is_empty() {
        return Err("The run subcommand requires --r1".into());
    }
    if opt.index.is_some() {
        for (label, p) in [
            ("--host-fa", &opt.host_fa),
            ("--target-fa", &opt.target_fa),
            ("--decoy-fa", &opt.decoy_fa),
            ("--contam-fa", &opt.contam_fa),
        ] {
            if p.is_some() {
                return Err(format!(
                    "--index cannot be combined with {label} (a loaded index obtains reference metadata from its manifest)"
                ));
            }
        }
    } else if opt.host_fa.is_none() || opt.target_fa.is_none() {
        return Err(
            "The run subcommand requires --host-fa and --target-fa when --index is omitted".into(),
        );
    }
    Ok(Command::Run(opt))
}

fn print_usage() {
    println!(
        "viroflash {} — viral candidate detection command-line tool\n\
\n\
Usage:\n\
  viroflash index \\\n\
    --host-fa <host.fa> --target-fa <target-virus.fa> \\\n\
    [--contam-fa <contaminant.fa>] [--decoy-fa <decoy.fa>] \\\n\
    [--decoy-ani 85] [--decoy-per-layer 1] [--decoy-seed 0] \\\n\
    --out <index-directory> [--k 21 (1..=31)] [--threads 8]\n\
\n\
  viroflash run \\\n\
    --r1 <reads_R1.fastq.gz> [--r2 <reads_R2.fastq.gz>] \\\n\
    (--index <index-directory> | --host-fa <host.fa> --target-fa <target-virus.fa> \\\n\
      [--contam-fa <contaminant.fa>] [--decoy-fa <decoy.fa>]) \\\n\
    [--threads 8 (pipeline compute budget)] [--out <output-prefix>] [--k 21 (1..=31)]\n\
\n\
  viroflash version\n\
\n\
index: build a reusable index directory (ref.mmi + bloom.bin + manifest.json + targets.fa).\n\
  Without --decoy-fa, generate decoys using --decoy-ani/--decoy-per-layer/--decoy-seed.\n\
run: reuse --index (mutually exclusive with FASTA options; --k must match), or build one\n\
  automatically when --index is omitted.\n\
Input: sequencing fastq.gz; omit --r2 for single-end mode. Output: <out>.json / <out>.tsv / <out>.html / <out>.csv\n\
Performance: index and run emit <out>.perf.json / <out>.perf.tsv automatically.",
        env!("CARGO_PKG_VERSION")
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn args(list: &[&str]) -> Vec<String> {
        list.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn run_index_conflicts_with_fasta() {
        let err = parse_args_from(&args(&[
            "run",
            "--r1",
            "r.fq.gz",
            "--index",
            "idx",
            "--host-fa",
            "h.fa",
            "--target-fa",
            "t.fa",
        ]))
        .unwrap_err();
        assert!(err.contains("cannot be combined"), "err={err}");
    }

    #[test]
    fn run_requires_r1_and_target_source() {
        let err = parse_args_from(&args(&["run", "--index", "idx"])).unwrap_err();
        assert!(err.contains("--r1"), "err={err}");
        let err = parse_args_from(&args(&["run", "--r1", "r.fq.gz"])).unwrap_err();
        assert!(err.contains("--host-fa"), "err={err}");
    }

    #[test]
    fn run_accepts_index_only() {
        match parse_args_from(&args(&["run", "--r1", "r.fq.gz", "--index", "idx"])).unwrap() {
            Command::Run(opt) => {
                assert_eq!(opt.index, Some(PathBuf::from("idx")));
                assert!(opt.host_fa.is_none());
            }
            _ => panic!("expected Run"),
        }
    }

    #[test]
    fn index_requires_host_target_out() {
        let err =
            parse_args_from(&args(&["index", "--target-fa", "t.fa", "--out", "idx"])).unwrap_err();
        assert!(err.contains("--host-fa"), "err={err}");
        let err = parse_args_from(&args(&[
            "index",
            "--host-fa",
            "h.fa",
            "--target-fa",
            "t.fa",
        ]))
        .unwrap_err();
        assert!(err.contains("--out"), "err={err}");
    }

    #[test]
    fn index_decoy_fa_conflicts_with_tuning() {
        let err = parse_args_from(&args(&[
            "index",
            "--host-fa",
            "h.fa",
            "--target-fa",
            "t.fa",
            "--out",
            "idx",
            "--decoy-fa",
            "d.fa",
            "--decoy-ani",
            "82,88",
        ]))
        .unwrap_err();
        assert!(err.contains("cannot be combined"), "err={err}");
    }

    #[test]
    fn index_parses_ani_and_defaults() {
        match parse_args_from(&args(&[
            "index",
            "--host-fa",
            "h.fa",
            "--target-fa",
            "t.fa",
            "--out",
            "idx",
            "--decoy-ani",
            "82, 88",
            "--decoy-per-layer",
            "3",
            "--decoy-seed",
            "7",
        ]))
        .unwrap()
        {
            Command::Index(opt) => {
                assert_eq!(opt.decoy_anis, vec![82, 88]);
                assert_eq!(opt.decoy_per_layer, 3);
                assert_eq!(opt.decoy_seed, 7);
                assert!(opt.decoy_fa.is_none());
                assert_eq!(opt.contam_fa, None);
            }
            _ => panic!("expected Index"),
        }
    }

    #[test]
    fn invalid_values_and_unknown_args() {
        let err = parse_args_from(&args(&[
            "index",
            "--host-fa",
            "h",
            "--target-fa",
            "t",
            "--out",
            "o",
            "--k",
            "abc",
        ]))
        .unwrap_err();
        assert!(err.contains("Invalid"), "err={err}");
        let err = parse_args_from(&args(&[
            "index",
            "--host-fa",
            "h",
            "--target-fa",
            "t",
            "--out",
            "o",
            "--decoy-ani",
            "82x",
        ]))
        .unwrap_err();
        assert!(err.contains("Invalid"), "err={err}");
        let err = parse_args_from(&args(&["frobnicate"])).unwrap_err();
        assert!(err.contains("Unknown subcommand"), "err={err}");
        let err = parse_args_from(&args(&["run", "--r1", "r", "--wat", "x"])).unwrap_err();
        assert!(err.contains("Unknown argument"), "err={err}");
    }
}
