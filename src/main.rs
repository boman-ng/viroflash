//! 命令行入口：`viroflash index` / `viroflash run` / `viroflash version`。
//! 手工解析参数，保持零额外依赖（轻量约束）。

use std::path::PathBuf;
use std::process::ExitCode;

use viroflash::index::{self, IndexOptions};
use viroflash::reference::Role;
use viroflash::{run_pipeline, Options};

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("错误: {e}");
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
                "索引已构建: {}（contig: host={} target={} decoy={} contam={}；k={}）",
                opt.out_dir.display(),
                n(Role::Host),
                n(Role::Target),
                n(Role::Decoy),
                n(Role::Contaminant),
                opt.k
            );
            eprintln!(
                "  ref.mmi / bloom.bin / manifest.json（manifest_blake3={}）",
                built.manifest_blake3
            );
            Ok(())
        }
        Command::Run(opt) => {
            let summary = run_pipeline(&opt)?;
            eprintln!("输入 pairs: {}", summary.input_pairs);
            eprintln!("预筛通过 pairs: {}", summary.prescreen_pairs);
            if summary.map_errors > 0 {
                eprintln!(
                    "警告: {} 个 fragment 比对失败（详见报告 run.map_errors）",
                    summary.map_errors
                );
            }
            eprintln!("候选数: {}", summary.candidates.len());
            for c in &summary.candidates {
                eprintln!(
                    "  {}  q={:.4}  {}  覆盖={:.2}  reads={}  split={}  discordant={}",
                    c.contig,
                    c.q_value,
                    c.confidence(),
                    c.covered_frac,
                    c.reads,
                    c.split_events,
                    c.discordant
                );
            }
            eprintln!(
                "报告: {} / {}",
                summary.result_json.display(),
                summary.result_tsv.display()
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
        other => Err(format!("未知子命令: {other}")),
    }
}

fn parse_num<T: std::str::FromStr>(flag: &str, value: &str) -> Result<T, String> {
    value
        .parse()
        .map_err(|_| format!("{flag} 取值非法: {value}"))
}

fn parse_anis(value: &str) -> Result<Vec<u8>, String> {
    value
        .split(',')
        .map(|s| {
            s.trim()
                .parse()
                .map_err(|_| format!("--decoy-ani 取值非法: {s}"))
        })
        .collect()
}

/// `viroflash index`：构建可复用索引目录；未提供 --decoy-fa 时自动生成诱饵。
fn parse_index(args: &[String]) -> Result<Command, String> {
    let mut opt = IndexOptions::default();
    let mut decoy_tuning = false;
    let mut i = 0;
    while i < args.len() {
        let name = args[i].as_str();
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("参数 {name} 缺少值"))?;
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
            other => return Err(format!("未知参数: {other}")),
        }
        i += 2;
    }
    if opt.decoy_fa.is_some() && decoy_tuning {
        return Err(
            "--decoy-fa 与 --decoy-ani/--decoy-per-layer/--decoy-seed 不能同时使用（诱饵调参仅对自动生成有效）"
                .into(),
        );
    }
    if opt.host_fa.as_os_str().is_empty() || opt.target_fa.as_os_str().is_empty() {
        return Err("index 子命令缺少 --host-fa 或 --target-fa".into());
    }
    if opt.out_dir.as_os_str().is_empty() {
        return Err("index 子命令缺少 --out".into());
    }
    Ok(Command::Index(opt))
}

/// `viroflash run`：--index 复用索引，或直接给 FASTA 自动构建（互斥）。
fn parse_run(args: &[String]) -> Result<Command, String> {
    let mut opt = Options::default();
    let mut i = 0;
    while i < args.len() {
        let name = args[i].as_str();
        let value = args
            .get(i + 1)
            .ok_or_else(|| format!("参数 {name} 缺少值"))?;
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
            other => return Err(format!("未知参数: {other}")),
        }
        i += 2;
    }
    if opt.r1.as_os_str().is_empty() {
        return Err("run 子命令缺少 --r1".into());
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
                    "--index 与 {label} 不能同时使用（加载索引时参考信息取自 manifest）"
                ));
            }
        }
    } else if opt.host_fa.is_none() || opt.target_fa.is_none() {
        return Err(
            "run 子命令缺少 --host-fa/--target-fa（未提供 --index 时自动构建索引需要它们）".into(),
        );
    }
    Ok(Command::Run(opt))
}

fn print_usage() {
    println!(
        "viroflash {} — 病毒候选检测命令行工具\n\
\n\
用法:\n\
  viroflash index \\\n\
    --host-fa <宿主.fa> --target-fa <目标病毒.fa> \\\n\
    [--contam-fa <污染.fa>] [--decoy-fa <诱饵.fa>] \\\n\
    [--decoy-ani 82,85,88] [--decoy-per-layer 4] [--decoy-seed 0] \\\n\
    --out <索引目录> [--k 21（1..=31）] [--threads 8]\n\
\n\
  viroflash run \\\n\
    --r1 <reads_R1.fastq.gz> [--r2 <reads_R2.fastq.gz>] \\\n\
    (--index <索引目录> | --host-fa <宿主.fa> --target-fa <目标病毒.fa> \\\n\
      [--contam-fa <污染.fa>] [--decoy-fa <诱饵.fa>]) \\\n\
    [--threads 8（进程总线程数）] [--out <输出前缀>] [--k 21（1..=31）]\n\
\n\
  viroflash version\n\
\n\
index: 构建可复用索引目录（ref.mmi + bloom.bin + manifest.json）；未提供\n\
  --decoy-fa 时按 --decoy-ani/--decoy-per-layer/--decoy-seed 从目标自动生成诱饵。\n\
run: --index 复用索引（与 FASTA 参数互斥，--k 需与索引一致）；未提供 --index\n\
  时自动构建（诱饵未提供时按默认参数自动生成）。\n\
输入: 测序 fastq.gz（省略 --r2 即单端模式）；输出: <out>.json / <out>.tsv",
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
        assert!(err.contains("不能同时使用"), "err={err}");
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
            _ => panic!("应为 Run"),
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
        assert!(err.contains("不能同时使用"), "err={err}");
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
            _ => panic!("应为 Index"),
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
        assert!(err.contains("非法"), "err={err}");
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
        assert!(err.contains("非法"), "err={err}");
        let err = parse_args_from(&args(&["frobnicate"])).unwrap_err();
        assert!(err.contains("未知子命令"), "err={err}");
        let err = parse_args_from(&args(&["run", "--r1", "r", "--wat", "x"])).unwrap_err();
        assert!(err.contains("未知参数"), "err={err}");
    }
}
