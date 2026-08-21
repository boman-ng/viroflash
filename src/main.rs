//! 命令行入口：`viroflash run` / `viroflash version`。
//! 手工解析参数，保持零额外依赖（轻量约束）。

use std::path::PathBuf;
use std::process::ExitCode;

use viroflash::decoy::{self, DecoyOptions};
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

enum Command {
    Help,
    Version,
    Run(Options),
    Decoy(DecoyOptions),
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
        Command::Decoy(opt) => {
            let msg = decoy::generate(&opt)?;
            eprintln!("{msg}");
            Ok(())
        }
    }
}

fn parse_args() -> Result<Command, String> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    if args.is_empty() || args[0] == "-h" || args[0] == "--help" {
        return Ok(Command::Help);
    }
    match args[0].as_str() {
        "version" | "--version" => Ok(Command::Version),
        "run" => {
            let mut opt = Options::default();
            let mut i = 1;
            while i < args.len() {
                let name = args[i].as_str();
                let value = args
                    .get(i + 1)
                    .ok_or_else(|| format!("参数 {name} 缺少值"))?;
                match name {
                    "--r1" => opt.r1 = PathBuf::from(value),
                    "--r2" => opt.r2 = Some(PathBuf::from(value)),
                    "--host-fa" => opt.host_fa = PathBuf::from(value),
                    "--target-fa" => opt.target_fa = PathBuf::from(value),
                    "--decoy-fa" => opt.decoy_fa = PathBuf::from(value),
                    "--contam-fa" => opt.contam_fa = PathBuf::from(value),
                    "--threads" => {
                        opt.threads = value
                            .parse()
                            .map_err(|_| format!("--threads 取值非法: {value}"))?
                    }
                    "--out" => opt.out = PathBuf::from(value),
                    "--k" => {
                        opt.k = value
                            .parse()
                            .map_err(|_| format!("--k 取值非法: {value}"))?
                    }
                    other => return Err(format!("未知参数: {other}")),
                }
                i += 2;
            }
            for (label, p) in [
                ("--r1", &opt.r1),
                ("--host-fa", &opt.host_fa),
                ("--target-fa", &opt.target_fa),
                ("--decoy-fa", &opt.decoy_fa),
                ("--contam-fa", &opt.contam_fa),
            ] {
                if p.as_os_str().is_empty() {
                    return Err(format!("run 子命令缺少 {label}"));
                }
            }
            Ok(Command::Run(opt))
        }
        "decoy" => {
            let mut opt = DecoyOptions::default();
            let mut i = 1;
            while i < args.len() {
                let name = args[i].as_str();
                let value = args
                    .get(i + 1)
                    .ok_or_else(|| format!("参数 {name} 缺少值"))?;
                match name {
                    "--target-fa" => opt.target_fa = PathBuf::from(value),
                    "--out" => opt.out = PathBuf::from(value),
                    "--ani" => {
                        opt.anis = value
                            .split(',')
                            .map(|s| s.trim().parse().map_err(|_| format!("--ani 取值非法: {s}")))
                            .collect::<Result<Vec<u8>, _>>()?
                    }
                    "--per-layer" => {
                        opt.per_layer = value
                            .parse()
                            .map_err(|_| format!("--per-layer 取值非法: {value}"))?
                    }
                    "--seed" => {
                        opt.seed = value
                            .parse()
                            .map_err(|_| format!("--seed 取值非法: {value}"))?
                    }
                    "--report" => opt.report = Some(PathBuf::from(value)),
                    other => return Err(format!("未知参数: {other}")),
                }
                i += 2;
            }
            if opt.target_fa.as_os_str().is_empty() || opt.out.as_os_str().is_empty() {
                return Err("decoy 子命令缺少 --target-fa 或 --out".into());
            }
            Ok(Command::Decoy(opt))
        }
        other => Err(format!("未知子命令: {other}")),
    }
}

fn print_usage() {
    println!(
        "viroflash {} — 病毒候选检测命令行工具\n\
\n\
用法:\n\
  viroflash run \\\n\
    --r1 <reads_R1.fastq.gz> [--r2 <reads_R2.fastq.gz>] \\\n\
    --host-fa <宿主.fa> --target-fa <目标病毒.fa> \\\n\
    --decoy-fa <诱饵.fa> --contam-fa <污染.fa> \\\n\
    [--threads 8（进程总线程数）] [--out <输出前缀>] [--k 21（1..=31）]\n\
\n\
  viroflash decoy \\\n\
    --target-fa <目标病毒.fa> --out <诱饵.fa> \\\n\
    [--ani 82,85,88] [--per-layer 4] [--seed 0] [--report <decoy.tsv>]\n\
\n\
  viroflash version\n\
\n\
输入: 测序 fastq.gz（省略 --r2 即单端模式）+ 4 类基因组 FASTA；输出: <out>.json / <out>.tsv\n\
decoy: 从目标基因组离线生成 SNP-only 突变诱饵（ANI 分层，确定性 seed）",
        env!("CARGO_PKG_VERSION")
    );
}
