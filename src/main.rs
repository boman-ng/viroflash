use std::path::PathBuf;
use std::process::ExitCode;

use viroflash::{build_index, run_pipeline, IndexOptions, Precision, RunMode, RunOptions};

fn main() -> ExitCode {
    match run() {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("Error: {error}");
            ExitCode::FAILURE
        }
    }
}

enum Command {
    Help,
    Index(IndexOptions),
    Run(RunOptions),
}

fn run() -> Result<(), String> {
    match parse_args(&std::env::args().skip(1).collect::<Vec<_>>())? {
        Command::Help => {
            print_usage();
            Ok(())
        }
        Command::Index(options) => {
            let index = build_index(&options)?;
            eprintln!(
                "Index built: {} target groups; profile {}",
                index.target_groups.len(),
                index.profile_digest
            );
            Ok(())
        }
        Command::Run(options) => {
            let summary = run_pipeline(&options)?;
            eprintln!(
                "Completed {} fragments; selected {}; target rows {}; reports: {}",
                summary.input_fragments,
                summary.selected_fragments,
                summary.target_signal_rows,
                summary.report_dir.display()
            );
            Ok(())
        }
    }
}

fn parse_args(args: &[String]) -> Result<Command, String> {
    match args.first().map(String::as_str) {
        None | Some("--help" | "-h") => Ok(Command::Help),
        Some("index") => parse_index(&args[1..]).map(Command::Index),
        Some("run") => parse_run(&args[1..]).map(Command::Run),
        Some(command) => Err(format!("Unknown subcommand: {command}")),
    }
}

fn parse_index(args: &[String]) -> Result<IndexOptions, String> {
    let values = parse_flags(args, &["--host-fa", "--target-fa", "--out", "--threads"])?;
    Ok(IndexOptions {
        host_fa: required(&values, "--host-fa")?.into(),
        target_fa: required(&values, "--target-fa")?.into(),
        out_dir: required(&values, "--out")?.into(),
        threads: parse_threads(values.get("--threads"))?,
    })
}

fn parse_run(args: &[String]) -> Result<RunOptions, String> {
    let values = parse_flags(
        args,
        &[
            "--r1",
            "--r2",
            "--index",
            "--out",
            "--threads",
            "--precision",
            "--mode",
        ],
    )?;
    Ok(RunOptions {
        r1: required(&values, "--r1")?.into(),
        r2: values.get("--r2").map(PathBuf::from),
        index_dir: required(&values, "--index")?.into(),
        out_dir: required(&values, "--out")?.into(),
        threads: parse_threads(values.get("--threads"))?,
        precision: values
            .get("--precision")
            .map_or(Ok(Precision::Standard), |value| value.parse())?,
        mode: values
            .get("--mode")
            .map_or(Ok(RunMode::Screen), |value| value.parse())?,
    })
}

fn parse_flags<'a>(
    args: &'a [String],
    allowed: &[&str],
) -> Result<std::collections::HashMap<&'a str, &'a str>, String> {
    if args
        .iter()
        .any(|argument| argument == "--help" || argument == "-h")
    {
        return Err("Subcommand help is available with viroflash --help".into());
    }
    if !args.len().is_multiple_of(2) {
        return Err(format!("Missing value for {}", args.last().unwrap()));
    }
    let mut values = std::collections::HashMap::new();
    for pair in args.chunks_exact(2) {
        if !allowed.contains(&pair[0].as_str()) {
            return Err(format!("Unknown argument: {}", pair[0]));
        }
        if values.insert(pair[0].as_str(), pair[1].as_str()).is_some() {
            return Err(format!("Duplicate argument: {}", pair[0]));
        }
    }
    Ok(values)
}

fn required<'a>(
    values: &'a std::collections::HashMap<&str, &'a str>,
    flag: &str,
) -> Result<&'a str, String> {
    values
        .get(flag)
        .copied()
        .ok_or_else(|| format!("Missing required argument: {flag}"))
}
fn parse_threads(value: Option<&&str>) -> Result<usize, String> {
    value
        .map_or(Ok(1), |value| {
            value
                .parse::<usize>()
                .map_err(|_| format!("Invalid --threads value: {value}"))
        })
        .and_then(|threads| {
            if threads == 0 {
                Err("--threads must be greater than zero".into())
            } else {
                Ok(threads)
            }
        })
}

fn print_usage() {
    println!("Viroflash {}\n\nUSAGE:\n  viroflash index --host-fa HOST.fa --target-fa TARGET.fa --out INDEX_DIR [--threads N]\n  viroflash run --r1 SAMPLE_R1.fastq.gz [--r2 SAMPLE_R2.fastq.gz] --index INDEX_DIR --out SAMPLE_REPORT_DIR [--threads N] [--mode full|screen] [--precision fast|standard|sensitive]\n\nScreen mode (default): sample input fragments, then Bloom-screen and competitively align against HOST+TARGET.\nFull mode: Bloom-screen all fragments, then sample and competitively align candidates.\nBoth modes use deterministic bottom-k sampling with a capacity derived from precision and panel size. No candidate files are written.\nPrecision: fast 100 ppm, standard 10 ppm (default), sensitive 1 ppm. Reports show original-library abundance percentages, sampling intervals and a score relative to the selected ppm target.\nSuccessful runs write report.csv, report.html and perf.json. HTML shows the top 20 supported reference groups and embeds the complete CSV.", env!("CARGO_PKG_VERSION"));
}

#[cfg(test)]
mod tests {
    use super::*;
    fn strings(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn cli_accepts_only_current_contract() {
        assert!(matches!(
            parse_args(&strings(&[
                "index",
                "--host-fa",
                "h",
                "--target-fa",
                "t",
                "--out",
                "i"
            ])),
            Ok(Command::Index(_))
        ));
        assert!(matches!(
            parse_args(&strings(&[
                "run", "--r1", "r", "--index", "i", "--out", "o"
            ])),
            Ok(Command::Run(_))
        ));
        for (preset, expected) in [
            ("fast", Precision::Fast),
            ("standard", Precision::Standard),
            ("sensitive", Precision::Sensitive),
        ] {
            let Command::Run(options) = parse_args(&strings(&[
                "run",
                "--r1",
                "r",
                "--index",
                "i",
                "--out",
                "o",
                "--precision",
                preset,
            ]))
            .unwrap() else {
                panic!("run command")
            };
            assert_eq!(options.precision, expected);
        }
        assert!(parse_args(&strings(&[
            "run",
            "--r1",
            "r",
            "--index",
            "i",
            "--out",
            "o",
            "--precision",
            "0.05"
        ]))
        .is_err());
        for command in ["index", "run"] {
            let error = parse_args(&strings(&[command, "--unknown", "unknown.fa"]))
                .err()
                .unwrap();
            assert!(error.contains("Unknown argument: --unknown"), "{error}");
        }
        let Command::Run(defaults) = parse_args(&strings(&[
            "run", "--r1", "r", "--index", "i", "--out", "o",
        ]))
        .unwrap() else {
            panic!("run command")
        };
        assert_eq!(defaults.mode, RunMode::Screen);
        assert_eq!(defaults.precision, Precision::Standard);
        for mode in ["full", "screen"] {
            let Command::Run(options) = parse_args(&strings(&[
                "run",
                "--r1",
                "r",
                "--index",
                "i",
                "--out",
                "o",
                "--mode",
                mode,
                "--precision",
                "fast",
            ]))
            .unwrap() else {
                panic!("run command");
            };
            assert_eq!(options.mode.as_str(), mode);
            assert_eq!(options.precision, Precision::Fast);
        }
        assert!(parse_args(&strings(&[
            "run", "--r1", "r", "--index", "i", "--out", "o", "--mode", "standard"
        ]))
        .is_err());
    }
}
