use std::collections::BTreeSet;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use flate2::write::GzEncoder;
use flate2::Compression;

static NEXT: AtomicU64 = AtomicU64::new(0);

fn workspace(tag: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "viroflash-{tag}-{}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&path).unwrap();
    path
}

fn sequence(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[(state & 3) as usize]
        })
        .collect()
}

fn write_fasta(path: &Path, id: &str, sequence: &[u8]) {
    let mut file = File::create(path).unwrap();
    writeln!(file, ">{id}").unwrap();
    for chunk in sequence.chunks(60) {
        writeln!(file, "{}", String::from_utf8_lossy(chunk)).unwrap();
    }
}

fn write_fastq(path: &Path, target: &[u8], mate: Option<u8>) {
    let file = File::create(path).unwrap();
    let mut writer = GzEncoder::new(file, Compression::fast());
    for index in 0..20 {
        let start = index % 30;
        let read = &target[start..start + 120];
        let suffix = mate.map_or(String::new(), |mate| format!("/{mate}"));
        writeln!(
            writer,
            "@fragment-{index}{suffix}\n{}\n+\n{}",
            String::from_utf8_lossy(read),
            "I".repeat(read.len())
        )
        .unwrap();
    }
    writer.finish().unwrap();
}

fn command(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_viroflash"))
        .args(args)
        .output()
        .unwrap()
}

fn build_fixture(tag: &str) -> (PathBuf, Vec<u8>) {
    let root = workspace(tag);
    let host = sequence(17, 600);
    let target = sequence(91, 600);
    write_fasta(&root.join("host.fa"), "host", &host);
    write_fasta(&root.join("target.fa"), "target", &target);
    let output = command(&[
        "index",
        "--host-fa",
        root.join("host.fa").to_str().unwrap(),
        "--target-fa",
        root.join("target.fa").to_str().unwrap(),
        "--out",
        root.join("index").to_str().unwrap(),
        "--threads",
        "2",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    (root, target)
}

fn output_names(path: &Path) -> BTreeSet<String> {
    std::fs::read_dir(path)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .collect()
}

#[test]
fn cli_index_and_se_run_produce_three_source_consistent_files() {
    let (root, target) = build_fixture("se");
    write_fastq(&root.join("sample.fastq.gz"), &target, None);
    for (threads, output_name) in [(1, "out-one"), (4, "out-four")] {
        let output = command(&[
            "run",
            "--r1",
            root.join("sample.fastq.gz").to_str().unwrap(),
            "--index",
            root.join("index").to_str().unwrap(),
            "--out",
            root.join(output_name).to_str().unwrap(),
            "--threads",
            &threads.to_string(),
        ]);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(
            output_names(&root.join(output_name)),
            BTreeSet::from([
                "perf.json".into(),
                "report.csv".into(),
                "report.html".into()
            ])
        );
        let perf: serde_json::Value = serde_json::from_slice(
            &std::fs::read(root.join(output_name).join("perf.json")).unwrap(),
        )
        .unwrap();
        let actual_fields = perf
            .as_object()
            .unwrap()
            .keys()
            .cloned()
            .collect::<BTreeSet<_>>();
        let output_fields = std::fs::read_to_string("evaluation/phase0/output-fields.tsv").unwrap();
        let expected_fields = output_fields
            .lines()
            .skip(1)
            .filter_map(|line| {
                let columns = line.split('\t').collect::<Vec<_>>();
                (columns[0] == "perf.json").then(|| columns[2].to_string())
            })
            .collect::<BTreeSet<_>>();
        assert_eq!(actual_fields, expected_fields);
    }
    let csv = std::fs::read(root.join("out-one/report.csv")).unwrap();
    assert_eq!(
        csv,
        std::fs::read(root.join("out-four/report.csv")).unwrap()
    );
    let text = String::from_utf8(csv).unwrap();
    let lines = text.lines().collect::<Vec<_>>();
    assert_eq!(
        lines
            .iter()
            .filter(|line| line.contains("TARGET_SIGNAL"))
            .count(),
        1
    );
    let header = lines[0].split(',').collect::<Vec<_>>();
    let values = lines[2].split(',').collect::<Vec<_>>();
    let value = |field: &str| {
        values[header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap()]
    };
    let html = std::fs::read_to_string(root.join("out-one/report.html")).unwrap();
    assert!(html.contains(&format!(
        "data-fraction=\"{}\"",
        value("attributed_fragment_fraction")
    )));
    assert!(html.contains(&format!("data-lower=\"{}\"", value("interval_lower"))));
    assert!(html.contains(&format!("data-upper=\"{}\"", value("interval_upper"))));

    let output_fields = std::fs::read_to_string("evaluation/phase0/output-fields.tsv").unwrap();
    let expected = output_fields
        .lines()
        .skip(1)
        .filter_map(|line| {
            let columns = line.split('\t').collect::<Vec<_>>();
            (columns[0] == "report.csv").then(|| columns[2].to_string())
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(
        header
            .into_iter()
            .map(str::to_string)
            .collect::<BTreeSet<_>>(),
        expected
    );
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn cli_pe_run_counts_fragments_once() {
    let (root, target) = build_fixture("pe");
    write_fastq(&root.join("sample_R1.fastq.gz"), &target, Some(1));
    write_fastq(&root.join("sample_R2.fastq.gz"), &target, Some(2));
    let output = command(&[
        "run",
        "--r1",
        root.join("sample_R1.fastq.gz").to_str().unwrap(),
        "--r2",
        root.join("sample_R2.fastq.gz").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out").to_str().unwrap(),
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let csv = std::fs::read_to_string(root.join("out/report.csv")).unwrap();
    let header = csv.lines().next().unwrap().split(',').collect::<Vec<_>>();
    let run = csv.lines().nth(1).unwrap().split(',').collect::<Vec<_>>();
    let value = |field: &str| {
        run[header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap()]
    };
    assert_eq!(value("input_mode"), "PE");
    assert_eq!(value("input_fragments"), "20");
    assert_eq!(value("read_ends_per_fragment"), "2");
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn failed_run_leaves_only_error_perf_json() {
    let (root, _) = build_fixture("failure");
    std::fs::write(root.join("broken.fastq"), b"@broken\nACGT\n+\n").unwrap();
    let output = command(&[
        "run",
        "--r1",
        root.join("broken.fastq").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out").to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert_eq!(
        output_names(&root.join("out")),
        BTreeSet::from(["perf.json".into()])
    );
    let perf: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("out/perf.json")).unwrap()).unwrap();
    assert_eq!(perf["status"], "ERROR");
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn run_requires_current_reusable_index() {
    let root = workspace("invalid-index");
    std::fs::create_dir(root.join("index")).unwrap();
    std::fs::write(root.join("index/manifest.json"), b"{}\n").unwrap();
    std::fs::write(root.join("sample.fastq"), b"@x\nACGT\n+\nIIII\n").unwrap();
    let output = command(&[
        "run",
        "--r1",
        root.join("sample.fastq").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out").to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("rebuild"));
    assert!(
        !command(&["run", "--r1", "r", "--unexpected", "x", "--out", "o"])
            .status
            .success()
    );
    let _ = std::fs::remove_dir_all(root);
}
