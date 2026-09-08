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

fn fastq_bytes(sequence: &[u8], count: usize) -> Vec<u8> {
    let mut bytes = Vec::new();
    for index in 0..count {
        writeln!(
            bytes,
            "@fragment-{index}\n{}\n+\n{}",
            String::from_utf8_lossy(sequence),
            "I".repeat(sequence.len())
        )
        .unwrap();
    }
    bytes
}

fn csv_rows(text: &str) -> (Vec<&str>, Vec<Vec<&str>>) {
    let mut lines = text.lines();
    let header = lines.next().unwrap().split(',').collect();
    let rows = lines
        .map(|line| line.split(',').collect::<Vec<_>>())
        .collect();
    (header, rows)
}

fn assert_model_cell(field: &str, value: &serde_json::Value, cell: &str) {
    match value {
        serde_json::Value::String(value) => assert_eq!(cell, value),
        serde_json::Value::Number(value) if value.is_f64() => {
            let csv_value = cell.parse::<f64>().unwrap();
            let model_value = value.as_f64().unwrap();
            assert!(
                (csv_value - model_value).abs() <= f64::EPSILON,
                "field {field}: CSV={csv_value:?}, model={model_value:?}"
            )
        }
        serde_json::Value::Number(value) => assert_eq!(cell, value.to_string()),
        serde_json::Value::Array(values) => assert_eq!(
            cell,
            values
                .iter()
                .map(|value| value.as_str().unwrap())
                .collect::<Vec<_>>()
                .join(";")
        ),
        _ => panic!("unsupported report model value: {value}"),
    }
}

fn assert_csv_and_html_share_complete_model(csv: &str, html: &str) {
    let marker = r#"<script id="evidence-report" type="application/json">"#;
    let embedded = html
        .split_once(marker)
        .unwrap()
        .1
        .split_once("</script>")
        .unwrap()
        .0;
    let model: serde_json::Value = serde_json::from_str(embedded).unwrap();
    let (header, rows) = csv_rows(csv);
    let run = &rows[0];
    for (field, value) in model["run"].as_object().unwrap() {
        let column = header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap();
        assert_model_cell(field, value, run[column]);
    }
    let signals = model["target_signals"].as_array().unwrap();
    assert_eq!(signals.len(), rows.len() - 1);
    for (signal, row) in signals.iter().zip(&rows[1..]) {
        for (field, value) in signal.as_object().unwrap() {
            let column = header
                .iter()
                .position(|candidate| *candidate == field)
                .unwrap();
            assert_model_cell(field, value, row[column]);
        }
    }
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
    let html = std::fs::read_to_string(root.join("out-one/report.html")).unwrap();
    assert_csv_and_html_share_complete_model(&text, &html);
    let run = lines[1].split(',').collect::<Vec<_>>();
    let run_value = |field: &str| {
        run[header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap()]
    };
    assert_eq!(
        run_value("selected_fragments"),
        run_value("prescreen_passed_fragments")
    );
    assert_eq!(
        run_value("selected_fragments"),
        run_value("aligned_fragments")
    );

    let output_fields = std::fs::read_to_string("evaluation/phase0/output-fields.tsv").unwrap();
    let mut seen = BTreeSet::new();
    let expected = output_fields
        .lines()
        .skip(1)
        .filter_map(|line| {
            let columns = line.split('\t').collect::<Vec<_>>();
            (columns[0] == "report.csv" && seen.insert(columns[2])).then_some(columns[2])
        })
        .collect::<Vec<_>>();
    assert_eq!(header, expected);
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

#[test]
fn run_rejects_manifest_group_mirror_and_non_dense_contigs() {
    let (root, target) = build_fixture("manifest-tamper");
    write_fastq(&root.join("sample.fastq.gz"), &target, None);
    let manifest_path = root.join("index/manifest.json");
    let mut manifest: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&manifest_path).unwrap()).unwrap();
    manifest["target_groups"] = serde_json::json!([{"ordinal": 0}, {"ordinal": 1}]);
    std::fs::write(
        &manifest_path,
        serde_json::to_vec_pretty(&manifest).unwrap(),
    )
    .unwrap();
    let output = command(&[
        "run",
        "--r1",
        root.join("sample.fastq.gz").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out-mirror").to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("rebuild"));

    manifest.as_object_mut().unwrap().remove("target_groups");
    manifest["contigs"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|contig| contig["role"] == "Target")
        .unwrap()["target_group_ordinal"] = serde_json::json!(1);
    std::fs::write(
        &manifest_path,
        serde_json::to_vec_pretty(&manifest).unwrap(),
    )
    .unwrap();
    let output = command(&[
        "run",
        "--r1",
        root.join("sample.fastq.gz").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out-contig").to_str().unwrap(),
    ]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("dense"));
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn short_selected_fragments_are_reported_as_a_quantified_limitation() {
    let (root, target) = build_fixture("short-reads");
    std::fs::write(root.join("short.fastq"), fastq_bytes(&target[..20], 100)).unwrap();
    let output = command(&[
        "run",
        "--r1",
        root.join("short.fastq").to_str().unwrap(),
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
    let (header, rows) = csv_rows(&csv);
    assert_eq!(rows.len(), 1);
    let run = &rows[0];
    let value = |field: &str| {
        run[header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap()]
    };
    assert_eq!(value("analysis_status"), "CONFORMANT_WITH_LIMITATIONS");
    assert_eq!(value("prescreen_passed_fragments"), "0");
    let selected = value("selected_fragments").parse::<u64>().unwrap();
    assert!(selected > 0);
    assert_eq!(
        value("reason_codes"),
        format!("TARGET_KMER_NOT_EVALUABLE_SELECTED_FRAGMENTS={selected}")
    );
    let _ = std::fs::remove_dir_all(root);
}

#[cfg(target_os = "linux")]
#[test]
fn changed_fifo_bytes_between_passes_are_rejected() {
    let (root, target) = build_fixture("fifo-change");
    let fifo = root.join("sample.fastq");
    assert!(Command::new("mkfifo")
        .arg(&fifo)
        .status()
        .unwrap()
        .success());
    let host = sequence(17, 600);
    std::fs::write(root.join("pass1.fastq"), fastq_bytes(&host[..120], 20)).unwrap();
    std::fs::write(root.join("pass2.fastq"), fastq_bytes(&target[..120], 20)).unwrap();
    let writer = Command::new("python3")
        .arg("-c")
        .arg(
            r#"import errno, os, sys, time
fifo, first, second = sys.argv[1:]
deadline = time.monotonic() + 10
def open_writer():
    while time.monotonic() < deadline:
        try:
            return os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
        except OSError as error:
            if error.errno != errno.ENXIO:
                raise
            time.sleep(0.01)
    raise TimeoutError('FIFO reader did not open')
for index, path in enumerate((first, second)):
    descriptor = open_writer()
    with os.fdopen(descriptor, 'wb') as stream, open(path, 'rb') as source:
        stream.write(source.read())
    if index == 0:
        time.sleep(0.5)
"#,
        )
        .arg(&fifo)
        .arg(root.join("pass1.fastq"))
        .arg(root.join("pass2.fastq"))
        .spawn()
        .unwrap();
    let output = Command::new("timeout")
        .arg("15s")
        .arg(env!("CARGO_BIN_EXE_viroflash"))
        .args([
            "run",
            "--r1",
            fifo.to_str().unwrap(),
            "--index",
            root.join("index").to_str().unwrap(),
            "--out",
            root.join("out").to_str().unwrap(),
            "--threads",
            "2",
        ])
        .output()
        .unwrap();
    assert!(writer.wait_with_output().unwrap().status.success());
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("bytes changed between pass"));
    assert_eq!(
        output_names(&root.join("out")),
        BTreeSet::from(["perf.json".into()])
    );
    let _ = std::fs::remove_dir_all(root);
}
