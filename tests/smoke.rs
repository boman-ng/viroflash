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

fn write_gzip_members(path: &Path, members: &[&[u8]]) {
    let mut file = File::create(path).unwrap();
    for member in members {
        let mut writer = GzEncoder::new(Vec::new(), Compression::fast());
        writer.write_all(member).unwrap();
        file.write_all(&writer.finish().unwrap()).unwrap();
    }
}

fn csv_rows(text: &str) -> (Vec<&str>, Vec<Vec<&str>>) {
    let mut lines = text.lines();
    let header = lines.next().unwrap().split(',').collect();
    let rows = lines
        .map(|line| line.split(',').collect::<Vec<_>>())
        .collect();
    (header, rows)
}

fn visible_fields(fragment: &str) -> Vec<(String, String)> {
    fragment
        .split(r#"<tr data-field=""#)
        .skip(1)
        .map(|row| {
            let (field, remainder) = row.split_once(r#""><th>"#).unwrap();
            let value = remainder
                .split_once("<td>")
                .unwrap()
                .1
                .split_once("</td>")
                .unwrap()
                .0;
            (field.to_string(), value.to_string())
        })
        .collect()
}

fn contract_fields(record_type: &str) -> Vec<String> {
    std::fs::read_to_string("tests/fixtures/output-fields.tsv")
        .unwrap()
        .lines()
        .skip(1)
        .filter_map(|line| {
            let columns = line.split('\t').collect::<Vec<_>>();
            (columns[0] == "report.csv" && columns[1] == record_type)
                .then(|| columns[2].to_string())
        })
        .collect()
}

fn assert_csv_and_visible_html_share_all_fields(csv: &str, html: &str) {
    assert!(!html.contains("application/json"));
    let (header, rows) = csv_rows(csv);
    let run_section = html
        .split_once(r#"<section id="run-integrity">"#)
        .unwrap()
        .1
        .split_once(r#"<section id="observed-target-signals">"#)
        .unwrap()
        .0;
    let visible_run = visible_fields(run_section);
    assert_eq!(
        visible_run
            .iter()
            .map(|(field, _)| field.clone())
            .collect::<Vec<_>>(),
        contract_fields("RUN")
    );
    for (field, value) in visible_run {
        let column = header
            .iter()
            .position(|candidate| **candidate == field)
            .unwrap();
        assert_eq!(value, rows[0][column], "RUN field {field}");
    }

    let target_area = html
        .split_once(r#"<section id="observed-target-signals">"#)
        .unwrap()
        .1
        .split_once("<section><h2>Evidence detail</h2>")
        .unwrap()
        .0;
    let visible_targets = target_area
        .split(r#"<section class="target-signal""#)
        .skip(1)
        .map(|signal| visible_fields(signal.split_once("</section>").unwrap().0))
        .collect::<Vec<_>>();
    assert_eq!(visible_targets.len(), rows.len() - 1);
    for (visible, row) in visible_targets.iter().zip(&rows[1..]) {
        assert_eq!(
            visible
                .iter()
                .map(|(field, _)| field.clone())
                .collect::<Vec<_>>(),
            contract_fields("TARGET_SIGNAL")
        );
        for (field, value) in visible {
            let column = header
                .iter()
                .position(|candidate| *candidate == field)
                .unwrap();
            assert_eq!(value, row[column], "TARGET_SIGNAL field {field}");
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
        let output_fields = std::fs::read_to_string("tests/fixtures/output-fields.tsv").unwrap();
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
    assert_csv_and_visible_html_share_all_fields(&text, &html);
    let run = lines[1].split(',').collect::<Vec<_>>();
    let run_value = |field: &str| {
        run[header
            .iter()
            .position(|candidate| *candidate == field)
            .unwrap()]
    };
    assert!(
        run_value("selected_fragments").parse::<u64>().unwrap()
            <= run_value("prescreen_passed_fragments")
                .parse::<u64>()
                .unwrap()
    );
    assert_eq!(
        run_value("selected_fragments"),
        run_value("aligned_fragments")
    );

    let output_fields = std::fs::read_to_string("tests/fixtures/output-fields.tsv").unwrap();
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
fn phase5_report_is_invariant_across_threads_reloads_and_gzip_segmentation() {
    let (root, target) = build_fixture("phase5-invariance");
    let fastq = fastq_bytes(&target[..120], 40);
    let split = fastq
        .windows(2)
        .enumerate()
        .filter_map(|(index, bytes)| (bytes == b"\n@").then_some(index + 1))
        .nth(19)
        .unwrap();

    let plain = root.join("plain/sample.fastq");
    let single = root.join("single/sample.fastq.gz");
    let multiple = root.join("multiple/sample.fastq.gz");
    for path in [&plain, &single, &multiple] {
        std::fs::create_dir(path.parent().unwrap()).unwrap();
    }
    std::fs::write(&plain, &fastq).unwrap();
    write_gzip_members(&single, &[&fastq]);
    write_gzip_members(&multiple, &[&fastq[..split], &fastq[split..]]);

    let mut expected = None;
    for (encoding, input) in [
        ("plain", &plain),
        ("single", &single),
        ("multiple", &multiple),
    ] {
        for threads in [1, 2, 4, 8] {
            for repeat in 1..=2 {
                let output_dir = root.join(format!("out-{encoding}-{threads}-{repeat}"));
                let output = command(&[
                    "run",
                    "--r1",
                    input.to_str().unwrap(),
                    "--index",
                    root.join("index").to_str().unwrap(),
                    "--out",
                    output_dir.to_str().unwrap(),
                    "--threads",
                    &threads.to_string(),
                ]);
                assert!(
                    output.status.success(),
                    "{}",
                    String::from_utf8_lossy(&output.stderr)
                );
                let report = std::fs::read(output_dir.join("report.csv")).unwrap();
                if let Some(expected) = &expected {
                    assert_eq!(&report, expected);
                } else {
                    expected = Some(report);
                }
            }
        }
    }
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
    let html = std::fs::read_to_string(root.join("out/report.html")).unwrap();
    assert_csv_and_visible_html_share_all_fields(&csv, &html);
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
    assert_eq!(selected, 0);
    assert_eq!(
        value("reason_codes"),
        "TARGET_KMER_NOT_EVALUABLE_INPUT_FRAGMENTS=100"
    );
    let _ = std::fs::remove_dir_all(root);
}
