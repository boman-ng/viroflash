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

fn csv_rows(text: &str) -> (Vec<String>, Vec<Vec<String>>) {
    let mut reader = csv::Reader::from_reader(text.trim_start_matches('\u{feff}').as_bytes());
    let header = reader
        .headers()
        .unwrap()
        .iter()
        .map(str::to_string)
        .collect();
    let rows = reader
        .records()
        .map(|row| row.unwrap().iter().map(str::to_string).collect())
        .collect();
    (header, rows)
}

fn assert_csv_and_visible_html_share_all_fields(csv: &str, html: &str) {
    use scraper::{Html, Selector};
    let (header, rows) = csv_rows(csv);
    let expected = include_str!("fixtures/output-fields.tsv")
        .lines()
        .skip(1)
        .filter_map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            (fields[0] == "report.csv").then_some(fields[2])
        })
        .collect::<Vec<_>>();
    assert_eq!(header, expected);
    let document = Html::parse_document(html);
    let row_selector = Selector::parse("dl[data-research-row]").unwrap();
    let field_selector = Selector::parse("dd[data-field]").unwrap();
    let html_rows = document
        .select(&row_selector)
        .map(|row| {
            let fields = row.select(&field_selector).collect::<Vec<_>>();
            assert_eq!(
                fields
                    .iter()
                    .map(|field| field.value().attr("data-field").unwrap())
                    .collect::<Vec<_>>(),
                header
            );
            fields
                .iter()
                .map(|field| field.text().collect::<String>())
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(html_rows, rows);
    let schema = Selector::parse("meta[name=viroflash-report-schema]").unwrap();
    assert_eq!(
        document
            .select(&schema)
            .next()
            .unwrap()
            .value()
            .attr("content"),
        Some("viroflash.evidence-report.v1")
    );
    assert!(document
        .select(&Selector::parse("script").unwrap())
        .next()
        .is_none());
    let download = document
        .select(&Selector::parse("a[download]").unwrap())
        .next()
        .unwrap();
    let encoded = download
        .value()
        .attr("href")
        .unwrap()
        .strip_prefix("data:text/csv;charset=utf-8,")
        .unwrap();
    let embedded = encoded
        .as_bytes()
        .chunks_exact(3)
        .map(|chunk| {
            assert_eq!(chunk[0], b'%');
            u8::from_str_radix(std::str::from_utf8(&chunk[1..]).unwrap(), 16).unwrap()
        })
        .collect::<Vec<_>>();
    assert_eq!(embedded, csv.as_bytes());
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
    let output = command(&[
        "run",
        "--r1",
        root.join("sample.fastq.gz").to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out").to_str().unwrap(),
        "--threads",
        "1",
    ]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output_names(&root.join("out")),
        BTreeSet::from([
            "perf.json".into(),
            "report.csv".into(),
            "report.html".into()
        ])
    );
    let perf: serde_json::Value =
        serde_json::from_slice(&std::fs::read(root.join("out").join("perf.json")).unwrap())
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
    let csv = std::fs::read(root.join("out/report.csv")).unwrap();
    let text = String::from_utf8(csv).unwrap();
    let (header, rows) = csv_rows(&text);
    assert_eq!(rows.len(), 1);
    let html = std::fs::read_to_string(root.join("out/report.html")).unwrap();
    assert_csv_and_visible_html_share_all_fields(&text, &html);
    let value = |field: &str| {
        rows[0][header
            .iter()
            .position(|candidate| candidate == field)
            .unwrap()]
        .as_str()
    };
    assert_eq!(value("sample_id"), "sample");
    assert_eq!(value("support_fragments"), value("selected_fragments"));
    assert!((1..=20).contains(&value("selected_fragments").parse::<u64>().unwrap()));
    assert_eq!(value("support_ppm"), "1000000");
    assert_eq!(value("target_support_share_pct"), "100");

    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn report_is_invariant_across_threads_reloads_and_gzip_segmentation() {
    let (root, target) = build_fixture("invariance");
    // Cross several input batches so thread invariance exercises parallel dispatch.
    let mut fastq = Vec::new();
    for ordinal in 0..4099 {
        // Repeated IDs and non-candidates produce sparse original ordinals.
        let read = if ordinal % 3 == 0 {
            b"AAAAAAAAAAAAAAAAAAAA".as_slice()
        } else {
            &target[..120]
        };
        fastq.extend(fastq_bytes(read, 1));
    }
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

    // Vary workers on identical input, then encoding at a fixed worker count.
    let mut expected = None;
    for (case, input, threads) in [
        ("plain-1", &plain, 1),
        ("plain-2", &plain, 2),
        ("plain-4", &plain, 4),
        ("plain-8", &plain, 8),
        ("gzip", &single, 4),
        ("multimember", &multiple, 4),
    ] {
        let output_dir = root.join(case);
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
        let report = (
            std::fs::read(output_dir.join("report.csv")).unwrap(),
            std::fs::read(output_dir.join("report.html")).unwrap(),
        );
        if let Some(expected) = &expected {
            assert_eq!(&report, expected, "{case}");
        } else {
            expected = Some(report);
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
    let (header, rows) = csv_rows(&csv);
    let value = |field: &str| {
        rows[0][header
            .iter()
            .position(|candidate| candidate == field)
            .unwrap()]
        .as_str()
    };
    assert_eq!(value("input_mode"), "PE");
    assert_eq!(value("input_fragments"), "20");
    assert_eq!(value("support_fragments"), value("selected_fragments"));
    assert_eq!(
        value("total_target_support_fragments"),
        value("selected_fragments")
    );
    assert!((1..=20).contains(&value("selected_fragments").parse::<u64>().unwrap()));
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn research_report_ranks_groups_and_preserves_indexed_descriptions() {
    let root = workspace("research-report");
    let host = sequence(17, 600);
    let alpha = sequence(91, 600);
    let zeta = sequence(321, 600);
    write_fasta(&root.join("host.fa"), "host", &host);
    let description = "Virus \"alpha\", <script>window.injected=true</script> {{AUDIT_SIGNALS}}";
    let fasta = format!(
        ">zeta Zeta virus\n{}\n>alpha |{description}\n{}\n>alpha-copy Equivalent alpha\n{}\n",
        String::from_utf8_lossy(&zeta),
        String::from_utf8_lossy(&alpha),
        String::from_utf8_lossy(&alpha)
    );
    std::fs::write(root.join("target.fa"), fasta).unwrap();
    let built = command(&[
        "index",
        "--host-fa",
        root.join("host.fa").to_str().unwrap(),
        "--target-fa",
        root.join("target.fa").to_str().unwrap(),
        "--out",
        root.join("index").to_str().unwrap(),
    ]);
    assert!(
        built.status.success(),
        "{}",
        String::from_utf8_lossy(&built.stderr)
    );
    // Reusable indexes must carry descriptions without consulting the original files.
    std::fs::remove_file(root.join("target.fa")).unwrap();
    std::fs::remove_file(root.join("host.fa")).unwrap();
    let fastq = [
        fastq_bytes(&zeta[..120], 3),
        fastq_bytes(&alpha[..120], 6),
        fastq_bytes(&host[..120], 2),
    ]
    .concat();
    let input = root.join("sample_{{RUN_FIELDS}}.fastq");
    std::fs::write(&input, fastq).unwrap();
    let result = command(&[
        "run",
        "--r1",
        input.to_str().unwrap(),
        "--index",
        root.join("index").to_str().unwrap(),
        "--out",
        root.join("out").to_str().unwrap(),
        "--threads",
        "3",
    ]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let csv = std::fs::read_to_string(root.join("out/report.csv")).unwrap();
    let html = std::fs::read_to_string(root.join("out/report.html")).unwrap();
    assert_csv_and_visible_html_share_all_fields(&csv, &html);
    let (header, rows) = csv_rows(&csv);
    let value = |row: usize, field: &str| {
        rows[row][header.iter().position(|key| key == field).unwrap()].as_str()
    };
    assert_eq!(rows.len(), 2);
    assert_eq!(value(0, "reference_id"), "alpha");
    assert_eq!(value(1, "reference_id"), "zeta");
    assert_eq!(value(0, "reference_description"), description);
    assert_eq!(value(0, "reference_member_ids"), "alpha;alpha-copy");
    assert_eq!(value(0, "reference_member_count"), "2");
    assert_eq!(value(0, "support_fragments"), "6");
    assert_eq!(value(1, "support_fragments"), "3");
    assert_eq!(value(0, "support_rank"), "1");
    assert_eq!(value(1, "support_rank"), "2");
    assert_eq!(value(0, "input_fragments"), "11");
    assert_eq!(value(0, "selected_fragments"), "9");
    assert_eq!(value(0, "total_target_support_fragments"), "9");
    assert_eq!(value(0, "supported_reference_groups"), "2");
    assert_eq!(value(0, "target_support_share_pct"), "66.666667");
    assert_eq!(value(1, "target_support_share_pct"), "33.333333");
    assert_eq!(value(0, "support_ppm"), "545454.545455");
    // Every candidate was selected: six supports among eleven original fragments.
    let lower: f64 = value(0, "support_ci_lower_ppm").parse().unwrap();
    let upper: f64 = value(0, "support_ci_upper_ppm").parse().unwrap();
    assert!(lower <= 6.0 / 11.0 * 1_000_000.0 && upper >= 6.0 / 11.0 * 1_000_000.0);
    assert!(upper - lower < 0.00001);
    assert!(html.contains("sample_{{RUN_FIELDS}}"));
    assert!(html.contains("{{AUDIT_SIGNALS}}"));
    assert!(!html.contains("<script>window.injected"));
    assert!(html.contains("&lt;script&gt;window.injected"));
    assert!(html.contains("lang=\"en\""));
    let _ = std::fs::remove_dir_all(root);
}

#[test]
fn invalid_index_metadata_is_rejected_without_success_reports() {
    let (root, target) = build_fixture("invalid-index");
    write_fastq(&root.join("sample.fastq.gz"), &target, None);
    let path = root.join("index/manifest.json");
    let original: serde_json::Value =
        serde_json::from_slice(&std::fs::read(&path).unwrap()).unwrap();
    for (case, message) in [
        ("missing-description", "rebuild"),
        ("duplicate-group-owner", "rebuild"),
        ("non-dense-contigs", "dense"),
    ] {
        let mut manifest = original.clone();
        match case {
            "missing-description" => {
                manifest
                    .as_object_mut()
                    .unwrap()
                    .remove("target_descriptions");
            }
            "duplicate-group-owner" => {
                manifest["target_groups"] = serde_json::json!([{"ordinal": 0}])
            }
            "non-dense-contigs" => {
                manifest["contigs"]
                    .as_array_mut()
                    .unwrap()
                    .iter_mut()
                    .find(|contig| contig["role"] == "Target")
                    .unwrap()["target_group_ordinal"] = serde_json::json!(1);
            }
            _ => unreachable!(),
        }
        std::fs::write(&path, serde_json::to_vec(&manifest).unwrap()).unwrap();
        let result = command(&[
            "run",
            "--r1",
            root.join("sample.fastq.gz").to_str().unwrap(),
            "--index",
            root.join("index").to_str().unwrap(),
            "--out",
            root.join(case).to_str().unwrap(),
        ]);
        assert!(!result.status.success(), "{case}");
        assert!(
            String::from_utf8_lossy(&result.stderr).contains(message),
            "{case}"
        );
        assert_eq!(
            output_names(&root.join(case)),
            BTreeSet::from(["perf.json".into()])
        );
        let perf: serde_json::Value =
            serde_json::from_slice(&std::fs::read(root.join(case).join("perf.json")).unwrap())
                .unwrap();
        assert_eq!(perf["status"], "ERROR");
    }
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
fn zero_candidates_report_full_input_prescreen_limitation() {
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
        .as_str()
    };
    assert_eq!(value("reference_id"), "");
    assert_eq!(value("support_fragments"), "");
    assert_eq!(value("total_target_support_fragments"), "0");
    assert_eq!(value("supported_reference_groups"), "0");
    let selected = value("selected_fragments").parse::<u64>().unwrap();
    assert_eq!(selected, 0);
    assert!(html.contains("CONFORMANT_WITH_LIMITATIONS"));
    assert!(html.contains("TARGET_KMER_NOT_EVALUABLE_INPUT_FRAGMENTS=100"));
    let _ = std::fs::remove_dir_all(root);
}
