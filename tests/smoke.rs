//! End-to-end smoke tests using random host, target, decoy, and contaminant references plus
//! synthetic read pairs. Covers reusable-index/autobuild equivalence, target reporting through the
//! exact-rate and distribution gates, automatic decoys, k mismatch, and conflicting input modes.

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use flate2::write::GzEncoder;
use flate2::Compression;

use viroflash::index::{self, IndexOptions};
use viroflash::{run_pipeline, Options};

fn parse_csv_record(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let mut field = String::new();
    let mut chars = line.chars().peekable();
    let mut quoted = false;
    while let Some(ch) = chars.next() {
        match ch {
            '"' if quoted && chars.peek() == Some(&'"') => {
                field.push('"');
                chars.next();
            }
            '"' => quoted = !quoted,
            ',' if !quoted => {
                fields.push(std::mem::take(&mut field));
            }
            _ => field.push(ch),
        }
    }
    assert!(!quoted, "unterminated quoted CSV field: {line}");
    fields.push(field);
    fields
}

fn assert_only_perf_json(out_prefix: &std::path::Path) {
    let parent = out_prefix.parent().unwrap();
    let label = out_prefix.file_name().unwrap().to_string_lossy();
    let mut reports: Vec<String> = std::fs::read_dir(parent)
        .unwrap()
        .map(|entry| entry.unwrap().file_name().to_string_lossy().into_owned())
        .filter(|name| name.starts_with(&format!("{label}.perf.")))
        .collect();
    reports.sort();
    assert_eq!(reports, vec![format!("{label}.perf.json")]);
}

/// Deterministic xorshift64 random source.
struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    fn seq(&mut self, n: usize) -> Vec<u8> {
        (0..n)
            .map(|_| b"ACGT"[(self.next() % 4) as usize])
            .collect()
    }

    /// Random sequence with exactly 50% GC, shuffled after generating equal AT and CG halves.
    /// This keeps the target and every decoy in the same GC stratum.
    fn half_gc_seq(&mut self, n: usize) -> Vec<u8> {
        let mut half: Vec<u8> = (0..n)
            .map(|i| {
                if i % 2 == 0 {
                    b"AT"[(self.next() % 2) as usize]
                } else {
                    b"CG"[(self.next() % 2) as usize]
                }
            })
            .collect();
        // Fisher-Yates shuffle using the same RNG.
        for i in (1..half.len()).rev() {
            let j = (self.next() % (i as u64 + 1)) as usize;
            half.swap(i, j);
        }
        half
    }
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            other => other,
        })
        .collect()
}

fn write_fasta(path: &std::path::Path, records: &[(&str, &[u8])]) {
    let mut w = BufWriter::new(File::create(path).unwrap());
    for (header, seq) in records {
        writeln!(w, ">{header}").unwrap();
        w.write_all(seq).unwrap();
        writeln!(w).unwrap();
    }
}

fn write_fastq_gz(path: &std::path::Path, records: &[(&str, &[u8])]) {
    let file = File::create(path).unwrap();
    let mut w = GzEncoder::new(file, Compression::default());
    for (id, seq) in records {
        writeln!(w, "@{id}").unwrap();
        w.write_all(seq).unwrap();
        writeln!(w, "\n+").unwrap();
        w.write_all(&vec![b'I'; seq.len()]).unwrap();
        writeln!(w).unwrap();
    }
    w.finish().unwrap();
}

fn synthetic_workspace(tag: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("viroflash_smoke_{tag}_{}", std::process::id()));
    if dir.exists() {
        let _ = std::fs::remove_dir_all(&dir);
    }
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

/// Write all four reference roles and synthetic read pairs, returning the work directory.
fn setup_synthetic(tag: &str) -> PathBuf {
    let dir = synthetic_workspace(tag);
    let mut rng = Rng(0x9e3779b97f4a7c15);

    let host = rng.seq(5000);
    let target = rng.half_gc_seq(1000);
    let decoys: Vec<Vec<u8>> = (0..20).map(|_| rng.half_gc_seq(1000)).collect();
    let contam = rng.seq(1000);

    write_fasta(&dir.join("host.fa"), &[("chrH", &host)]);
    write_fasta(&dir.join("target.fa"), &[("TESTVIR", &target)]);
    let decoy_records: Vec<(String, Vec<u8>)> = decoys
        .iter()
        .enumerate()
        .map(|(i, d)| (format!("dec{i}"), d.clone()))
        .collect();
    write_fasta(
        &dir.join("decoy.fa"),
        &decoy_records
            .iter()
            .map(|(h, s)| (h.as_str(), s.as_slice()))
            .collect::<Vec<_>>(),
    );
    write_fasta(&dir.join("contam.fa"), &[("mycoplasma", &contam)]);

    // Thirty target pairs across three distributed regions plus thirty host pairs exercise the
    // distributed-window PASS path without relying on split evidence.
    let mut r1: Vec<(String, Vec<u8>)> = Vec::new();
    let mut r2: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..30 {
        r1.push((format!("t{i}/1"), target[100..250].to_vec()));
        let r2_start = if i % 2 == 0 { 400 } else { 700 };
        r2.push((
            format!("t{i}/2"),
            reverse_complement(&target[r2_start..r2_start + 150]),
        ));
    }
    for i in 0..30 {
        let off = i * 10;
        r1.push((format!("h{i}/1"), host[off..off + 150].to_vec()));
        r2.push((
            format!("h{i}/2"),
            reverse_complement(&host[1000 + off..1000 + off + 150]),
        ));
    }
    write_fastq_gz(
        &dir.join("reads_R1.fq.gz"),
        &r1.iter()
            .map(|(h, s)| (h.as_str(), s.as_slice()))
            .collect::<Vec<_>>(),
    );
    write_fastq_gz(
        &dir.join("reads_R2.fq.gz"),
        &r2.iter()
            .map(|(h, s)| (h.as_str(), s.as_slice()))
            .collect::<Vec<_>>(),
    );
    dir
}

/// Assert target reporting, controlled background, and generated output files.
fn assert_detected(summary: &viroflash::RunSummary) {
    assert_eq!(summary.input_pairs, 60);
    assert!(
        summary.prescreen_pairs >= 30,
        "prescreen must retain at least all target pairs: {}",
        summary.prescreen_pairs
    );

    let cand = summary
        .candidates
        .iter()
        .find(|c| c.contig == "target_0")
        .expect("target_0 candidate should be reported");
    assert_eq!(summary.test_family_size, 1);
    assert!(
        cand.covered_frac >= 0.2,
        "breadth is too low: {:.3}",
        cand.covered_frac
    );
    assert!(
        cand.q_value <= 0.2,
        "model adjusted p-value should pass the exploratory threshold with no same-stratum decoy coverage: {:.4}",
        cand.q_value
    );
    assert_eq!(cand.decision, "PASS");
    assert!(cand.distinct_windows >= 3);
    assert_eq!(cand.confidence(), "UNVALIDATED");
    assert_eq!(
        cand.n_plain, cand.reads,
        "split evidence must not change general detection counts"
    );
    assert_eq!(
        cand.background_status, "SYNTHETIC_DECOY_UNVALIDATED",
        "the report must disclose that the synthetic-decoy null is uncalibrated"
    );
    // Synthetic reads are not chimeric and should have no split or discordant evidence.
    assert_eq!(cand.split_events, 0);
    assert_eq!(cand.discordant, 0);
    assert_eq!(cand.integration_evidence, "NONE");

    assert!(summary.result_json.exists());
    assert!(summary.result_tsv.exists());
    assert!(summary.result_html.exists());
    assert!(summary.result_csv.exists());
    let result_json = std::fs::read_to_string(&summary.result_json).unwrap();
    let result_tsv = std::fs::read_to_string(&summary.result_tsv).unwrap();
    assert!(result_json.contains("\"schema\": \"viroflash.result.v1\""));
    assert!(result_json.contains("\"sample_conclusion\": \"not_computed\""));
    assert!(result_json.contains("\"fdr_control_validated\": false"));
    assert!(result_json.contains("\"test\": \"exact_conditional_two_poisson_rates\""));
    assert!(result_tsv.contains("row_type=candidate;"));
    assert!(result_tsv.contains("result_schema=viroflash.result.v1"));
    assert!(result_tsv.contains("\treference_group\t"));
    assert!(result_tsv.contains("qc_status=NOT_EVALUATED"));
    assert!(result_tsv.contains("member_attribution=not_resolved"));
    let result_html = std::fs::read_to_string(&summary.result_html).unwrap();
    let result_csv = std::fs::read_to_string(&summary.result_csv).unwrap();
    assert!(result_html.contains("viroflash evidence report"));
    assert!(result_html.contains("Interpretation boundary"));
    assert!(result_html.contains("FDR not validated"));
    assert!(result_html.contains("data-decision=\"PASS\""));
    assert!(result_html.contains("Download core CSV"));
    assert!(result_csv
        .starts_with("csv_schema,result_schema,sample_id,sample_conclusion,qc_status,qc_issues"));
    assert!(result_csv.contains("viroflash.candidates.csv.v1"));
    assert!(result_csv.contains("not_computed,NOT_EVALUATED"));
    assert!(result_csv.contains("model_adjusted_p_max"));
    assert!(result_csv.contains("coverage_min"));
    assert!(result_csv.contains("min_distributed_windows"));
    assert!(result_csv.contains("manifest_blake3"));
    let json: serde_json::Value = serde_json::from_str(&result_json).unwrap();
    let mut csv_lines = result_csv.lines();
    let csv_header = parse_csv_record(csv_lines.next().unwrap());
    let csv_record = parse_csv_record(csv_lines.next().unwrap());
    assert_eq!(csv_header.len(), csv_record.len());
    assert!(csv_lines.next().is_none());
    let csv: std::collections::HashMap<_, _> = csv_header
        .iter()
        .zip(csv_record.iter())
        .map(|(key, value)| (key.as_str(), value.as_str()))
        .collect();
    let json_candidate = &json["candidates"][0];
    assert_eq!(csv["member_attribution"], "not_resolved");
    assert_eq!(csv["sample_id"], json["run"]["sample"].as_str().unwrap());
    assert_eq!(csv["qc_status"], json["quality_control"]["status"]);
    assert_eq!(
        serde_json::from_str::<serde_json::Value>(csv["qc_issues"]).unwrap(),
        json["quality_control"]["issues"]
    );
    assert_eq!(
        csv["input_pairs"].parse::<u64>().unwrap(),
        json["run"]["input_pairs"].as_u64().unwrap()
    );
    assert_eq!(csv["candidate_id"], json_candidate["contig"]);
    assert_eq!(csv["decision"], json_candidate["decision"]);
    assert_eq!(
        csv["validation_read_ends"].parse::<u64>().unwrap(),
        json_candidate["evidence"]["reads"].as_u64().unwrap()
    );
    assert_eq!(
        csv["coverage_breadth"].parse::<f64>().unwrap(),
        json_candidate["evidence"]["covered_frac"].as_f64().unwrap()
    );
    assert_eq!(
        csv["manifest_blake3"],
        json["index"]["manifest_blake3"].as_str().unwrap()
    );
    for line in result_tsv.lines() {
        assert_eq!(
            line.split('\t').count(),
            22,
            "fixed 22-column TSV contract: {line}"
        );
    }
    let perf_json = summary.result_json.with_extension("perf.json");
    assert!(perf_json.exists());
    assert_only_perf_json(&summary.result_json.with_extension(""));
    let perf = std::fs::read_to_string(perf_json).unwrap();
    assert!(perf.contains("\"schema\": \"viroflash.perf.v1\""));
    assert!(perf.contains("\"status\": \"success\""));
    assert!(
        !perf.contains("reads_R1.fq.gz"),
        "performance reports must not expose input paths"
    );
}

#[test]
fn e2e_index_and_autobuild_paths_agree() {
    let dir = setup_synthetic("equiv");

    // Path A: build an index directory, then load it with --index.
    let index_opts = IndexOptions {
        host_fa: dir.join("host.fa"),
        target_fa: dir.join("target.fa"),
        contam_fa: Some(dir.join("contam.fa")),
        decoy_fa: Some(dir.join("decoy.fa")),
        out_dir: dir.join("idx"),
        k: 21,
        threads: 4,
        ..IndexOptions::default()
    };
    let built = index::build_index(&index_opts).unwrap();
    assert!(dir.join("idx.perf.json").is_file());
    assert_only_perf_json(&dir.join("idx"));
    assert_eq!(
        built
            .contigs
            .iter()
            .filter(|c| c.role == viroflash::reference::Role::Decoy)
            .count(),
        20
    );
    let summary_loaded = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        r2: Some(dir.join("reads_R2.fq.gz")),
        index: Some(dir.join("idx")),
        threads: 4,
        out: dir.join("out_idx"),
        k: 21,
        ..Options::default()
    })
    .unwrap();

    // Path B: pass FASTA files directly for automatic construction through the shared builder.
    let summary_built = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        r2: Some(dir.join("reads_R2.fq.gz")),
        host_fa: Some(dir.join("host.fa")),
        target_fa: Some(dir.join("target.fa")),
        decoy_fa: Some(dir.join("decoy.fa")),
        contam_fa: Some(dir.join("contam.fa")),
        threads: 4,
        out: dir.join("out_fa"),
        k: 21,
        ..Options::default()
    })
    .unwrap();

    assert_detected(&summary_loaded);
    assert_detected(&summary_built);
    assert_eq!(summary_loaded.candidates, summary_built.candidates);
    // Automatic construction writes the index under <out>.work/index/.
    assert!(dir
        .join("out_fa.work")
        .join("index")
        .join("manifest.json")
        .is_file());

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn e2e_auto_decoy_when_decoy_fa_absent() {
    let dir = setup_synthetic("autodecoy");

    // Without --decoy-fa, index construction generates decoys with default ANI and seed.
    let summary = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        r2: Some(dir.join("reads_R2.fq.gz")),
        host_fa: Some(dir.join("host.fa")),
        target_fa: Some(dir.join("target.fa")),
        threads: 4,
        out: dir.join("out"),
        k: 21,
        ..Options::default()
    })
    .unwrap();
    assert_detected(&summary);

    // Generated decoy artifacts and their manifest provenance.
    let idx = dir.join("out.work").join("index");
    assert!(idx.join("decoys.fa").is_file());
    assert!(idx.join("decoys.tsv").is_file());
    let manifest = std::fs::read_to_string(idx.join("manifest.json")).unwrap();
    assert!(
        manifest.contains("\"generated\""),
        "manifest should record generated decoys: {manifest}"
    );
    assert!(
        manifest.contains("\"anis\""),
        "manifest should record ANI layers"
    );

    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn e2e_rejects_k_mismatch_and_index_fasta_conflict() {
    let dir = setup_synthetic("errors");

    let index_opts = IndexOptions {
        host_fa: dir.join("host.fa"),
        target_fa: dir.join("target.fa"),
        decoy_fa: Some(dir.join("decoy.fa")),
        out_dir: dir.join("idx"),
        k: 21,
        threads: 4,
        ..IndexOptions::default()
    };
    index::build_index(&index_opts).unwrap();

    // Loading must reject a k mismatch because both Bloom and MMI are built for one k.
    let err = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        index: Some(dir.join("idx")),
        out: dir.join("k_mismatch"),
        k: 6,
        ..Options::default()
    })
    .unwrap_err();
    assert!(err.contains("does not match"), "err={err}");

    // --index and FASTA inputs are mutually exclusive at both library and CLI boundaries.
    let err = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        index: Some(dir.join("idx")),
        host_fa: Some(dir.join("host.fa")),
        out: dir.join("index_fasta_conflict"),
        ..Options::default()
    })
    .unwrap_err();
    assert!(err.contains("cannot be combined"), "err={err}");

    // Automatic construction rejects missing required references.
    let err = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        target_fa: Some(dir.join("target.fa")),
        out: dir.join("missing_host"),
        ..Options::default()
    })
    .unwrap_err();
    assert!(err.contains("--host-fa"), "err={err}");

    let _ = std::fs::remove_dir_all(&dir);
}
