//! 合成数据端到端冒烟测试：随机宿主/目标/诱饵/污染参考 + 合成读对。
//! 覆盖：
//! - 「--index 加载」与「FASTA 自动构建」两条路径结果逐候选一致（等价性契约）；
//! - 目标候选通过精确率检验 adjusted-p 与分布门（输出文件生成）；
//! - 未提供 --decoy-fa 时自动生成诱饵并写入索引；
//! - 索引 k 不一致、--index 与 FASTA 互斥的错误路径。

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use flate2::write::GzEncoder;
use flate2::Compression;

use viroflash::index::{self, IndexOptions};
use viroflash::{run_pipeline, Options};

/// xorshift64 确定性随机源。
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

    /// 精确 50% GC 的随机序列：一半 AT、一半 CG 后洗牌，
    /// 保证目标与全部诱饵落入同一 gc 分层（gc:0.50-0.60）。
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
        // Fisher-Yates 洗牌（用同一 rng）
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

/// 写四类参考与合成读对，返回工作目录。
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

    // 30 个目标对（来自 target 三段分布区域）+ 30 个宿主对。三段使该夹具明确
    // 覆盖通用检测的 distributed-windows PASS 路径，而不依赖 split 证据。
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

/// 基础断言：目标候选被检出、背景受控、输出文件生成。
fn assert_detected(summary: &viroflash::RunSummary) {
    assert_eq!(summary.input_pairs, 60);
    assert!(
        summary.prescreen_pairs >= 30,
        "预筛至少保留全部目标对: {}",
        summary.prescreen_pairs
    );

    let cand = summary
        .candidates
        .iter()
        .find(|c| c.contig == "target_0")
        .expect("应检出 target_0 候选");
    assert_eq!(summary.test_family_size, 1);
    assert!(
        cand.covered_frac >= 0.2,
        "覆盖比例过低: {:.3}",
        cand.covered_frac
    );
    assert!(
        cand.q_value <= 0.2,
        "模型 adjusted-p 应通过探索阈值（同层 decoy 无覆盖）: {:.4}",
        cand.q_value
    );
    assert_eq!(cand.decision, "PASS");
    assert!(cand.distinct_windows >= 3);
    assert_eq!(cand.confidence(), "UNVALIDATED");
    assert_eq!(cand.n_plain, cand.reads, "split 不得改变通用检测计数");
    assert_eq!(
        cand.background_status, "SYNTHETIC_DECOY_UNVALIDATED",
        "报告必须披露 synthetic-decoy null 未校准"
    );
    // 合成读无嵌合，不应有 split/discordant 证据
    assert_eq!(cand.split_events, 0);
    assert_eq!(cand.discordant, 0);
    assert_eq!(cand.integration_evidence, "NONE");

    assert!(summary.result_json.exists());
    assert!(summary.result_tsv.exists());
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
    for line in result_tsv.lines() {
        assert_eq!(line.split('\t').count(), 22, "TSV 固定列契约: {line}");
    }
    let perf_json = summary.result_json.with_extension("perf.json");
    let perf_tsv = summary.result_tsv.with_extension("perf.tsv");
    assert!(perf_json.exists());
    assert!(perf_tsv.exists());
    let perf = std::fs::read_to_string(perf_json).unwrap();
    assert!(perf.contains("\"schema\": \"viroflash.perf.v1\""));
    assert!(perf.contains("\"status\": \"success\""));
    assert!(!perf.contains("reads_R1.fq.gz"), "性能报告不得泄露输入路径");
}

#[test]
fn e2e_index_and_autobuild_paths_agree() {
    let dir = setup_synthetic("equiv");

    // 路径 A：先构建索引目录，再 --index 加载运行。
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
    assert!(dir.join("idx.perf.tsv").is_file());
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

    // 路径 B：直接给 FASTA，自动构建（与 index 命令共用同一构建函数）。
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
    // 自动构建的索引落于 <out>.work/index/
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

    // 未提供 --decoy-fa：索引构建阶段按默认 ANI/seed 从目标自动生成诱饵。
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

    // 自动诱饵产物与 manifest 的 generated 来源记录
    let idx = dir.join("out.work").join("index");
    assert!(idx.join("decoys.fa").is_file());
    assert!(idx.join("decoys.tsv").is_file());
    let manifest = std::fs::read_to_string(idx.join("manifest.json")).unwrap();
    assert!(
        manifest.contains("\"generated\""),
        "manifest 应记录自动生成诱饵: {manifest}"
    );
    assert!(manifest.contains("\"anis\""), "manifest 应记录 ANI 层");

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

    // k 不一致：加载必须拒绝（Bloom 与索引均按 k 构建）。
    let err = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        index: Some(dir.join("idx")),
        out: dir.join("k_mismatch"),
        k: 6,
        ..Options::default()
    })
    .unwrap_err();
    assert!(err.contains("不一致"), "err={err}");

    // --index 与 FASTA 互斥（库层防御，与 CLI 解析双重校验）。
    let err = run_pipeline(&Options {
        r1: dir.join("reads_R1.fq.gz"),
        index: Some(dir.join("idx")),
        host_fa: Some(dir.join("host.fa")),
        out: dir.join("index_fasta_conflict"),
        ..Options::default()
    })
    .unwrap_err();
    assert!(err.contains("不能同时使用"), "err={err}");

    // 自动构建缺失必填参考。
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
