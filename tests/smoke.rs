//! 合成数据端到端冒烟测试：随机宿主/目标/诱饵/污染参考 + 合成读对，
//! 验证「目标候选被检出且 q 值显著」。
//! 测试同时检查背景控制与输出文件生成。

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;

use flate2::write::GzEncoder;
use flate2::Compression;

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

fn synthetic_workspace() -> PathBuf {
    let dir = std::env::temp_dir().join(format!("viroflash_smoke_{}", std::process::id()));
    if dir.exists() {
        let _ = std::fs::remove_dir_all(&dir);
    }
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

#[test]
fn e2e_detects_target_and_controls_background() {
    let dir = synthetic_workspace();
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

    // 30 个目标对（来自 target 两段区域）+ 30 个宿主对
    let mut r1: Vec<(String, Vec<u8>)> = Vec::new();
    let mut r2: Vec<(String, Vec<u8>)> = Vec::new();
    for i in 0..30 {
        r1.push((format!("t{i}/1"), target[100..250].to_vec()));
        r2.push((format!("t{i}/2"), reverse_complement(&target[400..550])));
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

    let opt = Options {
        r1: dir.join("reads_R1.fq.gz"),
        r2: Some(dir.join("reads_R2.fq.gz")),
        host_fa: dir.join("host.fa"),
        target_fa: dir.join("target.fa"),
        decoy_fa: dir.join("decoy.fa"),
        contam_fa: dir.join("contam.fa"),
        threads: 4,
        out: dir.join("out"),
        k: 21,
    };

    let summary = run_pipeline(&opt).unwrap();

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
    assert!(
        cand.covered_frac >= 0.2,
        "覆盖比例过低: {:.3}",
        cand.covered_frac
    );
    assert!(
        cand.q_value <= 0.2,
        "q 值应显著（同层 decoy 无覆盖）: {:.4}",
        cand.q_value
    );
    assert_eq!(cand.stratum_decoy_count, 20);
    // 合成读无嵌合，不应有 split/discordant 证据
    assert_eq!(cand.split_events, 0);
    assert_eq!(cand.discordant, 0);

    assert!(summary.result_json.exists());
    assert!(summary.result_tsv.exists());

    let _ = std::fs::remove_dir_all(&dir);
}
