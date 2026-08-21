//! 病毒候选检测管线编排：
//! 参考构建 → 预筛+竞争比对融合（单趟并行，k-mer 门在工作线程内）→ 候选聚类 → decoy 统计判定 → 报告。

pub mod align;
pub mod cluster;
pub mod decoy;
pub mod fastq;
pub mod prescreen;
pub mod reference;
pub mod report;
pub mod stats;

use std::collections::{HashMap, HashSet};
use std::path::{Path, PathBuf};

use reference::Role;

#[derive(Debug, Clone)]
pub struct Options {
    pub r1: PathBuf,
    /// None = 单端模式（454/单端测序），每 read 一个 fragment。
    pub r2: Option<PathBuf>,
    pub host_fa: PathBuf,
    pub target_fa: PathBuf,
    pub decoy_fa: PathBuf,
    pub contam_fa: PathBuf,
    pub threads: usize,
    pub out: PathBuf,
    pub k: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            r1: PathBuf::new(),
            r2: None,
            host_fa: PathBuf::new(),
            target_fa: PathBuf::new(),
            decoy_fa: PathBuf::new(),
            contam_fa: PathBuf::new(),
            threads: 8,
            out: PathBuf::from("viroflash_out"),
            k: prescreen::DEFAULT_K,
        }
    }
}

/// 证据层组合门槛（Mourik 2024 J Clin Microbiol 62(6):e00345-24，
/// doi:10.1128/jcm.00345-24）：q<0.2 候选须同时满足 500 RPM 深度与 10% 覆盖
/// 两门槛才进入 PASS/WEAK，否则 BELOW_THRESHOLD。RPM = reads per million input
/// fragments（PE pair = 1 fragment，SE read = 1 fragment，分母即 input_pairs）。
/// 门槛为 AND 判定，不乘入 p/q。
pub const DEPTH_RPM_MIN: f64 = 500.0;
pub const COVERAGE_MIN: f64 = 0.10;

#[derive(Debug, Clone)]
pub struct RunSummary {
    pub input_pairs: u64,
    pub prescreen_pairs: u64,
    /// 比对失败的 fragment 数（minimap2 返回错误；不中止运行但显式披露）。
    pub map_errors: u64,
    pub candidates: Vec<report::Candidate>,
    pub result_json: PathBuf,
    pub result_tsv: PathBuf,
}

pub fn run_pipeline(opt: &Options) -> Result<RunSummary, String> {
    if opt.threads == 0 {
        return Err("--threads 必须大于 0".into());
    }
    if !(1..=prescreen::K_MAX).contains(&opt.k) {
        return Err(format!(
            "--k 必须在 1..={} 之间（2-bit 编码上限），得到 {}",
            prescreen::K_MAX,
            opt.k
        ));
    }
    for (label, p) in [
        ("r1", &opt.r1),
        ("host-fa", &opt.host_fa),
        ("target-fa", &opt.target_fa),
        ("decoy-fa", &opt.decoy_fa),
        ("contam-fa", &opt.contam_fa),
    ] {
        if p.as_os_str().is_empty() {
            return Err(format!("缺少 --{label}"));
        }
    }

    let work_dir = PathBuf::from(format!("{}.work", opt.out.display()));
    let fastas = vec![
        (Role::Host, opt.host_fa.clone()),
        (Role::Target, opt.target_fa.clone()),
        (Role::Decoy, opt.decoy_fa.clone()),
        (Role::Contaminant, opt.contam_fa.clone()),
    ];
    let mut t0 = std::time::Instant::now();
    let (mmi_path, contigs) = reference::build_reference(&fastas, &work_dir, opt.threads)?;
    eprintln!("[阶段] 参考构建+索引 {:.1}s", t0.elapsed().as_secs_f64());
    let roles: HashMap<String, Role> = contigs.iter().map(|c| (c.name.clone(), c.role)).collect();

    // 1. 预筛（宽松保留）
    let target_refs: Vec<&reference::Contig> =
        contigs.iter().filter(|c| c.role == Role::Target).collect();
    if target_refs.is_empty() {
        return Err("目标 FASTA 中没有序列".into());
    }
    // Bloom 覆盖目标与诱饵的 k-mer 并集，使 82–88% ANI 的诱饵 reads 能通过
    // 比例门进入竞争比对；否则诱饵层 reads 恒为 0，导致 λ̂=0、p=0。
    let bloom_refs: Vec<&reference::Contig> = contigs
        .iter()
        .filter(|c| matches!(c.role, Role::Target | Role::Decoy))
        .collect();
    let bloom =
        prescreen::KmerBloom::build(&bloom_refs, opt.k, prescreen::GATE_FPR).ok_or_else(|| {
            format!(
                "目标+诱饵 FASTA 中没有 ≥k={} 的有效 k-mer（全部序列过短或含非 ACGT 字符）",
                opt.k
            )
        })?;
    // N 总线程 = 1 主 + D 解压 + W worker，D:W 约为 1:3。
    let (decomp_threads, workers) = thread_budget(opt.threads);
    let aligner = align::CompetitiveAligner::open(&mmi_path, workers)?;
    t0 = std::time::Instant::now();
    let (input_pairs, prescreen_pairs, map_errors, evidences) = prescreen_and_align(
        &opt.r1,
        opt.r2.as_deref(),
        opt.k,
        &bloom,
        &aligner,
        &roles,
        decomp_threads,
    )?;
    eprintln!("[阶段] 预筛+竞争比对 {:.1}s", t0.elapsed().as_secs_f64());

    // 3. 聚合 + 候选聚类
    t0 = std::time::Instant::now();
    let aggs = aggregate_evidence(&contigs, &evidences, opt.threads)?;
    eprintln!("[阶段] 证据聚合 {:.1}s", t0.elapsed().as_secs_f64());

    // 4. decoy 统计判定（N_total = 输入 read 边数：PE 每 pair 2 边，SE 每 read 1 边）
    t0 = std::time::Instant::now();
    let total_read_sides = input_pairs * if opt.r2.is_some() { 2 } else { 1 };
    let candidates = decide_candidates(&contigs, &aggs, total_read_sides, input_pairs)?;
    eprintln!("[阶段] 统计判定 {:.1}s", t0.elapsed().as_secs_f64());

    // 5. 报告
    let sample = sample_name(&opt.r1);
    t0 = std::time::Instant::now();
    let (result_json, result_tsv) = report::write_report(
        &opt.out,
        &sample,
        opt.threads,
        opt.k,
        input_pairs,
        prescreen_pairs,
        map_errors,
        &candidates,
    )?;
    eprintln!("[阶段] 报告 {:.1}s", t0.elapsed().as_secs_f64());

    Ok(RunSummary {
        input_pairs,
        prescreen_pairs,
        map_errors,
        candidates,
        result_json,
        result_tsv,
    })
}

fn sample_name(r1: &Path) -> String {
    let stem = r1
        .file_stem()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_default();
    for suffix in [".fastq", ".fq"] {
        if let Some(s) = stem.strip_suffix(suffix) {
            return s.to_string();
        }
    }
    stem
}

/// N 总线程 → (D 解压线程, W worker 线程)：D = ⌈(N−1)/4⌉，W = N−1−D。
/// N=8 → (2, 5)；N=1 → (0, 0)（解压与处理全内联主线程）。
fn thread_budget(total: usize) -> (usize, usize) {
    let rest = total.saturating_sub(1);
    let d = rest.div_ceil(4);
    (d, rest.saturating_sub(d))
}

/// 预筛 + 竞争比对融合：单趟流式读入原始 FASTQ，chunk 内先 k-mer 门再比对。
/// 融合理由：pass 中间文件的串行写/读是预筛阶段的串行瓶颈；比对本身
/// 并行吞吐足够，把门控移入工作线程可整段消除该 IO，且减少一个中间产物
/// （无过度设计，反是简化）。证据仍按 chunk 序确定性输出。
/// `decomp_threads`：成员级并行解压预算；为 0 或数据无 IG 索引时使用单线程
/// `MultiGzDecoder`，解析结果与并行路径一致。
fn prescreen_and_align(
    r1: &Path,
    r2: Option<&Path>,
    k: usize,
    bloom: &prescreen::KmerBloom,
    aligner: &align::CompetitiveAligner,
    roles: &HashMap<String, Role>,
    decomp_threads: usize,
) -> Result<(u64, u64, u64, Vec<align::FragmentEvidence>), String> {
    let items: Box<dyn Iterator<Item = Result<fastq::SpanPair, String>>> = match r2 {
        Some(r2_path) => Box::new(fastq::PairSpanIter::open_parallel(
            r1,
            r2_path,
            decomp_threads,
        )?),
        None => Box::new(fastq::SingleSpanIter::open_parallel(r1, decomp_threads)?),
    };
    let input = std::sync::atomic::AtomicU64::new(0);
    let passed = std::sync::atomic::AtomicU64::new(0);
    let map_errors = std::sync::atomic::AtomicU64::new(0);
    let items = items.inspect(|r| {
        if r.is_ok() {
            input.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
        }
    });
    let mut evidences: Vec<align::FragmentEvidence> = Vec::new();
    align::par_map_chunks(
        items,
        aligner.threads,
        1024,
        |chunk| {
            let mut out = Vec::with_capacity(chunk.len());
            let mut scratch = prescreen::GateScratch::default();
            for pair in chunk {
                // k-mer 门：双端任一 read 通过才比对；单端 r2 为空自然只看 s1。
                // 热路径只借共享缓冲切片，零每记录分配。
                let s1 = pair.r1.seq();
                let s2 = pair.r2.seq();
                if !prescreen::read_passes_gate_scratched(s1, k, bloom, &mut scratch)
                    && !prescreen::read_passes_gate_scratched(s2, k, bloom, &mut scratch)
                {
                    continue;
                }
                passed.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                // ID 仅在通过门后物化（每样本仅 ~20 对），热路径保持零分配。
                let key = String::from_utf8_lossy(pair.r1.id()).into_owned();
                let mapped = if pair.r2.has_seq() {
                    aligner.map_pair(s1, s2, &key)
                } else {
                    // 单端：r2 无序列，仅裁决 r1
                    aligner.map_seq(s1, &key).map(|m| (m, Vec::new()))
                };
                match mapped {
                    Ok((m1, m2)) => {
                        let hits1 = align::hits_of(&m1, roles);
                        let hits2 = align::hits_of(&m2, roles);
                        out.push(align::adjudicate_fragment(&key, &hits1, &hits2));
                    }
                    Err(_) => {
                        // 不中止整趟运行，但计数并在报告中披露（避免静默假阴性）
                        map_errors.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                    }
                }
            }
            out
        },
        |evs| {
            evidences.extend(evs);
            Ok(())
        },
    )?;
    Ok((
        input.into_inner(),
        passed.into_inner(),
        map_errors.into_inner(),
        evidences,
    ))
}

/// 每个 contig 的证据聚合。
#[derive(Default)]
struct ContigAgg {
    /// 覆盖区间（半开 [start, end)）。
    intervals: Vec<(i32, i32)>,
    split_events: Vec<cluster::SiteEvent>,
    discordant: u64,
    /// 高置信 reads 数：target 计 confident_target，诱饵计 best_overall
    /// （λ̂ 零分布依赖诱饵层 reads 计数）。
    reads: u64,
    plus: u64,
    minus: u64,
}

/// 把单个 fragment 的证据并入（局部）聚合表。
fn aggregate_one(ev: &align::FragmentEvidence, aggs: &mut HashMap<String, ContigAgg>) {
    for rd in [&ev.r1, &ev.r2] {
        // 覆盖区间 = 裁决命中 ∪ split 半段；split reads 的目标侧片段同样计入覆盖。
        // 重叠区间合并后取并集长度，天然去重。
        let mut coverage_hits: Vec<&align::Hit> = Vec::new();
        let primary = rd.confident_target.as_ref().or(rd.best_overall.as_ref());
        if let Some(h) = primary {
            coverage_hits.push(h);
        }
        if let Some(h) = rd.split_target.as_ref() {
            coverage_hits.push(h);
        }
        if let Some(h) = rd.split_host.as_ref() {
            coverage_hits.push(h);
        }
        let mut seen = std::collections::HashSet::new();
        for hit in coverage_hits {
            let key = (&hit.contig, hit.tstart, hit.tend);
            if !seen.insert(key) {
                continue;
            }
            let agg = aggs.entry(hit.contig.clone()).or_default();
            agg.intervals.push((hit.tstart, hit.tend));
        }
        if let Some(hit) = rd.confident_target.as_ref() {
            let agg = aggs.entry(hit.contig.clone()).or_default();
            agg.reads += 1;
            if hit.strand == '+' {
                agg.plus += 1;
            } else {
                agg.minus += 1;
            }
        } else if let Some(hit) = rd.best_overall.as_ref() {
            // 诱饵 reads 计数：最优命中为诱饵时计入其 reads；判定链 λ̂ 的
            // 零分布依赖该计数，否则诱饵层 reads 恒为 0，导致 λ̂=0、p=0。
            if hit.role == Role::Decoy {
                let agg = aggs.entry(hit.contig.clone()).or_default();
                agg.reads += 1;
                if hit.strand == '+' {
                    agg.plus += 1;
                } else {
                    agg.minus += 1;
                }
            }
        }
    }
    for se in &ev.split_events {
        aggs.entry(se.target_contig.clone())
            .or_default()
            .split_events
            .push(cluster::SiteEvent {
                contig: se.target_contig.clone(),
                pos: se.target_pos as i64,
                host_contig: se.host_contig.clone(),
                host_pos: se.host_pos as i64,
                direction: se.direction.clone(),
                qname: se.qname.clone(),
            });
    }
    if let Some(target) = &ev.discordant {
        aggs.entry(target.clone()).or_default().discordant += 1;
    }
}

/// 每个 contig 的证据聚合：chunk 并行出局部聚合表，主线程有序合并。
/// 合并只做区间拼接与计数相加（交换、可结合），与顺序无关，不影响确定性结果。
fn aggregate_evidence(
    contigs: &[reference::Contig],
    evidences: &[align::FragmentEvidence],
    threads: usize,
) -> Result<HashMap<String, ContigAgg>, String> {
    let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
    for c in contigs {
        aggs.insert(c.name.clone(), ContigAgg::default());
    }
    let items = evidences.iter().map(Ok);
    // --threads N 语义：N 为总线程数，本阶段 worker 取 N−1（主线程占 1）。
    let workers = threads.saturating_sub(1);
    align::par_map_chunks(
        items,
        workers,
        1024,
        |chunk| {
            let mut partial: HashMap<String, ContigAgg> = HashMap::new();
            for ev in chunk {
                aggregate_one(ev, &mut partial);
            }
            partial.into_iter().collect()
        },
        |partial: Vec<(String, ContigAgg)>| {
            for (name, p) in partial {
                let agg = aggs
                    .get_mut(&name)
                    .ok_or_else(|| format!("未知 contig: {name}"))?;
                agg.intervals.extend(p.intervals);
                agg.split_events.extend(p.split_events);
                agg.discordant += p.discordant;
                agg.reads += p.reads;
                agg.plus += p.plus;
                agg.minus += p.minus;
            }
            Ok(())
        },
    )?;
    Ok(aggs)
}

/// (size_bin, gc_bin) 分层键；边界由本模块集中定义。
fn size_bin(len: u64) -> &'static str {
    if len < 1_000 {
        "<1kb"
    } else if len < 5_000 {
        "1-5kb"
    } else if len < 50_000 {
        "5-50kb"
    } else {
        ">=50kb"
    }
}

fn gc_bin(gc: f64) -> &'static str {
    if gc < 0.40 {
        "<0.40"
    } else if gc < 0.50 {
        "0.40-0.50"
    } else if gc < 0.60 {
        "0.50-0.60"
    } else {
        ">=0.60"
    }
}

fn strata(len: u64, gc: f64) -> String {
    format!("sz:{},gc:{}", size_bin(len), gc_bin(gc))
}

/// 区间合并后覆盖碱基数；仅当 start < current_end 时合并。
fn merged_covered_bases(intervals: &[(i32, i32)]) -> u64 {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut bases = 0u64;
    let mut cur_start = 0i32;
    let mut cur_end = 0i32;
    for (s, e) in sorted {
        if s < cur_end {
            cur_end = cur_end.max(e);
        } else {
            bases += (cur_end - cur_start).max(0) as u64;
            cur_start = s;
            cur_end = e.max(s);
        }
    }
    bases += (cur_end - cur_start).max(0) as u64;
    bases
}

/// 合并后的区间数（与 merged_covered_bases 同合并规则，只数正长区间）。
fn merged_interval_count(intervals: &[(i32, i32)]) -> usize {
    let mut sorted = intervals.to_vec();
    sorted.sort_unstable();
    let mut count = 0usize;
    let mut cur_end = 0i32;
    let mut open_start = 0i32;
    let mut open = false;
    for (s, e) in sorted {
        if s < cur_end {
            cur_end = cur_end.max(e);
        } else {
            if open && cur_end > open_start {
                count += 1;
            }
            open_start = s;
            cur_end = e.max(s);
            open = true;
        }
    }
    if open && cur_end > open_start {
        count += 1;
    }
    count
}

/// 每目标 Poisson 计数链 + 全样本离散感知 Storey q（CDH2018 π0）+ 决策。
///
/// 层背景率 λ̂_ℓ = trimmed_mean_20({n_d/L_d})·(10^7/N_total)，深度归一化到 10M
/// reads 参考深度；与暴露量 E_t = L_t·N_total/10^7 相乘时 N 抵消，Poisson 期望
/// = 长度校正的诱饵 reads 期望（跨层/跨样本可比）。
/// n_plain = reads − 唯一 split qname 数（同一 read 双侧 softclip 只扣一次）。
/// NB 对照 α̂ 取同层诱饵原始 reads 计数（与 n_plain 同计数尺度；
/// <20 层 MoM，≥20 层 MoM 初值 + Newton MLE，不收敛回退 MoM）。
/// 决策：q≥0.2 → NOT_SIGNIFICANT；q<0.2 但不满足深度×覆盖组合门槛
/// （DEPTH_RPM_MIN / COVERAGE_MIN）→ BELOW_THRESHOLD；q<0.2 且过门槛：
/// split_events>0 → PASS，无 split 确认 → WEAK。
/// 当前 split 确认层使用 split_events>0；未计算 LLR_split，因为 μ_s⁺/μ_s⁻
/// 的估计方式未定义。
fn decide_candidates(
    contigs: &[reference::Contig],
    aggs: &HashMap<String, ContigAgg>,
    total_read_sides: u64,
    input_pairs: u64,
) -> Result<Vec<report::Candidate>, String> {
    let n_total = total_read_sides.max(1) as f64;
    let depth_norm = 1e7 / n_total;

    // 各分层 decoy 背景：归一化率（λ̂ 口径）+ 原始 reads 计数（NB α̂ 口径）
    let mut decoy_strata: HashMap<String, Vec<f64>> = HashMap::new();
    let mut decoy_counts: HashMap<String, Vec<f64>> = HashMap::new();
    let mut all_rates: Vec<f64> = Vec::new();
    for c in contigs.iter().filter(|c| c.role == Role::Decoy) {
        let agg = &aggs[&c.name];
        let stratum = strata(c.len() as u64, c.gc_frac);
        let rate = agg.reads as f64 / c.len() as f64 * depth_norm;
        decoy_strata.entry(stratum.clone()).or_default().push(rate);
        decoy_counts
            .entry(stratum)
            .or_default()
            .push(agg.reads as f64);
        all_rates.push(rate);
    }
    let lambda_glob = stats::trimmed_mean_20(&all_rates).unwrap_or(0.0);

    struct RawCandidate {
        contig: String,
        len: u64,
        bases: u64,
        windows: u64,
        frac: f64,
        reads: u64,
        n_plain: u64,
        split_events: Vec<cluster::SiteEvent>,
        discordant: u64,
        plus: u64,
        minus: u64,
        stratum: String,
        expected: f64,
        lambda_layer: f64,
        p: f64,
        /// reads per million input fragments（Mourik 2024 RPM 口径）。
        rpm: f64,
        stratum_decoy_count: usize,
    }

    let mut raw: Vec<RawCandidate> = Vec::new();
    for c in contigs.iter().filter(|c| c.role == Role::Target) {
        let agg = &aggs[&c.name];
        let bases = merged_covered_bases(&agg.intervals);
        let has_evidence = bases > 0 || agg.discordant > 0 || !agg.split_events.is_empty();
        if !has_evidence {
            continue;
        }
        let stratum = strata(c.len() as u64, c.gc_frac);
        let layer = decoy_strata.get(&stratum);
        let lambda_layer = layer.and_then(|v| stats::trimmed_mean_20(v));
        // EB 收缩（κ=2）：层空 → 全局（退化）；层大 → 逼近层均值
        let lambda_star = match lambda_layer {
            Some(l) => stats::eb_shrink_layer(l, layer.map_or(0, Vec::len), lambda_glob, 2.0),
            None => lambda_glob,
        };
        let expected = lambda_star * stats::expected_hits(c.len() as u64, total_read_sides.max(1));
        // n_plain：reads 扣唯一 split qname 数（避免 split read 双计）
        let mut split_qnames: HashSet<&str> = HashSet::new();
        for se in &agg.split_events {
            split_qnames.insert(se.qname.as_str());
        }
        let n_plain = agg.reads.saturating_sub(split_qnames.len() as u64);
        let p = stats::poisson_upper_tail(n_plain, expected);
        let rpm = agg.reads as f64 / input_pairs.max(1) as f64 * 1e6;
        raw.push(RawCandidate {
            contig: c.name.clone(),
            len: c.len() as u64,
            bases,
            windows: merged_interval_count(&agg.intervals) as u64,
            frac: bases as f64 / c.len() as f64,
            reads: agg.reads,
            n_plain,
            split_events: agg.split_events.clone(),
            discordant: agg.discordant,
            plus: agg.plus,
            minus: agg.minus,
            stratum,
            expected,
            lambda_layer: lambda_star,
            p,
            rpm,
            stratum_decoy_count: layer.map_or(0, Vec::len),
        });
    }

    // 全样本池化离散感知 Storey q（CDH2018 π0 + aBH 单调 q）
    let ps: Vec<f64> = raw.iter().map(|r| r.p).collect();
    let rates: Vec<f64> = raw.iter().map(|r| r.expected).collect();
    let (pi0, qs) = stats::discrete_q_values(&ps, &rates);

    let mut candidates = Vec::with_capacity(raw.len());
    for (i, r) in raw.into_iter().enumerate() {
        // NB 过度离散对照（Robinson & Smyth 2008 α̂）
        let nb_p = if let Some(counts) = decoy_counts.get(&r.stratum) {
            if counts.len() < 2 {
                None
            } else {
                let mean = counts.iter().sum::<f64>() / counts.len() as f64;
                let var = counts
                    .iter()
                    .map(|x| {
                        let d = x - mean;
                        d * d
                    })
                    .sum::<f64>()
                    / (counts.len() - 1) as f64;
                let alpha = if counts.len() >= 20 {
                    let mom = stats::nb_alpha_mom(mean, var).max(1e-4);
                    stats::nb_alpha_mle(counts, mom).unwrap_or(mom)
                } else {
                    stats::nb_alpha_mom(mean, var)
                };
                Some(stats::nb_upper_tail(
                    r.n_plain,
                    r.expected,
                    alpha.max(1e-10),
                ))
            }
        } else {
            None
        };
        // p 地板固定为 1/6：仅披露 flag，不进 q（p 已是连续化精确上尾）。
        let p_floor_flag = r.p <= 1.0 / 6.0;
        // 证据层组合门槛：AND 判定，不乘入 p/q。
        let decision = if qs[i] >= 0.2 {
            "NOT_SIGNIFICANT"
        } else if r.rpm < DEPTH_RPM_MIN || r.frac < COVERAGE_MIN {
            "BELOW_THRESHOLD"
        } else if !r.split_events.is_empty() {
            "PASS"
        } else {
            "WEAK"
        };
        let sites = cluster::cluster_sites(&r.split_events, cluster::MIN_SITE_SUPPORT);
        candidates.push(report::Candidate {
            contig: r.contig,
            contig_len: r.len,
            covered_bases: r.bases,
            covered_frac: r.frac,
            reads: r.reads,
            split_events: r.split_events.len() as u64,
            discordant: r.discordant,
            plus_strand: r.plus,
            minus_strand: r.minus,
            sites,
            p_value: r.p,
            q_value: qs[i],
            p_resolution_floor: 1.0 / (r.stratum_decoy_count as f64 + 1.0),
            poisson_p: Some(r.p),
            nb_p,
            stratum: r.stratum,
            stratum_decoy_count: r.stratum_decoy_count,
            n_plain: r.n_plain,
            expected_hits: r.expected,
            lambda_bg_layer: r.lambda_layer,
            depth_fold: if r.expected > 0.0 {
                Some(r.n_plain as f64 / r.expected)
            } else {
                None
            },
            p_floor_flag,
            pi0,
            decision,
            distinct_windows: r.windows,
            depth_rpm: r.rpm,
        });
    }
    // 按 q 升序输出
    candidates.sort_by(|a, b| {
        a.q_value
            .partial_cmp(&b.q_value)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    Ok(candidates)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn merge_intervals_counts_union_length() {
        assert_eq!(merged_covered_bases(&[(0, 100), (50, 150)]), 150);
        assert_eq!(merged_covered_bases(&[(0, 100), (200, 250)]), 150);
        assert_eq!(merged_covered_bases(&[(100, 50)]), 0); // 非法区间防御
        assert_eq!(merged_covered_bases(&[]), 0);
    }

    #[test]
    fn strata_binning() {
        assert_eq!(strata(800, 0.45), "sz:<1kb,gc:0.40-0.50");
        assert_eq!(strata(3000, 0.55), "sz:1-5kb,gc:0.50-0.60");
        assert_eq!(strata(100_000, 0.65), "sz:>=50kb,gc:>=0.60");
    }

    #[test]
    fn merged_interval_count_matches_merge_rule() {
        assert_eq!(merged_interval_count(&[(0, 100), (50, 150)]), 1);
        assert_eq!(merged_interval_count(&[(0, 100), (200, 250)]), 2);
        assert_eq!(merged_interval_count(&[(100, 50)]), 0);
        assert_eq!(merged_interval_count(&[]), 0);
    }

    #[test]
    fn decide_chain_pass_and_null_with_pooled_pi0() {
        // 2 目标各占一层、各 1 条零率诱饵：
        // target_0 reads=10、split 唯一 qname 1 条 → n_plain=9；λ̂=0 → 期望 0 → p=0。
        // 组合门槛：input_pairs=10_000 → rpm=1000≥500，intervals 覆盖 300/3000=10%≥10%，
        // split 非空 → PASS。
        // target_1 仅 discordant=1 → n_plain=0；p=P(Pois(0)≥0)=1 → q=1 → NOT_SIGNIFICANT。
        // 池化 π0=1.0（已手算：CDH2018 各 τ 估计均被钳制到 1）。
        let contigs = vec![
            reference::Contig {
                name: "target_0".into(),
                role: Role::Target,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
            reference::Contig {
                name: "target_1".into(),
                role: Role::Target,
                seq: vec![b'A'; 30000],
                gc_frac: 0.55,
            },
            reference::Contig {
                name: "decoy_0".into(),
                role: Role::Decoy,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
            reference::Contig {
                name: "decoy_1".into(),
                role: Role::Decoy,
                seq: vec![b'A'; 30000],
                gc_frac: 0.55,
            },
        ];
        let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
        for c in &contigs {
            aggs.insert(c.name.clone(), ContigAgg::default());
        }
        aggs.get_mut("target_0").unwrap().reads = 10;
        aggs.get_mut("target_0").unwrap().intervals = vec![(0, 300)];
        aggs.get_mut("target_0").unwrap().split_events = vec![cluster::SiteEvent {
            contig: "target_0".into(),
            pos: 1,
            host_contig: "host".into(),
            host_pos: 1,
            direction: "+".into(),
            qname: "a".into(),
        }];
        aggs.get_mut("target_1").unwrap().discordant = 1;
        let cands = decide_candidates(&contigs, &aggs, 20_000, 10_000).unwrap();
        assert_eq!(cands.len(), 2);
        // 按 q 升序：target_0 在前
        let t0 = &cands[0];
        assert_eq!(t0.contig, "target_0");
        assert_eq!(t0.n_plain, 9);
        assert!(t0.p_value == 0.0, "p={}", t0.p_value);
        assert!(t0.q_value == 0.0, "q={}", t0.q_value);
        assert_eq!(t0.decision, "PASS");
        assert!((t0.depth_rpm - 1000.0).abs() < 1e-9, "rpm={}", t0.depth_rpm);
        assert!((t0.pi0 - 1.0).abs() < 1e-12, "π̂0={}", t0.pi0);
        assert_eq!(t0.stratum_decoy_count, 1);
        assert!((t0.p_resolution_floor - 0.5).abs() < 1e-9);
        let t1 = &cands[1];
        assert_eq!(t1.contig, "target_1");
        assert_eq!(t1.n_plain, 0);
        assert!((t1.p_value - 1.0).abs() < 1e-12, "p={}", t1.p_value);
        assert!((t1.q_value - 1.0).abs() < 1e-12, "q={}", t1.q_value);
        assert_eq!(t1.decision, "NOT_SIGNIFICANT");
        assert_eq!(t1.stratum_decoy_count, 1);
    }

    #[test]
    fn n_plain_dedups_split_qnames() {
        // reads=3、split qnames ["a","a","b"] → 唯一 2 → n_plain=1。
        let contigs = vec![
            reference::Contig {
                name: "target_0".into(),
                role: Role::Target,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
            reference::Contig {
                name: "decoy_0".into(),
                role: Role::Decoy,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
        ];
        let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
        for c in &contigs {
            aggs.insert(c.name.clone(), ContigAgg::default());
        }
        aggs.get_mut("target_0").unwrap().reads = 3;
        let ev = |q: &str| cluster::SiteEvent {
            contig: "target_0".into(),
            pos: 1,
            host_contig: "host".into(),
            host_pos: 1,
            direction: "+".into(),
            qname: q.into(),
        };
        aggs.get_mut("target_0").unwrap().split_events = vec![ev("a"), ev("a"), ev("b")];
        let cands = decide_candidates(&contigs, &aggs, 1_000_000, 500_000).unwrap();
        assert_eq!(cands.len(), 1);
        assert_eq!(cands[0].n_plain, 1);
        assert_eq!(cands[0].split_events, 3);
    }

    #[test]
    fn depth_coverage_gate_belows_threshold() {
        // 组合门槛为 AND 判定：q=0 但深度或覆盖任一不过 → BELOW_THRESHOLD。
        // deep 目标：rpm=5 < 500（深度不过）；wide 目标：frac=0.02 < 0.10（覆盖不过）。
        // 两目标分层相同、诱饵零 reads → λ̂=0 → p=0 → q=0（无门槛时本应为 PASS/WEAK）。
        let contigs = vec![
            reference::Contig {
                name: "deep".into(),
                role: Role::Target,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
            reference::Contig {
                name: "wide".into(),
                role: Role::Target,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
            reference::Contig {
                name: "decoy_0".into(),
                role: Role::Decoy,
                seq: vec![b'A'; 3000],
                gc_frac: 0.5,
            },
        ];
        let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
        for c in &contigs {
            aggs.insert(c.name.clone(), ContigAgg::default());
        }
        aggs.get_mut("deep").unwrap().reads = 5; // rpm = 5/1e6*1e6 = 5
        aggs.get_mut("deep").unwrap().intervals = vec![(0, 300)]; // frac=0.10 过覆盖门槛
        aggs.get_mut("wide").unwrap().reads = 10_000; // rpm = 10000 过深度门槛
        aggs.get_mut("wide").unwrap().intervals = vec![(0, 60)]; // frac=0.02 不过覆盖门槛
        let cands = decide_candidates(&contigs, &aggs, 2_000_000, 1_000_000).unwrap();
        assert_eq!(cands.len(), 2);
        for c in &cands {
            assert!(c.q_value == 0.0, "{} q={}", c.contig, c.q_value);
            assert_eq!(c.decision, "BELOW_THRESHOLD", "{}", c.contig);
        }
        let deep = cands.iter().find(|c| c.contig == "deep").unwrap();
        let wide = cands.iter().find(|c| c.contig == "wide").unwrap();
        assert!((deep.depth_rpm - 5.0).abs() < 1e-9);
        assert!((wide.depth_rpm - 10_000.0).abs() < 1e-6);
    }

    #[test]
    fn run_pipeline_rejects_out_of_range_k() {
        let mut opt = Options {
            k: 0,
            ..Options::default()
        };
        let err = run_pipeline(&opt).unwrap_err();
        assert!(err.contains("--k"), "err={err}");
        opt.k = 32;
        let err = run_pipeline(&opt).unwrap_err();
        assert!(err.contains("--k"), "err={err}");
    }
}
