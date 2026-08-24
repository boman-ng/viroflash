//! 病毒候选检测管线编排：
//! 索引加载/构建 → 确定性 bottom-k 采样 → Bloom 门控与 ALL_CHAINS 竞争审计 →
//! discovery 兼容集构造 → validation 独立计数与 decoy 统计判定 → 报告。

pub mod align;
pub mod cluster;
pub mod decoy;
pub mod equivalence;
pub mod fastq;
pub mod group;
pub mod hash;
pub mod index;
pub mod perf;
pub mod prescreen;
pub mod reference;
pub mod report;
pub mod sampling;
pub mod stats;

use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::path::{Path, PathBuf};

use reference::Role;

#[derive(Debug, Clone)]
pub struct Options {
    pub r1: PathBuf,
    /// None = 单端模式（454/单端测序），每 read 一个 fragment。
    pub r2: Option<PathBuf>,
    /// 已构建索引目录；与四类 FASTA 参数互斥（main.rs 解析与 run_pipeline 双重校验）。
    /// None = 自动构建（索引产物落 `<out>.work/index/`）。
    pub index: Option<PathBuf>,
    pub host_fa: Option<PathBuf>,
    pub target_fa: Option<PathBuf>,
    pub decoy_fa: Option<PathBuf>,
    pub contam_fa: Option<PathBuf>,
    pub threads: usize,
    pub out: PathBuf,
    pub k: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            r1: PathBuf::new(),
            r2: None,
            index: None,
            host_fa: None,
            target_fa: None,
            decoy_fa: None,
            contam_fa: None,
            threads: 8,
            out: PathBuf::from("viroflash_out"),
            k: prescreen::DEFAULT_K,
        }
    }
}

/// 证据层 corroboration：Mourik 2024 在特定 pan-viral capture 流程中评估了
/// 10% genome breadth；Miller 2019 在 CSF mNGS 中采用至少 3 个非重叠区域。
/// 两者都没有验证本项目的跨样本通用阈值。本实现公开采用代表序列 10 等分，以每条
/// validation read 对齐中点占据的不同分区数衡量；它们是待外部验证的报告门，
/// 不乘入模型 p/adjusted-p。
pub const COVERAGE_MIN: f64 = 0.10;
pub const MIN_DISTRIBUTED_WINDOWS: u64 = 3;
pub const DISTRIBUTED_WINDOW_BINS: u64 = 10;
/// 模型 adjusted-p 的探索性候选报告阈值；尚未由独立验证集校准。
pub const MODEL_ADJUSTED_P_MAX: f64 = 0.20;

#[derive(Debug, Clone)]
pub struct RunSummary {
    pub input_pairs: u64,
    /// 在 bottom-k 入选样本中通过 k-mer 门的 pair 数；`input_pairs <= K` 时即全量
    /// 精确值，否则作用域由 `sampling` 显式披露为 selected sample。
    pub prescreen_pairs: u64,
    /// 比对失败的 fragment 数（minimap2 返回错误；不中止运行但显式披露）。
    pub map_errors: u64,
    pub sampling: report::SamplingReport,
    /// discovery 固定的统计检验族大小，包含 validation 零证据假设。
    pub test_family_size: usize,
    pub candidates: Vec<report::Candidate>,
    pub result_json: PathBuf,
    pub result_tsv: PathBuf,
}

pub fn run_pipeline(opt: &Options) -> Result<RunSummary, String> {
    validate_run_options(opt)?;
    let sample = sample_name(&opt.r1);
    let monitor = perf::PerfMonitor::start("run", sample.clone(), &opt.out, opt.threads)?;
    let result = run_pipeline_inner(opt, &sample, &monitor);
    monitor.complete(result)
}

fn validate_run_options(opt: &Options) -> Result<(), String> {
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
    if opt.r1.as_os_str().is_empty() {
        return Err("缺少 --r1".into());
    }
    // --index 与四类 FASTA 互斥：加载路径的角色/参考信息全部来自 manifest，
    // 同时给 FASTA 会造成两个矛盾的参考来源。
    let fasta_args: Vec<&str> = [
        ("--host-fa", &opt.host_fa),
        ("--target-fa", &opt.target_fa),
        ("--decoy-fa", &opt.decoy_fa),
        ("--contam-fa", &opt.contam_fa),
    ]
    .into_iter()
    .filter(|(_, p)| p.is_some())
    .map(|(label, _)| label)
    .collect();
    if opt.index.is_some() && !fasta_args.is_empty() {
        return Err(format!(
            "--index 与 {} 不能同时使用（加载索引时参考信息取自 manifest）",
            fasta_args.join("/")
        ));
    }
    Ok(())
}

fn run_pipeline_inner(
    opt: &Options,
    sample: &str,
    monitor: &perf::PerfMonitor,
) -> Result<RunSummary, String> {
    let work_dir = PathBuf::from(format!("{}.work", opt.out.display()));
    let mut t0 = std::time::Instant::now();
    // 索引解析：--index 直接加载并校验；否则自动构建（诱饵未提供时按固定默认
    // 参数自动生成），与 `viroflash index` 共用同一构建入口，两路径结果一致。
    let (built, index_source) = match &opt.index {
        Some(dir) => {
            monitor.stage("index_load");
            let loaded = index::load_index(dir, opt.k)?;
            eprintln!("[阶段] 索引加载 {:.1}s", t0.elapsed().as_secs_f64());
            (loaded, "loaded")
        }
        None => {
            let host_fa = require_fa(opt.host_fa.as_deref(), "--host-fa")?;
            let target_fa = require_fa(opt.target_fa.as_deref(), "--target-fa")?;
            let spec = index::IndexOptions {
                host_fa: host_fa.to_path_buf(),
                target_fa: target_fa.to_path_buf(),
                contam_fa: opt.contam_fa.clone(),
                decoy_fa: opt.decoy_fa.clone(),
                out_dir: work_dir.join("index"),
                k: opt.k,
                threads: opt.threads,
                ..index::IndexOptions::default()
            };
            let built = index::build_index_with_monitor(&spec, monitor)?;
            eprintln!("[阶段] 参考构建+索引 {:.1}s", t0.elapsed().as_secs_f64());
            (built, "built")
        }
    };

    // 1. 抽样 + 宽松预筛 + 全链直接证据审计。
    if !built.contigs.iter().any(|c| c.role == Role::Target) {
        return Err("目标参考中没有序列（索引未含目标或目标 FASTA 为空）".into());
    }
    // Bloom 覆盖全部原始目标与进入索引的诱饵 k-mer 并集。
    // Bloom 在索引构建阶段生成并持久化（bloom.bin），加载路径直接复用。
    // --threads 保持为数据管线预算；低频 telemetry sampler 不占 mapper/decompressor。
    monitor.stage("index_open_mmi");
    let runtime_threads = monitor.work_thread_budget();
    let (decomp_threads, workers) = thread_budget(runtime_threads);
    let aligner = align::CompetitiveAligner::open(&built.mmi_path, workers)?;
    let role_members = RoleMembers::from_contigs(&built.contigs);
    monitor.stage("sample_prescreen_audit");
    t0 = std::time::Instant::now();
    let audited = sample_and_audit(
        &opt.r1,
        opt.r2.as_deref(),
        opt.k,
        &built.bloom,
        &aligner,
        &built.roles,
        &role_members,
        decomp_threads,
    )?;
    eprintln!(
        "[阶段] 抽样+预筛+全链审计 {:.1}s",
        t0.elapsed().as_secs_f64()
    );

    // 2. discovery 直接集合构造 OR 假设；validation 只验证既有假设。
    monitor.stage("equivalence_aggregate");
    t0 = std::time::Instant::now();
    let composite = build_composite_evidence(
        &built.contigs,
        &built.target_groups,
        &role_members,
        &audited.evidence,
    )?;
    eprintln!("[阶段] 等价类聚合 {:.1}s", t0.elapsed().as_secs_f64());

    // 3. synthetic-decoy 统计判定。暴露量只取独立 validation 折的实际 read-end；
    // discovery 只定义假设，不能再进入检验统计量。
    monitor.stage("statistics");
    t0 = std::time::Instant::now();
    let sides_per_pair = if opt.r2.is_some() { 2 } else { 1 };
    let validation_read_sides = audited
        .validation_pairs
        .checked_mul(sides_per_pair)
        .ok_or_else(|| "validation read-end 计数溢出".to_string())?;
    let decision = decide_candidates(
        &built.contigs,
        &composite.aggs,
        &composite.decoy_background,
        validation_read_sides,
        Some(&composite.details),
    )?;
    let candidates = decision.candidates;
    let test_family_size = decision.test_family_size;
    eprintln!("[阶段] 统计判定 {:.1}s", t0.elapsed().as_secs_f64());

    let input_pairs = audited.input_pairs;
    let prescreen_pairs = audited.selected_passed_pairs;
    let map_errors = audited.evidence.map_error_pairs;
    let sampling = report::SamplingReport {
        method: "deterministic_bottom_k_blake3",
        capacity_pairs: sampling::DEFAULT_SAMPLE_PAIRS,
        selected_pairs: audited.selected_pairs,
        selected_passed_pairs: audited.selected_passed_pairs,
        discovery_pairs: audited.discovery_pairs,
        validation_pairs: audited.validation_pairs,
        validation_read_sides,
        inclusion_probability: audited.inclusion_probability,
        duplicate_selected_qnames: audited.duplicate_selected_qnames,
        audit_reads: audited.evidence.audit_reads,
        audit_no_evidence: audited.evidence.audit_no_evidence,
        audit_overflows: audited.evidence.audit_overflows,
        audit_map_error_reads: audited.evidence.audit_map_error_reads,
        target_validation_unassigned: composite.target_validation_unassigned,
        decoy_validation_unassigned: composite.decoy_validation_unassigned,
    };

    // 4. 报告
    monitor.stage("report");
    t0 = std::time::Instant::now();
    let (result_json, result_tsv) = report::write_report(
        &opt.out,
        sample,
        opt.threads,
        opt.k,
        input_pairs,
        prescreen_pairs,
        map_errors,
        index_source,
        built.format_version,
        &built.manifest_blake3,
        &sampling,
        test_family_size,
        &candidates,
    )?;
    eprintln!("[阶段] 报告 {:.1}s", t0.elapsed().as_secs_f64());

    // minimap2 索引与 Bloom 占据绝大多数常驻内存；显式释放并单独计时，避免进程
    // 退出前的大对象回收被误记为报告阶段。
    monitor.stage("cleanup");
    t0 = std::time::Instant::now();
    drop(aligner);
    drop(built);
    eprintln!("[阶段] 资源释放 {:.1}s", t0.elapsed().as_secs_f64());

    Ok(RunSummary {
        input_pairs,
        prescreen_pairs,
        map_errors,
        sampling,
        test_family_size,
        candidates,
        result_json,
        result_tsv,
    })
}

/// run 未提供 --index 时自动构建索引所需的必填 FASTA 校验。
fn require_fa<'a>(p: Option<&'a Path>, label: &str) -> Result<&'a Path, String> {
    match p {
        Some(path) if !path.as_os_str().is_empty() => Ok(path),
        _ => Err(format!(
            "缺少 {label}（未提供 --index 时自动构建索引需要它）"
        )),
    }
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

/// N 个管线计算线程 → (D 辅助解压线程, W worker 线程)。真实 IG 采样显示
/// DEFLATE helper 各只消耗约 5% CPU，不应从持续满载的计算预算中一比一扣除；
/// W=N−1，使 main caller-runs + workers 最多占满 N 个核。D 沿用既有约 1:3
/// 配比，但只有输入带 IG 成员索引时才实际创建；普通 gzip 自动回退且不创建。
/// CLI `--threads 8` → main 1 + worker 7，必要时另有 2 个有界低占用 helper；
/// telemetry sampler 也是额外的低频休眠线程。`--threads 1` → 全内联。
fn thread_budget(total: usize) -> (usize, usize) {
    let workers = total.saturating_sub(1);
    let auxiliary_decompressors = workers.div_ceil(4);
    (auxiliary_decompressors, workers)
}

#[derive(Debug, Clone)]
struct DirectPlacement {
    member: usize,
    start: i32,
    end: i32,
    strand: char,
}

/// 单个 read-end 的直接最高分兼容集合。成员数不改变该观测的一票权重。
#[derive(Debug, Clone)]
struct DirectObservation {
    members: Vec<usize>,
    placements: Vec<DirectPlacement>,
}

#[derive(Default)]
struct AuditChunk {
    target_discovery: Vec<Vec<usize>>,
    target_validation: Vec<DirectObservation>,
    decoy_discovery: Vec<Vec<usize>>,
    decoy_validation: Vec<DirectObservation>,
    validation_fragments: Vec<align::FragmentEvidence>,
    audit_reads: u64,
    audit_no_evidence: u64,
    audit_overflows: u64,
    audit_map_error_reads: u64,
    map_error_pairs: u64,
}

impl AuditChunk {
    fn merge(&mut self, mut other: Self) {
        self.target_discovery.append(&mut other.target_discovery);
        self.target_validation.append(&mut other.target_validation);
        self.decoy_discovery.append(&mut other.decoy_discovery);
        self.decoy_validation.append(&mut other.decoy_validation);
        self.validation_fragments
            .append(&mut other.validation_fragments);
        self.audit_reads += other.audit_reads;
        self.audit_no_evidence += other.audit_no_evidence;
        self.audit_overflows += other.audit_overflows;
        self.audit_map_error_reads += other.audit_map_error_reads;
        self.map_error_pairs += other.map_error_pairs;
    }
}

struct AuditRun {
    input_pairs: u64,
    selected_pairs: u64,
    selected_passed_pairs: u64,
    discovery_pairs: u64,
    validation_pairs: u64,
    inclusion_probability: f64,
    duplicate_selected_qnames: u64,
    evidence: AuditChunk,
}

struct RoleMembers {
    lookup: HashMap<String, (Role, usize)>,
    target_names: Vec<String>,
    decoy_names: Vec<String>,
}

impl RoleMembers {
    fn from_contigs(contigs: &[reference::ContigMeta]) -> Self {
        let mut lookup = HashMap::new();
        let mut target_names = Vec::new();
        let mut decoy_names = Vec::new();
        for contig in contigs {
            let names = match contig.role {
                Role::Target => &mut target_names,
                Role::Decoy => &mut decoy_names,
                _ => continue,
            };
            let member = names.len();
            names.push(contig.name.clone());
            lookup.insert(contig.name.clone(), (contig.role, member));
        }
        Self {
            lookup,
            target_names,
            decoy_names,
        }
    }

    fn names(&self, role: Role) -> &[String] {
        match role {
            Role::Target => &self.target_names,
            Role::Decoy => &self.decoy_names,
            _ => &[],
        }
    }
}

#[derive(Default)]
struct MappedRead {
    hits: Vec<align::Hit>,
    direct: Option<(Role, DirectObservation)>,
    map_error: bool,
}

fn audit_read(
    sequence: &[u8],
    qname: &str,
    aligner: &align::CompetitiveAligner,
    roles: &HashMap<String, Role>,
    members: &RoleMembers,
    output: &mut AuditChunk,
) -> MappedRead {
    if sequence.is_empty() {
        return MappedRead::default();
    }
    output.audit_reads += 1;
    let mappings = match aligner.audit_map_seq(sequence, qname) {
        Ok(mappings) => mappings,
        Err(_) => {
            output.audit_map_error_reads += 1;
            return MappedRead {
                map_error: true,
                ..MappedRead::default()
            };
        }
    };
    // 在物化 Hit 前检查 mapper 返回量，避免超限 read 先克隆数千个结构再丢弃。
    if mappings.len() > align::MAX_AUDIT_HITS {
        output.audit_overflows += 1;
        return MappedRead::default();
    }
    let hits = align::hits_of(&mappings, roles);
    let direct = match align::adjudicate_audit(&hits) {
        align::AuditDecision::Resolved { role, hits } => {
            let mut placements = hits
                .into_iter()
                .filter_map(|hit| {
                    let &(hit_role, member) = members.lookup.get(&hit.contig)?;
                    (hit_role == role).then_some(DirectPlacement {
                        member,
                        start: hit.tstart,
                        end: hit.tend,
                        strand: hit.strand,
                    })
                })
                .collect::<Vec<_>>();
            placements.sort_unstable_by_key(|placement| placement.member);
            placements.dedup_by_key(|placement| placement.member);
            let members = placements
                .iter()
                .map(|placement| placement.member)
                .collect::<Vec<_>>();
            (!members.is_empty()).then_some((
                role,
                DirectObservation {
                    members,
                    placements,
                },
            ))
        }
        align::AuditDecision::NoEvidence => {
            output.audit_no_evidence += 1;
            None
        }
        align::AuditDecision::Overflow { .. } => {
            // 上面的 mapping 数门已先拦截；保留该分支使接口未来调整仍显式计数。
            output.audit_overflows += 1;
            None
        }
    };
    MappedRead {
        hits,
        direct,
        map_error: false,
    }
}

fn store_direct(
    fold: sampling::EvidenceFold,
    direct: Option<(Role, DirectObservation)>,
    output: &mut AuditChunk,
) {
    let Some((role, observation)) = direct else {
        return;
    };
    match (role, fold) {
        (Role::Target, sampling::EvidenceFold::Discovery) => {
            output.target_discovery.push(observation.members)
        }
        (Role::Target, sampling::EvidenceFold::Validation) => {
            output.target_validation.push(observation)
        }
        (Role::Decoy, sampling::EvidenceFold::Discovery) => {
            output.decoy_discovery.push(observation.members)
        }
        (Role::Decoy, sampling::EvidenceFold::Validation) => {
            output.decoy_validation.push(observation)
        }
        _ => {}
    }
}

/// 全输入先做确定性 bottom-k，再只对最终样本访问 Bloom 和 ALL_CHAINS mapper。
/// 发现/验证按 fragment 折分，两端始终同折；直接兼容集合按 read-end 各计一票。
#[allow(clippy::too_many_arguments)]
fn sample_and_audit(
    r1: &Path,
    r2: Option<&Path>,
    k: usize,
    bloom: &prescreen::KmerBloom,
    aligner: &align::CompetitiveAligner,
    roles: &HashMap<String, Role>,
    role_members: &RoleMembers,
    decomp_threads: usize,
) -> Result<AuditRun, String> {
    let items: Box<dyn Iterator<Item = Result<fastq::SpanPair, String>>> = match r2 {
        Some(r2_path) => Box::new(fastq::PairSpanIter::open_parallel(
            r1,
            r2_path,
            decomp_threads,
        )?),
        None => Box::new(fastq::SingleSpanIter::open_parallel(r1, decomp_threads)?),
    };
    let mut reservoir = sampling::PairReservoir::new(sampling::DEFAULT_SAMPLE_PAIRS)?;
    let mut ordinal = 0u64;
    for pair in items {
        let pair = pair?;
        reservoir.observe(ordinal, pair.r1.id(), pair.r1.seq(), pair.r2.seq())?;
        ordinal = ordinal
            .checked_add(1)
            .ok_or_else(|| "fragment 序号溢出".to_string())?;
    }
    let sampled = reservoir.finish();
    let mut seen_qnames = HashSet::with_capacity(sampled.pairs.len());
    let duplicate_selected_qnames = sampled
        .pairs
        .iter()
        .filter(|pair| !seen_qnames.insert(pair.qname.as_str()))
        .count() as u64;
    let (discovery_pairs, validation_pairs) =
        sampled.pairs.iter().fold((0u64, 0u64), |counts, pair| {
            match sampling::evidence_fold(&pair.qname, pair.ordinal) {
                sampling::EvidenceFold::Discovery => (counts.0 + 1, counts.1),
                sampling::EvidenceFold::Validation => (counts.0, counts.1 + 1),
            }
        });

    let mut scratch = prescreen::GateScratch::default();
    let passed_pairs = sampled
        .pairs
        .into_iter()
        .filter(|pair| {
            prescreen::read_passes_gate_scratched(&pair.r1_seq, k, bloom, &mut scratch)
                || prescreen::read_passes_gate_scratched(&pair.r2_seq, k, bloom, &mut scratch)
        })
        .collect::<Vec<_>>();
    let selected_passed_pairs = passed_pairs.len() as u64;
    let mut evidence = AuditChunk::default();
    align::par_map_chunks(
        passed_pairs.into_iter().map(Ok),
        aligner.threads,
        64,
        |chunk| {
            let mut out = AuditChunk::default();
            for pair in chunk {
                let fold = sampling::evidence_fold(&pair.qname, pair.ordinal);
                let r1 = audit_read(
                    &pair.r1_seq,
                    &pair.qname,
                    aligner,
                    roles,
                    role_members,
                    &mut out,
                );
                let r2 = audit_read(
                    &pair.r2_seq,
                    &pair.qname,
                    aligner,
                    roles,
                    role_members,
                    &mut out,
                );
                if r1.map_error || r2.map_error {
                    out.map_error_pairs += 1;
                }
                store_direct(fold, r1.direct, &mut out);
                store_direct(fold, r2.direct, &mut out);
                if fold == sampling::EvidenceFold::Validation {
                    let evidence_key = format!("{}#{}", pair.qname, pair.ordinal);
                    out.validation_fragments.push(align::adjudicate_fragment(
                        &evidence_key,
                        &r1.hits,
                        &r2.hits,
                    ));
                }
            }
            vec![out]
        },
        |chunks| {
            for chunk in chunks {
                evidence.merge(chunk);
            }
            Ok(())
        },
    )?;

    Ok(AuditRun {
        input_pairs: sampled.seen_pairs,
        selected_pairs: sampled.selected_pairs,
        selected_passed_pairs,
        discovery_pairs,
        validation_pairs,
        inclusion_probability: sampled.inclusion_probability,
        duplicate_selected_qnames,
        evidence,
    })
}

/// 每个内部代表 contig 的 validation 证据聚合。
#[derive(Default)]
struct ContigAgg {
    /// 覆盖区间（半开 [start, end)）。
    intervals: Vec<(i32, i32)>,
    split_events: Vec<cluster::SiteEvent>,
    discordant: u64,
    /// 直接兼容 read-end 数；每个观测无论兼容成员多少只计一次。
    reads: u64,
    plus: u64,
    minus: u64,
}

/// 固定 decoy 参考集在 validation 折中的直接证据。与 discovery 是否观察到该
/// decoy 无关，避免把“未被 discovery 选中”误当成 validation 的结构性零计数。
#[derive(Debug, Clone, Copy, Default)]
struct DecoyLayerBackground {
    reads: u64,
    reference_bases: u64,
    contigs: usize,
}

#[derive(Default)]
struct DecoyBackground {
    layers: HashMap<String, DecoyLayerBackground>,
    global: DecoyLayerBackground,
    /// 一个 read-end 的直接最高分兼容 decoy 跨越多个 size/GC 层；该 read 在每个
    /// 兼容层中至多计一次，使对应目标检验保持保守而不在单层重复计票。
    cross_stratum_reads: u64,
}

#[derive(Default)]
struct CompositeDetails {
    representative: String,
    members: Vec<String>,
    explanation: Vec<String>,
    discovery_reads: u64,
    /// discovery OR 假设实际包含的内部索引 target contig 数；与展开后的原始
    /// accession 数不同。
    index_member_count: usize,
    /// 上述内部索引 contig 的长度总和，是聚合 target read 计数对应的 exposure。
    target_exposure_bases: u64,
}

struct ResolvedHypothesis {
    hypothesis: equivalence::CompositeHypothesis,
    representative: usize,
}

struct RoleAggregation {
    member_to_representative: Vec<Option<String>>,
    hypotheses: Vec<ResolvedHypothesis>,
    validation_unassigned: u64,
}

fn aggregate_role_hypotheses(
    role: Role,
    discovery: &[Vec<usize>],
    validation: &[DirectObservation],
    role_members: &RoleMembers,
    aggs: &mut HashMap<String, ContigAgg>,
) -> Result<RoleAggregation, String> {
    let hypotheses = equivalence::build_composite_hypotheses(discovery);
    let mut support = vec![BTreeMap::<usize, u64>::new(); hypotheses.len()];
    for members in discovery {
        if let Some(index) = equivalence::validation_hypothesis_index(&hypotheses, members) {
            for &member in members {
                *support[index].entry(member).or_default() += 1;
            }
        }
    }
    let names = role_members.names(role);
    let mut member_to_representative = vec![None; names.len()];
    let mut resolved = Vec::with_capacity(hypotheses.len());
    for (index, hypothesis) in hypotheses.into_iter().enumerate() {
        let representative = hypothesis
            .members
            .iter()
            .copied()
            .max_by(|left, right| {
                support[index]
                    .get(left)
                    .copied()
                    .unwrap_or(0)
                    .cmp(&support[index].get(right).copied().unwrap_or(0))
                    .then_with(|| right.cmp(left))
            })
            .ok_or_else(|| "复合假设没有成员".to_string())?;
        let representative_name = names
            .get(representative)
            .ok_or_else(|| format!("{role:?} 代表成员下标越界: {representative}"))?
            .clone();
        for &member in &hypothesis.members {
            let slot = member_to_representative
                .get_mut(member)
                .ok_or_else(|| format!("{role:?} 假设成员下标越界: {member}"))?;
            *slot = Some(representative_name.clone());
        }
        resolved.push(ResolvedHypothesis {
            hypothesis,
            representative,
        });
    }

    let hypotheses = resolved
        .iter()
        .map(|resolved| resolved.hypothesis.clone())
        .collect::<Vec<_>>();
    let mut unassigned = 0u64;
    for observation in validation {
        let Some(hypothesis_index) =
            equivalence::validation_hypothesis_index(&hypotheses, &observation.members)
        else {
            unassigned += 1;
            continue;
        };
        let resolved_hypothesis = &resolved[hypothesis_index];
        let representative_name = names
            .get(resolved_hypothesis.representative)
            .ok_or_else(|| "validation 代表成员下标越界".to_string())?;
        let agg = aggs
            .get_mut(representative_name)
            .ok_or_else(|| format!("未知代表 contig: {representative_name}"))?;
        agg.reads += 1;
        if let Some(placement) = observation
            .placements
            .iter()
            .find(|placement| placement.member == resolved_hypothesis.representative)
        {
            agg.intervals.push((placement.start, placement.end));
            if placement.strand == '+' {
                agg.plus += 1;
            } else {
                agg.minus += 1;
            }
        }
    }
    Ok(RoleAggregation {
        member_to_representative,
        hypotheses: resolved,
        validation_unassigned: unassigned,
    })
}

struct CompositeEvidence {
    aggs: HashMap<String, ContigAgg>,
    details: HashMap<String, CompositeDetails>,
    decoy_background: DecoyBackground,
    target_validation_unassigned: u64,
    decoy_validation_unassigned: u64,
}

fn build_decoy_background(
    contigs: &[reference::ContigMeta],
    role_members: &RoleMembers,
    validation: &[DirectObservation],
) -> Result<DecoyBackground, String> {
    let mut background = DecoyBackground::default();
    let metadata = contigs
        .iter()
        .map(|contig| (contig.name.as_str(), contig))
        .collect::<HashMap<_, _>>();
    let mut member_strata = Vec::with_capacity(role_members.decoy_names.len());
    for name in role_members.names(Role::Decoy) {
        let contig = metadata
            .get(name.as_str())
            .ok_or_else(|| format!("decoy 成员缺少 metadata: {name}"))?;
        let stratum = strata(contig.len, contig.gc_frac);
        let layer = background.layers.entry(stratum.clone()).or_default();
        layer.reference_bases = layer
            .reference_bases
            .checked_add(contig.len)
            .ok_or_else(|| format!("decoy 层 {stratum} reference bases 溢出"))?;
        layer.contigs = layer
            .contigs
            .checked_add(1)
            .ok_or_else(|| format!("decoy 层 {stratum} contig 数溢出"))?;
        background.global.reference_bases = background
            .global
            .reference_bases
            .checked_add(contig.len)
            .ok_or_else(|| "全局 decoy reference bases 溢出".to_string())?;
        background.global.contigs = background
            .global
            .contigs
            .checked_add(1)
            .ok_or_else(|| "全局 decoy contig 数溢出".to_string())?;
        member_strata.push(stratum);
    }

    for observation in validation {
        let mut observed_strata = BTreeSet::new();
        for &member in &observation.members {
            let stratum = member_strata
                .get(member)
                .ok_or_else(|| format!("decoy validation 成员下标越界: {member}"))?;
            observed_strata.insert(stratum.as_str());
        }
        if observed_strata.len() > 1 {
            background.cross_stratum_reads += 1;
        }
        if !observed_strata.is_empty() {
            background.global.reads += 1;
        }
        for stratum in observed_strata {
            background
                .layers
                .get_mut(stratum)
                .ok_or_else(|| format!("decoy validation 指向未知层: {stratum}"))?
                .reads += 1;
        }
    }
    Ok(background)
}

fn build_composite_evidence(
    contigs: &[reference::ContigMeta],
    target_groups: &[index::TargetGroupMeta],
    role_members: &RoleMembers,
    audit: &AuditChunk,
) -> Result<CompositeEvidence, String> {
    let decoy_background = build_decoy_background(contigs, role_members, &audit.decoy_validation)?;
    let mut aggs = contigs
        .iter()
        .map(|contig| (contig.name.clone(), ContigAgg::default()))
        .collect::<HashMap<_, _>>();
    let target = aggregate_role_hypotheses(
        Role::Target,
        &audit.target_discovery,
        &audit.target_validation,
        role_members,
        &mut aggs,
    )?;
    let decoy = aggregate_role_hypotheses(
        Role::Decoy,
        &audit.decoy_discovery,
        &audit.decoy_validation,
        role_members,
        &mut aggs,
    )?;

    // split/discordant 只使用 validation 折，并归入包含其直接 target 成员的 OR 假设。
    for fragment in &audit.validation_fragments {
        for event in &fragment.split_events {
            let Some(&(Role::Target, member)) = role_members.lookup.get(&event.target_contig)
            else {
                continue;
            };
            let Some(representative) = target
                .member_to_representative
                .get(member)
                .and_then(|value| value.as_ref())
            else {
                continue;
            };
            aggs.get_mut(representative)
                .ok_or_else(|| format!("未知 split 代表 contig: {representative}"))?
                .split_events
                .push(cluster::SiteEvent {
                    contig: representative.clone(),
                    pos: event.target_pos as i64,
                    host_contig: event.host_contig.clone(),
                    host_pos: event.host_pos as i64,
                    direction: event.direction.clone(),
                    qname: event.qname.clone(),
                });
        }
        if let Some(contig) = &fragment.discordant {
            let Some(&(Role::Target, member)) = role_members.lookup.get(contig) else {
                continue;
            };
            if let Some(representative) = target
                .member_to_representative
                .get(member)
                .and_then(|value| value.as_ref())
            {
                aggs.get_mut(representative)
                    .ok_or_else(|| format!("未知 discordant 代表 contig: {representative}"))?
                    .discordant += 1;
            }
        }
    }

    let groups = target_groups
        .iter()
        .map(|group| (group.contig.as_str(), group))
        .collect::<HashMap<_, _>>();
    let metadata = contigs
        .iter()
        .map(|contig| (contig.name.as_str(), contig))
        .collect::<HashMap<_, _>>();
    let target_names = role_members.names(Role::Target);
    let mut details = HashMap::new();
    for resolved in target.hypotheses {
        let representative_contig = target_names
            .get(resolved.representative)
            .ok_or_else(|| "target 代表成员下标越界".to_string())?;
        let representative_group = groups
            .get(representative_contig.as_str())
            .ok_or_else(|| format!("target_groups 缺少 {representative_contig}"))?;
        let mut original_members = BTreeSet::new();
        let mut target_exposure_bases = 0u64;
        for &member in &resolved.hypothesis.members {
            let contig = target_names
                .get(member)
                .ok_or_else(|| format!("target 假设成员下标越界: {member}"))?;
            let contig_meta = metadata
                .get(contig.as_str())
                .ok_or_else(|| format!("target 假设成员缺少 metadata: {contig}"))?;
            target_exposure_bases = target_exposure_bases
                .checked_add(contig_meta.len)
                .ok_or_else(|| format!("target 假设 {representative_contig} exposure 溢出"))?;
            let group = groups
                .get(contig.as_str())
                .ok_or_else(|| format!("target_groups 缺少 {contig}"))?;
            original_members.extend(group.members.iter().cloned());
        }
        let explanation = resolved
            .hypothesis
            .explanation
            .iter()
            .map(|&member| {
                let contig = target_names
                    .get(member)
                    .ok_or_else(|| format!("target explanation 下标越界: {member}"))?;
                groups
                    .get(contig.as_str())
                    .map(|group| group.representative.clone())
                    .ok_or_else(|| format!("target_groups 缺少 {contig}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        details.insert(
            representative_contig.clone(),
            CompositeDetails {
                representative: representative_group.representative.clone(),
                members: original_members.into_iter().collect(),
                explanation,
                discovery_reads: resolved.hypothesis.discovery_observations,
                index_member_count: resolved.hypothesis.members.len(),
                target_exposure_bases,
            },
        );
    }

    Ok(CompositeEvidence {
        aggs,
        details,
        decoy_background,
        target_validation_unassigned: target.validation_unassigned,
        decoy_validation_unassigned: decoy.validation_unassigned,
    })
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

/// 代表序列按固定数量等分；每条直接 validation read 只按对齐中点占据一个位置
/// 分区。这样一条长 read 不能伪造三个区域，连续高覆盖也不会被合并成“一个窗口”。
fn distributed_window_count(intervals: &[(i32, i32)], contig_len: u64) -> u64 {
    if contig_len == 0 {
        return 0;
    }
    let mut occupied = BTreeSet::new();
    for &(start, end) in intervals {
        let start = i64::from(start).max(0) as u64;
        let end = (i64::from(end).max(0) as u64).min(contig_len);
        if start >= end {
            continue;
        }
        let midpoint = start + (end - start) / 2;
        let bin = ((u128::from(midpoint) * u128::from(DISTRIBUTED_WINDOW_BINS))
            / u128::from(contig_len))
        .min(u128::from(DISTRIBUTED_WINDOW_BINS - 1)) as u64;
        occupied.insert(bin);
    }
    occupied.len() as u64
}

struct CandidateDecision {
    candidates: Vec<report::Candidate>,
    test_family_size: usize,
}

/// discovery 固定检验族 + validation 精确两 Poisson 率检验 + BH adjusted-p +
/// coverage/distribution 报告门。
///
/// target 计数的 exposure 是 discovery OR 假设内全部内部索引 contig 的长度总和。
/// 单内部成员假设与同 size/GC 层的固定 synthetic-decoy 比较；多内部成员假设使用
/// 全局固定 decoy，双方 read-end 均只计一次、reference bases 均各计一次。精确
/// 条件检验在背景观测为 0 时仍保留有限样本不确定性。全部 discovery 假设（包括
/// validation 为 0）进入 BH，避免用 validation 再筛检验族。synthetic decoy
/// 可交换性尚未独立验证，因而 adjusted-p 是研究性模型证据，不宣称经典 FDR 控制。
///
/// 通用病毒检出不再依赖宿主-病毒 split；split/site 作为独立 integration evidence
/// 报告。`PASS` 仅表示候选通过当前探索性模型与分布门，不是样本级或临床结论。
fn decide_candidates(
    contigs: &[reference::ContigMeta],
    aggs: &HashMap<String, ContigAgg>,
    decoy_background: &DecoyBackground,
    total_read_sides: u64,
    composite_details: Option<&HashMap<String, CompositeDetails>>,
) -> Result<CandidateDecision, String> {
    let n_total = total_read_sides.max(1) as f64;
    let depth_norm = 1e7 / n_total;

    struct RawCandidate {
        contig: String,
        len: u64,
        index_member_count: usize,
        target_exposure_bases: u64,
        has_evidence: bool,
        bases: u64,
        windows: u64,
        frac: f64,
        reads: u64,
        split_events: Vec<cluster::SiteEvent>,
        discordant: u64,
        plus: u64,
        minus: u64,
        stratum: String,
        expected: f64,
        lambda_layer: f64,
        p: stats::TailProbability,
        rpm: f64,
        background: DecoyLayerBackground,
        background_scope: &'static str,
        background_status: &'static str,
        p_resolution_reference: f64,
    }

    let mut raw = Vec::new();
    for contig in contigs.iter().filter(|contig| contig.role == Role::Target) {
        // 生产路径只检验 discovery 定义的假设；测试传 None 时保留显式传入目标族。
        if composite_details.is_some_and(|details| !details.contains_key(&contig.name)) {
            continue;
        }
        let agg = &aggs[&contig.name];
        let details = composite_details.and_then(|details| details.get(&contig.name));
        let index_member_count = details
            .map(|details| details.index_member_count)
            .unwrap_or(1);
        let target_exposure_bases = details
            .map(|details| details.target_exposure_bases)
            .unwrap_or(contig.len);
        if index_member_count == 0 || target_exposure_bases == 0 {
            return Err(format!("target 假设 {} exposure 为空", contig.name));
        }
        let bases = merged_covered_bases(&agg.intervals);
        let has_evidence =
            agg.reads > 0 || bases > 0 || agg.discordant > 0 || !agg.split_events.is_empty();
        let stratum = strata(contig.len, contig.gc_frac);
        let (background, background_scope) = if index_member_count > 1 {
            if decoy_background.global.reference_bases > 0 {
                (decoy_background.global, "global_composite_hypothesis")
            } else {
                (DecoyLayerBackground::default(), "none")
            }
        } else {
            decoy_background
                .layers
                .get(&stratum)
                .filter(|layer| layer.reference_bases > 0)
                .copied()
                .map(|layer| (layer, "size_gc_stratum"))
                .or_else(|| {
                    (decoy_background.global.reference_bases > 0)
                        .then_some((decoy_background.global, "global_fallback"))
                })
                .unwrap_or((DecoyLayerBackground::default(), "none"))
        };
        let background_status = if background.reference_bases > 0 {
            "SYNTHETIC_DECOY_UNVALIDATED"
        } else {
            "UNAVAILABLE"
        };
        let p = if background.reference_bases > 0 {
            stats::exact_poisson_rate_upper_tail(
                agg.reads,
                target_exposure_bases as f64,
                background.reads,
                background.reference_bases as f64,
            )?
        } else {
            stats::TailProbability {
                probability: 1.0,
                ln_probability: 0.0,
                underflow: false,
            }
        };
        let expected = if background.reference_bases > 0 {
            background.reads as f64 * target_exposure_bases as f64
                / background.reference_bases as f64
        } else {
            0.0
        };
        let lambda_layer = if background.reference_bases > 0 {
            background.reads as f64 / background.reference_bases as f64 * depth_norm
        } else {
            0.0
        };
        let p_resolution_reference = if background.reference_bases > 0 {
            target_exposure_bases as f64
                / (target_exposure_bases as f64 + background.reference_bases as f64)
        } else {
            1.0
        };
        raw.push(RawCandidate {
            contig: contig.name.clone(),
            len: contig.len,
            index_member_count,
            target_exposure_bases,
            has_evidence,
            bases,
            windows: distributed_window_count(&agg.intervals, contig.len),
            frac: bases as f64 / contig.len as f64,
            reads: agg.reads,
            split_events: agg.split_events.clone(),
            discordant: agg.discordant,
            plus: agg.plus,
            minus: agg.minus,
            stratum,
            expected,
            lambda_layer,
            p,
            rpm: agg.reads as f64 / total_read_sides.max(1) as f64 * 1e6,
            background,
            background_scope,
            background_status,
            p_resolution_reference,
        });
    }

    // validation 计数为 0 的 discovery 假设也以 p=1 进入固定检验族。
    let ln_p_values = raw
        .iter()
        .map(|candidate| candidate.p.ln_probability)
        .collect::<Vec<_>>();
    let adjusted = stats::benjamini_hochberg_from_log(&ln_p_values)?;
    let family_size = raw.len();
    let mut candidates = Vec::new();
    for (index, candidate) in raw.into_iter().enumerate() {
        if !candidate.has_evidence {
            continue;
        }
        let q = adjusted[index];
        let mut decision_reasons = Vec::new();
        let decision = if candidate.background_status == "UNAVAILABLE" {
            decision_reasons.push("background_unavailable");
            "NOT_SIGNIFICANT"
        } else if q.probability >= MODEL_ADJUSTED_P_MAX {
            decision_reasons.push("adjusted_p_at_or_above_cutoff");
            "NOT_SIGNIFICANT"
        } else {
            if candidate.frac < COVERAGE_MIN {
                decision_reasons.push("coverage_below_min");
            }
            if candidate.windows < MIN_DISTRIBUTED_WINDOWS {
                decision_reasons.push("distributed_windows_below_min");
            }
            if decision_reasons.is_empty() {
                decision_reasons.push("model_and_distribution_gates_met");
                "PASS"
            } else {
                "BELOW_THRESHOLD"
            }
        };
        let sites = cluster::cluster_sites(&candidate.split_events, cluster::MIN_SITE_SUPPORT);
        let integration_evidence = if !sites.is_empty() {
            "SUPPORTED_SITE"
        } else if !candidate.split_events.is_empty() {
            "UNCLUSTERED_SPLIT"
        } else {
            "NONE"
        };
        let details = composite_details.and_then(|details| details.get(&candidate.contig));
        let poisson_p = Some(stats::poisson_upper_tail(
            candidate.reads,
            candidate.expected,
        ));
        candidates.push(report::Candidate {
            representative: details
                .map(|details| details.representative.clone())
                .unwrap_or_else(|| candidate.contig.clone()),
            hypothesis_members: details
                .map(|details| details.members.clone())
                .unwrap_or_else(|| vec![candidate.contig.clone()]),
            hypothesis_explanation: details
                .map(|details| details.explanation.clone())
                .unwrap_or_else(|| vec![candidate.contig.clone()]),
            discovery_reads: details.map(|details| details.discovery_reads).unwrap_or(0),
            contig: candidate.contig,
            contig_len: candidate.len,
            index_member_count: candidate.index_member_count,
            target_exposure_bases: candidate.target_exposure_bases,
            covered_bases: candidate.bases,
            covered_frac: candidate.frac,
            reads: candidate.reads,
            split_events: candidate.split_events.len() as u64,
            discordant: candidate.discordant,
            plus_strand: candidate.plus,
            minus_strand: candidate.minus,
            sites,
            p_value: candidate.p.probability,
            q_value: q.probability,
            ln_p_value: candidate.p.ln_probability,
            ln_q_value: q.ln_probability,
            p_underflow: candidate.p.underflow,
            q_underflow: q.underflow,
            p_resolution_floor: candidate.p_resolution_reference,
            poisson_p,
            nb_p: None,
            stratum: candidate.stratum,
            stratum_decoy_count: candidate.background.contigs,
            n_plain: candidate.reads,
            expected_hits: candidate.expected,
            lambda_bg_layer: candidate.lambda_layer,
            depth_fold: if candidate.expected > 0.0 {
                Some(candidate.reads as f64 / candidate.expected)
            } else {
                None
            },
            depth_rpm: candidate.rpm,
            p_floor_flag: candidate.p.probability <= candidate.p_resolution_reference,
            pi0: 1.0,
            decision,
            decision_reasons,
            integration_evidence,
            distinct_windows: candidate.windows,
            test_family_size: family_size,
            background_status: candidate.background_status,
            background_scope: candidate.background_scope,
            background_reads: candidate.background.reads,
            background_reference_bases: candidate.background.reference_bases,
            background_cross_stratum_reads: decoy_background.cross_stratum_reads,
        });
    }
    candidates.sort_by(|left, right| left.ln_q_value.total_cmp(&right.ln_q_value));
    Ok(CandidateDecision {
        candidates,
        test_family_size: family_size,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_decoy_background(
        contigs: &[reference::ContigMeta],
        aggs: &HashMap<String, ContigAgg>,
    ) -> DecoyBackground {
        let mut background = DecoyBackground::default();
        for contig in contigs.iter().filter(|contig| contig.role == Role::Decoy) {
            let reads = aggs[&contig.name].reads;
            let layer = background
                .layers
                .entry(strata(contig.len, contig.gc_frac))
                .or_default();
            layer.reads += reads;
            layer.reference_bases += contig.len;
            layer.contigs += 1;
            background.global.reads += reads;
            background.global.reference_bases += contig.len;
            background.global.contigs += 1;
        }
        background
    }

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
    fn fixed_decoy_validation_background_does_not_require_discovery_selection() {
        let contigs = vec![
            reference::ContigMeta {
                name: "decoy_short".into(),
                role: Role::Decoy,
                len: 3_000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_long".into(),
                role: Role::Decoy,
                len: 30_000,
                gc_frac: 0.5,
            },
        ];
        let members = RoleMembers::from_contigs(&contigs);
        // 一个直接 validation read 同时兼容两个跨层 decoy；全局只计一票，每个
        // 相关层各计一票作为该层目标检验的保守背景。
        let observations = vec![DirectObservation {
            members: vec![0, 1],
            placements: Vec::new(),
        }];
        let background = build_decoy_background(&contigs, &members, &observations).unwrap();
        assert_eq!(background.global.reads, 1);
        assert_eq!(background.global.reference_bases, 33_000);
        assert_eq!(background.global.contigs, 2);
        assert_eq!(background.cross_stratum_reads, 1);
        assert_eq!(background.layers["sz:1-5kb,gc:0.50-0.60"].reads, 1);
        assert_eq!(background.layers["sz:5-50kb,gc:0.50-0.60"].reads, 1);
    }

    #[test]
    fn composite_build_sums_internal_index_member_exposure() {
        let contigs = vec![
            reference::ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 1_000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "target_1".into(),
                role: Role::Target,
                len: 3_000,
                gc_frac: 0.5,
            },
        ];
        let target_groups = vec![
            index::TargetGroupMeta {
                contig: "target_0".into(),
                representative: "ACC0".into(),
                members: vec!["ACC0".into()],
            },
            index::TargetGroupMeta {
                contig: "target_1".into(),
                representative: "ACC1".into(),
                members: vec!["ACC1".into()],
            },
        ];
        let role_members = RoleMembers::from_contigs(&contigs);
        let audit = AuditChunk {
            target_discovery: vec![vec![0, 1]],
            target_validation: vec![DirectObservation {
                members: vec![0, 1],
                placements: vec![
                    DirectPlacement {
                        member: 0,
                        start: 0,
                        end: 100,
                        strand: '+',
                    },
                    DirectPlacement {
                        member: 1,
                        start: 0,
                        end: 100,
                        strand: '+',
                    },
                ],
            }],
            ..AuditChunk::default()
        };

        let composite =
            build_composite_evidence(&contigs, &target_groups, &role_members, &audit).unwrap();
        assert_eq!(composite.details.len(), 1);
        let details = composite.details.values().next().unwrap();
        assert_eq!(details.index_member_count, 2);
        assert_eq!(details.target_exposure_bases, 4_000);
        assert_eq!(details.members, ["ACC0", "ACC1"]);
        assert_eq!(
            composite
                .aggs
                .values()
                .map(|aggregate| aggregate.reads)
                .sum::<u64>(),
            1
        );
    }

    #[test]
    fn thread_budget_keeps_aux_decompression_outside_compute_budget() {
        assert_eq!(thread_budget(1), (0, 0));
        assert_eq!(thread_budget(8), (2, 7));
    }

    #[test]
    fn distributed_windows_use_alignment_midpoints() {
        assert_eq!(distributed_window_count(&[(0, 100), (20, 120)], 1000), 1);
        assert_eq!(
            distributed_window_count(&[(0, 100), (450, 550), (900, 1000)], 1000),
            3
        );
        assert_eq!(distributed_window_count(&[(0, 1000)], 1000), 1);
        assert_eq!(distributed_window_count(&[(100, 50)], 1000), 0);
        assert_eq!(distributed_window_count(&[], 1000), 0);
    }

    #[test]
    fn decide_chain_uses_exact_rate_test_and_fixed_family_bh() {
        // 2 目标各占一层、各 1 条零计数诱饵：
        // target_0 reads=10、target/decoy exposure 相等 → exact p=(1/2)^10；
        // 两个 discovery 假设都进入 BH，adjusted-p=2p。
        // validation 暴露 20_000 read-ends → end-RPM=500；三个分布窗口合计
        // 300/3000=10% → PASS；split 只形成独立 integration evidence。
        // target_1 validation 完全无证据，仍以 p=1 进入固定检验族，但不生成候选行。
        let contigs = vec![
            reference::ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 3000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "target_1".into(),
                role: Role::Target,
                len: 30000,
                gc_frac: 0.55,
            },
            reference::ContigMeta {
                name: "decoy_0".into(),
                role: Role::Decoy,
                len: 3000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_1".into(),
                role: Role::Decoy,
                len: 30000,
                gc_frac: 0.55,
            },
        ];
        let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
        for c in &contigs {
            aggs.insert(c.name.clone(), ContigAgg::default());
        }
        aggs.get_mut("target_0").unwrap().reads = 10;
        aggs.get_mut("target_0").unwrap().intervals = vec![(0, 100), (1000, 1100), (2000, 2100)];
        aggs.get_mut("target_0").unwrap().split_events = vec![cluster::SiteEvent {
            contig: "target_0".into(),
            pos: 1,
            host_contig: "host".into(),
            host_pos: 1,
            direction: "+".into(),
            qname: "a".into(),
        }];
        let details = ["target_0", "target_1"]
            .into_iter()
            .map(|name| {
                (
                    name.to_string(),
                    CompositeDetails {
                        representative: name.to_string(),
                        members: vec![name.to_string()],
                        explanation: vec![name.to_string()],
                        discovery_reads: 1,
                        index_member_count: 1,
                        target_exposure_bases: if name == "target_0" { 3_000 } else { 30_000 },
                    },
                )
            })
            .collect::<HashMap<_, _>>();
        let background = test_decoy_background(&contigs, &aggs);
        let cands = decide_candidates(&contigs, &aggs, &background, 20_000, Some(&details))
            .unwrap()
            .candidates;
        assert_eq!(cands.len(), 1);
        let t0 = &cands[0];
        assert_eq!(t0.contig, "target_0");
        assert_eq!(t0.n_plain, 10);
        assert!((t0.p_value - 0.5f64.powi(10)).abs() < 1e-12);
        assert!((t0.q_value - 2.0 * 0.5f64.powi(10)).abs() < 1e-12);
        assert_eq!(t0.decision, "PASS");
        assert_eq!(t0.integration_evidence, "UNCLUSTERED_SPLIT");
        assert!((t0.depth_rpm - 500.0).abs() < 1e-9, "rpm={}", t0.depth_rpm);
        assert!((t0.pi0 - 1.0).abs() < 1e-12);
        assert_eq!(t0.test_family_size, 2);
        assert_eq!(t0.stratum_decoy_count, 1);
        assert!((t0.p_resolution_floor - 0.5).abs() < 1e-9);
    }

    #[test]
    fn missing_decoy_background_fails_closed() {
        let contigs = vec![reference::ContigMeta {
            name: "target_0".into(),
            role: Role::Target,
            len: 3_000,
            gc_frac: 0.5,
        }];
        let mut aggs = HashMap::from([("target_0".into(), ContigAgg::default())]);
        aggs.get_mut("target_0").unwrap().reads = 10;
        aggs.get_mut("target_0").unwrap().intervals =
            vec![(0, 100), (1_000, 1_100), (2_000, 2_100)];

        let cands = decide_candidates(&contigs, &aggs, &DecoyBackground::default(), 20_000, None)
            .unwrap()
            .candidates;
        assert_eq!(cands.len(), 1);
        assert_eq!(cands[0].background_status, "UNAVAILABLE");
        assert_eq!(cands[0].p_value, 1.0);
        assert_eq!(cands[0].q_value, 1.0);
        assert_eq!(cands[0].decision, "NOT_SIGNIFICANT");
        assert_eq!(cands[0].decision_reasons, ["background_unavailable"]);
    }

    #[test]
    fn composite_hypothesis_uses_all_index_member_exposure_and_global_decoys() {
        let contigs = vec![
            reference::ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 1_000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "target_1".into(),
                role: Role::Target,
                len: 3_000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_short".into(),
                role: Role::Decoy,
                len: 1_000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_long".into(),
                role: Role::Decoy,
                len: 9_000,
                gc_frac: 0.5,
            },
        ];
        let mut aggs = contigs
            .iter()
            .map(|contig| (contig.name.clone(), ContigAgg::default()))
            .collect::<HashMap<_, _>>();
        aggs.get_mut("target_0").unwrap().reads = 2;
        aggs.get_mut("target_0").unwrap().intervals = vec![(0, 100), (400, 500), (800, 900)];
        aggs.get_mut("decoy_long").unwrap().reads = 1;
        let details = HashMap::from([(
            "target_0".into(),
            CompositeDetails {
                representative: "target_0".into(),
                members: vec!["target_0".into(), "target_1".into()],
                explanation: vec!["target_0".into()],
                discovery_reads: 2,
                index_member_count: 2,
                target_exposure_bases: 4_000,
            },
        )]);
        let background = test_decoy_background(&contigs, &aggs);

        let cands = decide_candidates(&contigs, &aggs, &background, 20_000, Some(&details))
            .unwrap()
            .candidates;
        assert_eq!(cands.len(), 1);
        let candidate = &cands[0];
        // Conditional total=3, target exposure probability=4000/(4000+10000)=2/7:
        // P[Binom(3, 2/7) >= 2] = 68/343.
        assert!((candidate.p_value - 68.0 / 343.0).abs() < 1e-12);
        assert_eq!(candidate.index_member_count, 2);
        assert_eq!(candidate.target_exposure_bases, 4_000);
        assert_eq!(candidate.background_scope, "global_composite_hypothesis");
        assert_eq!(candidate.background_reads, 1);
        assert_eq!(candidate.background_reference_bases, 10_000);
        assert_eq!(candidate.stratum_decoy_count, 2);
        assert!((candidate.expected_hits - 0.4).abs() < 1e-12);
    }

    #[test]
    fn split_evidence_does_not_remove_direct_detection_counts() {
        // reads=3、split qnames ["a","a","b"]；通用检测仍只按 3 个直接 read-end
        // 计数，split 不再从统计量中扣除或提升候选判定。
        let contigs = vec![
            reference::ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 3000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_0".into(),
                role: Role::Decoy,
                len: 3000,
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
        let background = test_decoy_background(&contigs, &aggs);
        let cands = decide_candidates(&contigs, &aggs, &background, 1_000_000, None)
            .unwrap()
            .candidates;
        assert_eq!(cands.len(), 1);
        assert_eq!(cands[0].n_plain, 3);
        assert_eq!(cands[0].split_events, 3);
    }

    #[test]
    fn low_end_rpm_is_not_a_hard_gate_but_breadth_and_windows_are() {
        // 两目标分层相同、诱饵零 reads；精确率检验仍给出有限 p。low_rpm 只有 5 个 read-end，
        // 但跨 3 个分区且 breadth=10%，因此不因未校准 RPM 被拒；narrow 即使深度
        // 很高，breadth=2% 仍 BELOW_THRESHOLD。
        let contigs = vec![
            reference::ContigMeta {
                name: "low_rpm".into(),
                role: Role::Target,
                len: 3000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "narrow".into(),
                role: Role::Target,
                len: 3000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_0".into(),
                role: Role::Decoy,
                len: 3000,
                gc_frac: 0.5,
            },
        ];
        let mut aggs: HashMap<String, ContigAgg> = HashMap::new();
        for c in &contigs {
            aggs.insert(c.name.clone(), ContigAgg::default());
        }
        aggs.get_mut("low_rpm").unwrap().reads = 5;
        aggs.get_mut("low_rpm").unwrap().intervals = vec![(0, 100), (1000, 1100), (2000, 2100)];
        aggs.get_mut("narrow").unwrap().reads = 10_000;
        aggs.get_mut("narrow").unwrap().intervals = vec![(0, 60)];
        let background = test_decoy_background(&contigs, &aggs);
        let cands = decide_candidates(&contigs, &aggs, &background, 2_000_000, None)
            .unwrap()
            .candidates;
        assert_eq!(cands.len(), 2);
        let low_rpm = cands.iter().find(|c| c.contig == "low_rpm").unwrap();
        let narrow = cands.iter().find(|c| c.contig == "narrow").unwrap();
        assert_eq!(low_rpm.decision, "PASS");
        assert_eq!(low_rpm.distinct_windows, 3);
        assert!((low_rpm.depth_rpm - 2.5).abs() < 1e-9);
        assert_eq!(narrow.decision, "BELOW_THRESHOLD");
        assert!((narrow.depth_rpm - 5_000.0).abs() < 1e-6);
    }

    #[test]
    fn one_validation_read_counts_once_for_a_multi_member_hypothesis() {
        let contigs = (0..3)
            .map(|index| reference::ContigMeta {
                name: format!("target_{index}"),
                role: Role::Target,
                len: 3000,
                gc_frac: 0.5,
            })
            .collect::<Vec<_>>();
        let role_members = RoleMembers::from_contigs(&contigs);
        let mut aggs = contigs
            .iter()
            .map(|contig| (contig.name.clone(), ContigAgg::default()))
            .collect::<HashMap<_, _>>();
        let discovery = vec![vec![0, 1], vec![1, 2]];
        let validation = vec![DirectObservation {
            members: vec![0, 1, 2],
            placements: vec![
                DirectPlacement {
                    member: 0,
                    start: 10,
                    end: 110,
                    strand: '+',
                },
                DirectPlacement {
                    member: 1,
                    start: 20,
                    end: 120,
                    strand: '-',
                },
                DirectPlacement {
                    member: 2,
                    start: 30,
                    end: 130,
                    strand: '+',
                },
            ],
        }];

        let aggregation = aggregate_role_hypotheses(
            Role::Target,
            &discovery,
            &validation,
            &role_members,
            &mut aggs,
        )
        .unwrap();

        assert_eq!(aggregation.hypotheses.len(), 1);
        assert_eq!(aggregation.hypotheses[0].representative, 1);
        assert_eq!(aggregation.validation_unassigned, 0);
        assert!(aggregation
            .member_to_representative
            .iter()
            .all(|representative| representative.as_deref() == Some("target_1")));
        assert_eq!(aggs["target_1"].reads, 1);
        assert_eq!(aggs["target_1"].minus, 1);
        assert_eq!(aggs["target_1"].intervals, vec![(20, 120)]);
        assert_eq!(aggs["target_0"].reads + aggs["target_2"].reads, 0);
    }

    #[test]
    fn target_and_decoy_validation_use_the_same_subset_rule() {
        let contigs = vec![
            reference::ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 1000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "target_1".into(),
                role: Role::Target,
                len: 1000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_0".into(),
                role: Role::Decoy,
                len: 1000,
                gc_frac: 0.5,
            },
            reference::ContigMeta {
                name: "decoy_1".into(),
                role: Role::Decoy,
                len: 1000,
                gc_frac: 0.5,
            },
        ];
        let role_members = RoleMembers::from_contigs(&contigs);
        let observation = DirectObservation {
            members: vec![1],
            placements: vec![DirectPlacement {
                member: 1,
                start: 0,
                end: 100,
                strand: '+',
            }],
        };
        for role in [Role::Target, Role::Decoy] {
            let mut aggs = contigs
                .iter()
                .map(|contig| (contig.name.clone(), ContigAgg::default()))
                .collect::<HashMap<_, _>>();
            let aggregation = aggregate_role_hypotheses(
                role,
                &[vec![0]],
                std::slice::from_ref(&observation),
                &role_members,
                &mut aggs,
            )
            .unwrap();
            assert_eq!(aggregation.validation_unassigned, 1, "role={role:?}");
            assert_eq!(
                aggs.values().map(|aggregate| aggregate.reads).sum::<u64>(),
                0,
                "role={role:?}"
            );
        }
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
