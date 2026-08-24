//! 竞争比对 + 候选证据：minimap2 sr 预设（Rust 绑定），fragment 级竞争裁决，
//! split-read / discordant 证据提取，以及有界队列 + 有序聚合的并行处理。
//! minimap2 通过 Rust 原生绑定调用，不启动外部 CLI。

use std::collections::{BTreeMap, HashMap};
use std::path::Path;
use std::sync::mpsc::{sync_channel, TrySendError};
use std::sync::Arc;
use std::thread;

use minimap2::{Aligner, Built, Mapping};

use crate::reference::Role;

/// 竞争裁决阈值。
pub const MIN_MAPQ: u32 = 20;
pub const MIN_AS_DIFF: i32 = 12;
pub const MAX_NM: i32 = 8;
/// split-read 阈值（ViFi 社区惯例）。
pub const SPLIT_SOFTCLIP: i32 = 20;
pub const SPLIT_MAPQ: u32 = 10;
/// minimap2 二级命中保留数量。
pub const BEST_N: i32 = 15;
/// 单 read ALL_CHAINS 审计最多接受的命中数。4096 相对 fast `BEST_N` 保留超过
/// 两个数量级的审计余量，同时把后续克隆、按 contig 去重和稳定排序限制在线性
/// 内存与可控 CPU 内；超限是显式审计状态，不能当作空证据。
pub const MAX_AUDIT_HITS: usize = 4096;
/// 全进程 kalloc 预算 1GB，按线程均分。
const TOTAL_KALLOC_BUDGET: i64 = 1_000_000_000;

pub struct CompetitiveAligner {
    aligner: Arc<Aligner<Built>>,
    /// 从 fast aligner cheap clone；底层 `idx`/`idx_parts` 由 `Arc` 共享。
    audit_aligner: Arc<Aligner<Built>>,
    /// 比对 worker 线程数 = N−1（N 为计算线程预算；低占用解压 helper 不从中
    /// 扣减，见 lib::thread_budget）；0 表示 `--threads 1`：走主线程内联路径。
    pub threads: usize,
}

impl CompetitiveAligner {
    /// `workers` 为比对 worker 线程数（= N−1，见 lib::thread_budget）；
    /// 0 表示 `--threads 1`：走主线程内联路径。
    pub fn open(mmi_path: &Path, workers: usize) -> Result<Self, &'static str> {
        let mut builder = Aligner::builder().sr().with_cigar();
        builder.mapopt.best_n = BEST_N;
        // 队列饱和时主线程也会 caller-runs，最大并发 mapper = workers + 1；
        // 内联路径 workers=0 时仍为 1。1GB kalloc 预算按该上界均分。
        let mapper_concurrency = workers.saturating_add(1);
        let mapper_concurrency = i64::try_from(mapper_concurrency).unwrap_or(i64::MAX);
        builder.mapopt.cap_kalloc = (TOTAL_KALLOC_BUDGET / mapper_concurrency).max(1);
        let aligner = builder.with_index(mmi_path, None)?;
        crate::reference::ensure_single_part_index(aligner.idx_parts.len())?;
        // `Aligner::clone` 只复制 mapopt 等轻量状态，索引由内部 Arc 共享；MMI 不再加载。
        let mut audit_aligner = aligner.clone();
        audit_aligner.mapopt.flag |= minimap2::ffi::MM_F_ALL_CHAINS as i64;
        Ok(Self {
            aligner: Arc::new(aligner),
            audit_aligner: Arc::new(audit_aligner),
            threads: workers,
        })
    }

    pub fn map_pair(
        &self,
        r1: &[u8],
        r2: &[u8],
        name: &str,
    ) -> Result<(Vec<Mapping>, Vec<Mapping>), &'static str> {
        self.aligner
            .map_pair(r1, r2, false, false, None, None, Some(name.as_bytes()))
    }

    /// 单端模式：只比对一条 read。
    pub fn map_seq(&self, seq: &[u8], name: &str) -> Result<Vec<Mapping>, &'static str> {
        self.aligner
            .map(seq, false, false, None, None, Some(name.as_bytes()))
    }

    /// 单 read 歧义审计：保留 minimap2 的全部 chains，供 [`adjudicate_audit`] 使用。
    /// 生产检测对 bottom-k 样本中通过 Bloom 的每条 read 直接调用此路径；普通
    /// `best_n=15` mapper 仅保留给公开的常规映射/证据工具，不参与主检测裁决。
    pub fn audit_map_seq(&self, seq: &[u8], name: &str) -> Result<Vec<Mapping>, &'static str> {
        self.audit_aligner
            .map(seq, false, false, None, None, Some(name.as_bytes()))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Hit {
    pub role: Role,
    pub contig: String,
    pub as_score: i32,
    pub mapq: u32,
    pub is_primary: bool,
    pub is_supplementary: bool,
    pub nm: i32,
    pub qstart: i32,
    pub qend: i32,
    pub qlen: i32,
    pub tstart: i32,
    pub tend: i32,
    pub strand: char,
}

/// 提取一条 read 的全部命中（含 supplementary，即 SA 样式的 split 另一半）。
pub fn hits_of(mappings: &[Mapping], roles: &HashMap<String, Role>) -> Vec<Hit> {
    mappings
        .iter()
        .filter_map(|m| {
            let contig = m.target_name.as_ref()?.to_string();
            let role = *roles.get(&contig)?;
            let aln = m.alignment.as_ref()?;
            Some(Hit {
                role,
                contig,
                as_score: aln.alignment_score.unwrap_or(0),
                mapq: m.mapq,
                is_primary: m.is_primary,
                is_supplementary: m.is_supplementary,
                nm: aln.nm,
                qstart: m.query_start,
                qend: m.query_end,
                qlen: m.query_len.map(|l| l.get()).unwrap_or(m.query_end.max(1)),
                tstart: m.target_start,
                tend: m.target_end,
                strand: if matches!(m.strand, minimap2::Strand::Forward) {
                    '+'
                } else {
                    '-'
                },
            })
        })
        .collect()
}

/// 单 read ALL_CHAINS 歧义审计结果。`Overflow` 与 `NoEvidence` 分离，防止资源
/// 超限被误当成成功形状的空证据。
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AuditDecision {
    /// 仅产出 Target 或 Decoy 的直接、精确最高 AS 命中。
    Resolved { role: Role, hits: Vec<Hit> },
    /// 审计完整执行，但没有角色满足保守胜出条件。
    NoEvidence,
    /// ALL_CHAINS 命中数超过单 read 资源上限，裁决未执行。
    Overflow { hit_count: usize, limit: usize },
}

fn audit_hit_is_better(candidate: &Hit, current: &Hit) -> bool {
    candidate
        .as_score
        .cmp(&current.as_score)
        .then_with(|| current.nm.cmp(&candidate.nm))
        .then_with(|| candidate.mapq.cmp(&current.mapq))
        .then_with(|| candidate.is_primary.cmp(&current.is_primary))
        .then_with(|| current.is_supplementary.cmp(&candidate.is_supplementary))
        .is_gt()
}

fn resolve_audit_role(hits: &[Hit], role: Role) -> Option<Vec<Hit>> {
    let role_best = hits
        .iter()
        .filter(|hit| hit.role == role)
        .map(|hit| hit.as_score)
        .max()?;
    let other_best = hits
        .iter()
        .filter(|hit| hit.role != role)
        .map(|hit| hit.as_score)
        .max()
        .unwrap_or(0);
    if i64::from(role_best) - i64::from(other_best) < i64::from(MIN_AS_DIFF) {
        return None;
    }

    // 第一版只接受角色内精确最高 AS ties；ΔAS>0 的兼容窗口延期到匹配校准后。
    let mut by_contig: BTreeMap<&str, Hit> = BTreeMap::new();
    for hit in hits
        .iter()
        .filter(|hit| hit.role == role && hit.nm <= MAX_NM && hit.as_score == role_best)
    {
        match by_contig.get_mut(hit.contig.as_str()) {
            Some(current) if audit_hit_is_better(hit, current) => *current = hit.clone(),
            Some(_) => {}
            None => {
                by_contig.insert(hit.contig.as_str(), hit.clone());
            }
        }
    }
    let mut resolved = by_contig.into_values().collect::<Vec<_>>();
    resolved.sort_by(|a, b| {
        b.as_score
            .cmp(&a.as_score)
            .then_with(|| a.contig.cmp(&b.contig))
            .then_with(|| a.nm.cmp(&b.nm))
            .then_with(|| b.mapq.cmp(&a.mapq))
    });
    (!resolved.is_empty()).then_some(resolved)
}

/// 对单 read 的 ALL_CHAINS 命中作保守审计裁决。
///
/// Target 必须以至少 [`MIN_AS_DIFF`] 严格跨过 Host/Contaminant/Decoy；Decoy
/// 对所有非 Decoy 角色应用同一门槛。角色内只返回 `AS == S_best` 且
/// `NM <= MAX_NM` 的直接命中，按 contig 去重；不维护静态 component，也不做
/// A-B-C 传递。`ΔAS>0` 的角色内兼容窗口明确延期到匹配校准后。
pub fn adjudicate_audit(hits: &[Hit]) -> AuditDecision {
    if hits.len() > MAX_AUDIT_HITS {
        return AuditDecision::Overflow {
            hit_count: hits.len(),
            limit: MAX_AUDIT_HITS,
        };
    }
    if let Some(hits) = resolve_audit_role(hits, Role::Target) {
        return AuditDecision::Resolved {
            role: Role::Target,
            hits,
        };
    }
    if let Some(hits) = resolve_audit_role(hits, Role::Decoy) {
        return AuditDecision::Resolved {
            role: Role::Decoy,
            hits,
        };
    }
    AuditDecision::NoEvidence
}

fn best_hit(hits: &[Hit], role: Role) -> Option<&Hit> {
    hits.iter()
        .filter(|h| h.role == role)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
}

fn best_non_target(hits: &[Hit]) -> Option<&Hit> {
    hits.iter()
        .filter(|h| h.role != Role::Target)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
}

/// 单 read 裁决结果。
#[derive(Debug, Clone, Default)]
pub struct ReadDecision {
    /// 高置信目标命中（mapq≥20、nm≤8、AS 差 ≥12）。
    pub confident_target: Option<Hit>,
    /// 最优宿主命中（不限 mapq，discordant 判定时再查阈值）。
    pub best_host: Option<Hit>,
    /// split 检测用：mapq ≥ SPLIT_MAPQ 的最优目标命中。
    pub split_target: Option<Hit>,
    /// split 检测用：mapq ≥ SPLIT_MAPQ 的最优宿主命中。
    pub split_host: Option<Hit>,
    /// 全局最优命中（mapq≥20，纯 AS 排序），覆盖度统计用。
    pub best_overall: Option<Hit>,
}

pub fn adjudicate_read(hits: &[Hit]) -> ReadDecision {
    let best_other = best_non_target(hits).cloned();
    let best_host = best_hit(hits, Role::Host).cloned();
    let best_overall = hits
        .iter()
        .filter(|h| h.mapq >= MIN_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
        .cloned();

    // 高置信目标：只在 mapq≥20 的目标命中里选最优（secondary 命中 mapq=0 会因 AS
    // 更高而挤掉 primary，故先按 mapq 过滤再比 AS），并做 AS 竞争差 + NM 检查。
    let confident_target = {
        let best_target_conf = hits
            .iter()
            .filter(|h| h.role == Role::Target && h.mapq >= MIN_MAPQ)
            .max_by(|a, b| a.as_score.cmp(&b.as_score).then(a.mapq.cmp(&b.mapq)))
            .cloned();
        best_target_conf.as_ref().and_then(|t| {
            let other_as = best_other.as_ref().map(|o| o.as_score).unwrap_or(0);
            if t.nm <= MAX_NM && t.as_score - other_as >= MIN_AS_DIFF {
                Some(t.clone())
            } else {
                None
            }
        })
    };

    let split_target = hits
        .iter()
        .filter(|h| h.role == Role::Target && h.mapq >= SPLIT_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score))
        .cloned();
    let split_host = hits
        .iter()
        .filter(|h| h.role == Role::Host && h.mapq >= SPLIT_MAPQ)
        .max_by(|a, b| a.as_score.cmp(&b.as_score))
        .cloned();

    ReadDecision {
        confident_target,
        best_host,
        split_target,
        split_host,
        best_overall,
    }
}

/// 普通 `best_n` 映射是否需要升级为单 read ALL_CHAINS 审计。
///
/// 触发条件集中为：primary 低 MAPQ、secondary 数达到 [`BEST_N`]（结果可能被
/// 截断），或普通结果出现 Target 却未产出 `confident_target`。明确、高 MAPQ、
/// 仅 Host 的正常路径不会无条件进入审计。
pub fn needs_ambiguity_audit(hits: &[Hit], ordinary: &ReadDecision) -> bool {
    let has_low_mapq_primary = hits.iter().any(|hit| hit.is_primary && hit.mapq < MIN_MAPQ);
    let secondary_saturated = hits
        .iter()
        .filter(|hit| !hit.is_primary && !hit.is_supplementary)
        .count()
        >= BEST_N as usize;
    let has_unresolved_target =
        ordinary.confident_target.is_none() && hits.iter().any(|hit| hit.role == Role::Target);
    has_low_mapq_primary || secondary_saturated || has_unresolved_target
}

/// split-read 断点事件。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplitEvent {
    pub target_contig: String,
    pub target_pos: i32,
    pub host_contig: String,
    pub host_pos: i32,
    pub direction: String,
    /// 来源 fragment 的 read ID（位点支持度按唯一 qname 计数）。
    pub qname: String,
}

/// 双侧 softclip 达到门槛的跨类 split read。
/// 取 AS 最高的 target/host 命中配对；候选位点支持度在聚类阶段按唯一 qname 计算。
fn split_event(read: &ReadDecision, qname: &str) -> Option<SplitEvent> {
    let t = read.split_target.as_ref()?;
    let h = read.split_host.as_ref()?;
    let t_left = t.qstart >= SPLIT_SOFTCLIP;
    let t_right = t.qlen - t.qend >= SPLIT_SOFTCLIP;
    let h_left = h.qstart >= SPLIT_SOFTCLIP;
    let h_right = h.qlen - h.qend >= SPLIT_SOFTCLIP;
    if !(t_left || t_right) || !(h_left || h_right) {
        return None;
    }
    let target_pos = if t_left { t.tstart } else { t.tend };
    let host_pos = if h_left { h.tstart } else { h.tend };
    Some(SplitEvent {
        target_contig: t.contig.clone(),
        target_pos,
        host_contig: h.contig.clone(),
        host_pos,
        direction: format!("host{}_{}", h.strand, t.strand),
        qname: qname.to_string(),
    })
}

/// fragment 级证据。
#[derive(Debug, Clone, Default)]
pub struct FragmentEvidence {
    pub r1: ReadDecision,
    pub r2: ReadDecision,
    pub split_events: Vec<SplitEvent>,
    /// 一端高置信目标、另一端高置信宿主（mapq≥20）的目标 contig 名。
    pub discordant: Option<String>,
}

pub fn adjudicate_fragment(key: &str, hits1: &[Hit], hits2: &[Hit]) -> FragmentEvidence {
    let r1 = adjudicate_read(hits1);
    let r2 = adjudicate_read(hits2);
    let mut split_events = Vec::new();
    split_events.extend(split_event(&r1, key));
    split_events.extend(split_event(&r2, key));
    let discordant = if r1.confident_target.is_some() && r2.confident_target.is_some() {
        None
    } else if let (Some(t), Some(h)) = (&r1.confident_target, &r2.best_host) {
        if h.mapq >= MIN_MAPQ {
            Some(t.contig.clone())
        } else {
            None
        }
    } else if let (Some(h), Some(t)) = (&r1.best_host, &r2.confident_target) {
        if h.mapq >= MIN_MAPQ {
            Some(t.contig.clone())
        } else {
            None
        }
    } else {
        None
    };
    FragmentEvidence {
        r1,
        r2,
        split_events,
        discordant,
    }
}

fn accept_ordered<R, G>(
    pending: &mut BTreeMap<usize, Vec<R>>,
    next: &mut usize,
    sink_err: &mut Option<String>,
    sink: &mut G,
    allow_sink: bool,
    idx: usize,
    results: Vec<R>,
) where
    G: FnMut(Vec<R>) -> Result<(), String>,
{
    pending.insert(idx, results);
    while let Some(ready) = pending.remove(next) {
        if allow_sink && sink_err.is_none() {
            if let Err(err) = sink(ready) {
                *sink_err = Some(err);
            }
        }
        *next += 1;
    }
}

/// 有界队列 + 有序聚合的并行处理。
/// 主线程投递 chunk，队列饱和时 caller-runs；结果按 chunk 序确定性回调 `sink`。
/// 输入流中的 Err 会中止处理并返回。
/// 两个使用方：预筛+竞争比对（融合单趟，收集证据）与证据聚合（局部聚合表有序合并）。
///
/// 并发契约（勿破坏）：
/// 1. worker 持有 tx_done 的克隆；主线程投递结束后 drop 自己的 tx_done，
///    最终 recv 才会以 Err 结束（否则永不返回）。
/// 2. 生产期每发一个 chunk 后机会式 try_recv 排空结果，避免两个有界通道互填挂死。
/// 3. caller-runs 结果也必须经 pending/next 有序提交，不得直接调用 sink。
pub(crate) fn par_map_chunks<T, R, F, G>(
    items: impl Iterator<Item = Result<T, String>>,
    workers: usize,
    chunk_size: usize,
    process: F,
    mut sink: G,
) -> Result<(), String>
where
    T: Send,
    R: Send,
    F: Fn(&[T]) -> Vec<R> + Sync + Send,
    G: FnMut(Vec<R>) -> Result<(), String>,
{
    let chunk_size = chunk_size.max(1);
    // workers == 0（--threads 1）：不创建 worker 线程，主线程内联处理所有 chunk，
    // 保证"总线程数 = N"的语义在 N=1 时同样成立。
    if workers == 0 {
        let mut first_err: Option<String> = None;
        let mut chunk: Vec<T> = Vec::with_capacity(chunk_size);
        for item in items {
            match item {
                Ok(t) => {
                    chunk.push(t);
                    if chunk.len() >= chunk_size {
                        let results = process(&std::mem::take(&mut chunk));
                        if let Err(e) = sink(results) {
                            first_err = Some(e);
                            break;
                        }
                    }
                }
                Err(e) => {
                    first_err = Some(e);
                    break;
                }
            }
        }
        if !chunk.is_empty() && first_err.is_none() {
            if let Err(e) = sink(process(&chunk)) {
                first_err = Some(e);
            }
        }
        return if let Some(e) = first_err {
            Err(e)
        } else {
            Ok(())
        };
    }
    let (tx_chunk, rx_chunk) = sync_channel::<(usize, Vec<T>)>(workers * 2);
    let (tx_done, rx_done) = sync_channel::<(usize, Vec<R>)>(workers * 2);
    let rx_chunk = std::sync::Mutex::new(rx_chunk);
    let mut first_err: Option<String> = None;

    thread::scope(|scope| {
        for _ in 0..workers {
            let rx = &rx_chunk;
            // 每个 worker 持有 Sender 克隆；主线程结束后 drop 自己的 Sender，
            // 通道才能在对侧全部退出后关闭（否则最终 recv 永不返回 Err）。
            let tx = tx_done.clone();
            let process = &process;
            scope.spawn(move || loop {
                let item = {
                    let guard = match rx.lock() {
                        Ok(g) => g,
                        Err(_) => return,
                    };
                    guard.recv()
                };
                let (idx, chunk) = match item {
                    Ok(x) => x,
                    Err(_) => return,
                };
                let results: Vec<R> = process(&chunk);
                if tx.send((idx, results)).is_err() {
                    return;
                }
            });
        }

        let mut idx = 0usize;
        let mut chunk: Vec<T> = Vec::with_capacity(chunk_size);
        let mut pending: BTreeMap<usize, Vec<R>> = BTreeMap::new();
        let mut next = 0usize;
        let mut sink_err: Option<String> = None;
        // 机会式排空结果通道：生产期也要消费结果，避免两个有界通道互填死锁。
        let drain = |pending: &mut BTreeMap<usize, Vec<R>>,
                     next: &mut usize,
                     sink_err: &mut Option<String>,
                     sink: &mut G| {
            while let Ok((i, results)) = rx_done.try_recv() {
                accept_ordered(pending, next, sink_err, sink, true, i, results);
            }
        };
        let submit = |i: usize,
                      chunk: Vec<T>,
                      pending: &mut BTreeMap<usize, Vec<R>>,
                      next: &mut usize,
                      sink_err: &mut Option<String>,
                      sink: &mut G| {
            match tx_chunk.try_send((i, chunk)) {
                Ok(()) => Ok(()),
                Err(TrySendError::Full((i, chunk))) => {
                    let results = process(&chunk);
                    accept_ordered(pending, next, sink_err, sink, true, i, results);
                    Ok(())
                }
                Err(TrySendError::Disconnected(_)) => Err("处理 worker 意外退出".to_string()),
            }
        };
        for item in items {
            match item {
                Ok(t) => {
                    chunk.push(t);
                    if chunk.len() >= chunk_size {
                        // 队列满时由主线程处理当前 chunk；仍经 pending 按序提交。
                        drain(&mut pending, &mut next, &mut sink_err, &mut sink);
                        let full = std::mem::replace(&mut chunk, Vec::with_capacity(chunk_size));
                        if let Err(err) =
                            submit(idx, full, &mut pending, &mut next, &mut sink_err, &mut sink)
                        {
                            first_err = Some(err);
                            break;
                        }
                        idx += 1;
                        drain(&mut pending, &mut next, &mut sink_err, &mut sink);
                    }
                }
                Err(e) => {
                    first_err = Some(e);
                    break;
                }
            }
        }
        if !chunk.is_empty() && first_err.is_none() {
            drain(&mut pending, &mut next, &mut sink_err, &mut sink);
            if let Err(err) = submit(
                idx,
                chunk,
                &mut pending,
                &mut next,
                &mut sink_err,
                &mut sink,
            ) {
                first_err = Some(err);
            }
            drain(&mut pending, &mut next, &mut sink_err, &mut sink);
        }
        drop(tx_chunk);
        // 主线程释放自己的 Sender：worker 全部退出后 rx_done 才会关闭，
        // 最终 recv 循环才能以 Err 结束。
        drop(tx_done);
        // 最终全量排空（worker 全部退出后 tx_done 关闭，recv 返回 Err 结束）。
        // 无论是否有错都必须消费完结果，否则 worker 会阻塞在 tx_done.send 上。
        while let Ok((i, results)) = rx_done.recv() {
            accept_ordered(
                &mut pending,
                &mut next,
                &mut sink_err,
                &mut sink,
                first_err.is_none(),
                i,
                results,
            );
        }
        if let Some(e) = first_err.take() {
            Err(e)
        } else if let Some(e) = sink_err.take() {
            Err(e)
        } else {
            Ok(())
        }
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn hit(
        role: Role,
        as_score: i32,
        mapq: u32,
        qstart: i32,
        qend: i32,
        tstart: i32,
        tend: i32,
    ) -> Hit {
        Hit {
            role,
            contig: format!("{}_{}", role.prefix(), 0),
            as_score,
            mapq,
            is_primary: true,
            is_supplementary: false,
            nm: 0,
            qstart,
            qend,
            qlen: 150,
            tstart,
            tend,
            strand: '+',
        }
    }

    fn named_hit(role: Role, contig: &str, as_score: i32, mapq: u32, nm: i32) -> Hit {
        let mut hit = hit(role, as_score, mapq, 0, 150, 0, 150);
        hit.contig = contig.to_string();
        hit.nm = nm;
        hit
    }

    #[test]
    fn adjudicates_confident_target_with_as_margin() {
        let hits = vec![
            hit(Role::Target, 300, 60, 0, 150, 0, 150),
            hit(Role::Host, 280, 60, 0, 150, 0, 150), // AS 差 20 ≥ 12
        ];
        let d = adjudicate_read(&hits);
        assert!(d.confident_target.is_some());
    }

    #[test]
    fn rejects_target_without_as_margin() {
        let hits = vec![
            hit(Role::Target, 300, 60, 0, 150, 0, 150),
            hit(Role::Host, 295, 60, 0, 150, 0, 150), // AS 差 5 < 12
        ];
        let d = adjudicate_read(&hits);
        assert!(d.confident_target.is_none());
    }

    #[test]
    fn audit_recovers_compatible_mapq_zero_targets_over_weak_non_target() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 300, 0, 1),
            named_hit(Role::Host, "host_0", 280, 0, 0),
        ];

        match adjudicate_audit(&hits) {
            AuditDecision::Resolved { role, hits } => {
                assert_eq!(role, Role::Target);
                assert_eq!(
                    hits.iter()
                        .map(|hit| hit.contig.as_str())
                        .collect::<Vec<_>>(),
                    ["target_a", "target_b"]
                );
            }
            decision => panic!("目标应从 ALL_CHAINS 结果中恢复，实际为 {decision:?}"),
        }
    }

    #[test]
    fn audit_excludes_non_top_target_until_delta_as_is_calibrated() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 299, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&hits) else {
            panic!("最优目标应恢复");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(hits.len(), 1, "ΔAS>0 窗口须延期到匹配校准后");
        assert_eq!(hits[0].contig, "target_a");
    }

    #[test]
    fn audit_tied_host_or_near_contaminant_blocks_target() {
        for (competitor, competitor_score) in [(Role::Host, 300), (Role::Contaminant, 289)] {
            let hits = vec![
                named_hit(Role::Target, "target_a", 300, 0, 0),
                named_hit(competitor, "non_target", competitor_score, 0, 0),
            ];
            assert!(matches!(adjudicate_audit(&hits), AuditDecision::NoEvidence));
        }
    }

    #[test]
    fn audit_returns_only_targets_hit_directly_by_this_read() {
        let prior_read = vec![
            named_hit(Role::Target, "target_b", 300, 0, 0),
            named_hit(Role::Target, "target_c", 300, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        assert!(matches!(
            adjudicate_audit(&prior_read),
            AuditDecision::Resolved { .. }
        ));

        let this_read = vec![
            named_hit(Role::Target, "target_a", 300, 0, 0),
            named_hit(Role::Target, "target_b", 300, 0, 0),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&this_read) else {
            panic!("A/B 直接命中应恢复");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(
            hits.iter()
                .map(|hit| hit.contig.as_str())
                .collect::<Vec<_>>(),
            ["target_a", "target_b"]
        );
    }

    #[test]
    fn audit_deduplicates_contigs_and_keeps_the_best_hit() {
        let hits = vec![
            named_hit(Role::Target, "target_a", 294, 0, 0),
            named_hit(Role::Target, "target_a", 300, 0, 2),
            named_hit(Role::Target, "target_b", 300, 0, 1),
            named_hit(Role::Host, "host_0", 270, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&hits) else {
            panic!("目标直接命中应恢复");
        };
        assert_eq!(role, Role::Target);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0].contig, "target_a");
        assert_eq!(hits[0].as_score, 300);
        assert_eq!(hits[1].contig, "target_b");
    }

    #[test]
    fn audit_decoy_resolution_is_symmetric_and_conservative() {
        let winning = vec![
            named_hit(Role::Decoy, "decoy_a", 300, 0, 0),
            named_hit(Role::Decoy, "decoy_b", 300, 0, 1),
            named_hit(Role::Decoy, "decoy_c", 299, 0, 0),
            named_hit(Role::Host, "host_0", 280, 0, 0),
            named_hit(Role::Target, "target_a", 279, 0, 0),
        ];
        let AuditDecision::Resolved { role, hits } = adjudicate_audit(&winning) else {
            panic!("明确胜出的诱饵应产出同构证据");
        };
        assert_eq!(role, Role::Decoy);
        assert_eq!(
            hits.iter()
                .map(|hit| hit.contig.as_str())
                .collect::<Vec<_>>(),
            ["decoy_a", "decoy_b"],
            "Decoy 的 ΔAS>0 窗口同样须延期到匹配校准后"
        );

        let near_tie = vec![
            named_hit(Role::Decoy, "decoy_a", 300, 0, 0),
            named_hit(Role::Contaminant, "contam_0", 289, 0, 0),
        ];
        assert!(matches!(
            adjudicate_audit(&near_tie),
            AuditDecision::NoEvidence
        ));
    }

    #[test]
    fn audit_overflow_is_explicit() {
        let hits = (0..=MAX_AUDIT_HITS)
            .map(|i| named_hit(Role::Target, &format!("target_{i}"), 300, 0, 0))
            .collect::<Vec<_>>();
        assert_eq!(
            adjudicate_audit(&hits),
            AuditDecision::Overflow {
                hit_count: MAX_AUDIT_HITS + 1,
                limit: MAX_AUDIT_HITS,
            }
        );
    }

    #[test]
    fn best_n_saturation_requests_audit_before_a_hidden_host_can_be_missed() {
        let mut ordinary_hits = vec![named_hit(Role::Target, "target_primary", 300, 60, 0)];
        ordinary_hits[0].is_primary = true;
        for i in 0..BEST_N {
            let mut secondary = named_hit(
                Role::Target,
                &format!("target_secondary_{i}"),
                299 - i,
                0,
                0,
            );
            secondary.is_primary = false;
            ordinary_hits.push(secondary);
        }
        let ordinary = adjudicate_read(&ordinary_hits);
        assert!(ordinary.confident_target.is_some());
        assert!(needs_ambiguity_audit(&ordinary_hits, &ordinary));

        let mut all_chains = ordinary_hits;
        all_chains.push(named_hit(Role::Target, "target_after_best_n", 280, 0, 0));
        all_chains.push(named_hit(Role::Host, "host_after_best_n", 295, 0, 0));
        assert!(matches!(
            adjudicate_audit(&all_chains),
            AuditDecision::NoEvidence
        ));
    }

    #[test]
    fn audit_trigger_covers_low_mapq_primary_and_unresolved_target_but_not_clear_host() {
        let mut low_mapq_primary = named_hit(Role::Host, "host_0", 300, MIN_MAPQ - 1, 0);
        low_mapq_primary.is_primary = true;
        let hits = vec![low_mapq_primary];
        assert!(needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));

        let mut unresolved_target = named_hit(Role::Target, "target_a", 300, 0, 0);
        unresolved_target.is_primary = false;
        let hits = vec![unresolved_target];
        assert!(needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));

        let mut clear_host = named_hit(Role::Host, "host_0", 300, 60, 0);
        clear_host.is_primary = true;
        let hits = vec![clear_host];
        assert!(!needs_ambiguity_audit(&hits, &adjudicate_read(&hits)));
    }

    #[test]
    fn hits_preserve_primary_and_supplementary_flags() {
        let mapping = Mapping {
            target_name: Some(Arc::new("target_a".to_string())),
            is_primary: false,
            is_supplementary: true,
            alignment: Some(minimap2::Alignment {
                nm: 0,
                cigar: None,
                cigar_str: None,
                md: None,
                cs: None,
                alignment_score: Some(300),
            }),
            ..Mapping::default()
        };
        let roles = HashMap::from([("target_a".to_string(), Role::Target)]);
        let hits = hits_of(&[mapping], &roles);
        assert_eq!(hits.len(), 1);
        assert!(!hits[0].is_primary);
        assert!(hits[0].is_supplementary);
    }

    #[test]
    fn split_event_requires_both_side_softclips() {
        let hits1 = vec![
            hit(Role::Target, 200, 40, 0, 100, 0, 100),   // 右 clip 50
            hit(Role::Host, 180, 40, 100, 150, 500, 550), // 左 clip 100
        ];
        let hits2 = vec![
            hit(Role::Target, 200, 40, 0, 150, 0, 150), // 无 clip
            hit(Role::Host, 180, 40, 100, 150, 500, 550),
        ];
        let r1 = adjudicate_read(&hits1);
        let r2 = adjudicate_read(&hits2);
        assert!(split_event(&r1, "q1").is_some());
        assert!(split_event(&r2, "q2").is_none());
    }

    #[test]
    fn parallel_map_is_deterministic() {
        let items: Vec<Result<u64, String>> = (0..1000).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            4,
            64,
            |chunk| chunk.iter().map(|&v| v * 2).collect(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(out.len(), 1000);
        for (i, v) in out.iter().enumerate() {
            assert_eq!(*v, (i * 2) as u64);
        }
    }

    #[test]
    fn saturated_queue_uses_caller_and_preserves_order() {
        let caller = std::thread::current().id();
        let caller_used = std::sync::atomic::AtomicBool::new(false);
        let items: Vec<Result<u64, String>> = (0..32).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            1,
            1,
            |chunk| {
                let value = chunk[0];
                if std::thread::current().id() == caller {
                    caller_used.store(true, std::sync::atomic::Ordering::Release);
                } else if value == 0 {
                    // 固定阻塞首个 worker chunk，迫使两槽任务队列饱和并触发 caller-runs。
                    let deadline = std::time::Instant::now() + std::time::Duration::from_secs(1);
                    while !caller_used.load(std::sync::atomic::Ordering::Acquire) {
                        if std::time::Instant::now() >= deadline {
                            break;
                        }
                        std::thread::yield_now();
                    }
                }
                vec![value]
            },
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert!(caller_used.load(std::sync::atomic::Ordering::Acquire));
        assert_eq!(out, (0..32).collect::<Vec<_>>());
    }

    #[test]
    fn parallel_map_propagates_stream_error() {
        let mut items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        items.push(Err("流错误".into()));
        let mut out = Vec::new();
        let result = par_map_chunks(
            items.into_iter(),
            2,
            4,
            |chunk| chunk.to_vec(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        );
        assert_eq!(result.unwrap_err(), "流错误");
    }

    #[test]
    fn parallel_map_propagates_sink_error() {
        let items: Vec<Result<u64, String>> = (0..100).map(Ok).collect();
        let mut calls = 0u64;
        let result = par_map_chunks(
            items.into_iter(),
            4,
            8,
            |chunk| chunk.to_vec(),
            |rs| {
                calls += rs.len() as u64;
                if calls > 30 {
                    Err("sink 失败".into())
                } else {
                    Ok(())
                }
            },
        );
        assert_eq!(result.unwrap_err(), "sink 失败");
        assert_eq!(calls, 32, "sink 失败后不得继续回调后续 chunk");
    }

    #[test]
    fn inline_map_zero_workers_is_deterministic() {
        // --threads 1 语义：0 worker，主线程内联处理，顺序与并行路径一致。
        let items: Vec<Result<u64, String>> = (0..1000).map(Ok).collect();
        let mut out = Vec::new();
        par_map_chunks(
            items.into_iter(),
            0,
            64,
            |chunk| chunk.iter().map(|&v| v * 2).collect(),
            |rs| {
                out.extend(rs);
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(out.len(), 1000);
        for (i, v) in out.iter().enumerate() {
            assert_eq!(*v, (i * 2) as u64);
        }
    }

    #[test]
    fn inline_map_propagates_errors() {
        let items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        assert!(
            par_map_chunks(items.into_iter(), 0, 4, |chunk| chunk.to_vec(), |_| Ok(()),).is_ok()
        );
        let mut items: Vec<Result<u64, String>> = (0..10).map(Ok).collect();
        items.push(Err("流错误".into()));
        assert!(
            par_map_chunks(items.into_iter(), 0, 4, |chunk| chunk.to_vec(), |_| Ok(()),).is_err()
        );
        let items: Vec<Result<u64, String>> = (0..100).map(Ok).collect();
        let mut calls = 0u64;
        assert!(par_map_chunks(
            items.into_iter(),
            0,
            8,
            |chunk| chunk.to_vec(),
            |rs| {
                calls += rs.len() as u64;
                if calls > 30 {
                    Err("sink 失败".into())
                } else {
                    Ok(())
                }
            },
        )
        .is_err());
    }

    #[test]
    fn open_rejects_multi_part_index() {
        let dir = std::env::temp_dir().join(format!("vf_multipart_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("multi.fa");
        let mmi = dir.join("multi.mmi");
        let mut file = std::fs::File::create(&fasta).unwrap();
        let seq = "ACGT".repeat(200);
        for i in 0..4 {
            writeln!(file, ">contig_{i}\n{seq}").unwrap();
        }
        drop(file);

        let mut builder = Aligner::builder().sr().with_index_threads(1);
        builder.idxopt.mini_batch_size = 1;
        builder.idxopt.batch_size = 1;
        let indexed = builder
            .with_index(&fasta, Some(mmi.to_str().unwrap()))
            .unwrap();
        assert!(indexed.idx_parts.len() > 1);
        drop(indexed);

        let err = match CompetitiveAligner::open(&mmi, 1) {
            Ok(_) => panic!("多分片索引不应被接受"),
            Err(err) => err,
        };
        assert!(err.contains("单分片"), "err={err}");
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn open_shares_one_index_with_the_all_chains_audit_aligner() {
        let dir = std::env::temp_dir().join(format!("vf_audit_aligner_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let fasta = dir.join("single.fa");
        let mmi = dir.join("single.mmi");
        let mut file = std::fs::File::create(&fasta).unwrap();
        writeln!(file, ">target_0\n{}", "ACGT".repeat(200)).unwrap();
        drop(file);

        let indexed = Aligner::builder()
            .sr()
            .with_index_threads(1)
            .with_index(&fasta, Some(mmi.to_str().unwrap()))
            .unwrap();
        assert_eq!(indexed.idx_parts.len(), 1);
        drop(indexed);

        let aligners = CompetitiveAligner::open(&mmi, 1).unwrap();
        let all_chains = minimap2::ffi::MM_F_ALL_CHAINS as i64;
        assert_eq!(aligners.aligner.mapopt.best_n, BEST_N);
        assert_eq!(aligners.aligner.mapopt.flag & all_chains, 0);
        assert_ne!(aligners.audit_aligner.mapopt.flag & all_chains, 0);
        assert!(Arc::ptr_eq(
            &aligners.aligner.idx_parts[0],
            &aligners.audit_aligner.idx_parts[0]
        ));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
