//! 竞争比对 + 候选证据：minimap2 sr 预设（Rust 绑定），fragment 级竞争裁决，
//! split-read / discordant 证据提取，以及有界队列 + 有序聚合的并行处理。
//! minimap2 通过 Rust 原生绑定调用，不启动外部 CLI。

use std::collections::HashMap;
use std::path::Path;
use std::sync::mpsc::sync_channel;
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
/// 全进程 kalloc 预算 1GB，按线程均分。
const TOTAL_KALLOC_BUDGET: i64 = 1_000_000_000;

pub struct CompetitiveAligner {
    aligner: Arc<Aligner<Built>>,
    /// 比对 worker 线程数 = N−1−D（N 总线程、D 解压线程，见 lib::thread_budget）；
    /// 0 表示 `--threads 1`：走主线程内联路径，进程保持单线程。
    pub threads: usize,
}

impl CompetitiveAligner {
    /// `workers` 为比对 worker 线程数（= N−1−D，见 lib::thread_budget）；
    /// 0 表示 `--threads 1`：走主线程内联路径，进程保持单线程。
    pub fn open(mmi_path: &Path, workers: usize) -> Result<Self, &'static str> {
        let mut builder = Aligner::builder().sr().with_cigar();
        builder.mapopt.best_n = BEST_N;
        // 并发 mapper 数 = workers（内联时为 1），1GB 预算按并发 mapper 均分。
        builder.mapopt.cap_kalloc = (TOTAL_KALLOC_BUDGET / workers.max(1) as i64).max(1);
        let aligner = builder.with_index(mmi_path, None)?;
        Ok(Self {
            aligner: Arc::new(aligner),
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
}

#[derive(Debug, Clone)]
pub struct Hit {
    pub role: Role,
    pub contig: String,
    pub as_score: i32,
    pub mapq: u32,
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

/// 有界队列 + 有序聚合的并行处理。
/// 主线程投递 chunk，N 个 worker 并行处理，结果按 chunk 序确定性回调 `sink`；
/// 输入流中的 Err 会中止处理并返回。
/// 两个使用方：预筛+竞争比对（融合单趟，收集证据）与证据聚合（局部聚合表有序合并）。
///
/// 并发契约（勿破坏）：
/// 1. worker 持有 tx_done 的克隆；主线程投递结束后 drop 自己的 tx_done，
///    最终 recv 才会以 Err 结束（否则永不返回）。
/// 2. 生产期每发一个 chunk 后机会式 try_recv 排空结果，避免两个有界通道互填挂死。
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
        let mut pending: std::collections::BTreeMap<usize, Vec<R>> =
            std::collections::BTreeMap::new();
        let mut next = 0usize;
        let mut sink_err: Option<String> = None;
        // 机会式排空结果通道：生产期也要消费结果，避免两个有界通道互填死锁。
        let drain = |pending: &mut std::collections::BTreeMap<usize, Vec<R>>,
                     next: &mut usize,
                     sink_err: &mut Option<String>,
                     sink: &mut G| {
            while let Ok((i, results)) = rx_done.try_recv() {
                pending.insert(i, results);
                while let Some(v) = pending.remove(next) {
                    if sink_err.is_none() {
                        if let Err(e) = sink(v) {
                            *sink_err = Some(e);
                        }
                    }
                    *next += 1;
                }
            }
        };
        for item in items {
            match item {
                Ok(t) => {
                    chunk.push(t);
                    if chunk.len() >= chunk_size {
                        // send 前先排空结果：主线程在 send 上阻塞时（chunk 通道满）
                        // 不排空会让结果通道积满、worker 停在结果 send 上，损失吞吐。
                        drain(&mut pending, &mut next, &mut sink_err, &mut sink);
                        if tx_chunk.send((idx, std::mem::take(&mut chunk))).is_err() {
                            first_err = Some("处理 worker 意外退出".to_string());
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
        if !chunk.is_empty() && first_err.is_none() && tx_chunk.send((idx, chunk)).is_err() {
            first_err = Some("处理 worker 意外退出".to_string());
        }
        drop(tx_chunk);
        // 主线程释放自己的 Sender：worker 全部退出后 rx_done 才会关闭，
        // 最终 recv 循环才能以 Err 结束。
        drop(tx_done);
        // 最终全量排空（worker 全部退出后 tx_done 关闭，recv 返回 Err 结束）。
        // 无论是否有错都必须消费完结果，否则 worker 会阻塞在 tx_done.send 上。
        while let Ok((i, results)) = rx_done.recv() {
            pending.insert(i, results);
            while let Some(v) = pending.remove(&next) {
                if first_err.is_none() && sink_err.is_none() {
                    if let Err(e) = sink(v) {
                        sink_err = Some(e);
                    }
                }
                next += 1;
            }
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
            nm: 0,
            qstart,
            qend,
            qlen: 150,
            tstart,
            tend,
            strand: '+',
        }
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
        assert!(result.is_err());
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
        assert!(result.is_err());
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
}
