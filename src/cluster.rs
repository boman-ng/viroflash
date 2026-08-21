//! 候选聚类：位点 2bp 锚点聚类（唯一 qname 支持度 ≥2）→ 断裂点 5bp 锚点去重（均值合并）。
//! 聚类、去重和最小支持度阈值集中定义在本模块。

use std::collections::HashSet;

/// 位点聚类容差（bp）。
pub const CLUSTER_TOLERANCE: i64 = 2;
/// 断裂点去重容差（bp）。
pub const DEDUP_TOLERANCE: i64 = 5;
/// 出位点的最小支持 read（唯一 qname）数。
pub const MIN_SITE_SUPPORT: usize = 2;

fn mean_round(values: &[i64]) -> i64 {
    let sum: i64 = values.iter().sum();
    let n = values.len() as i64;
    (sum as f64 / n as f64).round() as i64
}

/// 断裂点事件：目标 contig 上的一个 split-read 断点。
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SiteEvent {
    pub contig: String,
    pub pos: i64,
    pub host_contig: String,
    pub host_pos: i64,
    pub direction: String,
    /// 来源 read ID：支持度按唯一 qname 计数。
    pub qname: String,
}

/// 聚类后的整合位点。
#[derive(Debug, Clone, PartialEq)]
pub struct Site {
    pub contig: String,
    pub pos: i64,
    /// 唯一 qname 支持数（5bp 去重合并时相加）。
    pub support: usize,
    pub host_contig: String,
    pub host_pos: i64,
    pub direction: String,
}

/// 位点聚类：同 (target contig, host contig) 内按 2bp 锚点聚类（支持度 = 唯一 qname
/// ≥ min_support），再对出站位点做 5bp 锚点去重（位置均值合并，支持度相加）。
pub fn cluster_sites(events: &[SiteEvent], min_support: usize) -> Vec<Site> {
    let mut sorted: Vec<&SiteEvent> = events.iter().collect();
    sorted.sort_by(|a, b| {
        (&a.contig, &a.host_contig, a.pos).cmp(&(&b.contig, &b.host_contig, b.pos))
    });
    let mut sites = Vec::new();
    let mut group: Vec<&SiteEvent> = Vec::new();
    for ev in sorted {
        if let Some(last) = group.last() {
            if last.contig != ev.contig || last.host_contig != ev.host_contig {
                push_clustered(&mut sites, &group, min_support);
                group.clear();
            }
        }
        group.push(ev);
    }
    push_clustered(&mut sites, &group, min_support);
    dedup_sites(&mut sites);
    sites
}

/// 2bp 锚点聚类（viroflash SummarizeIntegrationSite:1111）：同 (contig × host) 组内
/// 窗口起点锚定，positions[j] − positions[i] ≤ 2 归入同簇；支持度 = 簇内唯一 qname 数。
fn push_clustered(sites: &mut Vec<Site>, group: &[&SiteEvent], min_support: usize) {
    let mut i = 0;
    while i < group.len() {
        let anchor = group[i];
        let mut j = i;
        let mut qnames: HashSet<&str> = HashSet::new();
        while j < group.len() && group[j].pos - anchor.pos <= CLUSTER_TOLERANCE {
            qnames.insert(group[j].qname.as_str());
            j += 1;
        }
        let cluster = &group[i..j];
        i = j;
        if qnames.len() < min_support {
            continue;
        }
        let pos = mean_round(&cluster.iter().map(|e| e.pos).collect::<Vec<_>>());
        let host_pos = mean_round(&cluster.iter().map(|e| e.host_pos).collect::<Vec<_>>());
        sites.push(Site {
            contig: anchor.contig.clone(),
            pos,
            support: qnames.len(),
            host_contig: anchor.host_contig.clone(),
            host_pos,
            direction: anchor.direction.clone(),
        });
    }
}

/// 5bp 锚点去重（viroflash DeduplicateBreakpointPositions:1314）：位点位置差 ≤5 合并为
/// 均值；合并组的支持度相加（各 2bp 簇的 qname 已去重且位置不重叠）。
fn dedup_sites(sites: &mut Vec<Site>) {
    sites.sort_by(|a, b| {
        (&a.contig, &a.host_contig, a.pos).cmp(&(&b.contig, &b.host_contig, b.pos))
    });
    let mut out: Vec<Site> = Vec::new();
    let mut i = 0;
    while i < sites.len() {
        let anchor_pos = sites[i].pos;
        let anchor_host = sites[i].host_contig.clone();
        let anchor_contig = sites[i].contig.clone();
        let anchor_dir = sites[i].direction.clone();
        let mut j = i;
        let mut pos_sum = 0i64;
        let mut host_sum = 0i64;
        let mut support = 0usize;
        while j < sites.len()
            && sites[j].contig == anchor_contig
            && sites[j].host_contig == anchor_host
            && sites[j].pos - anchor_pos <= DEDUP_TOLERANCE
        {
            pos_sum += sites[j].pos;
            host_sum += sites[j].host_pos;
            support += sites[j].support;
            j += 1;
        }
        let n = (j - i) as i64;
        out.push(Site {
            contig: anchor_contig,
            pos: (pos_sum as f64 / n as f64).round() as i64,
            support,
            host_contig: anchor_host,
            host_pos: (host_sum as f64 / n as f64).round() as i64,
            direction: anchor_dir,
        });
        i = j;
    }
    *sites = out;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ev(contig: &str, pos: i64, host_pos: i64, qname: &str) -> SiteEvent {
        SiteEvent {
            contig: contig.into(),
            pos,
            host_contig: "host_0".into(),
            host_pos,
            direction: "host+_target+".into(),
            qname: qname.into(),
        }
    }

    #[test]
    fn cluster_sites_clusters_2bp_and_counts_qnames() {
        // 100/101 同 2bp 簇 → 位点 (100+101)/2 → 101；105 单独成簇支持 1 被丢弃
        let events = vec![
            ev("target_0", 100, 1000, "r1"),
            ev("target_0", 101, 1001, "r2"),
            ev("target_0", 105, 1005, "r3"),
        ];
        let sites = cluster_sites(&events, 2);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].contig, "target_0");
        assert_eq!(sites[0].pos, 101);
        assert_eq!(sites[0].support, 2);
        assert_eq!(sites[0].host_pos, 1001);
    }

    #[test]
    fn cluster_sites_dedups_nearby_sites_by_5bp() {
        // 两个 2bp 簇（101 与 106），5bp 去重合并为 (101+106)/2 → 104，支持度相加
        let events = vec![
            ev("target_0", 100, 1000, "r1"),
            ev("target_0", 101, 1001, "r2"),
            ev("target_0", 105, 1005, "r3"),
            ev("target_0", 106, 1006, "r4"),
        ];
        let sites = cluster_sites(&events, 2);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].pos, 104);
        assert_eq!(sites[0].support, 4);
    }

    #[test]
    fn cluster_sites_splits_beyond_5bp() {
        let events = vec![
            ev("target_0", 100, 1000, "r1"),
            ev("target_0", 101, 1001, "r2"),
            ev("target_0", 107, 1007, "r3"),
            ev("target_0", 108, 1008, "r4"),
        ];
        let sites = cluster_sites(&events, 2);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].pos, 101);
        assert_eq!(sites[0].support, 2);
        assert_eq!(sites[1].pos, 108);
        assert_eq!(sites[1].support, 2);
    }

    #[test]
    fn cluster_sites_dedups_qnames_for_support() {
        // 同一 fragment（r1/r2 各一条事件，同 qname）只计 1 个支持
        let events = vec![
            ev("target_0", 100, 1000, "same_read"),
            ev("target_0", 100, 1000, "same_read"),
        ];
        assert!(cluster_sites(&events, 2).is_empty());
        let mut with_other = events.clone();
        with_other.push(ev("target_0", 100, 1000, "other_read"));
        let sites = cluster_sites(&with_other, 2);
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].support, 2);
    }

    #[test]
    fn cluster_sites_groups_by_host_contig() {
        let events = vec![
            ev("target_0", 100, 1000, "r1"),
            ev("target_0", 101, 1001, "r2"),
            {
                let mut e = ev("target_0", 100, 1000, "r3");
                e.host_contig = "host_1".into();
                e
            },
            {
                let mut e = ev("target_0", 101, 1001, "r4");
                e.host_contig = "host_1".into();
                e
            },
        ];
        let sites = cluster_sites(&events, 2);
        assert_eq!(sites.len(), 2);
        assert_eq!(sites[0].host_contig, "host_0");
        assert_eq!(sites[1].host_contig, "host_1");
    }

    #[test]
    fn cluster_sites_respects_min_support() {
        let events = vec![ev("target_0", 100, 1000, "r1")];
        assert!(cluster_sites(&events, 2).is_empty());
    }
}
