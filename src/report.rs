//! 报告输出：候选级 JSON + TSV（手工序列化，避免引入 serde 依赖）。
//! 输出保持候选级结论，不生成样本级结论。

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cluster::Site;

/// 基于 q 值的置信度档位。
pub fn confidence_for_q(q: f64) -> &'static str {
    if q < 0.01 {
        "HIGH"
    } else if q < 0.05 {
        "MEDIUM"
    } else if q < 0.2 {
        "LOW"
    } else {
        "NOT_SIGNIFICANT"
    }
}

#[derive(Debug, Clone)]
pub struct Candidate {
    pub contig: String,
    pub contig_len: u64,
    pub covered_bases: u64,
    pub covered_frac: f64,
    pub reads: u64,
    /// split 断点事件数（同一 read 双侧 softclip 可产生 2 条事件）。
    pub split_events: u64,
    pub discordant: u64,
    pub plus_strand: u64,
    pub minus_strand: u64,
    pub sites: Vec<Site>,
    pub p_value: f64,
    pub q_value: f64,
    /// 本层经验 p 的分辨率下限 = 1/(层内 decoy 数 + 1)，供审计（层内 decoy 少时
    /// HIGH/MEDIUM 档在数学上不可达）。
    pub p_resolution_floor: f64,
    pub poisson_p: Option<f64>,
    pub nb_p: Option<f64>,
    pub stratum: String,
    pub stratum_decoy_count: usize,
    /// 非 split 接受 reads（reads − 唯一 split qname 数，避免 split read 双计）。
    pub n_plain: u64,
    /// Poisson 期望 λ̂_ℓ*·E_t（长度校正的诱饵 reads 期望）。
    pub expected_hits: f64,
    /// 层背景率 λ̂_ℓ*（10M reads 参考深度的每碱基率，EB 收缩后）。
    pub lambda_bg_layer: f64,
    /// n_plain / expected_hits；期望 = 0 时 None（TSV 输出 "-"）。
    pub depth_fold: Option<f64>,
    /// reads per million input fragments（Mourik 2024 RPM 口径；组合门槛主维）。
    pub depth_rpm: f64,
    /// p ≤ 1/6 地板 flag；仅披露，不进 q。
    pub p_floor_flag: bool,
    /// 全样本池化 π̂0（CDH2018 离散感知估计）。
    pub pi0: f64,
    /// 决策：q≥0.2 → NOT_SIGNIFICANT；q<0.2 但不满足深度×覆盖组合门槛 → BELOW_THRESHOLD；
    /// 过门槛：split 确认 → PASS，无 split → WEAK。
    pub decision: &'static str,
    /// 合并后 distinct 窗口数，定义为合并区间数。
    pub distinct_windows: u64,
}

impl Candidate {
    pub fn confidence(&self) -> &'static str {
        confidence_for_q(self.q_value)
    }
}

pub fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 8);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        Some(x) => format!("{x:.6}"),
        None => "null".to_string(),
    }
}

/// 写 `{out}.json` 与 `{out}.tsv`，返回两个路径。
// 参数直接对应稳定的运行元数据与候选输出，聚合封装会增加无必要层级。
#[allow(clippy::too_many_arguments)]
pub fn write_report(
    out_prefix: &Path,
    sample: &str,
    threads: usize,
    k: usize,
    input_pairs: u64,
    prescreen_pairs: u64,
    map_errors: u64,
    candidates: &[Candidate],
) -> Result<(PathBuf, PathBuf), String> {
    // 前缀拼接（与 work_dir 的 `{out}.work` 约定一致；with_extension 会剥掉 .out 等后缀）
    let json_path = PathBuf::from(format!("{}.json", out_prefix.display()));
    let tsv_path = PathBuf::from(format!("{}.tsv", out_prefix.display()));

    {
        let mut body = String::new();
        body.push_str("{\n");
        body.push_str("  \"schema\": \"viroflash.result.v0\",\n");
        body.push_str(&format!(
            "  \"run\": {{\"sample\": \"{}\", \"threads\": {threads}, \"k\": {k}, \"input_pairs\": {input_pairs}, \"prescreen_pairs\": {prescreen_pairs}, \"map_errors\": {map_errors}}},\n",
            json_escape(sample)
        ));
        // 阈值 source 分类：literature = 文献公式；community-convention = 社区惯例；
        // uncalibrated = 尚无正式校准依据。
        body.push_str(&format!(
            "  \"thresholds\": {{\n    \"k\": {{\"value\": {k}, \"source\": \"community-convention\"}},\n    \"min_mapq\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"max_nm\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"min_as_diff\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"best_n\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"split_softclip\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"split_mapq\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"site_cluster_bp\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"site_dedup_bp\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"min_site_support\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"decoy_p_plus_one\": {{\"source\": \"literature\"}},\n    \"bh_q\": {{\"source\": \"literature\"}},\n    \"nb_poisson_tail\": {{\"source\": \"literature\"}},\n    \"depth_rpm_min\": {{\"value\": {:.0}, \"source\": \"literature\"}},\n    \"coverage_min\": {{\"value\": {:.2}, \"source\": \"literature\"}},\n    \"confidence_bands\": {{\"source\": \"uncalibrated\"}}\n  }},\n",
            crate::align::MIN_MAPQ,
            crate::align::MAX_NM,
            crate::align::MIN_AS_DIFF,
            crate::align::BEST_N,
            crate::align::SPLIT_SOFTCLIP,
            crate::align::SPLIT_MAPQ,
            crate::cluster::CLUSTER_TOLERANCE,
            crate::cluster::DEDUP_TOLERANCE,
            crate::cluster::MIN_SITE_SUPPORT,
            crate::DEPTH_RPM_MIN,
            crate::COVERAGE_MIN,
        ));
        body.push_str("  \"candidates\": [\n");
        for (i, c) in candidates.iter().enumerate() {
            body.push_str("    {\n");
            body.push_str(&format!(
                "      \"contig\": \"{}\",\n",
                json_escape(&c.contig)
            ));
            body.push_str(&format!("      \"confidence\": \"{}\",\n", c.confidence()));
            body.push_str(&format!(
                "      \"q_value\": {:.6},\n      \"p_value\": {:.6},\n",
                c.q_value, c.p_value
            ));
            body.push_str(&format!(
                "      \"p_resolution_floor\": {:.6},\n",
                c.p_resolution_floor
            ));
            body.push_str(&format!(
                "      \"poisson_p\": {},\n      \"nb_p\": {},\n",
                fmt_opt(c.poisson_p),
                fmt_opt(c.nb_p)
            ));
            body.push_str(&format!(
                "      \"stratum\": \"{}\",\n      \"stratum_decoy_count\": {},\n",
                json_escape(&c.stratum),
                c.stratum_decoy_count
            ));
            body.push_str(&format!(
                "      \"pi0\": {:.6},\n      \"decision\": \"{}\",\n",
                c.pi0, c.decision
            ));
            body.push_str(&format!(
                "      \"n_plain\": {},\n      \"expected_hits\": {:.3},\n      \"lambda_bg_layer\": {:.6e},\n      \"depth_fold\": {},\n      \"p_floor_flag\": {},\n      \"distinct_windows\": {},\n",
                c.n_plain,
                c.expected_hits,
                c.lambda_bg_layer,
                fmt_opt(c.depth_fold),
                c.p_floor_flag,
                c.distinct_windows
            ));
            body.push_str(&format!(
                "      \"evidence\": {{\"reads\": {}, \"covered_bases\": {}, \"contig_len\": {}, \"covered_frac\": {:.6}, \"depth_rpm\": {:.3}, \"split_events\": {}, \"discordant\": {}, \"plus_strand\": {}, \"minus_strand\": {}}},\n",
                c.reads, c.covered_bases, c.contig_len, c.covered_frac, c.depth_rpm, c.split_events, c.discordant, c.plus_strand, c.minus_strand
            ));
            body.push_str("      \"sites\": [");
            for (j, s) in c.sites.iter().enumerate() {
                if j > 0 {
                    body.push_str(", ");
                }
                body.push_str(&format!(
                    "{{\"pos\": {}, \"support\": {}, \"host\": \"{}\", \"host_pos\": {}, \"direction\": \"{}\"}}",
                    s.pos,
                    s.support,
                    json_escape(&s.host_contig),
                    s.host_pos,
                    json_escape(&s.direction)
                ));
            }
            body.push_str("]\n");
            body.push_str(if i + 1 < candidates.len() {
                "    },\n"
            } else {
                "    }\n"
            });
        }
        body.push_str("  ]\n}\n");

        let file = File::create(&json_path)
            .map_err(|e| format!("无法创建 {}: {e}", json_path.display()))?;
        let mut w = BufWriter::new(file);
        w.write_all(body.as_bytes())
            .map_err(|e| format!("写报告失败: {e}"))?;
        w.flush().map_err(|e| format!("写报告失败: {e}"))?;
    }

    {
        // TSV：固定 22 列的行级契约。
        // q<0.2 候选逐行明细；其余候选合并为每样本一行 NOT_SIGNIFICANT 汇总（供审计）。
        let file =
            File::create(&tsv_path).map_err(|e| format!("无法创建 {}: {e}", tsv_path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(
            w,
            "sample_id\tvirus_id\ttaxid\tresolution_level\treads_plain\tsplit_reads\tdiscordant_pairs\tdistinct_windows\taligned_bases\texpected_hits\tlambda_bg_layer\tdepth_fold\tp_value\tp_floor_flag\tpi0\tq_value\tdecision\tevidence_strength\tcoverage_breadth\tbccp\tdecaf_grade\tnotes"
        )
        .map_err(|e| format!("写报告失败: {e}"))?;
        let fmt_depth = |v: Option<f64>| match v {
            Some(x) => format!("{x:.3}"),
            None => "-".to_string(),
        };
        let significant: Vec<&Candidate> = candidates.iter().filter(|c| c.q_value < 0.2).collect();
        let merged: Vec<&Candidate> = candidates.iter().filter(|c| c.q_value >= 0.2).collect();
        for c in &significant {
            // 明细行 notes 包含 stratum、decoy 数和 depth_rpm，供审计。
            let notes = format!(
                "stratum={};decoy_n={};depth_rpm={:.3}",
                c.stratum, c.stratum_decoy_count, c.depth_rpm
            );
            writeln!(
                w,
                "{}\t{}\t\tcontig\t{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.6e}\t{}\t{:.6}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{:.6}\t\t\t{}",
                sample,
                c.contig,
                c.n_plain,
                c.split_events,
                c.discordant,
                c.distinct_windows,
                c.covered_bases,
                c.expected_hits,
                c.lambda_bg_layer,
                fmt_depth(c.depth_fold),
                c.p_value,
                if c.p_floor_flag { "1" } else { "0" },
                c.pi0,
                c.q_value,
                c.decision,
                c.confidence(),
                c.covered_frac,
                notes
            )
            .map_err(|e| format!("写报告失败: {e}"))?;
        }
        // NOT_SIGNIFICANT 合并行：无任何候选时也输出一行（targets_merged=0），保持每样本至少一行。
        if significant.is_empty() || !merged.is_empty() {
            let max_q = merged
                .iter()
                .map(|c| c.q_value)
                .fold(f64::NEG_INFINITY, f64::max);
            let max_q_s = if merged.is_empty() {
                "-".to_string()
            } else {
                format!("{max_q:.6}")
            };
            writeln!(
                w,
                "{}\t-\t\t\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNOT_SIGNIFICANT\tNOT_SIGNIFICANT\t-\t\t\ttargets_merged={}; max_q={}",
                sample,
                merged.len(),
                max_q_s
            )
            .map_err(|e| format!("写报告失败: {e}"))?;
        }
        w.flush().map_err(|e| format!("写报告失败: {e}"))?;
    }

    Ok((json_path, tsv_path))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn confidence_bands() {
        assert_eq!(confidence_for_q(0.001), "HIGH");
        assert_eq!(confidence_for_q(0.02), "MEDIUM");
        assert_eq!(confidence_for_q(0.1), "LOW");
        assert_eq!(confidence_for_q(0.5), "NOT_SIGNIFICANT");
    }

    #[test]
    fn escape_quotes_and_controls() {
        assert_eq!(json_escape("a\"b"), "a\\\"b");
        assert_eq!(json_escape("a\nb"), "a\\nb");
        assert_eq!(json_escape("正常中文"), "正常中文");
    }
}
