//! 报告输出：候选级 JSON + TSV（手工序列化，避免引入 serde 依赖）。
//! 输出保持候选级结论，不生成样本级结论。

use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::cluster::Site;

/// 当前候选结果契约。v1 明确区分模型 adjusted-p、候选报告门和 integration 证据；
/// 这些语义与旧 v0 不兼容，不能静默沿用同一版本号。
pub const RESULT_SCHEMA: &str = "viroflash.result.v1";

/// 兼容字段 `confidence` 的输出。当前 synthetic-decoy null 尚未独立校准，因此
/// adjusted-p 不能转换成 HIGH/MEDIUM/LOW 置信度。
pub fn confidence_for_q(q: f64) -> &'static str {
    if q < crate::MODEL_ADJUSTED_P_MAX {
        "UNVALIDATED"
    } else {
        "NOT_SIGNIFICANT"
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub contig: String,
    /// 检测组代表在用户 target FASTA 中的原始名称。
    pub representative: String,
    /// 复合 OR 假设覆盖的原始 target 名称；不能据此归因到其中任一成员。
    pub hypothesis_members: Vec<String>,
    /// discovery 证据的最小 hitting-set 解释，仅用于压缩展示。
    pub hypothesis_explanation: Vec<String>,
    /// 构造该假设的 discovery read-end 数；不参与 validation 统计量。
    pub discovery_reads: u64,
    pub contig_len: u64,
    /// discovery OR 假设包含的内部索引 target contig 数。
    pub index_member_count: usize,
    /// target 直接 read 计数对应的内部索引参考总长度。
    pub target_exposure_bases: u64,
    pub covered_bases: u64,
    pub covered_frac: f64,
    pub reads: u64,
    /// split 断点事件数（同一 read 双侧 softclip 可产生 2 条事件）。
    pub split_events: u64,
    pub discordant: u64,
    pub plus_strand: u64,
    pub minus_strand: u64,
    pub sites: Vec<Site>,
    /// 精确条件 Poisson 率检验的普通 f64 p；若下溢为 0，`ln_p_value` 仍保留信息。
    pub p_value: f64,
    /// 固定 discovery 检验族内的 BH adjusted-p；不是已验证的 classical-FDR q。
    pub q_value: f64,
    pub ln_p_value: f64,
    pub ln_q_value: f64,
    pub p_underflow: bool,
    pub q_underflow: bool,
    /// 零 background read、一个 target read 时精确条件检验的 exposure 参考概率。
    /// 字段名为 v0 兼容保留；它不是多 read 检验的数学硬地板。
    pub p_resolution_floor: f64,
    /// 旧 plug-in Poisson 上尾，仅作诊断，不参与主 adjusted-p 或判定。
    pub poisson_p: Option<f64>,
    /// v0 兼容字段；当前主路径不拟合未经验证的 NB dispersion。
    pub nb_p: Option<f64>,
    pub stratum: String,
    pub stratum_decoy_count: usize,
    /// v0 兼容字段；现与 `reads` 相同，所有直接 validation read-end 只计一次。
    /// split 是独立 integration evidence，不再从通用检测计数中扣除。
    pub n_plain: u64,
    /// 同层诱饵直接 reads 按 target/background 参考长度比换算的 target 期望。
    pub expected_hits: f64,
    /// 同层诱饵直接 reads 率，归一化为每 10M validation read-ends 的每碱基率；
    /// 不含未经验证的 EB 收缩。
    pub lambda_bg_layer: f64,
    /// n_plain / expected_hits；期望 = 0 时 None（TSV 输出 "-"）。
    pub depth_fold: Option<f64>,
    /// validation read-ends / million selected validation read-ends；只作丰度描述，
    /// 不作为未经校准的硬判定门槛。
    pub depth_rpm: f64,
    /// p 是否不高于单 target event 的 exposure 参考概率；仅披露，不进 adjusted-p。
    pub p_floor_flag: bool,
    /// v0 兼容字段；固定族 BH 使用 π0=1。
    pub pi0: f64,
    /// 候选报告门：模型 adjusted-p、breadth 和分布窗口；不是样本级/临床结论。
    pub decision: &'static str,
    pub decision_reasons: Vec<&'static str>,
    /// 与通用检测正交：NONE / UNCLUSTERED_SPLIT / SUPPORTED_SITE。
    pub integration_evidence: &'static str,
    /// validation 对齐中点占据的代表序列固定位置分区数。
    pub distinct_windows: u64,
    /// discovery 固定的 BH 检验族大小，包含 validation 为零的假设。
    pub test_family_size: usize,
    pub background_status: &'static str,
    pub background_scope: &'static str,
    pub background_reads: u64,
    pub background_reference_bases: u64,
    pub background_cross_stratum_reads: u64,
}

/// 抽样、折分和歧义审计元数据。所有计数均为实际观测值；不把 Horvitz–Thompson
/// 估计伪装成全量精确计数。
#[derive(Debug, Clone, PartialEq)]
pub struct SamplingReport {
    pub method: &'static str,
    pub capacity_pairs: usize,
    pub selected_pairs: u64,
    pub selected_passed_pairs: u64,
    pub discovery_pairs: u64,
    pub validation_pairs: u64,
    pub validation_read_sides: u64,
    pub inclusion_probability: f64,
    pub duplicate_selected_qnames: u64,
    pub audit_reads: u64,
    pub audit_no_evidence: u64,
    pub audit_overflows: u64,
    pub audit_map_error_reads: u64,
    pub target_validation_unassigned: u64,
    pub decoy_validation_unassigned: u64,
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
        Some(x) => x.to_string(),
        None => "null".to_string(),
    }
}

fn fmt_sci(v: f64) -> String {
    format!("{v:.17e}")
}

fn probability_status(underflow: bool) -> &'static str {
    if underflow {
        "NUMERIC_UNDERFLOW_LOG_AVAILABLE"
    } else {
        "FINITE"
    }
}

/// 写 `{out}.json` 与 `{out}.tsv`，返回两个路径。
/// `index_*`：本次运行所用索引的溯源摘要（source=loaded|built，manifest 的
/// BLAKE3），写入 JSON 的顶层 `index` 块；
/// TSV 列契约不变。
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
    index_source: &str,
    index_format_version: u32,
    manifest_blake3: &str,
    sampling: &SamplingReport,
    test_family_size: usize,
    candidates: &[Candidate],
) -> Result<(PathBuf, PathBuf), String> {
    if candidates.len() > test_family_size
        || candidates
            .iter()
            .any(|candidate| candidate.test_family_size != test_family_size)
    {
        return Err("报告候选与固定检验族大小不一致".to_string());
    }
    for candidate in candidates {
        let required = [
            ("p_value", candidate.p_value),
            ("q_value", candidate.q_value),
            ("ln_p_value", candidate.ln_p_value),
            ("ln_q_value", candidate.ln_q_value),
            ("covered_frac", candidate.covered_frac),
            ("p_resolution_floor", candidate.p_resolution_floor),
            ("expected_hits", candidate.expected_hits),
            ("lambda_bg_layer", candidate.lambda_bg_layer),
            ("depth_rpm", candidate.depth_rpm),
            ("pi0", candidate.pi0),
        ];
        for (name, value) in required {
            if !value.is_finite() {
                return Err(format!("候选 {} 的 {name} 非有限", candidate.contig));
            }
        }
        for (name, value) in [
            ("poisson_p", candidate.poisson_p),
            ("nb_p", candidate.nb_p),
            ("depth_fold", candidate.depth_fold),
        ] {
            if value.is_some_and(|value| !value.is_finite()) {
                return Err(format!("候选 {} 的 {name} 非有限", candidate.contig));
            }
        }
    }
    // 前缀拼接（与 work_dir 的 `{out}.work` 约定一致；with_extension 会剥掉 .out 等后缀）
    let json_path = PathBuf::from(format!("{}.json", out_prefix.display()));
    let tsv_path = PathBuf::from(format!("{}.tsv", out_prefix.display()));

    {
        let mut body = String::new();
        body.push_str("{\n");
        body.push_str(&format!("  \"schema\": \"{RESULT_SCHEMA}\",\n"));
        body.push_str(&format!(
            "  \"run\": {{\"sample\": \"{}\", \"threads\": {threads}, \"k\": {k}, \"input_pairs\": {input_pairs}, \"prescreen_pairs\": {prescreen_pairs}, \"map_errors\": {map_errors}}},\n",
            json_escape(sample)
        ));
        body.push_str(&format!(
            "  \"sampling\": {{\"method\": \"{}\", \"capacity_pairs\": {}, \"selected_pairs\": {}, \"selected_passed_pairs\": {}, \"discovery_pairs\": {}, \"validation_pairs\": {}, \"validation_read_sides\": {}, \"inclusion_probability\": {:.9}, \"duplicate_selected_qnames\": {}, \"audit_reads\": {}, \"audit_no_evidence\": {}, \"audit_overflows\": {}, \"audit_map_error_reads\": {}, \"target_validation_unassigned\": {}, \"decoy_validation_unassigned\": {}, \"prescreen_scope\": \"selected_sample\"}},\n",
            json_escape(sampling.method),
            sampling.capacity_pairs,
            sampling.selected_pairs,
            sampling.selected_passed_pairs,
            sampling.discovery_pairs,
            sampling.validation_pairs,
            sampling.validation_read_sides,
            sampling.inclusion_probability,
            sampling.duplicate_selected_qnames,
            sampling.audit_reads,
            sampling.audit_no_evidence,
            sampling.audit_overflows,
            sampling.audit_map_error_reads,
            sampling.target_validation_unassigned,
            sampling.decoy_validation_unassigned,
        ));
        body.push_str("  \"report_semantics\": {\"scope\": \"candidate\", \"sample_conclusion\": \"not_computed\", \"empty_candidates\": \"no_candidates_reported_not_equivalent_to_not_detected\", \"decision_scope\": \"research_candidate_reporting_gate_not_clinical_diagnosis\", \"confidence_calibration\": \"not_validated\", \"q_value_semantics\": \"model_adjusted_p_compatibility_alias_not_validated_fdr\", \"p_resolution_floor_semantics\": \"one_target_event_zero_background_exposure_reference_legacy_field_name\", \"legacy_probability_fields\": \"poisson_p_diagnostic_only_nb_p_disabled\"},\n");
        body.push_str(&format!(
            "  \"candidate_testing\": {{\"family\": \"all_discovery_hypotheses_including_zero_validation_evidence\", \"test_family_size\": {test_family_size}, \"reported_candidates\": {}, \"unreported_zero_evidence_hypotheses\": {}}},\n",
            candidates.len(),
            test_family_size - candidates.len(),
        ));
        let mut qc_issues = Vec::new();
        if map_errors > 0 || sampling.audit_map_error_reads > 0 {
            qc_issues.push("mapping_errors_observed");
        }
        if sampling.audit_overflows > 0 {
            qc_issues.push("ambiguity_audit_overflow_observed");
        }
        if sampling.target_validation_unassigned > 0 {
            qc_issues.push("target_validation_unassigned_observed");
        }
        if sampling.decoy_validation_unassigned > 0 {
            qc_issues.push("decoy_validation_unassigned_observed");
        }
        body.push_str("  \"quality_control\": {\"status\": \"NOT_EVALUATED\", \"method\": \"viroflash_internal_observations\", \"issues\": [");
        for (index, issue) in qc_issues.iter().enumerate() {
            if index > 0 {
                body.push_str(", ");
            }
            body.push_str(&format!("\"{}\"", json_escape(issue)));
        }
        body.push_str(&format!(
            "], \"details\": {{\"map_errors\": {map_errors}, \"audit_overflows\": {}, \"target_validation_unassigned\": {}, \"decoy_validation_unassigned\": {}}}, \"interpretation\": \"observations_only_no_validated_qc_pass_fail_threshold\"}},\n",
            sampling.audit_overflows,
            sampling.target_validation_unassigned,
            sampling.decoy_validation_unassigned,
        ));
        // 阈值 source 分类：literature = 文献公式；community-convention = 社区惯例；
        // literature-adapted = 文献原则的公开算法化实现；uncalibrated = 尚无正式校准依据。
        body.push_str(&format!(
            "  \"thresholds\": {{\n    \"k\": {{\"value\": {k}, \"source\": \"community-convention\"}},\n    \"min_mapq\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"max_nm\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"min_as_diff\": {{\"value\": {}, \"source\": \"uncalibrated\"}},\n    \"best_n\": {{\"value\": {}, \"source\": \"uncalibrated\", \"scope\": \"ordinary_mapper_only\"}},\n    \"ambiguity_audit\": {{\"mode\": \"all_chains\", \"max_hits_per_read\": {}}},\n    \"split_softclip\": {{\"value\": {}, \"source\": \"community-convention\", \"decision_role\": \"integration_evidence_only\"}},\n    \"split_mapq\": {{\"value\": {}, \"source\": \"community-convention\", \"decision_role\": \"integration_evidence_only\"}},\n    \"site_cluster_bp\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"site_dedup_bp\": {{\"value\": {}, \"source\": \"community-convention\"}},\n    \"min_site_support\": {{\"value\": {}, \"source\": \"community-convention\", \"calibration_status\": \"unvalidated\"}},\n    \"synthetic_decoy_calibration\": {{\"source\": \"uncalibrated\", \"classical_fdr_guarantee\": false, \"role\": \"target_derived_null_stress_reference\"}},\n    \"model_p\": {{\"method\": \"exact_conditional_two_poisson_rates\", \"reference_id\": \"R_stats_poisson.test\", \"target_exposure\": \"sum_internal_index_contig_bases\", \"background_pool\": \"size_gc_for_single_index_member_global_for_composite\", \"null_calibration_status\": \"unvalidated_synthetic_decoy_exchangeability\"}},\n    \"multiple_testing\": {{\"method\": \"benjamini_hochberg\", \"family\": \"all_discovery_hypotheses_including_zero_validation_evidence\", \"classical_fdr_guarantee\": false}},\n    \"model_adjusted_p_max\": {{\"value\": {:.2}, \"source\": \"uncalibrated\", \"prespecified\": true}},\n    \"coverage_min\": {{\"value\": {:.2}, \"source\": \"literature-adapted\", \"reference_id\": \"doi:10.1128/jcm.00345-24\", \"applicability\": \"not_validated_for_all_viroflash_inputs\"}},\n    \"min_distributed_windows\": {{\"value\": {}, \"bins\": {}, \"source\": \"uncalibrated\", \"rationale\": \"operationalizes_distributed_nonoverlapping_evidence\"}},\n    \"confidence_bands\": {{\"status\": \"withdrawn_until_calibrated\"}}\n  }},\n",
            crate::align::MIN_MAPQ,
            crate::align::MAX_NM,
            crate::align::MIN_AS_DIFF,
            crate::align::BEST_N,
            crate::align::MAX_AUDIT_HITS,
            crate::align::SPLIT_SOFTCLIP,
            crate::align::SPLIT_MAPQ,
            crate::cluster::CLUSTER_TOLERANCE,
            crate::cluster::DEDUP_TOLERANCE,
            crate::cluster::MIN_SITE_SUPPORT,
            crate::MODEL_ADJUSTED_P_MAX,
            crate::COVERAGE_MIN,
            crate::MIN_DISTRIBUTED_WINDOWS,
            crate::DISTRIBUTED_WINDOW_BINS,
        ));
        // 索引溯源：本次运行所用索引的来源、格式版本与 manifest 校验和
        //（可审计的参数来源约定；旧字段不受影响）。
        body.push_str(&format!(
            "  \"index\": {{\"source\": \"{}\", \"format_version\": {}, \"manifest_blake3\": \"{}\"}},\n",
            index_source, index_format_version, manifest_blake3
        ));
        body.push_str("  \"candidates\": [\n");
        for (i, c) in candidates.iter().enumerate() {
            let stratum_scope = if c.background_scope == "global_composite_hypothesis" {
                "representative_only_not_test_background"
            } else {
                "test_background_stratum"
            };
            body.push_str("    {\n");
            body.push_str(&format!(
                "      \"contig\": \"{}\",\n",
                json_escape(&c.contig)
            ));
            body.push_str(&format!(
                "      \"representative\": \"{}\",\n",
                json_escape(&c.representative)
            ));
            body.push_str("      \"hypothesis\": {\"semantics\": \"at_least_one_member_present\", \"resolution\": \"reference_group\", \"member_attribution\": \"not_resolved\", \"members\": [");
            for (member_index, member) in c.hypothesis_members.iter().enumerate() {
                if member_index > 0 {
                    body.push_str(", ");
                }
                body.push_str(&format!("\"{}\"", json_escape(member)));
            }
            body.push_str("], \"explanation\": [");
            for (member_index, member) in c.hypothesis_explanation.iter().enumerate() {
                if member_index > 0 {
                    body.push_str(", ");
                }
                body.push_str(&format!("\"{}\"", json_escape(member)));
            }
            body.push_str(&format!(
                "], \"discovery_reads\": {}, \"index_member_count\": {}}},\n",
                c.discovery_reads, c.index_member_count
            ));
            body.push_str(&format!(
                "      \"confidence\": \"{}\",\n      \"confidence_basis\": \"not_calibrated_model_adjusted_p\",\n",
                c.confidence()
            ));
            body.push_str(&format!(
                "      \"q_value\": {},\n      \"p_value\": {},\n",
                c.q_value, c.p_value
            ));
            body.push_str(&format!(
                "      \"p_resolution_floor\": {},\n",
                c.p_resolution_floor
            ));
            body.push_str(&format!(
                "      \"poisson_p\": {},\n      \"nb_p\": {},\n",
                fmt_opt(c.poisson_p),
                fmt_opt(c.nb_p)
            ));
            body.push_str(&format!(
                "      \"stratum\": \"{}\",\n      \"stratum_scope\": \"{}\",\n      \"stratum_decoy_count\": {},\n",
                json_escape(&c.stratum),
                stratum_scope,
                c.stratum_decoy_count
            ));
            body.push_str(&format!(
                "      \"statistical_evidence\": {{\"test\": \"exact_conditional_two_poisson_rates\", \"alternative\": \"target_rate_greater\", \"target_exposure_bases\": {}, \"p_value\": {}, \"log10_p_value\": {}, \"p_value_status\": \"{}\", \"adjustment\": \"benjamini_hochberg_fixed_discovery_family\", \"adjusted_p_value\": {}, \"log10_adjusted_p_value\": {}, \"adjusted_p_status\": \"{}\", \"family_size\": {}, \"fdr_control_validated\": false, \"background\": {{\"status\": \"{}\", \"scope\": \"{}\", \"reads\": {}, \"reference_bases\": {}, \"cross_stratum_reads\": {}}}}},\n",
                c.target_exposure_bases,
                c.p_value,
                c.ln_p_value / std::f64::consts::LN_10,
                probability_status(c.p_underflow),
                c.q_value,
                c.ln_q_value / std::f64::consts::LN_10,
                probability_status(c.q_underflow),
                c.test_family_size,
                c.background_status,
                c.background_scope,
                c.background_reads,
                c.background_reference_bases,
                c.background_cross_stratum_reads,
            ));
            body.push_str(&format!(
                "      \"pi0\": {},\n      \"decision\": \"{}\",\n",
                c.pi0, c.decision
            ));
            body.push_str("      \"decision_reasons\": [");
            for (reason_index, reason) in c.decision_reasons.iter().enumerate() {
                if reason_index > 0 {
                    body.push_str(", ");
                }
                body.push_str(&format!("\"{}\"", json_escape(reason)));
            }
            body.push_str("],\n");
            body.push_str(&format!(
                "      \"integration_evidence\": {{\"status\": \"{}\", \"decision_role\": \"separate_from_virus_detection\", \"split_events\": {}, \"supported_sites\": {}}},\n",
                c.integration_evidence,
                c.split_events,
                c.sites.len(),
            ));
            body.push_str(&format!(
                "      \"n_plain\": {},\n      \"expected_hits\": {},\n      \"lambda_bg_layer\": {},\n      \"depth_fold\": {},\n      \"p_floor_flag\": {},\n      \"distinct_windows\": {},\n",
                c.n_plain,
                c.expected_hits,
                c.lambda_bg_layer,
                fmt_opt(c.depth_fold),
                c.p_floor_flag,
                c.distinct_windows
            ));
            body.push_str(&format!(
                "      \"evidence\": {{\"reads\": {}, \"unit\": \"validation_read_end\", \"covered_bases\": {}, \"contig_len\": {}, \"covered_frac\": {}, \"breadth_scope\": \"validation_sample_lower_bound\", \"depth_rpm\": {}, \"rpm_unit\": \"read_end_per_million_validation_read_ends\", \"split_events\": {}, \"discordant\": {}, \"plus_strand\": {}, \"minus_strand\": {}}},\n",
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
        // adjusted-p<探索阈值的候选逐行明细；其余候选合并为审计汇总。空候选汇总
        // 明确声明不等价于样本级 NOT_DETECTED。
        let file =
            File::create(&tsv_path).map_err(|e| format!("无法创建 {}: {e}", tsv_path.display()))?;
        let mut w = BufWriter::new(file);
        writeln!(
            w,
            "sample_id\tvirus_id\ttaxid\tresolution_level\treads_plain\tsplit_reads\tdiscordant_pairs\tdistinct_windows\taligned_bases\texpected_hits\tlambda_bg_layer\tdepth_fold\tp_value\tp_floor_flag\tpi0\tq_value\tdecision\tevidence_strength\tcoverage_breadth\tbccp\tdecaf_grade\tnotes"
        )
        .map_err(|e| format!("写报告失败: {e}"))?;
        let fmt_depth = |v: Option<f64>| match v {
            Some(x) => fmt_sci(x),
            None => "-".to_string(),
        };
        let significant: Vec<&Candidate> = candidates
            .iter()
            .filter(|c| c.q_value < crate::MODEL_ADJUSTED_P_MAX)
            .collect();
        let merged: Vec<&Candidate> = candidates
            .iter()
            .filter(|c| c.q_value >= crate::MODEL_ADJUSTED_P_MAX)
            .collect();
        for c in &significant {
            let reasons = c.decision_reasons.join(",");
            let stratum_scope = if c.background_scope == "global_composite_hypothesis" {
                "representative_only_not_test_background"
            } else {
                "test_background_stratum"
            };
            let notes = format!(
                "row_type=candidate;result_schema={};decision_scope=research_candidate_reporting_gate;sample_conclusion=not_computed;qc_status=NOT_EVALUATED;representative={};hypothesis_semantics=at_least_one_member_present;member_attribution=not_resolved;or_members={};index_members={};target_exposure_bases={};discovery_reads={};stratum={};stratum_scope={};decoy_n={};end_rpm={};breadth_scope=validation_sample_lower_bound;p_method=exact_conditional_two_poisson_rates;p_value_status={};log10_p_value={};adjustment=benjamini_hochberg;family_size={};adjusted_p_status={};log10_adjusted_p={};q_value_semantics=model_adjusted_p_compatibility_alias_not_validated_fdr;fdr_control_validated=false;background_status={};background_scope={};background_reads={};decision_reasons={};integration_evidence={};confidence_basis=not_calibrated_model_adjusted_p",
                RESULT_SCHEMA,
                c.representative,
                c.hypothesis_members.len(),
                c.index_member_count,
                c.target_exposure_bases,
                c.discovery_reads,
                c.stratum,
                stratum_scope,
                c.stratum_decoy_count,
                fmt_sci(c.depth_rpm),
                probability_status(c.p_underflow),
                fmt_sci(c.ln_p_value / std::f64::consts::LN_10),
                c.test_family_size,
                probability_status(c.q_underflow),
                fmt_sci(c.ln_q_value / std::f64::consts::LN_10),
                c.background_status,
                c.background_scope,
                c.background_reads,
                reasons,
                c.integration_evidence,
            );
            writeln!(
                w,
                "{}\t{}\t\treference_group\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t\t\t{}",
                sample,
                c.contig,
                c.n_plain,
                c.split_events,
                c.discordant,
                c.distinct_windows,
                c.covered_bases,
                fmt_sci(c.expected_hits),
                fmt_sci(c.lambda_bg_layer),
                fmt_depth(c.depth_fold),
                fmt_sci(c.p_value),
                if c.p_floor_flag { "1" } else { "0" },
                fmt_sci(c.pi0),
                fmt_sci(c.q_value),
                c.decision,
                c.confidence(),
                fmt_sci(c.covered_frac),
                notes
            )
            .map_err(|e| format!("写报告失败: {e}"))?;
        }
        // 仅合并实际存在但 model adjusted-p 未过门的候选。真正空候选 TSV 只保留
        // 表头，避免伪候选行被下游误解为样本级阴性结论。
        if !candidates.is_empty() && (significant.is_empty() || !merged.is_empty()) {
            let max_q = merged
                .iter()
                .map(|c| c.q_value)
                .fold(f64::NEG_INFINITY, f64::max);
            let max_q_s = if merged.is_empty() {
                "-".to_string()
            } else {
                fmt_sci(max_q)
            };
            writeln!(
                w,
                "{}\t-\t\t\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\tNOT_SIGNIFICANT\tNOT_SIGNIFICANT\t-\t\t\trow_type=candidate_summary;result_schema={};decision_scope=research_candidate_reporting_gate;sample_conclusion=not_computed;qc_status=NOT_EVALUATED;targets_merged={};max_adjusted_p={};q_value_semantics=model_adjusted_p_compatibility_alias_not_validated_fdr;fdr_control_validated=false;no_candidates_reported=false;not_equivalent_to_not_detected=true",
                sample,
                RESULT_SCHEMA,
                merged.len(),
                max_q_s,
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

    fn sampling_report() -> SamplingReport {
        SamplingReport {
            method: "test",
            capacity_pairs: 10,
            selected_pairs: 10,
            selected_passed_pairs: 2,
            discovery_pairs: 5,
            validation_pairs: 5,
            validation_read_sides: 10,
            inclusion_probability: 1.0,
            duplicate_selected_qnames: 0,
            audit_reads: 4,
            audit_no_evidence: 0,
            audit_overflows: 1,
            audit_map_error_reads: 0,
            target_validation_unassigned: 0,
            decoy_validation_unassigned: 1,
        }
    }

    fn candidate() -> Candidate {
        Candidate {
            contig: "target_0".into(),
            representative: "VIRUS_A".into(),
            hypothesis_members: vec!["VIRUS_A".into(), "VIRUS_B".into()],
            hypothesis_explanation: vec!["VIRUS_A".into()],
            discovery_reads: 4,
            contig_len: 1_000,
            index_member_count: 2,
            target_exposure_bases: 2_000,
            covered_bases: 300,
            covered_frac: 0.3,
            reads: 3,
            split_events: 0,
            discordant: 0,
            plus_strand: 2,
            minus_strand: 1,
            sites: Vec::new(),
            p_value: 0.0,
            q_value: 0.0,
            ln_p_value: -1_000.0,
            ln_q_value: -999.0,
            p_underflow: true,
            q_underflow: true,
            p_resolution_floor: 0.01,
            poisson_p: Some(0.0),
            nb_p: None,
            stratum: "sz:<1kb,gc:0.40-0.50".into(),
            stratum_decoy_count: 100,
            n_plain: 3,
            expected_hits: 0.0,
            lambda_bg_layer: 0.0,
            depth_fold: None,
            depth_rpm: 300_000.0,
            p_floor_flag: true,
            pi0: 1.0,
            decision: "PASS",
            decision_reasons: vec!["model_and_distribution_gates_met"],
            integration_evidence: "NONE",
            distinct_windows: 3,
            test_family_size: 2,
            background_status: "SYNTHETIC_DECOY_UNVALIDATED",
            background_scope: "global_composite_hypothesis",
            background_reads: 0,
            background_reference_bases: 100_000,
            background_cross_stratum_reads: 0,
        }
    }

    #[test]
    fn confidence_is_not_inferred_from_uncalibrated_adjusted_p() {
        assert_eq!(confidence_for_q(0.001), "UNVALIDATED");
        assert_eq!(confidence_for_q(0.02), "UNVALIDATED");
        assert_eq!(confidence_for_q(0.1), "UNVALIDATED");
        assert_eq!(confidence_for_q(0.5), "NOT_SIGNIFICANT");
    }

    #[test]
    fn escape_quotes_and_controls() {
        assert_eq!(json_escape("a\"b"), "a\\\"b");
        assert_eq!(json_escape("a\nb"), "a\\nb");
        assert_eq!(json_escape("正常中文"), "正常中文");
    }

    #[test]
    fn report_exposes_candidate_semantics_underflow_and_empty_summary() {
        let unique = format!(
            "viroflash_report_test_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("unnamed")
        );
        let prefix = std::env::temp_dir().join(unique);
        let (json_path, tsv_path) = write_report(
            &prefix,
            "sample",
            1,
            21,
            10,
            2,
            0,
            "loaded",
            2,
            "abc",
            &sampling_report(),
            2,
            &[candidate()],
        )
        .unwrap();
        let json = std::fs::read_to_string(&json_path).unwrap();
        let tsv = std::fs::read_to_string(&tsv_path).unwrap();
        assert!(json.contains("\"schema\": \"viroflash.result.v1\""));
        assert!(json.contains("\"sample_conclusion\": \"not_computed\""));
        assert!(json.contains(
            "\"q_value_semantics\": \"model_adjusted_p_compatibility_alias_not_validated_fdr\""
        ));
        assert!(json.contains(
            "\"legacy_probability_fields\": \"poisson_p_diagnostic_only_nb_p_disabled\""
        ));
        assert!(json.contains("\"confidence\": \"UNVALIDATED\""));
        assert!(json.contains("\"p_value_status\": \"NUMERIC_UNDERFLOW_LOG_AVAILABLE\""));
        assert!(json.contains("\"fdr_control_validated\": false"));
        assert!(json.contains("\"integration_evidence\": {\"status\": \"NONE\""));
        assert!(json.contains("\"member_attribution\": \"not_resolved\""));
        assert!(json.contains("\"index_member_count\": 2"));
        assert!(json.contains("\"target_exposure_bases\": 2000"));
        assert!(json.contains("\"stratum_scope\": \"representative_only_not_test_background\""));
        assert!(json.contains("\"test_family_size\": 2"));
        assert!(json.contains("\"unreported_zero_evidence_hypotheses\": 1"));
        for line in tsv.lines() {
            assert_eq!(line.split('\t').count(), 22, "TSV 列数错误: {line}");
        }
        assert!(tsv.contains("row_type=candidate;"));
        assert!(tsv.contains("result_schema=viroflash.result.v1"));
        assert!(tsv.contains("\treference_group\t"));
        assert!(tsv.contains("decision_scope=research_candidate_reporting_gate"));
        assert!(tsv.contains("qc_status=NOT_EVALUATED"));
        assert!(tsv.contains("stratum_scope=representative_only_not_test_background"));
        assert!(tsv.contains("fdr_control_validated=false"));
        assert!(tsv.contains("adjusted_p_status=NUMERIC_UNDERFLOW_LOG_AVAILABLE"));
        assert!(tsv
            .contains("q_value_semantics=model_adjusted_p_compatibility_alias_not_validated_fdr"));
        let _ = std::fs::remove_file(json_path);
        let _ = std::fs::remove_file(tsv_path);

        let empty_prefix = prefix.with_extension("empty");
        let (json_path, tsv_path) = write_report(
            &empty_prefix,
            "empty",
            1,
            21,
            10,
            0,
            0,
            "loaded",
            2,
            "abc",
            &sampling_report(),
            2,
            &[],
        )
        .unwrap();
        let json = std::fs::read_to_string(&json_path).unwrap();
        let tsv = std::fs::read_to_string(&tsv_path).unwrap();
        assert!(json.contains("\"test_family_size\": 2"));
        assert!(json.contains("\"reported_candidates\": 0"));
        assert!(json.contains("\"unreported_zero_evidence_hypotheses\": 2"));
        assert_eq!(tsv.lines().count(), 1, "空候选 TSV 只能包含表头");
        assert!(!tsv.contains("candidate_summary"));
        let _ = std::fs::remove_file(json_path);
        let _ = std::fs::remove_file(tsv_path);
    }
}
