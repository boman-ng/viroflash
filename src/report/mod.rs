use crate::output::atomic_write;
use std::path::Path;

use crate::evidence::{
    finite_population_interval, AttributionStatus, EvidenceAccumulator, EvidenceStatus,
    INTERVAL_METHOD,
};
use crate::fastq::InputCensus;
use crate::profile::AnalysisProfile;
use crate::sampling::{RunMode, SamplingDesign};

pub const REPORT_SCHEMA: &str = "viroflash.evidence-report.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AnalysisStatus {
    ConformantComplete,
    ConformantWithLimitations,
}
impl AnalysisStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ConformantComplete => "CONFORMANT_COMPLETE",
            Self::ConformantWithLimitations => "CONFORMANT_WITH_LIMITATIONS",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct RunEvidence {
    pub sample_id: String,
    pub analysis_status: AnalysisStatus,
    pub reason_codes: Vec<String>,
    pub input_mode: String,
    pub input_fragments: u64,
    pub selected_fragments: u64,
    pub precision: &'static str,
    pub sampling_population: &'static str,
    pub selection_probability: f64,
    pub minimum_relevant_fraction: f64,
    pub familywise_miss_probability: f64,
    pub interval_level: f64,
    pub target_family_size: usize,
    pub prescreen_passed_fragments: u64,
    pub aligned_fragments: u64,
    pub unassigned_fragments: u64,
    pub profile_digest: String,
    pub index_digest: String,
    pub input_digest: String,
    pub read_ends_per_fragment: u8,
    pub mode: &'static str,
    pub sampling_population_fragments: u64,
    pub sample_capacity: u64,
    pub sampling_expected_support: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct TargetSignal {
    pub target_group_id: String,
    pub representative_id: String,
    pub representative_description: String,
    pub member_ids: Vec<String>,
    pub evidence_status: EvidenceStatus,
    pub attribution_status: AttributionStatus,
    pub supporting_selected_fragments: u64,
    pub selected_fragment_denominator: u64,
    pub attributed_fragment_fraction: f64,
    pub interval_lower: f64,
    pub interval_upper: f64,
    pub interval_level: f64,
    pub interval_method: &'static str,
    pub estimated_input_supporting_fragments: f64,
    pub covered_bases: u64,
    pub representative_length: u64,
    pub coverage_fraction: f64,
    pub occupied_windows: u64,
    pub host_confounded_fragments: u64,
    pub cross_group_ambiguous_fragments: u64,
    pub integration_status: &'static str,
    pub split_events: u64,
    pub discordant_fragments: u64,
    pub limitation_codes: Vec<String>,
    pub sampling_target_score: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceReport {
    pub schema_id: &'static str,
    pub run: RunEvidence,
    pub target_signals: Vec<TargetSignal>,
    pub research_rows: Vec<[String; 24]>,
}

pub struct ReportInputs<'a> {
    pub sample_id: String,
    pub census: &'a InputCensus,
    pub design: SamplingDesign,
    pub selected_fragments: u64,
    pub prescreen_passed_fragments: u64,
    pub profile: AnalysisProfile,
    pub index_digest: String,
    pub unevaluable_fragments: u64,
    pub mode: RunMode,
}

pub fn build_evidence_report(
    inputs: ReportInputs<'_>,
    accumulator: EvidenceAccumulator,
) -> Result<EvidenceReport, String> {
    let family_size = accumulator.groups.len();
    let per_group_alpha = inputs.profile.familywise_interval_error / family_size as f64;
    let mut target_signals = Vec::new();
    for evidence in accumulator
        .groups
        .into_iter()
        .filter(|evidence| evidence.observed_or_indeterminate())
    {
        let interval = finite_population_interval(
            inputs.design.population_fragments,
            inputs.selected_fragments,
            evidence.supporting_selected_fragments,
            per_group_alpha,
        )?;
        let fraction = evidence.supporting_selected_fragments as f64
            / inputs.selected_fragments as f64
            * inputs.design.population_fragments as f64
            / inputs.census.fragments as f64;
        let covered_bases = evidence.covered_bases();
        target_signals.push(TargetSignal {
            target_group_id: evidence.group.target_group_id.clone(),
            representative_id: evidence.group.representative_id.clone(),
            representative_description: evidence.group.representative_description.clone(),
            member_ids: evidence.group.member_ids.clone(),
            evidence_status: evidence.evidence_status(),
            attribution_status: evidence.attribution_status(),
            supporting_selected_fragments: evidence.supporting_selected_fragments,
            selected_fragment_denominator: inputs.selected_fragments,
            attributed_fragment_fraction: fraction,
            interval_lower: interval.lower_count as f64 / inputs.census.fragments as f64,
            interval_upper: interval.upper_count as f64 / inputs.census.fragments as f64,
            interval_level: interval.level,
            interval_method: INTERVAL_METHOD,
            estimated_input_supporting_fragments: fraction * inputs.census.fragments as f64,
            covered_bases,
            representative_length: evidence.group.representative_length,
            coverage_fraction: covered_bases as f64 / evidence.group.representative_length as f64,
            occupied_windows: evidence.occupied_windows(inputs.profile.occupied_window_bins),
            host_confounded_fragments: evidence.host_confounded_fragments,
            cross_group_ambiguous_fragments: evidence.cross_group_ambiguous_fragments,
            integration_status: evidence.integration_status().as_str(),
            split_events: evidence.split_events,
            discordant_fragments: evidence.discordant_fragments,
            limitation_codes: Vec::new(),
            sampling_target_score: sampling_target_score(
                fraction,
                inputs.design.precision.minimum_fraction(),
            ),
        });
    }
    let reason_codes = (inputs.unevaluable_fragments > 0)
        .then(|| {
            format!(
                "TARGET_KMER_NOT_EVALUABLE_INPUT_FRAGMENTS={}",
                inputs.unevaluable_fragments
            )
        })
        .into_iter()
        .collect();
    let run = RunEvidence {
        sample_id: inputs.sample_id,
        analysis_status: if inputs.unevaluable_fragments == 0 {
            AnalysisStatus::ConformantComplete
        } else {
            AnalysisStatus::ConformantWithLimitations
        },
        reason_codes,
        input_mode: inputs.census.input_mode.into(),
        input_fragments: inputs.census.fragments,
        selected_fragments: inputs.selected_fragments,
        precision: inputs.design.precision.as_str(),
        sampling_population: inputs.mode.population(),
        selection_probability: inputs.design.selection_probability,
        minimum_relevant_fraction: inputs.design.precision.minimum_fraction(),
        familywise_miss_probability: inputs.profile.familywise_miss_probability,
        interval_level: 1.0 - inputs.profile.familywise_interval_error,
        target_family_size: family_size,
        prescreen_passed_fragments: inputs.prescreen_passed_fragments,
        aligned_fragments: accumulator.aligned_fragments,
        unassigned_fragments: accumulator.unassigned_fragments,
        profile_digest: inputs.profile.digest(),
        index_digest: inputs.index_digest,
        input_digest: inputs.census.input_digest.clone(),
        read_ends_per_fragment: inputs.census.read_ends_per_fragment,
        mode: inputs.mode.as_str(),
        sampling_population_fragments: inputs.design.population_fragments,
        sample_capacity: inputs.design.sample_capacity,
        sampling_expected_support: (inputs.selected_fragments > 0
            && inputs.design.population_fragments > 0)
            .then(|| {
                inputs.selected_fragments as f64
                    * inputs.design.precision.minimum_fraction()
                    * inputs.census.fragments as f64
                    / inputs.design.population_fragments as f64
            }),
    };
    let research_rows =
        build_research_rows(&run, &target_signals, inputs.profile.occupied_window_bins);
    Ok(EvidenceReport {
        schema_id: REPORT_SCHEMA,
        run,
        target_signals,
        research_rows,
    })
}

pub const RESEARCH_HEADER: [&str; 24] = [
    "sample_id",
    "reference_description",
    "reference_id",
    "support_rank",
    "support_fragments",
    "support_pct",
    "target_support_share_pct",
    "coverage_pct",
    "covered_bp",
    "reference_length_bp",
    "occupied_windows",
    "window_count",
    "support_ci_lower_pct",
    "support_ci_upper_pct",
    "split_support_fragments",
    "discordant_fragments",
    "reference_member_count",
    "reference_member_ids",
    "input_mode",
    "input_fragments",
    "selected_fragments",
    "total_target_support_fragments",
    "supported_reference_groups",
    "sampling_target_score",
];

const RESEARCH_LABELS: [&str; 24] = [
    "Sample ID",
    "Reference description",
    "Reference accession",
    "Support rank",
    "Supporting fragments",
    "Library abundance · %",
    "Target share · %",
    "Coverage breadth · %",
    "Covered bases · bp",
    "Reference length · bp",
    "Occupied windows",
    "Total windows",
    "Interval lower bound · %",
    "Interval upper bound · %",
    "Split-support fragments",
    "Discordant fragments",
    "Equivalent reference count",
    "Equivalent reference IDs",
    "Input mode",
    "Input fragments",
    "Selected fragments",
    "Total target-support fragments",
    "Supported reference groups",
    "Target score",
];

fn build_research_rows(
    run: &RunEvidence,
    signals: &[TargetSignal],
    windows: usize,
) -> Vec<[String; 24]> {
    let mut supported = signals
        .iter()
        .filter(|signal| signal.supporting_selected_fragments > 0)
        .collect::<Vec<_>>();
    supported.sort_by(|left, right| {
        right
            .supporting_selected_fragments
            .cmp(&left.supporting_selected_fragments)
            .then_with(|| left.representative_id.cmp(&right.representative_id))
    });
    let total: u64 = supported
        .iter()
        .map(|signal| signal.supporting_selected_fragments)
        .sum();
    let common = || {
        let mut row = std::array::from_fn(|_| String::new());
        row[0] = run.sample_id.clone();
        row[18] = run.input_mode.clone();
        row[19] = run.input_fragments.to_string();
        row[20] = run.selected_fragments.to_string();
        row[21] = total.to_string();
        row[22] = supported.len().to_string();
        row
    };
    if supported.is_empty() {
        return vec![common()];
    }
    supported
        .iter()
        .enumerate()
        .map(|(rank, signal)| {
            let mut row = common();
            row[1] = signal.representative_description.clone();
            row[2] = signal.representative_id.clone();
            row[3] = (rank + 1).to_string();
            row[4] = signal.supporting_selected_fragments.to_string();
            row[5] = research_decimal_precise(signal.attributed_fragment_fraction * 100.0);
            row[6] = research_decimal(
                signal.supporting_selected_fragments as f64 / total as f64 * 100.0,
            );
            row[7] = research_decimal(signal.coverage_fraction * 100.0);
            row[8] = signal.covered_bases.to_string();
            row[9] = signal.representative_length.to_string();
            row[10] = signal.occupied_windows.to_string();
            row[11] = windows.to_string();
            row[12] = interval_pct(signal.interval_lower, false);
            row[13] = interval_pct(signal.interval_upper, true);
            row[14] = signal.split_events.to_string();
            row[15] = signal.discordant_fragments.to_string();
            row[16] = signal.member_ids.len().to_string();
            row[17] = signal.member_ids.join(";");
            row[23] = signal.sampling_target_score.to_string();
            row
        })
        .collect()
}

fn research_decimal(value: f64) -> String {
    format!("{value:.6}")
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn research_decimal_precise(value: f64) -> String {
    format!("{value:.10}")
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}

fn interval_pct(fraction: f64, upper: bool) -> String {
    // One displayed unit is 1e-10 percent. Include multiplication error
    // before rounding outward so formatting cannot narrow the original interval.
    let scaled = fraction * 1_000_000_000_000.0;
    let units = if upper {
        scaled.next_up().ceil()
    } else {
        scaled.next_down().floor()
    };
    research_decimal_precise(units.clamp(0.0, 1_000_000_000_000.0) / 10_000_000_000.0)
}

fn csv_text(report: &EvidenceReport) -> String {
    let mut output = String::from("\u{feff}");
    output.push_str(&RESEARCH_HEADER.join(","));
    output.push_str("\r\n");
    for row in &report.research_rows {
        push_csv_row(&mut output, row);
    }
    output
}

pub fn write_report_csv(path: &Path, report: &EvidenceReport) -> Result<(), String> {
    atomic_write(path, csv_text(report).as_bytes())
}

fn csv_data_uri(report: &EvidenceReport) -> String {
    const HEX: &[u8; 16] = b"0123456789ABCDEF";
    let csv = csv_text(report);
    let mut uri = String::with_capacity(28 + csv.len() * 3);
    uri.push_str("data:text/csv;charset=utf-8,");
    for byte in csv.bytes() {
        uri.push('%');
        uri.push(HEX[(byte >> 4) as usize] as char);
        uri.push(HEX[(byte & 15) as usize] as char);
    }
    uri
}

pub fn write_report_html(path: &Path, report: &EvidenceReport) -> Result<(), String> {
    let rows = &report.research_rows;
    let visible = |row: &&[String; 24]| !row[2].is_empty();
    let count = rows.iter().filter(visible).take(20).count();
    let all_count = rows.iter().filter(|row| !row[2].is_empty()).count();
    let total: u64 = report
        .target_signals
        .iter()
        .map(|s| s.supporting_selected_fragments)
        .sum();
    let mut result_list = String::new();
    for row in rows.iter().filter(visible).take(20) {
        result_list.push_str(&research_html(row));
    }
    if count == 0 {
        let explanation = if report.run.selected_fragments == 0 {
            "No target candidates were retained."
        } else if report.run.prescreen_passed_fragments == 0 {
            "No screened fragments passed the target Bloom filter."
        } else {
            "No selected fragments were attributed to a target reference group."
        };
        let sample_fields = if all_count == 0 {
            research_field_list(&rows[0])
        } else {
            String::new()
        };
        result_list = format!("<div class=\"empty-state\"><h3>No attributed target signals</h3><p>{explanation}</p>{sample_fields}</div>");
    }
    let mut audit = String::new();
    for signal in report.target_signals.iter().filter(|s| {
        rows.iter()
            .filter(|row| !row[2].is_empty())
            .take(20)
            .any(|row| row[2] == s.representative_id)
    }) {
        audit.push_str(&format!(
            "<section class=\"target-signal\" data-group=\"{}\"><h4>{}</h4><table><tbody>{}</tbody></table></section>",
            html(&signal.target_group_id), html(&signal.representative_id),
            visible_field_rows(&target_fields(signal))
        ));
    }
    let values = [
        ("STYLE", include_str!("report.css").to_string()),
        ("CSV_DOWNLOAD", csv_data_uri(report)),
        ("SCHEMA", report.schema_id.to_string()),
        ("SAMPLE", html(&report.run.sample_id)),
        ("INPUT_MODE", html(&report.run.input_mode)),
        ("INPUT_COUNT", report.run.input_fragments.to_string()),
        ("SELECTED_COUNT", report.run.selected_fragments.to_string()),
        ("TOTAL_SUPPORT", total.to_string()),
        ("GROUP_COUNT", count.to_string()),
        ("ALL_GROUP_COUNT", all_count.to_string()),
        ("SAMPLED_LABEL", if report.run.mode == "full" { "Sampled candidates" } else { "Sampled input fragments" }.into()),
        ("METHOD", if report.run.mode == "full" {
            "All input fragments undergo a target 21-mer Bloom prescreen. Candidate fragments are selected by deterministic bottom-k sampling before HOST+TARGET competitive alignment."
        } else {
            "Input fragments are selected by deterministic bottom-k sampling before the target 21-mer Bloom prescreen and HOST+TARGET competitive alignment."
        }.into()),
        ("RESULT_LIST", result_list),
        ("RUN_FIELDS", visible_field_rows(&run_fields(report))),
        ("AUDIT_SIGNALS", audit),
    ];
    // Only parse the static template: data containing template delimiters stays data.
    let mut parts = include_str!("report.html").split("{{");
    let mut document = parts.next().unwrap_or_default().to_string();
    for part in parts {
        let (key, rest) = part
            .split_once("}}")
            .ok_or_else(|| "Unclosed report template field".to_string())?;
        let (_, value) = values
            .iter()
            .find(|(name, _)| *name == key)
            .ok_or_else(|| format!("Unknown report template field: {key}"))?;
        document.push_str(value);
        document.push_str(rest);
    }
    atomic_write(path, document.as_bytes())
}

fn research_field_list(row: &[String; 24]) -> String {
    let mut fields = String::from("<dl class=\"research-fields\" data-research-row>");
    for ((key, label), value) in RESEARCH_HEADER.iter().zip(RESEARCH_LABELS).zip(row) {
        fields.push_str(&format!(
            "<div><dt>{label}</dt><dd data-field=\"{key}\">{}</dd></div>",
            html(value)
        ));
    }
    fields.push_str("</dl>");
    fields
}

fn research_html(row: &[String; 24]) -> String {
    let title = if row[1].is_empty() { &row[2] } else { &row[1] };
    format!(
        r#"<details class="signal">
<summary><span class="rank"><span class="sr-only">Rank </span>{rank}</span>
<span class="reference"><span class="reference-name">{title}</span><span class="accession">{reference}</span></span>
<span class="metric"><span class="mobile-label">Fragments</span><strong>{support}</strong></span>
<span class="metric"><span class="mobile-label">Library interval · %</span><strong class="abundance-interval">[{lower}%, {upper}%]</strong><span class="point-estimate">Estimate {pct}%</span></span>
<span class="metric coverage"><span class="mobile-label">Coverage</span><strong>{coverage}%</strong><span class="coverage-track" aria-hidden="true"><span style="width:{coverage}%"></span></span></span>
<span class="metric target-score"><span class="mobile-label">Target score</span><strong>{score:+.6}</strong></span>
<svg class="chevron" viewBox="0 0 16 16" aria-hidden="true"><path d="m4 6 4 4 4-4"/></svg></summary>
<div class="signal-detail"><h3>Research data <span>{reference}</span></h3>{fields}</div></details>"#,
        rank = html(&row[3]),
        title = html(title),
        reference = html(&row[2]),
        support = html(&row[4]),
        pct = html(&row[5]),
        lower = html(&row[12]),
        upper = html(&row[13]),
        coverage = html(&row[7]),
        score = row[23].parse::<f64>().unwrap(),
        fields = research_field_list(row)
    )
}

fn sampling_target_score(fraction: f64, target: f64) -> f64 {
    (fraction - target) / (fraction + target)
}

fn run_fields(report: &EvidenceReport) -> Vec<(&'static str, String)> {
    let run = &report.run;
    vec![
        ("schema_id", report.schema_id.into()),
        ("record_type", "RUN".into()),
        ("sample_id", run.sample_id.clone()),
        ("analysis_status", run.analysis_status.as_str().into()),
        ("reason_codes", run.reason_codes.join(";")),
        ("input_mode", run.input_mode.clone()),
        ("input_fragments", run.input_fragments.to_string()),
        ("selected_fragments", run.selected_fragments.to_string()),
        ("precision", run.precision.into()),
        ("mode", run.mode.into()),
        ("sample_capacity", run.sample_capacity.to_string()),
        (
            "sampling_population_fragments",
            run.sampling_population_fragments.to_string(),
        ),
        (
            "sampling_expected_support",
            run.sampling_expected_support
                .map(format_float)
                .unwrap_or_default(),
        ),
        ("sampling_population", run.sampling_population.into()),
        (
            "selection_probability",
            format_float(run.selection_probability),
        ),
        (
            "minimum_relevant_fraction",
            format_float(run.minimum_relevant_fraction),
        ),
        (
            "familywise_miss_probability",
            format_float(run.familywise_miss_probability),
        ),
        ("interval_level", format_float(run.interval_level)),
        ("target_family_size", run.target_family_size.to_string()),
        (
            "prescreen_passed_fragments",
            run.prescreen_passed_fragments.to_string(),
        ),
        ("aligned_fragments", run.aligned_fragments.to_string()),
        ("unassigned_fragments", run.unassigned_fragments.to_string()),
        ("profile_digest", run.profile_digest.clone()),
        ("index_digest", run.index_digest.clone()),
        ("input_digest", run.input_digest.clone()),
        (
            "read_ends_per_fragment",
            run.read_ends_per_fragment.to_string(),
        ),
    ]
}

fn target_fields(signal: &TargetSignal) -> Vec<(&'static str, String)> {
    vec![
        ("target_group_id", signal.target_group_id.clone()),
        ("record_type", "TARGET_SIGNAL".into()),
        ("representative_id", signal.representative_id.clone()),
        (
            "representative_description",
            signal.representative_description.clone(),
        ),
        ("member_ids", signal.member_ids.join(";")),
        ("evidence_status", signal.evidence_status.as_str().into()),
        (
            "attribution_status",
            signal.attribution_status.as_str().into(),
        ),
        (
            "supporting_selected_fragments",
            signal.supporting_selected_fragments.to_string(),
        ),
        (
            "selected_fragment_denominator",
            signal.selected_fragment_denominator.to_string(),
        ),
        (
            "attributed_fragment_fraction",
            format_float(signal.attributed_fragment_fraction),
        ),
        ("interval_lower", format_float(signal.interval_lower)),
        ("interval_upper", format_float(signal.interval_upper)),
        ("interval_level", format_float(signal.interval_level)),
        ("interval_method", signal.interval_method.into()),
        (
            "estimated_input_supporting_fragments",
            format_float(signal.estimated_input_supporting_fragments),
        ),
        ("covered_bases", signal.covered_bases.to_string()),
        (
            "representative_length",
            signal.representative_length.to_string(),
        ),
        ("coverage_fraction", format_float(signal.coverage_fraction)),
        ("occupied_windows", signal.occupied_windows.to_string()),
        (
            "host_confounded_fragments",
            signal.host_confounded_fragments.to_string(),
        ),
        (
            "cross_group_ambiguous_fragments",
            signal.cross_group_ambiguous_fragments.to_string(),
        ),
        ("integration_status", signal.integration_status.into()),
        ("split_events", signal.split_events.to_string()),
        (
            "discordant_fragments",
            signal.discordant_fragments.to_string(),
        ),
        ("limitation_codes", signal.limitation_codes.join(";")),
        (
            "sampling_target_score",
            signal.sampling_target_score.to_string(),
        ),
    ]
}

fn visible_field_rows(fields: &[(&str, String)]) -> String {
    fields
        .iter()
        .map(|(field, value)| {
            format!(
                "<tr data-field=\"{}\"><th>{}</th><td>{}</td></tr>",
                field,
                field,
                html(value)
            )
        })
        .collect()
}

fn push_csv_row(output: &mut String, row: &[String]) {
    output.push_str(
        &row.iter()
            .map(|cell| csv(cell))
            .collect::<Vec<_>>()
            .join(","),
    );
    output.push_str("\r\n");
}
fn csv(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}
fn html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}
fn format_float(value: f64) -> String {
    format!("{value:.17}")
        .trim_end_matches('0')
        .trim_end_matches('.')
        .to_string()
}
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn candidate_sampling_estimates_the_original_library() {
        use crate::index::reference::{build_reference_groups, FastaRecord};
        use crate::sampling::Precision;
        let groups = build_reference_groups(&[FastaRecord {
            id: "virus".into(),
            description: "Virus reference".into(),
            sequence: b"ACGT".repeat(100),
        }]);
        let census = InputCensus {
            input_mode: "PE",
            fragments: 10_000,
            input_digest: "input".into(),
            read_ends_per_fragment: 2,
        };
        let build = |selected: u64, support: u64| {
            let mut accumulator = EvidenceAccumulator::new(&groups);
            accumulator.groups[0].supporting_selected_fragments = support;
            build_evidence_report(
                ReportInputs {
                    sample_id: "sample".into(),
                    census: &census,
                    design: SamplingDesign {
                        precision: Precision::Standard,
                        population_fragments: 1_000,
                        selection_probability: 0.1,
                        sample_capacity: selected,
                    },
                    selected_fragments: selected,
                    prescreen_passed_fragments: 1_000,
                    profile: AnalysisProfile::FROZEN,
                    index_digest: "index".into(),
                    unevaluable_fragments: 0,
                    mode: RunMode::Full,
                },
                accumulator,
            )
            .unwrap()
        };
        let report = build(100, 10);
        let signal = &report.target_signals[0];
        assert!((signal.attributed_fragment_fraction - 0.01).abs() < 1e-15);
        assert_eq!(signal.estimated_input_supporting_fragments, 100.0);
        let interval = finite_population_interval(1_000, 100, 10, 0.05).unwrap();
        assert_eq!(
            signal.interval_lower,
            interval.lower_count as f64 / 10_000.0
        );
        assert_eq!(
            signal.interval_upper,
            interval.upper_count as f64 / 10_000.0
        );
        assert_eq!(report.research_rows[0][5], "1");
        let empty = build(0, 0);
        assert!(empty.target_signals.is_empty());
        assert_eq!(empty.research_rows.len(), 1);
        assert!(empty.research_rows[0][12].is_empty());
    }

    #[test]
    fn target_score_tracks_library_target_in_both_populations() {
        let delta = 1e-5;
        assert_eq!(sampling_target_score(0.0, delta), -1.0);
        assert_eq!(sampling_target_score(delta, delta), 0.0);
        for ratio in [0.5, 2.0, 10.0] {
            let actual = sampling_target_score(delta * ratio, delta);
            assert!((actual - (ratio - 1.0) / (ratio + 1.0)).abs() < 1e-15);
        }
        let mut previous = -1.0;
        for x in 1..100 {
            let score = sampling_target_score(x as f64 / 1_000_000.0, delta);
            assert!(score > previous && score < 1.0);
            previous = score;
        }
        let screen = sampling_target_score(20.0 / 1_000_000.0, delta);
        let full = sampling_target_score(
            (1_000_000.0 / 100_000_000.0) * (2000.0 / 1_000_000.0),
            delta,
        );
        assert!((screen - full).abs() < 1e-15);
    }

    #[test]
    fn top_twenty_includes_negative_scores_and_embeds_every_csv_row() {
        use crate::index::reference::{build_reference_groups, FastaRecord};
        use crate::sampling::Precision;
        use scraper::{Html, Selector};
        let groups = build_reference_groups(
            &(0..25)
                .map(|i| FastaRecord {
                    id: format!("virus-{i}"),
                    description: format!("Reference {i}"),
                    sequence: b"ACGT".repeat(100 + i),
                })
                .collect::<Vec<_>>(),
        );
        let census = InputCensus {
            input_mode: "SE",
            fragments: 100_000_000,
            input_digest: "input".into(),
            read_ends_per_fragment: 1,
        };
        let mut accumulator = EvidenceAccumulator::new(&groups);
        for (i, group) in accumulator.groups.iter_mut().enumerate() {
            group.supporting_selected_fragments = (i + 1) as u64;
        }
        let report = build_evidence_report(
            ReportInputs {
                sample_id: "top20".into(),
                census: &census,
                design: SamplingDesign {
                    precision: Precision::Standard,
                    population_fragments: census.fragments,
                    selection_probability: 0.01,
                    sample_capacity: 1_000_000,
                },
                selected_fragments: 1_000_000,
                prescreen_passed_fragments: 325,
                profile: AnalysisProfile::FROZEN,
                index_digest: "index".into(),
                unevaluable_fragments: 0,
                mode: RunMode::Screen,
            },
            accumulator,
        )
        .unwrap();
        let path =
            std::env::temp_dir().join(format!("viroflash-top20-{}.html", std::process::id()));
        write_report_html(&path, &report).unwrap();
        let text = std::fs::read_to_string(&path).unwrap();
        std::fs::remove_file(path).unwrap();
        let html = Html::parse_document(&text);
        assert_eq!(report.research_rows.len(), 25);
        let rows = html
            .select(&Selector::parse("dl[data-research-row]").unwrap())
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 20);
        for (actual, expected) in rows.iter().zip(report.research_rows.iter()) {
            let values = actual
                .select(&Selector::parse("dd[data-field]").unwrap())
                .map(|v| v.text().collect::<String>())
                .collect::<Vec<_>>();
            assert_eq!(values, expected.as_slice());
        }
        assert!(report.research_rows[..20]
            .iter()
            .any(|r| r[23].parse::<f64>().unwrap() < 0.0));
        let href = html
            .select(&Selector::parse("a[download]").unwrap())
            .next()
            .unwrap()
            .value()
            .attr("href")
            .unwrap();
        assert_eq!(href, csv_data_uri(&report));
    }

    #[test]
    fn displayed_intervals_round_outward_at_decimal_and_float_boundaries() {
        for value in [
            0.0_f64,
            1e-18,
            1e-12,
            1e-6,
            0.000_000_811_627_703_13,
            0.1,
            0.5,
            1.0,
        ] {
            for fraction in [value.next_down().max(0.0), value, value.next_up().min(1.0)] {
                let lower = interval_pct(fraction, false).parse::<f64>().unwrap() / 100.0;
                let upper = interval_pct(fraction, true).parse::<f64>().unwrap() / 100.0;
                assert!(
                    lower <= fraction && fraction <= upper,
                    "{lower} <= {fraction} <= {upper}"
                );
                assert!(upper - lower <= 2.1e-12);
            }
        }
    }
}
