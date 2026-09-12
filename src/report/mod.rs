use crate::output::atomic_write;
use std::path::Path;

use crate::evidence::{
    finite_population_interval, AttributionStatus, EvidenceAccumulator, EvidenceStatus,
    INTERVAL_METHOD,
};
use crate::fastq::InputCensus;
use crate::profile::AnalysisProfile;
use crate::sampling::SamplingDesign;

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
}

#[derive(Debug, Clone, PartialEq)]
pub struct TargetSignal {
    pub target_group_id: String,
    pub representative_id: String,
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
}

#[derive(Debug, Clone, PartialEq)]
pub struct EvidenceReport {
    pub schema_id: &'static str,
    pub run: RunEvidence,
    pub target_signals: Vec<TargetSignal>,
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
            inputs.census.fragments,
            inputs.selected_fragments,
            evidence.supporting_selected_fragments,
            per_group_alpha,
        )?;
        let fraction =
            evidence.supporting_selected_fragments as f64 / inputs.selected_fragments as f64;
        let covered_bases = evidence.covered_bases();
        target_signals.push(TargetSignal {
            target_group_id: evidence.group.target_group_id.clone(),
            representative_id: evidence.group.representative_id.clone(),
            member_ids: evidence.group.member_ids.clone(),
            evidence_status: evidence.evidence_status(),
            attribution_status: evidence.attribution_status(),
            supporting_selected_fragments: evidence.supporting_selected_fragments,
            selected_fragment_denominator: inputs.selected_fragments,
            attributed_fragment_fraction: fraction,
            interval_lower: interval.lower,
            interval_upper: interval.upper,
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
        });
    }
    let reason_codes = (inputs.unevaluable_fragments > 0)
        .then(|| {
            format!(
                "TARGET_KMER_NOT_EVALUABLE_SELECTED_FRAGMENTS={}",
                inputs.unevaluable_fragments
            )
        })
        .into_iter()
        .collect();
    Ok(EvidenceReport {
        schema_id: REPORT_SCHEMA,
        run: RunEvidence {
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
            selection_probability: inputs.design.selection_probability,
            minimum_relevant_fraction: inputs.profile.minimum_relevant_fraction,
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
        },
        target_signals,
    })
}

const HEADER: &[&str] = &[
    "schema_id",
    "record_type",
    "sample_id",
    "analysis_status",
    "reason_codes",
    "input_mode",
    "input_fragments",
    "selected_fragments",
    "selection_probability",
    "minimum_relevant_fraction",
    "familywise_miss_probability",
    "interval_level",
    "target_family_size",
    "prescreen_passed_fragments",
    "aligned_fragments",
    "unassigned_fragments",
    "profile_digest",
    "index_digest",
    "input_digest",
    "read_ends_per_fragment",
    "target_group_id",
    "representative_id",
    "member_ids",
    "evidence_status",
    "attribution_status",
    "supporting_selected_fragments",
    "selected_fragment_denominator",
    "attributed_fragment_fraction",
    "interval_lower",
    "interval_upper",
    "interval_method",
    "estimated_input_supporting_fragments",
    "covered_bases",
    "representative_length",
    "coverage_fraction",
    "occupied_windows",
    "host_confounded_fragments",
    "cross_group_ambiguous_fragments",
    "integration_status",
    "split_events",
    "discordant_fragments",
    "limitation_codes",
];

pub fn write_report_csv(path: &Path, report: &EvidenceReport) -> Result<(), String> {
    let mut output = String::new();
    output.push_str(&HEADER.join(","));
    output.push('\n');
    let mut row = vec![String::new(); HEADER.len()];
    for (field, value) in run_fields(report) {
        set(&mut row, field, &value);
    }
    push_csv_row(&mut output, &row);
    for signal in &report.target_signals {
        let mut row = vec![String::new(); HEADER.len()];
        for (field, value) in target_fields(signal) {
            set(&mut row, field, &value);
        }
        push_csv_row(&mut output, &row);
    }
    atomic_write(path, output.as_bytes())
}

pub fn write_report_html(path: &Path, report: &EvidenceReport) -> Result<(), String> {
    let run_rows = visible_field_rows(&run_fields(report));
    let mut signals = String::new();
    for signal in &report.target_signals {
        signals.push_str(&format!(
            "<section class=\"target-signal\" data-group=\"{}\"><h3>{}</h3><table><tbody>{}</tbody></table></section>",
            html(&signal.target_group_id),
            html(&signal.target_group_id),
            visible_field_rows(&target_fields(signal))
        ));
    }
    let html_document = format!(
        r#"<!doctype html><html><head><meta charset="utf-8"><title>Viroflash evidence report</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17202a}}table{{border-collapse:collapse;width:100%;margin-bottom:1.5rem}}th,td{{border:1px solid #ccd1d1;padding:.45rem;text-align:left;overflow-wrap:anywhere}}th{{width:22rem}}.boundary{{background:#fff4d6;padding:1rem;border-left:4px solid #d68910}}</style></head><body>
<h1>Viroflash Evidence Report</h1><section><h2>Interpretation boundary</h2><p class="boundary">Profile-attributed fragment fraction measures input-library fragments attributed under the frozen Viroflash profile. It is not viral load, absolute quantitation, or a clinical positive/negative conclusion.</p></section>
<section id="run-integrity"><h2>Run integrity</h2><table><tbody>{}</tbody></table></section>
<section id="observed-target-signals"><h2>Observed target signals</h2>{}</section>
<section><h2>Evidence detail</h2><p>The table reports every model value for host conflict, cross-group ambiguity, coverage, ten occupied windows, split events, discordance, and integration. These diagnostics do not create a decision cutoff.</p></section>
<section><h2>Methods and limitations</h2><p>FASTQ is fully validated and counted in pass 1. Pass 2 applies deterministic BLAKE3 Bernoulli fragment inclusion, a target-only 21-mer Bloom workload gate, and HOST+TARGET competitive minimap2 alignment. Intervals use simultaneous equal-tailed exact hypergeometric inversion. References absent from the index remain unmodeled.</p></section></body></html>"#,
        run_rows, signals
    );
    atomic_write(path, html_document.as_bytes())
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

fn set(row: &mut [String], field: &str, value: &str) {
    row[HEADER
        .iter()
        .position(|candidate| *candidate == field)
        .unwrap()] = value.to_string();
}
fn push_csv_row(output: &mut String, row: &[String]) {
    output.push_str(
        &row.iter()
            .map(|cell| csv(cell))
            .collect::<Vec<_>>()
            .join(","),
    );
    output.push('\n');
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
