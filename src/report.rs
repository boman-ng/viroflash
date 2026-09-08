use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

use crate::analysis_profile::AnalysisProfile;
use crate::evidence::{
    finite_population_interval, AttributionStatus, EvidenceAccumulator, EvidenceStatus,
    INTERVAL_METHOD,
};
use crate::fastq_input::InputCensus;
use crate::sampling_design::SamplingDesign;
use serde::Serialize;

pub const REPORT_SCHEMA: &str = "viroflash.evidence-report.v1";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "SCREAMING_SNAKE_CASE")]
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

#[derive(Debug, Clone, Serialize)]
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

#[derive(Debug, Clone, Serialize)]
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

#[derive(Debug, Clone, Serialize)]
pub struct EvidenceReport {
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
    let run = &report.run;
    let mut row = vec![String::new(); HEADER.len()];
    set(&mut row, "schema_id", REPORT_SCHEMA);
    set(&mut row, "record_type", "RUN");
    set(&mut row, "sample_id", &run.sample_id);
    set(&mut row, "analysis_status", run.analysis_status.as_str());
    set(&mut row, "reason_codes", &run.reason_codes.join(";"));
    set(&mut row, "input_mode", &run.input_mode);
    for (field, value) in [
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
        (
            "read_ends_per_fragment",
            run.read_ends_per_fragment.to_string(),
        ),
    ] {
        set(&mut row, field, &value);
    }
    set(&mut row, "profile_digest", &run.profile_digest);
    set(&mut row, "index_digest", &run.index_digest);
    set(&mut row, "input_digest", &run.input_digest);
    push_csv_row(&mut output, &row);
    for signal in &report.target_signals {
        let mut row = vec![String::new(); HEADER.len()];
        set(&mut row, "record_type", "TARGET_SIGNAL");
        set(&mut row, "target_group_id", &signal.target_group_id);
        set(&mut row, "representative_id", &signal.representative_id);
        set(&mut row, "member_ids", &signal.member_ids.join(";"));
        set(&mut row, "evidence_status", signal.evidence_status.as_str());
        set(
            &mut row,
            "attribution_status",
            signal.attribution_status.as_str(),
        );
        for (field, value) in [
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
            ("split_events", signal.split_events.to_string()),
            (
                "discordant_fragments",
                signal.discordant_fragments.to_string(),
            ),
        ] {
            set(&mut row, field, &value);
        }
        set(&mut row, "interval_method", INTERVAL_METHOD);
        set(&mut row, "integration_status", signal.integration_status);
        set(
            &mut row,
            "limitation_codes",
            &signal.limitation_codes.join(";"),
        );
        push_csv_row(&mut output, &row);
    }
    atomic_write(path, output.as_bytes())
}

pub fn write_report_html(path: &Path, report: &EvidenceReport) -> Result<(), String> {
    let mut rows = String::new();
    for signal in &report.target_signals {
        rows.push_str(&format!("<tr data-group=\"{}\"><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td data-fraction=\"{}\">{}</td><td data-lower=\"{}\" data-upper=\"{}\">[{}, {}]</td><td>{}</td><td>{}</td><td>{}/{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
            html(&signal.target_group_id), html(&signal.target_group_id), html(&signal.representative_id), html(&signal.member_ids.join(";")), signal.supporting_selected_fragments,
            format_float(signal.attributed_fragment_fraction), format_float(signal.attributed_fragment_fraction), format_float(signal.interval_lower), format_float(signal.interval_upper),
            format_float(signal.interval_lower), format_float(signal.interval_upper), signal.evidence_status.as_str(), signal.attribution_status.as_str(), signal.covered_bases,
            signal.representative_length, signal.occupied_windows, signal.host_confounded_fragments, signal.cross_group_ambiguous_fragments, signal.split_events,
            signal.discordant_fragments, signal.integration_status));
    }
    let embedded_report = serde_json::to_string(report)
        .map_err(|error| format!("Cannot serialize HTML EvidenceReport model: {error}"))?
        .replace('&', "\\u0026")
        .replace('<', "\\u003c")
        .replace('>', "\\u003e");
    let html_document = format!(
        r#"<!doctype html><html><head><meta charset="utf-8"><title>Viroflash evidence report</title><style>body{{font:16px system-ui;max-width:1100px;margin:2rem auto;padding:0 1rem;color:#17202a}}table{{border-collapse:collapse;width:100%}}th,td{{border:1px solid #ccd1d1;padding:.45rem;text-align:left}}code{{overflow-wrap:anywhere}}.boundary{{background:#fff4d6;padding:1rem;border-left:4px solid #d68910}}</style></head><body>
<h1>Viroflash Evidence Report</h1><section><h2>Interpretation boundary</h2><p class="boundary">Profile-attributed fragment fraction measures input-library fragments attributed under the frozen Viroflash profile. It is not viral load, absolute quantitation, or a clinical positive/negative conclusion.</p></section>
<section><h2>Run integrity</h2><dl><dt>Status</dt><dd>{}</dd><dt>Reason codes</dt><dd>{}</dd><dt>Sample</dt><dd>{}</dd><dt>Input</dt><dd>{} fragments; {} selected at probability {}</dd><dt>Profile digest</dt><dd><code>{}</code></dd><dt>Index digest</dt><dd><code>{}</code></dd><dt>Input digest</dt><dd><code>{}</code></dd></dl></section>
<section><h2>Observed target signals</h2><table><thead><tr><th>Group</th><th>Representative</th><th>Members</th><th>Supporting fragments</th><th>Fraction</th><th>Simultaneous interval</th><th>Evidence</th><th>Attribution</th><th>Coverage</th><th>Windows</th><th>Host-confounded</th><th>Cross-group ambiguous</th><th>Split</th><th>Discordant</th><th>Integration</th></tr></thead><tbody>{}</tbody></table></section>
<section><h2>Evidence detail</h2><p>The table reports every model value for host conflict, cross-group ambiguity, coverage, ten occupied windows, split events, discordance, and integration. These diagnostics do not create a decision cutoff.</p></section>
<section><h2>Methods and limitations</h2><p>FASTQ is fully validated and counted in pass 1. Pass 2 applies deterministic BLAKE3 Bernoulli fragment inclusion, a target-only 21-mer Bloom workload gate, and HOST+TARGET competitive minimap2 alignment. Intervals use simultaneous equal-tailed exact hypergeometric inversion. References absent from the index remain unmodeled.</p></section><script id="evidence-report" type="application/json">{}</script></body></html>"#,
        report.run.analysis_status.as_str(),
        html(&report.run.reason_codes.join(";")),
        html(&report.run.sample_id),
        report.run.input_fragments,
        report.run.selected_fragments,
        format_float(report.run.selection_probability),
        report.run.profile_digest,
        report.run.index_digest,
        report.run.input_digest,
        rows,
        embedded_report
    );
    atomic_write(path, html_document.as_bytes())
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
fn atomic_write(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let temporary = path.with_extension(format!(
        "{}.part.{}",
        path.extension()
            .and_then(|value| value.to_str())
            .unwrap_or("tmp"),
        std::process::id()
    ));
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| format!("Cannot create {}: {error}", temporary.display()))?;
    let result = file
        .write_all(bytes)
        .and_then(|_| file.flush())
        .map_err(|error| format!("Cannot write {}: {error}", temporary.display()))
        .and_then(|_| {
            std::fs::rename(&temporary, path)
                .map_err(|error| format!("Cannot finalize {}: {error}", path.display()))
        });
    if result.is_err() {
        let _ = std::fs::remove_file(temporary);
    }
    result
}
