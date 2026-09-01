//! Self-contained HTML and flat CSV views for human review of candidate evidence.
//!
//! These files are presentation layers over the stable JSON contract. They do not recompute
//! candidates, alter decisions, or imply sample-level or clinical conclusions.

use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::report::{Candidate, SamplingReport, RESULT_SCHEMA};

pub const CSV_SCHEMA: &str = "viroflash.candidates.csv.v1";

const CSV_HEADER: &str = "csv_schema,result_schema,sample_id,sample_conclusion,qc_status,qc_issues,input_pairs,prescreen_pairs,selected_pairs,inclusion_probability,discovery_pairs,validation_pairs,validation_read_sides,map_errors,audit_map_error_reads,audit_overflows,target_validation_unassigned,decoy_validation_unassigned,test_family_size,reported_candidates,threads,k,index_source,index_format_version,manifest_blake3,candidate_id,representative,resolution_level,hypothesis_members,member_attribution,discovery_reads,index_member_count,target_exposure_bases,decision,decision_reasons,adjusted_p_value,log10_adjusted_p_value,adjusted_p_status,model_adjusted_p_max,p_value,log10_p_value,p_value_status,p_resolution_floor,p_resolution_floor_reached,validation_read_ends,covered_bases,contig_length,coverage_breadth,coverage_min,distributed_windows,total_windows,min_distributed_windows,end_rpm,expected_hits,depth_fold,stratum,stratum_decoy_count,background_status,background_scope,background_reads,background_reference_bases,background_cross_stratum_reads,integration_status,split_events,supported_sites,discordant_pairs,plus_strand,minus_strand";

pub struct HumanReportInput<'a> {
    pub sample: &'a str,
    pub threads: usize,
    pub k: usize,
    pub input_pairs: u64,
    pub prescreen_pairs: u64,
    pub map_errors: u64,
    pub index_source: &'a str,
    pub index_format_version: u32,
    pub manifest_blake3: &'a str,
    pub sampling: &'a SamplingReport,
    pub test_family_size: usize,
    pub candidates: &'a [Candidate],
}

pub fn write_human_report(
    out_prefix: &Path,
    input: &HumanReportInput<'_>,
) -> Result<(PathBuf, PathBuf), String> {
    let html_path = PathBuf::from(format!("{}.html", out_prefix.display()));
    let csv_path = PathBuf::from(format!("{}.csv", out_prefix.display()));
    write_csv(&csv_path, input)?;
    write_html(&html_path, &csv_path, input)?;
    Ok((html_path, csv_path))
}

fn write_csv(path: &Path, input: &HumanReportInput<'_>) -> Result<(), String> {
    let file = File::create(path).map_err(|e| format!("Cannot create {}: {e}", path.display()))?;
    let mut out = BufWriter::new(file);
    writeln!(out, "{CSV_HEADER}")
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;

    let qc_issues = qc_issue_codes(input);
    for candidate in input.candidates {
        let values = [
            csv_cell(CSV_SCHEMA),
            csv_cell(RESULT_SCHEMA),
            csv_cell(input.sample),
            "not_computed".to_string(),
            "NOT_EVALUATED".to_string(),
            csv_cell(&json_str_array(&qc_issues)),
            input.input_pairs.to_string(),
            input.prescreen_pairs.to_string(),
            input.sampling.selected_pairs.to_string(),
            input.sampling.inclusion_probability.to_string(),
            input.sampling.discovery_pairs.to_string(),
            input.sampling.validation_pairs.to_string(),
            input.sampling.validation_read_sides.to_string(),
            input.map_errors.to_string(),
            input.sampling.audit_map_error_reads.to_string(),
            input.sampling.audit_overflows.to_string(),
            input.sampling.target_validation_unassigned.to_string(),
            input.sampling.decoy_validation_unassigned.to_string(),
            input.test_family_size.to_string(),
            input.candidates.len().to_string(),
            input.threads.to_string(),
            input.k.to_string(),
            csv_cell(input.index_source),
            input.index_format_version.to_string(),
            csv_cell(input.manifest_blake3),
            csv_cell(&candidate.contig),
            csv_cell(&candidate.representative),
            "reference_group".to_string(),
            csv_cell(&json_string_array(&candidate.hypothesis_members)),
            "not_resolved".to_string(),
            candidate.discovery_reads.to_string(),
            candidate.index_member_count.to_string(),
            candidate.target_exposure_bases.to_string(),
            candidate.decision.to_string(),
            csv_cell(&json_str_array(&candidate.decision_reasons)),
            candidate.q_value.to_string(),
            (candidate.ln_q_value / std::f64::consts::LN_10).to_string(),
            probability_status(candidate.q_underflow).to_string(),
            crate::MODEL_ADJUSTED_P_MAX.to_string(),
            candidate.p_value.to_string(),
            (candidate.ln_p_value / std::f64::consts::LN_10).to_string(),
            probability_status(candidate.p_underflow).to_string(),
            candidate.p_resolution_floor.to_string(),
            candidate.p_floor_flag.to_string(),
            candidate.reads.to_string(),
            candidate.covered_bases.to_string(),
            candidate.contig_len.to_string(),
            candidate.covered_frac.to_string(),
            crate::COVERAGE_MIN.to_string(),
            candidate.distinct_windows.to_string(),
            crate::DISTRIBUTED_WINDOW_BINS.to_string(),
            crate::MIN_DISTRIBUTED_WINDOWS.to_string(),
            candidate.depth_rpm.to_string(),
            candidate.expected_hits.to_string(),
            candidate
                .depth_fold
                .map_or_else(String::new, |value| value.to_string()),
            csv_cell(&candidate.stratum),
            candidate.stratum_decoy_count.to_string(),
            candidate.background_status.to_string(),
            candidate.background_scope.to_string(),
            candidate.background_reads.to_string(),
            candidate.background_reference_bases.to_string(),
            candidate.background_cross_stratum_reads.to_string(),
            candidate.integration_evidence.to_string(),
            candidate.split_events.to_string(),
            candidate.sites.len().to_string(),
            candidate.discordant.to_string(),
            candidate.plus_strand.to_string(),
            candidate.minus_strand.to_string(),
        ];
        writeln!(out, "{}", values.join(","))
            .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
    }
    out.flush()
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))
}

fn write_html(path: &Path, csv_path: &Path, input: &HumanReportInput<'_>) -> Result<(), String> {
    let mut body = String::with_capacity(48_000 + input.candidates.len() * 5_000);
    let pass_count = input
        .candidates
        .iter()
        .filter(|candidate| candidate.decision == "PASS")
        .count();
    let below_count = input
        .candidates
        .iter()
        .filter(|candidate| candidate.decision == "BELOW_THRESHOLD")
        .count();
    let not_significant_count = input
        .candidates
        .iter()
        .filter(|candidate| candidate.decision == "NOT_SIGNIFICANT")
        .count();
    let integration_count = input
        .candidates
        .iter()
        .filter(|candidate| candidate.integration_evidence != "NONE")
        .count();
    let qc_issues = qc_issue_codes(input);
    let csv_name = csv_path
        .file_name()
        .and_then(|name| name.to_str())
        .unwrap_or("candidates.csv");

    body.push_str("<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">");
    body.push_str("<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">");
    body.push_str("<meta name=\"color-scheme\" content=\"light\">");
    let _ = write!(
        body,
        "<title>{} · viroflash evidence report</title>",
        html_escape(input.sample)
    );
    body.push_str("<style>");
    body.push_str(STYLES);
    body.push_str("</style></head><body>");
    body.push_str(
        "<!--\nTHESIS: A viral candidate report is an evidence ledger, not a diagnostic scorecard.\nOWN-WORLD: Paper-white assay field, slate ink, hairline rules, and narrow mint/rose/violet evidence bands.\nSTORY: Orient to scope, triage the review queue, inspect candidate evidence, then audit method and provenance.\nFIRST VIEWPORT: Sample identity and interpretation boundary lead; a compact review queue and evidence matrix follow without hero metrics.\nFORM: Evidence Ledger, grounded direction 6, seed 40c95b97.\nFINISH: unreviewed and undocumented is unfinished; this build ends with the finish review, the verdict, DESIGN.md, and every shipping raster carrying its provenance\n-->",
    );

    body.push_str("<a class=\"skip-link\" href=\"#candidates\">Skip to candidates</a>");
    body.push_str("<header class=\"report-head\"><div class=\"head-main\">");
    body.push_str("<div class=\"wordmark\">viroflash <span>evidence report</span></div>");
    let _ = write!(body, "<h1>{}</h1>", html_escape(input.sample));
    body.push_str("<p class=\"scope\">Viral candidate evidence from the validation sample. This report does not compute a sample-level conclusion.</p>");
    body.push_str("</div><div class=\"head-actions\">");
    let _ = write!(
        body,
        "<a class=\"button\" href=\"{}\" download>Download core CSV</a>",
        html_escape(csv_name)
    );
    body.push_str("<button class=\"button secondary\" id=\"print-report\" type=\"button\">Print report</button>");
    body.push_str("</div></header>");

    body.push_str("<main><section class=\"orientation\" aria-labelledby=\"orientation-title\">");
    body.push_str(
        "<div class=\"interpretation\"><h2 id=\"orientation-title\">Interpretation boundary</h2>",
    );
    body.push_str("<p><strong>PASS</strong> means that a candidate crossed the current exploratory statistical, breadth, and distribution gates. It is not a clinical positive call. <strong>No candidates reported</strong> is not equivalent to NOT_DETECTED.</p>");
    body.push_str("<div class=\"semantic-tags\"><span>Candidate scope</span><span>QC not evaluated</span><span>FDR not validated</span><span>Member attribution unresolved</span></div></div>");
    body.push_str(
        "<div class=\"review-queue\" aria-label=\"Review queue summary\"><h2>Review queue</h2><dl>",
    );
    summary_row(&mut body, "PASS", pass_count, "pass");
    summary_row(&mut body, "Below threshold", below_count, "below");
    summary_row(
        &mut body,
        "Not significant",
        not_significant_count,
        "not-significant",
    );
    summary_row(
        &mut body,
        "Integration evidence",
        integration_count,
        "integration",
    );
    body.push_str("</dl></div></section>");

    body.push_str("<section class=\"run-strip\" aria-label=\"Run overview\">");
    metric(&mut body, "Input pairs", &format_integer(input.input_pairs));
    metric(
        &mut body,
        "Sampled pairs",
        &format_integer(input.sampling.selected_pairs),
    );
    metric(
        &mut body,
        "Validation read ends",
        &format_integer(input.sampling.validation_read_sides),
    );
    metric(
        &mut body,
        "Discovery family",
        &format_integer(input.test_family_size as u64),
    );
    metric(
        &mut body,
        "Reported candidates",
        &format_integer(input.candidates.len() as u64),
    );
    metric(
        &mut body,
        "Sampling inclusion",
        &format!("{:.1}%", input.sampling.inclusion_probability * 100.0),
    );
    body.push_str("</section>");

    body.push_str("<section class=\"candidate-section\" id=\"candidates\" aria-labelledby=\"candidate-title\">");
    body.push_str("<div class=\"section-head\"><div><h2 id=\"candidate-title\">Candidate evidence</h2><p>Sorted by adjusted p-value. Search includes candidate, representative, and unresolved members.</p></div><p class=\"result-count\" id=\"result-count\" aria-live=\"polite\"></p></div>");
    body.push_str("<div class=\"controls\" aria-label=\"Candidate controls\">");
    body.push_str("<label>Search<input id=\"candidate-search\" type=\"search\" placeholder=\"Candidate or member\" autocomplete=\"off\"></label>");
    body.push_str("<label>Decision<select id=\"decision-filter\"><option value=\"all\">All decisions</option><option value=\"PASS\">PASS</option><option value=\"BELOW_THRESHOLD\">Below threshold</option><option value=\"NOT_SIGNIFICANT\">Not significant</option></select></label>");
    body.push_str("<label>Sort<select id=\"candidate-sort\"><option value=\"q\">Adjusted p-value</option><option value=\"reads\">Validation reads</option><option value=\"breadth\">Coverage breadth</option><option value=\"windows\">Distributed windows</option></select></label>");
    body.push_str("<button class=\"text-button\" id=\"toggle-details\" type=\"button\">Expand all details</button></div>");

    body.push_str("<div class=\"candidate-list\" id=\"candidate-list\">");
    if input.candidates.is_empty() {
        body.push_str("<div class=\"empty-state\"><h3>No candidates reported</h3><p>No discovery hypothesis had reportable validation evidence. This is not a sample-level negative conclusion. Review sampling, QC observations, reference scope, and independent controls before interpretation.</p></div>");
    } else {
        for (index, candidate) in input.candidates.iter().enumerate() {
            candidate_html(&mut body, candidate, index);
        }
        body.push_str("<div class=\"no-filter-results\" id=\"no-filter-results\" hidden><h3>No candidates match these filters</h3><p>Clear the search or select another decision state.</p></div>");
    }
    body.push_str("</div></section>");

    method_html(&mut body, input, &qc_issues);
    let _ = write!(
        body,
        "</main><footer><p>Generated offline by viroflash. Machine contract: <code>{RESULT_SCHEMA}</code>. Human CSV contract: <code>{CSV_SCHEMA}</code>.</p><p>Review the JSON report for complete structured metadata and integration-site details.</p></footer>"
    );
    body.push_str("<script>");
    body.push_str(SCRIPT);
    body.push_str("</script></body></html>");

    let file = File::create(path).map_err(|e| format!("Cannot create {}: {e}", path.display()))?;
    let mut out = BufWriter::new(file);
    out.write_all(body.as_bytes())
        .and_then(|_| out.flush())
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))
}

fn candidate_html(out: &mut String, candidate: &Candidate, index: usize) {
    let searchable = format!(
        "{} {} {}",
        candidate.contig,
        candidate.representative,
        candidate.hypothesis_members.join(" ")
    )
    .to_lowercase();
    let decision_class = match candidate.decision {
        "PASS" => "pass",
        "BELOW_THRESHOLD" => "below",
        _ => "not-significant",
    };
    let log_q = candidate.ln_q_value / std::f64::consts::LN_10;
    let _ = write!(
        out,
        "<article class=\"candidate\" data-decision=\"{}\" data-search=\"{}\" data-q=\"{}\" data-reads=\"{}\" data-breadth=\"{}\" data-windows=\"{}\" data-order=\"{}\">",
        html_escape(candidate.decision),
        html_escape(&searchable),
        log_q,
        candidate.reads,
        candidate.covered_frac,
        candidate.distinct_windows,
        index
    );
    out.push_str("<div class=\"candidate-summary\"><div class=\"candidate-identity\">");
    let _ = write!(
        out,
        "<span class=\"status {}\">{}</span><h3>{}</h3><p>Representative: <strong>{}</strong></p>",
        decision_class,
        human_decision(candidate.decision),
        html_escape(&candidate.contig),
        html_escape(&candidate.representative)
    );
    out.push_str("</div><div class=\"primary-evidence\">");
    evidence_value(
        out,
        "Adjusted p",
        &probability_label(candidate.q_value, log_q, candidate.q_underflow),
        "BH-adjusted over the fixed discovery family; not a validated FDR.",
    );
    evidence_value(
        out,
        "Validation read ends",
        &format_integer(candidate.reads),
        "Each directly compatible validation read end is counted once per hypothesis.",
    );
    evidence_value(
        out,
        "Coverage breadth",
        &format!("{:.1}%", candidate.covered_frac * 100.0),
        "Observed representative breadth in the validation sample.",
    );
    evidence_value(
        out,
        "Distributed windows",
        &format!(
            "{} / {}",
            candidate.distinct_windows,
            crate::DISTRIBUTED_WINDOW_BINS
        ),
        "Count of fixed positional windows with direct evidence.",
    );
    out.push_str("</div></div>");

    out.push_str("<div class=\"evidence-rails\">");
    let breadth = (candidate.covered_frac * 100.0).clamp(0.0, 100.0);
    let _ = write!(
        out,
        "<div class=\"rail-block\"><div class=\"rail-label\"><span>Representative breadth</span><span>{:.1}% · gate ≥ {:.0}%</span></div><div class=\"breadth-rail\" role=\"img\" aria-label=\"Coverage breadth {:.1} percent; gate {:.0} percent\"><span style=\"width:{:.3}%\"></span><i style=\"left:{}%\"></i></div></div>",
        breadth,
        crate::COVERAGE_MIN * 100.0,
        breadth,
        crate::COVERAGE_MIN * 100.0,
        breadth,
        crate::COVERAGE_MIN * 100.0
    );
    out.push_str(
        "<div class=\"rail-block\"><div class=\"rail-label\"><span>Evidence distribution</span>",
    );
    let _ = write!(
        out,
        "<span>{} windows · gate ≥ {}</span></div><div class=\"window-rail\" role=\"img\" aria-label=\"Evidence observed in {} of {} positional windows; window locations are not retained in this summary\">",
        candidate.distinct_windows,
        crate::MIN_DISTRIBUTED_WINDOWS,
        candidate.distinct_windows,
        crate::DISTRIBUTED_WINDOW_BINS
    );
    for window in 0..crate::DISTRIBUTED_WINDOW_BINS {
        out.push_str(if window < candidate.distinct_windows {
            "<span class=\"observed\"></span>"
        } else {
            "<span></span>"
        });
    }
    out.push_str("</div><small>Quantity indicator only; segment positions do not encode observed genomic locations.</small></div></div>");

    out.push_str("<details><summary>Inspect hypothesis, background, and integration evidence</summary><div class=\"detail-grid\">");
    out.push_str("<section><h4>Reference-group hypothesis</h4><p class=\"definition\">At least one listed member is supported under the current reference model. Individual member attribution is unresolved.</p><ul class=\"member-list\">");
    for member in &candidate.hypothesis_members {
        let _ = write!(out, "<li>{}</li>", html_escape(member));
    }
    out.push_str("</ul></section>");

    out.push_str("<section><h4>Statistical background</h4><dl class=\"fact-list\">");
    fact(
        out,
        "Raw p-value",
        &format_probability(
            candidate.p_value,
            candidate.ln_p_value,
            candidate.p_underflow,
        ),
    );
    fact(
        out,
        "Expected hits",
        &format_compact(candidate.expected_hits),
    );
    fact(
        out,
        "Observed / expected",
        &candidate
            .depth_fold
            .map_or_else(|| "Not defined".to_string(), format_compact),
    );
    fact(out, "Background status", candidate.background_status);
    fact(out, "Background scope", candidate.background_scope);
    fact(
        out,
        "Background reads",
        &format_integer(candidate.background_reads),
    );
    out.push_str("</dl></section>");

    out.push_str("<section><h4>Evidence audit</h4><dl class=\"fact-list\">");
    fact(
        out,
        "Covered bases",
        &format!(
            "{} / {}",
            format_integer(candidate.covered_bases),
            format_integer(candidate.contig_len)
        ),
    );
    fact(out, "End RPM", &format_compact(candidate.depth_rpm));
    fact(
        out,
        "Discovery reads",
        &format_integer(candidate.discovery_reads),
    );
    fact(
        out,
        "Strand support",
        &format!("+{} / −{}", candidate.plus_strand, candidate.minus_strand),
    );
    fact(
        out,
        "Discordant pairs",
        &format_integer(candidate.discordant),
    );
    out.push_str("</dl></section>");

    out.push_str("<section><h4>Integration evidence</h4>");
    let _ = write!(
        out,
        "<p><span class=\"integration-state\">{}</span></p><p class=\"definition\">Separate from general viral candidate detection. {} split events; {} supported sites.</p>",
        human_token(candidate.integration_evidence),
        candidate.split_events,
        candidate.sites.len()
    );
    if !candidate.sites.is_empty() {
        out.push_str("<div class=\"table-wrap\"><table><thead><tr><th>Virus position</th><th>Host</th><th>Host position</th><th>Direction</th><th>Support</th></tr></thead><tbody>");
        for site in &candidate.sites {
            let _ = write!(
                out,
                "<tr><td>{}</td><td>{}</td><td>{}</td><td>{}</td><td>{}</td></tr>",
                site.pos,
                html_escape(&site.host_contig),
                site.host_pos,
                html_escape(&site.direction),
                site.support
            );
        }
        out.push_str("</tbody></table></div>");
    }
    out.push_str(
        "</section></div><div class=\"decision-reasons\"><strong>Decision basis</strong><ul>",
    );
    for reason in &candidate.decision_reasons {
        let _ = write!(out, "<li>{}</li>", human_token(reason));
    }
    out.push_str("</ul></div></details></article>");
}

fn method_html(out: &mut String, input: &HumanReportInput<'_>, qc_issues: &[&str]) {
    out.push_str("<section class=\"audit-section\" id=\"audit\" aria-labelledby=\"audit-title\"><div class=\"section-head\"><div><h2 id=\"audit-title\">Run audit and method</h2><p>Observed metadata and prespecified gates. QC status is not evaluated.</p></div><span class=\"status below\">NOT EVALUATED</span></div>");
    if qc_issues.is_empty() {
        out.push_str("<p class=\"audit-note\"><strong>No flagged internal observations.</strong> This does not establish QC pass; validated pass/fail thresholds are not available.</p>");
    } else {
        out.push_str("<div class=\"audit-note warning\"><strong>Internal observations requiring review</strong><ul>");
        for issue in qc_issues {
            let _ = write!(out, "<li>{}</li>", human_token(issue));
        }
        out.push_str("</ul></div>");
    }
    out.push_str("<div class=\"audit-columns\"><section><h3>Sampling</h3><dl class=\"fact-list\">");
    fact(out, "Method", input.sampling.method);
    fact(
        out,
        "Capacity pairs",
        &format_integer(input.sampling.capacity_pairs as u64),
    );
    fact(
        out,
        "Selected pairs",
        &format_integer(input.sampling.selected_pairs),
    );
    fact(
        out,
        "Discovery pairs",
        &format_integer(input.sampling.discovery_pairs),
    );
    fact(
        out,
        "Validation pairs",
        &format_integer(input.sampling.validation_pairs),
    );
    fact(
        out,
        "Inclusion probability",
        &format!("{:.6}", input.sampling.inclusion_probability),
    );
    out.push_str("</dl></section><section><h3>Audit observations</h3><dl class=\"fact-list\">");
    fact(out, "Map errors", &format_integer(input.map_errors));
    fact(
        out,
        "Audit overflows",
        &format_integer(input.sampling.audit_overflows),
    );
    fact(
        out,
        "Target unassigned",
        &format_integer(input.sampling.target_validation_unassigned),
    );
    fact(
        out,
        "Decoy unassigned",
        &format_integer(input.sampling.decoy_validation_unassigned),
    );
    fact(
        out,
        "Prescreened sample pairs",
        &format_integer(input.prescreen_pairs),
    );
    out.push_str("</dl></section><section><h3>Candidate gates</h3><dl class=\"fact-list\">");
    fact(
        out,
        "Adjusted p-value",
        &format!("< {:.2}", crate::MODEL_ADJUSTED_P_MAX),
    );
    fact(
        out,
        "Coverage breadth",
        &format!("≥ {:.0}%", crate::COVERAGE_MIN * 100.0),
    );
    fact(
        out,
        "Distributed windows",
        &format!(
            "≥ {} of {}",
            crate::MIN_DISTRIBUTED_WINDOWS,
            crate::DISTRIBUTED_WINDOW_BINS
        ),
    );
    fact(out, "Primary test", "Exact conditional two-Poisson rates");
    fact(
        out,
        "Adjustment",
        "Benjamini–Hochberg fixed discovery family",
    );
    out.push_str("</dl></section><section><h3>Provenance</h3><dl class=\"fact-list\">");
    fact(out, "Result schema", RESULT_SCHEMA);
    fact(out, "Index source", input.index_source);
    fact(out, "Index format", &input.index_format_version.to_string());
    fact(out, "Manifest BLAKE3", input.manifest_blake3);
    fact(out, "Threads", &input.threads.to_string());
    fact(out, "k-mer length", &input.k.to_string());
    out.push_str("</dl></section></div>");
    out.push_str("<details class=\"limitations\"><summary>Interpretation limitations</summary><ul><li>Synthetic decoys are unvalidated null-model stress references and do not substitute for extraction blanks, batch-matched controls, or laboratory contamination controls.</li><li>The adjusted p-value does not carry a validated classical FDR guarantee or clinical confidence calibration.</li><li>Sampling capacity, alignment gates, breadth and window thresholds require validation for the intended matrix, viral classes, limit of detection, and near-neighbor interference.</li><li>An unresolved OR hypothesis supports at least one member under the current reference model; it does not identify an individual accession.</li></ul></details></section>");
}

fn summary_row(out: &mut String, label: &str, value: usize, class_name: &str) {
    let _ = write!(
        out,
        "<div><dt><span class=\"queue-mark {}\"></span>{}</dt><dd>{}</dd></div>",
        class_name,
        html_escape(label),
        value
    );
}

fn metric(out: &mut String, label: &str, value: &str) {
    let _ = write!(
        out,
        "<div><span>{}</span><strong>{}</strong></div>",
        html_escape(label),
        html_escape(value)
    );
}

fn evidence_value(out: &mut String, label: &str, value: &str, help: &str) {
    let _ = write!(
        out,
        "<div><span>{}</span><strong>{}</strong><small>{}</small></div>",
        html_escape(label),
        html_escape(value),
        html_escape(help)
    );
}

fn fact(out: &mut String, label: &str, value: &str) {
    let _ = write!(
        out,
        "<div><dt>{}</dt><dd>{}</dd></div>",
        html_escape(label),
        html_escape(value)
    );
}

fn qc_issue_codes(input: &HumanReportInput<'_>) -> Vec<&'static str> {
    let mut issues = Vec::new();
    if input.map_errors > 0 || input.sampling.audit_map_error_reads > 0 {
        issues.push("mapping_errors_observed");
    }
    if input.sampling.audit_overflows > 0 {
        issues.push("ambiguity_audit_overflow_observed");
    }
    if input.sampling.target_validation_unassigned > 0 {
        issues.push("target_validation_unassigned_observed");
    }
    if input.sampling.decoy_validation_unassigned > 0 {
        issues.push("decoy_validation_unassigned_observed");
    }
    issues
}

fn probability_status(underflow: bool) -> &'static str {
    if underflow {
        "NUMERIC_UNDERFLOW_LOG_AVAILABLE"
    } else {
        "FINITE"
    }
}

fn probability_label(value: f64, log10_value: f64, underflow: bool) -> String {
    if underflow {
        format!("underflow · log₁₀ {log10_value:.2}")
    } else if value == 0.0 {
        "0".to_string()
    } else if value < 0.001 {
        format!("{value:.2e}")
    } else {
        format!("{value:.4}")
    }
}

fn format_probability(value: f64, ln_value: f64, underflow: bool) -> String {
    probability_label(value, ln_value / std::f64::consts::LN_10, underflow)
}

fn format_compact(value: f64) -> String {
    if value == 0.0 {
        "0".to_string()
    } else if value.abs() >= 100_000.0 || value.abs() < 0.001 {
        format!("{value:.3e}")
    } else if value.fract().abs() < f64::EPSILON {
        format!("{value:.0}")
    } else {
        format!("{value:.3}")
    }
}

fn format_integer(value: u64) -> String {
    let digits = value.to_string();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    for (index, byte) in digits.bytes().enumerate() {
        if index > 0 && (digits.len() - index).is_multiple_of(3) {
            out.push(',');
        }
        out.push(byte as char);
    }
    out
}

fn human_decision(value: &str) -> &'static str {
    match value {
        "PASS" => "PASS",
        "BELOW_THRESHOLD" => "BELOW THRESHOLD",
        _ => "NOT SIGNIFICANT",
    }
}

fn human_token(value: &str) -> String {
    value.replace('_', " ").to_lowercase()
}

fn html_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\'' => out.push_str("&#39;"),
            _ => out.push(character),
        }
    }
    out
}

fn csv_cell(value: &str) -> String {
    let protected;
    let value = if value
        .chars()
        .next()
        .is_some_and(|character| matches!(character, '=' | '+' | '-' | '@'))
    {
        protected = format!("'{value}");
        &protected
    } else {
        value
    };
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn json_string_array(values: &[String]) -> String {
    let values = values
        .iter()
        .map(|value| format!("\"{}\"", crate::report::json_escape(value)))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{values}]")
}

fn json_str_array(values: &[&str]) -> String {
    let values = values
        .iter()
        .map(|value| format!("\"{}\"", crate::report::json_escape(value)))
        .collect::<Vec<_>>()
        .join(",");
    format!("[{values}]")
}

const STYLES: &str = r#"
:root{--paper:#f7f8f5;--surface:#fff;--ink:#182029;--muted:#59646f;--rule:#d8ddd9;--rule-strong:#aeb8b3;--mint:#1c775f;--mint-soft:#dff2e9;--rose:#a63d59;--rose-soft:#f8e5ea;--violet:#6550a5;--violet-soft:#ece8f8;--amber:#8a5b12;--amber-soft:#f7ecd5;--focus:#3659c9;--shadow:0 14px 34px rgba(31,42,38,.08)}
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:var(--paper);color:var(--ink);font-family:ui-sans-serif,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;font-size:16px;line-height:1.55}::selection{background:var(--violet-soft);color:var(--ink)}button,input,select{font:inherit}button,a,input,select,summary{outline-offset:3px}:focus-visible{outline:3px solid var(--focus)}a{color:var(--violet);text-underline-offset:.18em}code{font-family:"SFMono-Regular",Consolas,monospace;font-size:.9em}.skip-link{position:fixed;left:1rem;top:1rem;z-index:20;transform:translateY(-180%);background:var(--ink);color:#fff;padding:.65rem 1rem}.skip-link:focus{transform:none}.report-head{max-width:1480px;margin:0 auto;padding:2.5rem 3rem 2rem;display:flex;justify-content:space-between;gap:2rem;align-items:flex-end;border-bottom:1px solid var(--rule-strong)}.wordmark{font-weight:760;letter-spacing:-.02em}.wordmark span{font-weight:450;color:var(--muted)}h1{font-size:clamp(2.2rem,4vw,3.8rem);line-height:.98;letter-spacing:-.035em;margin:2.2rem 0 .85rem;max-width:18ch;overflow-wrap:anywhere}.scope{max-width:67ch;color:var(--muted);font-size:1.08rem;margin:0}.head-actions{display:flex;gap:.7rem;flex-wrap:wrap}.button{display:inline-flex;align-items:center;justify-content:center;border:1px solid var(--ink);border-radius:12px;background:var(--ink);color:#fff;text-decoration:none;padding:.72rem 1rem;font-weight:700;cursor:pointer;transition:transform .18s ease-out,box-shadow .18s ease-out}.button:hover{transform:translateY(-1px);box-shadow:0 7px 18px rgba(24,32,41,.14)}.button.secondary{background:transparent;color:var(--ink);border-color:var(--rule-strong)}main{max-width:1480px;margin:0 auto;padding:0 3rem 5rem}.orientation{display:grid;grid-template-columns:minmax(0,1.65fr) minmax(280px,.75fr);gap:0;border-bottom:1px solid var(--rule-strong)}.interpretation,.review-queue{padding:2.3rem 0}.interpretation{padding-right:4rem}.review-queue{border-left:1px solid var(--rule);padding-left:2.2rem}h2{font-size:1.45rem;line-height:1.2;letter-spacing:-.02em;margin:0 0 .7rem}.interpretation p,.section-head p{max-width:72ch;color:var(--muted);margin:.5rem 0}.semantic-tags{display:flex;gap:.5rem;flex-wrap:wrap;margin-top:1.4rem}.semantic-tags span,.integration-state{border:1px solid var(--rule-strong);border-radius:999px;padding:.28rem .62rem;font-size:.78rem;font-weight:720;letter-spacing:.03em;text-transform:uppercase}.review-queue dl{margin:1.2rem 0 0}.review-queue dl>div{display:flex;justify-content:space-between;align-items:center;padding:.48rem 0;border-bottom:1px solid var(--rule)}.review-queue dt{display:flex;align-items:center;gap:.65rem;color:var(--muted)}.review-queue dd{font-variant-numeric:tabular-nums;font-weight:780;font-size:1.2rem}.queue-mark{width:.75rem;height:.75rem;border-radius:2px;background:var(--rule-strong)}.queue-mark.pass{background:var(--mint)}.queue-mark.below{background:var(--amber)}.queue-mark.not-significant{background:var(--rose)}.queue-mark.integration{background:var(--violet)}.run-strip{display:grid;grid-template-columns:repeat(6,1fr);border-bottom:1px solid var(--rule-strong)}.run-strip>div{padding:1.1rem 1.1rem 1.15rem 0;border-right:1px solid var(--rule);margin-right:1.1rem}.run-strip>div:last-child{border-right:0}.run-strip span{display:block;color:var(--muted);font-size:.78rem}.run-strip strong{display:block;margin-top:.2rem;font-size:1.08rem;font-variant-numeric:tabular-nums}.candidate-section,.audit-section{padding-top:4.5rem}.section-head{display:flex;align-items:flex-end;justify-content:space-between;gap:2rem;margin-bottom:1.35rem}.result-count{font-variant-numeric:tabular-nums;white-space:nowrap}.controls{display:grid;grid-template-columns:minmax(220px,1.3fr) minmax(170px,.7fr) minmax(180px,.7fr) auto;gap:.75rem;align-items:end;padding:1rem 0 1.2rem;border-top:1px solid var(--rule-strong);border-bottom:1px solid var(--rule-strong)}.controls label{font-size:.75rem;font-weight:720;color:var(--muted);letter-spacing:.02em}.controls input,.controls select{display:block;width:100%;margin-top:.35rem;border:1px solid var(--rule-strong);border-radius:10px;background:var(--surface);color:var(--ink);padding:.65rem .75rem}.text-button{border:0;background:transparent;color:var(--violet);font-weight:720;padding:.68rem;cursor:pointer;text-decoration:underline;text-underline-offset:.18em}.candidate{background:var(--surface);border-bottom:1px solid var(--rule-strong);position:relative}.candidate::before{content:"";position:absolute;inset:0 0 auto 0;height:3px;background:var(--rule)}.candidate[data-decision="PASS"]::before{background:linear-gradient(var(--mint),var(--violet))}.candidate[data-decision="BELOW_THRESHOLD"]::before{background:linear-gradient(var(--amber),var(--rose))}.candidate-summary{display:grid;grid-template-columns:minmax(250px,.8fr) minmax(520px,1.5fr);gap:2.4rem;padding:1.7rem 1.8rem 1.3rem 2rem}.candidate-identity h3{font-size:1.35rem;line-height:1.25;margin:.75rem 0 .25rem;overflow-wrap:anywhere}.candidate-identity p{color:var(--muted);font-size:.85rem;margin:0;overflow-wrap:anywhere}.status{display:inline-flex;border-radius:999px;padding:.28rem .6rem;font-size:.72rem;font-weight:800;letter-spacing:.05em}.status.pass{color:#105c48;background:var(--mint-soft)}.status.below{color:#714708;background:var(--amber-soft)}.status.not-significant{color:#873049;background:var(--rose-soft)}.primary-evidence{display:grid;grid-template-columns:repeat(4,1fr);gap:1.5rem}.primary-evidence span,.primary-evidence strong,.primary-evidence small{display:block}.primary-evidence span{font-size:.75rem;color:var(--muted)}.primary-evidence strong{font-size:1.12rem;margin:.18rem 0;font-variant-numeric:tabular-nums}.primary-evidence small{font-size:.7rem;color:var(--muted);line-height:1.35}.evidence-rails{display:grid;grid-template-columns:1fr 1fr;gap:2.5rem;padding:0 1.8rem 1.55rem 2rem}.rail-label{display:flex;justify-content:space-between;gap:1rem;font-size:.75rem;color:var(--muted);margin-bottom:.45rem}.breadth-rail{height:10px;background:#e4e8e5;position:relative;overflow:visible}.breadth-rail span{display:block;height:100%;background:var(--violet)}.breadth-rail i{position:absolute;top:-4px;height:18px;width:2px;background:var(--ink)}.window-rail{display:grid;grid-template-columns:repeat(10,1fr);gap:4px}.window-rail span{height:10px;background:#e4e8e5}.window-rail span.observed{background:var(--mint)}.rail-block small{display:block;color:var(--muted);font-size:.67rem;margin-top:.38rem}.candidate details{border-top:1px solid var(--rule);margin-left:2rem}.candidate summary,.limitations summary{cursor:pointer;padding:1rem 1.8rem 1rem 0;font-weight:720;color:var(--violet)}.candidate summary::marker,.limitations summary::marker{color:var(--violet)}.detail-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:2.5rem;padding:1.2rem 1.8rem 1.8rem 0}.detail-grid section{min-width:0}.detail-grid h4{font-size:.9rem;margin:0 0 .55rem}.definition{color:var(--muted);font-size:.82rem;max-width:70ch}.member-list{max-height:11rem;overflow:auto;margin:.75rem 0 0;padding:0;list-style:none;border-top:1px solid var(--rule)}.member-list li{padding:.42rem 0;border-bottom:1px solid var(--rule);overflow-wrap:anywhere;font-size:.83rem}.fact-list{margin:0}.fact-list>div{display:grid;grid-template-columns:minmax(120px,.85fr) minmax(0,1.15fr);gap:1rem;border-bottom:1px solid var(--rule);padding:.38rem 0}.fact-list dt{color:var(--muted);font-size:.78rem}.fact-list dd{margin:0;text-align:right;overflow-wrap:anywhere;font-size:.8rem;font-variant-numeric:tabular-nums}.decision-reasons{margin:0 1.8rem 1.8rem 0;padding:1rem 1.1rem;background:var(--violet-soft)}.decision-reasons ul{display:flex;flex-wrap:wrap;gap:.4rem 1.5rem;margin:.45rem 0 0;padding-left:1.1rem;font-size:.8rem}.table-wrap{overflow:auto}.table-wrap table{width:100%;border-collapse:collapse;font-size:.75rem}.table-wrap th,.table-wrap td{text-align:left;padding:.4rem;border-bottom:1px solid var(--rule);font-variant-numeric:tabular-nums}.empty-state,.no-filter-results{padding:3rem 0;border-bottom:1px solid var(--rule-strong);max-width:72ch}.empty-state h3,.no-filter-results h3{margin:0}.empty-state p,.no-filter-results p{color:var(--muted)}.audit-note{padding:1rem 0;border-top:1px solid var(--rule-strong);border-bottom:1px solid var(--rule);color:var(--muted)}.audit-note.warning{color:var(--ink);background:var(--amber-soft);padding:1rem}.audit-columns{display:grid;grid-template-columns:repeat(4,1fr);gap:2.5rem;padding:2rem 0}.audit-columns h3{font-size:.95rem;margin:0 0 .8rem}.limitations{border-top:1px solid var(--rule-strong);border-bottom:1px solid var(--rule-strong)}.limitations ul{max-width:88ch;padding:0 0 1.5rem 1.2rem;color:var(--muted)}footer{max-width:1480px;margin:0 auto;padding:2rem 3rem 4rem;border-top:1px solid var(--rule-strong);display:flex;justify-content:space-between;gap:2rem;color:var(--muted);font-size:.78rem}footer p{margin:0;max-width:70ch}
@media(max-width:1000px){.report-head,.orientation{display:block}.head-actions{margin-top:1.5rem}.review-queue{border-left:0;border-top:1px solid var(--rule);padding-left:0}.interpretation{padding-right:0}.run-strip{grid-template-columns:repeat(3,1fr)}.candidate-summary{grid-template-columns:1fr}.primary-evidence{grid-template-columns:repeat(2,1fr)}.controls{grid-template-columns:1fr 1fr}.audit-columns{grid-template-columns:repeat(2,1fr)}}
@media(max-width:640px){body{font-size:15px}.report-head,main,footer{padding-left:1rem;padding-right:1rem}.report-head{padding-top:1.5rem}.head-actions .button{flex:1}.orientation{border-bottom:0}.run-strip{grid-template-columns:repeat(2,1fr);border-top:1px solid var(--rule-strong)}.run-strip>div{margin-right:.6rem}.candidate-section,.audit-section{padding-top:3rem}.section-head{display:block}.controls{grid-template-columns:1fr}.candidate-summary{padding:1.4rem 1rem 1.2rem 1.3rem}.primary-evidence{gap:1rem}.evidence-rails{grid-template-columns:1fr;gap:1.2rem;padding:0 1rem 1.3rem 1.3rem}.candidate details{margin-left:1.3rem}.detail-grid{grid-template-columns:1fr;padding-right:1rem}.audit-columns{grid-template-columns:1fr;gap:2rem}.fact-list>div{grid-template-columns:1fr 1fr}footer{display:block}footer p+p{margin-top:.7rem}}
@media(prefers-reduced-motion:reduce){html{scroll-behavior:auto}.button{transition:none}}
@media print{body{background:#fff;font-size:10pt}.skip-link,.head-actions,.controls{display:none!important}.report-head,main,footer{max-width:none;padding-left:0;padding-right:0}.report-head{padding-top:0}.orientation{grid-template-columns:1.5fr 1fr}.candidate-section,.audit-section{padding-top:2rem}.candidate{break-inside:avoid}.candidate details:not([open])>:not(summary){display:block}.candidate summary{display:none}.detail-grid{grid-template-columns:1fr 1fr}.status{border:1px solid currentColor}.audit-columns{grid-template-columns:1fr 1fr}.run-strip{grid-template-columns:repeat(3,1fr)}footer{margin-top:2rem}}
"#;

const SCRIPT: &str = r#"
(()=>{const list=document.querySelector('#candidate-list');if(!list)return;const rows=[...list.querySelectorAll('.candidate')];const search=document.querySelector('#candidate-search');const decision=document.querySelector('#decision-filter');const sort=document.querySelector('#candidate-sort');const count=document.querySelector('#result-count');const empty=document.querySelector('#no-filter-results');const detailsButton=document.querySelector('#toggle-details');const apply=()=>{const term=(search?.value||'').trim().toLowerCase();const state=decision?.value||'all';let visible=0;rows.forEach(row=>{const show=(!term||row.dataset.search.includes(term))&&(state==='all'||row.dataset.decision===state);row.hidden=!show;if(show)visible++});const mode=sort?.value||'q';rows.sort((a,b)=>{if(mode==='q')return Number(a.dataset.q)-Number(b.dataset.q);if(mode==='reads')return Number(b.dataset.reads)-Number(a.dataset.reads);if(mode==='breadth')return Number(b.dataset.breadth)-Number(a.dataset.breadth);if(mode==='windows')return Number(b.dataset.windows)-Number(a.dataset.windows);return Number(a.dataset.order)-Number(b.dataset.order)}).forEach(row=>list.appendChild(row));if(empty){empty.hidden=visible!==0;list.appendChild(empty)}if(count)count.textContent=`${visible} of ${rows.length} candidates shown`};[search,decision,sort].forEach(control=>control?.addEventListener('input',apply));let expanded=false;detailsButton?.addEventListener('click',()=>{expanded=!expanded;rows.forEach(row=>{const details=row.querySelector('details');if(details)details.open=expanded});detailsButton.textContent=expanded?'Collapse all details':'Expand all details'});document.querySelector('#print-report')?.addEventListener('click',()=>window.print());apply()})();
"#;

#[cfg(test)]
mod tests {
    use super::*;

    fn sampling_report() -> SamplingReport {
        SamplingReport {
            method: "test",
            capacity_pairs: 10,
            selected_pairs: 8,
            selected_passed_pairs: 2,
            discovery_pairs: 4,
            validation_pairs: 4,
            validation_read_sides: 8,
            inclusion_probability: 0.8,
            duplicate_selected_qnames: 0,
            audit_reads: 4,
            audit_no_evidence: 4,
            audit_overflows: 0,
            audit_map_error_reads: 0,
            target_validation_unassigned: 0,
            decoy_validation_unassigned: 0,
        }
    }

    #[test]
    fn escapes_html_and_csv_cells() {
        assert_eq!(
            html_escape("A&B <x> \"q\""),
            "A&amp;B &lt;x&gt; &quot;q&quot;"
        );
        assert_eq!(csv_cell("plain"), "plain");
        assert_eq!(csv_cell("a,\"b\""), "\"a,\"\"b\"\"\"");
        assert_eq!(csv_cell("=SUM(A1:A2)"), "'=SUM(A1:A2)");
    }

    #[test]
    fn formats_integer_groups() {
        assert_eq!(format_integer(0), "0");
        assert_eq!(format_integer(999), "999");
        assert_eq!(format_integer(12_345_678), "12,345,678");
    }

    #[test]
    fn empty_report_preserves_non_diagnostic_semantics() {
        let prefix = std::env::temp_dir().join(format!(
            "viroflash_human_report_test_{}_{}",
            std::process::id(),
            std::thread::current().name().unwrap_or("unnamed")
        ));
        let sampling = sampling_report();
        let input = HumanReportInput {
            sample: "sample<script>",
            threads: 1,
            k: 21,
            input_pairs: 10,
            prescreen_pairs: 2,
            map_errors: 0,
            index_source: "loaded",
            index_format_version: 2,
            manifest_blake3: "abc",
            sampling: &sampling,
            test_family_size: 2,
            candidates: &[],
        };
        let (html_path, csv_path) = write_human_report(&prefix, &input).unwrap();
        let html = std::fs::read_to_string(&html_path).unwrap();
        let csv = std::fs::read_to_string(&csv_path).unwrap();
        assert!(html.contains("sample&lt;script&gt;"));
        assert!(!html.contains("<h1>sample<script></h1>"));
        assert!(html.contains("No candidates reported"));
        assert!(html.contains("not equivalent to NOT_DETECTED"));
        assert!(csv.starts_with(
            "csv_schema,result_schema,sample_id,sample_conclusion,qc_status,qc_issues"
        ));
        assert!(csv.contains("background_reference_bases"));
        assert_eq!(csv.lines().count(), 1);
        let _ = std::fs::remove_file(html_path);
        let _ = std::fs::remove_file(csv_path);
    }
}
