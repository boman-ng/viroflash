use std::path::{Path, PathBuf};

use crate::analysis_profile::AnalysisProfile;
use crate::competitive_alignment::CompetitiveAligner;
use crate::evidence::EvidenceAccumulator;
use crate::fastq_input::{census_fastq, compute_input_digest, FragmentReader};
use crate::performance_report::{stage_start, write_perf_json, PerformanceMonitor, StageTimes};
use crate::reference_index::load_index;
use crate::report::{build_evidence_report, write_report_csv, write_report_html, ReportInputs};
use crate::sampling_design::{derive_sampling_design, fragment_selection_key, include_fragment};

#[derive(Debug, Clone)]
pub struct RunOptions {
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
    pub index_dir: PathBuf,
    pub out_dir: PathBuf,
    pub threads: usize,
}

#[derive(Debug, Clone)]
pub struct RunSummary {
    pub input_fragments: u64,
    pub selected_fragments: u64,
    pub target_signal_rows: usize,
    pub report_dir: PathBuf,
}

pub fn run_pipeline(options: &RunOptions) -> Result<RunSummary, String> {
    validate_options(options)?;
    let sample_id = sample_id(&options.r1);
    let monitor = PerformanceMonitor::start();
    let staging = PathBuf::from(format!(
        "{}.part.{}",
        options.out_dir.display(),
        std::process::id()
    ));
    if staging.exists() {
        return Err(format!(
            "Report staging path already exists: {}",
            staging.display()
        ));
    }
    let mut stages = StageTimes::default();
    let mut counts = (0, 0, 0, 0);
    let result = (|| {
        std::fs::create_dir(&staging)
            .map_err(|error| format!("Cannot create {}: {error}", staging.display()))?;
        let profile = AnalysisProfile::FROZEN;
        let index = load_index(&options.index_dir)?;

        let started = stage_start();
        let census = census_fastq(&options.r1, options.r2.as_deref())?;
        stages.pass1_count = started.elapsed().as_millis() as u64;
        counts.0 = census.fragments;
        let design = derive_sampling_design(
            &census,
            index.target_groups.len(),
            profile.minimum_relevant_fraction,
            profile.familywise_miss_probability,
        )?;

        let started = stage_start();
        let aligner = CompetitiveAligner::open(&index.mmi_path)?;
        let mut fragments = FragmentReader::open(&options.r1, options.r2.as_deref())?;
        let mut accumulator = EvidenceAccumulator::new(&index.target_groups);
        let mut pass2_fragments = 0;
        while let Some(fragment) = fragments.next_fragment()? {
            pass2_fragments += 1;
            let key = fragment_selection_key(
                &index.profile_digest,
                &census.input_digest,
                &fragment.id,
                fragment.ordinal,
            );
            if !include_fragment(key, design.selection_probability) {
                continue;
            }
            counts.1 += 1;
            if !index
                .bloom
                .passes_target_kmer_gate(&fragment.r1, fragment.r2.as_deref())
            {
                continue;
            }
            counts.2 += 1;
            let evidence = aligner.align_fragment_competitively(&fragment, &index.contigs)?;
            accumulator.accumulate_group_evidence(evidence);
        }
        if pass2_fragments != census.fragments {
            return Err(format!("FASTQ changed between passes: pass 1 counted {} fragments, pass 2 counted {pass2_fragments}", census.fragments));
        }
        if compute_input_digest(&options.r1, options.r2.as_deref())? != census.input_digest {
            return Err("FASTQ bytes changed between pass 1 and pass 2".into());
        }
        counts.3 = accumulator.aligned_fragments;
        stages.pass2_sample_prescreen_align = started.elapsed().as_millis() as u64;

        let report = build_evidence_report(
            ReportInputs {
                sample_id: sample_id.clone(),
                census: &census,
                design,
                selected_fragments: counts.1,
                prescreen_passed_fragments: counts.2,
                profile,
                index_digest: index.index_digest,
            },
            accumulator,
        )?;
        let started = stage_start();
        write_report_csv(&staging.join("report.csv"), &report)?;
        write_report_html(&staging.join("report.html"), &report)?;
        stages.report_write = started.elapsed().as_millis() as u64;
        Ok((census, report))
    })();

    match result {
        Ok((census, report)) => {
            let target_signal_rows = report.target_signals.len();
            let perf = monitor.finish(
                "SUCCESS",
                sample_id,
                options.threads,
                stages,
                counts,
                Vec::new(),
            );
            write_perf_json(&staging.join("perf.json"), &perf)?;
            if let Err(error) = std::fs::rename(&staging, &options.out_dir) {
                let _ = std::fs::remove_dir_all(&staging);
                return Err(format!(
                    "Cannot finalize report directory {}: {error}",
                    options.out_dir.display()
                ));
            }
            Ok(RunSummary {
                input_fragments: census.fragments,
                selected_fragments: counts.1,
                target_signal_rows,
                report_dir: options.out_dir.clone(),
            })
        }
        Err(error) => {
            let _ = std::fs::remove_dir_all(&staging);
            std::fs::create_dir(&options.out_dir).map_err(|create_error| {
                format!("{error}; also cannot create error report directory: {create_error}")
            })?;
            let perf = monitor.finish(
                "ERROR",
                sample_id,
                options.threads,
                stages,
                counts,
                vec![error.clone()],
            );
            write_perf_json(&options.out_dir.join("perf.json"), &perf).map_err(|perf_error| {
                format!("{error}; also failed to write perf.json: {perf_error}")
            })?;
            Err(error)
        }
    }
}

fn validate_options(options: &RunOptions) -> Result<(), String> {
    if options.threads == 0 {
        return Err("--threads must be greater than zero".into());
    }
    if options.out_dir.exists() {
        return Err(format!(
            "Output directory already exists: {}",
            options.out_dir.display()
        ));
    }
    if options.r1.as_os_str().is_empty()
        || options.index_dir.as_os_str().is_empty()
        || options.out_dir.as_os_str().is_empty()
    {
        return Err("run requires --r1, --index, and --out".into());
    }
    Ok(())
}

fn sample_id(path: &Path) -> String {
    let mut name = path
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    for suffix in [".gz", ".fastq", ".fq", "_R1", "_1"] {
        if let Some(stripped) = name.strip_suffix(suffix) {
            name = stripped.to_string();
        }
    }
    name
}
