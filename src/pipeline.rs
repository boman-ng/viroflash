use std::path::{Path, PathBuf};

use crate::analysis_profile::AnalysisProfile;
use crate::competitive_alignment::align_fragments_bounded;
use crate::evidence::EvidenceAccumulator;
use crate::fastq_input::{census_fastq, FragmentReader};
use crate::kmer_gate::GateEvaluation;
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
    let mut unevaluable_fragments = 0;
    let result = (|| {
        std::fs::create_dir(&staging)
            .map_err(|error| format!("Cannot create {}: {error}", staging.display()))?;
        let profile = AnalysisProfile::FROZEN;
        let index = load_index(&options.index_dir)?;

        let started = stage_start();
        let census = census_fastq(&options.r1, options.r2.as_deref())?;
        stages.pass1_count = started.elapsed().as_millis() as u64;
        counts.0 = census.fragments;
        let design = derive_sampling_design(&census, index.target_groups.len(), profile)?;

        let started = stage_start();
        let mut fragments = FragmentReader::open(&options.r1, options.r2.as_deref())?;
        let mut accumulator = EvidenceAccumulator::new(&index.target_groups);
        let mut pass2_fragments = 0;
        align_fragments_bounded(
            &index.mmi_path,
            &index.contigs,
            options.threads,
            || loop {
                let Some(fragment) = fragments.next_fragment()? else {
                    return Ok(None);
                };
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
                match index
                    .bloom
                    .evaluate_fragment(&fragment.r1, fragment.r2.as_deref())
                {
                    GateEvaluation::Pass => {
                        counts.2 += 1;
                        return Ok(Some(fragment));
                    }
                    GateEvaluation::Negative => continue,
                    GateEvaluation::NotEvaluable => {
                        unevaluable_fragments += 1;
                        continue;
                    }
                }
            },
            |evidence| accumulator.accumulate_group_evidence(evidence),
        )?;
        if pass2_fragments != census.fragments {
            return Err(format!("FASTQ changed between passes: pass 1 counted {} fragments, pass 2 counted {pass2_fragments}", census.fragments));
        }
        if fragments.input_digest() != census.input_digest {
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
                unevaluable_fragments,
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

#[cfg(test)]
mod tests {
    use std::fs::File;
    use std::io::Write;
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;
    use crate::competitive_alignment::{CompetitiveAligner, FragmentAdjudication};
    use crate::fastq_input::FragmentReader;
    use crate::reference_index::{build_index, IndexOptions};

    static NEXT: AtomicU64 = AtomicU64::new(0);

    fn sequence(seed: u64, length: usize) -> Vec<u8> {
        let mut state = seed;
        (0..length)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn write_fasta(path: &Path, id: &str, sequence: &[u8]) {
        let mut file = File::create(path).unwrap();
        writeln!(file, ">{id}\n{}", String::from_utf8_lossy(sequence)).unwrap();
    }

    #[test]
    fn phase5_gate_counterfactual_quantifies_exhaustive_evidence_difference() {
        let root = std::env::temp_dir().join(format!(
            "viroflash-phase5-gate-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir(&root).unwrap();
        let host = sequence(17, 600);
        let target = sequence(91, 600);
        write_fasta(&root.join("host.fa"), "host", &host);
        write_fasta(&root.join("target.fa"), "target", &target);
        let index = build_index(&IndexOptions {
            host_fa: root.join("host.fa"),
            target_fa: root.join("target.fa"),
            out_dir: root.join("index"),
            threads: 1,
        })
        .unwrap();

        let exact = target[100..220].to_vec();
        let mut approximate = exact.clone();
        for offset in (10..approximate.len()).step_by(20) {
            approximate[offset] = match approximate[offset] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
        }
        let reads = [exact, approximate, host[100..220].to_vec()];
        let mut fastq = File::create(root.join("small.fastq")).unwrap();
        for (ordinal, read) in reads.iter().enumerate() {
            writeln!(
                fastq,
                "@fragment-{ordinal}\n{}\n+\n{}",
                String::from_utf8_lossy(read),
                "I".repeat(read.len())
            )
            .unwrap();
        }
        drop(fastq);

        let aligner = CompetitiveAligner::open(&index.mmi_path).unwrap();
        let mut reader = FragmentReader::open(&root.join("small.fastq"), None).unwrap();
        let mut gate_evaluations = Vec::new();
        let mut exhaustive = Vec::new();
        while let Some(fragment) = reader.next_fragment().unwrap() {
            gate_evaluations.push(index.bloom.evaluate_fragment(&fragment.r1, None));
            exhaustive.push(
                aligner
                    .align_fragment_competitively(&fragment, &index.contigs)
                    .unwrap()
                    .adjudication,
            );
        }

        assert_eq!(
            gate_evaluations,
            [
                GateEvaluation::Pass,
                GateEvaluation::Negative,
                GateEvaluation::Negative
            ]
        );
        assert!(matches!(exhaustive[0], FragmentAdjudication::Supporting(0)));
        assert!(matches!(
            exhaustive[1],
            FragmentAdjudication::NoTargetEvidence
        ));
        assert!(matches!(
            exhaustive[2],
            FragmentAdjudication::NoTargetEvidence
        ));
        let exhaustive_support = exhaustive
            .iter()
            .filter(|result| matches!(result, FragmentAdjudication::Supporting(0)))
            .count();
        let gated_support = exhaustive
            .iter()
            .zip(&gate_evaluations)
            .filter(|(result, gate)| {
                matches!(result, FragmentAdjudication::Supporting(0))
                    && **gate == GateEvaluation::Pass
            })
            .count();
        assert_eq!((gated_support, exhaustive_support), (1, 1));
        assert_eq!(
            gate_evaluations
                .iter()
                .filter(|evaluation| **evaluation == GateEvaluation::Pass)
                .count(),
            1
        );
        assert_eq!(exhaustive.len(), 3);
        let _ = std::fs::remove_dir_all(root);
    }
}
