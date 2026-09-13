use std::path::{Path, PathBuf};
use std::time::Instant;

use crate::alignment::{align_fragments_bounded, AnalysisWorkerConfig, CompetitiveAligner};
use crate::candidates::{sample_input, SamplingConfig};
use crate::evidence::EvidenceAccumulator;
use crate::fastq::FragmentReader;
use crate::index::load_index;
use crate::profile::AnalysisProfile;
use crate::report::{build_evidence_report, write_report_csv, write_report_html, ReportInputs};
use crate::sampling::{sample_capacity, FragmentKeys, Precision, RunMode, SamplingDesign};
use crate::telemetry::{write_perf_json, PerformanceMonitor, StageTimes};

#[derive(Debug, Clone)]
pub struct RunOptions {
    pub r1: PathBuf,
    pub r2: Option<PathBuf>,
    pub index_dir: PathBuf,
    pub out_dir: PathBuf,
    pub threads: usize,
    pub precision: Precision,
    pub mode: RunMode,
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
    let result: Result<_, String> = (|| {
        std::fs::create_dir(&staging)
            .map_err(|error| format!("Cannot create {}: {error}", staging.display()))?;
        let profile = AnalysisProfile::FROZEN;
        let started = Instant::now();
        let index = load_index(&options.index_dir, options.threads)?;
        stages.index_load_ms = started.elapsed().as_millis() as u64;
        let capacity = sample_capacity(index.target_groups.len(), profile, options.precision)?;
        let (minimum_hits, minimum_covered_bases) =
            CompetitiveAligner::short_read_chain_requirements(profile.kmer_length)?;
        let keys = FragmentKeys::new(&index.profile_digest, &index.index_digest);
        let started = Instant::now();
        let mut reader = FragmentReader::open(&options.r1, options.r2.as_deref())?;
        let mut sampled = sample_input(
            &mut reader,
            SamplingConfig {
                mode: options.mode,
                bloom: &index.bloom,
                keys: &keys,
                minimum_hits,
                minimum_covered_bases,
                threads: options.threads,
                paired: options.r2.is_some(),
                capacity,
            },
        )?;
        drop(reader);
        stages.scan_sample_ms = started.elapsed().as_millis() as u64;
        stages.sampling_peak_buffered_fragments = sampled.selected_fragments;
        let census = sampled.input;
        counts.0 = census.fragments;
        counts.1 = sampled.selected_fragments;
        unevaluable_fragments = sampled.unevaluable;
        let design = SamplingDesign {
            precision: options.precision,
            population_fragments: sampled.population_fragments,
            selection_probability: if sampled.population_fragments == 0 {
                0.0
            } else {
                counts.1 as f64 / sampled.population_fragments as f64
            },
            sample_capacity: capacity as u64,
        };
        let mut accumulator = EvidenceAccumulator::new(&index.target_groups);
        let started = Instant::now();
        let bloom = if options.mode == RunMode::Screen {
            Some(index.bloom)
        } else {
            drop(index.bloom);
            None
        };
        let analyzed = align_fragments_bounded(
            AnalysisWorkerConfig {
                index_path: &index.mmi_path,
                contigs: &index.contigs,
                threads: options.threads,
                bloom: bloom.as_ref(),
                minimum_hits,
                minimum_covered_bases,
            },
            || Ok(sampled.selected.next_batch()),
            |evidence| accumulator.accumulate_group_evidence(evidence),
        )?;
        drop(bloom);
        unevaluable_fragments += analyzed.unevaluable_fragments;
        counts.2 = if options.mode == RunMode::Full {
            sampled.population_fragments
        } else {
            analyzed.passed_fragments
        };
        counts.3 = accumulator.aligned_fragments;
        stages.selected_analysis_ms = started.elapsed().as_millis() as u64;
        let started = Instant::now();
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
                mode: options.mode,
            },
            accumulator,
        )?;
        write_report_csv(&staging.join("report.csv"), &report)?;
        write_report_html(&staging.join("report.html"), &report)?;
        stages.report_write = started.elapsed().as_millis() as u64;
        Ok((census, report))
    })();

    match result {
        Ok((census, report)) => {
            let target_signal_rows = report
                .research_rows
                .iter()
                .filter(|row| !row[2].is_empty())
                .count();
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
    use crate::alignment::CompetitiveAligner;
    use crate::evidence::EvidenceAccumulator;
    use crate::fastq::{census_fastq, FragmentReader};
    use crate::gate::{GateEvaluation, GateScratch};
    use crate::index::{build_index, IndexOptions};
    use crate::report::{build_evidence_report, EvidenceReport, ReportInputs};
    use crate::sampling::SamplingDesign;

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

    fn gate_counterfactual_report(
        fastq: &Path,
        index: &crate::index::ReferenceIndex,
        exhaustive: bool,
    ) -> EvidenceReport {
        let census = census_fastq(fastq, None).unwrap();
        let aligner = CompetitiveAligner::open(&index.mmi_path).unwrap();
        let mut reader = FragmentReader::open(fastq, None).unwrap();
        let mut accumulator = EvidenceAccumulator::new(&index.target_groups);
        let mut gate_scratch = GateScratch::default();
        let (minimum_hits, minimum_covered_bases) =
            CompetitiveAligner::short_read_chain_requirements(AnalysisProfile::FROZEN.kmer_length)
                .unwrap();
        let mut submitted = 0;
        while let Some(batch) = reader.next_batch().unwrap() {
            for fragment in batch.fragments() {
                let fragment = fragment.unwrap();
                if !exhaustive
                    && index.bloom.evaluate_fragment(
                        fragment.r1,
                        None,
                        minimum_hits,
                        minimum_covered_bases,
                        &mut gate_scratch,
                    ) != GateEvaluation::Pass
                {
                    continue;
                }
                submitted += 1;
                accumulator.accumulate_group_evidence(
                    aligner
                        .align_fragment_competitively(&fragment, &index.contigs)
                        .unwrap(),
                );
            }
        }
        build_evidence_report(
            ReportInputs {
                sample_id: "gate-counterfactual".into(),
                census: &census,
                design: SamplingDesign {
                    precision: Precision::Standard,
                    population_fragments: census.fragments,
                    selection_probability: 1.0,
                    sample_capacity: census.fragments,
                },
                selected_fragments: census.fragments,
                prescreen_passed_fragments: submitted,
                profile: AnalysisProfile::FROZEN,
                index_digest: index.index_digest.clone(),
                unevaluable_fragments: 0,
                mode: RunMode::Full,
            },
            accumulator,
        )
        .unwrap()
    }

    #[test]
    fn gate_counterfactual_quantifies_exhaustive_evidence_difference() {
        let root = std::env::temp_dir().join(format!(
            "viroflash-gate-{}-{}",
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

        let gated = gate_counterfactual_report(&root.join("small.fastq"), &index, false);
        let exhaustive = gate_counterfactual_report(&root.join("small.fastq"), &index, true);
        assert_eq!(gated.schema_id, exhaustive.schema_id);
        assert_eq!(gated.target_signals, exhaustive.target_signals);

        let mut normalized_exhaustive_run = exhaustive.run.clone();
        normalized_exhaustive_run.prescreen_passed_fragments = gated.run.prescreen_passed_fragments;
        normalized_exhaustive_run.aligned_fragments = gated.run.aligned_fragments;
        normalized_exhaustive_run.unassigned_fragments = gated.run.unassigned_fragments;
        assert_eq!(gated.run, normalized_exhaustive_run);
        assert_eq!(
            (
                gated.run.prescreen_passed_fragments,
                gated.run.aligned_fragments,
                gated.run.unassigned_fragments,
            ),
            (1, 1, 0)
        );
        assert_eq!(
            (
                exhaustive.run.prescreen_passed_fragments,
                exhaustive.run.aligned_fragments,
                exhaustive.run.unassigned_fragments,
            ),
            (3, 2, 1)
        );
        let _ = std::fs::remove_dir_all(root);
    }
}
