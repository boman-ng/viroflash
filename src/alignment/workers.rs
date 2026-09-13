use super::{CompetitiveAligner, FragmentAlignmentEvidence};
use crate::fastq::FragmentBatch;
#[cfg(test)]
use crate::fastq::{Fragment, FASTQ_BATCH_RECORDS};
use crate::gate::{GateEvaluation, GateScratch, TargetKmerBloom};
use crate::index::ReferenceContig;
use crate::workers::process_batches_bounded;
#[cfg(test)]
use crate::workers::AlignmentRetention;
use std::path::Path;
use std::sync::OnceLock;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct ParallelAnalysisCounts {
    pub passed_fragments: u64,
    pub unevaluable_fragments: u64,
}

pub(crate) struct AnalysisWorkerConfig<'a> {
    pub index_path: &'a Path,
    pub contigs: &'a std::collections::HashMap<String, ReferenceContig>,
    pub threads: usize,
    pub bloom: Option<&'a TargetKmerBloom>,
    pub minimum_hits: usize,
    pub minimum_covered_bases: usize,
}

pub(crate) fn align_fragments_bounded<N, S>(
    config: AnalysisWorkerConfig<'_>,
    mut next_batch: N,
    mut sink: S,
) -> Result<ParallelAnalysisCounts, String>
where
    N: FnMut() -> Result<Option<FragmentBatch>, String>,
    S: FnMut(FragmentAlignmentEvidence),
{
    align_fragments_bounded_inner(
        &config,
        &mut next_batch,
        &mut sink,
        #[cfg(test)]
        |_| {},
    )
}

fn align_fragments_bounded_inner<N, S>(
    config: &AnalysisWorkerConfig<'_>,
    next_batch: &mut N,
    sink: &mut S,
    #[cfg(test)] observe_retention: impl FnMut(AlignmentRetention),
) -> Result<ParallelAnalysisCounts, String>
where
    N: FnMut() -> Result<Option<FragmentBatch>, String>,
    S: FnMut(FragmentAlignmentEvidence),
{
    let aligner = OnceLock::new();
    let mut counts = ParallelAnalysisCounts::default();
    process_batches_bounded(
        config.threads,
        next_batch,
        || {
            let aligner = &aligner;
            let mut scratch = GateScratch::default();
            move |batch: FragmentBatch| {
                let mut alignments = Vec::new();
                let mut unevaluable = 0;
                for fragment in batch.fragments() {
                    let f = fragment?;
                    if let Some(bloom) = config.bloom {
                        match bloom.evaluate_fragment(
                            f.r1,
                            f.r2,
                            config.minimum_hits,
                            config.minimum_covered_bases,
                            &mut scratch,
                        ) {
                            GateEvaluation::Pass => {}
                            GateEvaluation::Negative => continue,
                            GateEvaluation::NotEvaluable => {
                                unevaluable += 1;
                                continue;
                            }
                        }
                    }
                    let mapper = aligner
                        .get_or_init(|| CompetitiveAligner::open(config.index_path))
                        .as_ref()
                        .map_err(Clone::clone)?;
                    alignments.push(mapper.align_fragment_competitively(&f, config.contigs)?);
                }
                Ok((alignments, unevaluable))
            }
        },
        &mut |(alignments, unevaluable): (Vec<FragmentAlignmentEvidence>, u64)| {
            counts.passed_fragments += alignments.len() as u64;
            counts.unevaluable_fragments += unevaluable;
            for evidence in alignments {
                sink(evidence);
            }
            Ok(())
        },
        #[cfg(test)]
        observe_retention,
    )?;
    Ok(counts)
}

#[cfg(test)]
#[path = "worker_tests.rs"]
mod tests;
