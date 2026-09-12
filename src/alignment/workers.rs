use super::{CompetitiveAligner, FragmentAlignmentEvidence};
use crate::fastq::FragmentBatch;
#[cfg(test)]
use crate::fastq::{Fragment, FASTQ_BATCH_RECORDS};
use crate::index::ReferenceContig;
use crate::sampling::FragmentSelector;
use crate::workers::process_batches_bounded;
#[cfg(test)]
use crate::workers::AlignmentRetention;
use std::path::Path;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(crate) struct ParallelAnalysisCounts {
    pub input_fragments: u64,
    pub selected_fragments: u64,
}

pub(crate) struct AnalysisWorkerConfig<'a> {
    pub index_path: &'a Path,
    pub contigs: &'a std::collections::HashMap<String, ReferenceContig>,
    pub threads: usize,
    pub selector: &'a FragmentSelector,
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
    let aligner = CompetitiveAligner::open(config.index_path)?;
    let mut counts = ParallelAnalysisCounts::default();
    process_batches_bounded(
        config.threads,
        next_batch,
        || {
            let aligner = &aligner;
            move |batch: &FragmentBatch| {
                let mut selected = Vec::new();
                for fragment in batch.fragments() {
                    let fragment = fragment?;
                    if config.selector.includes(fragment.id, fragment.ordinal) {
                        selected
                            .push(aligner.align_fragment_competitively(&fragment, config.contigs)?);
                    }
                }
                Ok((batch.len() as u64, selected))
            }
        },
        &mut |(input, selected): (u64, Vec<FragmentAlignmentEvidence>)| {
            counts.input_fragments += input;
            counts.selected_fragments += selected.len() as u64;
            for evidence in selected {
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
