use super::{CompetitiveAligner, FragmentAlignmentEvidence};
use crate::index::ReferenceContig;
use crate::sampling::SelectedBatch;
use crate::workers::process_batches_bounded;
use std::path::Path;

pub(crate) struct AnalysisWorkerConfig<'a> {
    pub index_path: &'a Path,
    pub contigs: &'a std::collections::HashMap<String, ReferenceContig>,
    pub threads: usize,
}

pub(crate) fn align_fragments_bounded<N, S>(
    config: AnalysisWorkerConfig<'_>,
    mut next_batch: N,
    mut sink: S,
) -> Result<(), String>
where
    N: FnMut() -> Result<Option<SelectedBatch>, String>,
    S: FnMut(FragmentAlignmentEvidence),
{
    let Some(first) = next_batch()? else {
        return Ok(());
    };
    let aligner = CompetitiveAligner::open(config.index_path)?;
    let mut first = Some(first);
    process_batches_bounded(
        config.threads,
        &mut || match first.take() {
            Some(batch) => Ok(Some(batch)),
            None => next_batch(),
        },
        || {
            let aligner = &aligner;
            let contigs = config.contigs;
            move |batch: SelectedBatch| {
                batch
                    .fragments()
                    .map(|f| aligner.align_fragment_competitively(&f, contigs))
                    .collect()
            }
        },
        &mut |alignments: Vec<FragmentAlignmentEvidence>| {
            for evidence in alignments {
                sink(evidence);
            }
            Ok(())
        },
    )
}
