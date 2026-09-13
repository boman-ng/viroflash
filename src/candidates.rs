//! One-pass input census and an in-memory bottom-k sample.
use crate::fastq::{FragmentBatch, FragmentReader, InputCensus};
use crate::gate::{GateEvaluation, GateScratch, TargetKmerBloom};
use crate::sampling::{BottomKSampler, FragmentKeys, RunMode, SelectedFragments};
use crate::workers::process_batches_bounded;

pub(crate) struct SamplingConfig<'a> {
    pub mode: RunMode,
    pub bloom: &'a TargetKmerBloom,
    pub keys: &'a FragmentKeys,
    pub minimum_hits: usize,
    pub minimum_covered_bases: usize,
    pub threads: usize,
    pub paired: bool,
    pub capacity: usize,
}

pub(crate) struct SampledInput {
    pub input: InputCensus,
    pub selected: SelectedFragments,
    pub population_fragments: u64,
    pub selected_fragments: u64,
    pub unevaluable: u64,
}

struct KeyedBatch {
    batch: FragmentBatch,
    keys: Vec<Option<u128>>,
    unevaluable: u64,
}

pub(crate) fn sample_input(
    reader: &mut FragmentReader,
    config: SamplingConfig<'_>,
) -> Result<SampledInput, String> {
    let mut sampler = BottomKSampler::new(config.capacity, config.paired);
    let mut input = 0;
    let mut population = 0;
    let mut unevaluable = 0;
    process_batches_bounded(
        config.threads,
        &mut || reader.next_batch(),
        || {
            let mut scratch = GateScratch::default();
            let config = &config;
            move |batch: FragmentBatch| {
                let mut keys = Vec::with_capacity(batch.len());
                let mut unevaluable = 0;
                for fragment in batch.fragments() {
                    let f = fragment?;
                    let retained = if config.mode == RunMode::Full {
                        match config.bloom.evaluate_fragment(
                            f.r1,
                            f.r2,
                            config.minimum_hits,
                            config.minimum_covered_bases,
                            &mut scratch,
                        ) {
                            GateEvaluation::Pass => true,
                            GateEvaluation::Negative => false,
                            GateEvaluation::NotEvaluable => {
                                unevaluable += 1;
                                false
                            }
                        }
                    } else {
                        true
                    };
                    keys.push(retained.then(|| config.keys.key(f.id, f.ordinal)));
                }
                Ok(KeyedBatch {
                    batch,
                    keys,
                    unevaluable,
                })
            }
        },
        &mut |result: KeyedBatch| {
            input += result.batch.len() as u64;
            unevaluable += result.unevaluable;
            for (f, key) in result.batch.fragments().zip(result.keys) {
                let f = f?;
                if let Some(key) = key {
                    population += 1;
                    sampler.consider(key, f);
                }
            }
            Ok(())
        },
        #[cfg(test)]
        |_| {},
    )?;
    if input == 0 {
        return Err("FASTQ contains no fragments".into());
    }
    let selected_fragments = sampler.len() as u64;
    Ok(SampledInput {
        input: InputCensus {
            input_mode: if config.paired { "PE" } else { "SE" },
            fragments: input,
            input_digest: reader.input_digest()?,
            read_ends_per_fragment: if config.paired { 2 } else { 1 },
        },
        selected: sampler.finish(),
        population_fragments: population,
        selected_fragments,
        unevaluable,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::index::reference::FastaRecord;
    use std::collections::BTreeSet;
    use std::io::Write;

    #[test]
    fn parallel_sampling_preserves_population_sequences_and_mode_containment() {
        let sequence = |mut state: u64| -> Vec<u8> {
            (0..600)
                .map(|_| {
                    state ^= state << 13;
                    state ^= state >> 7;
                    state ^= state << 17;
                    b"ACGT"[(state & 3) as usize]
                })
                .collect()
        };
        let target = sequence(91);
        let background = sequence(17);
        let bloom = TargetKmerBloom::build(
            &[FastaRecord {
                id: "target".into(),
                description: String::new(),
                sequence: target.clone(),
            }],
            21,
        );
        let path =
            std::env::temp_dir().join(format!("viroflash-sample-{}.fastq", std::process::id()));
        let mut file = std::fs::File::create(&path).unwrap();
        for i in 0..3000 {
            let read = if i % 7 == 0 {
                &target[..120]
            } else {
                &background[..120]
            };
            writeln!(
                file,
                "@read-{i}\n{}\n+\n{}",
                String::from_utf8_lossy(read),
                "I".repeat(120)
            )
            .unwrap();
        }
        drop(file);
        let keys = FragmentKeys::new("profile", "index");
        let (minimum_hits, minimum_covered_bases) =
            crate::alignment::CompetitiveAligner::short_read_chain_requirements(21).unwrap();
        let mut modes = Vec::new();
        for mode in [RunMode::Screen, RunMode::Full] {
            let mut expected = None;
            for threads in [1, 4] {
                let mut reader = FragmentReader::open(&path, None).unwrap();
                let mut result = sample_input(
                    &mut reader,
                    SamplingConfig {
                        mode,
                        bloom: &bloom,
                        keys: &keys,
                        minimum_hits,
                        minimum_covered_bases,
                        threads,
                        paired: false,
                        capacity: 31,
                    },
                )
                .unwrap();
                assert_eq!(result.input.fragments, 3000);
                assert_eq!(
                    result.population_fragments,
                    if mode == RunMode::Full { 429 } else { 3000 }
                );
                assert_eq!(result.selected_fragments, 31);
                let mut selected = Vec::new();
                let mut passing = BTreeSet::new();
                let mut scratch = GateScratch::default();
                while let Some(batch) = result.selected.next_batch() {
                    for f in batch.fragments() {
                        let f = f.unwrap();
                        let sequence = if f.ordinal % 7 == 0 {
                            &target[..120]
                        } else {
                            &background[..120]
                        };
                        assert_eq!(f.r1, sequence);
                        assert_eq!(f.id, format!("read-{}", f.ordinal));
                        selected.push(f.ordinal);
                        if bloom.evaluate_fragment(
                            f.r1,
                            None,
                            minimum_hits,
                            minimum_covered_bases,
                            &mut scratch,
                        ) == GateEvaluation::Pass
                        {
                            passing.insert(f.ordinal);
                        }
                    }
                }
                if let Some((previous_selected, previous_passing)) = &expected {
                    assert_eq!(&selected, previous_selected);
                    assert_eq!(&passing, previous_passing);
                } else {
                    expected = Some((selected, passing));
                }
            }
            modes.push(expected.unwrap());
        }
        assert!(modes[0].1.len() < 31);
        assert_eq!(modes[1].1.len(), 31);
        assert!(modes[0].1.is_subset(&modes[1].1));
        std::fs::remove_file(path).unwrap();
    }
}
