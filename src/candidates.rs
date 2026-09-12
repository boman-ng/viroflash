//! Run-local candidate records, consumed only after the full input census.
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;

use crate::fastq::{Fragment, FragmentBatch, FragmentReader, InputCensus, FASTQ_BATCH_RECORDS};
use crate::gate::{GateEvaluation, GateScratch, TargetKmerBloom};
use crate::workers::process_batches_bounded;

pub(crate) struct CandidateCensus {
    pub input: InputCensus,
    pub candidates: u64,
    pub unevaluable: u64,
    pub spool_bytes: u64,
}

pub(crate) fn prescreen(
    reader: &mut FragmentReader,
    path: &Path,
    bloom: &TargetKmerBloom,
    minimum_hits: usize,
    minimum_covered_bases: usize,
    threads: usize,
    paired: bool,
) -> Result<CandidateCensus, String> {
    let mut writer = BufWriter::with_capacity(
        1 << 20,
        File::create(path).map_err(|e| format!("Cannot create candidate spool: {e}"))?,
    );
    let mut input = 0;
    let mut candidates = 0;
    let mut unevaluable = 0;
    let mut spool_bytes = 0;
    process_batches_bounded(
        threads,
        &mut || reader.next_batch(),
        || {
            let mut scratch = GateScratch::default();
            move |batch: &FragmentBatch| {
                let mut bytes = Vec::new();
                let mut retained = 0;
                let mut unscorable = 0;
                for fragment in batch.fragments() {
                    let fragment = fragment?;
                    match bloom.evaluate_fragment(
                        fragment.r1,
                        fragment.r2,
                        minimum_hits,
                        minimum_covered_bases,
                        &mut scratch,
                    ) {
                        GateEvaluation::Pass => {
                            retained += 1;
                            encode(&mut bytes, &fragment);
                        }
                        GateEvaluation::NotEvaluable => unscorable += 1,
                        GateEvaluation::Negative => {}
                    }
                }
                Ok((batch.len() as u64, retained, unscorable, bytes))
            }
        },
        &mut |(count, retained, unscorable, bytes): (u64, u64, u64, Vec<u8>)| {
            input += count;
            candidates += retained;
            unevaluable += unscorable;
            spool_bytes += bytes.len() as u64;
            writer
                .write_all(&bytes)
                .map_err(|e| format!("Cannot write candidate spool: {e}"))
        },
        #[cfg(test)]
        |_| {},
    )?;
    writer
        .flush()
        .map_err(|e| format!("Cannot flush candidate spool: {e}"))?;
    if input == 0 {
        return Err("FASTQ contains no fragments".into());
    }
    Ok(CandidateCensus {
        input: InputCensus {
            input_mode: if paired { "PE" } else { "SE" },
            fragments: input,
            input_digest: reader.input_digest()?,
            read_ends_per_fragment: if paired { 2 } else { 1 },
        },
        candidates,
        unevaluable,
        spool_bytes,
    })
}

fn encode(bytes: &mut Vec<u8>, fragment: &Fragment<'_>) {
    for value in [
        fragment.ordinal,
        fragment.id.len() as u64,
        fragment.r1.len() as u64,
        fragment.r2.map_or(0, |r| r.len()) as u64,
    ] {
        bytes.extend_from_slice(&value.to_le_bytes());
    }
    bytes.extend_from_slice(fragment.id.as_bytes());
    bytes.extend_from_slice(fragment.r1);
    if let Some(right) = fragment.r2 {
        bytes.extend_from_slice(right);
    }
}

pub(crate) struct CandidateReader {
    reader: BufReader<File>,
    remaining: u64,
    paired: bool,
}

impl CandidateReader {
    pub fn open(path: &Path, count: u64, paired: bool) -> Result<Self, String> {
        Ok(Self {
            reader: BufReader::with_capacity(
                1 << 20,
                File::open(path).map_err(|e| format!("Cannot open candidate spool: {e}"))?,
            ),
            remaining: count,
            paired,
        })
    }

    pub fn next_batch(&mut self) -> Result<Option<FragmentBatch>, String> {
        if self.remaining == 0 {
            return Ok(None);
        }
        let count = self.remaining.min(FASTQ_BATCH_RECORDS as u64);
        let mut batch = FragmentBatch::new(self.paired);
        let mut bytes = Vec::new();
        for _ in 0..count {
            let mut header = [0; 32];
            self.reader
                .read_exact(&mut header)
                .map_err(|e| format!("Cannot read candidate header: {e}"))?;
            let values: [u64; 4] = std::array::from_fn(|i| {
                u64::from_le_bytes(header[i * 8..i * 8 + 8].try_into().unwrap())
            });
            let id_len = usize::try_from(values[1]).map_err(|e| e.to_string())?;
            let left_len = usize::try_from(values[2]).map_err(|e| e.to_string())?;
            let right_len = usize::try_from(values[3]).map_err(|e| e.to_string())?;
            let size = id_len
                .checked_add(left_len)
                .and_then(|n| n.checked_add(right_len))
                .ok_or("Candidate length overflow")?;
            bytes.resize(size, 0);
            self.reader
                .read_exact(&mut bytes)
                .map_err(|e| format!("Cannot read candidate sequence: {e}"))?;
            batch.push(Fragment {
                ordinal: values[0],
                id: std::str::from_utf8(&bytes[..id_len])
                    .map_err(|e| format!("Invalid candidate ID: {e}"))?,
                r1: &bytes[id_len..id_len + left_len],
                r2: self.paired.then_some(&bytes[id_len + left_len..]),
            });
        }
        self.remaining -= count;
        Ok(Some(batch))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn spool_preserves_sparse_identity_and_rejects_truncation() {
        let path = std::env::temp_dir().join(format!("vf-candidates-{}", std::process::id()));
        for paired in [false, true] {
            let fragments = [99, 2, 1_000_003].map(|ordinal| Fragment {
                ordinal,
                id: "repeat-α/1",
                r1: b"ACGTN",
                r2: paired.then_some(b"TGCA".as_slice()),
            });
            let mut bytes = Vec::new();
            for fragment in &fragments {
                encode(&mut bytes, fragment);
            }
            std::fs::write(&path, &bytes).unwrap();
            let mut reader = CandidateReader::open(&path, 3, paired).unwrap();
            let batch = reader.next_batch().unwrap().unwrap();
            assert_eq!(
                batch.fragments().collect::<Result<Vec<_>, _>>().unwrap(),
                fragments
            );
            assert!(reader.next_batch().unwrap().is_none());
            std::fs::write(&path, &bytes[..bytes.len() - 1]).unwrap();
            assert!(CandidateReader::open(&path, 3, paired)
                .unwrap()
                .next_batch()
                .unwrap_err()
                .contains("Cannot read candidate sequence"));
        }
        std::fs::remove_file(path).unwrap();
    }
}
