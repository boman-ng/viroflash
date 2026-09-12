use crate::fastq::FragmentBatch;
#[cfg(test)]
use crate::fastq::FASTQ_BATCH_RECORDS;
use std::sync::mpsc::sync_channel;

#[cfg(test)]
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct AlignmentRetention {
    pub pending_fragments: usize,
    pub pending_sequence_bytes: usize,
    pub maximum_pending_fragments: usize,
    pub maximum_pending_sequence_bytes: usize,
    pub maximum_fragment_bytes: usize,
    pub total_sequence_bytes: usize,
}

#[cfg(test)]
impl AlignmentRetention {
    fn claim(
        &mut self,
        fragment_count: usize,
        sequence_bytes: usize,
        maximum_fragment_bytes: usize,
        threads: usize,
    ) {
        self.pending_fragments += fragment_count;
        self.pending_sequence_bytes += sequence_bytes;
        self.maximum_pending_fragments = self.maximum_pending_fragments.max(self.pending_fragments);
        self.maximum_pending_sequence_bytes = self
            .maximum_pending_sequence_bytes
            .max(self.pending_sequence_bytes);
        self.maximum_fragment_bytes = self.maximum_fragment_bytes.max(maximum_fragment_bytes);
        self.total_sequence_bytes += sequence_bytes;
        let maximum_pending_fragments = threads.saturating_mul(FASTQ_BATCH_RECORDS);
        assert!(self.pending_fragments <= maximum_pending_fragments);
        assert!(
            self.pending_sequence_bytes
                <= maximum_pending_fragments.saturating_mul(self.maximum_fragment_bytes)
        );
    }

    fn release(&mut self, fragment_count: usize, sequence_bytes: usize) {
        self.pending_fragments -= fragment_count;
        self.pending_sequence_bytes -= sequence_bytes;
    }
}

struct CompletedBatch<R> {
    worker_index: usize,
    #[cfg(test)]
    fragment_count: usize,
    #[cfg(test)]
    sequence_bytes: usize,
    result: Result<R, String>,
}

pub(crate) fn process_batches_bounded<N, F, W, S, R>(
    threads: usize,
    next_batch: &mut N,
    make_worker: F,
    sink: &mut S,
    #[cfg(test)] mut observe_retention: impl FnMut(AlignmentRetention),
) -> Result<(), String>
where
    N: FnMut() -> Result<Option<FragmentBatch>, String>,
    F: Fn() -> W,
    W: FnMut(&FragmentBatch) -> Result<R, String> + Send,
    S: FnMut(R) -> Result<(), String>,
    R: Send,
{
    let queue_capacity = threads;
    let (result_sender, result_receiver) = sync_channel(queue_capacity);
    std::thread::scope(|scope| -> Result<(), String> {
        let mut handles = Vec::new();
        let mut task_senders = Vec::with_capacity(threads);
        for worker_index in 0..threads {
            let (task_sender, task_receiver) = sync_channel::<FragmentBatch>(1);
            task_senders.push(task_sender);
            let mut process = make_worker();
            let result_sender = result_sender.clone();
            handles.push(scope.spawn(move || loop {
                let task = task_receiver.recv();
                let Ok(batch) = task else { break };
                #[cfg(test)]
                let fragment_count = batch.len();
                #[cfg(test)]
                let sequence_bytes = batch.sequence_bytes();
                let result =
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| process(&batch)));
                drop(batch);
                let panicked = result.is_err();
                let result = result.unwrap_or_else(|_| Err("Batch worker panicked".into()));
                if result_sender
                    .send(CompletedBatch {
                        worker_index,
                        #[cfg(test)]
                        fragment_count,
                        #[cfg(test)]
                        sequence_bytes,
                        result,
                    })
                    .is_err()
                {
                    break;
                }
                if panicked {
                    break;
                }
            }));
        }
        drop(result_sender);
        let processing = (|| -> Result<(), String> {
            #[cfg(test)]
            let mut retention = AlignmentRetention::default();

            let mut active = 0;
            for sender in &task_senders {
                let Some(batch) = next_batch()? else {
                    break;
                };
                #[cfg(test)]
                {
                    let sequence_bytes = batch.sequence_bytes();
                    let maximum_fragment_bytes = batch.maximum_fragment_bytes();
                    retention.claim(batch.len(), sequence_bytes, maximum_fragment_bytes, threads);
                    observe_retention(retention);
                }
                sender
                    .send(batch)
                    .map_err(|_| "Batch worker queue closed".to_string())?;
                active += 1;
            }
            while active > 0 {
                let completed = result_receiver
                    .recv()
                    .map_err(|_| "Batch worker result queue closed".to_string())?;
                active -= 1;
                #[cfg(test)]
                {
                    retention.release(completed.fragment_count, completed.sequence_bytes);
                    observe_retention(retention);
                }
                let analysis = completed.result?;
                sink(analysis)?;
                if let Some(batch) = next_batch()? {
                    #[cfg(test)]
                    {
                        let sequence_bytes = batch.sequence_bytes();
                        let maximum_fragment_bytes = batch.maximum_fragment_bytes();
                        retention.claim(
                            batch.len(),
                            sequence_bytes,
                            maximum_fragment_bytes,
                            threads,
                        );
                        observe_retention(retention);
                    }
                    task_senders[completed.worker_index]
                        .send(batch)
                        .map_err(|_| "Batch worker queue closed".to_string())?;
                    active += 1;
                }
            }
            #[cfg(test)]
            {
                assert_eq!(retention.pending_fragments, 0);
                assert_eq!(retention.pending_sequence_bytes, 0);
            }
            Ok(())
        })();
        drop(task_senders);
        let mut join_error = None;
        for handle in handles {
            if handle.join().is_err() {
                join_error = Some("Batch worker panicked".to_string());
            }
        }
        processing.and_then(|()| join_error.map_or(Ok(()), Err))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::Fragment;

    #[test]
    fn sink_failure_joins_workers_without_waiting_for_more_input() {
        let mut emitted = 0;
        let error = process_batches_bounded(
            4,
            &mut || {
                emitted += 1;
                Ok(Some(FragmentBatch::from_fragments(&[Fragment {
                    ordinal: emitted,
                    id: "read",
                    r1: b"ACGT",
                    r2: None,
                }])))
            },
            || |_: &FragmentBatch| Ok(()),
            &mut |_| Err("candidate write failed".into()),
            |_| {},
        )
        .unwrap_err();
        assert_eq!(error, "candidate write failed");
        assert_eq!(emitted, 4);
    }
}
