use std::sync::mpsc::sync_channel;

struct CompletedBatch<R> {
    worker_index: usize,
    result: Result<R, String>,
}

pub(crate) fn process_batches_bounded<B, N, F, W, S, R>(
    threads: usize,
    next_batch: &mut N,
    make_worker: F,
    sink: &mut S,
) -> Result<(), String>
where
    B: Send,
    N: FnMut() -> Result<Option<B>, String>,
    F: Fn() -> W,
    W: FnMut(B) -> Result<R, String> + Send,
    S: FnMut(R) -> Result<(), String>,
    R: Send,
{
    let queue_capacity = threads;
    let (result_sender, result_receiver) = sync_channel(queue_capacity);
    std::thread::scope(|scope| -> Result<(), String> {
        let mut handles = Vec::new();
        let mut task_senders = Vec::with_capacity(threads);
        for worker_index in 0..threads {
            let (task_sender, task_receiver) = sync_channel::<B>(1);
            task_senders.push(task_sender);
            let mut process = make_worker();
            let result_sender = result_sender.clone();
            handles.push(scope.spawn(move || loop {
                let task = task_receiver.recv();
                let Ok(batch) = task else { break };
                let result =
                    std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| process(batch)));
                let panicked = result.is_err();
                let result = result.unwrap_or_else(|_| Err("Batch worker panicked".into()));
                if result_sender
                    .send(CompletedBatch {
                        worker_index,
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
            let mut active = 0;
            for sender in &task_senders {
                let Some(batch) = next_batch()? else {
                    break;
                };
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
                let analysis = completed.result?;
                if let Some(batch) = next_batch()? {
                    task_senders[completed.worker_index]
                        .send(batch)
                        .map_err(|_| "Batch worker queue closed".to_string())?;
                    active += 1;
                }
                sink(analysis)?;
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
    use std::sync::atomic::{AtomicUsize, Ordering};
    use std::sync::Barrier;

    #[test]
    fn parallel_work_keeps_live_batches_bounded() {
        struct Batch<'a>(&'a AtomicUsize);
        impl Drop for Batch<'_> {
            fn drop(&mut self) {
                self.0.fetch_sub(1, Ordering::SeqCst);
            }
        }
        let threads = 4;
        let barrier = Barrier::new(threads);
        let live = AtomicUsize::new(0);
        let mut emitted = 0;
        let mut completed = 0;
        process_batches_bounded(
            threads,
            &mut || {
                if emitted == 100 {
                    return Ok(None);
                }
                emitted += 1;
                assert!(live.fetch_add(1, Ordering::SeqCst) < threads + 1);
                Ok(Some(Batch(&live)))
            },
            || {
                let barrier = &barrier;
                let mut first = true;
                move |batch| {
                    if first {
                        barrier.wait();
                        first = false;
                    }
                    Ok(batch)
                }
            },
            &mut |_| {
                completed += 1;
                Ok(())
            },
        )
        .unwrap();
        assert_eq!(completed, 100);
        assert_eq!(live.load(Ordering::SeqCst), 0);
    }

    #[test]
    fn sink_failure_stops_after_at_most_one_refill() {
        let mut emitted = 0;
        let error = process_batches_bounded(
            4,
            &mut || {
                emitted += 1;
                Ok(Some(()))
            },
            || |_| Ok(()),
            &mut |_| Err("aggregation failed".into()),
        )
        .unwrap_err();
        assert_eq!(error, "aggregation failed");
        assert_eq!(emitted, 5);
    }

    #[test]
    fn worker_errors_and_panics_propagate() {
        for panics in [false, true] {
            let error = process_batches_bounded(
                4,
                &mut || Ok(Some(())),
                || {
                    move |_| {
                        assert!(!panics, "worker panic");
                        Err::<(), _>("worker failure".into())
                    }
                },
                &mut |_| Ok(()),
            )
            .unwrap_err();
            assert!(error.contains(if panics { "panicked" } else { "worker failure" }));
        }
    }
}
