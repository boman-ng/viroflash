use super::GateScratch;
use crate::profile::{SDUST_THRESHOLD, SDUST_WINDOW};

#[derive(Debug)]
pub(super) struct PerfectInterval {
    start: usize,
    finish: usize,
    score: i64,
    length: i64,
}

fn save_masked_regions(
    result: &mut Vec<(usize, usize)>,
    perfect: &mut Vec<PerfectInterval>,
    start: usize,
) {
    if perfect.is_empty() || perfect[perfect.len() - 1].start >= start {
        return;
    }
    let last = &perfect[perfect.len() - 1];
    if let Some(previous) = result.last_mut() {
        if last.start <= previous.1 {
            previous.1 = previous.1.max(last.finish);
        } else {
            result.push((last.start, last.finish));
        }
    } else {
        result.push((last.start, last.finish));
    }
    let mut keep = 0;
    for index in (0..perfect.len()).rev() {
        if perfect[index].start >= start {
            keep = index + 1;
            break;
        }
    }
    perfect.truncate(keep);
}

const SDUST_RING_CAPACITY: usize = 128;
const SDUST_RING_MASK: usize = SDUST_RING_CAPACITY - 1;

struct SdustRing {
    values: [u32; SDUST_RING_CAPACITY],
    head: usize,
    length: usize,
}

impl SdustRing {
    fn new() -> Self {
        Self {
            values: [0; SDUST_RING_CAPACITY],
            head: 0,
            length: 0,
        }
    }

    fn push(&mut self, value: u32) {
        self.values[(self.head + self.length) & SDUST_RING_MASK] = value;
        self.length += 1;
    }

    fn pop_front(&mut self) -> u32 {
        debug_assert!(self.length > 0);
        let value = self.values[self.head];
        self.head = (self.head + 1) & SDUST_RING_MASK;
        self.length -= 1;
        value
    }

    fn at(&self, index: usize) -> u32 {
        self.values[(self.head + index) & SDUST_RING_MASK]
    }
}

#[allow(clippy::too_many_arguments)]
fn shift_sdust_window(
    queue: &mut SdustRing,
    triplet: u32,
    active_length: &mut usize,
    window_score: &mut i64,
    suffix_score: &mut i64,
    window_counts: &mut [i64; 64],
    suffix_counts: &mut [i64; 64],
) {
    if queue.length >= SDUST_WINDOW - 2 {
        let symbol = (queue.pop_front() as usize) & 63;
        window_counts[symbol] -= 1;
        *window_score -= window_counts[symbol];
        if *active_length > queue.length {
            *active_length -= 1;
            suffix_counts[symbol] -= 1;
            *suffix_score -= suffix_counts[symbol];
        }
    }
    queue.push(triplet);
    *active_length += 1;
    let symbol = triplet as usize;
    *window_score += window_counts[symbol];
    window_counts[symbol] += 1;
    *suffix_score += suffix_counts[symbol];
    suffix_counts[symbol] += 1;
    if suffix_counts[symbol] * 10 > SDUST_THRESHOLD << 1 {
        loop {
            let removed = (queue.at(queue.length - *active_length) as usize) & 63;
            suffix_counts[removed] -= 1;
            *suffix_score -= suffix_counts[removed];
            *active_length -= 1;
            if removed == symbol {
                break;
            }
        }
    }
}

fn find_perfect_intervals(
    perfect: &mut Vec<PerfectInterval>,
    queue: &SdustRing,
    start: usize,
    active_length: usize,
    suffix_score: i64,
    suffix_counts: &[i64; 64],
) {
    let mut counts = *suffix_counts;
    let mut score = suffix_score;
    let mut maximum_score = 0;
    let mut maximum_length = 0;
    for index in (0..queue.length.saturating_sub(active_length)).rev() {
        let triplet = (queue.at(index) as usize) & 63;
        score += counts[triplet];
        counts[triplet] += 1;
        let new_score = score;
        let new_length = (queue.length - index - 1) as i64;
        if new_score * 10 > SDUST_THRESHOLD * new_length {
            let mut insertion = 0;
            while insertion < perfect.len() && perfect[insertion].start >= index + start {
                if maximum_score == 0
                    || perfect[insertion].score * maximum_length
                        > maximum_score * perfect[insertion].length
                {
                    maximum_score = perfect[insertion].score;
                    maximum_length = perfect[insertion].length;
                }
                insertion += 1;
            }
            if maximum_score == 0 || new_score * maximum_length >= maximum_score * new_length {
                maximum_score = new_score;
                maximum_length = new_length;
                perfect.insert(
                    insertion,
                    PerfectInterval {
                        start: index + start,
                        finish: queue.length + 2 + start,
                        score: new_score,
                        length: new_length,
                    },
                );
            }
        }
    }
}

fn sdust_base(symbol: u8) -> Option<u8> {
    match symbol.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

pub(super) fn sdust_intervals_into(sequence: &[u8], scratch: &mut GateScratch) {
    scratch.intervals.clear();
    scratch.perfect_intervals.clear();
    let mut queue = SdustRing::new();
    let mut window_counts = [0_i64; 64];
    let mut suffix_counts = [0_i64; 64];
    let mut suffix_score = 0;
    let mut window_score = 0;
    let mut active_length = 0;
    let mut contiguous_length = 0_usize;
    let mut triplet = 0_u32;

    for index in 0..=sequence.len() {
        if let Some(base) = sequence.get(index).copied().and_then(sdust_base) {
            contiguous_length += 1;
            triplet = ((triplet << 2) | u32::from(base)) & 0x3f;
            if contiguous_length >= 3 {
                let start = contiguous_length.saturating_sub(SDUST_WINDOW)
                    + (index + 1 - contiguous_length);
                save_masked_regions(
                    &mut scratch.intervals,
                    &mut scratch.perfect_intervals,
                    start,
                );
                shift_sdust_window(
                    &mut queue,
                    triplet,
                    &mut active_length,
                    &mut window_score,
                    &mut suffix_score,
                    &mut window_counts,
                    &mut suffix_counts,
                );
                if window_score * 10 > active_length as i64 * SDUST_THRESHOLD {
                    find_perfect_intervals(
                        &mut scratch.perfect_intervals,
                        &queue,
                        start,
                        active_length,
                        suffix_score,
                        &suffix_counts,
                    );
                }
            }
        } else {
            let mut start = contiguous_length.saturating_sub(SDUST_WINDOW.saturating_sub(1))
                + (index + 1 - contiguous_length);
            while !scratch.perfect_intervals.is_empty() {
                save_masked_regions(
                    &mut scratch.intervals,
                    &mut scratch.perfect_intervals,
                    start,
                );
                start += 1;
            }
            contiguous_length = 0;
            triplet = 0;
        }
    }
}
