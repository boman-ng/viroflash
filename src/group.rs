//! Sequence-only detection groups using fixed candidate screening and representative-star confirmation.
//!
//! Edges only identify members that one representative may cover; connected components are not
//! treated as equivalence classes.

use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::sync::atomic::{AtomicUsize, Ordering as AtomicOrdering};
use std::sync::mpsc::sync_channel;
use std::thread;

use minimap2::{Aligner, Built, Mapping, Strand};

const K: usize = 21;
const FRAC_SCALE: u64 = 64;
const BOTTOM_K: usize = 256;
const JACCARD_NUMERATOR: usize = 1;
const JACCARD_DENOMINATOR: usize = 8;
const MASH_ANI_NUMERATOR: usize = 96;
const MASH_ANI_DENOMINATOR: usize = 100;
const ALIGN_IDENTITY_NUMERATOR: usize = 97;
const ALIGN_IDENTITY_DENOMINATOR: usize = 100;
const COVERAGE_NUMERATOR: usize = 95;
const COVERAGE_DENOMINATOR: usize = 100;
const POSITION_BINS: usize = 16;
const MIN_POSITION_BINS: usize = 14;

/// One detection group. `representative` and `members` index the caller's input slice.
/// `members` includes the representative, and every input index appears exactly once.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DetectionGroup {
    pub representative: usize,
    pub members: Vec<usize>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Token {
    // splitmix64 is bijective over u64 while canonical k=21 codes occupy 42 bits. Retaining the
    // code keeps token identity and global `(df, canonical code)` order independent of that fact.
    hash: u64,
    code: u64,
}

#[derive(Clone, Copy, Debug)]
struct SketchEntry {
    token: Token,
    position: usize,
    forward_canonical: bool,
    occurrences: u32,
}

#[derive(Debug)]
struct ExactClass {
    source_indices: Vec<usize>,
    sequence: Vec<u8>,
    sample: Vec<Token>,
    sketch: Vec<SketchEntry>,
}

type ClassRank<'a> = (usize, std::cmp::Reverse<usize>, [u8; 32], &'a [u8]);

/// Build groups with fixed k=21, scaled=64 FracMinHash, bottom-k Mash, and minimap2 `asm5`.
///
/// `threads == 0` runs inline; positive values bound precomputation and confirmation workers.
pub fn construct_detection_groups(
    records: &[(String, Vec<u8>)],
    threads: usize,
) -> Result<Vec<DetectionGroup>, String> {
    let mut classes = exact_classes(records);
    if classes.is_empty() {
        return Ok(Vec::new());
    }
    populate_signatures(&mut classes, threads)?;

    let lengths: Vec<usize> = classes.iter().map(|class| class.sequence.len()).collect();
    // PPJoin is the sample's final consumer. Transfer ownership to avoid retaining three copies of
    // canonical FracMinHash values across classes, the call site, and this function.
    let samples: Vec<Vec<Token>> = classes
        .iter_mut()
        .map(|class| std::mem::take(&mut class.sample))
        .collect();
    let candidate_pairs = ppjoin_pairs(samples);
    let accepted_pairs = filter_candidate_pairs(&classes, &lengths, &candidate_pairs, threads)?;
    let mut adjacency = vec![Vec::new(); classes.len()];
    for (left, right) in accepted_pairs {
        adjacency[left].push(right);
        adjacency[right].push(left);
    }
    for neighbors in &mut adjacency {
        neighbors.sort_unstable();
    }
    confirm_stars(&classes, &adjacency, threads)
}

fn exact_classes(records: &[(String, Vec<u8>)]) -> Vec<ExactClass> {
    // Digests only select candidate buckets; full sequences are still compared within each bucket.
    // This avoids copying large genomes into BTreeMap keys and prevents collision-based collapsing.
    let mut by_digest = BTreeMap::<[u8; 32], Vec<usize>>::new();
    let mut classes = Vec::<ExactClass>::new();
    for (source_index, (_, raw_sequence)) in records.iter().enumerate() {
        let sequence = normalized_sequence(raw_sequence);
        let digest = *blake3::hash(&sequence).as_bytes();
        let existing = by_digest.get(&digest).and_then(|candidates| {
            candidates
                .iter()
                .copied()
                .find(|&index| classes[index].sequence == sequence)
        });
        if let Some(class_index) = existing {
            classes[class_index].source_indices.push(source_index);
            continue;
        }
        let class_index = classes.len();
        classes.push(ExactClass {
            source_indices: vec![source_index],
            sequence,
            sample: Vec::new(),
            sketch: Vec::new(),
        });
        by_digest.entry(digest).or_default().push(class_index);
    }
    classes
}

fn normalized_sequence(sequence: &[u8]) -> Vec<u8> {
    sequence.iter().map(u8::to_ascii_uppercase).collect()
}

fn populate_signatures(classes: &mut [ExactClass], threads: usize) -> Result<(), String> {
    if threads == 0 || classes.len() == 1 {
        for class in classes {
            class.sample = frac_sample(&class.sequence);
            class.sketch = bottom_sketch(&class.sequence);
        }
        return Ok(());
    }
    let workers = threads.min(classes.len());
    let next = AtomicUsize::new(0);
    let (sender, receiver) = sync_channel::<(usize, Vec<Token>, Vec<SketchEntry>)>(workers * 2);
    let sequences: Vec<_> = classes
        .iter()
        .map(|class| class.sequence.as_slice())
        .collect();
    let mut signatures = vec![None; classes.len()];
    thread::scope(|scope| {
        for _ in 0..workers {
            let sender = sender.clone();
            let next = &next;
            let sequences = &sequences;
            scope.spawn(move || loop {
                let index = next.fetch_add(1, AtomicOrdering::Relaxed);
                let Some(sequence) = sequences.get(index) else {
                    break;
                };
                if sender
                    .send((index, frac_sample(sequence), bottom_sketch(sequence)))
                    .is_err()
                {
                    break;
                }
            });
        }
        drop(sender);
        for (index, sample, sketch) in receiver {
            let Some(signature) = signatures.get_mut(index) else {
                return Err("Detection-group signature index out of bounds".to_string());
            };
            *signature = Some((sample, sketch));
        }
        Ok(())
    })?;
    for (class, signature) in classes.iter_mut().zip(signatures) {
        let Some((sample, sketch)) = signature else {
            return Err("Detection-group signature is incomplete".to_string());
        };
        class.sample = sample;
        class.sketch = sketch;
    }
    Ok(())
}

fn dna_bits(base: u8) -> Option<u64> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

fn splitmix64(mut value: u64) -> u64 {
    value = value.wrapping_add(0x9e37_79b9_7f4a_7c15);
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    value ^ (value >> 31)
}

fn frac_sample(sequence: &[u8]) -> Vec<Token> {
    if sequence.len() < K {
        return Vec::new();
    }
    let intervals = crate::prescreen::sdust_intervals(
        sequence,
        crate::prescreen::SDUST_W,
        crate::prescreen::SDUST_T,
    );
    let mask = (1_u64 << (2 * K)) - 1;
    let threshold = u64::MAX / FRAC_SCALE;
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mut filled = 0_usize;
    let mut interval_index = 0_usize;
    let mut tokens = Vec::with_capacity(sequence.len() / FRAC_SCALE as usize + 1);
    // Add k-1 wrapped starts so rotations of one circular sequence share candidates; final alignment
    // still performs strict confirmation.
    for end in 0..sequence.len() + K - 1 {
        let base = sequence[end % sequence.len()];
        let Some(bits) = dna_bits(base) else {
            forward = 0;
            reverse = 0;
            filled = 0;
            continue;
        };
        forward = ((forward << 2) | bits) & mask;
        reverse = (reverse >> 2) | ((3 - bits) << (2 * K - 2));
        filled = (filled + 1).min(K);
        if filled != K {
            continue;
        }
        let start = end + 1 - K;
        if start >= sequence.len() {
            break;
        }
        while interval_index < intervals.len() && intervals[interval_index].1 <= start {
            interval_index += 1;
        }
        if interval_index < intervals.len() && intervals[interval_index].0 <= start {
            continue;
        }
        let code = forward.min(reverse);
        let hash = splitmix64(code);
        if hash < threshold {
            tokens.push(Token { hash, code });
        }
    }
    tokens.sort_unstable();
    tokens.dedup();
    tokens
}

fn bottom_sketch(sequence: &[u8]) -> Vec<SketchEntry> {
    if sequence.len() < K {
        return Vec::new();
    }
    let mask = (1_u64 << (2 * K)) - 1;
    let mut forward = 0_u64;
    let mut reverse = 0_u64;
    let mut filled = 0_usize;
    let mut entries = BTreeMap::<Token, SketchEntry>::new();
    for (end, &base) in sequence.iter().enumerate() {
        let Some(bits) = dna_bits(base) else {
            forward = 0;
            reverse = 0;
            filled = 0;
            continue;
        };
        forward = ((forward << 2) | bits) & mask;
        reverse = (reverse >> 2) | ((3 - bits) << (2 * K - 2));
        filled = (filled + 1).min(K);
        if filled != K {
            continue;
        }
        let code = forward.min(reverse);
        let token = Token {
            hash: splitmix64(code),
            code,
        };
        if let Some(entry) = entries.get_mut(&token) {
            entry.occurrences = entry.occurrences.saturating_add(1);
        } else {
            entries.insert(
                token,
                SketchEntry {
                    token,
                    position: end + 1 - K,
                    forward_canonical: forward <= reverse,
                    occurrences: 1,
                },
            );
        }
    }
    let mut sketch: Vec<_> = entries.into_values().collect();
    sketch.sort_unstable_by_key(|entry| entry.token);
    sketch.truncate(BOTTOM_K);
    sketch
}

/// AllPairs/PPJoin prefix join. Candidates arise only from rare-token prefixes without enumerating
/// all global posting pairs; the final intersection lower bound is `ceil((|A| + |B|) / 9)`.
fn ppjoin_pairs(mut samples: Vec<Vec<Token>>) -> Vec<(usize, usize)> {
    for sample in &mut samples {
        sample.sort_unstable();
        sample.dedup();
    }
    let mut frequencies = BTreeMap::<Token, usize>::new();
    for sample in &samples {
        for &token in sample {
            *frequencies.entry(token).or_default() += 1;
        }
    }
    let mut ordered = samples.to_vec();
    for sample in &mut ordered {
        sample.sort_unstable_by_key(|token| (frequencies[token], token.code));
    }
    let mut sequence_order: Vec<usize> = (0..samples.len()).collect();
    sequence_order.sort_unstable_by_key(|&index| (samples[index].len(), index));
    let mut inverted = BTreeMap::<Token, Vec<usize>>::new();
    let mut seen = vec![usize::MAX; samples.len()];
    let mut pairs = Vec::new();
    for (epoch, index) in sequence_order.into_iter().enumerate() {
        let current_len = samples[index].len();
        if current_len == 0 {
            continue;
        }
        let needed = ceil_div(JACCARD_NUMERATOR * current_len, JACCARD_DENOMINATOR).max(1);
        let prefix_len = current_len.saturating_sub(needed) + 1;
        let mut candidates = Vec::new();
        for token in ordered[index].iter().take(prefix_len) {
            if let Some(previous) = inverted.get(token) {
                for &candidate in previous {
                    let length_lower_bound =
                        ceil_div(JACCARD_NUMERATOR * current_len, JACCARD_DENOMINATOR);
                    if samples[candidate].len() >= length_lower_bound && seen[candidate] != epoch {
                        seen[candidate] = epoch;
                        candidates.push(candidate);
                    }
                }
            }
        }
        candidates.sort_unstable();
        for candidate in candidates {
            if sampled_jaccard_at_least(&samples[index], &samples[candidate]) {
                pairs.push((candidate.min(index), candidate.max(index)));
            }
        }
        for token in ordered[index].iter().take(prefix_len) {
            inverted.entry(*token).or_default().push(index);
        }
    }
    pairs.sort_unstable();
    pairs.dedup();
    pairs
}

fn ceil_div(numerator: usize, denominator: usize) -> usize {
    numerator / denominator + usize::from(!numerator.is_multiple_of(denominator))
}

fn sampled_jaccard_at_least(left: &[Token], right: &[Token]) -> bool {
    let mut left_index = 0_usize;
    let mut right_index = 0_usize;
    let mut intersection = 0_usize;
    while left_index < left.len() && right_index < right.len() {
        match left[left_index].cmp(&right[right_index]) {
            Ordering::Less => left_index += 1,
            Ordering::Greater => right_index += 1,
            Ordering::Equal => {
                intersection += 1;
                left_index += 1;
                right_index += 1;
            }
        }
    }
    let alpha = ceil_div_u128(
        JACCARD_NUMERATOR as u128 * (left.len() as u128 + right.len() as u128),
        (JACCARD_NUMERATOR + JACCARD_DENOMINATOR) as u128,
    );
    intersection as u128 >= alpha
}

fn ceil_div_u128(numerator: u128, denominator: u128) -> u128 {
    numerator / denominator + u128::from(!numerator.is_multiple_of(denominator))
}

fn compatible_lengths(left: usize, right: usize) -> bool {
    let shorter = left.min(right);
    let longer = left.max(right);
    shorter.saturating_mul(COVERAGE_DENOMINATOR) >= longer.saturating_mul(COVERAGE_NUMERATOR)
}

fn filter_candidate_pairs(
    classes: &[ExactClass],
    lengths: &[usize],
    pairs: &[(usize, usize)],
    threads: usize,
) -> Result<Vec<(usize, usize)>, String> {
    if pairs.is_empty() {
        return Ok(Vec::new());
    }
    if threads == 0 {
        return Ok(pairs
            .iter()
            .copied()
            .filter(|&(left, right)| {
                compatible_lengths(lengths[left], lengths[right])
                    && mash_ani_at_least(&classes[left].sketch, &classes[right].sketch)
                    && collinear_coverage(
                        &classes[left].sketch,
                        &classes[right].sketch,
                        lengths[left],
                        lengths[right],
                    )
            })
            .collect());
    }
    let workers = threads.max(1).min(pairs.len());
    let next = AtomicUsize::new(0);
    let (sender, receiver) = sync_channel::<Vec<(usize, usize)>>(workers * 2);
    thread::scope(|scope| {
        for _ in 0..workers {
            let sender = sender.clone();
            let next = &next;
            scope.spawn(move || {
                let mut accepted = Vec::new();
                loop {
                    let start = next.fetch_add(1_024, AtomicOrdering::Relaxed);
                    if start >= pairs.len() {
                        break;
                    }
                    for &(left, right) in &pairs[start..(start + 1_024).min(pairs.len())] {
                        if compatible_lengths(lengths[left], lengths[right])
                            && mash_ani_at_least(&classes[left].sketch, &classes[right].sketch)
                            && collinear_coverage(
                                &classes[left].sketch,
                                &classes[right].sketch,
                                lengths[left],
                                lengths[right],
                            )
                        {
                            accepted.push((left, right));
                        }
                    }
                }
                let _ = sender.send(accepted);
            });
        }
        drop(sender);
        let mut accepted = Vec::new();
        for mut local in receiver {
            accepted.append(&mut local);
        }
        accepted.sort_unstable();
        Ok(accepted)
    })
}

fn mash_ani_at_least(left: &[SketchEntry], right: &[SketchEntry]) -> bool {
    let Some(left_last) = left.last() else {
        return false;
    };
    let Some(right_last) = right.last() else {
        return false;
    };
    let cutoff = left_last.token.min(right_last.token);
    let left_count = left.partition_point(|entry| entry.token <= cutoff);
    let right_count = right.partition_point(|entry| entry.token <= cutoff);
    let mut left_index = 0_usize;
    let mut right_index = 0_usize;
    let mut shared = 0_usize;
    while left_index < left_count && right_index < right_count {
        match left[left_index].token.cmp(&right[right_index].token) {
            Ordering::Less => left_index += 1,
            Ordering::Greater => right_index += 1,
            Ordering::Equal => {
                shared += 1;
                left_index += 1;
                right_index += 1;
            }
        }
    }
    let union = left_count + right_count - shared;
    if union == 0 {
        return false;
    }
    let jaccard = shared as f64 / union as f64;
    let retained = (2.0 * jaccard) / (1.0 + jaccard);
    if retained <= 0.0 {
        return false;
    }
    let ani = 1.0 - (-retained.ln() / K as f64);
    ani >= MASH_ANI_NUMERATOR as f64 / MASH_ANI_DENOMINATOR as f64
}

fn collinear_coverage(
    left: &[SketchEntry],
    right: &[SketchEntry],
    left_length: usize,
    right_length: usize,
) -> bool {
    let mut forward = Vec::new();
    let mut reverse = Vec::new();
    let mut left_index = 0_usize;
    let mut right_index = 0_usize;
    while left_index < left.len() && right_index < right.len() {
        match left[left_index].token.cmp(&right[right_index].token) {
            Ordering::Less => left_index += 1,
            Ordering::Greater => right_index += 1,
            Ordering::Equal => {
                let a = left[left_index];
                let b = right[right_index];
                if a.occurrences == 1 && b.occurrences == 1 {
                    if a.forward_canonical == b.forward_canonical {
                        forward.push((a.position, b.position));
                    } else {
                        reverse.push((a.position, right_length - (b.position + K)));
                    }
                }
                left_index += 1;
                right_index += 1;
            }
        }
    }
    ordered_once_with_coverage(&mut forward, left_length, right_length)
        || ordered_once_with_coverage(&mut reverse, left_length, right_length)
}

fn ordered_once_with_coverage(
    anchors: &mut [(usize, usize)],
    left_length: usize,
    right_length: usize,
) -> bool {
    if anchors.len() < MIN_POSITION_BINS {
        return false;
    }
    let mut left_bins = [false; POSITION_BINS];
    let mut right_bins = [false; POSITION_BINS];
    for &(left, right) in anchors.iter() {
        left_bins
            [(left.saturating_mul(POSITION_BINS) / left_length.max(1)).min(POSITION_BINS - 1)] =
            true;
        right_bins
            [(right.saturating_mul(POSITION_BINS) / right_length.max(1)).min(POSITION_BINS - 1)] =
            true;
    }
    if left_bins.into_iter().filter(|covered| *covered).count() < MIN_POSITION_BINS
        || right_bins.into_iter().filter(|covered| *covered).count() < MIN_POSITION_BINS
    {
        return false;
    }
    anchors.sort_unstable();
    let breaks = anchors
        .windows(2)
        .filter(|pair| pair[1].1 <= pair[0].1)
        .count();
    breaks <= 1
}

fn confirm_stars(
    classes: &[ExactClass],
    adjacency: &[Vec<usize>],
    threads: usize,
) -> Result<Vec<DetectionGroup>, String> {
    let ranks: Vec<_> = classes.iter().map(class_rank).collect();
    let mut remaining = vec![true; classes.len()];
    let mut remaining_count = classes.len();
    let mut groups = Vec::new();
    while remaining_count != 0 {
        let center = choose_center(&remaining, adjacency, &ranks)
            .ok_or_else(|| "Inconsistent remaining detection-group state".to_string())?;
        let mut member_classes = vec![center];
        let candidates: Vec<_> = adjacency[center]
            .iter()
            .copied()
            .filter(|&index| remaining[index])
            .collect();
        let confirmed = confirm_members(&classes[center].sequence, classes, &candidates, threads)?;
        member_classes.extend(confirmed);
        member_classes.sort_unstable_by(|left, right| ranks[*left].cmp(&ranks[*right]));
        member_classes.dedup();
        let representative_class = center;
        let mut members = Vec::new();
        for class_index in &member_classes {
            if remaining[*class_index] {
                remaining[*class_index] = false;
                remaining_count -= 1;
            }
            members.extend_from_slice(&classes[*class_index].source_indices);
        }
        members.sort_unstable();
        let representative = classes[representative_class]
            .source_indices
            .first()
            .copied()
            .ok_or_else(|| "Detection group has no representative sequence".to_string())?;
        groups.push((
            representative_class,
            DetectionGroup {
                representative,
                members,
            },
        ));
    }
    groups.sort_unstable_by(|(left, _), (right, _)| ranks[*left].cmp(&ranks[*right]));
    Ok(groups.into_iter().map(|(_, group)| group).collect())
}

fn class_rank(class: &ExactClass) -> ClassRank<'_> {
    let invalid = class
        .sequence
        .iter()
        .filter(|&&base| dna_bits(base).is_none())
        .count();
    (
        invalid,
        std::cmp::Reverse(class.sequence.len()),
        *blake3::hash(&class.sequence).as_bytes(),
        class.sequence.as_slice(),
    )
}

fn choose_center(
    remaining: &[bool],
    adjacency: &[Vec<usize>],
    ranks: &[ClassRank<'_>],
) -> Option<usize> {
    (0..remaining.len())
        .filter(|&index| remaining[index])
        .max_by(|&left, &right| {
            let left_degree = adjacency[left].iter().filter(|&&n| remaining[n]).count();
            let right_degree = adjacency[right].iter().filter(|&&n| remaining[n]).count();
            left_degree
                .cmp(&right_degree)
                .then_with(|| ranks[right].cmp(&ranks[left]))
        })
}

fn confirm_members(
    center: &[u8],
    classes: &[ExactClass],
    candidates: &[usize],
    threads: usize,
) -> Result<Vec<usize>, String> {
    if candidates.is_empty() {
        return Ok(Vec::new());
    }
    let mut builder = Aligner::builder().asm5().with_cigar();
    builder.mapopt.best_n = 8;
    let worker_count = threads.max(1).min(candidates.len()).min(16);
    builder.mapopt.cap_kalloc = (800_000_000_i64 / worker_count as i64).max(1);
    let aligner = builder
        .with_seq_and_id(center, b"group-center")
        .map_err(|error| format!("Failed to build minimap2 detection-group index: {error}"))?;
    if worker_count == 1 {
        let mut accepted = Vec::new();
        for &candidate in candidates {
            if strict_confirmation(&aligner, center, &classes[candidate].sequence)? {
                accepted.push(candidate);
            }
        }
        return Ok(accepted);
    }
    let next = AtomicUsize::new(0);
    let (sender, receiver) = sync_channel::<Result<(usize, bool), String>>(worker_count * 2);
    thread::scope(|scope| {
        for _ in 0..worker_count {
            let sender = sender.clone();
            let next = &next;
            let aligner = &aligner;
            scope.spawn(move || loop {
                let position = next.fetch_add(1, AtomicOrdering::Relaxed);
                let Some(&candidate) = candidates.get(position) else {
                    break;
                };
                let result = strict_confirmation(aligner, center, &classes[candidate].sequence)
                    .map(|matches| (candidate, matches));
                if sender.send(result).is_err() {
                    break;
                }
            });
        }
        drop(sender);
        let mut accepted = Vec::new();
        for result in receiver {
            let (candidate, matches) = result?;
            if matches {
                accepted.push(candidate);
            }
        }
        accepted.sort_unstable();
        Ok(accepted)
    })
}

fn strict_confirmation(
    aligner: &Aligner<Built>,
    center: &[u8],
    member: &[u8],
) -> Result<bool, String> {
    if center.is_empty() || member.is_empty() {
        return Ok(false);
    }
    let mappings = aligner
        .map(member, false, false, None, None, Some(b"group-member"))
        .map_err(|error| format!("Failed to confirm minimap2 detection group: {error}"))?;
    Ok(mappings
        .iter()
        .any(|mapping| linear_match(mapping, center.len(), member.len()))
        || circular_two_segment_match(&mappings, center.len(), member.len()))
}

fn linear_match(mapping: &Mapping, target_length: usize, query_length: usize) -> bool {
    mapping.block_len > 0
        && mapping.match_len.max(0) as usize * ALIGN_IDENTITY_DENOMINATOR
            >= mapping.block_len as usize * ALIGN_IDENTITY_NUMERATOR
        && span(mapping.query_start, mapping.query_end) * COVERAGE_DENOMINATOR
            >= query_length * COVERAGE_NUMERATOR
        && span(mapping.target_start, mapping.target_end) * COVERAGE_DENOMINATOR
            >= target_length * COVERAGE_NUMERATOR
}

fn circular_two_segment_match(
    mappings: &[Mapping],
    target_length: usize,
    query_length: usize,
) -> bool {
    for (first_index, first) in mappings.iter().enumerate() {
        for second in mappings.iter().skip(first_index + 1) {
            let (first, second) = if first.query_start <= second.query_start {
                (first, second)
            } else {
                (second, first)
            };
            if first.strand != second.strand
                || first.target_id != second.target_id
                || span(first.query_start, first.query_end) == 0
                || span(second.query_start, second.query_end) == 0
                || first.query_end > second.query_start
                || overlaps(
                    first.target_start,
                    first.target_end,
                    second.target_start,
                    second.target_end,
                )
            {
                continue;
            }
            let wraps_once = match first.strand {
                Strand::Forward => first.target_start > second.target_start,
                Strand::Reverse => first.target_start < second.target_start,
            };
            if !wraps_once {
                continue;
            }
            let matches = first.match_len.max(0) as usize + second.match_len.max(0) as usize;
            let blocks = first.block_len.max(0) as usize + second.block_len.max(0) as usize;
            if blocks != 0
                && matches * ALIGN_IDENTITY_DENOMINATOR >= blocks * ALIGN_IDENTITY_NUMERATOR
                && (span(first.query_start, first.query_end)
                    + span(second.query_start, second.query_end))
                    * COVERAGE_DENOMINATOR
                    >= query_length * COVERAGE_NUMERATOR
                && (span(first.target_start, first.target_end)
                    + span(second.target_start, second.target_end))
                    * COVERAGE_DENOMINATOR
                    >= target_length * COVERAGE_NUMERATOR
            {
                return true;
            }
        }
    }
    false
}

fn span(start: i32, end: i32) -> usize {
    end.saturating_sub(start).max(0) as usize
}

fn overlaps(a_start: i32, a_end: i32, b_start: i32, b_end: i32) -> bool {
    a_start < b_end && b_start < a_end
}

#[cfg(test)]
mod tests {
    use super::*;

    fn token(code: u64) -> Token {
        Token {
            hash: splitmix64(code),
            code,
        }
    }

    fn random_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = splitmix64(state);
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn mutate_substitutions(sequence: &[u8], count: usize) -> Vec<u8> {
        let mut result = sequence.to_vec();
        let mut selected = vec![false; result.len()];
        let mut state = 0x8d26_1c47_b1e9_53f1_u64;
        let mut changed = 0_usize;
        while changed < count.min(result.len()) {
            state = splitmix64(state);
            let position = state as usize % result.len();
            if selected[position] {
                continue;
            }
            selected[position] = true;
            result[position] = match result[position] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
            changed += 1;
        }
        result
    }

    fn groups(records: Vec<(&str, Vec<u8>)>) -> Vec<DetectionGroup> {
        let records: Vec<_> = records
            .into_iter()
            .map(|(name, sequence)| (name.to_string(), sequence))
            .collect();
        construct_detection_groups(&records, 2).unwrap()
    }

    fn brute_pairs(samples: &[Vec<Token>]) -> Vec<(usize, usize)> {
        let mut canonical = samples.to_vec();
        for sample in &mut canonical {
            sample.sort_unstable();
            sample.dedup();
        }
        let mut pairs = Vec::new();
        for left in 0..canonical.len() {
            for right in left + 1..canonical.len() {
                if sampled_jaccard_at_least(&canonical[left], &canonical[right]) {
                    pairs.push((left, right));
                }
            }
        }
        pairs
    }

    #[test]
    fn ppjoin_matches_bruteforce_at_boundaries_df_ties_and_permutations() {
        let samples = vec![
            vec![token(1), token(2), token(3), token(4)],
            vec![token(1), token(5), token(6), token(7), token(8)], // 1 / 8: accept
            vec![token(2), token(9), token(10), token(11)],
            vec![token(3), token(12), token(13), token(14)],
            vec![token(4), token(15), token(16), token(17)],
        ];
        assert_eq!(ppjoin_pairs(samples.clone()), brute_pairs(&samples));
        let permutation = [3_usize, 1, 4, 0, 2];
        let permuted: Vec<_> = permutation
            .iter()
            .map(|&index| samples[index].clone())
            .collect();
        assert_eq!(ppjoin_pairs(permuted.clone()), brute_pairs(&permuted));
    }

    #[test]
    fn grouping_is_identical_inline_and_parallel() {
        let original = random_dna(3_000, 101);
        let records = vec![
            ("original".to_string(), original.clone()),
            (
                "near".to_string(),
                mutate_substitutions(&original, original.len() * 2 / 100),
            ),
            ("duplicate".to_string(), original),
            ("unrelated".to_string(), random_dna(3_000, 202)),
        ];
        let inline = construct_detection_groups(&records, 0).unwrap();
        let parallel = construct_detection_groups(&records, 4).unwrap();
        assert_eq!(parallel, inline);
    }

    #[test]
    fn exact_duplicates_stay_together() {
        let sequence = random_dna(3_000, 1);
        let result = groups(vec![
            ("one", sequence.clone()),
            ("two", sequence),
            ("other", random_dna(3_000, 2)),
        ]);
        assert_eq!(
            result.iter().map(|group| group.members.len()).max(),
            Some(2)
        );
        let mut members: Vec<_> = result.into_iter().flat_map(|group| group.members).collect();
        members.sort_unstable();
        assert_eq!(members, vec![0, 1, 2]);
    }

    #[test]
    fn accepts_97_percent_substitutions() {
        let sequence = random_dna(3_000, 3);
        let member = mutate_substitutions(&sequence, 90);
        let result = groups(vec![("a", sequence.clone()), ("b", member)]);
        assert_eq!(result.len(), 1);
    }

    #[test]
    fn accepts_98_percent_with_indel() {
        let sequence = random_dna(3_000, 4);
        let mut member = sequence.clone();
        member.drain(1_000..1_030);
        member = mutate_substitutions(&member, 30);
        assert_eq!(groups(vec![("a", sequence), ("b", member)]).len(), 1);
    }

    #[test]
    fn accepts_reverse_complement() {
        let sequence = random_dna(3_000, 5);
        assert_eq!(
            groups(vec![
                ("a", sequence.clone()),
                ("b", crate::prescreen::reverse_complement(&sequence))
            ])
            .len(),
            1
        );
    }

    #[test]
    fn accepts_one_circular_rotation() {
        let sequence = random_dna(3_000, 6);
        let mut rotated = sequence[1_000..].to_vec();
        rotated.extend_from_slice(&sequence[..1_000]);
        assert_eq!(groups(vec![("a", sequence), ("b", rotated)]).len(), 1);
    }

    #[test]
    fn rejects_80_percent_local_only_match() {
        let sequence = random_dna(3_000, 7);
        let mut member = sequence[..2_400].to_vec();
        member.extend_from_slice(&random_dna(600, 8));
        assert_eq!(groups(vec![("a", sequence), ("b", member)]).len(), 2);
    }

    #[test]
    fn accepts_explicit_four_percent_nonhomologous_tail() {
        let sequence = random_dna(3_000, 9);
        let mut member = sequence[..2_880].to_vec();
        member.extend_from_slice(&random_dna(120, 10));
        assert_eq!(groups(vec![("a", sequence), ("b", member)]).len(), 1);
    }

    #[test]
    fn rejects_three_block_rearrangement() {
        let sequence = random_dna(3_000, 11);
        let mut member = sequence[1_000..2_000].to_vec();
        member.extend_from_slice(&sequence[..1_000]);
        member.extend_from_slice(&sequence[2_000..]);
        assert_eq!(groups(vec![("a", sequence), ("b", member)]).len(), 2);
    }
}
