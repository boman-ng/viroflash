# Analysis methods

## References and attribution

HOST is the background competitor. TARGET sequences share a reference group when their complete normalized IUPAC sequences are identical or reverse complements. The smallest member ID represents each group; FASTA descriptions are retained. Indexes contain the composite reference, minimap2 index, target Bloom filter and reference metadata.

The target-only 21-mer Bloom filter selects work. Passing fragments have an unmasked target hit and satisfy minimap2 short-read chaining hit-count and union-coverage requirements. Low-complexity masking uses symmetric DUST. Both ends compete against HOST and TARGET with minimap2 `sr`.

Alignment scores are compared as score/query length using integer cross multiplication. Host ties or advantages give no target support. Paired target sets are intersected; an unmapped end is neutral. Unresolved cross-group ties remain unassigned. One SE read or PE pair is one fragment and supports at most one reference group.

Coverage is the union of 0-based half-open reference intervals. Ten equal reference windows summarize spatial occupancy. Split and discordant counts describe alignment patterns.

## Sampling

Both modes make one validated pass over plain or gzip FASTQ and retain a bounded sample in memory. Screen samples the original input before Bloom. Full applies Bloom to every fragment and samples the resulting candidates. Final selected fragments are analyzed in bounded parallel batches.

Let N be the original fragment count, C the candidate count, m the number of reference groups and β=0.05. Precision sets δ to 1e-4 (fast), 1e-5 (standard) or 1e-6 (sensitive). The common capacity is obtained by inverting the binomial zero-hit probability:

```text
n0 = ceil(log(β/m) / log(1-δ))
```

The implementation rounds conservatively and verifies the miss-probability bound. A bottom-k reservoir selects min(N,n0) input fragments in screen or min(C,n0) candidates in full. Full uses δ as a conservative lower bound on the candidate fraction corresponding to a whole-library target δ. It retains at least the Bloom-passing screen sample under matching input, index and precision.

BLAKE3 priorities use profile/index identity, normalized fragment ID and original ordinal. The ordinal also breaks hash ties. Selection is independent of worker completion order and FASTQ compression. Under the ideal uniform-priority model conditioned on distinct keys, bottom-k is a simple random sample without replacement. The zero-hit bound and a union bound over the reference family control sampling loss for sufficiently abundant detectable signals. The design protects relative abundance, not a fixed absolute copy count.

For each input end, decoding and digesting overlap FASTQ parsing through bounded byte chunks. Readers prefetch bounded record batches. Workers calculate keys and, in full mode, evaluate Bloom. A single global reservoir retains selected sequences without quality strings. Screen applies Bloom to its finalized sample. Both modes release Bloom before loading the shared mapping index and transfer selected batches without copying sequences. No candidate files are written.

## Abundance and target score

Let P be N for screen or C for full, n the actual selected count and x a group's attributed support. The estimated original-library fraction is:

```text
p = (P/N) * (x/n)
support_pct = 100 * p
```

Screen's n includes selected Bloom negatives. Exact equal-tailed hypergeometric inversion uses P,n,x with per-group error 0.05/m. Integer population-count endpoints are divided by N, converted to percentages and rounded outward to 1e-10 percentage units. Intervals have simultaneous 95% sampling coverage under the sampling model; they do not include prescreen or classification error.

At the whole-library target δ, expected selected support is μ=n*δ*N/P. The reported score is:

```text
sampling_target_score = (x-μ)/(x+μ) = (p-δ)/(p+δ)
```

The score is negative below the target, zero at the target and positive above it. Its range is [-1,1); it measures relative abundance, not confidence. The baseline uses the actual selected count. Sample-only rows have an empty score and interval.

## Outputs

Successful runs write `report.csv`, `report.html` and `perf.json` atomically. CSV has 24 ordered fields with BOM and CRLF, sorted by support count then reference ID. HTML displays the top 20 supported groups and embeds the complete CSV. Percentages use original input fragments; target share uses all attributed target fragments.

Run information records the mode, precision, sample capacity, population size, actual selection fraction and expected support at the target. Performance records index loading, scan/sample, selected-analysis and report-writing times, CPU time, peak memory, storage I/O and peak reservoir size. Scan/sample includes input, sampling and Bloom in both modes; selected-analysis includes mapper loading and competitive alignment.

The embedded base profile bytes define index identity. Runtime mode, precision and sample population are recorded separately. Production and tests use Rust and the existing dependencies.
