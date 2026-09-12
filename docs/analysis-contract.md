# Analysis methods

Viroflash measures profile-attributed fragment fraction for a fixed family of reference groups. The production profile is embedded from `src/analysis-profile.json`; its exact bytes define the SHA-256 profile identity. It fixes δ = 1e-5, β = 0.05, α = 0.05, target k = 21, minimap2 `sr` with all-chain enumeration, and ten diagnostic windows. The frozen bytes remain the base profile and index identity. Run sampling uses the explicit precision preset and Bloom candidate population.

## Index and sampling

The index combines HOST and TARGET references; HOST is the sole background competitor. Target records share a group only if their complete normalized IUPAC sequences are identical or reverse complements. The representative is the lexicographically smallest member ID. Group membership is fixed before analysis; names and taxonomy do not merge sequences. Representative descriptions are stored in the index, so analysis does not need the original FASTA files.

A single original-input pass validates FASTQ records and paired identifiers, counts all fragments N, hashes decoded content, and applies the Bloom prescreen to every fragment. Readers reject file replacement, size or modification-time changes observed during reading. Passing candidates are written to a run-local binary spool containing original ordinal, normalized ID and read sequences; quality strings are validated but not stored. After EOF fixes the candidate count C and input identity, deterministic BLAKE3 Bernoulli selection operates on that spool. A paired-end fragment shares one selection key and contributes at most one target assignment.

For C candidate fragments and m target groups in the complete index:

```text
M_min = ceil(δC)
π = min(1, 1 - (β/m)^(1/M_min))
```

The presets `fast`, `standard` (default), and `sensitive` set δ to 1e-4, 1e-5, and 1e-6, respectively; β and α remain 0.05. The guarantee concerns retaining at least one candidate supporting a group with M_min or more such fragments, under the ideal-hash sampling model. It does not guarantee a fixed confidence-interval width. The same hash keys give nested selections across presets. The implementation rounds inclusion thresholds outward to preserve the familywise miss budget. Readers and analysis workers exchange bounded batches rather than retaining the input library.

## Prescreen and attribution

A target-only 21-mer Bloom gate selects alignment work. A passing fragment needs a non-SDUST-masked target hit on either end, plus a read end with sufficient hit count and union query coverage for the active minimap2 short-read chaining requirements. Both ends of passing fragments undergo competitive HOST+TARGET alignment. The gate is a workload filter, not a positive call.

Alignment scores are compared as `alignment_score / query_length` using exact integer cross multiplication. Viroflash adds no MAPQ, score-margin, or coverage cutoff.

- Host ties or advantages on a target-evidenced fragment are confounded and contribute no target support.
- Paired target sets are intersected; an unmapped end is neutral. Disjoint sets remain unresolved.
- A single surviving group receives one supporting fragment. Exact-equivalent members remain indistinguishable within that group.
- Cross-group ties remain unresolved; no fractional assignment is made.

Input fragments without an evaluable target k-mer are counted and reported as a limitation; they do not enter the candidate pool. HOST–TARGET split evidence requires disjoint same-end chains with supplementary geometry; split and discordant counts are diagnostic fragment counts, not confirmed integration events.

## Estimates and reports

For n sampled candidates and x uniquely attributed supporting fragments, estimated support in the original library is C*x/n, and abundance is (C/N)*(x/n). Simultaneous intervals use the existing equal-tailed exact hypergeometric inversion on population C, sample n and count x, conditional on the realized sample size, with Bonferroni α/m across the complete target family. Integer candidate-total endpoints L and U become [L/N, U/N], with outward ppm rounding. Sampling all candidates collapses the interval to x/N. The estimand is signal detectable under the current prescreen and attribution rules; intervals cover sampling uncertainty, not gate losses, classification accuracy or biological variation. With no candidates or no selected candidates, the report contains only the sample overview and no fabricated per-virus interval.

Coverage is the union of supporting reference intervals in 0-based half-open coordinates. It measures breadth, not depth; overlapping chains and paired ends do not double-count bases. Occupied windows summarize distribution over the representative sequence.

`report.csv` contains the 23 research fields in `tests/fixtures/output-fields.tsv`, sorted by support descending then reference ID ascending. A sample without attributed targets has one sample-only row. `support_ppm` and its interval estimate detectable support per million original fragments; `support_fragments` remains the observed sampled support, and `selected_fragments` is the sampled candidate count; `target_support_share_pct` uses total target support. HTML embeds its own CSV download and presents the same values plus full evidence and limitations; telemetry is confined to `perf.json`.

The report identity remains `viroflash.evidence-report.v1`. The 23-column layout replaced the former mixed RUN/TARGET_SIGNAL CSV in place, so consumers must validate the actual header. Older indexes without `target_descriptions` require rebuilding.

## Implementation references

The reference Bloom screening approach follows [BioBloom Tools (Chu et al., 2014)](https://doi.org/10.1093/bioinformatics/btu558); its matches remain candidates requiring competitive attribution. Sampling and exact finite-population intervals reuse Viroflash's existing implementations. Input decoding retains [flate2 with the zlib-rs backend](https://github.com/rust-lang/flate2-rs). No digital normalization, adaptive stopping, or additional classifier is introduced.
