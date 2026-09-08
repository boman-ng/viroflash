# Phase 0 analysis contract

This document freezes Phase 0 evidence and the explicitly authorized scientific decisions. It does not alter the production path, change a production schema, or authorize inspecting a v0.5 result.

## Frozen AnalysisProfile

The user explicitly authorized δ=1e-5, β=.05, and α=.05 on 2026-09-08. `evaluation/phase0/analysis-profile.json` records the actor as `user`, the date, and exactly those three parameters with decision status `FROZEN`; no person or approval system is invented. The 68 internal and 59 external runs were not used to select or revise them.

| Field | Symbol | Value | Product decision |
|---|---:|---:|---|
| `minimum_relevant_fraction` | δ | 0.00001 | Design for a group at 10 fragments per million input fragments; not a clinical LoD. |
| `familywise_interval_error` | α | 0.05 | Require 95% simultaneous target-family interval coverage. |
| `familywise_miss_probability` | β | 0.05 | Permit at most 5% planned familywise probability of sampling no fragment at δ. |

There is no minimum-read, MAPQ, edit-distance, alignment-margin, coverage, window, p/q, or PASS cutoff. `k=21`, minimap2 `sr`, all-chain enumeration, and ten equal-width diagnostic windows are profile inputs, not CLI choices. Changing any profile byte creates a new profile digest, requires new indexes, and restarts all 127 evaluations.

The values are authorized product/scientific choices, not estimates or calibration results from these datasets. α is the familywise error for simultaneous target-family intervals. Horvitz and Thompson support explicit inclusion-probability accounting; exact hypergeometric bounds support finite-population intervals; Dunn supports Bonferroni simultaneous coverage. These references support the methods, not empirical calibration of δ, β, or α: <https://doi.org/10.1080/01621459.1952.10483446>, <https://doi.org/10.1007/978-1-4612-3140-0>, <https://doi.org/10.1080/01621459.1961.10482090>.

The point estimate and `interval_lower`, `interval_upper`, `interval_level`, and `interval_method` are reported continuously and unchanged. The interval is the only precision expression. No interval width produces a binary category, and there is no “trusted”, “gray-zone”, or precision-met mapping.

## Digest and analysis identity

All digests are lowercase SHA-256. A file digest covers exact file bytes. `profile_digest` covers the exact committed bytes of `analysis-profile.json`; JSON semantic reserialization is therefore a profile change.

Composite identities use this unambiguous framing: the ASCII domain followed by NUL, then each ordered item as `u16_be(tag_byte_length) || tag_utf8 || 32_digest_bytes`. Paths, sample labels, output directories, thread counts, jobs, and timestamps are excluded.

- `input_digest`: domain `viroflash-input-v1`; ordered items are `mode`, then `r1` and, for PE, `r2`. The mode item is the SHA-256 of literal `SE` or `PE`; each read item is the exact compressed-file digest. PE order is R1 then R2.
- `reference_set_digest`: domain `viroflash-reference-set-v1`; ordered items are `host_fasta`, `target_fasta`, and `reference_group_ledger`.
All input-file digests are frozen, but the v0.5 binary, index artifacts, and their final artifact contract do not exist. Phase 0 therefore does not define or claim a computable index digest, artifact binding, or complete analysis identity.

## ReferenceGroup ledger

`evaluation/phase0/reference-group-ledger-contract.tsv` is the executable ledger contract. `evaluation/phase0/reference-groups.json` freezes actual ledgers for the internal DNA panel, external respiratory panel, and external HPV panel. No historical 97%-identity or 95%-coverage group is imported.

Every target FASTA record appears in exactly one group. The record ID is the first whitespace-delimited header token and must be non-empty UTF-8 and unique within the FASTA. Sequence lines are joined after removing ASCII whitespace and uppercasing. Accepted DNA symbols are exactly `ACGTMRWSYKVHDBN`; complement pairs are `A-T`, `C-G`, `M-K`, `R-Y`, `W-W`, `S-S`, `V-B`, `H-D`, and `N-N`. Any other character rejects the file; it is never silently removed or replaced. This follows the NCBI FASTA requirement for a unique SeqID and IUPAC symbols and the INSDC standard nucleotide-code table: <https://www.ncbi.nlm.nih.gov/genbank/fastaformat>, <https://www.insdc.org/submitting-standards/feature-table/>.

Two records may share a group only when their complete normalized sequence is byte-identical to the other's complete sequence or reverse complement. Length must match. Circular rotations, high identity, partial coverage, shared taxonomy, names, and runtime reads never merge groups. `target_group_id` is `sha256:` plus the SHA-256 of the lexicographically smaller byte string of the normalized sequence and its reverse complement. The representative is the lexicographically smallest member ID; membership is OR-only and never supports a member-level claim.

The canonical ledger is UTF-8 TSV with one row per FASTA record, the contract header order excluding the manifest-owned `ledger_digest`, LF endings, groups ordered by `target_group_id`, members ordered bytewise by ID, and a final LF. The index manifest's `ledger_digest` covers those exact ledger bytes; it is not embedded in the hashed file. The frozen counts are 21,790 records/20,560 groups for the internal panel and 2 records/2 groups for each external panel. Runtime FASTQ cannot create, split, merge, or reorder groups; distinct group count fixes the multiplicity denominator.

## Competitive adjudication

The candidate boundary is every positive-query-length CIGAR mapping returned by minimap2 `sr` with `MM_F_ALL_CHAINS`; Viroflash adds no acceptance cutoff. This is profile-relative evidence because minimap2's preset has its own documented chaining and alignment settings.

For each read end and role/group, compare `alignment_score / query_length` exactly by integer cross multiplication. Higher is better. MAPQ, NM, primary/supplementary flags, aligned length, contig order, and IDs are audit metadata, never tie breakers. Equality remains a tie; there is no epsilon or score margin.

Fragment adjudication is executable in this order:

1. Build the equal-best host and target-group sets independently for each end. An unmapped end is neutral.
2. If an end has host equal to or better than its best target, that end is host-confounded. A fragment with both target evidence and any host-confounded end is `CONFOUNDED_WITH_HOST` and contributes no target numerator.
3. For non-confounded PE, intersect the non-empty target sets from R1 and R2. Disjoint sets are unresolved across their union; an unmapped end does not erase the mapped end.
4. One surviving group contributes one supporting fragment. A singleton group is `RESOLVED_TO_REFERENCE_GROUP`; a multi-member exact-equivalence group is `AMBIGUOUS_WITHIN_GROUP`.
5. Multiple surviving groups are `UNRESOLVED_ACROSS_GROUPS`: no group receives numerator support, the RUN `unassigned_fragments` count increases once, and each involved group receives one `cross_group_ambiguous_fragments` diagnostic incidence and an indeterminate row.

All groups use the same selected-fragment denominator. Cross-group incidence may sum above the number of unresolved fragments and is never a second denominator or fractional allocation.

Distribution is fragment-owned and representative-coordinate based. `covered_bases` is the union of half-open intervals from supporting fragment chains; overlap, split chains, and two PE ends do not double count bases. `coverage_fraction` is that union size divided by representative length. `occupied_windows` counts the ten equal-width representative bins touched by the union. `split_events` and `discordant_fragments` count supporting fragments at most once per group. These fields are diagnostics only and create no cutoff or confidence score.

The minimap2 manual defines `sr`, all-chain behavior, pairing, and preset-internal filters; the SAM specification defines paired, secondary, and supplementary records; Salmon's equivalence-class model supports preserving ambiguous sets without forcing allocation. Viroflash adopts only those candidate and ambiguity semantics, not RNA abundance estimation: <https://github.com/lh3/minimap2/blob/master/minimap2.1>, <https://samtools.github.io/hts-specs/SAMv1.pdf>, <https://pmc.ncbi.nlm.nih.gov/articles/PMC5600148/>.

## Phase ownership switch

Phase 1 introduces one internal `EvidenceReport` owner for fragment evidence while leaving current serialization untouched. It must not adapt new fragment fields into the old read-end model, shadow-write, alias, or maintain two accumulators. Phase 4 deletes the old report owner and atomically switches the sample directory to the new three-file output. A failed run may leave only `perf.json` with error status; successful `report.csv` and `report.html` appear together from the same immutable `EvidenceReport`.

## CLI and outputs

The complete normal CLI is:

```text
viroflash index --host-fa HOST.fa --target-fa TARGET.fa --out INDEX_DIR [--threads N]
viroflash run --r1 SAMPLE_R1.fastq.gz [--r2 SAMPLE_R2.fastq.gz] --index INDEX_DIR --out SAMPLE_REPORT_DIR [--threads N]
```

There is no automatic index build, scientific parameter option, legacy alias, converter, symlink, or dual write. A successful sample directory contains exactly `report.csv`, `report.html`, and `perf.json`. `evaluation/phase0/output-fields.tsv` freezes every field/section and its record applicability. CSV and HTML serialize the same `EvidenceReport`; HTML does not recompute it, and `perf.json` contains no biological result.

## Evaluation and historical baseline

`evaluation/phase0/evaluation-manifest.json` freezes 68 internal and 59 external run IDs, cohorts, SE/PE layouts, absolute input/reference paths, byte sizes, truth labels, provenance, checksums, known label/history discrepancies, and planned resources. It is evaluation metadata, not production configuration.

The external 82 FASTQ digests correspond exactly to `metadata/FASTQ_SHA256SUMS`; external reference digests correspond to `refs/SHA256SUMS`. Internal truth and target-reference digests correspond to the existing `run_config.json`. The 136 internal FASTQ digests were computed on 2026-09-08 by parallel streaming SHA-256 over the existing gzip bytes without decompression and are independently frozen in `internal-fastq-sha256.tsv` with status `PHASE0_FROZEN_COMPRESSED_SHA256`; this is explicit Phase 0 provenance, not an upstream authoritative source. The internal host digest remains `PHASE0_COMPUTED` metadata.

`PRE_EXECUTION_FROZEN` means sample, truth, checksum, ledger, and the authorized profile bytes are fixed before execution. It does not mean any binary/index artifact is bound. Phase 6 owns artifact validation: the final E2E runner must stream SHA-256 from the actual binary and cohort-index paths under the then-current concrete artifact contract before running results. Git review, not this verifier, controls any attempted rebinding. No Phase 0 artifact gate, registry, migration, alternate schema, or binding claim exists.

`evaluation/phase0/historical/` is an immutable snapshot of old-system observations. It is comparison evidence only: not v0.5 acceptance truth, not a pristine holdout, not a scorer input, and not permission to tune δ/α/β.

## Verification

Run `python3 evaluation/phase0/verify.py`. It directly checks all 68 internal runs against `panel_truth.json`, all 59 external runs against `label_comparison.tsv`, run-to-FASTQ names, the independent 136-file checksum ledger, the exact frozen three-parameter profile, and every target FASTA record against the three committed ReferenceGroup ledgers. It also checks output fields, historical snapshots, and Cargo version `0.3.0`. Its success means checksum metadata consistency; it does not read or recompute FASTQ bytes or validate future execution artifacts.

Phase 6 Input Pass 1 owns streaming recomputation of the 136 internal FASTQ digests while reading the exact compressed inputs. Any mismatch must stop that run before analysis. No separate default Phase 0 byte scan is implied.
