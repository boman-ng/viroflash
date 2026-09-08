# Phase 0 analysis contract

This document freezes only Phase 0 decisions and audit assets. It does not alter the production path, change a production schema, or authorize inspecting a v0.5 result.

## Frozen AnalysisProfile

`evaluation/phase0/analysis-profile.json` was frozen before a v0.5 implementation or result existed. The 68 internal and 59 external runs must not be used to revise it.

| Field | Symbol | Value | Product decision |
|---|---:|---:|---|
| `minimum_relevant_fraction` | δ | 0.00001 | Design for a group at 10 fragments per million input fragments; not a clinical LoD. |
| `familywise_interval_error` | α | 0.05 | Require 95% simultaneous target-family interval coverage. |
| `familywise_miss_probability` | β | 0.05 | Permit at most 5% planned familywise probability of sampling no fragment at δ. |
| `maximum_interval_width` | w | 0.00001 | Call precision met only at interval width no greater than 10 fragments per million. |

There is no minimum-read, MAPQ, edit-distance, alignment-margin, coverage, window, p/q, or PASS cutoff. `k=21`, minimap2 `sr`, all-chain enumeration, and ten equal-width diagnostic windows are profile inputs, not CLI choices. Changing any profile byte creates a new profile digest, requires new indexes, and restarts all 127 evaluations.

## Digest and analysis identity

All digests are lowercase SHA-256. A file digest covers exact file bytes. `profile_digest` covers the exact committed bytes of `analysis-profile.json`; JSON semantic reserialization is therefore a profile change.

Composite identities use this unambiguous framing: the ASCII domain followed by NUL, then each ordered item as `u16_be(tag_byte_length) || tag_utf8 || 32_digest_bytes`. Paths, sample labels, output directories, thread counts, jobs, and timestamps are excluded.

- `input_digest`: domain `viroflash-input-v1`; ordered items are `mode`, then `r1` and, for PE, `r2`. The mode item is the SHA-256 of literal `SE` or `PE`; each read item is the exact compressed-file digest. PE order is R1 then R2.
- `reference_set_digest`: domain `viroflash-reference-set-v1`; ordered items are `host_fasta`, `target_fasta`, and `reference_group_ledger`.
- `index_digest`: domain `viroflash-index-v1`; ordered items are `profile`, `reference_set`, and every required index artifact digest in manifest order.
- `analysis_identity`: domain `viroflash-analysis-v1`; ordered items are `binary`, `profile`, `index`, and `input`.

Missing component digests leave the composite identity unavailable; null is required and no placeholder is hashed. This is the Phase 0 gap for the 136 internal FASTQ files. Phase 1/2 may implement these formulas but may not add a second identity path.

## ReferenceGroup ledger

`evaluation/phase0/reference-group-ledger-contract.tsv` is the executable ledger contract. Actual ledgers are index outputs and cannot be materialized in Phase 0 without reading reference sequence content; no historical 97%-identity group is imported.

Every target FASTA record appears in exactly one group. Two records may share a group only when their complete ASCII-uppercased sequence is byte-identical, allowing complete reverse-complement equality because alignment is strand-independent. Length must match. Circular rotations, high identity, partial coverage, shared taxonomy, names, and runtime reads never merge groups. `target_group_id` is `sha256:` plus the SHA-256 of the lexicographically smaller byte string of the normalized sequence and its reverse complement. The representative is the lexicographically smallest member ID; membership is OR-only and never supports a member-level claim.

The canonical ledger is UTF-8 TSV with the contract header order excluding the manifest-owned `ledger_digest`, LF endings, rows in `group_ordinal` order, member IDs sorted bytewise and joined by `;`, no quoting, and a final LF. The index manifest's `ledger_digest` covers those exact ledger bytes; it is not embedded in the hashed file. Runtime FASTQ cannot create, split, merge, or reorder groups; ledger row count fixes the multiplicity denominator.

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

`evaluation/phase0/evaluation-manifest.json` freezes 68 internal and 59 external run IDs, cohorts, SE/PE layouts, absolute input/reference paths, byte sizes, truth labels, provenance, available checksums, known label/history discrepancies, and planned resources. It is read-only evaluation metadata, not production configuration.

The external 82 FASTQ digests correspond exactly to `metadata/FASTQ_SHA256SUMS`; external reference digests correspond to `refs/SHA256SUMS`. Internal truth and target-reference digests correspond to the existing `run_config.json`. The internal host digest is retained as non-authoritative `PHASE0_COMPUTED` metadata and is not eligible for `reference_set_digest`. No authoritative checksum source exists for the 136 internal FASTQs, so all are null `MISSING_AUTHORITATIVE_CHECKSUM`; this evidence gap is not repaired or disguised in Phase 0. Binary/index digests remain null because v0.5 artifacts do not exist.

`evaluation/phase0/historical/` is an immutable snapshot of old-system observations. It is comparison evidence only: not v0.5 acceptance truth, not a pristine holdout, not a scorer input, and not permission to tune δ/α/β/w.

## Verification

Run `python3 evaluation/phase0/verify.py`. It reads only committed assets and small truth/label/checksum/config metadata, checks 68+59 counts, unique IDs, exact fields, paths, sizes, checksum syntax and authoritative path correspondence, historical counts, the frozen profile digest, and Cargo version `0.3.0`. It has no FASTQ reader or checksum mode.
