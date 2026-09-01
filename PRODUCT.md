# Product

<!-- impeccable:product-schema 1 -->

## Platform

web

## Users

Bioinformatics analysts reviewing per-sample viral candidate evidence from sequencing runs.

## Product Purpose

viroflash detects and reports viral candidate hypotheses from FASTQ sequencing data. Success means that analysts can inspect evidence, uncertainty, quality observations, provenance, and integration evidence without mistaking candidate-level research gates for sample-level or clinical conclusions.

## Positioning

The report preserves viroflash's discovery/validation split, set-valued reference-group hypotheses, synthetic-decoy background model, and fixed-family statistical evidence as auditable first-class concepts.

## Operating Context

The Rust CLI runs offline and emits per-run reports. JSON and TSV are stable machine-readable outputs. The human-readable layer consists of a self-contained interactive HTML report and a companion CSV containing core candidate fields.

## Capabilities and Constraints

- Report language is English.
- HTML must work offline without CDN, remote fonts, or a frontend runtime.
- JSON and TSV schemas, deterministic ordering, sampling semantics, candidate decisions, and scientific limitations remain unchanged.
- A reference-group OR hypothesis must not be expanded into individual positive calls.
- Integration evidence remains separate from general viral candidate detection.
- Empty candidates means no candidates were reported, not `NOT_DETECTED`.
- The report must not imply validated QC, classical FDR control, clinical confidence, or diagnosis.

## Evidence on Hand

The implementation provides run metadata, sampling and audit observations, index provenance, thresholds, complete candidate evidence, decision reasons, statistical background details, and integration sites. No externally validated clinical thresholds, sensitivity/specificity claims, institutional branding, or taxonomy lookup are available and none may be fabricated.

## Product Principles

- Put interpretation boundaries next to the values they constrain.
- Support rapid triage first and auditable evidence review second.
- Prefer direct labels and explicit units over unexplained acronyms.
- Preserve unresolved hypotheses and uncertainty instead of forcing specificity.
- Keep human and machine outputs traceable to the same candidate records.

## Accessibility & Inclusion

The report should support keyboard navigation, reduced motion, high-contrast text, color-independent status encoding, responsive layouts, and printable review records.
