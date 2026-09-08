# Viroflash Evidence Report Design

## Purpose

The HTML report is a read-only rendering of the same `EvidenceReport` model written to
`report.csv`. It presents measured fragment evidence and explicit limitations. It does not
produce a clinical conclusion, quantitation status, confidence tier, or acceptance decision.

## Source Of Truth

- `EvidenceReport` is the sole biological report model.
- CSV serializes the frozen field order in `evaluation/phase0/output-fields.tsv`.
- HTML renders the supplied model without recomputing fractions, intervals, coverage, windows,
  ambiguity, host conflict, or integration diagnostics.
- `perf.json` contains execution telemetry only.

## Information Order

1. **Interpretation boundary** states what the attributed fragment fraction measures and excludes.
2. **Run integrity** shows status, reason codes, counts, sampling probability, and all digests.
3. **Observed target signals** shows every emitted ReferenceGroup row and its complete evidence.
4. **Evidence detail** explains ambiguity, host conflict, coverage, windows, and integration fields.
5. **Methods and limitations** describes the two FASTQ passes, target k-mer workload gate,
   competitive HOST+TARGET alignment, and finite-population interval.

## Visual Rules

- Use a flat, high-contrast ledger with tabular values and explicit text labels.
- Keep uncertainty and interpretation boundaries adjacent to the values they qualify.
- Never encode scientific meaning by color alone.
- Preserve readable responsive and print layouts without hiding report fields.
- Diagnostics remain descriptive; visual prominence must not turn them into decision rules.

## Analysis States

- `CONFORMANT_COMPLETE` means every selected fragment was evaluable by the target 21-mer gate.
- `CONFORMANT_WITH_LIMITATIONS` names the unevaluable selected-fragment count in `reason_codes`.
- Target rows independently report evidence and attribution states; no run or target state is a
  positive/negative verdict.

## Output Boundary

A successful output directory contains exactly `report.csv`, `report.html`, and `perf.json`.
Failures contain only an error-shaped `perf.json`. Report files are written atomically before the
directory is finalized.
