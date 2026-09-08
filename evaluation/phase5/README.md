# Phase 5 Statistical Oracles

Phase 5 adds deterministic counterfactual tests only. It does not tune the frozen profile, add a
runtime cutoff, or use these results as training data.

Run the independent finite-population enumeration and deterministic simulation:

```bash
python3 evaluation/phase5/statistical_oracles.py
```

Run the gate/exhaustive counterfactual and report invariance checks:

```bash
cargo test pipeline::tests::phase5_gate_counterfactual_quantifies_exhaustive_evidence_difference --locked
cargo test phase5_report_is_invariant_across_threads_reloads_and_gzip_segmentation --locked
```

Build the release binary, then record the small-to-large memory curve:

```bash
cargo build --release --locked
python3 evaluation/phase5/memory_curve.py --binary target/release/viroflash
```

The memory command reports raw measurements and adjacent slopes without a pass tolerance. Its
structural claim is narrow: the FASTQ reader has a 1 MiB buffer per input end, alignment task and
result channels each hold at most one fragment per configured thread, and no selected-sequence
reservoir exists. This is not a universal byte bound because a single FASTQ record has no byte cap,
minimap2 owns internal allocations, retained merged intervals depend on the fixed index, and Linux
`VmHWM` is documented as potentially imprecise.

Evidence sources used to resolve implementation ambiguities:

- [RFC 1952](https://www.rfc-editor.org/rfc/rfc1952.html) defines a gzip file as a series of
  consecutive members and requires bounded intermediate decompression storage. Report invariance
  therefore compares the same decoded FASTQ stream across plain, one-member gzip, and differently
  segmented multi-member gzip; it does not reorder decoded FASTQ records.
- [minimap2 manual](https://lh3.github.io/minimap2/minimap2.html) documents the `sr` preset and
  retention of all chains. The Bloom test consequently treats exhaustive minimap2 evidence as a
  workload-gate counterfactual, not biological truth, and reports exact target-evidence and
  submitted-fragment differences for the fixture.
- [Exact hypergeometric interval construction](https://www.samplingbook.manitz.org/articles/samplingbook-Sprop.html)
  states the finite-population tail inversion used by the independent rational-arithmetic oracle.
- [Linux `/proc` documentation](https://docs.kernel.org/filesystems/proc.html) defines `VmHWM` as
  peak resident set size and warns that RSS-related values may be imprecise; the curve is therefore
  measurement evidence rather than a proof by itself.
