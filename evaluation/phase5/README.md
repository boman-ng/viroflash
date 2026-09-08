# Phase 5 Statistical Oracles

Phase 5 adds deterministic counterfactual tests only. It does not tune the frozen profile, add a
runtime cutoff, or use these results as training data.

Run the production-connected exact finite-population oracles. The Python command is only an
orchestrator: test-only Rust boundaries call the production sampling selector and interval function,
then compare their outputs with independent integer/rational enumeration. Coverage uses the frozen
20,560-group family and conditions estimator bias and intervals on each realized sample size.

```bash
python3 evaluation/phase5/statistical_oracles.py
```

Run the gate/exhaustive counterfactual and report invariance checks:

```bash
cargo test pipeline::tests::phase5_gate_counterfactual_quantifies_exhaustive_evidence_difference --locked
cargo test phase5_report_is_invariant_across_threads_reloads_and_gzip_segmentation --locked
cargo test phase5_alignment_retention_is_threads_times_max_fragment_not_total_selected --locked
```

Build the release binary, then record the small-to-large memory curve:

```bash
cargo build --release --locked
python3 evaluation/phase5/memory_curve.py --binary target/release/viroflash
```

The memory command first runs exact instrumentation over 40-, 120-, and 4,096-base reads with
1/2/4/8 workers. It requires retained alignment sequence bytes to remain bounded by configured
threads times the largest current fragment and explicitly verifies that total-selected linear
retention is rejected. It then reports raw RSS measurements and adjacent slopes without a subjective
tolerance. This is not a universal constant-byte bound: maximum fragment length is an unavoidable
term, minimap2 owns internal allocations, retained merged intervals depend on the fixed index, and
Linux `VmHWM` is documented as potentially imprecise.

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
- [Rust `sync_channel` documentation](https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html)
  defines its pending-message buffer as fixed-size. Viroflash narrows the owned sequence bound to
  one current batch of at most one fragment per worker and tests bytes, not only object counts.
- [Rust `f64::next_up` documentation](https://doc.rust-lang.org/std/primitive.f64.html#method.next_up)
  defines the least representable successor. Sampling uses outward-rounded multiplication to bound
  miss probability monotonically, while a separate exact 128-bit-selector oracle verifies the
  production threshold without an epsilon.
