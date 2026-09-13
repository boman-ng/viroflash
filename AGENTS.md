# Working on Viroflash

Viroflash is a Rust CLI for panel-based viral screening. Users install a prebuilt executable, build a HOST+TARGET index, run a sample and open its report.

## Code ownership

- `main.rs` defines the CLI; `pipeline.rs` orchestrates analysis.
- `fastq.rs`, `candidates.rs`, `sampling.rs` and `workers.rs` handle input, in-memory bottom-k selection and bounded parallel work.
- `gate/` handles Bloom filtering and symmetric DUST masking; `alignment/` handles competitive mapping; `evidence.rs` aggregates attribution and coverage.
- `index/` owns reusable references and metadata. `profile.rs` and `analysis-profile.json` define the frozen base profile; preserve its bytes and digest.
- `report/` owns CSV and standalone HTML; `telemetry.rs` records performance; `output.rs` owns atomic writes. `tests/fixtures/output-fields.tsv` records output fields.

## Current behavior

- One SE read or PE pair is one fragment. Selection and attribution agree across thread counts and FASTQ compression formats.
- Screen is the default: sample input, Bloom, align. Full: Bloom all input, sample candidates, align. Both use the same precision-derived capacity and one global reservoir, without candidate files.
- Precision fast/standard/sensitive sets 100/10/1 ppm. Full's analyzed fragment set contains screen's for matching input, index and precision.
- Each fragment supports at most one exact reference group. HOST is the sole background competitor. Coordinates are 0-based and half-open.
- Library abundance and intervals use original input fragments; target share uses all attributed target fragments. Keep exact sampling intervals and outward rounding.
- `sampling_target_score` is (estimated library fraction - ppm target)/(estimated library fraction + ppm target).
- Success writes exactly `report.csv`, `report.html` and `perf.json`. CSV has 24 columns, BOM and CRLF. HTML is English, script-free, shows Top 20 and embeds the complete CSV.
- Invalid input or indexes fail explicitly. Keep bounded streaming, shared indexes and existing dependencies.

## Verification

Use the pinned toolchain and lockfile. For Rust or pipeline changes:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
cargo test --release --test smoke --locked
```

Keep distinct checks for input errors, index reuse, attribution, sampling, intervals, report parity and bounded concurrency. Tests are Rust and independent of local datasets. CI runs formatting, Clippy and core tests. Release checks exercise the installed static binary and compare container reports. UI changes need fresh desktop/mobile reports and interaction checks. Documentation changes need link checks and `git diff --check`.

## Local work

`evaluation/`, `.local/`, `.tmp/`, `target/`, `dist/`, `*.work/` and generated images are local artifacts, excluded from Git, default searches and container contexts. Preserve existing results. Use scoped `rg --hidden --no-ignore` to inspect them.

Preserve unrelated changes and review the complete diff before a focused Conventional Commit. Pushes, tags, releases, history/remote changes and version changes require explicit authorization.
