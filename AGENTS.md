# Repository Guidelines

## Project Structure

`viroflash` is a Rust 2021 command-line application. `src/main.rs` owns the strict `index` and
`run` CLI, `src/lib.rs` exposes only supported entry points, and `src/pipeline.rs` orchestrates
analysis. `analysis_profile.rs`, `fastq_input.rs`, and `sampling_design.rs` own the frozen profile
and two FASTQ passes. `reference_group.rs` and `reference_index.rs` build and load reusable
HOST+TARGET indexes. `kmer_gate.rs`, `competitive_alignment.rs`, `evidence.rs`, and
`integration_evidence.rs` own fragment evidence. `report.rs` writes `report.csv` and visible
`report.html`; `performance_report.rs` writes telemetry-only `perf.json`.

Integration coverage is in `tests/smoke.rs`; focused unit tests live beside their modules. Never
commit local datasets, indexes, reports, `.local/`, `target/`, `.tmp/`, or `*.work/` directories.

## Build and Test

Use the committed lockfile:

```bash
cargo fmt --all -- --check
cargo check --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
cargo test --release --locked
cargo build --release --locked
python3 evaluation/phase0/verify.py
```

Run the narrowest relevant test first, such as `cargo test sampling_design::tests --locked`, then
broaden checks for shared pipeline changes.

## Coding Conventions

Follow default `rustfmt`: four-space indentation, `snake_case` functions and modules, and
`CamelCase` types. Stable scientific constants belong to `analysis_profile.rs`; do not create
unsupported runtime choices. Production input paths return `Result<T, String>` and must not panic
or hide errors. Prefer the standard library and existing dependencies. Preserve deterministic
fragment selection, report column order, atomic output writes, and 0-based half-open coordinates.

## Testing Contract

Add unit tests near changed logic and end-to-end behavior to `tests/smoke.rs`. Cover success and
error paths for CLI or parser changes. Reporting changes require strict CSV column-order checks and
comparison of every visible HTML field against CSV. Concurrency and sampling changes must verify
deterministic output across thread counts. Index changes must exercise reusable `index` followed by
`run --index`; analysis always consumes a reusable index.

Successful runs must produce exactly `report.csv`, `report.html`, and `perf.json`. Failed runs may
leave only an error-shaped `perf.json`.

## Change Control

Use focused Conventional Commits such as `fix(evidence): ...`, `fix(index): ...`, or
`docs: align release guidance`. Review the complete staged diff before committing. Do not push, tag, rewrite
history, modify remotes, or change package/schema versions without explicit authorization.
