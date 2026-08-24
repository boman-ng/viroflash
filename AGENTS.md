# Repository Guidelines

## Project Structure & Module Organization

`viroflash` is a Rust 2021 command-line application. `src/main.rs` defines the CLI and exit behavior, while `src/lib.rs` orchestrates the detection pipeline. Domain modules own individual stages: `index.rs` and `reference.rs` build reusable indexes, `prescreen.rs` performs k-mer filtering, `align.rs` wraps minimap2 alignment, `stats.rs` evaluates candidates, and `report.rs` writes JSON/TSV output. Performance reporting lives in `perf.rs`. Integration coverage is in `tests/smoke.rs`; focused unit tests live beside their modules. `scripts/` contains portable synthetic-data utilities. Never commit local datasets, indexes, reports, `.local/`, `target/`, or `*.work/` directories.

## Build, Test, and Development Commands

Use the committed lockfile for reproducible builds:

```bash
cargo fmt --all -- --check                         # verify formatting
cargo check --locked                               # type-check quickly
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked                  # unit and smoke tests
cargo test --release --locked                      # optimized-path tests
cargo build --release --locked                     # production binary
cargo run --locked -- --help                       # inspect the CLI
```

Run the narrowest relevant test first, for example `cargo test stats::tests --locked`, then broaden checks for shared pipeline changes.

## Coding Style & Naming Conventions

Follow default `rustfmt`: four-space indentation, `snake_case` functions and modules, and `CamelCase` types. Keep thresholds in the module that owns the corresponding stage. Production input paths return `Result<T, String>` and must not panic or hide errors. Prefer the standard library and existing dependencies. Preserve deterministic ordering, fixed seeds, CLI options, JSON schemas, TSV column order, and 0-based half-open coordinates unless a contract change is explicitly required.

## Testing Guidelines

Add unit tests near changed logic and end-to-end behavior to `tests/smoke.rs`. Cover success and error paths for CLI or parsing changes. Reporting changes require JSON parsing plus JSON/TSV field-consistency checks. Concurrency, sampling, and index changes must verify deterministic output and equivalence between reusable-index and automatic-build execution.

## Commit & Pull Request Guidelines

Use Conventional Commits seen in history, such as `feat(index): ...`, `fix(report): ...`, or `test(smoke): ...`. Keep each commit limited to one self-contained purpose and review its staged diff before committing. Pull requests should explain the user-visible effect, affected contracts, verification commands, and any compatibility or performance impact. Link relevant issues and include representative output when report formats change. Do not push, tag, rewrite history, or modify remotes without explicit authorization.
