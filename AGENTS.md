# Working on Viroflash

Viroflash is a Rust CLI. Keep the user path simple: install a prebuilt package, build a reference index, run a sample, open its report.

## Where changes belong

- `src/main.rs`: CLI; `pipeline.rs`: orchestration.
- `profile.rs` and `analysis-profile.json`: frozen base profile and index identity. Run precision and sampling population are recorded separately. Preserve the JSON bytes and digest when relocating files.
- `fastq.rs`, `candidates.rs`, `sampling.rs`, and `evidence.rs`: input, candidate spooling, selection, and evidence aggregation.
- `alignment/`: mapper and attribution; `src/workers.rs`: shared bounded execution for prescreen and alignment. `gate/`: Bloom prescreen; `sdust.rs`: low-complexity masking.
- `index/`: reusable HOST+TARGET indexes and reference grouping. HOST is the only background competitor; there is no separate decoy input.
- `report/`: research values and embedded HTML/CSS; `telemetry.rs`: performance; `output.rs`: shared atomic writes. `tests/fixtures/output-fields.tsv` records the output contract.
- `.github/workflows/`: Rust CI and prebuilt releases; `.github/scripts/create-smoke-fixture.rs`: release fixture.

## Contracts to preserve

- One SE read or PE pair is one fragment; deterministic selection and attribution must agree across worker counts and FASTQ compression formats.
- A fragment supports at most one reference group. Library ppm estimates candidate support over all original input fragments; target share uses all attributed target fragments. Preserve exact intervals and outward rounding.
- Analysis consumes a reusable index, including its reference descriptions. Invalid indexes fail explicitly; do not infer missing metadata or rebuild automatically.
- Successful output is exactly `report.csv`, `report.html`, and `perf.json`; failure may leave only error telemetry. Preserve atomic writes.
- CSV has 23 ordered columns with BOM and CRLF. HTML shares those research values and retains full evidence. Keep the interface English, script-free, and standalone, including embedded CSV download; preserve source descriptions.
- Use existing dependencies, bounded streaming, explicit errors, and 0-based half-open coordinates. Only `--precision fast|standard|sensitive` changes sampling sensitivity; presets use the candidate population. Keep other scientific constants fixed.

## Verify the changed behavior

Use the pinned toolchain and lockfile. Tests are Rust and need no Python, JavaScript, or local datasets. Run focused tests first; for Rust or shared pipeline changes:

```bash
cargo fmt --all -- --check
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
cargo test --release --test smoke --locked
```

Keep tests that catch distinct failures: malformed or changed input, index reuse, attribution, sampling/interval correctness, report parity, and bounded concurrency. Avoid duplicate assertions and Cartesian test matrices. CI runs formatting, Clippy, and core tests; release preparation tests the installed static binary and compares SIF/OCI reports. UI changes need fresh desktop/mobile reports and interaction checks. Documentation-only changes need link/contract checks and `git diff --check`.

## Local work and change control

`evaluation/`, `.local/`, `.tmp/`, `target/`, `dist/`, `*.work/`, and generated images are local artifacts. Keep them out of Git, default searches, and container build contexts. Build and test must work without `evaluation/`; preserve existing local results. Use scoped `rg --hidden --no-ignore` when inspecting them.

Preserve unrelated work. Review the complete diff before a focused Conventional Commit. Pushes, tags, releases, history/remote changes, and package/schema version changes require explicit authorization.
