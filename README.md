# Viroflash

Viroflash is a Rust command-line tool that measures **profile-attributed fragment fraction** for a frozen family of target reference groups. It is a library-fragment measurement under a defined software profile, not viral load, absolute quantitation, or a clinical positive/negative result.

## Build

```bash
cargo build --release --locked
```

The binary remains version `0.3.0`.

## Commands

Build a reusable HOST+TARGET index:

```bash
viroflash index \
  --host-fa HOST.fa \
  --target-fa TARGET.fa \
  --out INDEX_DIR \
  [--threads N]
```

Analyze single-end or paired-end FASTQ:

```bash
viroflash run \
  --r1 SAMPLE_R1.fastq.gz \
  [--r2 SAMPLE_R2.fastq.gz] \
  --index INDEX_DIR \
  --out SAMPLE_REPORT_DIR \
  [--threads N]
```

`--out` must not already exist. A successful run creates exactly:

```text
report.csv
report.html
perf.json
```

`report.csv` is the machine-readable scientific artifact. `report.html` renders the same immutable `EvidenceReport`; it does not recompute statistics. `perf.json` contains execution telemetry, not biological conclusions. A failed analysis may create only an error-shaped `perf.json`.

## Analysis Contract

The only production `AnalysisProfile` freezes:

- minimum relevant fraction δ = `1e-5`;
- familywise miss probability β = `0.05`;
- familywise interval error α = `0.05`;
- target k-mer length `21`, minimap2 `sr` with all-chain enumeration, and ten diagnostic windows.

The profile digest is the lowercase SHA-256 of the exact committed `evaluation/phase0/analysis-profile.json` bytes. Scientific parameters are not CLI options.

Target records are grouped during indexing only when their complete uppercase IUPAC sequence is byte-identical, directly or after reverse complementation. Invalid symbols reject the FASTA. Each target record belongs to exactly one fixed `ReferenceGroup`.

FASTQ analysis has two streaming passes:

1. validate every record and PE identifier, count fragments exactly, and compute the compressed-byte input identity from that same stream;
2. repeat validation and digesting on the analysis stream, reject any cross-pass change, apply deterministic BLAKE3 Bernoulli inclusion, then run a target-only 21-mer Bloom workload gate and HOST+TARGET competitive alignment through a bounded worker queue.

For `N` input fragments and `m` fixed target groups:

```text
M_min = ceil(δN)
π = min(1, 1 - (β/m)^(1/M_min))
```

PE ends share one selection key and one attribution. A fragment contributes at most once to one group. Host ties or advantages are confounded; cross-group ties remain unresolved; multiple exact-equivalent members within one group retain group-level support. A split diagnostic requires disjoint HOST and TARGET chains on one read end with a supplementary alignment; TARGET alternatives alone are not split evidence.

Selected fragments shorter than 21 bases or without an encodable 21-mer are counted explicitly and produce `CONFORMANT_WITH_LIMITATIONS`; they are not silently treated as target-negative.

The point estimate is `x/n`. Simultaneous intervals use Bonferroni `α/m` and equal-tailed exact hypergeometric inversion conditional on the realized `n`. Census intervals collapse to the exact fraction. The implementation follows the finite-population inversion in [samplingbook `Sprop`](https://rdrr.io/cran/samplingbook/src/R/Sprop.R); distribution evaluation uses the maintained Rust `statrs` implementation. Coverage, windows, split, discordant, and integration fields are diagnostics only.

## Verification

```bash
cargo fmt --all -- --check
cargo check --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
cargo test --release --locked
cargo build --release --locked
python3 evaluation/phase0/verify.py
```
