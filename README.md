# viroflash

`viroflash` detects viral candidates in FASTQ/FASTQ.gz sequencing data and produces an offline, interactive HTML report plus analysis-ready CSV, JSON, and TSV files.

> Research use only. A `PASS` candidate crossed the current exploratory evidence gates; it is not a clinical diagnosis.

## Quick start

### 1. Download and unpack

Download the Linux amd64 archive and its `.sha256` file from [GitHub Releases](https://github.com/boman-ng/viroflash/releases). No Rust, minimap2, container runtime, or system library installation is required.

```bash
VERSION=0.2.0
ARCHIVE="viroflash-${VERSION}-linux-x86_64.tar.gz"

sha256sum --check "${ARCHIVE}.sha256"
tar -xzf "${ARCHIVE}"
cd "${ARCHIVE%.tar.gz}"
./viroflash version
```

The downloadable binary currently supports 64-bit Linux on x86_64/amd64 CPUs.

### 2. Run a sample

The shortest path builds a temporary index from your references and runs detection in one command:

```bash
./viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --host-fa host.fa \
  --target-fa viruses.fa \
  --threads 8 \
  --out result
```

For single-end data, omit `--r2`. `--host-fa` and `--target-fa` are required; contaminant and decoy references are optional:

```bash
  --contam-fa contaminants.fa \
  --decoy-fa decoys.fa
```

When `--decoy-fa` is omitted, viroflash generates deterministic synthetic decoys. These are statistical stress references, not replacements for laboratory negative controls.

### 3. Review the result

Open `result.html` in a web browser. It is self-contained, works offline, and supports candidate search, filtering, sorting, expandable evidence, printing, and CSV download.

| File | Use |
| --- | --- |
| `result.html` | Primary human-readable report |
| `result.csv` | One row per reported candidate for spreadsheets or downstream analysis |
| `result.json` | Complete structured result and audit metadata |
| `result.tsv` | Fixed-column compatibility table |
| `result.perf.json`, `result.perf.tsv` | Runtime and resource metrics |

## Reuse an index

For repeated analysis against the same references, build the index once:

```bash
./viroflash index \
  --host-fa host.fa \
  --target-fa viruses.fa \
  --out virus-index \
  --threads 8
```

Then reuse it for each sample:

```bash
./viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --index virus-index \
  --threads 8 \
  --out result
```

An index contains its reference metadata and cannot be combined with `--host-fa`, `--target-fa`, `--contam-fa`, or `--decoy-fa` at run time.

## Interpret results

- Results are candidate-scoped. No reported candidates does not mean `NOT_DETECTED`.
- `PASS` means a candidate passed the configured statistical, breadth, and positional-distribution gates. It does not mean clinically positive.
- `q_value` is a model-adjusted p-value, not a validated classical FDR or confidence score.
- An OR hypothesis lists unresolved alternatives; do not attribute evidence to one member without additional analysis.
- `quality_control.status=NOT_EVALUATED` means audit metrics are present but validated QC thresholds have not been applied.
- Validate reference scope, sampling capacity, synthetic decoys, gates, LoD/LoB, near-neighbor interference, and independent controls for the intended assay before operational use.

The HTML and CSV schemas are presentation layers over `viroflash.result.v1`; JSON remains the complete source of truth. Runtime operation is offline and does not query taxonomy or external databases.

## Containers

Use the standalone archive for ordinary Linux use. Versioned Docker/OCI and Apptainer images are also published for managed or HPC environments:

```bash
docker run --rm -v "$PWD:/work" ghcr.io/boman-ng/viroflash:0.2.0 version
apptainer run --bind "$PWD:/work" --pwd /work viroflash-0.2.0-x86_64.sif version
```

See the [release guide](https://github.com/boman-ng/viroflash/blob/master/docs/releasing.md) for image details, supported platforms, and release verification.

## Build from source

Source builds require the pinned Rust toolchain and a C compiler:

```bash
cargo build --release --locked
cargo test --all-targets --locked
```

Run `./viroflash --help` from an unpacked release, or `cargo run --locked -- --help` from source, for the complete CLI reference.
