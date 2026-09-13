# Viroflash

Viroflash is a Rust CLI for fast detection of viral reference signals in FASTQ data. Supply a host reference and a viral panel, build an index, and analyze samples locally.

## Install

Extract the Linux x86_64 prebuilt package supplied by your administrator or available from [Releases](https://github.com/boman-ng/viroflash/releases):

```bash
tar -xzf viroflash-VERSION-linux-x86_64.tar.gz
cd viroflash-VERSION-linux-x86_64
install -Dm755 viroflash "$HOME/.local/bin/viroflash"
export PATH="$HOME/.local/bin:$PATH"
viroflash --help
```

The executable also runs directly as `./viroflash`. Installation and analysis work offline and need no compiler. GitHub downloads require access to this private repository.

## Analyze a sample

Build a reusable index from host and target FASTA files:

```bash
viroflash index --host-fa host.fa --target-fa viruses.fa --out virus-index --threads 8
```

Run a paired-end sample:

```bash
viroflash run --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --index virus-index --out sample-results --threads 8
```

Omit `--r2` for single-end input. Plain and gzip FASTQ are supported. Output directories must be new. `--threads` defaults to `1` and sets the number of analysis workers; input readers also use background threads.

| Mode | Workflow |
|---|---|
| `screen` (default) | Sample input fragments, apply Bloom, then competitively align against HOST+TARGET. |
| `full` | Apply Bloom to all input fragments, sample candidates, then competitively align. |

Both modes read FASTQ once and retain an in-memory deterministic bottom-k sample. They use the same sample capacity; full selects from the candidate pool and therefore retains at least the candidates analyzed by screen for the same input, index and precision.

`--precision` controls sample capacity in either mode:

| Precision | Relative-abundance target |
|---|---:|
| `fast` | 100 ppm |
| `standard` (default) | 10 ppm |
| `sensitive` | 1 ppm |

Capacity is derived from the target, panel size and a 5% familywise sampling-miss budget. With 20,560 reference groups, standard uses a capacity of 1,292,678 fragments. Smaller populations are fully selected. These settings control sampling loss; they do not specify interval width.

```bash
viroflash run --r1 sample.fastq.gz --index virus-index --out full-results \
  --mode full --precision standard --threads 8
```

## Read the report

| File | Contents |
|---|---|
| `report.html` | Top 20 supported reference groups, expandable data and an embedded download of the complete CSV. |
| `report.csv` | All supported groups, 24 columns, ranked by supporting fragments. |
| `perf.json` | Wall-time, CPU time, peak memory, I/O, stage timings and fragment counts. |

Start with the reference description, support count, whole-library abundance interval and coverage. One SE read or PE pair is one fragment. Each row represents identical or reverse-complement reference sequences grouped together.

`support_pct`, `support_ci_lower_pct` and `support_ci_upper_pct` are numeric percentages: `0.013` means `0.013%`. They estimate detectable, attributed fragments over the original input library. Intervals describe sampling uncertainty. `target_support_share_pct` uses all attributed target fragments as its denominator.

`sampling_target_score` compares the estimated whole-library fraction **p** with the selected ppm target **δ**:

```text
score = (p - δ) / (p + δ)
```

Negative values are below the target, zero is at the target, and positive values are above it. Half the target gives −1/3; twice the target gives +1/3. The score describes relative abundance, not confidence. Coverage and alignment diagnostics remain separate measurements.

The English HTML is standalone and script-free. Viewing it and downloading its CSV require no neighboring files. Samples without attributed groups retain a sample-only CSV row.

## Containers and development

Run a supplied SIF with `apptainer run viroflash-VERSION-x86_64.sif --help`, using the same `index` and `run` arguments and the appropriate input/output binds.

See [analysis methods](docs/analysis-contract.md), [contributor instructions](AGENTS.md), and [packaging](docs/releasing.md).
