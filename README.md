# Viroflash

Viroflash analyzes FASTQ files for viral reference signals and reports their abundance, coverage, and sampling intervals. It estimates detectable viral fragments in the original library; it does not report viral load or clinical positive/negative status.

## Install

Use the **Linux x86_64 prebuilt package** supplied by your administrator or downloaded from [Releases](https://github.com/boman-ng/viroflash/releases). No Rust, compiler, or root access is needed. This repository is private; GitHub downloads require access.

From the download directory, replace `VERSION` with the package version:

```bash
tar -xzf viroflash-VERSION-linux-x86_64.tar.gz
cd viroflash-VERSION-linux-x86_64
install -Dm755 viroflash "$HOME/.local/bin/viroflash"
export PATH="$HOME/.local/bin:$PATH"
viroflash --help
```

The executable can also run directly after extraction as `./viroflash`. Add the PATH line to your shell startup file once. Installation and analysis work offline. Older releases may use a different CLI; use a package built from this revision for the commands below.

## Run

Create an index once for your host and viral reference FASTA files. The host is the sole background competitor; no separate decoy genome is needed:

```bash
viroflash index --host-fa host.fa --target-fa viruses.fa --out virus-index --threads 8
```

Analyze a sample using that index:

```bash
viroflash run --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --index virus-index --out sample-results --threads 8
```

Omit `--r2` for single-end data. Plain FASTQ and gzip are supported. Output directories must be new. `--threads` defaults to `1`; choose a worker count appropriate for your machine. Indexes from before reference descriptions were stored must be rebuilt.

All input fragments undergo Bloom prescreening once. Viroflash then samples the retained candidates deterministically and aligns them against the host and viral references. Temporary candidate data are removed when processing finishes; allow disk space for retained IDs and sequences beside the output directory.

Use `--precision fast|standard|sensitive` to select the minimum relevant **candidate-pool** fraction: 100, 10 (default), or 1 ppm. Each preset limits the probability of completely missing a signal at or above that fraction to 5% across the reference family. More sensitive presets usually align more candidates. This controls sampling loss, not clinical sensitivity or interval width.

## Read the results

Open `sample-results/report.html` in a browser. The English report has expandable details and an embedded CSV download. Use your browser to find text or print; signals are ranked by support.

| File | Purpose |
|---|---|
| `report.html` | Browse reference signals and review methods and full evidence. |
| `report.csv` | Analyze the same research values in a spreadsheet: 23 columns, ranked by supporting fragments. |
| `perf.json` | Inspect elapsed time, CPU usage, memory, and I/O. |

Start with the reference name, supporting fragments, **library abundance (ppm)**, coverage, and **target share (%)**. Library abundance scales candidate-sample support back to the original input and displays a simultaneous 95% interval; target share divides observed support by all attributed target fragments. One paired-end read pair counts as one fragment. Reference groups are not species counts.

A sample without attributed targets retains a sample-only CSV row. Interval bounds describe candidate sampling uncertainty and do not correct prescreen or attribution losses; coverage and split/discordant counts provide supporting context. The HTML file is fully standalone: viewing and CSV download need no neighboring files or JavaScript.

## Containers and development

HPC users can run a supplied SIF with `apptainer run viroflash-VERSION-x86_64.sif --help`; pass the same `index` and `run` arguments and bind input/output directories as needed.

See [analysis methods](docs/analysis-contract.md), [contributor instructions](AGENTS.md), and [packaging and releases](docs/releasing.md) for details beyond normal use.
