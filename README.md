# viroflash

`viroflash` is a Rust command-line tool that produces viral candidate reports from FASTQ or FASTQ.gz sequencing data. Its pipeline combines deterministic bottom-k sampling, k-mer prescreening, minimap2 whole-chain competitive alignment, discovery/validation evidence hypotheses, synthetic-decoy background statistics, and JSON/TSV reporting. Both `index` and `run` also emit standalone performance reports. Runtime operation is fully offline: it does not query NCBI, taxonomy services, or panel-specific prior knowledge.

## Build and test

```bash
cargo build --release --locked
cargo test --all-targets --locked
cargo run --locked -- --help
```

## Containers

Published releases provide an amd64 Docker/OCI image at `ghcr.io/boman-ng/viroflash` and an amd64 Apptainer SIF file on the GitHub release. Input datasets and indexes are mounted at runtime; they are never embedded in an image.

```bash
docker run --rm -v "$PWD:/work" ghcr.io/boman-ng/viroflash:0.1.0 version
apptainer run --bind "$PWD:/work" --pwd /work viroflash-0.1.0-x86_64.sif version
```

Build either image locally with `Dockerfile` or `Apptainer.def`. See [the CI and release strategy](docs/releasing.md) for supported architectures, release artifacts, bind-mount permissions, and checksum verification.

## Build an index

```bash
viroflash index \
  --host-fa host.fa \
  --target-fa target.fa \
  --contam-fa contaminant.fa \
  --decoy-fa decoy.fa \
  --out index_dir \
  --k 21 \
  --threads 8
```

The reusable index directory contains:

- `ref.mmi`: a single-shard minimap2 `sr` index, keeping host, target, decoy, and contaminant hits in one competitive MAPQ space.
- `bloom.bin`: a canonical k-mer Bloom filter built from all original targets and the decoys used by the index.
- `targets.fa`: sequence-derived detection-group representatives included in the MMI.
- `manifest.json`: format version, k, role-to-contig mappings, representative-to-original-member mappings, BLAKE3 source checksums, and decoy provenance.

Detection groups use sequence evidence only—never names, taxonomy, or database lookups. Exact duplicates are collapsed first; fixed k=21 FracMinHash/PPJoin and Mash screening proposes candidates; minimap2 then confirms representative-centered groups at identity >=97% and bidirectional coverage >=95%. Grouping compresses the index without adding member counts. Runtime OR hypotheses still arise only from each read's directly compatible top-scoring target set.

`--contam-fa` and `--decoy-fa` are optional. Without a decoy FASTA, the index command deterministically generates SNP-only synthetic decoys from each group representative using `--decoy-ani` (default `85`), `--decoy-per-layer` (default `1`), and `--decoy-seed` (default `0`). Synthetic decoys are competition and null-model stress references; they are not substitutes for extraction blanks, batch-matched NTCs, or laboratory contamination controls. Their exchangeability has not been independently validated, so reports do not claim classical FDR control.

The command writes `<out>.perf.json` and `<out>.perf.tsv` beside the index. Existing multi-shard MMI indexes are rejected because cross-shard MAPQ does not satisfy the competitive-alignment contract.

## Run detection

Use a reusable index:

```bash
viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --index index_dir \
  --threads 8 \
  --out result
```

Single-end input is supported by omitting `--r2`. `--index` is mutually exclusive with reference FASTA options, and `--k` must match the index.

Without `--index`, `run` builds the same index under `<out>.work/index/`:

```bash
viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --host-fa host.fa \
  --target-fa target.fa \
  --threads 8 \
  --out result
```

Both paths share one index builder and produce equivalent candidate results. `--k` accepts values from 1 through 31. `--threads` is the compute budget for the data pipeline; bounded decompression helpers and one low-frequency telemetry sampler may run in addition to compute workers.

All fragments first enter a fixed-capacity, domain-separated BLAKE3 bottom-k sample. Selected pairs are then prescreened and aligned. Fragments are deterministically split into discovery and validation folds, and both ends of a pair always stay in the same fold. A read end compatible with several top-scoring targets contributes one set-valued observation, not one vote per member. A reported OR hypothesis means that at least one member is present under the current reference model; it does not resolve an individual accession.

Discovery fixes the target hypotheses to test. Validation counts each direct read end once per hypothesis. Single-member hypotheses are compared with fixed synthetic decoys in the same size/GC stratum; multi-member hypotheses use the fixed global decoy set. The primary test is a one-sided exact conditional two-Poisson rate test with reference length as exposure. All discovery hypotheses, including those with zero validation counts, enter a fixed-family Benjamini-Hochberg adjustment. The compatibility field `q_value` therefore contains a model adjusted p-value, not a validated classical FDR or clinical confidence score.

The exploratory candidate gate requires model adjusted p-value `<0.2`, representative observed breadth `>=10%`, and direct evidence in at least 3 of 10 fixed positional windows. `PASS` and `BELOW_THRESHOLD` are candidate-level research states. Host-virus split/site evidence is reported separately as `integration_evidence`; it neither subtracts detection reads nor determines general viral detection.

Generate deterministic synthetic inputs with:

```bash
python3 scripts/make_virtual_data.py synthetic.work --scenario all
```

## Output

`run` writes:

- `<out>.json`: run summary, sampling/fold/audit metadata, index provenance, and compound candidate details.
- `<out>.tsv`: the fixed 22-column candidate table.
- `<out>.perf.json`: process/system CPU, RSS, virtual memory, cumulative CPU time, I/O, sampling completeness, and stage-level timing/CPU/RSS.
- `<out>.perf.tsv`: a fixed-column performance summary for batch aggregation.

The current candidate JSON contract is `viroflash.result.v1`. Performance reports never include input or index paths. RSS includes shared pages and is not PSS; telemetry does not claim NUMA, cache-miss, or memory-bandwidth measurements.

## Interpretation and limitations

- Output is candidate-scoped, not a sample-level diagnosis. `candidates: []` means that no candidate was reported; it does not mean `NOT_DETECTED`.
- `decision=PASS` means the candidate passed the current exploratory gate, not that the sample is clinically positive. `confidence=UNVALIDATED` records the lack of external null/FDR calibration.
- `p_value` is the exact conditional rate-test result. `q_value` is its fixed-discovery-family BH adjustment. Log-scale fields preserve meaning when ordinary floating-point values underflow to zero.
- `candidate_testing` exposes the complete discovery family, reported candidates, and omitted zero-evidence hypotheses so the multiplicity boundary remains auditable.
- Legacy `poisson_p` is diagnostic only; `nb_p` is currently disabled and emitted as `null`.
- OR `members` are unresolved alternatives and must not be attributed individually. An empty `taxid` confirms that viroflash did not query taxonomy.
- `quality_control.status=NOT_EVALUATED` means that observed audit metrics are reported without validated QC pass/fail thresholds.
- Sampling capacity, reference scope, synthetic decoys, breadth/window gates, and alignment gates require validation against the intended sample matrix, viral classes, LoD/LoB, near-neighbor interference, and independent negative controls before supporting clinical performance or general sensitivity/specificity claims.
