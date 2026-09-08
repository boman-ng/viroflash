# Viroflash Product Contract

## Product Value

Viroflash reports the fraction of input-library fragments attributed to each frozen target `ReferenceGroup` under one auditable analysis profile. The measurement supports reproducible review of reference-relative signal while preserving host conflict and ambiguity.

It does not identify absolute viral quantity, copies per volume, clinical infection, target absence, or member-level composition inside an exact-equivalence group.

## Inputs and Index

The production index accepts one HOST FASTA and one TARGET FASTA. FASTA IDs are the first non-empty whitespace-delimited header token and must be unique within each file. Sequence normalization removes ASCII whitespace, uppercases symbols, and accepts exactly `ACGTMRWSYKVHDBN`; any other symbol rejects the file.

TARGET records share a group only when complete normalized sequences are identical directly or by reverse complement. Group membership, order, representative, and multiplicity are frozen in the index ledger before FASTQ is read. Every run requires a current reusable profile-bound index.

## Measurement

The observation unit is a fragment. SE reads are fragments; PE mates are one fragment and are always selected and counted together.

Pass 1 validates the complete FASTQ input, checks PE counts and normalized IDs, computes exact fragment count `N`, and establishes the compressed-input digest. The frozen δ, β, and fixed target-family size determine Bernoulli inclusion probability `π`. Pass 2 deterministically selects fragments with a domain-separated BLAKE3 key based on profile identity, input identity, fragment ID, and ordinal.

The target-only 21-mer Bloom filter is solely a workload gate. Gate-passing fragments enter one HOST+TARGET competitive minimap2 index. Host evidence equal to or better than target evidence confounds the fragment. Equal evidence across target groups remains unresolved. Exact-equivalent members within one group preserve group attribution but not member attribution.

For each observed or indeterminate group, `report.csv` contains `x/n`, a simultaneous exact finite-population interval, and separate evidence and attribution states. Coverage, occupied windows, split alignments, discordance, and integration evidence remain diagnostics and never change the point estimate or create a binary decision.

## Outputs

`EvidenceReport` is the sole scientific result owner. Successful sample output contains exactly `report.csv`, `report.html`, and `perf.json`.

- `report.csv` contains one `RUN` row and only observed or indeterminate `TARGET_SIGNAL` rows.
- `report.html` renders the same report model and states the interpretation boundary.
- `perf.json` contains execution status, resources, stage times, and pipeline counts without biological conclusions.

The report directory becomes visible atomically only after all three successful artifacts are complete. An analysis error produces no success-shaped CSV or HTML.
