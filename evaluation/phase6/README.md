# Phase 6 real E2E

These frozen, standard-library-only tools execute and score the 127-run manifest without using
historical results during production execution.

```bash
python3 evaluation/phase6/phase6_runner.py prepare --campaign-dir .tmp/phase6-real-e2e
python3 evaluation/phase6/phase6_runner.py run --campaign-dir .tmp/phase6-real-e2e
python3 evaluation/phase6/phase6_runner.py repro --campaign-dir .tmp/phase6-real-e2e
python3 evaluation/phase6/phase6_scorer.py --campaign-dir .tmp/phase6-real-e2e
```

`prepare` builds new reusable indexes and freezes actual binary, profile, reference, and every
index-artifact SHA-256 before any sample launch. It also writes `scoring-provenance.json`, binding
the scorer, shared evaluator, config, manifest identity, references, and historical scoring inputs.
`run`, `repro`, and scoring reject manifest or scoring-input drift. `run` verifies all original
FASTQ artifact sizes and SHA-256 values, then launches each manifest run once with its frozen cohort
resource settings.
The append-only ledger preserves the first terminal evidence and the runner refuses an existing
ledger or output root. `phase6_scorer.py` reads historical results only for post-run paired
evaluation. The frozen reproducibility matrix runs with at most eight concurrent jobs. Campaign
indexes, logs, reports, ledgers, and evidence remain under the ignored
`.tmp/` campaign directory.

The completed 2026 Phase 6 campaign predates pre-execution scoring provenance. Before its single
reviewer-requested rescore, create the separate retrospective ledger without changing
`digest-ledger.json`:

```bash
python3 evaluation/phase6/phase6_scorer.py \
  --campaign-dir .tmp/phase6-real-e2e \
  --write-retrospective-provenance
```

That ledger explicitly records the limitation and the four post-run evaluation-only corrections:
all exact-sequence ReferenceGroups matching an internal truth label participate in expected-signal
adjudication (`853677e`), and internal adjudication is limited to the frozen EBV/HBV/HPV16/HPV18
label scope while retaining all evidence (`26aeab2`). Taxon matching starts at the FASTA
description, excludes non-human `Heron hepatitis B virus`, and is accompanied by independent public
report semantic validation and scoring-provenance binding (`e55dbd2`). It does not claim any of
these scorer inputs or corrections were frozen before production execution. A bounded
metadata-prefix grammar then restores legitimate human-virus labels without admitting Heron HBV or
HPV numeric lookalikes, while exact enum/status validation and an independent cached
hypergeometric inversion oracle preserve interval mismatches as a false report-validity gate
(`5fc4efb`).
