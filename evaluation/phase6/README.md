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
index-artifact SHA-256 before any sample launch. `run` verifies all original FASTQ artifact sizes
and SHA-256 values, then launches each manifest run once with its frozen cohort resource settings.
The append-only ledger preserves the first terminal evidence and the runner refuses an existing
ledger or output root. `phase6_scorer.py` reads historical results only for post-run paired
evaluation. The frozen reproducibility matrix runs with at most eight concurrent jobs. Campaign
indexes, logs, reports, ledgers, and evidence remain under the ignored
`.tmp/` campaign directory.
