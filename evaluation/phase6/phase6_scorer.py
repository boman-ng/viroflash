#!/usr/bin/env python3
import argparse
import csv
import json
import math
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import NormalDist

from phase6_common import (
    CONFIG_PATH, MANIFEST_PATH, atomic_json, atomic_text, index_id_for, load_contract,
    load_digest_ledger, load_json, output_dir, parse_time_verbose, require, sha256_file,
    validate_error_directory, validate_success_directory, verify_campaign_digests,
)


def wilson(successes, total, confidence=0.95):
    if total == 0:
        return {"successes": 0, "total": 0, "rate": None, "lower": None, "upper": None, "confidence": confidence}
    z = NormalDist().inv_cdf(0.5 + confidence / 2)
    rate = successes / total
    denominator = 1 + z * z / total
    center = (rate + z * z / (2 * total)) / denominator
    radius = z * math.sqrt(rate * (1 - rate) / total + z * z / (4 * total * total)) / denominator
    return {"successes": successes, "total": total, "rate": rate, "lower": max(0.0, center - radius), "upper": min(1.0, center + radius), "confidence": confidence}


def percentile(values, fraction):
    require(values, "cannot summarize an empty performance series")
    ordered = sorted(values)
    rank = (len(ordered) - 1) * fraction
    lower = math.floor(rank)
    upper = math.ceil(rank)
    return ordered[lower] if lower == upper else ordered[lower] + (ordered[upper] - ordered[lower]) * (rank - lower)


def summary(values):
    return {"median": statistics.median(values), "p95": percentile(values, 0.95), "max": max(values)}


def fasta_labels(path, terms):
    labels = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            member = line[1:].split(maxsplit=1)[0]
            header = line[1:]
            matched = [label for label, patterns in terms.items() if any(pattern.casefold() in header.casefold() for pattern in patterns)]
            require(len(matched) <= 1, f"reference header maps to multiple expected labels: {member}")
            if matched:
                labels[member] = matched[0]
    return labels


def expected_signal(run, targets, internal_labels):
    expected = run.get("expected_group_key")
    if expected is None:
        return None
    for target in targets:
        members = set(target["member_ids"].split(";")) | {target["representative_id"]}
        if run["dataset_id"] == "internal-68":
            if any(internal_labels.get(member) == expected for member in members):
                return target
        elif expected in members:
            return target
    return None


def observed(target):
    return target["evidence_status"] == "REFERENCE_SIGNAL_OBSERVED"


def adjudicate(run, parsed, internal_labels):
    targets = parsed["targets"]
    expected = expected_signal(run, targets, internal_labels)
    observed_targets = [target for target in targets if observed(target)]
    expected_observed = expected is not None and observed(expected)
    if run["expectation_kind"] == "MOCK":
        classification = "LABEL_CONCORDANT_NO_SIGNAL" if not observed_targets else "LABEL_DISCORDANT_WITH_SEQUENCE_EVIDENCE"
    elif expected_observed:
        classification = "LABEL_CONCORDANT_SIGNAL"
    elif observed_targets:
        classification = "LABEL_DISCORDANT_WITH_SEQUENCE_EVIDENCE"
    else:
        classification = "LABEL_DISCORDANT_UNRESOLVED"
    wrong_resolved = any(
        target is not expected and target["attribution_status"] == "RESOLVED_TO_REFERENCE_GROUP"
        for target in observed_targets
    )
    ambiguous = any(target["attribution_status"] in {"AMBIGUOUS_WITHIN_GROUP", "UNRESOLVED_ACROSS_GROUPS"} for target in targets)
    evidence_summary = [
        {
            key: target[key]
            for key in (
                "target_group_id", "representative_id", "member_ids", "evidence_status",
                "attribution_status", "supporting_selected_fragments",
                "selected_fragment_denominator", "attributed_fragment_fraction",
                "interval_lower", "interval_upper", "covered_bases", "coverage_fraction",
                "occupied_windows", "host_confounded_fragments",
                "cross_group_ambiguous_fragments", "integration_status", "split_events",
                "discordant_fragments",
            )
        }
        for target in targets
    ]
    return {
        "classification": classification,
        "expected_observed": expected_observed,
        "expected_signal": expected,
        "observed_signal_count": len(observed_targets),
        "wrong_group_resolved": wrong_resolved,
        "ambiguous_attribution": ambiguous,
        "evidence_summary": evidence_summary,
    }


def historical_path(run, config):
    if run["dataset_id"] == "internal-68":
        return Path(config["historical"]["internal_results"]) / f"{run['run_id']}.json"
    return Path(config["historical"]["external_results"]) / run["cohort"] / f"{run['run_id']}.json"


def historical_perf_path(run, config):
    return historical_path(run, config).with_suffix(".perf.json")


def old_expected_observed(run, old, internal_labels):
    expected = run.get("expected_group_key")
    if expected is None:
        return False
    for candidate in old.get("candidates", []):
        if candidate.get("decision") != "PASS":
            continue
        members = set(candidate.get("hypothesis", {}).get("members", [])) | {candidate.get("representative")}
        if run["dataset_id"] == "internal-68":
            if any(internal_labels.get(member) == expected for member in members):
                return True
        elif expected in members:
            return True
    return False


def old_perf(perf):
    return {
        "wall_time_ms": perf["run"]["wall_time_ms"],
        "peak_rss_bytes": perf["process"]["rss_bytes_peak"],
        "stages": {stage["name"]: stage["wall_time_ms"] for stage in perf.get("stages", [])},
    }


def terminal_events(path):
    starts = {}
    terminals = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            event = json.loads(line)
            key = (event["dataset_id"], event["run_id"])
            if event["event"] == "STARTED":
                require(key not in starts, f"duplicate start event: {key}")
                starts[key] = event
            elif event["event"] == "TERMINAL":
                require(key not in terminals, f"duplicate terminal event: {key}")
                terminals[key] = event
    return starts, terminals


def reproducibility_summary(path, config):
    starts = set()
    terminals = {}
    with Path(path).open(encoding="utf-8") as handle:
        for line in handle:
            event = json.loads(line)
            key = (event["run_id"], event.get("threads"), event.get("repetition"), event.get("encoding"))
            if event["event"] == "STARTED":
                require(key not in starts, f"duplicate reproducibility start: {key}")
                starts.add(key)
            elif event["event"] == "TERMINAL":
                require(key not in terminals, f"duplicate reproducibility terminal: {key}")
                terminals[key] = event
    expected = {
        (item["run_id"], threads, repetition, None)
        for item in config["reproducibility_runs"]
        for threads in config["reproducibility_threads"]
        for repetition in range(1, config["reproducibility_repetitions"] + 1)
    }
    expected |= {
        (config["gzip_equivalence_run_id"], None, None, "plain"),
        (config["gzip_equivalence_run_id"], None, None, "multimember"),
    }
    require(starts == expected and set(terminals) == expected, "reproducibility ledger is incomplete")
    require(all(event["status"] == "SUCCESS" and event["report_csv_byte_identical"] for event in terminals.values()), "reproducibility report.csv mismatch or failure")
    return {"runs": len(terminals), "report_csv_byte_identical": len(terminals), "threads": config["reproducibility_threads"], "repetitions": config["reproducibility_repetitions"], "encoding_variants": ["plain", "multimember-gzip"]}


def score(campaign_dir):
    manifest, config = load_contract()
    campaign_dir = Path(campaign_dir).resolve()
    evidence_dir = campaign_dir / "evidence"
    require(not evidence_dir.exists(), "campaign evidence already exists; frozen scoring outputs are write-once")
    digest_ledger = verify_campaign_digests(campaign_dir)
    starts, terminals = terminal_events(campaign_dir / "run-ledger.jsonl")
    expected_keys = {(run["dataset_id"], run["run_id"]) for run in manifest["runs"]}
    require(set(starts) == expected_keys and set(terminals) == expected_keys, "campaign ledger does not contain exactly 127 started and terminal runs")
    reproducibility = reproducibility_summary(campaign_dir / "reproducibility/repro-ledger.jsonl", config)
    internal_target = next(run["target_reference"]["path"] for run in manifest["runs"] if run["dataset_id"] == "internal-68")
    internal_labels = fasta_labels(internal_target, config["internal_label_header_terms"])
    require(all(any(label == expected for label in internal_labels.values()) for expected in config["internal_label_header_terms"]), "an internal expected label maps to no frozen reference")
    records = []
    failures = []
    for run in manifest["runs"]:
        key = (run["dataset_id"], run["run_id"])
        terminal = terminals[key]
        directory = output_dir(campaign_dir, run)
        if terminal["status"] == "SUCCESS":
            index = digest_ledger["indexes"][index_id_for(run, config)]
            parsed = validate_success_directory(directory, run, digest_ledger["profile"]["sha256"], index["index_digest"])
            decision = adjudicate(run, parsed, internal_labels)
            records.append({"manifest": run, "parsed": parsed, "decision": decision, "terminal": terminal})
        else:
            failure = {"dataset_id": run["dataset_id"], "cohort": run["cohort"], "run_id": run["run_id"], "terminal": terminal}
            if terminal["status"] == "ERROR":
                failure["perf"] = validate_error_directory(directory)
            failures.append(failure)
    scientific = scientific_metrics(records, manifest)
    performance, paired = performance_metrics(records, config, digest_ledger, starts, terminals, internal_labels)
    adjudication_rows = [record for record in records if record["decision"]["classification"].startswith("LABEL_DISCORDANT")]
    evidence_dir.mkdir()
    write_adjudication(evidence_dir / "adjudication.tsv", adjudication_rows)
    write_run_table(evidence_dir / "run-results.tsv", records, failures)
    write_paired(evidence_dir / "old-system-paired.tsv", paired)
    summary_value = {
        "contract_id": config["contract_id"],
        "scored_at": datetime.now(timezone.utc).isoformat(),
        "terminal": {"total": 127, "success": len(records), "failure": len(failures)},
        "hard_gates": {
            "exactly_127_started_once": len(starts) == 127,
            "exactly_127_terminal": len(terminals) == 127,
            "all_success_directories_exact": not failures,
            "binary_profile_index_digests_reverified": True,
            "pass1_pass2_counts_ids_input_digest": "ENFORCED_BY_VERIFIED_BINARY_BEFORE_SUCCESS_REPORT",
            "reproducibility_report_csv_byte_identical": reproducibility["runs"] == reproducibility["report_csv_byte_identical"],
            "cargo_version_unchanged_0_3_0": cargo_version_is_frozen(),
        },
        "reproducibility": reproducibility,
        "scientific": scientific,
        "performance": performance,
        "failures": failures,
        "evidence_limits": [
            "Internal 68 are development/regression evidence, not an independent blind validation set.",
            "External 59 have historical results and manual review, so they are generalization regression rather than a pristine holdout.",
            "Labels cover only mapped groups and do not establish sensitivity or specificity across all indexed references.",
            "Real samples do not establish abundance bias, interval coverage, LoD/LoQ, absolute quantitation, or clinical validity.",
            "Pass1/Pass2 pair-ID and decoded-byte equality are binary-enforced invariants; reports expose the terminal digest/count consequences, not separate per-pass pair-ID ledgers.",
        ],
    }
    atomic_json(evidence_dir / "summary.json", summary_value)
    artifact_digests = {entry.name: sha256_file(entry) for entry in sorted(evidence_dir.iterdir()) if entry.is_file()}
    atomic_json(evidence_dir / "evidence-digests.json", artifact_digests)
    print(json.dumps(summary_value["terminal"], sort_keys=True))


def cargo_version_is_frozen():
    for line in (Path(__file__).resolve().parents[2] / "Cargo.toml").read_text().splitlines():
        if line.startswith("version = "):
            return line == 'version = "0.3.0"'
    return False


def scientific_metrics(records, manifest):
    expected = [record for record in records if record["manifest"]["expectation_kind"] == "EXPECTED_GROUP"]
    mocks = [record for record in records if record["manifest"]["expectation_kind"] == "MOCK"]
    result = {
        "expected_group_observed": wilson(sum(record["decision"]["expected_observed"] for record in expected), len(expected)),
        "mock_unexpected_signal": wilson(sum(record["decision"]["observed_signal_count"] > 0 for record in mocks), len(mocks)),
        "wrong_group_resolved_signal": wilson(sum(record["decision"]["wrong_group_resolved"] for record in records), len(records)),
        "ambiguous_attribution": wilson(sum(record["decision"]["ambiguous_attribution"] for record in records), len(records)),
        "not_evaluable": wilson(sum(run["evaluability_status"] == "NOT_EVALUABLE" for run in manifest["runs"]), len(manifest["runs"])),
    }
    interval_fields = ("interval_lower", "interval_upper")
    expected_signals = [record["decision"]["expected_signal"] for record in expected if record["decision"]["expected_observed"]]
    result["expected_group_interval_summaries"] = {
        field: summary([float(signal[field]) for signal in expected_signals]) if expected_signals else None
        for field in interval_fields
    }
    internal = [record for record in records if record["manifest"]["dataset_id"] == "internal-68"]
    internal_without_ambiguity = [record for record in internal if record["manifest"]["evaluability_status"] != "LABEL_SEQUENCE_AMBIGUITY"]
    concordant = lambda record: record["decision"]["classification"].startswith("LABEL_CONCORDANT")
    result["internal_label_concordance_including_hpv18_ambiguity"] = wilson(sum(concordant(record) for record in internal), len(internal))
    result["internal_label_concordance_excluding_hpv18_ambiguity"] = wilson(sum(concordant(record) for record in internal_without_ambiguity), len(internal_without_ambiguity))
    result["external_by_cohort"] = {}
    for cohort in ("gse147507_sars", "gse147507_rsv", "gse91065_hpv"):
        cohort_records = [record for record in records if record["manifest"]["cohort"] == cohort]
        cohort_expected = [record for record in cohort_records if record["manifest"]["expectation_kind"] == "EXPECTED_GROUP"]
        cohort_mocks = [record for record in cohort_records if record["manifest"]["expectation_kind"] == "MOCK"]
        result["external_by_cohort"][cohort] = {
            "expected_group_observed": wilson(sum(record["decision"]["expected_observed"] for record in cohort_expected), len(cohort_expected)),
            "mock_unexpected_signal": wilson(sum(record["decision"]["observed_signal_count"] > 0 for record in cohort_mocks), len(cohort_mocks)),
            "wrong_group_resolved_signal": wilson(sum(record["decision"]["wrong_group_resolved"] for record in cohort_records), len(cohort_records)),
            "ambiguous_attribution": wilson(sum(record["decision"]["ambiguous_attribution"] for record in cohort_records), len(cohort_records)),
        }
    return result


def performance_metrics(records, config, digest_ledger, starts, terminals, internal_labels):
    paired = []
    for record in records:
        run = record["manifest"]
        old_result_path = historical_path(run, config)
        old_perf_path = historical_perf_path(run, config)
        require(old_result_path.is_file() and old_perf_path.is_file(), f"historical evaluation artifact missing: {run['run_id']}")
        old_result = load_json(old_result_path)
        old = old_perf(load_json(old_perf_path))
        new = record["parsed"]["perf"]
        old_observed = old_expected_observed(run, old_result, internal_labels)
        paired.append({
            "dataset_id": run["dataset_id"], "cohort": run["cohort"], "run_id": run["run_id"],
            "old_expected_observed": old_observed,
            "new_expected_observed": record["decision"]["expected_observed"],
            "expected_observed_difference": int(record["decision"]["expected_observed"]) - int(old_observed),
            "old_wall_time_ms": old["wall_time_ms"], "new_wall_time_ms": new["wall_time_ms"],
            "wall_ratio": new["wall_time_ms"] / old["wall_time_ms"],
            "old_peak_rss_bytes": old["peak_rss_bytes"], "new_peak_rss_bytes": new["peak_rss_bytes"],
            "rss_ratio": new["peak_rss_bytes"] / old["peak_rss_bytes"],
            "new_pass1_ms": new["pass1_count"], "new_pass2_bloom_alignment_ms": new["pass2_sample_prescreen_align"],
            "old_index_load_ms": old["stages"].get("index_load", 0) + old["stages"].get("index_open_mmi", 0),
            "old_sample_prescreen_alignment_ms": old["stages"].get("sample_prescreen_audit", 0),
        })
        perf = record["parsed"]["perf"]
        compressed = run["r1"]["compressed_bytes"] + (run["r2"]["compressed_bytes"] if run.get("r2") else 0)
        record["throughput"] = {
            "pass1_compressed_bytes_per_second": compressed / (perf["pass1_count"] / 1000) if perf["pass1_count"] else None,
            "pass2_combined_compressed_bytes_per_second": compressed / (perf["pass2_sample_prescreen_align"] / 1000) if perf["pass2_sample_prescreen_align"] else None,
            "bloom_pass_rate": perf["prescreen_passed_fragments"] / perf["selected_fragments"],
            "aligned_fragments": perf["aligned_fragments"],
        }
    result = {"indexes": {}, "datasets": {}, "cohorts": {}}
    for index_id, index in digest_ledger["indexes"].items():
        result["indexes"][index_id] = index["build"]
    for grouping, key_name in (("datasets", "dataset_id"), ("cohorts", "cohort")):
        keys = sorted({row[key_name] for row in paired})
        for key in keys:
            rows = [row for row in paired if row[key_name] == key]
            corresponding_old_max = max(row["old_peak_rss_bytes"] for row in rows)
            result[grouping][key] = {
                "sample_wall_time_ms": summary([row["new_wall_time_ms"] for row in rows]),
                "peak_rss_bytes": summary([row["new_peak_rss_bytes"] for row in rows]),
                "pass1_compressed_bytes_per_second": summary([record["throughput"]["pass1_compressed_bytes_per_second"] for record in records if record["manifest"][key_name] == key and record["throughput"]["pass1_compressed_bytes_per_second"] is not None]),
                "pass2_combined_compressed_bytes_per_second": summary([record["throughput"]["pass2_combined_compressed_bytes_per_second"] for record in records if record["manifest"][key_name] == key and record["throughput"]["pass2_combined_compressed_bytes_per_second"] is not None]),
                "bloom_pass_rate": summary([record["throughput"]["bloom_pass_rate"] for record in records if record["manifest"][key_name] == key]),
                "aligned_fragments": summary([record["throughput"]["aligned_fragments"] for record in records if record["manifest"][key_name] == key]),
                "paired_wall_ratio": summary([row["wall_ratio"] for row in rows]),
                "paired_rss_ratio": summary([row["rss_ratio"] for row in rows]),
                "historical_peak_rss_max_bytes": corresponding_old_max,
                "new_samples_over_historical_max": [row["run_id"] for row in rows if row["new_peak_rss_bytes"] > corresponding_old_max],
                "systematic_wall_regression": statistics.median(row["wall_ratio"] for row in rows) > 1.0,
                "individual_wall_regressions": [
                    {"run_id": row["run_id"], "wall_ratio": row["wall_ratio"], "new_pass1_ms": row["new_pass1_ms"], "new_pass2_bloom_alignment_ms": row["new_pass2_bloom_alignment_ms"], "old_index_load_ms": row["old_index_load_ms"], "old_sample_prescreen_alignment_ms": row["old_sample_prescreen_alignment_ms"]}
                    for row in rows if row["wall_ratio"] > 1.0
                ],
            }
            if grouping == "cohorts":
                cohort_runs = [record["manifest"] for record in records if record["manifest"]["cohort"] == key]
                cohort_keys = [(run["dataset_id"], run["run_id"]) for run in cohort_runs]
                started = min(datetime.fromisoformat(starts[run_key]["at"]) for run_key in cohort_keys)
                finished = max(datetime.fromisoformat(terminals[run_key]["at"]) for run_key in cohort_keys)
                result[grouping][key]["batch_wall_time_ms"] = round((finished - started).total_seconds() * 1000)
                result[grouping][key]["configured_jobs"] = cohort_runs[0]["planned_jobs"]
                result[grouping][key]["configured_threads_per_run"] = cohort_runs[0]["planned_threads"]
    return result, paired


def write_adjudication(path, records):
    fields = ("dataset_id", "cohort", "run_id", "sample_id", "truth_label", "expected_group_key", "evaluability_status", "classification", "ambiguity_notes", "input_fragments", "selected_fragments", "input_digest", "evidence_summary")
    lines = ["\t".join(fields)]
    for record in records:
        run = record["manifest"]
        csv_run = record["parsed"]["run"]
        values = {
            **run,
            "classification": record["decision"]["classification"],
            "ambiguity_notes": json.dumps(run["ambiguity_notes"], ensure_ascii=False, separators=(",", ":")),
            "input_fragments": csv_run["input_fragments"], "selected_fragments": csv_run["selected_fragments"], "input_digest": csv_run["input_digest"],
            "evidence_summary": json.dumps(record["decision"]["evidence_summary"], ensure_ascii=False, separators=(",", ":")),
        }
        lines.append("\t".join("" if values.get(field) is None else str(values[field]).replace("\t", " ").replace("\n", " ") for field in fields))
    atomic_text(path, "\n".join(lines) + "\n")


def write_run_table(path, records, failures):
    fields = ("dataset_id", "cohort", "run_id", "terminal_status", "classification", "expected_observed", "observed_signal_count", "wrong_group_resolved", "ambiguous_attribution", "input_fragments", "selected_fragments", "prescreen_passed_fragments", "aligned_fragments", "wall_time_ms", "peak_rss_bytes", "pass1_bytes_per_second", "pass2_bytes_per_second", "bloom_pass_rate")
    lines = ["\t".join(fields)]
    for record in records:
        run = record["manifest"]
        decision = record["decision"]
        perf = record["parsed"]["perf"]
        throughput = record["throughput"]
        values = {
            **run, "terminal_status": "SUCCESS", **decision, **perf,
            "pass1_bytes_per_second": throughput["pass1_compressed_bytes_per_second"],
            "pass2_bytes_per_second": throughput["pass2_combined_compressed_bytes_per_second"],
            "bloom_pass_rate": throughput["bloom_pass_rate"],
        }
        lines.append("\t".join(str(values[field]) for field in fields))
    for failure in failures:
        lines.append("\t".join(str(failure.get(field, failure["terminal"].get("status", "") if field == "terminal_status" else "")) for field in fields))
    atomic_text(path, "\n".join(lines) + "\n")


def write_paired(path, paired):
    fields = tuple(paired[0]) if paired else ("dataset_id", "cohort", "run_id")
    lines = ["\t".join(fields)]
    for row in paired:
        lines.append("\t".join(json.dumps(row[field], separators=(",", ":")) if isinstance(row[field], (dict, list)) else str(row[field]) for field in fields))
    atomic_text(path, "\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Frozen Phase 6 real-E2E campaign scorer")
    parser.add_argument("--campaign-dir", required=True)
    args = parser.parse_args()
    score(args.campaign_dir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"Phase 6 scorer failed: {error}")
