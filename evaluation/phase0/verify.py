#!/usr/bin/env python3
import csv
import hashlib
import json
import tomllib
from collections import Counter
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
PHASE0 = Path(__file__).resolve().parent
MANIFEST = PHASE0 / "evaluation-manifest.json"
EXTERNAL_FASTQ_SUMS = Path("/home/wubw/data/viroflash-benchmarks/metadata/FASTQ_SHA256SUMS")
EXTERNAL_REFERENCE_SUMS = Path("/home/wubw/data/viroflash-benchmarks/refs/SHA256SUMS")
EXTERNAL_RESULT_SUMS = Path("/home/wubw/data/viroflash-benchmarks/results_20260827/SHA256SUMS")
INTERNAL_RUN_CONFIG = Path("/home/wubw/data/viroflash/dna_virus_68_20260827/run_config.json")
SHA256_LENGTH = 64

RUN_FIELDS = {
    "dataset_id", "cohort", "run_id", "sample_id", "input_mode", "r1", "r2",
    "host_reference", "target_reference", "truth_label", "expectation_kind",
    "expected_group_key", "label_provenance", "evaluability_status", "ambiguity_notes",
    "profile_digest", "binary_digest", "index_digest", "planned_threads", "planned_jobs",
}
FILE_FIELDS = {"path", "compressed_bytes", "sha256", "checksum_status"}
PROFILE_VALUES = {
    "minimum_relevant_fraction": ("δ", 0.00001),
    "familywise_interval_error": ("α", 0.05),
    "familywise_miss_probability": ("β", 0.05),
    "maximum_interval_width": ("w", 0.00001),
}
OUTPUT_FIELDS = {
    ("report.csv", "RUN"): {
        "schema_id", "record_type", "sample_id", "analysis_status", "reason_codes",
        "input_mode", "input_fragments", "selected_fragments", "selection_probability",
        "minimum_relevant_fraction", "familywise_miss_probability", "interval_level",
        "target_family_size", "prescreen_passed_fragments", "aligned_fragments",
        "unassigned_fragments", "profile_digest", "index_digest", "input_digest",
        "read_ends_per_fragment",
    },
    ("report.csv", "TARGET_SIGNAL"): {
        "target_group_id", "representative_id", "member_ids", "evidence_status",
        "attribution_status", "quantitation_status", "supporting_selected_fragments",
        "selected_fragment_denominator", "attributed_fragment_fraction", "interval_lower",
        "interval_upper", "interval_level", "interval_method",
        "estimated_input_supporting_fragments", "covered_bases", "representative_length",
        "coverage_fraction", "occupied_windows", "host_confounded_fragments",
        "cross_group_ambiguous_fragments", "integration_status", "split_events",
        "discordant_fragments", "limitation_codes",
    },
    ("report.html", "DOCUMENT"): {
        "interpretation_boundary", "run_integrity", "observed_target_signals",
        "evidence_detail", "methods_and_limitations",
    },
    ("perf.json", "RUN"): {
        "schema_id", "status", "sample_id", "wall_time_ms", "configured_threads",
        "process_cpu_time_ms", "peak_rss_bytes", "read_bytes", "written_bytes",
        "pass1_count", "pass2_sample_prescreen_align", "report_write", "input_fragments",
        "selected_fragments", "prescreen_passed_fragments", "aligned_fragments",
        "telemetry_status", "sample_errors",
    },
}


def require(condition, message):
    if not condition:
        raise SystemExit(f"Phase 0 verification failed: {message}")


def is_sha256(value):
    return isinstance(value, str) and len(value) == SHA256_LENGTH and all(
        character in "0123456789abcdef" for character in value
    )


def small_file_sha256(path):
    require(path.stat().st_size <= 10_000_000, f"refusing to hash non-metadata file: {path}")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_checksum_list(path):
    require(path.is_file(), f"missing checksum list: {path}")
    records = {}
    for line in path.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        require(is_sha256(digest), f"invalid checksum in {path}: {digest}")
        records[str(Path(name))] = digest
    return records


def validate_file(record, context, authoritative_checksums):
    require(isinstance(record, dict) and set(record) == FILE_FIELDS, f"{context}: field mismatch")
    path = Path(record["path"])
    require(path.is_absolute() and path.is_file(), f"{context}: missing file {path}")
    require(path.stat().st_size == record["compressed_bytes"], f"{context}: size mismatch")
    status = record["checksum_status"]
    if status == "MISSING_AUTHORITATIVE_CHECKSUM":
        require(record["sha256"] is None, f"{context}: checksum gap must be null")
        return
    require(status in {"AUTHORITATIVE_RECORDED", "PHASE0_COMPUTED"}, f"{context}: invalid checksum status")
    require(is_sha256(record["sha256"]), f"{context}: malformed checksum")
    if status == "AUTHORITATIVE_RECORDED":
        require(authoritative_checksums.get(str(path)) == record["sha256"], f"{context}: checksum-list mismatch")


def load_truth_mapping():
    with (PHASE0 / "truth-mapping.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(row["dataset_id"], row["cohort"], row["source_label"]) for row in rows]
    require(len(keys) == len(set(keys)), "truth mapping has duplicate keys")
    return dict(zip(keys, rows))


def validate_profile():
    profile_path = PHASE0 / "analysis-profile.json"
    profile = json.loads(profile_path.read_text())
    require(set(profile) == {"contract_id", "decision_status", "frozen_on", "index_and_alignment", "parameters"}, "profile fields changed")
    require(profile["decision_status"] == "FROZEN_BEFORE_FIRST_V0_5_RESULT", "profile is not frozen")
    require(set(profile["parameters"]) == set(PROFILE_VALUES), "profile parameter set changed")
    for field, (symbol, value) in PROFILE_VALUES.items():
        require(profile["parameters"][field] == {"symbol": symbol, "value": value}, f"{field} changed")
    require(profile["index_and_alignment"]["alignment_role_margin"] == 0, "extra alignment margin introduced")
    return small_file_sha256(profile_path)


def validate_manifest(profile_digest):
    manifest = json.loads(MANIFEST.read_text())
    require(set(manifest) == {"manifest_schema", "frozen_on", "profile_digest", "runs"}, "manifest fields changed")
    require(manifest["manifest_schema"] == "viroflash.phase0.evaluation-manifest", "manifest identity changed")
    require(manifest["profile_digest"] == profile_digest, "profile digest mismatch")
    runs = manifest["runs"]
    require(len(runs) == 127, "manifest must contain 127 runs")
    require(Counter(run.get("dataset_id") for run in runs) == {"internal-68": 68, "external-59": 59}, "cohort counts differ from 68+59")
    run_ids = [run.get("run_id") for run in runs]
    require(len(run_ids) == len(set(run_ids)), "run IDs are not unique")

    fastq_sums = read_checksum_list(EXTERNAL_FASTQ_SUMS)
    reference_sums = read_checksum_list(EXTERNAL_REFERENCE_SUMS)
    result_sums = read_checksum_list(EXTERNAL_RESULT_SUMS)
    authoritative = fastq_sums | reference_sums
    truth_mapping = load_truth_mapping()
    checksum_statuses = Counter()

    for run in runs:
        run_id = run.get("run_id", "<missing>")
        require(set(run) == RUN_FIELDS, f"{run_id}: run fields changed")
        require(run["input_mode"] in {"SE", "PE"}, f"{run_id}: invalid input mode")
        require((run["r2"] is None) == (run["input_mode"] == "SE"), f"{run_id}: input-mode mismatch")
        require(run["profile_digest"] == profile_digest, f"{run_id}: profile digest mismatch")
        for digest_name in ("binary_digest", "index_digest"):
            require(run[digest_name] is None or is_sha256(run[digest_name]), f"{run_id}: malformed {digest_name}")
        mapping = truth_mapping.get((run["dataset_id"], run["cohort"], run["truth_label"]))
        require(mapping is not None, f"{run_id}: truth mapping missing")
        require(run["expectation_kind"] == mapping["expectation_kind"], f"{run_id}: expectation mismatch")
        mapped_group = "" if mapping["expected_group_key"] == "-" else mapping["expected_group_key"]
        require((run["expected_group_key"] or "") == mapped_group, f"{run_id}: expected group mismatch")

        for field in ("r1", "r2", "host_reference", "target_reference"):
            record = run[field]
            if record is None:
                continue
            validate_file(record, f"{run_id}.{field}", authoritative)
            checksum_statuses[(run["dataset_id"], field, record["checksum_status"])] += 1

        provenance = run["label_provenance"]
        require(set(provenance) == {"path", "sha256"}, f"{run_id}: label provenance fields changed")
        path = Path(provenance["path"])
        require(path.is_file() and is_sha256(provenance["sha256"]), f"{run_id}: invalid label provenance")
        if run["dataset_id"] == "external-59":
            require(result_sums.get(str(path)) == provenance["sha256"], f"{run_id}: label checksum-list mismatch")

    expected_statuses = Counter({
        ("internal-68", "r1", "MISSING_AUTHORITATIVE_CHECKSUM"): 68,
        ("internal-68", "r2", "MISSING_AUTHORITATIVE_CHECKSUM"): 68,
        ("internal-68", "host_reference", "PHASE0_COMPUTED"): 68,
        ("internal-68", "target_reference", "PHASE0_COMPUTED"): 68,
        ("external-59", "r1", "AUTHORITATIVE_RECORDED"): 59,
        ("external-59", "r2", "AUTHORITATIVE_RECORDED"): 23,
        ("external-59", "host_reference", "AUTHORITATIVE_RECORDED"): 59,
        ("external-59", "target_reference", "AUTHORITATIVE_RECORDED"): 59,
    })
    require(checksum_statuses == expected_statuses, "checksum coverage/status matrix changed")

    run_config = json.loads(INTERNAL_RUN_CONFIG.read_text())
    internal = next(run for run in runs if run["dataset_id"] == "internal-68")
    require(internal["target_reference"]["sha256"] == run_config["target_fasta_sha256"], "internal target checksum metadata mismatch")
    require(internal["label_provenance"]["sha256"] == run_config["truth_sha256"], "internal truth checksum metadata mismatch")


def validate_contract_files():
    with (PHASE0 / "output-fields.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(row["artifact"], row["record_type"], row["field"]) for row in rows]
    require(len(keys) == len(set(keys)), "output field contract has duplicates")
    observed = {}
    for row in rows:
        observed.setdefault((row["artifact"], row["record_type"]), set()).add(row["field"])
    require(observed == OUTPUT_FIELDS, "output field contract changed")
    for name, count in (("internal-68.tsv", 68), ("external-59.tsv", 59)):
        with (PHASE0 / "historical" / name).open(newline="") as handle:
            require(sum(1 for _ in csv.DictReader(handle, delimiter="\t")) == count, f"{name}: row count changed")
    summary = json.loads((PHASE0 / "historical" / "summary.json").read_text())
    require(summary["interpretation"] == "Historical observations only; not v0.5 acceptance truth.", "historical boundary changed")
    require(summary["internal_68"]["source_sha256"] == small_file_sha256(Path(summary["internal_68"]["source_path"])), "internal baseline source mismatch")
    require(summary["external_59"]["source_sha256"] == small_file_sha256(Path(summary["external_59"]["source_path"])), "external baseline source mismatch")
    with (REPO / "Cargo.toml").open("rb") as handle:
        require(tomllib.load(handle)["package"]["version"] == "0.3.0", "Cargo version changed")


def main():
    profile_digest = validate_profile()
    validate_manifest(profile_digest)
    validate_contract_files()
    print("Phase 0 verification passed: 68 internal + 59 external runs; no FASTQ content read")


if __name__ == "__main__":
    main()
