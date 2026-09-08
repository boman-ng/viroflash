#!/usr/bin/env python3
import csv
import hashlib
import json
import sys
import tomllib
from collections import Counter, defaultdict
from pathlib import Path

from freeze_reference_groups import ALPHABET, COMPLEMENT, LEDGER_FIELDS, canonical_sequence, records


REPO = Path(__file__).resolve().parents[2]
PHASE0 = Path(__file__).resolve().parent
MANIFEST = PHASE0 / "evaluation-manifest.json"
EXTERNAL_FASTQ_SUMS = Path("/home/wubw/data/viroflash-benchmarks/metadata/FASTQ_SHA256SUMS")
EXTERNAL_REFERENCE_SUMS = Path("/home/wubw/data/viroflash-benchmarks/refs/SHA256SUMS")
EXTERNAL_RESULT_SUMS = Path("/home/wubw/data/viroflash-benchmarks/results_20260827/SHA256SUMS")
EXTERNAL_LABELS = Path("/home/wubw/data/viroflash-benchmarks/results_20260827/label_comparison.tsv")
EXTERNAL_COHORTS = Path("/home/wubw/data/viroflash-benchmarks/metadata/cohorts.tsv")
HPV_RUN_LABELS = Path("/home/wubw/data/viroflash-benchmarks/metadata/GSE91065_run_labels.tsv")
INTERNAL_RUN_CONFIG = Path("/home/wubw/data/viroflash/dna_virus_68_20260827/run_config.json")
INTERNAL_TRUTH = Path("/home/wubw/data/viroflash/dna_virus_68_20260827/panel_truth.json")
INTERNAL_CHECKSUM_STATUS = "PHASE0_FROZEN_COMPRESSED_SHA256"
SHA256_LENGTH = 64

RUN_FIELDS = {
    "dataset_id", "cohort", "run_id", "sample_id", "input_mode", "r1", "r2",
    "host_reference", "target_reference", "truth_label", "expectation_kind",
    "expected_group_key", "label_provenance", "evaluability_status", "ambiguity_notes",
    "profile_digest", "planned_threads", "planned_jobs",
}
FILE_FIELDS = {"path", "compressed_bytes", "sha256", "checksum_status"}
PROFILE_PARAMETERS = {
    "minimum_relevant_fraction": {"symbol": "δ", "value": 0.00001},
    "familywise_interval_error": {"symbol": "α", "value": 0.05},
    "familywise_miss_probability": {"symbol": "β", "value": 0.05},
    "maximum_interval_width": {"symbol": "w", "value": 0.00001},
}
INDEX_AND_ALIGNMENT = {
    "kmer_length": 21,
    "minimap2_preset": "sr",
    "reference_group_equivalence": "uppercase full sequence equals another member or its reverse complement",
    "fasta_record_id": "first whitespace-delimited header token, UTF-8 and unique within the file",
    "sequence_normalization": "remove ASCII whitespace between FASTA headers, then uppercase",
    "iupac_dna_alphabet": "ACGTMRWSYKVHDBN",
    "invalid_sequence_symbol": "reject the reference file without dropping or replacing the symbol",
    "alignment_acceptance": "all emitted chains with CIGAR and positive query length",
    "alignment_role_margin": 0,
    "occupied_window_bins": 10,
}
REFERENCE_GROUP_SOURCE_IDS = {"internal-dna-panel", "external-respiratory-panel", "external-hpv-panel"}
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
        "record_type", "target_group_id", "representative_id", "member_ids",
        "evidence_status", "attribution_status", "quantitation_status",
        "supporting_selected_fragments", "selected_fragment_denominator",
        "attributed_fragment_fraction", "interval_lower", "interval_upper", "interval_level",
        "interval_method", "estimated_input_supporting_fragments", "covered_bases",
        "representative_length", "coverage_fraction", "occupied_windows",
        "host_confounded_fragments", "cross_group_ambiguous_fragments", "integration_status",
        "split_events", "discordant_fragments", "limitation_codes",
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


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def frozen_evaluation_digest(manifest):
    frozen = {
        "profile_digest": manifest["profile_digest"],
        "checksum_provenance": manifest["checksum_provenance"],
        "reference_groups_sha256": manifest["reference_groups_sha256"],
        "runs": manifest["runs"],
    }
    encoded = json.dumps(frozen, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()
    return hashlib.sha256(encoded).hexdigest()


def read_checksum_list(path):
    require(path.is_file(), f"missing checksum list: {path}")
    checksums = {}
    for line in path.read_text().splitlines():
        digest, name = line.split(maxsplit=1)
        require(is_sha256(digest), f"invalid checksum in {path}: {digest}")
        require(name not in checksums, f"duplicate checksum path in {path}: {name}")
        checksums[name] = digest
    return checksums


def validate_profile():
    profile_path = PHASE0 / "analysis-profile.json"
    profile = json.loads(profile_path.read_text())
    require(set(profile) == {"contract_id", "decision_status", "frozen_on", "index_and_alignment", "parameters"}, "profile fields changed")
    require(profile["decision_status"] == "PENDING_PRODUCT_DECISION", "profile decision status changed")
    require(profile["parameters"] == PROFILE_PARAMETERS, "profile parameter contract changed")
    require(profile["index_and_alignment"] == INDEX_AND_ALIGNMENT, "index/alignment profile changed")
    require(ALPHABET == frozenset(b"ACGTMRWSYKVHDBN"), "parser IUPAC alphabet changed")
    require(b"ACGTMRWSYKVHDBN".translate(COMPLEMENT) == b"TGCAKYWSRMBDHVN", "parser IUPAC complement changed")
    return file_sha256(profile_path)


def validate_internal_checksum_freeze(runs):
    checksum_path = PHASE0 / "internal-fastq-sha256.tsv"
    with checksum_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    require(rows and tuple(rows[0]) == ("run_id", "read_end", "path", "compressed_bytes", "sha256", "provenance"), "internal checksum fields changed")
    require(len(rows) == 136, "internal checksum list must contain 136 files")
    row_keys = [(row["run_id"], row["read_end"]) for row in rows]
    require(len(row_keys) == len(set(row_keys)), "internal checksum run/read-end key duplicated")
    frozen = {(row["run_id"], row["read_end"]): row for row in rows}
    for run in runs:
        if run["dataset_id"] != "internal-68":
            continue
        for read_end in ("r1", "r2"):
            row = frozen.get((run["run_id"], read_end))
            require(row is not None, f"{run['run_id']}.{read_end}: frozen checksum missing")
            record = run[read_end]
            require(row["path"] == record["path"], f"{run['run_id']}.{read_end}: frozen path mismatch")
            require(int(row["compressed_bytes"]) == record["compressed_bytes"], f"{run['run_id']}.{read_end}: frozen size mismatch")
            require(row["sha256"] == record["sha256"], f"{run['run_id']}.{read_end}: frozen digest mismatch")
            require(row["provenance"] == INTERNAL_CHECKSUM_STATUS == record["checksum_status"], f"{run['run_id']}.{read_end}: provenance mismatch")


def validate_truth_sources(runs):
    by_key = {(run["dataset_id"], run["run_id"]): run for run in runs}
    internal_truth = json.loads(INTERNAL_TRUTH.read_text())
    require(len(internal_truth) == 68, "panel_truth.json must contain 68 samples")
    for truth in internal_truth:
        run_id = str(truth["data_id"])
        run = by_key.get(("internal-68", run_id))
        require(run is not None, f"internal truth run missing: {run_id}")
        require(run["sample_id"] == run_id and run["cohort"] == "internal_dna_virus_68", f"{run_id}: internal identity mismatch")
        require(run["input_mode"] == "PE", f"{run_id}: internal mode must be PE")
        require(Path(run["r1"]["path"]).name == truth["r1"], f"{run_id}: R1 does not match panel truth")
        require(Path(run["r2"]["path"]).name == truth["r2"], f"{run_id}: R2 does not match panel truth")
        require(run["truth_label"] == truth["positive"], f"{run_id}: truth label mismatch")
        require(run["expectation_kind"] == ("MOCK" if truth["positive"] is None else "EXPECTED_GROUP"), f"{run_id}: expectation kind mismatch")

    with EXTERNAL_LABELS.open(newline="") as handle:
        labels = list(csv.DictReader(handle, delimiter="\t"))
    require(len(labels) == 59, "label_comparison.tsv must contain 59 runs")
    with EXTERNAL_COHORTS.open(newline="") as handle:
        respiratory_metadata = {row["run_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with HPV_RUN_LABELS.open(newline="") as handle:
        hpv_metadata = {row["run_accession"]: row for row in csv.DictReader(handle, delimiter="\t")}
    for label in labels:
        run_id = label["run"]
        run = by_key.get(("external-59", run_id))
        require(run is not None, f"external label run missing: {run_id}")
        require(run["cohort"] == label["cohort"], f"{run_id}: cohort mismatch")
        require(run["truth_label"] == label["ground_truth_label"], f"{run_id}: truth label mismatch")
        expected_group = None if label["expected_target"] == "none" else label["expected_target"]
        require(run["expected_group_key"] == expected_group, f"{run_id}: expected target mismatch")
        if run["cohort"] == "gse91065_hpv":
            metadata = hpv_metadata.get(run_id)
            require(metadata is not None and run["sample_id"] == metadata["logical_sample"], f"{run_id}: HPV run/sample mapping mismatch")
        else:
            metadata = respiratory_metadata.get(run_id)
            require(metadata is not None and metadata["cohort"] == run["cohort"], f"{run_id}: respiratory run/cohort mapping mismatch")
            require(metadata["label"] == run["truth_label"] and run["sample_id"] == run_id, f"{run_id}: respiratory sample mapping mismatch")
        if run["input_mode"] == "SE":
            require(Path(run["r1"]["path"]).name == f"{run_id}.fastq.gz" and run["r2"] is None, f"{run_id}: SE FASTQ mismatch")
        else:
            require(Path(run["r1"]["path"]).name == f"{run_id}_1.fastq.gz", f"{run_id}: PE R1 mismatch")
            require(Path(run["r2"]["path"]).name == f"{run_id}_2.fastq.gz", f"{run_id}: PE R2 mismatch")

    with (PHASE0 / "truth-mapping.tsv").open(newline="") as handle:
        mappings = list(csv.DictReader(handle, delimiter="\t"))
    mapping_keys = []
    for mapping in mappings:
        source_label = None if mapping["source_label"] == "<NULL>" else mapping["source_label"]
        mapping_keys.append((mapping["dataset_id"], mapping["cohort"], source_label))
    require(len(mapping_keys) == len(set(mapping_keys)), "truth mapping keys are not unique")
    mapping_by_key = dict(zip(mapping_keys, mappings))
    run_mapping_keys = {(run["dataset_id"], run["cohort"], run["truth_label"]) for run in runs}
    require(set(mapping_by_key) == run_mapping_keys, "truth mapping contains missing or unused labels")
    for run in runs:
        mapping = mapping_by_key.get((run["dataset_id"], run["cohort"], run["truth_label"]))
        require(mapping is not None, f"{run['run_id']}: truth-mapping row missing")
        mapped_group = None if mapping["expected_group_key"] == "-" else mapping["expected_group_key"]
        require(mapping["expectation_kind"] == run["expectation_kind"], f"{run['run_id']}: mapped expectation mismatch")
        require(mapped_group == run["expected_group_key"], f"{run['run_id']}: mapped group mismatch")


def validate_reference_groups(profile_digest, manifest):
    group_manifest_path = PHASE0 / "reference-groups.json"
    require(file_sha256(group_manifest_path) == manifest["reference_groups_sha256"], "reference-group manifest digest mismatch")
    group_manifest = json.loads(group_manifest_path.read_text())
    require(set(group_manifest) == {"contract_id", "profile_digest", "sources"}, "reference-group manifest fields changed")
    require(group_manifest["contract_id"] == "viroflash.phase0.reference-groups", "reference-group manifest identity changed")
    require(group_manifest["profile_digest"] == profile_digest, "reference-group profile digest mismatch")
    require({source["source_id"] for source in group_manifest["sources"]} == REFERENCE_GROUP_SOURCE_IDS, "reference-group source set changed")
    external_reference_sums = read_checksum_list(EXTERNAL_REFERENCE_SUMS)
    internal_target_digest = json.loads(INTERNAL_RUN_CONFIG.read_text())["target_fasta_sha256"]

    for source in group_manifest["sources"]:
        require(set(source) == {"source_id", "target_fasta", "target_fasta_sha256", "ledger", "ledger_sha256", "record_count", "group_count"}, f"{source.get('source_id')}: reference-group source fields changed")
        fasta_path = Path(source["target_fasta"])
        expected_fasta_digest = (
            internal_target_digest
            if source["source_id"] == "internal-dna-panel"
            else external_reference_sums.get(str(fasta_path))
        )
        require(source["target_fasta_sha256"] == expected_fasta_digest, f"{source['source_id']}: target checksum metadata mismatch")
        ledger_path = PHASE0 / source["ledger"]
        require(file_sha256(ledger_path) == source["ledger_sha256"], f"{source['source_id']}: ledger digest mismatch")
        with ledger_path.open(newline="") as handle:
            ledger_rows = list(csv.DictReader(handle, delimiter="\t"))
        require(ledger_rows and tuple(ledger_rows[0]) == LEDGER_FIELDS, f"{source['source_id']}: ledger fields changed")
        require(len(ledger_rows) == source["record_count"], f"{source['source_id']}: ledger record count mismatch")
        ledger_by_member = {row["member_id"]: row for row in ledger_rows}
        require(len(ledger_by_member) == len(ledger_rows), f"{source['source_id']}: ledger member duplicated")
        ledger_group_sizes = Counter(row["target_group_id"] for row in ledger_rows)

        fasta_ids = set()
        observed_groups = defaultdict(list)
        duplicate_group_sequences = {}
        for member_id, sequence in records(fasta_path):
            require(member_id not in fasta_ids, f"{source['source_id']}: duplicate FASTA ID")
            fasta_ids.add(member_id)
            row = ledger_by_member.get(member_id)
            require(row is not None, f"{source['source_id']}: FASTA record absent from ledger: {member_id}")
            canonical = canonical_sequence(sequence)
            target_group_id = f"sha256:{hashlib.sha256(canonical).hexdigest()}"
            require(row["target_group_id"] == target_group_id, f"{source['source_id']}: group mismatch for {member_id}")
            if ledger_group_sizes[target_group_id] > 1:
                previous = duplicate_group_sequences.setdefault(target_group_id, canonical)
                require(previous == canonical, f"{source['source_id']}: non-identical records merged in {target_group_id}")
            require(int(row["representative_length"]) == len(sequence), f"{source['source_id']}: length mismatch for {member_id}")
            require(row["target_fasta_sha256"] == expected_fasta_digest, f"{source['source_id']}: row target digest mismatch")
            require(row["profile_digest"] == profile_digest, f"{source['source_id']}: row profile digest mismatch")
            observed_groups[target_group_id].append(member_id)
        require(fasta_ids == set(ledger_by_member), f"{source['source_id']}: ledger contains non-FASTA member")
        require(len(observed_groups) == source["group_count"], f"{source['source_id']}: group count mismatch")
        sorted_groups = sorted(observed_groups)
        for group_ordinal, target_group_id in enumerate(sorted_groups):
            members = sorted(observed_groups[target_group_id])
            rows = sorted((ledger_by_member[member] for member in members), key=lambda row: int(row["member_ordinal"]))
            require([int(row["member_ordinal"]) for row in rows] == list(range(len(rows))), f"{source['source_id']}: member ordinals invalid")
            require(all(int(row["group_ordinal"]) == group_ordinal for row in rows), f"{source['source_id']}: group ordinal invalid")
            require(all(row["representative_id"] == members[0] for row in rows), f"{source['source_id']}: representative invalid")


def validate_manifest(profile_digest):
    manifest = json.loads(MANIFEST.read_text())
    require(set(manifest) == {"manifest_schema", "frozen_on", "freeze_state", "profile_digest", "frozen_evaluation_digest", "checksum_provenance", "reference_groups_sha256", "runs"}, "manifest fields changed")
    require(manifest["manifest_schema"] == "viroflash.phase0.evaluation-manifest", "manifest identity changed")
    require(manifest["freeze_state"] == "PRE_EXECUTION_FROZEN", "manifest must remain pre-execution")
    require(manifest["profile_digest"] == profile_digest, "manifest profile digest mismatch")
    require(manifest["frozen_evaluation_digest"] == frozen_evaluation_digest(manifest), "frozen evaluation digest mismatch")
    require(manifest["checksum_provenance"] == {INTERNAL_CHECKSUM_STATUS: {"frozen_on": "2026-09-08", "method": "SHA-256 streamed over existing gzip file bytes without decompression", "list": "internal-fastq-sha256.tsv"}}, "internal checksum provenance changed")
    runs = manifest["runs"]
    require(len(runs) == 127, "manifest must contain 127 runs")
    require(Counter(run.get("dataset_id") for run in runs) == {"internal-68": 68, "external-59": 59}, "dataset counts differ from 68+59")
    run_ids = [run.get("run_id") for run in runs]
    require(len(run_ids) == len(set(run_ids)), "run IDs are not unique")
    require(Counter((run["dataset_id"], run["cohort"], run["input_mode"]) for run in runs) == Counter({("internal-68", "internal_dna_virus_68", "PE"): 68, ("external-59", "gse147507_rsv", "SE"): 6, ("external-59", "gse147507_sars", "SE"): 6, ("external-59", "gse91065_hpv", "SE"): 24, ("external-59", "gse91065_hpv", "PE"): 23}), "cohort/input-mode matrix changed")

    authoritative = read_checksum_list(EXTERNAL_FASTQ_SUMS) | read_checksum_list(EXTERNAL_REFERENCE_SUMS)
    result_sums = read_checksum_list(EXTERNAL_RESULT_SUMS)
    for run in runs:
        run_id = run.get("run_id", "<missing>")
        require(set(run) == RUN_FIELDS, f"{run_id}: run fields changed")
        require(run["profile_digest"] == profile_digest, f"{run_id}: profile identity mismatch")
        require((run["r2"] is None) == (run["input_mode"] == "SE"), f"{run_id}: input-mode mismatch")
        for field in ("r1", "r2", "host_reference", "target_reference"):
            record = run[field]
            if record is None:
                continue
            require(set(record) == FILE_FIELDS, f"{run_id}.{field}: fields changed")
            path = Path(record["path"])
            require(path.is_absolute() and path.is_file(), f"{run_id}.{field}: missing file")
            require(path.stat().st_size == record["compressed_bytes"], f"{run_id}.{field}: size mismatch")
            require(is_sha256(record["sha256"]), f"{run_id}.{field}: malformed checksum")
            if record["checksum_status"] == "AUTHORITATIVE_RECORDED":
                require(authoritative.get(str(path)) == record["sha256"], f"{run_id}.{field}: checksum-list mismatch")
            else:
                require(record["checksum_status"] in {"PHASE0_COMPUTED", INTERNAL_CHECKSUM_STATUS}, f"{run_id}.{field}: invalid checksum provenance")
        provenance = run["label_provenance"]
        require(set(provenance) == {"path", "sha256"} and is_sha256(provenance["sha256"]), f"{run_id}: label provenance invalid")
        if run["dataset_id"] == "external-59":
            require(result_sums.get(provenance["path"]) == provenance["sha256"], f"{run_id}: label checksum-list mismatch")

    validate_internal_checksum_freeze(runs)
    validate_truth_sources(runs)
    validate_reference_groups(profile_digest, manifest)


def validate_contract_files():
    with (PHASE0 / "output-fields.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    keys = [(row["artifact"], row["record_type"], row["field"]) for row in rows]
    require(len(keys) == len(set(keys)), "output field contract has duplicates")
    observed = defaultdict(set)
    for row in rows:
        observed[(row["artifact"], row["record_type"])].add(row["field"])
    require(dict(observed) == OUTPUT_FIELDS, "output field contract changed")
    for name, count in (("internal-68.tsv", 68), ("external-59.tsv", 59)):
        with (PHASE0 / "historical" / name).open(newline="") as handle:
            require(sum(1 for _ in csv.DictReader(handle, delimiter="\t")) == count, f"{name}: row count changed")
    summary = json.loads((PHASE0 / "historical" / "summary.json").read_text())
    require(summary["interpretation"] == "Historical observations only; not v0.5 acceptance truth.", "historical boundary changed")
    require(summary["internal_68"]["source_sha256"] == file_sha256(Path(summary["internal_68"]["source_path"])), "internal baseline source mismatch")
    require(summary["external_59"]["source_sha256"] == file_sha256(Path(summary["external_59"]["source_path"])), "external baseline source mismatch")
    with (REPO / "Cargo.toml").open("rb") as handle:
        require(tomllib.load(handle)["package"]["version"] == "0.3.0", "Cargo version changed")


def main():
    require(len(sys.argv) == 1, "this metadata verifier accepts no execution or E2E gate arguments")
    profile_digest = validate_profile()
    validate_manifest(profile_digest)
    validate_contract_files()
    print("Phase 0 checksum metadata consistency passed; FASTQ bytes were not read; analysis profile remains pending product decision")


if __name__ == "__main__":
    main()
