#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import os
from html.parser import HTMLParser
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO / "evaluation/phase0/evaluation-manifest.json"
PROFILE_PATH = REPO / "evaluation/phase0/analysis-profile.json"
CONFIG_PATH = Path(__file__).resolve().with_name("phase6_config.json")
BINARY_PATH = REPO / "target/release/viroflash"
INDEX_FILES = {"manifest.json", "ref.mmi", "bloom.bin", "reference-groups.tsv", "reference.fa"}
REPORT_FILES = {"report.csv", "report.html", "perf.json"}
REPORT_HEADER = (
    "schema_id", "record_type", "sample_id", "analysis_status", "reason_codes",
    "input_mode", "input_fragments", "selected_fragments", "selection_probability",
    "minimum_relevant_fraction", "familywise_miss_probability", "interval_level",
    "target_family_size", "prescreen_passed_fragments", "aligned_fragments",
    "unassigned_fragments", "profile_digest", "index_digest", "input_digest",
    "read_ends_per_fragment", "target_group_id", "representative_id", "member_ids",
    "evidence_status", "attribution_status", "supporting_selected_fragments",
    "selected_fragment_denominator", "attributed_fragment_fraction", "interval_lower",
    "interval_upper", "interval_method", "estimated_input_supporting_fragments",
    "covered_bases", "representative_length", "coverage_fraction", "occupied_windows",
    "host_confounded_fragments", "cross_group_ambiguous_fragments", "integration_status",
    "split_events", "discordant_fragments", "limitation_codes",
)
RUN_FIELDS = REPORT_HEADER[:20]
TARGET_FIELDS = (
    "record_type", "target_group_id", "representative_id", "member_ids",
    "evidence_status", "attribution_status", "supporting_selected_fragments",
    "selected_fragment_denominator", "attributed_fragment_fraction", "interval_lower",
    "interval_upper", "interval_level", "interval_method",
    "estimated_input_supporting_fragments", "covered_bases", "representative_length",
    "coverage_fraction", "occupied_windows", "host_confounded_fragments",
    "cross_group_ambiguous_fragments", "integration_status", "split_events",
    "discordant_fragments", "limitation_codes",
)
INTEGER_FIELDS = {
    "input_fragments", "selected_fragments", "target_family_size",
    "prescreen_passed_fragments", "aligned_fragments", "unassigned_fragments",
    "read_ends_per_fragment", "supporting_selected_fragments",
    "selected_fragment_denominator", "covered_bases", "representative_length",
    "occupied_windows", "host_confounded_fragments",
    "cross_group_ambiguous_fragments", "split_events", "discordant_fragments",
}
FLOAT_FIELDS = {
    "selection_probability", "minimum_relevant_fraction", "familywise_miss_probability",
    "interval_level", "attributed_fragment_fraction", "interval_lower", "interval_upper",
    "estimated_input_supporting_fragments", "coverage_fraction",
}
EMPTY_ALLOWED_FIELDS = {"reason_codes", "limitation_codes"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    with Path(path).open(encoding="utf-8") as handle:
        return json.load(handle)


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.part.{os.getpid()}")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=False)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    temporary.rename(path)


def atomic_text(path, value):
    path = Path(path)
    temporary = path.with_name(f".{path.name}.part.{os.getpid()}")
    with temporary.open("x", encoding="utf-8", newline="") as handle:
        handle.write(value)
        handle.flush()
        os.fsync(handle.fileno())
    temporary.rename(path)


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def is_sha256(value):
    return isinstance(value, str) and len(value) == 64 and all(c in "0123456789abcdef" for c in value)


def load_contract():
    manifest = load_json(MANIFEST_PATH)
    config = load_json(CONFIG_PATH)
    require(manifest["manifest_schema"] == "viroflash.phase0.evaluation-manifest", "unexpected evaluation manifest schema")
    require(manifest["freeze_state"] == "PRE_EXECUTION_FROZEN", "evaluation manifest is not pre-execution frozen")
    require(len(manifest["runs"]) == 127, "evaluation manifest does not contain 127 runs")
    require(len({(run["dataset_id"], run["run_id"]) for run in manifest["runs"]}) == 127, "manifest run identities are not unique")
    require({run["cohort"] for run in manifest["runs"]} == set(config["index_by_cohort"]), "cohort/index mapping is incomplete")
    require(all(len({run["planned_threads"] for run in manifest["runs"] if run["cohort"] == cohort}) == 1 for cohort in config["index_by_cohort"]), "threads vary within a cohort")
    require(all(len({run["planned_jobs"] for run in manifest["runs"] if run["cohort"] == cohort}) == 1 for cohort in config["index_by_cohort"]), "jobs vary within a cohort")
    repro_ids = [item["run_id"] for item in config["reproducibility_runs"]]
    require(len(repro_ids) == len(set(repro_ids)) and set(repro_ids) <= {run["run_id"] for run in manifest["runs"]}, "reproducibility run set is invalid")
    required_coverage = {
        "internal-largest", "internal-smallest", "EBV", "HBV", "HPV16",
        "HPV18-ambiguity", "negative", "SARS-positive", "SARS-mock",
        "RSV-positive", "RSV-mock", "external-HPV16", "external-HPV18",
        "SRR5090635", "plain-multimember-gzip",
    }
    require(set().union(*(set(item["covers"]) for item in config["reproducibility_runs"])) == required_coverage, "reproducibility coverage differs from frozen contract")
    require(config["reproducibility_threads"] == [1, 2, 4, 8] and config["reproducibility_repetitions"] == 2, "reproducibility matrix differs from frozen contract")
    require(sha256_file(PROFILE_PATH) == manifest["profile_digest"], "profile file digest differs from manifest")
    require(all(run["profile_digest"] == manifest["profile_digest"] for run in manifest["runs"]), "run profile digest differs from manifest")
    return manifest, config


def index_id_for(run, config):
    return config["index_by_cohort"][run["cohort"]]


def output_dir(campaign_dir, run):
    return Path(campaign_dir) / "runs" / run["dataset_id"] / run["run_id"]


def sample_id_from_path(path):
    name = Path(path).name
    for suffix in (".gz", ".fastq", ".fq", "_R1", "_1"):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


class ReportHtmlParser(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.section = None
        self.field = None
        self.capture = False
        self.value = []
        self.run = {}
        self.targets = []
        self.current_target = None

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if tag == "section":
            classes = set(attributes.get("class", "").split())
            if attributes.get("id") == "run-integrity":
                self.section = "run"
            elif "target-signal" in classes:
                self.section = "target"
                self.current_target = {"data_group": attributes.get("data-group", "")}
                self.targets.append(self.current_target)
        elif tag == "tr" and self.section in {"run", "target"}:
            self.field = attributes.get("data-field")
        elif tag == "td" and self.field:
            self.capture = True
            self.value = []

    def handle_data(self, data):
        if self.capture:
            self.value.append(data)

    def handle_endtag(self, tag):
        if tag == "td" and self.capture:
            value = "".join(self.value)
            destination = self.run if self.section == "run" else self.current_target
            require(self.field not in destination, f"duplicate HTML field: {self.field}")
            destination[self.field] = value
            self.capture = False
        elif tag == "tr":
            self.field = None
        elif tag == "section" and self.section in {"run", "target"}:
            self.section = None
            self.current_target = None


def parse_report_csv(path):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require(tuple(reader.fieldnames or ()) == REPORT_HEADER, "report.csv header differs from frozen contract")
        rows = list(reader)
    require(rows and rows[0]["record_type"] == "RUN", "report.csv must start with one RUN row")
    require(sum(row["record_type"] == "RUN" for row in rows) == 1, "report.csv RUN row count is not one")
    require(all(row["record_type"] == "TARGET_SIGNAL" for row in rows[1:]), "report.csv contains an unknown record type")
    targets = rows[1:]
    require(len({row["target_group_id"] for row in targets}) == len(targets), "duplicate target group in report.csv")
    for row in rows:
        fields = RUN_FIELDS if row["record_type"] == "RUN" else TARGET_FIELDS
        for field in fields:
            require(row[field] != "" or field in EMPTY_ALLOWED_FIELDS, f"empty required CSV field: {field}")
        for field in INTEGER_FIELDS & set(fields):
            require(int(row[field]) >= 0, f"invalid non-negative integer field: {field}")
        for field in FLOAT_FIELDS & set(fields):
            require(math.isfinite(float(row[field])), f"invalid finite float field: {field}")
    return rows[0], targets


def parse_report_html(path):
    parser = ReportHtmlParser()
    parser.feed(Path(path).read_text(encoding="utf-8"))
    parser.close()
    require(set(parser.run) == set(RUN_FIELDS), "HTML RUN fields differ from CSV-visible contract")
    for target in parser.targets:
        require(set(target) == set(TARGET_FIELDS) | {"data_group"}, "HTML target fields differ from CSV-visible contract")
        require(target["data_group"] == target["target_group_id"], "HTML target data-group mismatch")
    return parser.run, parser.targets


def parse_perf_json(path):
    perf = load_json(path)
    required = {
        "schema_id", "status", "sample_id", "wall_time_ms", "configured_threads",
        "process_cpu_time_ms", "peak_rss_bytes", "read_bytes", "written_bytes",
        "pass1_count", "pass2_sample_prescreen_align", "report_write", "input_fragments",
        "selected_fragments", "prescreen_passed_fragments", "aligned_fragments",
        "telemetry_status", "sample_errors",
    }
    require(set(perf) == required, "perf.json fields differ from frozen contract")
    require(perf["schema_id"] == "viroflash.perf.v1", "unexpected perf.json schema")
    require(perf["status"] in {"SUCCESS", "ERROR"}, "invalid perf terminal status")
    for field in required - {"schema_id", "status", "sample_id", "telemetry_status", "sample_errors"}:
        require(isinstance(perf[field], int) and perf[field] >= 0, f"invalid perf integer: {field}")
    require(isinstance(perf["sample_errors"], list), "perf sample_errors is not a list")
    return perf


def validate_success_directory(directory, run, expected_profile_digest, expected_index_digest):
    directory = Path(directory)
    require(directory.is_dir(), f"missing report directory: {directory}")
    require({entry.name for entry in directory.iterdir()} == REPORT_FILES, f"successful directory does not contain exactly three reports: {directory}")
    csv_run, csv_targets = parse_report_csv(directory / "report.csv")
    html_run, html_targets = parse_report_html(directory / "report.html")
    perf = parse_perf_json(directory / "perf.json")
    require(all(html_run[field] == csv_run[field] for field in RUN_FIELDS), "HTML and CSV RUN fields differ")
    html_by_group = {target["target_group_id"]: target for target in html_targets}
    require(len(html_by_group) == len(html_targets) == len(csv_targets), "HTML and CSV target counts differ")
    for target in csv_targets:
        html_target = html_by_group.get(target["target_group_id"])
        require(html_target is not None, f"target absent from HTML: {target['target_group_id']}")
        require(all(html_target[field] == target[field] for field in TARGET_FIELDS), f"HTML and CSV target fields differ: {target['target_group_id']}")
    expected_sample = sample_id_from_path(run["r1"]["path"])
    require(csv_run["sample_id"] == perf["sample_id"] == expected_sample, "reported sample identity mismatch")
    require(csv_run["analysis_status"] in {"CONFORMANT_COMPLETE", "CONFORMANT_WITH_LIMITATIONS"}, "invalid analysis status")
    require(csv_run["input_mode"] == run["input_mode"], "reported input mode mismatch")
    require(int(csv_run["read_ends_per_fragment"]) == (2 if run["input_mode"] == "PE" else 1), "read-end contract mismatch")
    require(csv_run["profile_digest"] == expected_profile_digest, "reported profile digest mismatch")
    require(csv_run["index_digest"] == expected_index_digest, "reported index digest mismatch")
    require(is_sha256(csv_run["input_digest"]), "reported decoded input digest is invalid")
    require(perf["status"] == "SUCCESS" and not perf["sample_errors"], "success report has error-shaped perf status")
    require(perf["configured_threads"] == run["planned_threads"], "reported thread count differs from frozen manifest")
    for field in ("input_fragments", "selected_fragments", "prescreen_passed_fragments", "aligned_fragments"):
        require(perf[field] == int(csv_run[field]), f"perf/CSV count mismatch: {field}")
    require(int(csv_run["input_fragments"]) > 0 and int(csv_run["selected_fragments"]) > 0, "successful run has no input or selected fragments")
    require(int(csv_run["prescreen_passed_fragments"]) <= int(csv_run["selected_fragments"]), "prescreen count exceeds selected count")
    require(int(csv_run["aligned_fragments"]) <= int(csv_run["prescreen_passed_fragments"]), "aligned count exceeds prescreen count")
    return {"run": csv_run, "targets": csv_targets, "perf": perf}


def validate_error_directory(directory):
    directory = Path(directory)
    require(directory.is_dir(), f"missing error directory: {directory}")
    names = {entry.name for entry in directory.iterdir()}
    require(names == {"perf.json"}, f"error directory is success-shaped or contains extra files: {directory}")
    perf = parse_perf_json(directory / "perf.json")
    require(perf["status"] == "ERROR" and perf["sample_errors"], "error perf.json lacks terminal error evidence")
    return perf


def verify_index_directory(index_dir, expected):
    index_dir = Path(index_dir)
    names = {entry.name for entry in index_dir.iterdir() if entry.is_file()}
    require(names == INDEX_FILES, f"index artifacts differ from frozen set: {index_dir}")
    actual = {name: sha256_file(index_dir / name) for name in sorted(names)}
    require(actual == expected["artifacts"], f"index artifact digest changed: {index_dir}")
    manifest = load_json(index_dir / "manifest.json")
    require(manifest["profile_digest"] == expected["profile_digest"], "index profile digest mismatch")
    require(manifest["host_fasta_sha256"] == expected["host_fasta_sha256"], "index host digest mismatch")
    require(manifest["target_fasta_sha256"] == expected["target_fasta_sha256"], "index target digest mismatch")
    require(manifest["mmi_digest"] == actual["ref.mmi"], "embedded MMI digest mismatch")
    require(manifest["bloom_digest"] == actual["bloom.bin"], "embedded Bloom digest mismatch")
    require(manifest["ledger_digest"] == actual["reference-groups.tsv"], "embedded ledger digest mismatch")
    return actual


def load_digest_ledger(campaign_dir):
    ledger = load_json(Path(campaign_dir) / "digest-ledger.json")
    require(ledger["binary"]["path"] == str(BINARY_PATH), "campaign binary path is not target/release/viroflash")
    return ledger


def verify_campaign_digests(campaign_dir, include_index_artifacts=True):
    ledger = load_digest_ledger(campaign_dir)
    require(sha256_file(BINARY_PATH) == ledger["binary"]["sha256"], "campaign binary digest changed")
    require(sha256_file(PROFILE_PATH) == ledger["profile"]["sha256"], "campaign profile digest changed")
    if include_index_artifacts:
        for index_id, expected in ledger["indexes"].items():
            verify_index_directory(Path(campaign_dir) / "indexes" / index_id, expected)
    return ledger


def parse_time_verbose(path):
    values = {}
    for line in Path(path).read_text(encoding="utf-8").splitlines():
        if ": " in line:
            key, value = line.strip().rsplit(": ", 1)
            values[key] = value
    elapsed = values.get("Elapsed (wall clock) time (h:mm:ss or m:ss)")
    require(elapsed is not None, f"missing elapsed time in {path}")
    parts = [float(part) for part in elapsed.split(":")]
    seconds = sum(value * (60 ** power) for power, value in enumerate(reversed(parts)))
    return {
        "wall_time_ms": round(seconds * 1000),
        "peak_rss_bytes": int(values["Maximum resident set size (kbytes)"]) * 1024,
        "exit_status": int(values["Exit status"]),
    }
