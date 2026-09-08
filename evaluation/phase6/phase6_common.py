#!/usr/bin/env python3
import csv
import hashlib
import json
import math
import os
import re
from functools import lru_cache
from html.parser import HTMLParser
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO / "evaluation/phase0/evaluation-manifest.json"
PROFILE_PATH = REPO / "evaluation/phase0/analysis-profile.json"
CONFIG_PATH = Path(__file__).resolve().with_name("phase6_config.json")
COMMON_PATH = Path(__file__).resolve()
SCORER_PATH = COMMON_PATH.with_name("phase6_scorer.py")
BINARY_PATH = REPO / "target/release/viroflash"
SCORING_PROVENANCE_NAME = "scoring-provenance.json"
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
INDEX_LEDGER_HEADER = (
    "group_ordinal", "target_group_id", "representative_id", "member_ordinal",
    "member_id", "representative_length", "target_fasta_sha256", "profile_digest",
)
SCORER_CORRECTIONS = (
    {
        "commit": "853677e",
        "timing": "PRE_SCORING",
        "discovery": "independent report audit",
        "correction": "Expected internal labels map to every matching exact-sequence ReferenceGroup; any observed matching row establishes expected-observed concordance.",
    },
    {
        "commit": "26aeab2",
        "timing": "POST_RUN_PRE_FINAL_SCORING",
        "discovery": "evaluation semantic review",
        "correction": "Internal adjudication is limited to the frozen EBV/HBV/HPV16/HPV18 label scope while preserving total observed evidence separately.",
    },
    {
        "commit": "e55dbd2",
        "timing": "POST_RUN_POST_INITIAL_SCORING_PRE_REVIEW_RESCORING",
        "discovery": "Phase 6 reviewer integrity review",
        "correction": "Internal taxon matching starts at the FASTA description and excludes non-human Heron hepatitis B virus; independent public report semantic validation and scoring provenance binding were added.",
    },
    {
        "commit": "5fc4efb",
        "timing": "POST_RUN_POST_REVIEW_RESCORING_PRE_SECOND_REVIEW_RESCORING",
        "discovery": "second Phase 6 reviewer integrity review",
        "correction": "A bounded metadata-aware species matcher restores legitimate scoped human-virus labels while excluding Heron HBV and HPV numeric lookalikes; exact enum/status validation and independent cached hypergeometric inversion expose invalid report intervals without suppressing evidence.",
    },
)
RETROSPECTIVE_PROVENANCE_LIMIT = (
    "The completed current campaign did not bind scorer/config/manifest inputs before execution. "
    "This ledger was created retrospectively before reviewer-requested rescoring and cannot establish pre-run scoring immutability."
)
EVIDENCE_STATUSES = {"REFERENCE_SIGNAL_OBSERVED", "INDETERMINATE_EVIDENCE"}
ATTRIBUTION_STATUSES = {
    "RESOLVED_TO_REFERENCE_GROUP", "AMBIGUOUS_WITHIN_GROUP",
    "UNRESOLVED_ACROSS_GROUPS", "CONFOUNDED_WITH_HOST",
}
INTEGRATION_STATUSES = {"NOT_OBSERVED", "DIAGNOSTIC_EVIDENCE_OBSERVED"}
RUN_LIMITATION = re.compile(r"TARGET_KMER_NOT_EVALUABLE_SELECTED_FRAGMENTS=([1-9][0-9]*)\Z")


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
    require(config["reproducibility_jobs"] == 8, "reproducibility job budget differs from frozen contract")
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


def close_float(actual, expected, absolute=1e-15, relative=1e-12):
    return math.isclose(actual, expected, abs_tol=absolute, rel_tol=relative)


def hypergeometric_tail_reaches(
        population, total_successes, sample, observed, upper_tail, threshold, inclusive=True):
    support_low = max(0, sample - (population - total_successes))
    support_high = min(sample, total_successes)
    if upper_tail:
        if observed <= support_low:
            return True
        if observed > support_high:
            return False
    else:
        if observed < support_low:
            return False
        if observed >= support_high:
            return True
    mode = ((sample + 1) * (total_successes + 1)) // (population + 2)
    if upper_tail and observed <= mode:
        return True
    if not upper_tail and observed >= mode:
        return True
    mode = min(max(mode, support_low), support_high)
    total_terms = [1.0]
    tail_terms = []
    weight = 1.0
    current = mode
    while current > support_low:
        numerator = current * (population - total_successes - sample + current)
        denominator = (total_successes - current + 1) * (sample - current + 1)
        require(numerator > 0 and denominator > 0, "invalid lower hypergeometric recurrence")
        weight *= numerator / denominator
        current -= 1
        if weight <= 1e-18:
            break
        total_terms.append(weight)
        if not upper_tail and current <= observed:
            tail_terms.append(weight)
    weight = 1.0
    current = mode
    while current < support_high:
        numerator = (total_successes - current) * (sample - current)
        current += 1
        denominator = current * (population - total_successes - sample + current)
        require(numerator > 0 and denominator > 0, "invalid upper hypergeometric recurrence")
        weight *= numerator / denominator
        if weight <= 1e-18:
            break
        total_terms.append(weight)
        if upper_tail and current >= observed:
            tail_terms.append(weight)
    total_weight = math.fsum(total_terms)
    tail_weight = math.fsum(tail_terms)
    scaled_threshold = threshold * total_weight
    return tail_weight >= scaled_threshold if inclusive else tail_weight > scaled_threshold


@lru_cache(maxsize=None)
def exact_hypergeometric_interval(population, sample, successes, level):
    require(0 < population and 0 < sample <= population and 0 <= successes <= sample, "invalid exact interval inputs")
    require(0 < level < 1, "invalid exact interval level")
    if sample == population:
        exact = successes / population
        return exact, exact
    tail = (1.0 - level) / 2.0
    feasible_low = successes
    feasible_high = population - (sample - successes)
    if successes == 0:
        lower = 0
    else:
        low, high = feasible_low, feasible_high
        while low < high:
            middle = low + (high - low) // 2
            if hypergeometric_tail_reaches(
                    population, middle, sample, successes, True, tail, inclusive=False):
                high = middle
            else:
                low = middle + 1
        lower = max(feasible_low, low - 1)
    if successes == sample:
        upper = population
    else:
        low, high = feasible_low, feasible_high
        while low < high:
            middle = low + (high - low + 1) // 2
            if hypergeometric_tail_reaches(population, middle, sample, successes, False, tail):
                low = middle
            else:
                high = middle - 1
        upper = min(feasible_high, low + 1)
    return lower / population, upper / population


def validate_report_invariants(csv_run, csv_targets, index_contract, interval_mismatches=None):
    input_fragments = int(csv_run["input_fragments"])
    selected_fragments = int(csv_run["selected_fragments"])
    prescreen_fragments = int(csv_run["prescreen_passed_fragments"])
    aligned_fragments = int(csv_run["aligned_fragments"])
    unassigned_fragments = int(csv_run["unassigned_fragments"])
    family_size = int(csv_run["target_family_size"])
    probability = float(csv_run["selection_probability"])
    minimum_fraction = float(csv_run["minimum_relevant_fraction"])
    miss_probability = float(csv_run["familywise_miss_probability"])
    family_interval_level = float(csv_run["interval_level"])
    require(family_size == index_contract["target_family_size"], "reported target family size differs from index ledger")
    require(0 < minimum_fraction <= 1 and 0 < miss_probability < 1, "reported sampling design probabilities are out of range")
    require(0 < probability <= 1, "reported selection probability is out of range")
    require(0 < selected_fragments <= input_fragments, "selected/input fragment counts are incoherent")
    require(unassigned_fragments <= aligned_fragments <= prescreen_fragments <= selected_fragments, "reported fragment counts are not ordered")
    minimum_relevant = math.ceil(minimum_fraction * input_fragments)
    require((1.0 - probability) ** minimum_relevant <= miss_probability / family_size, "selection probability exceeds the reported familywise miss budget")
    if probability == 1.0:
        require(selected_fragments == input_fragments, "census selection probability did not select every fragment")
    require(0 < family_interval_level < 1, "familywise interval level is out of range")
    expected_target_level = 1.0 - (1.0 - family_interval_level) / family_size
    previous_ordinal = -1
    for target in csv_targets:
        group = index_contract["groups"].get(target["target_group_id"])
        require(group is not None, f"reported target group is absent from index ledger: {target['target_group_id']}")
        require(group["ordinal"] > previous_ordinal, "reported target groups are not in index-ledger order")
        previous_ordinal = group["ordinal"]
        require(target["representative_id"] == group["representative_id"], "reported target representative differs from index ledger")
        require(target["member_ids"].split(";") == group["member_ids"], "reported target members differ from index ledger")
        require(int(target["representative_length"]) == group["representative_length"], "reported representative length differs from index ledger")
        supporting = int(target["supporting_selected_fragments"])
        denominator = int(target["selected_fragment_denominator"])
        attributed_fraction = float(target["attributed_fragment_fraction"])
        lower = float(target["interval_lower"])
        upper = float(target["interval_upper"])
        covered_bases = int(target["covered_bases"])
        representative_length = int(target["representative_length"])
        require(target["evidence_status"] in EVIDENCE_STATUSES, "unexpected target evidence status")
        require(target["attribution_status"] in ATTRIBUTION_STATUSES, "unexpected target attribution status")
        require(target["integration_status"] in INTEGRATION_STATUSES, "unexpected target integration status")
        require(denominator == selected_fragments and supporting <= denominator, "target support fraction denominator/count is incoherent")
        require(close_float(attributed_fraction, supporting / denominator), "target attributed fraction arithmetic mismatch")
        require(target["interval_method"] == "EQUAL_TAILED_EXACT_HYPERGEOMETRIC_INVERSION", "unexpected target interval method")
        require(close_float(float(target["interval_level"]), expected_target_level), "target interval level does not implement familywise allocation")
        require(0 <= lower <= attributed_fraction <= upper <= 1, "target interval bounds are unordered or out of range")
        expected_lower, expected_upper = exact_hypergeometric_interval(
            input_fragments, denominator, supporting, expected_target_level,
        )
        lower_matches = close_float(lower, expected_lower)
        upper_matches = close_float(upper, expected_upper)
        if interval_mismatches is not None and not (lower_matches and upper_matches):
            interval_mismatches.append({
                "target_group_id": target["target_group_id"],
                "input_fragments": input_fragments,
                "selected_fragments": denominator,
                "supporting_selected_fragments": supporting,
                "interval_level": expected_target_level,
                "reported_lower": lower,
                "reported_upper": upper,
                "expected_lower": expected_lower,
                "expected_upper": expected_upper,
                "reported_lower_population_count": round(lower * input_fragments),
                "reported_upper_population_count": round(upper * input_fragments),
                "expected_lower_population_count": round(expected_lower * input_fragments),
                "expected_upper_population_count": round(expected_upper * input_fragments),
            })
        else:
            require(lower_matches, "target interval lower endpoint differs from exact inversion")
            require(upper_matches, "target interval upper endpoint differs from exact inversion")
        require(covered_bases <= representative_length, "covered bases exceed representative length")
        require(close_float(float(target["coverage_fraction"]), covered_bases / representative_length), "target coverage fraction arithmetic mismatch")
        require(close_float(float(target["estimated_input_supporting_fragments"]), attributed_fraction * input_fragments, absolute=1e-9, relative=1e-11), "estimated input support arithmetic mismatch")
        require(int(target["occupied_windows"]) <= representative_length, "occupied windows exceed representative length")
        require(int(target["host_confounded_fragments"]) <= selected_fragments, "host-confounded count exceeds selected fragments")
        require(int(target["cross_group_ambiguous_fragments"]) <= selected_fragments, "cross-group ambiguous count exceeds selected fragments")
        require(int(target["split_events"]) <= supporting, "split-event count exceeds supporting fragments")
        require(int(target["discordant_fragments"]) <= supporting, "discordant count exceeds supporting fragments")
        observed = target["evidence_status"] == "REFERENCE_SIGNAL_OBSERVED"
        require(observed == (supporting > 0), "target evidence status contradicts supporting count")
        if observed:
            expected_attribution = "RESOLVED_TO_REFERENCE_GROUP" if len(group["member_ids"]) == 1 else "AMBIGUOUS_WITHIN_GROUP"
            require(target["attribution_status"] == expected_attribution, "observed target has incoherent attribution status")
        else:
            require(int(target["host_confounded_fragments"]) > 0 or int(target["cross_group_ambiguous_fragments"]) > 0, "indeterminate target has no indeterminate evidence")
            expected_attribution = "UNRESOLVED_ACROSS_GROUPS" if int(target["cross_group_ambiguous_fragments"]) > 0 else "CONFOUNDED_WITH_HOST"
            require(target["attribution_status"] == expected_attribution, "indeterminate target has incoherent attribution status")
        integration_observed = supporting > 0 and (
            int(target["split_events"]) > 0 or int(target["discordant_fragments"]) > 0
        )
        expected_integration = "DIAGNOSTIC_EVIDENCE_OBSERVED" if integration_observed else "NOT_OBSERVED"
        require(target["integration_status"] == expected_integration, "target integration status contradicts diagnostic counts")
        require(not target["limitation_codes"], "target limitation codes are not defined by the report contract")
    has_run_limitations = bool(csv_run["reason_codes"])
    require((csv_run["analysis_status"] == "CONFORMANT_WITH_LIMITATIONS") == has_run_limitations, "analysis status and run limitations are incoherent")
    require((csv_run["analysis_status"] == "CONFORMANT_COMPLETE") == (not has_run_limitations), "complete analysis status has limitations")
    if has_run_limitations:
        match = RUN_LIMITATION.fullmatch(csv_run["reason_codes"])
        require(match is not None, "unknown or malformed run limitation code")
        require(int(match.group(1)) <= selected_fragments, "run limitation count exceeds selected fragments")


def validate_success_directory(
    directory, run, expected_profile_digest, expected_index_digest, index_contract,
    interval_mismatches=None,
):
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
    validate_report_invariants(csv_run, csv_targets, index_contract, interval_mismatches)
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


def load_index_report_contract(index_dir, expected):
    index_dir = Path(index_dir)
    manifest = load_json(index_dir / "manifest.json")
    require(manifest["profile_digest"] == expected["profile_digest"], "index/report profile digest mismatch")
    require(manifest["target_fasta_sha256"] == expected["target_fasta_sha256"], "index/report target digest mismatch")
    groups = {}
    with (index_dir / "reference-groups.tsv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(tuple(reader.fieldnames or ()) == INDEX_LEDGER_HEADER, "reference-group ledger header differs from contract")
        for row in reader:
            group_ordinal = int(row["group_ordinal"])
            member_ordinal = int(row["member_ordinal"])
            group = groups.get(row["target_group_id"])
            if group is None:
                require(group_ordinal == len(groups), "reference groups are not contiguous and ordered")
                require(member_ordinal == 0, "reference-group member ordinals do not start at zero")
                group = {
                    "ordinal": group_ordinal,
                    "representative_id": row["representative_id"],
                    "representative_length": int(row["representative_length"]),
                    "member_ids": [],
                }
                groups[row["target_group_id"]] = group
            require(group["ordinal"] == group_ordinal, "target group appears in multiple ledger ordinals")
            require(group["representative_id"] == row["representative_id"], "group representative changes within ledger")
            require(group["representative_length"] == int(row["representative_length"]), "group representative length changes within ledger")
            require(member_ordinal == len(group["member_ids"]), "reference-group members are not contiguous and ordered")
            require(row["target_fasta_sha256"] == expected["target_fasta_sha256"], "ledger target digest mismatch")
            require(row["profile_digest"] == expected["profile_digest"], "ledger profile digest mismatch")
            group["member_ids"].append(row["member_id"])
    require(groups, "reference-group ledger contains no groups")
    return {"target_family_size": len(groups), "groups": groups}


def load_digest_ledger(campaign_dir):
    ledger = load_json(Path(campaign_dir) / "digest-ledger.json")
    require(ledger["binary"]["path"] == str(BINARY_PATH), "campaign binary path is not target/release/viroflash")
    return ledger


def verify_manifest_identity(manifest_path, prepared):
    manifest_path = Path(manifest_path)
    manifest = load_json(manifest_path)
    require(prepared["path"] == str(manifest_path), "campaign manifest path changed")
    require(sha256_file(manifest_path) == prepared["sha256"], "campaign manifest digest changed")
    require(manifest["frozen_evaluation_digest"] == prepared["frozen_evaluation_digest"], "campaign frozen evaluation identity changed")
    return manifest


def verify_prepared_digests(campaign_dir, include_index_artifacts=True):
    ledger = load_digest_ledger(campaign_dir)
    verify_manifest_identity(MANIFEST_PATH, ledger["manifest"])
    require(sha256_file(BINARY_PATH) == ledger["binary"]["sha256"], "campaign binary digest changed")
    require(sha256_file(PROFILE_PATH) == ledger["profile"]["sha256"], "campaign profile digest changed")
    if include_index_artifacts:
        for index_id, expected in ledger["indexes"].items():
            verify_index_directory(Path(campaign_dir) / "indexes" / index_id, expected)
    return ledger


def historical_path(run, config):
    if run["dataset_id"] == "internal-68":
        return Path(config["historical"]["internal_results"]) / f"{run['run_id']}.json"
    return Path(config["historical"]["external_results"]) / run["cohort"] / f"{run['run_id']}.json"


def historical_perf_path(run, config):
    return historical_path(run, config).with_suffix(".perf.json")


def digest_record(path):
    path = Path(path)
    require(path.is_file(), f"scoring input is missing: {path}")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def scoring_input_records(manifest, config):
    references = {}
    for run in manifest["runs"]:
        for role in ("host_reference", "target_reference"):
            record = run[role]
            key = (role, record["path"])
            if key in references:
                continue
            actual = digest_record(record["path"])
            require(actual["bytes"] == record["compressed_bytes"], f"frozen reference byte size changed: {record['path']}")
            require(actual["sha256"] == record["sha256"], f"frozen reference digest changed: {record['path']}")
            references[key] = {"role": role, **actual, "manifest_sha256": record["sha256"]}
    historical = []
    for run in manifest["runs"]:
        for kind, path in (("result", historical_path(run, config)), ("performance", historical_perf_path(run, config))):
            historical.append({"dataset_id": run["dataset_id"], "cohort": run["cohort"], "run_id": run["run_id"], "kind": kind, **digest_record(path)})
    historical.append({"kind": "internal_index_performance", **digest_record(config["historical"]["internal_index_perf"])})
    return {
        "evaluation_sources": [digest_record(COMMON_PATH), digest_record(SCORER_PATH)],
        "config": digest_record(CONFIG_PATH),
        "manifest": {**digest_record(MANIFEST_PATH), "frozen_evaluation_digest": manifest["frozen_evaluation_digest"]},
        "references": sorted(references.values(), key=lambda item: (item["role"], item["path"])),
        "historical": historical,
    }


def create_scoring_provenance(campaign_dir, manifest, config, mode, created_at):
    campaign_dir = Path(campaign_dir)
    destination = campaign_dir / SCORING_PROVENANCE_NAME
    require(not destination.exists(), "scoring provenance ledger already exists")
    verify_prepared_digests(campaign_dir)
    require(mode in {"PRE_EXECUTION_BOUND", "RETROSPECTIVE_PRE_RESCORING"}, "invalid scoring provenance mode")
    value = {
        "schema_id": "viroflash.phase6.scoring-provenance.v1",
        "contract_id": config["contract_id"],
        "mode": mode,
        "created_at": created_at,
        "inputs": scoring_input_records(manifest, config),
        "prior_evaluation_only_corrections": list(SCORER_CORRECTIONS),
        "provenance_limitation": None if mode == "PRE_EXECUTION_BOUND" else RETROSPECTIVE_PROVENANCE_LIMIT,
    }
    atomic_json(destination, value)
    return value


def verify_scoring_provenance(campaign_dir, manifest, config):
    path = Path(campaign_dir) / SCORING_PROVENANCE_NAME
    require(path.is_file(), "campaign scoring provenance ledger is missing")
    recorded = load_json(path)
    require(recorded["schema_id"] == "viroflash.phase6.scoring-provenance.v1", "unexpected scoring provenance schema")
    require(recorded["contract_id"] == config["contract_id"], "scoring provenance contract changed")
    require(recorded["mode"] in {"PRE_EXECUTION_BOUND", "RETROSPECTIVE_PRE_RESCORING"}, "invalid scoring provenance mode")
    require(recorded["prior_evaluation_only_corrections"] == list(SCORER_CORRECTIONS), "scoring correction disclosure changed")
    expected_limit = None if recorded["mode"] == "PRE_EXECUTION_BOUND" else RETROSPECTIVE_PROVENANCE_LIMIT
    require(recorded["provenance_limitation"] == expected_limit, "scoring provenance limitation changed")
    require(recorded["inputs"] == scoring_input_records(manifest, config), "scoring input digest ledger changed")
    return recorded


def verify_campaign_digests(campaign_dir, include_index_artifacts=True):
    manifest, config = load_contract()
    ledger = verify_prepared_digests(campaign_dir, include_index_artifacts)
    verify_scoring_provenance(campaign_dir, manifest, config)
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
