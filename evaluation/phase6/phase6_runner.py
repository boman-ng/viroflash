#!/usr/bin/env python3
import argparse
import concurrent.futures
import gzip
import json
import os
import shutil
import subprocess
import sys
import threading
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

from phase6_common import (
    BINARY_PATH, MANIFEST_PATH, PROFILE_PATH, atomic_json, atomic_text, index_id_for,
    load_contract, load_digest_ledger, load_json, output_dir, parse_time_verbose, require,
    sha256_file, validate_error_directory, validate_success_directory,
    verify_campaign_digests, verify_index_directory,
)


TIME = Path("/usr/bin/time")
LEDGER_LOCK = threading.Lock()


def now():
    return datetime.now(timezone.utc).isoformat()


def append_event(path, event):
    line = json.dumps(event, sort_keys=True, ensure_ascii=False) + "\n"
    with LEDGER_LOCK:
        with Path(path).open("a", encoding="utf-8") as handle:
            handle.write(line)
            handle.flush()
            os.fsync(handle.fileno())


def run_timed(command, stdout_path, stderr_path, time_path):
    full_command = [str(TIME), "-v", "-o", str(time_path), *map(str, command)]
    with Path(stdout_path).open("x", encoding="utf-8") as stdout, Path(stderr_path).open("x", encoding="utf-8") as stderr:
        completed = subprocess.run(full_command, stdout=stdout, stderr=stderr, check=False)
    timing = parse_time_verbose(time_path)
    require(timing["exit_status"] == completed.returncode, "time/subprocess exit status mismatch")
    return completed.returncode, timing


def reference_sets(manifest, config):
    indexes = {}
    for run in manifest["runs"]:
        index_id = index_id_for(run, config)
        value = {
            "host": run["host_reference"],
            "target": run["target_reference"],
            "threads": run["planned_threads"],
        }
        if index_id in indexes:
            require(indexes[index_id]["host"] == value["host"] and indexes[index_id]["target"] == value["target"], f"index source drift in {index_id}")
        else:
            indexes[index_id] = value
    return indexes


def prepare(campaign_dir):
    manifest, config = load_contract()
    campaign_dir = Path(campaign_dir).resolve()
    require(not campaign_dir.exists(), f"campaign path already exists: {campaign_dir}")
    require(BINARY_PATH.is_file() and os.access(BINARY_PATH, os.X_OK), f"missing executable release binary: {BINARY_PATH}")
    campaign_dir.mkdir(parents=True)
    (campaign_dir / "indexes").mkdir()
    (campaign_dir / "logs/indexes").mkdir(parents=True)
    binary_digest = sha256_file(BINARY_PATH)
    profile_digest = sha256_file(PROFILE_PATH)
    require(profile_digest == manifest["profile_digest"], "actual profile digest differs from frozen manifest")
    indexes = {}
    for index_id, source in reference_sets(manifest, config).items():
        host = Path(source["host"]["path"])
        target = Path(source["target"]["path"])
        require(host.stat().st_size == source["host"]["compressed_bytes"], f"host size mismatch: {host}")
        require(target.stat().st_size == source["target"]["compressed_bytes"], f"target size mismatch: {target}")
        host_digest = sha256_file(host)
        target_digest = sha256_file(target)
        require(host_digest == source["host"]["sha256"], f"host digest mismatch: {host}")
        require(target_digest == source["target"]["sha256"], f"target digest mismatch: {target}")
        index_dir = campaign_dir / "indexes" / index_id
        log_dir = campaign_dir / "logs/indexes"
        command = [BINARY_PATH, "index", "--host-fa", host, "--target-fa", target, "--out", index_dir, "--threads", str(source["threads"])]
        returncode, timing = run_timed(command, log_dir / f"{index_id}.stdout", log_dir / f"{index_id}.stderr", log_dir / f"{index_id}.time")
        require(returncode == 0, f"index build failed: {index_id}")
        artifacts = {entry.name: sha256_file(entry) for entry in sorted(index_dir.iterdir()) if entry.is_file()}
        index_manifest = load_json(index_dir / "manifest.json")
        entry = {
            "path": str(index_dir),
            "profile_digest": index_manifest["profile_digest"],
            "index_digest": artifacts["manifest.json"],
            "host_fasta_path": str(host),
            "host_fasta_sha256": host_digest,
            "target_fasta_path": str(target),
            "target_fasta_sha256": target_digest,
            "artifacts": artifacts,
            "build": timing,
            "threads": source["threads"],
        }
        require(entry["profile_digest"] == profile_digest, f"index profile digest mismatch: {index_id}")
        indexes[index_id] = entry
        verify_index_directory(index_dir, entry)
    ledger = {
        "contract_id": config["contract_id"],
        "prepared_at": now(),
        "manifest": {"path": str(MANIFEST_PATH), "sha256": sha256_file(MANIFEST_PATH), "frozen_evaluation_digest": manifest["frozen_evaluation_digest"]},
        "binary": {"path": str(BINARY_PATH), "sha256": binary_digest},
        "profile": {"path": str(PROFILE_PATH), "sha256": profile_digest},
        "indexes": indexes,
    }
    atomic_json(campaign_dir / "digest-ledger.json", ledger)
    verify_campaign_digests(campaign_dir)
    print(f"prepared {len(indexes)} indexes in {campaign_dir}")


def verify_fastq(run):
    result = {"dataset_id": run["dataset_id"], "run_id": run["run_id"], "files": {}}
    for read_end in ("r1", "r2"):
        record = run.get(read_end)
        if record is None:
            continue
        path = Path(record["path"])
        require(path.is_file(), f"missing FASTQ: {path}")
        require(path.stat().st_size == record["compressed_bytes"], f"FASTQ size mismatch: {path}")
        digest = sha256_file(path)
        require(digest == record["sha256"], f"FASTQ digest mismatch: {path}")
        result["files"][read_end] = {"path": str(path), "bytes": path.stat().st_size, "sha256": digest}
    require((run["input_mode"] == "PE") == ("r2" in result["files"]), f"input mode/path mismatch: {run['run_id']}")
    return result


def preflight_fastqs(campaign_dir, manifest):
    destination = Path(campaign_dir) / "fastq-digest-ledger.json"
    require(not destination.exists(), "FASTQ digest ledger already exists")
    print("verifying all 127 frozen FASTQ inputs before first launch", flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        rows = list(executor.map(verify_fastq, manifest["runs"]))
    atomic_json(destination, {"verified_at": now(), "runs": rows})


def run_one(campaign_dir, run, config, digest_ledger, ledger_path):
    index_id = index_id_for(run, config)
    index = digest_ledger["indexes"][index_id]
    directory = output_dir(campaign_dir, run)
    log_dir = Path(campaign_dir) / "logs/runs" / run["dataset_id"]
    log_dir.mkdir(parents=True, exist_ok=True)
    directory.parent.mkdir(parents=True, exist_ok=True)
    require(not directory.exists(), f"run output already exists: {directory}")
    require(sha256_file(BINARY_PATH) == digest_ledger["binary"]["sha256"], "binary changed before sample launch")
    command = [BINARY_PATH, "run", "--r1", run["r1"]["path"]]
    if run.get("r2"):
        command.extend(["--r2", run["r2"]["path"]])
    command.extend(["--index", index["path"], "--out", directory, "--threads", str(run["planned_threads"])])
    identity = {"dataset_id": run["dataset_id"], "cohort": run["cohort"], "run_id": run["run_id"], "attempt": 1}
    append_event(ledger_path, {**identity, "event": "STARTED", "at": now(), "command": [str(value) for value in command]})
    try:
        returncode, timing = run_timed(command, log_dir / f"{run['run_id']}.stdout", log_dir / f"{run['run_id']}.stderr", log_dir / f"{run['run_id']}.time")
        if returncode == 0:
            parsed = validate_success_directory(directory, run, digest_ledger["profile"]["sha256"], index["index_digest"])
            status = "SUCCESS"
            summary = {
                "analysis_status": parsed["run"]["analysis_status"],
                "input_fragments": int(parsed["run"]["input_fragments"]),
                "selected_fragments": int(parsed["run"]["selected_fragments"]),
                "target_signal_rows": len(parsed["targets"]),
                "input_digest": parsed["run"]["input_digest"],
                "pass1_pass2_contract": "BINARY_VERIFIED_COUNTS_IDS_AND_DECODED_INPUT_DIGEST",
            }
        else:
            perf = validate_error_directory(directory)
            status = "ERROR"
            summary = {"sample_errors": perf["sample_errors"]}
    except Exception as error:
        status = "HARNESS_ERROR"
        summary = {"error": f"{type(error).__name__}: {error}"}
        returncode = -1
        timing = {}
    append_event(ledger_path, {**identity, "event": "TERMINAL", "at": now(), "status": status, "returncode": returncode, "time": timing, "summary": summary})
    print(f"{run['dataset_id']} {run['run_id']} {status}", flush=True)
    return status


def execute(campaign_dir):
    manifest, config = load_contract()
    campaign_dir = Path(campaign_dir).resolve()
    digest_ledger = verify_campaign_digests(campaign_dir)
    ledger_path = campaign_dir / "run-ledger.jsonl"
    require(not ledger_path.exists(), "run ledger already exists; sample retries are forbidden")
    require(not (campaign_dir / "runs").exists(), "run output root already exists; sample retries are forbidden")
    preflight_fastqs(campaign_dir, manifest)
    verify_campaign_digests(campaign_dir)
    runs_by_cohort = defaultdict(list)
    for run in manifest["runs"]:
        runs_by_cohort[run["cohort"]].append(run)
    statuses = []
    for cohort in config["index_by_cohort"]:
        runs = runs_by_cohort[cohort]
        jobs = {run["planned_jobs"] for run in runs}
        require(len(jobs) == 1, f"cohort has inconsistent planned jobs: {cohort}")
        print(f"starting {cohort}: {len(runs)} runs, {next(iter(jobs))} jobs", flush=True)
        with concurrent.futures.ThreadPoolExecutor(max_workers=next(iter(jobs))) as executor:
            futures = [executor.submit(run_one, campaign_dir, run, config, digest_ledger, ledger_path) for run in runs]
            statuses.extend(future.result() for future in futures)
    verify_campaign_digests(campaign_dir)
    require(len(statuses) == 127, "campaign did not reach 127 terminal harness outcomes")
    print(json.dumps({status: statuses.count(status) for status in sorted(set(statuses))}, sort_keys=True))


def gzip_variants(campaign_dir, run):
    variants = Path(campaign_dir) / "reproducibility/input-variants" / run["run_id"]
    require(not variants.exists(), "gzip equivalence variants already exist")
    plain_dir = variants / "plain"
    multi_dir = variants / "multimember"
    plain_dir.mkdir(parents=True)
    multi_dir.mkdir(parents=True)
    source = Path(run["r1"]["path"])
    plain = plain_dir / source.name.removesuffix(".gz")
    multi = multi_dir / source.name
    midpoint = None
    with gzip.open(source, "rb") as decoded, plain.open("xb") as output:
        copied = 0
        while chunk := decoded.read(8 * 1024 * 1024):
            output.write(chunk)
            copied += len(chunk)
        midpoint = copied // 2
    with plain.open("rb") as decoded, multi.open("xb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as member:
            remaining = midpoint
            while remaining:
                chunk = decoded.read(min(8 * 1024 * 1024, remaining))
                if not chunk:
                    break
                member.write(chunk)
                remaining -= len(chunk)
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as member:
            shutil.copyfileobj(decoded, member, 8 * 1024 * 1024)
    return plain, multi


def execute_repro(campaign_dir):
    manifest, config = load_contract()
    campaign_dir = Path(campaign_dir).resolve()
    digest_ledger = verify_campaign_digests(campaign_dir)
    ledger_path = campaign_dir / "reproducibility/repro-ledger.jsonl"
    require(not ledger_path.exists(), "reproducibility ledger already exists; retries are forbidden")
    by_id = {run["run_id"]: run for run in manifest["runs"]}
    selected = [by_id[item["run_id"]] for item in config["reproducibility_runs"]]
    root = campaign_dir / "reproducibility/runs"
    require(not root.exists(), "reproducibility output root already exists")
    root.mkdir(parents=True)
    for run in selected:
        canonical = output_dir(campaign_dir, run) / "report.csv"
        require(canonical.is_file(), f"canonical campaign report missing: {run['run_id']}")
        for threads in config["reproducibility_threads"]:
            for repetition in range(1, config["reproducibility_repetitions"] + 1):
                repro_run = dict(run)
                repro_run["planned_threads"] = threads
                destination = root / run["run_id"] / f"threads-{threads}" / f"run-{repetition}"
                repro_run["_output"] = destination
                index_id = index_id_for(run, config)
                index = digest_ledger["indexes"][index_id]
                destination.parent.mkdir(parents=True, exist_ok=True)
                command = [BINARY_PATH, "run", "--r1", run["r1"]["path"]]
                if run.get("r2"):
                    command.extend(["--r2", run["r2"]["path"]])
                command.extend(["--index", index["path"], "--out", destination, "--threads", str(threads)])
                identity = {"run_id": run["run_id"], "threads": threads, "repetition": repetition}
                append_event(ledger_path, {**identity, "event": "STARTED", "at": now(), "command": [str(value) for value in command]})
                log = destination.parent / f"run-{repetition}"
                returncode, timing = run_timed(command, Path(str(log) + ".stdout"), Path(str(log) + ".stderr"), Path(str(log) + ".time"))
                require(returncode == 0, f"reproducibility run failed: {identity}")
                validate_success_directory(destination, repro_run, digest_ledger["profile"]["sha256"], index["index_digest"])
                identical = canonical.read_bytes() == (destination / "report.csv").read_bytes()
                append_event(ledger_path, {**identity, "event": "TERMINAL", "at": now(), "status": "SUCCESS", "report_csv_byte_identical": identical, "time": timing})
                require(identical, f"report.csv reproducibility failure: {identity}")
                print(f"repro {run['run_id']} t={threads} n={repetition} byte-identical", flush=True)
    variant_run = by_id[config["gzip_equivalence_run_id"]]
    plain, multi = gzip_variants(campaign_dir, variant_run)
    canonical = output_dir(campaign_dir, variant_run) / "report.csv"
    for variant, path in (("plain", plain), ("multimember", multi)):
        repro_run = dict(variant_run)
        repro_run["r1"] = dict(variant_run["r1"], path=str(path))
        repro_run["planned_threads"] = 8
        destination = root / variant_run["run_id"] / "encoding" / variant
        index = digest_ledger["indexes"][index_id_for(variant_run, config)]
        command = [BINARY_PATH, "run", "--r1", path, "--index", index["path"], "--out", destination, "--threads", "8"]
        identity = {"run_id": variant_run["run_id"], "encoding": variant}
        append_event(ledger_path, {**identity, "event": "STARTED", "at": now(), "command": [str(value) for value in command]})
        destination.parent.mkdir(parents=True, exist_ok=True)
        returncode, timing = run_timed(command, Path(str(destination) + ".stdout"), Path(str(destination) + ".stderr"), Path(str(destination) + ".time"))
        require(returncode == 0, f"encoding equivalence run failed: {variant}")
        validate_success_directory(destination, repro_run, digest_ledger["profile"]["sha256"], index["index_digest"])
        identical = canonical.read_bytes() == (destination / "report.csv").read_bytes()
        append_event(ledger_path, {**identity, "event": "TERMINAL", "at": now(), "status": "SUCCESS", "report_csv_byte_identical": identical, "time": timing})
        require(identical, f"report.csv encoding equivalence failure: {variant}")
    verify_campaign_digests(campaign_dir)


def main():
    parser = argparse.ArgumentParser(description="Frozen Phase 6 real-E2E campaign runner")
    parser.add_argument("action", choices=("prepare", "run", "repro"))
    parser.add_argument("--campaign-dir", required=True)
    args = parser.parse_args()
    {"prepare": prepare, "run": execute, "repro": execute_repro}[args.action](args.campaign_dir)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, json.JSONDecodeError) as error:
        raise SystemExit(f"Phase 6 runner failed: {error}")
