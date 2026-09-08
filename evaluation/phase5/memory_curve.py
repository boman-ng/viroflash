#!/usr/bin/env python3
import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def sequence(seed, length):
    state = seed
    bases = b"ACGT"
    result = bytearray()
    for _ in range(length):
        state ^= state << 13
        state ^= state >> 7
        state ^= state << 17
        state &= (1 << 64) - 1
        result.append(bases[state & 3])
    return bytes(result)


def write_fasta(path, name, bases):
    with path.open("wb") as handle:
        handle.write(b">" + name.encode() + b"\n")
        for start in range(0, len(bases), 60):
            handle.write(bases[start : start + 60] + b"\n")


def write_fastq(path, read, fragments):
    quality = b"I" * len(read)
    with path.open("wb") as handle:
        for ordinal in range(fragments):
            handle.write(
                b"@fragment-"
                + str(ordinal).encode()
                + b"\n"
                + read
                + b"\n+\n"
                + quality
                + b"\n"
            )


def run_checked(arguments):
    completed = subprocess.run(arguments, capture_output=True, text=True)
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(map(str, arguments))}\n"
            f"{completed.stderr}"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=Path("target/release/viroflash"))
    parser.add_argument("--sizes", type=int, nargs="+", default=[256, 4096, 65536])
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    sizes = sorted(set(args.sizes))
    if len(sizes) < 3 or sizes[0] <= 0 or args.threads <= 0:
        raise SystemExit("provide at least three positive sizes and positive threads")
    binary = args.binary.resolve()
    if not binary.is_file():
        raise SystemExit(f"binary does not exist: {binary}")

    measurements = []
    with tempfile.TemporaryDirectory(prefix="viroflash-phase5-memory-") as temporary:
        root = Path(temporary)
        host = sequence(17, 600)
        target = sequence(91, 600)
        write_fasta(root / "host.fa", "host", host)
        write_fasta(root / "target.fa", "target", target)
        run_checked(
            [
                binary,
                "index",
                "--host-fa",
                root / "host.fa",
                "--target-fa",
                root / "target.fa",
                "--out",
                root / "index",
                "--threads",
                str(args.threads),
            ]
        )
        profile_digest = None
        index_digest = None
        for fragments in sizes:
            fastq = root / f"sample-{fragments}.fastq"
            output = root / f"out-{fragments}"
            write_fastq(fastq, target[100:220], fragments)
            run_checked(
                [
                    binary,
                    "run",
                    "--r1",
                    fastq,
                    "--index",
                    root / "index",
                    "--out",
                    output,
                    "--threads",
                    str(args.threads),
                ]
            )
            if {path.name for path in output.iterdir()} != {
                "report.csv",
                "report.html",
                "perf.json",
            }:
                raise RuntimeError("successful run did not produce exactly three outputs")
            perf = json.loads((output / "perf.json").read_text(encoding="utf-8"))
            run_row = (output / "report.csv").read_text(encoding="utf-8").splitlines()
            header = run_row[0].split(",")
            values = run_row[1].split(",")
            report = dict(zip(header, values))
            profile_digest = profile_digest or report["profile_digest"]
            index_digest = index_digest or report["index_digest"]
            assert report["profile_digest"] == profile_digest
            assert report["index_digest"] == index_digest
            measurements.append(
                {
                    "input_fragments": fragments,
                    "fastq_bytes": fastq.stat().st_size,
                    "selected_fragments": int(report["selected_fragments"]),
                    "selected_sequence_bytes": int(report["selected_fragments"]) * 120,
                    "peak_rss_bytes": perf["peak_rss_bytes"],
                    "wall_time_ms": perf["wall_time_ms"],
                    "telemetry_status": perf["telemetry_status"],
                }
            )

    adjacent_slopes = []
    for left, right in zip(measurements, measurements[1:]):
        selected_delta = right["selected_sequence_bytes"] - left["selected_sequence_bytes"]
        rss_delta = right["peak_rss_bytes"] - left["peak_rss_bytes"]
        adjacent_slopes.append(
            {
                "from_input_fragments": left["input_fragments"],
                "to_input_fragments": right["input_fragments"],
                "rss_delta_bytes": rss_delta,
                "selected_sequence_delta_bytes": selected_delta,
                "rss_delta_per_selected_sequence_byte": rss_delta / selected_delta,
            }
        )
    result = {
        "measurements": measurements,
        "adjacent_slopes": adjacent_slopes,
        "structural_evidence": {
            "sequence_reservoir": "absent",
            "fastq_reader_buffer_bytes_per_input_end": 1048576,
            "alignment_task_queue_capacity_fragments": f"{args.threads} threads x 1",
            "alignment_result_queue_capacity_fragments": f"{args.threads} threads x 1",
            "retained_evidence_owner": "fixed reference groups plus merged reference intervals",
        },
        "evidence_limits": [
            "The curve is an observation, not a proof of a universal RSS bound.",
            "Linux VmHWM is a high-water estimate and may be imprecise.",
            "Minimap2 internal allocations are opaque to this structural audit.",
            "Queue item count is bounded, but one arbitrarily long FASTQ record is not byte-capped.",
            "The structural bound is independent of selected-fragment count only for a fixed index, thread count, and maximum fragment length.",
        ],
        "acceptance_role": "raw_measurement_no_subjective_tolerance",
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
