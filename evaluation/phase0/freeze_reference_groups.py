#!/usr/bin/env python3
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


PHASE0 = Path(__file__).resolve().parent
PROFILE = PHASE0 / "analysis-profile.json"
OUTPUT_DIR = PHASE0 / "reference-groups"
SOURCES = (
    (
        "internal-dna-panel",
        Path("/home/wubw/data/viroflash/dna_virus_68_20260827/refs/dna_virus_genome.fasta"),
    ),
    (
        "external-respiratory-panel",
        Path("/home/wubw/data/viroflash-benchmarks/refs/target_respiratory_panel.fa"),
    ),
    (
        "external-hpv-panel",
        Path("/home/wubw/data/viroflash-benchmarks/refs/target_hpv16_hpv18.fa"),
    ),
)
ALPHABET = frozenset(b"ACGTMRWSYKVHDBN")
COMPLEMENT = bytes.maketrans(b"ACGTMRWSYKVHDBN", b"TGCAKYWSRMBDHVN")
ASCII_WHITESPACE = b" \t\r\n\v\f"
LEDGER_FIELDS = (
    "group_ordinal",
    "target_group_id",
    "representative_id",
    "member_ordinal",
    "member_id",
    "representative_length",
    "target_fasta_sha256",
    "profile_digest",
)


def records(path):
    seen_ids = set()
    record_id = None
    sequence_parts = []
    with path.open("rb") as handle:
        for line_number, line in enumerate(handle, 1):
            if line.startswith(b">"):
                if record_id is not None:
                    if not sequence_parts:
                        raise ValueError(f"{path}:{line_number}: empty sequence for {record_id}")
                    yield record_id, b"".join(sequence_parts)
                tokens = line[1:].split()
                if not tokens:
                    raise ValueError(f"{path}:{line_number}: missing FASTA identifier")
                record_id = tokens[0].decode("utf-8")
                if record_id in seen_ids:
                    raise ValueError(f"{path}:{line_number}: duplicate FASTA identifier {record_id}")
                seen_ids.add(record_id)
                sequence_parts = []
                continue
            if record_id is None:
                if line.strip():
                    raise ValueError(f"{path}:{line_number}: sequence before first header")
                continue
            normalized = line.upper().translate(None, ASCII_WHITESPACE)
            invalid = sorted(set(normalized) - ALPHABET)
            if invalid:
                symbols = ", ".join(f"0x{byte:02x}" for byte in invalid)
                raise ValueError(f"{path}:{line_number}: invalid DNA symbol(s): {symbols}")
            sequence_parts.append(normalized)
    if record_id is None:
        raise ValueError(f"{path}: no FASTA records")
    if not sequence_parts:
        raise ValueError(f"{path}: empty sequence for {record_id}")
    yield record_id, b"".join(sequence_parts)


def canonical_sequence(sequence):
    reverse_complement = sequence.translate(COMPLEMENT)[::-1]
    return min(sequence, reverse_complement)


def profile_digest():
    return hashlib.sha256(PROFILE.read_bytes()).hexdigest()


def freeze_source(source_id, fasta_path, frozen_profile_digest):
    groups = defaultdict(list)
    digest_counts = Counter()
    fasta_digest = hashlib.sha256()
    with fasta_path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            fasta_digest.update(chunk)
    for member_id, sequence in records(fasta_path):
        sequence_digest = hashlib.sha256(canonical_sequence(sequence)).hexdigest()
        groups[sequence_digest].append((member_id, len(sequence)))
        digest_counts[sequence_digest] += 1

    duplicate_sequences = {}
    for member_id, sequence in records(fasta_path):
        canonical = canonical_sequence(sequence)
        sequence_digest = hashlib.sha256(canonical).hexdigest()
        if digest_counts[sequence_digest] > 1:
            previous = duplicate_sequences.setdefault(sequence_digest, canonical)
            if previous != canonical:
                raise ValueError(f"{fasta_path}: SHA-256 collision while grouping {member_id}")

    output_path = OUTPUT_DIR / f"{source_id}.tsv"
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LEDGER_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for group_ordinal, sequence_digest in enumerate(sorted(groups)):
            members = sorted(groups[sequence_digest])
            representative_id = members[0][0]
            representative_length = members[0][1]
            if any(length != representative_length for _, length in members):
                raise ValueError(f"{fasta_path}: unequal lengths in exact group {sequence_digest}")
            for member_ordinal, (member_id, _) in enumerate(members):
                writer.writerow({
                    "group_ordinal": group_ordinal,
                    "target_group_id": f"sha256:{sequence_digest}",
                    "representative_id": representative_id,
                    "member_ordinal": member_ordinal,
                    "member_id": member_id,
                    "representative_length": representative_length,
                    "target_fasta_sha256": fasta_digest.hexdigest(),
                    "profile_digest": frozen_profile_digest,
                })
    return {
        "source_id": source_id,
        "target_fasta": str(fasta_path),
        "target_fasta_sha256": fasta_digest.hexdigest(),
        "ledger": str(output_path.relative_to(PHASE0)),
        "ledger_sha256": hashlib.sha256(output_path.read_bytes()).hexdigest(),
        "record_count": sum(len(members) for members in groups.values()),
        "group_count": len(groups),
    }


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    frozen_profile_digest = profile_digest()
    frozen = [freeze_source(source_id, path, frozen_profile_digest) for source_id, path in SOURCES]
    manifest = {
        "contract_id": "viroflash.phase0.reference-groups",
        "profile_digest": frozen_profile_digest,
        "sources": frozen,
    }
    output = PHASE0 / "reference-groups.json"
    output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
