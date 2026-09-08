#!/usr/bin/env python3

import pathlib
import sys


def sequence(seed: int, length: int) -> str:
    state = seed
    bases = "ACGT"
    output = []
    for _ in range(length):
        state ^= state << 13
        state ^= state >> 7
        state ^= state << 17
        state &= (1 << 64) - 1
        output.append(bases[state & 3])
    return "".join(output)


def write_fasta(path: pathlib.Path, identifier: str, bases: str) -> None:
    with path.open("w", encoding="ascii") as output:
        output.write(f">{identifier}\n")
        for start in range(0, len(bases), 60):
            output.write(f"{bases[start:start + 60]}\n")


def main() -> None:
    output_dir = pathlib.Path(sys.argv[1])
    output_dir.mkdir(parents=True, exist_ok=True)
    host = sequence(17, 600)
    target = sequence(91, 600)
    write_fasta(output_dir / "host.fa", "host", host)
    write_fasta(output_dir / "target.fa", "target", target)
    with (output_dir / "sample.fastq").open("w", encoding="ascii") as output:
        for ordinal in range(20):
            read = target[ordinal % 30:ordinal % 30 + 120]
            output.write(f"@fragment-{ordinal}\n{read}\n+\n{'I' * len(read)}\n")


if __name__ == "__main__":
    main()
