#!/usr/bin/env python3
"""虚拟数据集生成器（合成 FASTQ.gz + 4 类参考 FASTA）。

场景（对应目标验证要求）：
  positive    阳性对照：reads 全部来自目标病毒（应检出）
  chimeric    嵌合：病毒整合入宿主，reads 跨越整合断点（应检出 + split 位点）
  contaminant 污染：reads 全部来自污染基因组（不应检出目标）
  decoy       诱饵噪声：reads 来自诱饵序列 + 宿主（不应检出目标）
  negative    阴性：reads 全部来自宿主（不应检出目标）

用法:
  scripts/make_virtual_data.py <输出目录> --scenario all [--seed 1] [--pairs 2000]
仅用 Python 标准库，确定性可复现。
"""

import argparse
import gzip
import os
import random
import sys

READ_LEN = 150
HOST_LEN = 100_000
TARGET_LEN = 5_000
DECOY_LEN = 5_000
CONTAM_LEN = 10_000
DECOY_COUNT = 20
INSERTION_SITE = 50_000  # 嵌合场景整合位点（宿主坐标）

COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def randseq(rng: random.Random, n: int, gc: float = 0.45) -> str:
    at = (1.0 - gc) / 2.0
    weights = (at, at, gc / 2.0, gc / 2.0)
    return "".join(rng.choices("ATCG", weights=weights, k=n))


def write_fasta(path: str, records: list[tuple[str, str]]) -> None:
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i : i + 60] + "\n")


def write_fastq_gz(path: str, records: list[tuple[str, str]]) -> None:
    with gzip.open(path, "wt") as f:
        for rid, seq in records:
            f.write(f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n")


def make_references(rng: random.Random, out_dir: str) -> dict[str, str]:
    host = randseq(rng, HOST_LEN, 0.42)
    target = randseq(rng, TARGET_LEN, 0.48)
    decoys = [randseq(rng, DECOY_LEN, 0.48) for _ in range(DECOY_COUNT)]
    contam = randseq(rng, CONTAM_LEN, 0.30)

    write_fasta(os.path.join(out_dir, "host.fa"), [("chrH", host)])
    write_fasta(os.path.join(out_dir, "target.fa"), [("TESTVIR", target)])
    write_fasta(
        os.path.join(out_dir, "decoy.fa"),
        [(f"dec{i}", d) for i, d in enumerate(decoys)],
    )
    write_fasta(os.path.join(out_dir, "contam.fa"), [("mycoplasma", contam)])
    return {"host": host, "target": target, "decoys": decoys, "contam": contam}


def pair(rng: random.Random, r1: str, r2: str, i: int) -> tuple[tuple[str, str], tuple[str, str]]:
    """构造一对读。r1/r2 为正链序列，r2 输出为反向互补（Illumina 惯例）。"""
    jitter = rng.choice([0, 0, 1, -1, 2])  # 少量坐标抖动
    return ((f"r{i}/1", r1), (f"r{i}/2", revcomp(r2)))


def scenario_positive(ref: dict, rng: random.Random, pairs: int):
    r1s, r2s = [], []
    t = ref["target"]
    for i in range(pairs):
        s = rng.randrange(0, TARGET_LEN - 2 * READ_LEN)
        r1s.append((f"p{i}/1", t[s : s + READ_LEN]))
        r2s.append((f"p{i}/2", revcomp(t[s + 200 : s + 200 + READ_LEN])))
    return r1s, r2s


def scenario_chimeric(ref: dict, rng: random.Random, pairs: int):
    """50% 跨断点读（宿主半 + 病毒半），50% 纯宿主读；断点两侧各半支持方向。"""
    host, t = ref["host"], ref["target"]
    r1s, r2s = [], []
    for i in range(pairs):
        if i % 2 == 0:
            half = READ_LEN // 2
            if i % 4 == 0:
                # 左跨断点：宿主[site-half, site) + 病毒[0, half)
                r1 = host[INSERTION_SITE - half : INSERTION_SITE] + t[0:half]
            else:
                # 右跨断点：病毒[-half:] + 宿主[site, site+half)
                r1 = t[-half:] + host[INSERTION_SITE : INSERTION_SITE + half]
            mate = rng.randrange(0, HOST_LEN - READ_LEN)
            r2 = revcomp(host[mate : mate + READ_LEN])
        else:
            s = rng.randrange(0, HOST_LEN - READ_LEN)
            r1 = host[s : s + READ_LEN]
            r2 = revcomp(host[(s + 200) % (HOST_LEN - READ_LEN) : (s + 200) % (HOST_LEN - READ_LEN) + READ_LEN])
        r1s.append((f"c{i}/1", r1))
        r2s.append((f"c{i}/2", r2))
    return r1s, r2s


def scenario_contaminant(ref: dict, rng: random.Random, pairs: int):
    c = ref["contam"]
    r1s, r2s = [], []
    for i in range(pairs):
        s = rng.randrange(0, CONTAM_LEN - 2 * READ_LEN)
        r1s.append((f"m{i}/1", c[s : s + READ_LEN]))
        r2s.append((f"m{i}/2", revcomp(c[s + 200 : s + 200 + READ_LEN])))
    return r1s, r2s


def scenario_decoy(ref: dict, rng: random.Random, pairs: int):
    host = ref["host"]
    r1s, r2s = [], []
    for i in range(pairs):
        if i % 3 == 0:
            d = ref["decoys"][i % DECOY_COUNT]
            s = rng.randrange(0, DECOY_LEN - 2 * READ_LEN)
            r1s.append((f"d{i}/1", d[s : s + READ_LEN]))
            r2s.append((f"d{i}/2", revcomp(d[s + 200 : s + 200 + READ_LEN])))
        else:
            s = rng.randrange(0, HOST_LEN - 2 * READ_LEN)
            r1s.append((f"d{i}/1", host[s : s + READ_LEN]))
            r2s.append((f"d{i}/2", revcomp(host[s + 200 : s + 200 + READ_LEN])))
    return r1s, r2s


def scenario_negative(ref: dict, rng: random.Random, pairs: int):
    host = ref["host"]
    r1s, r2s = [], []
    for i in range(pairs):
        s = rng.randrange(0, HOST_LEN - 2 * READ_LEN)
        r1s.append((f"n{i}/1", host[s : s + READ_LEN]))
        r2s.append((f"n{i}/2", revcomp(host[s + 200 : s + 200 + READ_LEN])))
    return r1s, r2s


SCENARIOS = {
    "positive": scenario_positive,
    "chimeric": scenario_chimeric,
    "contaminant": scenario_contaminant,
    "decoy": scenario_decoy,
    "negative": scenario_negative,
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out_dir")
    ap.add_argument("--scenario", default="all",
                    help="all / positive / chimeric / contaminant / decoy / negative")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--pairs", type=int, default=2000)
    args = ap.parse_args()

    names = list(SCENARIOS) if args.scenario == "all" else [args.scenario]
    for name in names:
        if name not in SCENARIOS:
            print(f"未知场景: {name}", file=sys.stderr)
            return 1
    os.makedirs(args.out_dir, exist_ok=True)

    for name in names:
        scenario_dir = os.path.join(args.out_dir, name)
        os.makedirs(scenario_dir, exist_ok=True)
        rng = random.Random(f"{args.seed}-{name}")
        ref = make_references(rng, scenario_dir)
        r1s, r2s = SCENARIOS[name](ref, rng, args.pairs)
        write_fastq_gz(os.path.join(scenario_dir, "reads_R1.fq.gz"), r1s)
        write_fastq_gz(os.path.join(scenario_dir, "reads_R2.fq.gz"), r2s)
        print(f"[{name}] 生成 {len(r1s)} pairs → {scenario_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
