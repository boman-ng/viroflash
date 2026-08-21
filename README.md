# viroflash

`viroflash` 是一个 Rust 命令行工具，用于从 FASTQ/FASTQ.gz 测序数据中生成病毒候选报告。管线：k-mer 预筛 → minimap2 竞争比对 → 候选聚合 → 诱饵零分布统计判定 → JSON/TSV 报告。

## 构建与测试

```bash
cargo build --release
cargo test
```

查看命令行帮助：

```bash
cargo run -- --help
```

## 输入

`viroflash run` 接受：

- 单端或双端 FASTQ/FASTQ.gz（省略 `--r2` 即单端模式）；
- 一个索引目录（`--index`），或用于自动构建索引的参考 FASTA（见下）。

## 索引

`viroflash index` 构建可复用索引目录：

```bash
viroflash index \
  --host-fa host.fa \
  --target-fa target.fa \
  --contam-fa contaminant.fa \
  --decoy-fa decoy.fa \
  --out index_dir \
  --k 21 \
  --threads 8
```

产物为一个目录：

- `ref.mmi`：minimap2 sr 索引；
- `bloom.bin`：目标+诱饵 canonical k-mer Bloom（预筛门）；
- `manifest.json`：版本、k、角色→contig 映射、来源 FASTA 的 BLAKE3 校验和、诱饵来源参数。

`--contam-fa` 可选；`--decoy-fa` 可选，省略时按 `--decoy-ani`（默认 `82,85,88`）、`--decoy-per-layer`（默认 `4`）、`--decoy-seed`（默认 `0`）从目标基因组确定性生成 SNP-only 诱饵（产物 `decoys.fa`/`decoys.tsv` 留在索引目录内）。

## 运行

复用已构建索引（推荐；`--index` 与 FASTA 参数互斥，`--k` 须与索引一致）：

```bash
viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --index index_dir \
  --threads 8 \
  --out result
```

未提供 `--index` 时自动构建索引（需 `--host-fa`/`--target-fa`，`--contam-fa`/`--decoy-fa` 可选，诱饵省略时按默认参数自动生成），构建产物落 `<out>.work/index/`：

```bash
viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --host-fa host.fa \
  --target-fa target.fa \
  --threads 8 \
  --out result
```

两条路径共用同一构建入口，结果一致。`--k` 可设置预筛 k-mer 长度，取值范围为 1 到 31。

生成确定性合成输入：

```bash
python3 scripts/make_virtual_data.py synthetic.work --scenario all
```

## 输出

`run` 根据 `--out` 前缀写出：

- `<out>.json`：运行摘要、所用索引溯源（顶层 `index` 块）和候选详情；
- `<out>.tsv`：候选表格。
