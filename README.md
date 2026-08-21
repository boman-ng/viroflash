# viroflash

`viroflash` 是一个 Rust 命令行工具，用于从 FASTQ/FASTQ.gz 测序数据和四类 FASTA 参考序列中生成病毒候选报告。

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

- 单端或双端 FASTQ/FASTQ.gz；
- 宿主 FASTA；
- 目标病毒 FASTA；
- 诱饵 FASTA；
- 污染参考 FASTA。

## 运行

```bash
viroflash run \
  --r1 reads_R1.fastq.gz \
  --r2 reads_R2.fastq.gz \
  --host-fa host.fa \
  --target-fa target.fa \
  --decoy-fa decoy.fa \
  --contam-fa contaminant.fa \
  --threads 8 \
  --out result
```

省略 `--r2` 即使用单端模式。`--k` 可设置预筛 k-mer 长度，取值范围为 1 到 31。

生成确定性合成输入：

```bash
python3 scripts/make_virtual_data.py synthetic.work --scenario all
```

生成分层诱饵：

```bash
viroflash decoy \
  --target-fa target.fa \
  --out decoy.fa \
  --ani 82,85,88 \
  --per-layer 4 \
  --seed 0
```

## 输出

`run` 根据 `--out` 前缀写出：

- `<out>.json`：运行摘要和候选详情；
- `<out>.tsv`：候选表格。

`decoy` 写出诱饵 FASTA；提供 `--report` 时同时写出诱饵明细 TSV。
