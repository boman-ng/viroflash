# viroflash

`viroflash` 是一个 Rust 命令行工具，用于从 FASTQ/FASTQ.gz 测序数据中生成病毒候选报告。管线：确定性 bottom-k → k-mer 预筛 → minimap2 全链竞争审计 → discovery/validation 直接证据等价类 → synthetic-decoy 背景统计 → JSON/TSV 报告；`index` 与 `run` 同时生成独立的性能报告。运行时完全离线，不查询 NCBI、taxonomy 或特定 Panel 先验。

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

- `ref.mmi`：minimap2 sr 单分片索引，使宿主/目标/诱饵/污染命中在同一全局竞争空间计算 MAPQ；
- `bloom.bin`：全部原始目标+实际诱饵的 canonical k-mer Bloom（预筛门；按实测填充率折叠到满足假阳性率上界的最小安全容量）；
- `targets.fa`：由序列相似度算法离线构造的检测组代表，仅这些代表进入 MMI；
- `manifest.json`：版本、k、角色→contig 映射、代表→原始目标成员、来源 FASTA 的 BLAKE3 校验和、诱饵来源参数。

检测组不使用名称、taxonomy 或数据库查询：先折叠完全相同序列，再用固定 k=21 的 FracMinHash/PPJoin 与 Mash 候选筛选，最后由 minimap2 全长确认 `identity≥97%` 且双向 `coverage≥95%` 的代表星形成组。该静态组只压缩索引并保留原始成员映射，不把成员 counts 相加；运行时的复合 OR 假设仍只由 read 的直接最高分兼容集合构造。

索引目录旁同时写出 `<out>.perf.json` 与 `<out>.perf.tsv`。

`--contam-fa` 可选；`--decoy-fa` 可选，省略时按 `--decoy-ani`（默认 `85`）、`--decoy-per-layer`（默认 `1`）、`--decoy-seed`（默认 `0`）从每个检测组代表确定性生成一条 SNP-only synthetic decoy（产物 `decoys.fa`/`decoys.tsv` 留在索引目录内）。额外 ANI 层和重复可显式用于独立压力测试，但默认统计不让索引按层数×重复数膨胀。target-derived synthetic decoy 是近邻竞争和模型 null 的压力参考，不等同于提取空白、同批 NTC 或实验室污染背景；其可交换性尚未独立验证，报告不宣称经典 FDR 保证。
加载时会拒绝多分片 MMI；这类索引的跨分片 MAPQ 不满足竞争比对契约，需用当前版本重建。

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
`--threads` 是数据管线的计算线程预算（主线程 caller-runs + 比对 worker）；带 IG 成员索引的 gzip 输入可另建少量有界、低占用解压 helper，普通 gzip 不创建。多线程运行另有 1 个低频、通常休眠的 telemetry 采样线程；这些辅助线程都不缩减计算预算。

每个样本先在全部 fragment 上做固定容量、域分离 BLAKE3 bottom-k（输入不超过容量时即全量），再只对入样 pair 做 Bloom 与比对。fragment 被确定性分到 discovery 或 validation，同一 pair 两端不会跨折。一个 read-end 对多个最高分目标只形成一个直接兼容集合，不按成员数重复计票；报告中的复合 OR 假设表示“当前参考模型下至少一个成员存在”，不等于已定位到某一 accession。

discovery 只固定待检验的目标 OR 假设。validation 中每个直接 read-end 在一个假设内只计一次；target exposure 是该假设实际包含的内部索引 contig 长度总和，而不是代表序列长度或展开后的原始 accession 数。单内部成员假设与同 size/GC 层的固定 synthetic-decoy 比较；跨多个内部成员的假设使用全局固定 decoy 集，双方 read-end 均去重、reference bases 均各计一次。decoy validation 计数不要求该 decoy 先在 discovery 出现，避免把未入选误作零背景。主检验是以参考长度为 exposure 的单侧两 Poisson 率精确条件检验；背景零命中时仍给出有限概率，而不是把背景率视作已知的精确零。全部 discovery 假设（包括 validation 零计数）进入固定族 Benjamini–Hochberg 调整。兼容字段 `q_value` 现表示该模型 adjusted-p；由于 synthetic-decoy null 尚未校准，它不能解释为已验证的 FDR 或临床置信度。

候选报告门依次使用模型 adjusted-p `<0.2`、代表序列 observed breadth `≥10%` 和 10 个固定位置分区中至少 3 个有直接 evidence。三个数值均在 JSON 中标明适用性或未校准状态；end-RPM 只作描述。通过门输出 `PASS`，覆盖/分布未过门输出 `BELOW_THRESHOLD`，但二者都是候选级研究状态。宿主–病毒 split/site 只写入独立 `integration_evidence`，不再扣除检测 reads，也不决定通用病毒检出。

生成确定性合成输入：

```bash
python3 scripts/make_virtual_data.py synthetic.work --scenario all
```

## 输出

`run` 根据 `--out` 前缀写出：

- `<out>.json`：运行摘要、抽样/折分/审计元数据、所用索引溯源和复合候选详情；
- `<out>.tsv`：候选表格；
- `<out>.perf.json`：进程/整机 CPU、RSS、虚拟内存、累计 CPU 时间、I/O、采样完整性及阶段级墙钟/CPU/RSS；
- `<out>.perf.tsv`：同一次运行的固定列性能摘要，便于批量汇总。

候选结果 JSON 的当前契约为 `viroflash.result.v1`。v1 从契约层明确引入模型
adjusted-p、`UNVALIDATED`、候选门与 integration 解耦等语义；它不是旧 v0 字段的
静默重解释。TSV 继续保持 22 列；非空数据行在 `notes.result_schema` 标明同一版本，
空 TSV 的版本、QC 和固定检验族元数据以配对 JSON 为准。

性能报告不写入输入或索引路径。单线程运行采用阶段边界采样；多线程运行默认每秒定向采样当前进程。RSS 包含共享页且不等同于 PSS，报告也不声称提供 NUMA、cache miss 或内存带宽指标。

### 结果解释与限制

- 输出作用域是候选，不计算样本级诊断结论；`candidates: []` 只表示“没有候选被报告”，不等价于 `NOT_DETECTED`。对应 TSV 只保留表头，不伪造 `NOT_SIGNIFICANT` 候选行。
- `decision=PASS` 表示通过当前探索性候选门，不代表临床阳性。`confidence=UNVALIDATED` 明确表示当前 adjusted-p 未完成外部 null/FDR 校准。
- `p_value` 是精确条件率检验结果；`q_value` 是固定 discovery 检验族的 BH adjusted-p。极小值若发生 f64 下溢，普通字段可为 `0`，但 `statistical_evidence.log10_*` 与状态字段保留其数值含义。
- 顶层 `candidate_testing` 披露完整 discovery 检验族、已报告候选和省略的零证据假设数；即使 `candidates: []` 也可审计多重校正边界。
- v0 兼容字段 `poisson_p` 仅保留旧 plug-in Poisson 诊断值，不参与 adjusted-p 或判定；`nb_p` 当前禁用并输出 `null`。
- OR `members` 是未解析的备选集合；不能把代表 accession 或任一成员单独归因。TSV 的 `virus_id` 是索引内部 hypothesis ID、`resolution_level=reference_group`，`notes` 明示 OR 语义和代表名称；空 `taxid` 表示本工具未查询 taxonomy。
- JSON 的 `quality_control.status=NOT_EVALUATED` 表示只报告 map error、audit overflow 和 unassigned 等观测量，尚无经验证的 QC 通过/失败阈值。
- TSV 保持既有 22 列兼容：`split_reads` 当前实际为 split event 数，`aligned_bases` 为代表参考上的 covered bases，`evidence_strength` 当前为 `UNVALIDATED/NOT_SIGNIFICANT`；精确单位和限制以 JSON 为准。
- 固定容量抽样、参考库范围、synthetic decoy、10% breadth、3/10 windows 和比对门均需在目标样本基质、病毒类别、LoD/LoB、近缘干扰及独立阴性对照上验证后，才能支持临床性能或通用敏感度/特异度主张。
