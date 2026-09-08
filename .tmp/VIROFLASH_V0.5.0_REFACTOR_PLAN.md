# Viroflash v0.5.0 优化重构方案

> 文档性质：实现与验收方案，不是代码变更。
> 版本约束：Cargo.toml 当前版本 0.3.0 保持不变；“v0.5.0”仅是本方案代号，不授权修改软件、索引或 schema 版本号。
> 输出约束：每个样本运行目录最终仅有 report.html、report.csv、perf.json。
> 数据约束：终验必须从原始 FASTQ 开始重跑 68 个内部样本与 59 个外部公开样本，不得用历史 JSON/TSV 代替。

## 1. 执行摘要

v0.5.0 的核心不是继续增加阈值或统计模型，而是把当前系统收敛为一条可解释、可复算、不可通过运行时调参改变结论的路径：

~~~text
HOST + TARGET genomes
  → 冻结 analysis profile 与 target groups
  → 第 1 遍完整校验并精确计数 FASTQ fragments
  → 由 N、最低相关比例和错误风险计算纳入概率
  → 第 2 遍确定性哈希等概率抽样
  → target k-mer Bloom 工作量门
  → HOST + TARGET 竞争比对
  → fragment 级 target-group 归属
  → profile-attributed fragment fraction + simultaneous interval
  → report.csv / report.html / perf.json
~~~

核心裁决：

1. 删除固定 2^19 fragment 容量；样本量由统计目标和真实输入规模计算，不设经验性最小/最大 reads。
2. 抽样单位统一为 fragment；PE 两端同入同出且一个 fragment 最多计一次。
3. 主估计量改为 **profile-attributed fragment fraction**，即“在冻结 Viroflash profile 下归属于某 target group 的输入文库 fragment 比例”。
4. 该比例不是原始样本中的绝对病毒载量、copies/mL、临床阳性概率或完整生物学丰度。
5. 删除 synthetic decoy、discovery/validation split、数据依赖 hypothesis、Poisson p 值、BH/q、PASS 和 coverage/windows cutoff 联合判定。
6. 不实现 EM、latent abundance、confidence sequence、插件化 profile、旧索引转换器或旧报告双写。
7. Bloom 只减少进入 minimap2 的工作量，不再被解释为“人类 reads 的统计排除器”或阳性证据。
8. 资源不足只产生明确的非结论状态；不得降低统计目标、改变抽样概率或产生成功形状的结果。

## 2. 背景与当前系统

### 2.1 当前真实数据流

源码显示当前流程为：

~~~text
全部 FASTQ
  → 固定 bottom-k 524,288 fragments
  → target + decoy Bloom（不含 host）
  → host + target + decoy + contaminant 联合竞争比对
  → discovery / validation 哈希拆分
  → discovery 生成 hypothesis
  → validation 条件双 Poisson 检验
  → BH adjusted p
  → adjusted p < 0.20、coverage ≥ 10%、windows ≥ 3
  → JSON / TSV / HTML / CSV / perf JSON
~~~

代码锚点：

- 固定样本容量：src/sampling.rs:8。
- reservoir 无条件使用固定容量：src/lib.rs:547。
- 三个生产判定阈值：src/lib.rs:62、src/lib.rs:63、src/lib.rs:66。
- 候选判定组合：src/lib.rs:1222。
- 当前比例实际为 validation read-end RPM：src/lib.rs:1198、src/report.rs:398。
- 当前报告仍输出五类文件：README.md:58、src/main.rs:278。

### 2.2 当前统计判定矩阵

| adjusted p | coverage | distributed windows | 当前 decision | 实际含义 |
|---:|---:|---:|---|---|
| ≥ 0.20 | 任意 | 任意 | NOT_SIGNIFICANT | 未通过未校准 synthetic-decoy 模型的 adjusted p 门 |
| < 0.20 | < 10% | 任意 | BELOW_THRESHOLD | 通过 p 门但未通过覆盖门 |
| < 0.20 | ≥ 10% | < 3 | BELOW_THRESHOLD | 通过 p/coverage 门但未通过分布门 |
| < 0.20 | ≥ 10% | ≥ 3 | PASS | 通过三个候选报告门 |
| 无 discovery hypothesis | 不适用 | 不适用 | 无候选 | 不能解释为病毒阴性 |

主要问题：

1. 0.20、10%、3 是三个独立自由度，缺少当前数据条件下的联合校准。
2. synthetic decoy 不代表宿主同源、未建模微生物、载体、试剂污染或 index hopping；其交换性未验证。
3. BH adjusted p 因 null proxy 未验证，不能自动获得经典 FDR 语义。
4. discovery/validation 拆分让固定目标族损失约一半证据，并增加动态 testing family。
5. read-end RPM 的分子与分母单位不一致于用户期望的 fragment 比例，且没有区间。
6. PASS 混合了统计、分布和产品阈值，不等于病毒存在、临床阳性或含量可信。

### 2.3 已确认的资源风险

- 当前 reservoir 保存最多 524,288 对完整序列，内存与 read 长度成正比。
- src/fastq.rs 的多成员 gzip 有序合并使用 pending BTreeMap；队列按对象数限制，但乱序缓存没有严格字节上界。
- 当前内部 68 样本基线单样本 RSS 中位数 25.26 GiB、最大 26.14 GiB；大索引构建峰值 45.53 GiB。
- 根据瞬时 RSS、线程调度或机器内存改变样本量，会导致相同数据在不同机器上产生不同科学结果，禁止采用。

## 3. 目标、边界与非目标

### 3.1 目标

1. 用统计目标反推抽样概率，而不是固定 reads 数或按文件大小线性缩放。
2. 在大体积差异样本上保持内存有界、结果确定、跨线程一致。
3. 输出 fragment 级 target-group 比例、区间、证据充分性和归属不确定性。
4. 将科学参数与资源参数严格分离，消除生产运行中的结果调参入口。
5. 将报告压缩为三个明确产物，并保证 HTML/CSV 同源。
6. 用 68 个内部样本和 59 个外部样本完成真实端到端评估。

### 3.2 可识别边界

仅有 HOST、TARGET 和 FASTQ 时，可以估计：

> 输入 FASTQ 文库中，被冻结 Viroflash profile 归属于某 target group 的 fragment 比例。

不能从这些输入单独识别：

- 原始样本中的病毒 copies/mL、绝对载量或总病毒颗粒数；
- 临床阳性/阴性、PPV、NPV、LoB、LoD、LoQ；
- 未提供污染物、近缘物种或未知序列导致的系统误差；
- accession 或 group 内成员的真实混合比例；
- PCR duplicates 对独立分子数的影响；
- Bloom、比对器和参考缺失造成的完整生物学 sensitivity。

报告必须使用上述限定语，不得把 profile-attributed fragment fraction 简写成 abundance 或 viral load。

### 3.3 明确非目标

- EM、VBEM 或其他 mixture deconvolution；
- accession 级硬分配；
- absolute quantitation；
- synthetic null、p/q 或 FDR 声明；
- evidence-driven 在线早停和 confidence sequence；
- 通用 profile registry、插件系统或多种统计后端；
- 自动学习阈值、训练模型或利用 68/59 样本调参；
- 旧 index、旧 JSON/TSV、旧字段名或旧 CLI 的兼容层；
- 自动迁移历史结果；
- 根据可用内存动态降低科学要求。

## 4. 科学合同

### 4.1 观测单位

- 总体单位是 fragment。
- SE：一条 read 是一个 fragment。
- PE：R1/R2 是同一 fragment，两端同入同出。
- 一个 fragment 对一个 target group 最多贡献一次。
- split、discordant、strand 和覆盖位置是同一 fragment 的诊断属性，不新增计数单位。

### 4.2 固定目标族

target groups 必须只由 HOST/TARGET references 和冻结的 grouping 规则在索引构建期产生。

- 运行时 FASTQ 不得创建、拆分、合并或重排确认性 group。
- 完全不可区分的 target members 在索引期合并为一个 ReferenceGroup。
- group 内只允许 OR 语义；不能声称具体 member 已解析。
- 跨 group 的等价证据标记为 UNRESOLVED_ACROSS_GROUPS，不强行分配。
- 删除运行期 discovery hypothesis，因此不再需要 discovery/validation fold。

### 4.3 主估计量

设输入有限总体含 N 个 fragments。对冻结 target group g，定义：

~~~text
Y_ig = 1，当 fragment i 通过冻结 k-mer gate，
           且竞争比对将其解析到 group g；
       0，否则。

P_g = (1 / N) × Σ_i Y_ig
~~~

P_g 是 profile-attributed fragment fraction。

该定义有意把 Bloom gate、比对器和归属规则纳入 profile；因此它是可复现的软件测量值，而不是参考无关的生物学真值。

### 4.4 最小不可避免的科学参数

v0.5 只保留一个编译期、公开、可审计的 AnalysisProfile，不提供运行时 profile 选择器。profile 至少冻结：

| 参数 | 含义 | 为什么无法由同一 FASTQ 自动推导 |
|---|---|---|
| minimum_relevant_fraction，δ | 产品希望有能力观察的最低文库比例 | 它表达科学/产品相关性，不是数据事实 |
| familywise_miss_probability，β | 对冻结 target family 的允许漏抽风险 | 它表达错误代价 |
| familywise_interval_error，α | 所有已报告 target intervals 的联合覆盖错误风险 | 它表达置信要求 |

不再保留 minimum reads、q cutoff、coverage cutoff、window cutoff 等互相独立的生产判定参数。

任何 profile 参数变更都必须：

1. 修改唯一源码所有者；
2. 生成新的 profile digest；
3. 重建索引；
4. 对 127 个样本完整重跑；
5. 不得把失败后的调参重跑称为同一次验证。

### 4.5 抗调参合同

- 生产 CLI 不暴露 k、sample size、阈值、seed、decoy 或 score margin。
- profile digest 在读取 FASTQ 结果前确定，并写入 index、CSV、HTML。
- 抽样 seed 由 domain-separated hash(profile digest, input digest) 派生，不允许用户指定。
- target family 顺序固定，所有 multiplicity denominator 来自索引 ledger。
- 线程数、chunk 大小和队列容量属于资源参数，不得改变 report.csv。
- 任何科学参数变化产生不同 analysis identity，不允许覆盖原结果。

## 5. 动态抽样设计

### 5.1 为什么不是“按文件体积成比例抽样”

对于最低比例 δ，当 N 远大于所需样本量时，达到既定检出概率所需的抽样量主要由 δ 和 β 决定，而不是与 N 线性增长。N 只通过有限总体修正、全量纳入和实际纳入概率参与计算。

因此“动态”应解释为：

> 使用每个样本的精确 fragment 总数计算统计设计，而不是使用固定 524,288、压缩字节数、主观上下限或机器资源。

### 5.2 两遍流式路径

**Pass 1：Input census**

- 完整读取并校验 FASTQ；
- 校验 PE 数量与 normalized pair ID；
- 精确计算 N、read length 摘要和 input digest；
- 不保存 read 序列、不运行 Bloom、不比对；
- gzip 并行缓冲必须受字节预算约束。

**Design**

- m 为冻结 ReferenceGroup 数。
- M_min = ceil(δN)。
- β_g = β / m，使用固定等权 union bound，不提供权重调节。
- 对 Bernoulli 纳入概率 π，选择满足下式的最小 π：

~~~text
Pr[Binomial(M_min, π) ≥ 1] ≥ 1 - β_g

因此：
π = min(1, 1 - β_g^(1 / M_min))
~~~

- 不再增加“至少 3 条”或“至少 10 条”的经验 cutoff；比例区间负责表达证据量是否足够。
- 当 π = 1 时执行全量 census，不把它当异常。

**Pass 2：Deterministic Bernoulli inclusion**

- 为每个 fragment 计算 128-bit domain-separated BLAKE3 selection key；
- key 映射到 [0,1)，小于 π 即纳入；
- PE 两端共享同一个 key；
- 纳入 fragment 立即进入 Bloom 与有界比对队列，不构建 sequence reservoir；
- 完成后记录实际 selected fragments n 和实际纳入率 n/N。

### 5.3 确定性与重现性

- 相同 FASTQ bytes、profile、index 在任意线程数下产生相同纳入集合和 CSV。
- fragment key 使用完整 read identifier 与输入 ordinal，保留重复记录的多重性。
- FASTQ 物理重排会改变 input digest 和 analysis identity；不声称重排后样本不变。
- 不使用瞬时 RSS、CPU、输入吞吐或 Bloom 命中数调整 π。
- 不允许观察到 target hits 后再扩大、缩小或停止抽样。

### 5.4 比例估计与区间

Pass 2 的 hash inclusion 对所有 fragments 使用相同 π。条件于实际样本量 n，纳入集合是有限总体的简单随机样本。

对 group g：

~~~text
x_g = selected fragments 中归属于 g 的数量
p_hat_g = x_g / n
~~~

通过反演 Hypergeometric(N, M_g, n) 得到 M_g 的有限总体区间，再除以 N 得到 P_g 区间。

- 区间对固定 target family 使用 α_g = α / m 的 simultaneous Bonferroni coverage。
- 不对 observed rows 重新定义 multiplicity denominator。
- 不用普通独立 read-end 二项区间。
- 不输出 p value、adjusted p value 或 q value。
- interval_lower、interval_upper、interval_level 和 interval_method 始终原样报告；不从区间宽度派生二元类别。

### 5.5 Bloom 的职责

Bloom index 只包含 TARGET k-mers；它不包含 HOST，也不是 host-removal index。

实际语义：

1. 抽中 fragments 才进入 Bloom；
2. Bloom-negative fragments 不进入昂贵比对；
3. Bloom-positive fragments 对 HOST + TARGET 联合索引竞争比对；
4. host 通过竞争得分抑制 target 错配。

需要保留并报告：

- k-mer length；
- Bloom fill fraction 和理论 k-mer false-positive rate；
- fragment gate pass count；
- 低复杂度/短 read 导致的不可评价计数；
- profile digest。

k=21 可作为当前算法候选值，但必须从 CLI 删除、进入唯一 profile，并通过“禁用 gate 的小规模 exhaustive alignment”反事实测试证明其工作量收益和可接受偏差。验证失败时修改 profile 并重新执行完整评估，不增加第二条 gate 路径。

## 6. 竞争比对与归属

### 6.1 索引角色

用户合同只允许 HOST 和 TARGET genomes，因此 v0.5 的生产索引只包含：

- Role::Host；
- Role::Target。

删除：

- Role::Decoy；
- Role::Contaminant；
- synthetic decoy 生成与外部 decoy 输入；
- contaminant CLI 输入；
- 对应 manifest、report 和测试字段。

未知污染和未提供近缘序列必须作为限制报告，不能由 synthetic decoy 伪装成已建模背景。

### 6.2 Fragment adjudication

每个 fragment 形成一个 FragmentAlignmentEvidence：

1. 收集 R1/R2 的 host 与 target chains；
2. 使用冻结、长度归一化的 alignment evidence 排序；
3. target 仅在最佳可接受证据不被 host 支配时进入 target support set；
4. 同一 ReferenceGroup 内多个 member 命中保留 group-level support；
5. 跨 group 同等可接受证据标记为 unresolved；
6. host 优势或 host-target 不可区分标记为 confounded；
7. 不再使用绝对 AS margin=12、NM≤8 或 accession 级硬分配。

对齐规则必须在 profile 中有一个所有者，不允许 report、evidence 和 align 模块分别重判。

### 6.3 分布证据

covered bases、breadth、occupied windows、strand、split 和 discordant 继续作为审阅诊断量，但：

- 不再参与 PASS/FAIL cutoff；
- 不合成为 confidence score；
- 不因“看起来更合理”而修改比例估计；
- 只在 target signal rows 中输出；
- integration evidence 与普通 target fraction 分开命名和解释。

## 7. 状态模型

不再使用“可信”“灰区”“PASS”“WEAK”“NOT_SIGNIFICANT”等混合语义。结果拆为 analysis、evidence、attribution 三个正交轴。

### 7.1 analysis_status

| 状态 | 含义 |
|---|---|
| CONFORMANT_COMPLETE | 两遍输入、抽样、比对和三个产物均完整 |
| CONFORMANT_WITH_LIMITATIONS | 运行完整，但存在已量化 exclusions/overflow |
| RESOURCE_TERMINATED_INCONCLUSIVE | 资源边界终止，不能解释为零信号 |
| INPUT_TERMINATED_INCONCLUSIVE | FASTQ 截断、配对错误或第二遍不一致 |
| INVALID | profile/index/input digest 或内部不变量失败 |

### 7.2 evidence_status

| 状态 | 含义 |
|---|---|
| REFERENCE_SIGNAL_OBSERVED | 至少一个纳入 fragment 对该 group 产生可接受 target evidence |
| REFERENCE_SIGNAL_NOT_OBSERVED | 完整运行中未观察到，不等于 target absent |
| INDETERMINATE_EVIDENCE | 证据存在但跨 group 或 host 冲突无法解析 |
| NOT_EVALUABLE | 该 group 或所需参考合同无法评价 |

report.csv 只写 observed/indeterminate target rows；零信号 groups 由 RUN 行中的 family size、profile digest 和 HTML 汇总表达。

### 7.3 attribution_status

- RESOLVED_TO_REFERENCE_GROUP
- AMBIGUOUS_WITHIN_GROUP
- UNRESOLVED_ACROSS_GROUPS
- CONFOUNDED_WITH_HOST

区间数值是唯一精度表达，不映射为“可信/灰区”或任何二元状态。

## 8. 输出合同

### 8.1 目录合同

viroflash run --out 指向一个新的样本运行目录。成功时原子地产生且仅产生：

~~~text
report.html
report.csv
perf.json
~~~

不再生成：

- result.json；
- result.tsv；
- perf.tsv；
- 旧前缀式 sample.html、sample.csv、sample.perf.json；
- compatibility alias、symlink、转换器或双写。

索引目录的 ref.mmi、bloom.bin、manifest 等不属于样本运行输出三文件约束。

### 8.2 单一 report model

EvidenceReport 是唯一科学结果对象：

- report.csv 直接序列化 EvidenceReport；
- report.html 从同一个 EvidenceReport 渲染；
- HTML 不读取 CSV 后重算；
- perf.json 不复制生物学结果；
- 任何 HTML/CSV 数值差异都是测试失败。

### 8.3 report.csv

采用单一 RFC 4180 宽表，record_type 只有两种：

**RUN：恰好一行**

- schema_id、record_type、sample_id；
- analysis_status、reason_codes；
- input_mode、input_fragments、selected_fragments、selection_probability；
- minimum_relevant_fraction、familywise_miss_probability；
- interval_level、target_family_size；
- prescreen_passed_fragments、aligned_fragments、unassigned_fragments；
- profile_digest、index_digest、input_digest；
- read_ends_per_fragment。

**TARGET_SIGNAL：每个 observed 或 indeterminate group 一行**

- target_group_id、representative_id、member_ids；
- evidence_status、attribution_status；
- supporting_selected_fragments；
- selected_fragment_denominator；
- attributed_fragment_fraction；
- interval_lower、interval_upper、interval_level、interval_method；
- estimated_input_supporting_fragments；
- covered_bases、representative_length、coverage_fraction、occupied_windows；
- host_confounded_fragments、cross_group_ambiguous_fragments；
- integration_status、split_events、discordant_fragments；
- limitation_codes。

删除旧字段：

- reads、n_plain、depth_rpm、end_rpm；
- discovery_reads、validation_read_sides 和 fold 字段；
- p_value、q_value、adjusted_p_value、pi0、poisson_p、nb_p；
- coverage_min、min_distributed_windows、model_adjusted_p_max；
- confidence、confidence_basis、evidence_strength；
- bccp、decaf_grade、空 taxid；
- decoy/background stratum 字段；
- notes key-value blob；
- 顶层与 candidate 内重复字段。

### 8.4 report.html

自包含、离线、无外部请求，只保留五区：

1. **Interpretation boundary**：明确这是文库中 profile-attributed fragment fraction，不是临床结论。
2. **Run integrity**：analysis status、profile/index/input digests、输入 fragment 和抽样设计。
3. **Observed target signals**：按 group 展示比例、simultaneous interval 和三轴状态。
4. **Evidence detail**：归属歧义、host conflict、coverage/windows 和 integration 诊断。
5. **Methods and limitations**：两遍抽样、Bloom 职责、竞争比对、未知污染与不可识别边界。

HTML 可以搜索、排序和过滤，但不得删除底层 EvidenceReport 中的 target rows 或创造新的统计值。

### 8.5 perf.json

只负责执行和资源遥测：

- schema_id、status、sample_id；
- wall_time_ms、configured_threads；
- process_cpu_time_ms、peak_rss_bytes、read_bytes、written_bytes；
- pass1_count、pass2_sample_prescreen_align 和 report_write 的 stage wall time；
- input_fragments、selected_fragments、prescreen_passed_fragments、aligned_fragments；
- telemetry_status 和 sample_errors。

删除 PID、virtual memory、系统总内存/可用内存、系统 CPU、重复 thread budget、低频采样推导的 stage CPU 均值等弱解释字段。

失败时 perf.json 可单独存在并标记 ERROR；report.csv/report.html 不得以成功形状部分写出。

## 9. 模块、文件和命名重构

### 9.1 目标模块

| 当前 | 目标 | 职责 |
|---|---|---|
| src/fastq.rs | src/fastq_input.rs | 两遍 FASTQ 校验、计数和稳定迭代 |
| src/sampling.rs | src/sampling_design.rs | π 计算、selection key、纳入审计 |
| src/prescreen.rs | src/kmer_gate.rs | TARGET Bloom 工作量门 |
| src/align.rs | src/competitive_alignment.rs | HOST/TARGET chains 与 fragment adjudication |
| src/group.rs | src/reference_group.rs | 索引期固定 target groups |
| src/index.rs | src/reference_index.rs | index/profile/reference digest |
| src/cluster.rs | src/integration_evidence.rs | integration-site 证据 |
| src/report.rs + src/visual_report.rs | src/report.rs | typed model、CSV、HTML 唯一序列化所有者 |
| src/perf.rs | src/performance_report.rs | perf.json |
| src/lib.rs | src/lib.rs + src/pipeline.rs | lib.rs 仅公开入口；pipeline 负责编排 |

新增：

- src/analysis_profile.rs：唯一编译期科学 profile 与 digest。
- src/evidence.rs：fragment evidence 聚合、区间和三轴状态。

删除：

- src/decoy.rs；
- src/equivalence.rs；
- src/stats.rs；
- src/hash.rs，其小型 digest helper 收归实际所有者。

不保留旧模块 re-export、deprecated alias 或 feature flag 双路径。

### 9.2 类型与函数命名

核心类型：

- AnalysisProfile
- ReferenceGroup
- InputCensus
- SamplingDesign
- FragmentSelection
- FragmentAlignmentEvidence
- TargetGroupEvidence
- FinitePopulationInterval
- EvidenceReport
- AnalysisStatus
- EvidenceStatus
- AttributionStatus

核心函数边界：

- census_fastq
- derive_sampling_design
- fragment_selection_key
- include_fragment
- passes_target_kmer_gate
- align_fragment_competitively
- adjudicate_fragment
- accumulate_group_evidence
- finite_population_interval
- build_evidence_report
- write_report_csv
- write_report_html
- write_perf_json

删除 Pair、read-side、depth RPM、fold、candidate confidence、q compatibility alias 等不准确概念。

### 9.3 CLI 收敛

目标正常路径：

~~~text
viroflash index \
  --host-fa HOST.fa \
  --target-fa TARGET.fa \
  --out INDEX_DIR \
  [--threads N]

viroflash run \
  --r1 SAMPLE_R1.fastq.gz \
  [--r2 SAMPLE_R2.fastq.gz] \
  --index INDEX_DIR \
  --out SAMPLE_REPORT_DIR \
  [--threads N]
~~~

删除：

- run 时自动临时建索引；
- --k；
- --decoy-fa、--decoy-ani、--decoy-per-layer、--decoy-seed；
- --contam-fa；
- 所有统计与报告阈值参数；
- 旧参数 alias。

threads 是结果不变量；它不进入 scientific profile。

## 10. 实施顺序

### Phase 0：冻结合同与基线

产出：

- AnalysisProfile 字段及符号 δ、α、β 的产品决策记录；
- ReferenceGroup ledger；
- 新 CLI 与三文件字段表；
- 127 样本只读 evaluation manifest；
- 当前 68/59 历史基线快照。

验收：

- 参数在首次查看 v0.5 结果前冻结；
- truth mapping、样本路径、reference paths 和 checksums 固定；
- Cargo.toml version 仍为 0.3.0。

此阶段不改生产路径。

### Phase 1：建立 fragment 级证据所有权

- 新建 analysis_profile.rs 与 evidence.rs；
- 固定 ReferenceGroup，统一 SE/PE fragment 计数；
- 将归属、歧义、host conflict 和分布诊断移入 evidence.rs；
- 报告模型先建立但暂不切换文件输出。

验证：

- PE 两端不双计；
- group 内 member ambiguity 不变成 accession 归属；
- 跨 group tie 不硬分配；
- 所有 evidence reason code 可复算。

### Phase 2：替换固定 reservoir

- Pass 1 精确 InputCensus；
- 按公式生成 SamplingDesign；
- Pass 2 确定性 Bernoulli inclusion；
- 流式 Bloom 和有界 alignment queue；
- 删除 PairReservoir、DEFAULT_SAMPLE_PAIRS 和 EvidenceFold。

验证：

- π 公式与枚举/高精度计算一致；
- empirical inclusion rate 与 π 相符；
- 相同输入跨线程字节级 CSV 一致；
- π=1 与全量路径一致；
- gzip pending 与 alignment queue 有严格字节上界；
- 无 sequence reservoir。

### Phase 3：删除旧统计路径

- 生产索引只保留 HOST/TARGET；
- 删除 synthetic decoy、contaminant、discovery/validation；
- 删除 conditional Poisson、BH/q 和三 cutoff 决策；
- 删除 decoy.rs、equivalence.rs、stats.rs；
- 删除普通未使用 aligner 路径和重复 hash 模块。

验证：

- 无 p/q/PASS 代码和字段；
- 缺少背景不会产生 positive-shaped 统计结论；
- 与旧结果的每个变化按 evidence ownership 解释，而不是追求阈值复现。

### Phase 4：切换三文件报告

- 建立唯一 EvidenceReport；
- 原子写入 report.csv/report.html；
- perf.json 与科学结果分离；
- 删除旧 JSON/TSV 和 prefix 输出；
- 重写 README.md、PRODUCT.md 和 CLI help。

验证：

- 成功目录恰好三个文件；
- CSV/HTML 对每个 target row 数值一致；
- 空 signal 时有 RUN 行、无伪 candidate；
- error 时不存在成功形状 CSV/HTML；
- 不存在兼容 alias、converter 或双写。

### Phase 5：统计反事实与规模测试

- 构造有限总体，枚举/模拟验证 sampling power、估计偏差和 simultaneous interval coverage；
- 对小 FASTQ 同时运行 gate-enabled 与 exhaustive alignment，量化 Bloom 引入的 profile 差异；
- 线程、chunk、gzip member 顺序只改变性能，不改变结果；
- 运行从小到大输入的内存曲线，证明 RSS 不随 selected sequence bytes 线性增长。

这是 127 个真实样本无法替代的统计 oracle，不属于训练或学习。

### Phase 6：127 样本真实 E2E

严格按第 11 节执行。任何 profile 或结果规则修改都使本轮作废，必须用新 analysis identity 从头重跑 127 个样本。

## 11. 68 内部 + 59 外部真实 E2E 方案

### 11.1 当前数据可用性预检

已在 2026-09-07 对 ~/data 做只读盘点：

**内部 68**

- truth：~/data/viroflash/dna_virus_68_20260827/panel_truth.json；
- FASTQ：~/data/samples/fastviro_testset；
- 68 个 PE 样本需要 136 个文件，已匹配 136/136；
- 压缩总量 954,071,909,453 bytes，即 888.55 GiB；
- references：
  - ~/data/viroflash/dna_virus_68_20260827/refs/host_grch38.fa；
  - ~/data/viroflash/dna_virus_68_20260827/refs/dna_virus_genome.fasta；
- 历史脚本默认 /home/wubw/__viroflash/.dataset/testset 已过期；新 harness 禁止保留该路径 fallback。

**外部 59**

- labels：~/data/viroflash-benchmarks/results_20260827/label_comparison.tsv；
- FASTQ：~/data/viroflash-benchmarks/reads；
- 59/59 runs 可用，共 82 个 SE/PE FASTQ；
- 82/82 已有 metadata/FASTQ_SHA256SUMS；
- cohorts：
  - GSE147507 SARS：3 mock + 3 SARS-CoV-2；
  - GSE147507 RSV：3 mock + 3 RSV；
  - GSE91065 HPV：47 runs。

数据已具备执行条件，但内部 136 个 FASTQ 尚需在 Phase 0 生成冻结 checksum manifest。

### 11.2 评价角色与偏倚边界

- 68 内部样本曾参与当前系统开发和复核，只能作为 development/regression evidence，不是独立盲法验证。
- 59 外部样本也已有历史结果和人工检查，只能作为 external generalization regression，不是 pristine holdout。
- 两套数据都不得用于选择 δ、α、β、k 或归属规则。
- 如结果不理想，修改 profile 后必须创建新 analysis identity 并全量重跑；不能只挑失败样本重跑。
- 真值标签只评价标签映射范围，不扩展为 21,790 references 的临床 sensitivity/specificity。

### 11.3 冻结 evaluation manifest

运行前生成开发期 evaluation manifest，至少包含：

- dataset_id、cohort、sample_id；
- SE/PE mode、R1/R2 absolute path；
- FASTQ sha256、compressed bytes；
- host/target FASTA path 与 sha256；
- expected group 或 mock；
- label provenance；
- evaluability status；
- 已知 ambiguity notes；
- profile digest、binary digest、index digest；
- planned threads/jobs。

manifest 和 scorer 在首次 v0.5 E2E 前冻结。它们是测试工件，不是新增的生产输出。

### 11.4 真实端到端边界

每个 cohort 必须执行：

1. 从 HOST/TARGET FASTA 构建新索引；
2. 验证 index digest 和 profile digest；
3. 从原始 FASTQ 执行 Pass 1；
4. 生成 SamplingDesign；
5. 从原始 FASTQ 执行 Pass 2、Bloom、竞争比对和 evidence；
6. 写入三个新文件；
7. 独立解析 CSV/HTML/perf；
8. 汇总标签一致性、证据状态和性能。

禁止：

- 读取历史 result JSON/TSV 作为 v0.5 输入；
- 复用旧索引；
- 跳过 Pass 1；
- 只测试报告转换；
- 根据样本标签改变 profile 或 threads；
- 失败后自动重试并覆盖第一次失败证据。

### 11.5 内部 68 评价矩阵

历史标签构成：

- EBV：17；
- HBV：3；
- HPV16：2；
- HPV18：2；
- negative：44。

已知边界：

- 历史结果为 expected target 22/24、negative 44/44 无 PASS；
- 两个 HPV18 标签样本在序列复核中 HPV18 为零、EBV 很强，不能机械记为算法 false negative；
- truth 只覆盖五类 panel labels，不覆盖整个病毒参考库；
- EBV 多 accession 必须在 ReferenceGroup 层计一次。

每个样本分类为：

- LABEL_CONCORDANT_SIGNAL；
- LABEL_CONCORDANT_NO_SIGNAL；
- LABEL_DISCORDANT_WITH_SEQUENCE_EVIDENCE；
- LABEL_DISCORDANT_UNRESOLVED；
- NOT_EVALUABLE。

不得把 label discordance 自动等同于算法错误；所有 discordance 必须保留 read-level evidence 摘要供人工复核。

### 11.6 外部 59 评价矩阵

历史基线：

- 36 PASS；
- 17 BELOW_THRESHOLD；
- 6 无候选；
- 6 mock 均无 PASS，但 4 个 mock 有 trace BELOW；
- 4 个 HPV16 not reported 对应 2 个 biological samples；全文件 40-mer probe 未发现 HPV16，存在标签与 RNA-seq 内容不一致；
- SRR5090635 有历史 SE/PE 执行不一致，v0.5 manifest 必须冻结正确 input mode。

评价维度：

- expected group 是否出现在 TARGET_SIGNAL；
- expected group 的 attribution/quantitation status；
- mock 中 observed unexpected signals 的数量、比例和区间；
- HPV16/HPV18 group discrimination；
- wrong-group resolved signal；
- ambiguous/cross-group signal；
- SE 与 PE input contract 完整性。

trace signal 不再被隐藏为 BELOW_THRESHOLD；它应以点估计和 simultaneous interval 原样展示。

### 11.7 硬性 E2E 验收

以下全部满足才算执行完整：

1. 68/68 内部样本和 59/59 外部样本都有终态，无静默跳过。
2. 所有成功样本目录恰好包含 report.html、report.csv、perf.json。
3. 三个文件均可由独立 parser 解析；HTML/CSV 数值逐 target 一致。
4. 所有样本使用同一 profile digest；各 cohort 使用冻结 index digest。
5. FASTQ Pass 1 与 Pass 2 的 fragment count、pair IDs 和 input digest 一致。
6. 没有 legacy JSON/TSV/perf TSV、兼容 alias 或旧 index。
7. 无 crash、OOM、panic、死锁、无界 pending growth 或 success-shaped partial report。
8. 所有 label discordance 都进入冻结的 adjudication table，无未审查变化。
9. 生产运行期间没有结果相关参数 override。
10. Cargo.toml version 保持 0.3.0。

### 11.8 科学评价指标

必须报告计数和 Wilson interval，而不是只给百分比：

- expected-group observed rate；
- expected-group interval lower/upper summaries；
- mock unexpected-signal rate；
- wrong-group resolved-signal rate；
- ambiguous attribution rate；
- not-evaluable rate；
- internal label concordance，含/不含已知 HPV18 ambiguity 两种口径；
- external 三 cohort 分层结果；
- 与旧系统逐样本 paired difference table。

这些指标用于描述变化，不用于反向调整 profile。

真实数据没有已知病毒 fragment fraction，因此不能用 68/59 验证：

- abundance bias；
- interval coverage；
- LoD/LoQ；
- absolute quantitation。

这些只能由 Phase 5 的已知有限总体、spike-in 或独立 wet-lab controls 验证。本版本只要求已知有限总体模拟，不凭空增加临床验证声明。

### 11.9 性能评价

内部历史基线：

- 68 样本总输入 888.55 GiB；
- batch wall 约 28.9 min；
- 单样本 wall median 147.3 s、P95 437.7 s、max 685.9 s；
- peak RSS median 25.26 GiB、P95 25.83 GiB、max 26.14 GiB；
- index build 1507.0 s、45.53 GiB。

外部历史基线：

- 单 run wall median 34.7 s、P95 53.7 s、max 106.7 s；
- peak RSS 约 10.7 GiB。

v0.5 必须报告：

- index build wall/peak RSS；
- 两遍读取各自 wall 和 read throughput；
- sample wall median/P95/max；
- peak RSS median/P95/max；
- Bloom pass rate；
- aligned fragments；
- batch wall 和并发配置；
- 相对旧基线的 paired ratio。

性能收敛规则：

- 结果正确性和内存有界是硬门；
- 相同 index 下单样本 peak RSS 不得高于对应旧基线最大值；
- 输入规模增长时，除 index 与有界队列外，RSS 不得与 selected sequence bytes 线性增长；
- wall time 若回退，必须定位到 Pass 1/Pass 2/Bloom/alignment 的具体 stage；
- 不设拍脑袋百分比容限：任何未被明确接受的系统性 wall-time 回退都视为未收敛，先优化结果不变量路径或由产品所有者记录接受；
- 不允许以减少 π、降低 δ、放宽 α/β 或跳过样本换取性能通过。

### 11.10 重现性子集

在完整 127 样本之外，冻结一个覆盖以下情况的最小子集：

- internal：最大文件、最小文件、EBV、HBV、HPV16、HPV18 ambiguity、negative；
- external：SARS positive/mock、RSV positive/mock、HPV16、HPV18、SRR5090635。

对子集执行：

- threads 1/2/4/8；
- reusable index 重复加载；
- 两次独立运行；
- plain/multi-member gzip 等价测试（适用时）。

report.csv 必须字节级一致；perf.json 允许时间和资源字段变化。

## 12. 明确收敛条件

只有以下全部满足，v0.5.0 方案才算实现完成：

### 合同收敛

- 一个 AnalysisProfile、一个 ReferenceGroup ledger、一个 EvidenceReport。
- 一个生产索引路径和一个样本运行路径。
- 生产 CLI 只保留输入、输出和资源参数。
- 旧名称、旧字段、旧输出和旧 index 路径已删除，没有 alias。

### 统计收敛

- sample selection power 公式通过枚举和模拟。
- finite-population estimator 无偏性和 simultaneous interval coverage 达到冻结合同。
- fragment 是唯一计数单位。
- 无 synthetic-decoy p/q、无 PASS、无数据后阈值调整。
- EM/latent abundance 明确为不可识别且未实现。

### 执行收敛

- gzip pending、reader、alignment 和 writer 队列均有字节上界。
- 相同 input/profile/index 跨线程 report.csv 一致。
- resource termination 只产生 inconclusive/error，不产生零信号成功报告。
- 成功目录恰好三个文件。

### 真实数据收敛

- 68 + 59 全部真实重跑完成。
- 所有 E2E 硬门通过。
- paired scientific/performance comparison 已生成。
- 所有 call/status 变化已有证据解释。
- 没有在看到结果后修改 profile；若修改则整轮作废并重新开始。

### 代码库收敛

- cargo fmt --all -- --check；
- cargo check --locked；
- cargo clippy --all-targets --all-features --locked -- -D warnings；
- cargo test --all-targets --locked；
- cargo test --release --locked；
- cargo build --release --locked；
- 完整 diff 审查确认无版本变更、兼容胶水、重复所有者、隐藏阈值和用户数据改动。

## 13. 删除/保留/延期清单

### 本版本删除

- 固定 524,288 reservoir；
- discovery/validation fold；
- runtime equivalence hypotheses；
- synthetic decoy 与 contaminant 角色；
- Poisson/BH/q/PASS 判定；
- coverage/windows cutoff；
- read-end RPM；
- JSON/TSV/perf TSV；
- auto-build run path；
- 所有旧 schema/CLI/module aliases。

### 本版本保留

- 严格 FASTQ 校验；
- deterministic hash selection；
- TARGET Bloom + SDUST 工作量门；
- minimap2 ALL_CHAINS HOST/TARGET 竞争；
- fixed ReferenceGroup；
- fragment-level coverage、split、discordant 和 integration diagnostics；
- reusable index；
- 有界并行执行；
- performance telemetry。

### 有证据后才重开

- confidence sequence：只有真实 streaming/early-stop 合同且固定设计两遍成本不可接受时。
- EM：只有 group 内丰度成为明确需求，并有可识别 mock mixtures 与独立验证时。
- HIBF/COBS/KMCP 分箱：只有 127 样本显示 alignment 是主导瓶颈且单 Bloom 不足时。
- contaminant reference：只有用户输入合同明确增加且有真实参考来源时。
- clinical LoD/LoQ：只有阴性对照、spike-in、重复梯度和 wet-lab protocol 时。

## 14. 证据来源与采用边界

- Horvitz–Thompson 设计估计说明纳入概率必须进入无偏估计；本方案采用等概率设计，不引入不必要的加权分层复杂度：
  https://doi.org/10.1080/01621459.1952.10483446
- Vitter reservoir 与 bottom-k 证明固定容量均匀抽样的统计性质；本方案保留 hash-randomization 思想，但删除保存完整序列的固定 reservoir：
  https://doi.org/10.1145/3147.3165
  https://arxiv.org/abs/1303.5479
- Waudby-Smith/Ramdas 与 Howard 等提供 anytime-valid confidence sequence；因当前没有边读边停止的生产合同，v0.5 明确延期而不是预先实现：
  https://proceedings.neurips.cc/paper/2020/hash/e96c7de8f6390b1e6c71556e4e0a4959-Abstract.html
  https://doi.org/10.1214/20-AOS1991
- Salmon/kallisto 的 equivalence-class 不确定性说明不能将多重归属硬分配；本方案只采用 set/group 语义，不移植 RNA effective-length 或 EM：
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5600148/
  https://pachterlab.github.io/kallisto/manual
- Kraken2、ganon、KMCP 说明 compact/minimizer/HIBF/COBS 可减少候选检索成本；没有当前瓶颈证据前不引入新索引层：
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6883579/
  https://pmc.ncbi.nlm.nih.gov/articles/PMC12267982/
  https://academic.oup.com/bioinformatics/article/39/1/btac845/6965021
- PathoScope/GRAMMy 说明竞争比对与软分配可用于混合样本；当前 unknown component 与 emission model 不可识别，因此拒绝生产 EM：
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4164323/
  https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0027992
- khmer digital normalization 的目标是压缩 assembly 工作量，不提供病原检出保证，因此不用于 target detection 抽样：
  https://pmc.ncbi.nlm.nih.gov/articles/PMC4608353/

## 15. 对既有文档的处理

- .local/UPGRADE_PLAN.md 保留为历史研究记录，但其 synthetic decoy、Storey q、EB shrinkage、Group-walk、固定 coverage/q gates 和推测性线程比例不进入 v0.5。
- PRODUCT.md 当前保留旧 JSON/TSV/synthetic-decoy 语义；实施 Phase 4 时必须直接重写，不加兼容段落。
- README.md 当前声明 JSON 是 source of truth；实施后改为 EvidenceReport 为内存 source of truth、report.csv 为唯一机器可读科学产物。
- 本方案不修改 Cargo.toml 版本号，也不授权 commit、tag、push 或 release。
