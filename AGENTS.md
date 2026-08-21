# AGENTS.md

## 项目概览

`viroflash` 是一个 Rust 2021 单 crate 命令行项目。它读取单端或双端
FASTQ/FASTQ.gz，以及宿主、目标、诱饵和污染四类 FASTA，执行 k-mer 预筛、
minimap2 竞争比对、候选聚合和统计判定，输出候选级 JSON 与 TSV。

主程序不调用外部 CLI；比对通过 `minimap2` Rust 绑定及其原生 FFI 完成。
保持依赖少、流式输入、有界并发、确定性聚合和可审计的参数来源。

## 目录与职责

```text
Cargo.toml / Cargo.lock    包配置与锁定依赖
src/main.rs                CLI 解析、命令分发和退出码
src/lib.rs                 管线编排、证据聚合和候选决策
src/reference.rs           FASTA 解析、复合参考和 minimap2 索引
src/fastq.rs               FASTQ/FASTQ.gz 解析与配对校验
src/prescreen.rs           canonical k-mer Bloom 与 SDUST 门控
src/align.rs               minimap2 竞争比对和有序并行
src/cluster.rs             位点聚类与去重
src/stats.rs               Poisson/NB 与 q 值计算
src/decoy.rs               确定性 SNP-only 诱饵生成
src/report.rs              JSON/TSV 序列化
tests/smoke.rs             确定性合成端到端测试
scripts/                   可移植的合成数据工具
```

`Cargo.toml`、`Cargo.lock`、当前 `src/` 和测试是构建及运行行为的事实来源。
阈值应由拥有该阶段的模块集中定义，不要散落在调用点。

## 仓库边界

- `.local/` 只保存本机真实测试、历史材料、运行结果、日志和缓存，不得提交。
- 不提交 FASTQ/FASTA 实际数据、样本真值、报告产物、索引、构建产物或工作目录。
- `target/`、`.local/` 和 `*.work/` 必须保持在 `.gitignore` 中。
- 公开脚本必须可在 clean checkout 中运行，不得依赖本机绝对路径或内部资料。
- 不把开发工具交互记录、过程性报告或内部基础设施信息写入版本库。

## 构建与验证

在仓库根目录使用锁文件：

```bash
cargo fmt --all -- --check
cargo check --locked
cargo clippy --all-targets --all-features --locked -- -D warnings
cargo test --all-targets --locked
cargo test --release --locked
cargo build --release --locked
```

查看 CLI：

```bash
cargo run --locked -- --help
cargo run --locked -- version
```

## 修改约定

- 使用默认 rustfmt，保持中文模块注释、代码注释和错误信息的现有风格。
- 标识符遵循 Rust `snake_case` 和 `CamelCase` 约定。
- 生产路径返回 `Result<T, String>`；输入错误不得 panic、吞错或返回成功形状的空结果。
- 测试中可以使用 `unwrap`/`expect`，生产输入路径避免使用。
- 优先标准库和现有依赖；新增依赖前先证明现有能力不足。
- 不用 `std::process::Command` 包装 minimap2、gzip 或其他外部工具替代现有库路径。
- 保持 CLI 参数、JSON schema、TSV 列顺序、fragment/read 口径和 0-based 半开坐标稳定。
- 保持固定 seed、稳定输入顺序、有序并行汇总和候选稳定排序。
- 修改并发通道、发送/接收顺序或排空逻辑前，先阅读 `src/align.rs` 的并发契约。
- 做最小、相关、可逆的修改；不要顺带重构、格式化无关文件或增加兼容层。

## 测试范围

- CLI/参数：相关单元测试、help/version 和至少一个错误路径。
- FASTA/FASTQ/gzip：相关单元测试与 smoke；并行解压还需验证单线程等价。
- 门控、比对、聚类或线程：相关模块测试和 smoke。
- 统计或候选决策：`stats`/`lib` 测试和 smoke。
- JSON/TSV：`report` 测试、JSON 解析和两种格式字段一致性。

## Git 约定

使用 Conventional Commits。每个 commit 只包含一个自相关目的，提交前检查
staged diff 并运行最窄验证。导入既有代码使用 `chore(import): ...`，不要伪造
开发历史。未经明确授权，不执行 push、tag、remote 修改或历史改写。

完成前检查 `git status`、tracked files、忽略规则和完整验证结果，确保内部资料、
构建产物及无关更改未进入版本库。
