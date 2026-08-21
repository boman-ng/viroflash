//! viroflash 索引目录：构建、加载与校验。
//!
//! 索引目录布局（`viroflash index --out <dir>` 产物）：
//! ```text
//! <dir>/
//!   ref.mmi        minimap2 sr 索引（minimap2 绑定 with_index(fa, Some(out)) 落盘）
//!   bloom.bin      目标+诱饵 canonical k-mer Bloom（版本化二进制）
//!   manifest.json  版本、k、角色→contig 元数据、诱饵来源参数、构建时间、来源 checksum
//!   decoys.fa      自动构建诱饵时的生成产物（+ decoys.tsv 元数据报告）
//! ```
//! 加载路径只读 manifest/bloom/ref.mmi，不需要序列本身；来源 FASTA 以路径 + BLAKE3
//! 记录在 manifest 中，保证索引可审计（STAR genomeParameters.txt / salmon
//! versionInfo.json 模式的先例）。构建采用「临时目录 + 原子改名」，避免并发竞态
//! 与半成品索引。
//!
//! JSON 为手工序列化/解析（延续 report.rs 零 serde 依赖约定），解析器仅支持本
//! manifest 需要的 JSON 子集。

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::decoy::{self, DecoyOptions};
use crate::hash;
use crate::prescreen::{self, KmerBloom};
use crate::reference::{self, Contig, ContigMeta, Role};
use crate::report::json_escape;

/// 索引格式版本（manifest 与 bloom.bin 共用同一版本号）。
pub const FORMAT_VERSION: u32 = 1;
pub const MANIFEST_NAME: &str = "manifest.json";
pub const MMI_NAME: &str = "ref.mmi";
pub const BLOOM_NAME: &str = "bloom.bin";
pub const DECOYS_FA_NAME: &str = "decoys.fa";
pub const DECOYS_TSV_NAME: &str = "decoys.tsv";

/// 索引构建参数（`viroflash index` 命令与 `run` 自动构建共用同一入口）。
#[derive(Debug, Clone)]
pub struct IndexOptions {
    pub host_fa: PathBuf,
    pub target_fa: PathBuf,
    pub contam_fa: Option<PathBuf>,
    pub decoy_fa: Option<PathBuf>,
    /// 未提供 decoy_fa 时自动生成诱饵的 ANI 层（百分比整数）。
    pub decoy_anis: Vec<u8>,
    pub decoy_per_layer: usize,
    pub decoy_seed: u64,
    pub out_dir: PathBuf,
    pub k: usize,
    /// minimap2 索引构建线程数。
    pub threads: usize,
}

impl Default for IndexOptions {
    fn default() -> Self {
        Self {
            host_fa: PathBuf::new(),
            target_fa: PathBuf::new(),
            contam_fa: None,
            decoy_fa: None,
            decoy_anis: decoy::DEFAULT_ANIS.to_vec(),
            decoy_per_layer: decoy::DEFAULT_PER_LAYER,
            decoy_seed: 0,
            out_dir: PathBuf::new(),
            k: prescreen::DEFAULT_K,
            threads: 8,
        }
    }
}

/// 诱饵来源（manifest 记录，保证零分布可审计）。
#[derive(Debug, Clone, PartialEq)]
pub enum DecoySource {
    File {
        path: String,
        blake3: String,
    },
    Generated {
        anis: Vec<u8>,
        per_layer: usize,
        seed: u64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct RefsInfo {
    pub host: FileInfo,
    pub target: FileInfo,
    pub contam: Option<FileInfo>,
    pub decoy: DecoySource,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FileInfo {
    pub path: String,
    pub blake3: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BloomInfo {
    pub n_inserted: u64,
    pub fill_frac: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndexManifest {
    pub format_version: u32,
    pub k: usize,
    pub created_at_unix: u64,
    pub refs: RefsInfo,
    pub contigs: Vec<ContigMeta>,
    pub bloom: BloomInfo,
}

/// 解析后的可用索引：构建与加载两条路径产出同一结构，
/// `run` 管线对其余逻辑完全一致（两路径结果等价的基础）。
#[derive(Debug)]
pub struct BuiltIndex {
    pub mmi_path: PathBuf,
    pub bloom: KmerBloom,
    pub roles: HashMap<String, Role>,
    pub contigs: Vec<ContigMeta>,
    pub manifest_blake3: String,
}

fn now_unix() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// 构建索引目录。`out_dir` 已存在时拒绝（不静默覆盖）；内部在临时目录构建
/// 完成后原子改名，失败时清理临时目录。
pub fn build_index(opt: &IndexOptions) -> Result<BuiltIndex, String> {
    validate_index_options(opt)?;
    if opt.out_dir.exists() {
        return Err(format!("索引目录已存在: {}", opt.out_dir.display()));
    }
    let part = PathBuf::from(format!(
        "{}.part.{}",
        opt.out_dir.display(),
        std::process::id()
    ));
    if part.exists() {
        // 上次同 pid 崩溃残留：仅清理自身专属的临时目录。
        let _ = std::fs::remove_dir_all(&part);
    }
    let result = build_index_into(&part, opt);
    if result.is_err() {
        let _ = std::fs::remove_dir_all(&part);
        return result;
    }
    std::fs::rename(&part, &opt.out_dir)
        .map_err(|e| format!("索引目录定稿失败 {}: {e}", opt.out_dir.display()))?;
    // part 目录已原子改名为最终目录：返回结构中的产物路径须指向最终位置。
    let mut built = result?;
    built.mmi_path = opt.out_dir.join(MMI_NAME);
    Ok(built)
}

fn build_index_into(part: &Path, opt: &IndexOptions) -> Result<BuiltIndex, String> {
    std::fs::create_dir_all(part)
        .map_err(|e| format!("无法创建索引目录 {}: {e}", part.display()))?;

    // 1. 诱饵：显式文件直接引用；否则按固定参数从目标自动生成（产物留在索引内，
    //    保证同一索引的零分布可复现、可审计）。
    let decoy_path: PathBuf;
    let decoy_source: DecoySource;
    match &opt.decoy_fa {
        Some(p) => {
            decoy_path = p.clone();
            decoy_source = DecoySource::File {
                path: p.display().to_string(),
                blake3: hash::blake3_file_hex(p)?,
            };
        }
        None => {
            let out = part.join(DECOYS_FA_NAME);
            let report = part.join(DECOYS_TSV_NAME);
            decoy::generate(&DecoyOptions {
                target_fa: opt.target_fa.clone(),
                out: out.clone(),
                anis: opt.decoy_anis.clone(),
                per_layer: opt.decoy_per_layer,
                seed: opt.decoy_seed,
                report: Some(report),
            })?;
            decoy_path = out;
            decoy_source = DecoySource::Generated {
                anis: opt.decoy_anis.clone(),
                per_layer: opt.decoy_per_layer,
                seed: opt.decoy_seed,
            };
        }
    }

    // 2. 复合参考 + minimap2 索引（构建在 build/ 子目录，落盘后移入根并丢弃 composite）。
    let build_dir = part.join("build");
    std::fs::create_dir_all(&build_dir)
        .map_err(|e| format!("无法创建构建目录 {}: {e}", build_dir.display()))?;
    let mut fastas = vec![
        (Role::Host, opt.host_fa.clone()),
        (Role::Target, opt.target_fa.clone()),
        (Role::Decoy, decoy_path.clone()),
    ];
    if let Some(c) = &opt.contam_fa {
        fastas.push((Role::Contaminant, c.clone()));
    }
    let (mmi_built, contigs) = reference::build_reference(&fastas, &build_dir, opt.threads)?;
    let mmi_path = part.join(MMI_NAME);
    std::fs::rename(&mmi_built, &mmi_path).map_err(|e| format!("移动索引文件失败: {e}"))?;
    let _ = std::fs::remove_file(build_dir.join("composite.fa"));
    let _ = std::fs::remove_dir(&build_dir);

    // 3. Bloom（目标+诱饵 canonical k-mer；与 run 的预筛门共用同一构建参数）。
    let bloom_refs: Vec<&Contig> = contigs
        .iter()
        .filter(|c| matches!(c.role, Role::Target | Role::Decoy))
        .collect();
    let bloom = KmerBloom::build(&bloom_refs, opt.k, prescreen::GATE_FPR).ok_or_else(|| {
        format!(
            "目标+诱饵 FASTA 中没有 ≥k={} 的有效 k-mer（全部序列过短或含非 ACGT 字符）",
            opt.k
        )
    })?;
    write_bloom(&part.join(BLOOM_NAME), &bloom)?;

    // 4. manifest（角色→contig 元数据 + 来源 checksum + 诱饵参数）。
    let roles: HashMap<String, Role> = contigs.iter().map(|c| (c.name.clone(), c.role)).collect();
    let metas: Vec<ContigMeta> = contigs.iter().map(ContigMeta::from).collect();
    let manifest = IndexManifest {
        format_version: FORMAT_VERSION,
        k: opt.k,
        created_at_unix: now_unix(),
        refs: RefsInfo {
            host: file_info(&opt.host_fa)?,
            target: file_info(&opt.target_fa)?,
            contam: opt.contam_fa.as_deref().map(file_info).transpose()?,
            decoy: decoy_source,
        },
        contigs: metas.clone(),
        bloom: BloomInfo {
            n_inserted: bloom.n_inserted,
            fill_frac: bloom.fill_frac(),
        },
    };
    let text = write_manifest(&manifest);
    let manifest_path = part.join(MANIFEST_NAME);
    std::fs::write(&manifest_path, &text).map_err(|e| format!("写 manifest.json 失败: {e}"))?;

    Ok(BuiltIndex {
        mmi_path,
        bloom,
        roles,
        contigs: metas,
        manifest_blake3: hash::blake3_hex(text.as_bytes()),
    })
}

fn file_info(p: &Path) -> Result<FileInfo, String> {
    Ok(FileInfo {
        path: p.display().to_string(),
        blake3: hash::blake3_file_hex(p)?,
    })
}

/// 加载既有索引目录并校验：存在性、格式版本、k 一致性、bloom 完整性。
pub fn load_index(dir: &Path, k: usize) -> Result<BuiltIndex, String> {
    if !dir.is_dir() {
        return Err(format!("索引目录不存在: {}", dir.display()));
    }
    let manifest_path = dir.join(MANIFEST_NAME);
    let text = std::fs::read_to_string(&manifest_path).map_err(|e| {
        format!(
            "读取 {} 失败（不是 viroflash 索引目录或缺少 manifest.json）: {e}",
            manifest_path.display()
        )
    })?;
    let manifest = parse_manifest(&text)?;
    if manifest.format_version > FORMAT_VERSION {
        return Err(format!(
            "索引格式版本 {} 高于当前支持的 {}，请用当前版本重新构建索引",
            manifest.format_version, FORMAT_VERSION
        ));
    }
    if manifest.k != k {
        return Err(format!(
            "索引 k={} 与 --k={} 不一致（Bloom 与索引均按 k 构建，需以同一 k 重建索引）",
            manifest.k, k
        ));
    }
    let mmi_path = dir.join(MMI_NAME);
    if !mmi_path.is_file() {
        return Err(format!("索引缺少 {}", MMI_NAME));
    }
    let bloom = read_bloom(&dir.join(BLOOM_NAME))?;
    if bloom.k != manifest.k {
        return Err(format!(
            "{} 与 {} 的 k 不一致（{} / {}），索引损坏，请重建",
            BLOOM_NAME, MANIFEST_NAME, bloom.k, manifest.k
        ));
    }
    let roles: HashMap<String, Role> = manifest
        .contigs
        .iter()
        .map(|c| (c.name.clone(), c.role))
        .collect();
    Ok(BuiltIndex {
        mmi_path,
        bloom,
        roles,
        contigs: manifest.contigs,
        manifest_blake3: hash::blake3_hex(text.as_bytes()),
    })
}

fn validate_index_options(opt: &IndexOptions) -> Result<(), String> {
    if !(1..=prescreen::K_MAX).contains(&opt.k) {
        return Err(format!(
            "--k 必须在 1..={} 之间（2-bit 编码上限），得到 {}",
            prescreen::K_MAX,
            opt.k
        ));
    }
    if opt.threads == 0 {
        return Err("--threads 必须大于 0".into());
    }
    if opt.host_fa.as_os_str().is_empty() || opt.target_fa.as_os_str().is_empty() {
        return Err("构建索引需要 --host-fa 与 --target-fa".into());
    }
    if opt.out_dir.as_os_str().is_empty() {
        return Err("构建索引需要 --out".into());
    }
    if opt.decoy_fa.is_none() {
        if opt.decoy_anis.is_empty() {
            return Err("--decoy-ani 至少一层".into());
        }
        for &a in &opt.decoy_anis {
            if a == 0 || a >= 100 {
                return Err(format!("--decoy-ani 层 {a} 非法（须在 1..99 之间）"));
            }
        }
        if opt.decoy_per_layer == 0 {
            return Err("--decoy-per-layer 必须大于 0".into());
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Bloom 二进制序列化（私有格式，带魔数与版本；跨平台确定性）
// ---------------------------------------------------------------------------

const BLOOM_MAGIC: [u8; 4] = *b"VFB1";
const BLOOM_VERSION: u32 = 1;

fn write_bloom(path: &Path, bloom: &KmerBloom) -> Result<(), String> {
    let (words, _mask) = bloom.to_raw();
    let mut out = Vec::with_capacity(40 + words.len() * 8);
    out.extend_from_slice(&BLOOM_MAGIC);
    out.extend_from_slice(&BLOOM_VERSION.to_le_bytes());
    out.extend_from_slice(&(bloom.k as u64).to_le_bytes());
    out.extend_from_slice(&bloom.n_inserted.to_le_bytes());
    out.extend_from_slice(&(words.len() as u64).to_le_bytes());
    for w in words {
        out.extend_from_slice(&w.to_le_bytes());
    }
    std::fs::write(path, out).map_err(|e| format!("写 {} 失败: {e}", path.display()))
}

fn read_bloom(path: &Path) -> Result<KmerBloom, String> {
    let bytes =
        std::fs::read(path).map_err(|e| format!("索引缺少或无法读取 {}: {e}", path.display()))?;
    let header = 4 + 4 + 8 + 8 + 8;
    if bytes.len() < header {
        return Err(format!("{} 损坏（长度不足）", path.display()));
    }
    if bytes[..4] != BLOOM_MAGIC {
        return Err(format!(
            "{} 魔数不匹配（版本不兼容或文件损坏）",
            path.display()
        ));
    }
    let version = u32::from_le_bytes(bytes[4..8].try_into().unwrap());
    if version > BLOOM_VERSION {
        return Err(format!(
            "{} 版本 {version} 高于当前支持的 {BLOOM_VERSION}，请重建索引",
            path.display()
        ));
    }
    let k = u64::from_le_bytes(bytes[8..16].try_into().unwrap()) as usize;
    if !(1..=prescreen::K_MAX).contains(&k) {
        return Err(format!("{} 损坏（k 非法: {k}）", path.display()));
    }
    let n_inserted = u64::from_le_bytes(bytes[16..24].try_into().unwrap());
    let words_len = u64::from_le_bytes(bytes[24..32].try_into().unwrap()) as usize;
    if words_len == 0 || !words_len.is_power_of_two() {
        return Err(format!(
            "{} 损坏（位容量非 2 的幂: {words_len}）",
            path.display()
        ));
    }
    if bytes.len() != header + words_len * 8 {
        return Err(format!(
            "{} 损坏（长度不一致: 期望 {}，实际 {}）",
            path.display(),
            header + words_len * 8,
            bytes.len()
        ));
    }
    let mut words = Vec::with_capacity(words_len);
    for i in 0..words_len {
        let off = header + i * 8;
        words.push(u64::from_le_bytes(bytes[off..off + 8].try_into().unwrap()));
    }
    Ok(KmerBloom::from_raw(
        words,
        (words_len * 64 - 1) as u64,
        k,
        n_inserted,
    ))
}

// ---------------------------------------------------------------------------
// manifest JSON 写/读（手工序列化，零依赖）
// ---------------------------------------------------------------------------

fn write_manifest(m: &IndexManifest) -> String {
    let mut s = String::new();
    s.push_str("{\n");
    s.push_str(&format!(
        "  \"format\": \"viroflash.index\",\n  \"version\": {},\n  \"k\": {},\n  \"created_at_unix\": {},\n",
        m.format_version, m.k, m.created_at_unix
    ));
    s.push_str("  \"refs\": {\n");
    s.push_str(&format!(
        "    \"host\": {},\n    \"target\": {},\n",
        file_info_json(&m.refs.host),
        file_info_json(&m.refs.target)
    ));
    match &m.refs.contam {
        Some(c) => s.push_str(&format!("    \"contam\": {},\n", file_info_json(c))),
        None => s.push_str("    \"contam\": null,\n"),
    }
    match &m.refs.decoy {
        DecoySource::File { path, blake3 } => s.push_str(&format!(
            "    \"decoy\": {{\"file\": {{\"path\": \"{}\", \"blake3\": \"{}\"}}}}\n",
            json_escape(path),
            blake3
        )),
        DecoySource::Generated {
            anis,
            per_layer,
            seed,
        } => {
            let anis_s: Vec<String> = anis.iter().map(|a| a.to_string()).collect();
            s.push_str(&format!(
                "    \"decoy\": {{\"generated\": {{\"anis\": [{}], \"per_layer\": {per_layer}, \"seed\": {seed}}}}}\n",
                anis_s.join(", ")
            ));
        }
    }
    s.push_str("  },\n");
    s.push_str("  \"contigs\": [\n");
    for (i, c) in m.contigs.iter().enumerate() {
        // gc 用 f64 最短往返表示（Display），解析后与构建时逐位一致，
        // 避免 6 位截断在分层边界（0.40/0.50/0.60）上翻转 stratum。
        s.push_str(&format!(
            "    {{\"name\": \"{}\", \"role\": \"{}\", \"len\": {}, \"gc\": {}}}{}\n",
            json_escape(&c.name),
            c.role.prefix(),
            c.len,
            c.gc_frac,
            if i + 1 < m.contigs.len() { "," } else { "" }
        ));
    }
    s.push_str("  ],\n");
    s.push_str(&format!(
        "  \"files\": {{\"mmi\": \"{}\", \"bloom\": \"{}\"}},\n",
        MMI_NAME, BLOOM_NAME
    ));
    s.push_str(&format!(
        "  \"bloom\": {{\"n_inserted\": {}, \"fill_frac\": {:.6}}}\n",
        m.bloom.n_inserted, m.bloom.fill_frac
    ));
    s.push_str("}\n");
    s
}

fn file_info_json(f: &FileInfo) -> String {
    format!(
        "{{\"path\": \"{}\", \"blake3\": \"{}\"}}",
        json_escape(&f.path),
        f.blake3
    )
}

// --- 最小 JSON 解析器（仅本 manifest 所需子集） ---

#[derive(Debug, PartialEq)]
enum Json {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

struct JsonParser<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> JsonParser<'a> {
    fn new(text: &'a str) -> Self {
        Self {
            bytes: text.as_bytes(),
            pos: 0,
        }
    }

    fn err(&self, msg: &str) -> String {
        format!("manifest.json 解析失败（偏移 {}）: {msg}", self.pos)
    }

    fn ws(&mut self) {
        while self.pos < self.bytes.len()
            && matches!(self.bytes[self.pos], b' ' | b'\t' | b'\n' | b'\r')
        {
            self.pos += 1;
        }
    }

    fn peek(&mut self) -> Result<u8, String> {
        self.ws();
        self.bytes
            .get(self.pos)
            .copied()
            .ok_or_else(|| self.err("意外结束"))
    }

    fn expect(&mut self, b: u8) -> Result<(), String> {
        if self.peek()? == b {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.err(&format!("期望 '{}'", b as char)))
        }
    }

    fn parse(&mut self) -> Result<Json, String> {
        let v = self.parse_value()?;
        self.ws();
        if self.pos != self.bytes.len() {
            return Err(self.err("多余内容"));
        }
        Ok(v)
    }

    fn parse_value(&mut self) -> Result<Json, String> {
        match self.peek()? {
            b'{' => self.parse_object(),
            b'[' => self.parse_array(),
            b'"' => Ok(Json::Str(self.parse_string()?)),
            b't' => {
                self.expect_lit(b"true")?;
                Ok(Json::Bool(true))
            }
            b'f' => {
                self.expect_lit(b"false")?;
                Ok(Json::Bool(false))
            }
            b'n' => {
                self.expect_lit(b"null")?;
                Ok(Json::Null)
            }
            b'-' | b'0'..=b'9' => self.parse_number(),
            other => Err(self.err(&format!("非法字符 '{}'", other as char))),
        }
    }

    fn expect_lit(&mut self, lit: &[u8]) -> Result<(), String> {
        for &b in lit {
            if self.peek()? == b {
                self.pos += 1;
            } else {
                return Err(self.err("非法字面量"));
            }
        }
        Ok(())
    }

    fn parse_object(&mut self) -> Result<Json, String> {
        self.expect(b'{')?;
        let mut pairs = Vec::new();
        if self.peek()? == b'}' {
            self.pos += 1;
            return Ok(Json::Obj(pairs));
        }
        loop {
            let key = self.parse_string()?;
            self.expect(b':')?;
            let value = self.parse_value()?;
            pairs.push((key, value));
            match self.peek()? {
                b',' => {
                    self.pos += 1;
                }
                b'}' => {
                    self.pos += 1;
                    break;
                }
                other => {
                    return Err(self.err(&format!("期望 ',' 或 '}}'，得到 '{}'", other as char)))
                }
            }
        }
        Ok(Json::Obj(pairs))
    }

    fn parse_array(&mut self) -> Result<Json, String> {
        self.expect(b'[')?;
        let mut items = Vec::new();
        if self.peek()? == b']' {
            self.pos += 1;
            return Ok(Json::Arr(items));
        }
        loop {
            items.push(self.parse_value()?);
            match self.peek()? {
                b',' => {
                    self.pos += 1;
                }
                b']' => {
                    self.pos += 1;
                    break;
                }
                other => {
                    return Err(self.err(&format!("期望 ',' 或 ']'，得到 '{}'", other as char)))
                }
            }
        }
        Ok(Json::Arr(items))
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut out = String::new();
        loop {
            let b = *self
                .bytes
                .get(self.pos)
                .ok_or_else(|| self.err("字符串未闭合"))?;
            self.pos += 1;
            match b {
                b'"' => break,
                b'\\' => {
                    let e = *self
                        .bytes
                        .get(self.pos)
                        .ok_or_else(|| self.err("转义未闭合"))?;
                    self.pos += 1;
                    match e {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'b' => out.push('\u{0008}'),
                        b'f' => out.push('\u{000c}'),
                        b'n' => out.push('\n'),
                        b'r' => out.push('\r'),
                        b't' => out.push('\t'),
                        b'u' => out.push(self.parse_unicode_escape()?),
                        other => return Err(self.err(&format!("非法转义 '\\{}'", other as char))),
                    }
                }
                0x20..=0x7e | 0x80..=0xff => out.push(b as char),
                other => return Err(self.err(&format!("字符串中非法控制字符 0x{other:02x}"))),
            }
        }
        Ok(out)
    }

    /// \uXXXX，含代理对组合。
    fn parse_unicode_escape(&mut self) -> Result<char, String> {
        let hi = self.parse_hex4()?;
        let cp = if (0xD800..=0xDBFF).contains(&hi) {
            // 高代理：必须跟随 \uDC00-\uDFFF。
            if self.bytes.get(self.pos) == Some(&b'\\')
                && self.bytes.get(self.pos + 1) == Some(&b'u')
            {
                self.pos += 2;
                let lo = self.parse_hex4()?;
                if !(0xDC00..=0xDFFF).contains(&lo) {
                    return Err(self.err("非法低代理"));
                }
                0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
            } else {
                return Err(self.err("高代理未配低代理"));
            }
        } else {
            hi
        };
        char::from_u32(cp).ok_or_else(|| self.err("非法 Unicode 码点"))
    }

    fn parse_hex4(&mut self) -> Result<u32, String> {
        let mut v = 0u32;
        for _ in 0..4 {
            let b = *self
                .bytes
                .get(self.pos)
                .ok_or_else(|| self.err("\\u 转义截断"))?;
            self.pos += 1;
            v = v * 16
                + match b {
                    b'0'..=b'9' => u32::from(b - b'0'),
                    b'a'..=b'f' => u32::from(b - b'a' + 10),
                    b'A'..=b'F' => u32::from(b - b'A' + 10),
                    _ => return Err(self.err("\\u 转义中非法十六进制")),
                };
        }
        Ok(v)
    }

    fn parse_number(&mut self) -> Result<Json, String> {
        let start = self.pos;
        while self.pos < self.bytes.len()
            && matches!(
                self.bytes[self.pos],
                b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E'
            )
        {
            self.pos += 1;
        }
        let text =
            std::str::from_utf8(&self.bytes[start..self.pos]).map_err(|_| self.err("非法数字"))?;
        text.parse::<f64>()
            .map(Json::Num)
            .map_err(|_| self.err(&format!("非法数字: {text}")))
    }
}

// --- 从解析结果提取 manifest 字段 ---

fn obj_get<'a>(obj: &'a [(String, Json)], key: &str) -> Option<&'a Json> {
    obj.iter().find(|(k, _)| k == key).map(|(_, v)| v)
}

fn as_obj<'a>(v: &'a Json, what: &str) -> Result<&'a [(String, Json)], String> {
    match v {
        Json::Obj(pairs) => Ok(pairs),
        _ => Err(format!("manifest.json 字段 {what} 应为对象")),
    }
}

fn as_str<'a>(v: &'a Json, what: &str) -> Result<&'a str, String> {
    match v {
        Json::Str(s) => Ok(s),
        _ => Err(format!("manifest.json 字段 {what} 应为字符串")),
    }
}

fn as_u64(v: &Json, what: &str) -> Result<u64, String> {
    match v {
        Json::Num(n) if *n >= 0.0 && n.fract() == 0.0 && *n <= u64::MAX as f64 => Ok(*n as u64),
        _ => Err(format!("manifest.json 字段 {what} 应为非负整数")),
    }
}

fn as_f64(v: &Json, what: &str) -> Result<f64, String> {
    match v {
        Json::Num(n) => Ok(*n),
        _ => Err(format!("manifest.json 字段 {what} 应为数字")),
    }
}

fn role_from_str(s: &str) -> Result<Role, String> {
    match s {
        "host" => Ok(Role::Host),
        "target" => Ok(Role::Target),
        "decoy" => Ok(Role::Decoy),
        "contam" => Ok(Role::Contaminant),
        other => Err(format!("manifest.json 中未知角色: {other}")),
    }
}

fn parse_manifest(text: &str) -> Result<IndexManifest, String> {
    let root = JsonParser::new(text).parse()?;
    let top = as_obj(&root, "根")?;
    let format = as_str(
        obj_get(top, "format").ok_or("manifest.json 缺少 format 字段")?,
        "format",
    )?;
    if format != "viroflash.index" {
        return Err(format!(
            "manifest.json format 应为 \"viroflash.index\"，得到 \"{format}\""
        ));
    }
    let format_version = as_u64(
        obj_get(top, "version").ok_or("manifest.json 缺少 version 字段")?,
        "version",
    )? as u32;
    let k = as_u64(obj_get(top, "k").ok_or("manifest.json 缺少 k 字段")?, "k")? as usize;
    if !(1..=prescreen::K_MAX).contains(&k) {
        return Err(format!("manifest.json 中 k 非法: {k}"));
    }
    let created_at_unix = as_u64(
        obj_get(top, "created_at_unix").ok_or("manifest.json 缺少 created_at_unix 字段")?,
        "created_at_unix",
    )?;
    let refs_root = as_obj(
        obj_get(top, "refs").ok_or("manifest.json 缺少 refs 字段")?,
        "refs",
    )?;
    let parse_file = |key: &str| -> Result<FileInfo, String> {
        let o = as_obj(
            obj_get(refs_root, key).ok_or_else(|| format!("manifest.json refs 缺少 {key}"))?,
            key,
        )?;
        Ok(FileInfo {
            path: as_str(
                obj_get(o, "path").ok_or_else(|| format!("refs.{key} 缺少 path"))?,
                "path",
            )?
            .to_string(),
            blake3: as_str(
                obj_get(o, "blake3").ok_or_else(|| format!("refs.{key} 缺少 blake3"))?,
                "blake3",
            )?
            .to_string(),
        })
    };
    let contam = match obj_get(refs_root, "contam") {
        None | Some(Json::Null) => None,
        Some(_) => Some(parse_file("contam")?),
    };
    let decoy = as_obj(
        obj_get(refs_root, "decoy").ok_or("manifest.json refs 缺少 decoy")?,
        "decoy",
    )?;
    let decoy_source = if let Some(v) = obj_get(decoy, "file") {
        let o = as_obj(v, "decoy.file")?;
        DecoySource::File {
            path: as_str(obj_get(o, "path").ok_or("decoy.file 缺少 path")?, "path")?.to_string(),
            blake3: as_str(
                obj_get(o, "blake3").ok_or("decoy.file 缺少 blake3")?,
                "blake3",
            )?
            .to_string(),
        }
    } else if let Some(v) = obj_get(decoy, "generated") {
        let o = as_obj(v, "decoy.generated")?;
        let anis = match obj_get(o, "anis").ok_or("decoy.generated 缺少 anis")? {
            Json::Arr(items) => items
                .iter()
                .map(|i| as_u64(i, "anis").map(|n| n as u8))
                .collect::<Result<Vec<u8>, _>>()?,
            _ => return Err("decoy.generated.anis 应为数组".into()),
        };
        DecoySource::Generated {
            anis,
            per_layer: as_u64(
                obj_get(o, "per_layer").ok_or("decoy.generated 缺少 per_layer")?,
                "per_layer",
            )? as usize,
            seed: as_u64(
                obj_get(o, "seed").ok_or("decoy.generated 缺少 seed")?,
                "seed",
            )?,
        }
    } else {
        return Err("manifest.json refs.decoy 应为 file 或 generated".into());
    };
    let contigs = match obj_get(top, "contigs").ok_or("manifest.json 缺少 contigs 字段")? {
        Json::Arr(items) => items
            .iter()
            .map(|v| {
                let o = as_obj(v, "contigs[]")?;
                Ok(ContigMeta {
                    name: as_str(obj_get(o, "name").ok_or("contig 缺少 name")?, "name")?
                        .to_string(),
                    role: role_from_str(as_str(
                        obj_get(o, "role").ok_or("contig 缺少 role")?,
                        "role",
                    )?)?,
                    len: as_u64(obj_get(o, "len").ok_or("contig 缺少 len")?, "len")?,
                    gc_frac: as_f64(obj_get(o, "gc").ok_or("contig 缺少 gc")?, "gc")?,
                })
            })
            .collect::<Result<Vec<ContigMeta>, String>>()?,
        _ => return Err("manifest.json contigs 应为数组".into()),
    };
    let bloom_root = as_obj(
        obj_get(top, "bloom").ok_or("manifest.json 缺少 bloom 字段")?,
        "bloom",
    )?;
    let bloom = BloomInfo {
        n_inserted: as_u64(
            obj_get(bloom_root, "n_inserted").ok_or("bloom 缺少 n_inserted")?,
            "n_inserted",
        )?,
        fill_frac: as_f64(
            obj_get(bloom_root, "fill_frac").ok_or("bloom 缺少 fill_frac")?,
            "fill_frac",
        )?,
    };
    Ok(IndexManifest {
        format_version,
        k,
        created_at_unix,
        refs: RefsInfo {
            host: parse_file("host")?,
            target: parse_file("target")?,
            contam,
            decoy: decoy_source,
        },
        contigs,
        bloom,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tmp_dir(tag: &str) -> PathBuf {
        let d =
            std::env::temp_dir().join(format!("viroflash_index_test_{tag}_{}", std::process::id()));
        if d.exists() {
            let _ = std::fs::remove_dir_all(&d);
        }
        std::fs::create_dir_all(&d).unwrap();
        d
    }

    fn write_fa(path: &Path, name: &str, seq: &[u8]) {
        let mut out = String::from(">");
        out.push_str(name);
        out.push('\n');
        out.push_str(std::str::from_utf8(seq).unwrap());
        out.push('\n');
        std::fs::write(path, out).unwrap();
    }

    fn default_index_opts(dir: &Path) -> IndexOptions {
        let refs = dir.join("refs");
        std::fs::create_dir_all(&refs).unwrap();
        write_fa(
            &refs.join("host.fa"),
            "chrH",
            b"ACGTACGTACGTACGTACGTACGTACGT",
        );
        write_fa(
            &refs.join("target.fa"),
            "tv",
            b"GATTACAGATTACAGATTACAGATTACA",
        );
        write_fa(
            &refs.join("decoy.fa"),
            "d0",
            b"TGGCTAGCTTGGCTAGCTTGGCTAGCTT",
        );
        write_fa(
            &refs.join("contam.fa"),
            "myco",
            b"CCTAGGCCTAGGCCTAGGCCTAGGCCTA",
        );
        IndexOptions {
            host_fa: refs.join("host.fa"),
            target_fa: refs.join("target.fa"),
            contam_fa: Some(refs.join("contam.fa")),
            decoy_fa: Some(refs.join("decoy.fa")),
            out_dir: dir.join("idx"),
            ..IndexOptions::default()
        }
    }

    #[test]
    fn build_and_load_roundtrip() {
        let dir = tmp_dir("roundtrip");
        let opt = default_index_opts(&dir);
        let built = build_index(&opt).unwrap();
        assert_eq!(built.contigs.len(), 4);
        assert!(built.mmi_path.is_file());
        assert!(opt.out_dir.join(MANIFEST_NAME).is_file());
        assert!(opt.out_dir.join(BLOOM_NAME).is_file());

        let loaded = load_index(&opt.out_dir, opt.k).unwrap();
        assert_eq!(loaded.contigs, built.contigs);
        assert_eq!(loaded.roles, built.roles);
        assert_eq!(loaded.manifest_blake3, built.manifest_blake3);
        assert_eq!(loaded.bloom.k, built.bloom.k);
        assert_eq!(loaded.bloom.n_inserted, built.bloom.n_inserted);
        // 角色齐全：host/target/decoy/contam 各 1。
        let counts = |c: &[ContigMeta], r: Role| c.iter().filter(|x| x.role == r).count();
        assert_eq!(counts(&loaded.contigs, Role::Host), 1);
        assert_eq!(counts(&loaded.contigs, Role::Target), 1);
        assert_eq!(counts(&loaded.contigs, Role::Decoy), 1);
        assert_eq!(counts(&loaded.contigs, Role::Contaminant), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn manifest_roundtrip_preserves_fields() {
        let m = IndexManifest {
            format_version: FORMAT_VERSION,
            k: 21,
            created_at_unix: 123456,
            refs: RefsInfo {
                host: FileInfo {
                    path: "/a/host.fa".into(),
                    blake3: "ab".repeat(32),
                },
                target: FileInfo {
                    path: "/a/target.fa".into(),
                    blake3: "cd".repeat(32),
                },
                contam: None,
                decoy: DecoySource::Generated {
                    anis: vec![82, 85, 88],
                    per_layer: 4,
                    seed: 7,
                },
            },
            contigs: vec![
                ContigMeta {
                    name: "host_0".into(),
                    role: Role::Host,
                    len: 1_000_000,
                    gc_frac: 0.40963,
                },
                ContigMeta {
                    name: "decoy:tv:ani82:i0".into(),
                    role: Role::Decoy,
                    len: 293116,
                    gc_frac: 0.52,
                },
            ],
            bloom: BloomInfo {
                n_inserted: 42,
                fill_frac: 0.4123,
            },
        };
        let parsed = parse_manifest(&write_manifest(&m)).unwrap();
        assert_eq!(parsed, m);
        // 字符串转义往返：路径含引号/反斜杠/控制字符。
        let mut m2 = m.clone();
        m2.refs.host.path = "a\"b\\c\nd".into();
        let parsed2 = parse_manifest(&write_manifest(&m2)).unwrap();
        assert_eq!(parsed2, m2);
    }

    #[test]
    fn load_rejects_k_mismatch_and_version() {
        let dir = tmp_dir("mismatch");
        let opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        let err = load_index(&opt.out_dir, opt.k + 1).unwrap_err();
        assert!(err.contains("不一致"), "err={err}");
        // 版本高于当前支持：手工改写 version 字段后必须拒绝。
        let mp = opt.out_dir.join(MANIFEST_NAME);
        let text =
            std::fs::read_to_string(&mp)
                .unwrap()
                .replacen("\"version\": 1", "\"version\": 99", 1);
        std::fs::write(&mp, text).unwrap();
        let err = load_index(&opt.out_dir, opt.k).unwrap_err();
        assert!(err.contains("高于"), "err={err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_rejects_missing_and_corrupt() {
        let dir = tmp_dir("corrupt");
        let opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        // 缺 manifest（空目录）
        let empty = dir.join("empty");
        std::fs::create_dir_all(&empty).unwrap();
        let err = load_index(&empty, opt.k).unwrap_err();
        assert!(err.contains("manifest.json"), "err={err}");
        // bloom 截断
        let bp = opt.out_dir.join(BLOOM_NAME);
        let bytes = std::fs::read(&bp).unwrap();
        std::fs::write(&bp, &bytes[..bytes.len() / 2]).unwrap();
        let err = load_index(&opt.out_dir, opt.k).unwrap_err();
        assert!(err.contains("损坏"), "err={err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn build_rejects_existing_dir_and_bad_params() {
        let dir = tmp_dir("exist");
        let mut opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        let err = build_index(&opt).unwrap_err();
        assert!(err.contains("已存在"), "err={err}");
        // 构建失败须清理临时目录。
        opt.out_dir = dir.join("idx2");
        opt.k = 0;
        let err = build_index(&opt).unwrap_err();
        assert!(err.contains("--k"), "err={err}");
        let leftovers: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().starts_with("idx2"))
            .collect();
        assert!(leftovers.is_empty(), "失败后残留临时目录");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn auto_decoy_generation_writes_into_index() {
        let dir = tmp_dir("autodecoy");
        let mut opt = default_index_opts(&dir);
        opt.decoy_fa = None;
        opt.decoy_anis = vec![82, 88];
        opt.decoy_per_layer = 2;
        opt.decoy_seed = 3;
        let built = build_index(&opt).unwrap();
        // host + target + 2 层 × 2 = 6 个 contig（无 contam 时省略污染）。
        let mut opt_nc = opt.clone();
        opt_nc.out_dir = dir.join("idx_nc");
        opt_nc.contam_fa = None;
        let built_nc = build_index(&opt_nc).unwrap();
        assert_eq!(built_nc.contigs.len(), 2 + 2 * 2);
        assert!(opt.out_dir.join(DECOYS_FA_NAME).is_file());
        assert!(opt.out_dir.join(DECOYS_TSV_NAME).is_file());
        let manifest =
            parse_manifest(&std::fs::read_to_string(opt.out_dir.join(MANIFEST_NAME)).unwrap())
                .unwrap();
        assert_eq!(
            manifest.refs.decoy,
            DecoySource::Generated {
                anis: vec![82, 88],
                per_layer: 2,
                seed: 3
            }
        );
        assert!(
            built
                .contigs
                .iter()
                .filter(|c| c.role == Role::Decoy)
                .count()
                == 4
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn bloom_binary_roundtrip() {
        let dir = tmp_dir("bloom");
        let opt = default_index_opts(&dir);
        let built = build_index(&opt).unwrap();
        let reloaded = read_bloom(&opt.out_dir.join(BLOOM_NAME)).unwrap();
        assert_eq!(reloaded.k, built.bloom.k);
        assert_eq!(reloaded.n_inserted, built.bloom.n_inserted);
        assert_eq!(reloaded.fill_frac(), built.bloom.fill_frac());
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn parser_rejects_malformed_json() {
        assert!(parse_manifest("").is_err());
        assert!(parse_manifest("{").is_err());
        assert!(parse_manifest("{}").unwrap_err().contains("format"));
        assert!(parse_manifest("{\"format\":\"viroflash.index\"}").is_err());
        assert!(parse_manifest("{\"format\":\"other\",\"version\":1,\"k\":21}").is_err());
    }
}
