//! 流式 FASTQ(.gz) 读取：4 行记录、配对 ID 校验、gzip 流式解压。
//! 不把整份数据载入内存。两条读取路径语义等价并逐记录输出一致：
//! - 单线程：`read_until` 字节级分行，复用行缓冲；1MB 读缓冲让 zlib-rs
//!   以较大块解压，减少小块调用开销。
//! - 成员级并行：解压线程按 gzip 成员并行解压并用 memchr 切行，主线程
//!   从行队列每 4 行构建一条记录；跨成员的行由队列自然拼接。

use std::collections::{BTreeMap, VecDeque};
use std::fs::File;
use std::io::{BufRead, BufReader, Read, Seek, SeekFrom};
use std::path::Path;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{sync_channel, Receiver};
use std::sync::Arc;
use std::thread::JoinHandle;

use flate2::read::{GzDecoder, MultiGzDecoder};
use memchr::memchr_iter;

#[derive(Debug, Clone)]
pub struct FastqRecord {
    pub id: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

/// 并行泵的零拷贝记录：成员内部完整 4 行组共享解压缓冲；
/// 跨成员边界 / EOF 尾记录（每成员至多 1 条）物化为 Owned。
/// 热路径（预筛）只消费 seq/id 切片，不产生每记录堆分配。
#[derive(Debug)]
pub enum SpanRecord {
    /// 共享缓冲上的字节区间（id 已去 @ 前缀与首尾 ASCII 空白）。
    Span(Arc<Vec<u8>>, (u32, u32), (u32, u32), (u32, u32)),
    Owned(FastqRecord),
}

impl SpanRecord {
    pub fn seq(&self) -> &[u8] {
        match self {
            SpanRecord::Span(buf, _, s, _) => &buf[s.0 as usize..s.1 as usize],
            SpanRecord::Owned(r) => &r.seq,
        }
    }
    pub fn id(&self) -> &[u8] {
        match self {
            SpanRecord::Span(buf, i, _, _) => &buf[i.0 as usize..i.1 as usize],
            SpanRecord::Owned(r) => r.id.as_bytes(),
        }
    }
    pub fn has_seq(&self) -> bool {
        !self.seq().is_empty()
    }
    pub fn qual(&self) -> &[u8] {
        match self {
            SpanRecord::Span(buf, _, _, q) => &buf[q.0 as usize..q.1 as usize],
            SpanRecord::Owned(r) => &r.qual,
        }
    }
}

/// 并行泵产出的配对（配对 ID 已按 normalize_pair_id 语义字节级校验一致）。
#[derive(Debug)]
pub struct SpanPair {
    pub r1: SpanRecord,
    pub r2: SpanRecord,
}

/// 单线程读取源：BufRead + 复用的行缓冲（现行为不变）。
struct StreamFeed {
    reader: Box<dyn BufRead>,
    line: Vec<u8>,
}

/// 读取源：单线程字节流 或 成员级 并行记录泵（见模块文档）。
enum Feed {
    Stream(StreamFeed),
    Pump(LinePump),
}

pub struct FastqReader {
    feed: Feed,
}

impl FastqReader {
    pub fn open(path: &Path) -> Result<Self, String> {
        let file = File::open(path).map_err(|e| format!("无法打开 {}: {e}", path.display()))?;
        let inner: Box<dyn Read> = if path.to_string_lossy().ends_with(".gz") {
            Box::new(MultiGzDecoder::new(file))
        } else {
            Box::new(file)
        };
        Ok(Self {
            feed: Feed::Stream(StreamFeed {
                reader: Box::new(BufReader::with_capacity(1 << 20, inner)),
                line: Vec::with_capacity(1024),
            }),
        })
    }

    /// 成员级并行解压打开（成员级）：扫描 gzip 成员索引（IG extra 子域），`d` 个
    /// 解压线程按成员并行解压+切行、按序供给行队列。无 IG 索引 / 非 gzip /
    /// d==0 时回退单线程 `open`（吞吐不降级为错）。输出与单线程逐记录一致。
    pub fn open_parallel(path: &Path, d: usize) -> Result<Self, String> {
        if d == 0 || !path.to_string_lossy().ends_with(".gz") {
            return Self::open(path);
        }
        let members = build_member_index(path)?;
        if members.len() < 2 {
            // 单成员/无索引：并行无收益，走单线程路径
            return Self::open(path);
        }
        let (rx, handles) = start_member_pumps(path, members.clone(), d);
        Ok(Self {
            feed: Feed::Pump(LinePump::new(rx, handles, members.len())),
        })
    }

    pub fn next_record(&mut self) -> Result<Option<FastqRecord>, String> {
        match &mut self.feed {
            Feed::Stream(f) => next_record_stream(f),
            Feed::Pump(p) => next_record_pump(p),
        }
    }

    /// 零拷贝变体：泵路径四行同缓冲时返回区间（无每记录分配），其余回退
    /// 物化。与 `next_record` 逐记录等价（同一校验顺序与文案）。
    pub fn next_record_span(&mut self) -> Result<Option<SpanRecord>, String> {
        match &mut self.feed {
            Feed::Stream(f) => next_record_stream(f).map(|r| r.map(SpanRecord::Owned)),
            Feed::Pump(p) => next_record_pump_span(p),
        }
    }
}

/// 单线程路径：读一行并剥去行尾 \n / \r\n，返回对内部缓冲的借用（下次调用前失效）。
fn read_line_raw(feed: &mut StreamFeed) -> Result<Option<&[u8]>, String> {
    feed.line.clear();
    let n = feed
        .reader
        .read_until(b'\n', &mut feed.line)
        .map_err(|e| format!("读取 FASTQ 失败: {e}"))?;
    if n == 0 {
        return Ok(None);
    }
    let mut end = feed.line.len();
    if end > 0 && feed.line[end - 1] == b'\n' {
        end -= 1;
    }
    if end > 0 && feed.line[end - 1] == b'\r' {
        end -= 1;
    }
    Ok(Some(&feed.line[..end]))
}

/// 单线程路径的逐记录解析。
fn next_record_stream(feed: &mut StreamFeed) -> Result<Option<FastqRecord>, String> {
    let name_line = match read_line_raw(feed)? {
        None => return Ok(None),
        Some(l) => l.to_vec(),
    };
    let seq = read_line_raw(feed)?
        .ok_or_else(|| "FASTQ 记录截断（缺序列行）".to_string())?
        .to_vec();
    let plus = read_line_raw(feed)?
        .ok_or_else(|| "FASTQ 记录截断（缺 + 行）".to_string())?
        .to_vec();
    let qual = read_line_raw(feed)?
        .ok_or_else(|| "FASTQ 记录截断（缺质量行）".to_string())?
        .to_vec();
    assemble_record(name_line, seq, plus, qual).map(Some)
}

/// 成员级 记录泵路径：行已由解压线程切好（Vec 所有权转移，无拷贝），
/// 检查顺序与错误文案和单线程路径完全一致。
fn next_record_pump(pump: &mut LinePump) -> Result<Option<FastqRecord>, String> {
    let name_line = match pump.pop_line()? {
        None => return Ok(None),
        Some(l) => l.as_bytes().to_vec(),
    };
    let seq = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺序列行）".to_string())?
        .as_bytes()
        .to_vec();
    let plus = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺 + 行）".to_string())?
        .as_bytes()
        .to_vec();
    let qual = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺质量行）".to_string())?
        .as_bytes()
        .to_vec();
    assemble_record(name_line, seq, plus, qual).map(Some)
}

/// 记录泵的零拷贝路径：四行同缓冲时返回共享缓冲上的区间（id 已剥 @ 与首尾
/// ASCII 空白），跨成员边界 / 非 ASCII ID 时物化回退。校验顺序与错误文案和
/// `assemble_record`（唯一事实来源）一致；非 ASCII ID 走物化路径，UTF-8 语义
/// 与单线程路径逐字节等价。
fn next_record_pump_span(pump: &mut LinePump) -> Result<Option<SpanRecord>, String> {
    let name = match pump.pop_line()? {
        None => return Ok(None),
        Some(l) => l,
    };
    let seq = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺序列行）".to_string())?;
    let plus = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺 + 行）".to_string())?;
    let qual = pump
        .pop_line()?
        .ok_or_else(|| "FASTQ 记录截断（缺质量行）".to_string())?;

    if plus.as_bytes().first() != Some(&b'+') {
        return Err(format!(
            "FASTQ + 行格式错误: {:?}",
            String::from_utf8_lossy(plus.as_bytes())
        ));
    }
    let name_bytes = name.as_bytes();
    if name_bytes.first() != Some(&b'@') {
        return Err(format!(
            "FASTQ 记录格式错误: {:?}",
            String::from_utf8_lossy(name_bytes)
        ));
    }
    // 非 ASCII ID：物化后由 assemble_record 执行 UTF-8 校验。
    if !name_bytes.is_ascii() {
        return assemble_record(
            name_bytes.to_vec(),
            seq.as_bytes().to_vec(),
            plus.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .map(SpanRecord::Owned)
        .map(Some);
    }
    // 剥全部前导 @（trim_start_matches('@') 等价）与首尾 ASCII 空白
    // （ASCII 输入上 trim() 等价）。
    let mut lo = 0usize;
    while lo < name_bytes.len() && name_bytes[lo] == b'@' {
        lo += 1;
    }
    while lo < name_bytes.len() && name_bytes[lo].is_ascii_whitespace() {
        lo += 1;
    }
    let mut hi = name_bytes.len();
    while hi > lo && name_bytes[hi - 1].is_ascii_whitespace() {
        hi -= 1;
    }
    if seq.as_bytes().len() != qual.as_bytes().len() {
        return Err(format!(
            "FASTQ 序列/质量长度不一致: {}",
            String::from_utf8_lossy(&name_bytes[lo..hi])
        ));
    }
    let same_buf = seq.buf.as_ptr() == name.buf.as_ptr()
        && plus.buf.as_ptr() == name.buf.as_ptr()
        && qual.buf.as_ptr() == name.buf.as_ptr();
    if same_buf {
        let base = name.r.0 as usize;
        Ok(Some(SpanRecord::Span(
            name.buf.clone(),
            ((base + lo) as u32, (base + hi) as u32),
            seq.r,
            qual.r,
        )))
    } else {
        assemble_record(
            name_bytes.to_vec(),
            seq.as_bytes().to_vec(),
            plus.as_bytes().to_vec(),
            qual.as_bytes().to_vec(),
        )
        .map(SpanRecord::Owned)
        .map(Some)
    }
}

/// 4 行组装 FastqRecord：两条读取路径共用的唯一事实来源。
/// 校验顺序：+ 前缀 → @ 前缀 → UTF-8 → seq/qual 等长；两条读取路径错误文案一致。
fn assemble_record(
    name_line: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    qual: Vec<u8>,
) -> Result<FastqRecord, String> {
    if plus.first() != Some(&b'+') {
        return Err(format!(
            "FASTQ + 行格式错误: {:?}",
            String::from_utf8_lossy(&plus)
        ));
    }
    if name_line.first() != Some(&b'@') {
        return Err(format!(
            "FASTQ 记录格式错误: {:?}",
            String::from_utf8_lossy(&name_line)
        ));
    }
    let id = String::from_utf8(name_line).map_err(|e| format!("FASTQ ID 非 UTF-8: {e}"))?;
    let id = id.trim_start_matches('@').trim().to_string();
    if seq.len() != qual.len() {
        return Err(format!("FASTQ 序列/质量长度不一致: {}", id));
    }
    Ok(FastqRecord { id, seq, qual })
}

/// 配对 ID 归一化：截断首个空白后的描述字段（Illumina CASAVA 1.8 `read 1:N:0:...`），
/// 再剥 /1 /2（CASAVA 1.7 后缀）——参照 viroflash ReadIdWithoutPrefixAndWhitespace。
pub fn normalize_pair_id(id: &str) -> String {
    let token = id.split_ascii_whitespace().next().unwrap_or(id);
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
        .to_string()
}

/// `normalize_pair_id` 的字节级等价实现（零分配）：跳前导 ASCII 空白取首个
/// token（split_ascii_whitespace 语义），再剥 /1 /2 后缀（CASAVA 1.7）。
/// 逐字节相等即判定配对一致。
pub fn pair_ids_eq(a: &[u8], b: &[u8]) -> bool {
    fn norm(id: &[u8]) -> &[u8] {
        let mut s = 0;
        while s < id.len() && id[s].is_ascii_whitespace() {
            s += 1;
        }
        if s == id.len() {
            // 全空白：split_ascii_whitespace 无 token，unwrap_or 原样返回
            return id;
        }
        let mut t = &id[s..];
        if let Some(pos) = t.iter().position(|c| c.is_ascii_whitespace()) {
            t = &t[..pos];
        }
        t.strip_suffix(b"/1")
            .or_else(|| t.strip_suffix(b"/2"))
            .unwrap_or(t)
    }
    norm(a) == norm(b)
}

/// 并行泵的配对迭代器（融合阶段消费）：两端共享成员级并行解压预算，
/// 逐对按 `pair_ids_eq` 字节级校验配对 ID，产零拷贝 `SpanPair`。
pub struct PairSpanIter {
    r1: FastqReader,
    r2: FastqReader,
}

impl PairSpanIter {
    /// `d_total` 为**总**解压线程预算（N=8 → 2，见 lib::thread_budget），
    /// R1/R2 各取一半（无索引自动回退单线程）。
    pub fn open_parallel(r1_path: &Path, r2_path: &Path, d_total: usize) -> Result<Self, String> {
        let d1 = d_total.div_ceil(2);
        let d2 = d_total - d1;
        Ok(Self {
            r1: FastqReader::open_parallel(r1_path, d1)?,
            r2: FastqReader::open_parallel(r2_path, d2)?,
        })
    }
}

impl Iterator for PairSpanIter {
    type Item = Result<SpanPair, String>;

    fn next(&mut self) -> Option<Self::Item> {
        let a = self.r1.next_record_span();
        let b = self.r2.next_record_span();
        match (a, b) {
            (Ok(None), Ok(None)) => None,
            (Ok(Some(_)), Ok(None)) | (Ok(None), Ok(Some(_))) => {
                Some(Err("R1/R2 记录数不一致".to_string()))
            }
            (Ok(Some(a)), Ok(Some(b))) => {
                if !pair_ids_eq(a.id(), b.id()) {
                    return Some(Err(format!(
                        "配对 ID 不一致: {} vs {}",
                        String::from_utf8_lossy(a.id()),
                        String::from_utf8_lossy(b.id())
                    )));
                }
                Some(Ok(SpanPair { r1: a, r2: b }))
            }
            (Err(e), _) | (_, Err(e)) => Some(Err(e)),
        }
    }
}

/// 单端流式迭代适配器：每个 read 一个 fragment（r2 为空 Owned，零分配）。
pub struct SingleSpanIter {
    inner: FastqReader,
}

impl SingleSpanIter {
    pub fn open_parallel(r1_path: &Path, d: usize) -> Result<Self, String> {
        Ok(Self {
            inner: FastqReader::open_parallel(r1_path, d)?,
        })
    }
}

impl Iterator for SingleSpanIter {
    type Item = Result<SpanPair, String>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.inner.next_record_span() {
            Ok(Some(r1)) => Some(Ok(SpanPair {
                r1,
                r2: SpanRecord::Owned(FastqRecord {
                    id: String::new(),
                    seq: Vec::new(),
                    qual: Vec::new(),
                }),
            })),
            Ok(None) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

// ---------------------------------------------------------------------------
// 成员级并行解压
//
// 解压线程并行解压并用 memchr 切行，主线程只做 4 行一组的轻组装。
// 多成员 gzip 的成员边界天然提供字节对齐的并行切分点（deflate 窗口在成员
// 边界重置，无需字典注入）；`IG` extra 子域携带成员总压缩长度（u32 LE），
// 可按成员链式跳转建索引。
// 参考：RFC 1952 §2.2/§2.3.1、bgzip/rapidgzip 对独立成员并行的先例。
// 未采用 pigz 式位级块扫描：deflate 块边界位对齐不可靠定位、dynamic 块需逐块
// 解 Huffman 头、LZ77 跨块 ≤32KB 需窗口传播（rapidgzip 两阶段），工程量与本
// 收益不匹配。无 IG 索引的普通 gzip 回退单线程路径。
// ---------------------------------------------------------------------------

/// 一个 gzip 成员：文件偏移 + 总压缩长（含 10B 头、extra/name/comment、数据、8B trailer）。
#[derive(Debug, Clone)]
struct Member {
    off: u64,
    clen: u64,
}

/// 头字段扫描上限（防御畸形输入；正常值远小于此）。
const MAX_XLEN: usize = 64 << 10;
const MAX_NUL_FIELD: usize = 64 << 10;
const MAX_MEMBERS: usize = 16 << 20;

/// 扫描 gzip 文件头链：FLG.FEXTRA 内找 `IG` 子域（u32 LE = 成员总压缩长），
/// 用 `next = off + ig` 链式跳转，并在每个候选落点校验 `1f 8b` magic。
/// 兼容两种 IG 编码：成员总长与仅压缩数据长；
/// 任何一步失败 → 返回空 Vec（调用方回退单线程，吞吐不降级为错）。
fn build_member_index(path: &Path) -> Result<Vec<Member>, String> {
    let mut file = File::open(path).map_err(|e| format!("无法打开 {}: {e}", path.display()))?;
    let size = file
        .metadata()
        .map_err(|e| format!("读取 {} 元数据失败: {e}", path.display()))?
        .len();
    let mut members = Vec::new();
    let mut off: u64 = 0;

    loop {
        if members.len() >= MAX_MEMBERS {
            return Ok(Vec::new()); // 成员过多（防御畸形输入）：回退单线程
        }
        // 任何头解析失败都整体回退：部分索引会让并行路径在最后一个已索引
        // 成员处提前 EOF，静默丢数据——宁回退不丢数据。
        let Some(hdr) = parse_member_header(&mut file, off)? else {
            return Ok(Vec::new());
        };
        let header_end = hdr.end;
        let Some(ig) = hdr.ig else {
            return Ok(Vec::new()); // 无 IG 子域 → 无索引，回退
        };
        if ig == 0 {
            return Ok(Vec::new());
        }
        // 两种解释分别计算落点；`== size` 视为合法末成员（magic 无从校验）。
        let next_total = off.checked_add(ig as u64);
        let next_data = header_end.checked_add(ig as u64);
        let next = match (next_total, next_data) {
            (Some(a), Some(_b)) if a == size || next_is_magic(&mut file, a)? => a,
            (_, Some(b)) if b == size || next_is_magic(&mut file, b)? => b,
            _ => return Ok(Vec::new()), // 两种解释都不成立 → 回退
        };
        members.push(Member {
            off,
            clen: next - off,
        });
        if next >= size {
            return Ok(members);
        }
        off = next;
    }
}

/// 单个成员头解析结果。
struct ParsedHeader {
    /// 压缩数据区（含 FHCRC 之后）起始偏移。
    end: u64,
    /// IG 子域值（u32 LE），无则 None。
    ig: Option<u32>,
}

/// 解析 `off` 处的 gzip 成员头（RFC 1952 §2.3）：magic/CM 校验，跳过
/// FEXTRA（其中找 `IG`）、FNAME、FCOMMENT、FHCRC。失败返回 Ok(None)（非致命，
/// 由调用方决定回退）；IO 错误返回 Err。
fn parse_member_header(file: &mut File, off: u64) -> Result<Option<ParsedHeader>, String> {
    let mut hdr = [0u8; 10];
    match file
        .seek(SeekFrom::Start(off))
        .and_then(|_| read_exact_at(file, &mut hdr))
    {
        Ok(()) => {}
        Err(_) => return Ok(None), // 头都不完整：视为无可索引结构
    }
    if hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 8 {
        return Ok(None);
    }
    let flg = hdr[3];
    let mut p = off + 10;
    let mut ig = None;

    if flg & 0x04 != 0 {
        // FEXTRA：XLEN + 子域链（SI1 SI2 LEN data），找 "IG"
        let mut xl = [0u8; 2];
        if read_exact_at(file, &mut xl).is_err() {
            return Ok(None);
        }
        let xlen = u16::from_le_bytes(xl) as usize;
        if xlen > MAX_XLEN {
            return Ok(None);
        }
        let mut extra = vec![0u8; xlen];
        if read_exact_at(file, &mut extra).is_err() {
            return Ok(None);
        }
        let mut q = 0usize;
        while q + 4 <= extra.len() {
            let sub_len = u16::from_le_bytes([extra[q + 2], extra[q + 3]]) as usize;
            let data_end = q + 4 + sub_len;
            if data_end > extra.len() {
                break;
            }
            if extra[q] == b'I' && extra[q + 1] == b'G' && sub_len == 4 {
                ig = Some(u32::from_le_bytes([
                    extra[q + 4],
                    extra[q + 5],
                    extra[q + 6],
                    extra[q + 7],
                ]));
            }
            q = data_end;
        }
        p += 2 + xlen as u64;
    }
    if flg & 0x08 != 0 {
        // FNAME：NUL 结尾
        let Some(n) = skip_to_nul(file, p)? else {
            return Ok(None);
        };
        p += n;
    }
    if flg & 0x10 != 0 {
        // FCOMMENT：NUL 结尾
        let Some(n) = skip_to_nul(file, p)? else {
            return Ok(None);
        };
        p += n;
    }
    if flg & 0x02 != 0 {
        p += 2; // FHCRC
    }
    Ok(Some(ParsedHeader { end: p, ig }))
}

fn read_exact_at(file: &mut File, buf: &mut [u8]) -> std::io::Result<()> {
    use std::io::Read;
    file.read_exact(buf)
}

/// 从 `p` 起扫描 NUL 字节，返回 NUL 之后的位置（相对 `p` 的字节数）。超限返回 None。
fn skip_to_nul(file: &mut File, p: u64) -> Result<Option<u64>, String> {
    use std::io::Read;
    file.seek(SeekFrom::Start(p))
        .map_err(|e| format!("gzip 头跳转失败: {e}"))?;
    let mut buf = [0u8; 256];
    let mut total = 0u64;
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|e| format!("gzip 头读取失败: {e}"))?;
        if n == 0 {
            return Ok(None);
        }
        if let Some(i) = buf[..n].iter().position(|&b| b == 0) {
            return Ok(Some(total + i as u64 + 1));
        }
        total += n as u64;
        if total > MAX_NUL_FIELD as u64 {
            return Ok(None);
        }
    }
}

/// `off` 处是否为下一个成员 magic（EOF 返回 false，交由调用方按 `== size` 判定末成员）。
fn next_is_magic(file: &mut File, off: u64) -> Result<bool, String> {
    let mut m = [0u8; 2];
    match file
        .seek(SeekFrom::Start(off))
        .and_then(|_| read_exact_at(file, &mut m))
    {
        Ok(()) => Ok(m == [0x1f, 0x8b]),
        Err(_) => Ok(false),
    }
}

/// 解压线程的一个成果：`lines` 为本成员内的完整行（已剥 \r\n/\n），
/// `tail` 为末尾无 \n 的部分行（可能被下一成员续上）。
struct MemberLines {
    /// 整成员解压缓冲（所有行区间共享，成员级单次分配）。
    buf: Arc<Vec<u8>>,
    /// 完整行区间（已剥 \r）；tail 为末段（无 \n 结尾，不剥 \r）。
    lines: Vec<(u32, u32)>,
    tail: Option<(u32, u32)>,
}

/// 解压线程侧切行（memchr）：完整行剥 \r\n/\n 记为区间；末尾无 \n 的段入 tail
/// （不剥 \r，跨成员拼接由主线程统一处理）。整缓冲共享、零行级分配。
fn split_lines(bytes: Vec<u8>) -> MemberLines {
    let mut lines = Vec::new();
    let mut start = 0usize;
    for pos in memchr_iter(b'\n', &bytes) {
        let mut end = pos;
        if end > start && bytes[end - 1] == b'\r' {
            end -= 1;
        }
        lines.push((start as u32, end as u32));
        start = pos + 1;
    }
    let tail = if start < bytes.len() {
        Some((start as u32, bytes.len() as u32))
    } else {
        None
    };
    MemberLines {
        buf: Arc::new(bytes),
        lines,
        tail,
    }
}

/// 解压+切行池：`d` 个线程原子取号，各自打开文件句柄 seek 到成员起点、
/// `take(clen)` 交给单成员 `GzDecoder` 解压（CRC32+ISIZE 由 GzDecoder 校验），
/// memchr 切行后按成员序回传 `(seq, Result<MemberLines>)`。有界通道（容量 2d）背压。
/// 非末成员解码失败为硬错误（与 MultiGzDecoder 语义一致）；仅截断 trailer 的
/// 退化末成员按干净 EOF 容错。
// 返回类型直接表达泵通道与线程句柄的绑定关系，拆分别名反而分散所有权。
#[allow(clippy::type_complexity)]
fn start_member_pumps(
    path: &Path,
    members: Vec<Member>,
    d: usize,
) -> (
    Receiver<(usize, Result<MemberLines, String>)>,
    Vec<JoinHandle<()>>,
) {
    let path = path.to_path_buf();
    let members = std::sync::Arc::new(members);
    let next = std::sync::Arc::new(AtomicUsize::new(0));
    let (tx, rx) = sync_channel::<(usize, Result<MemberLines, String>)>(d * 2);
    let mut handles = Vec::with_capacity(d);
    for _ in 0..d {
        let path = path.clone();
        let members = members.clone();
        let next = next.clone();
        let tx = tx.clone();
        handles.push(std::thread::spawn(move || {
            let mut file = match File::open(&path) {
                Ok(f) => f,
                Err(_e) => return, // 打开失败：reader 侧以 channel 关闭感知并报错
            };
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= members.len() {
                    return;
                }
                let m = &members[i];
                let result = decode_member(&mut file, m, i == members.len() - 1);
                if tx.send((i, result)).is_err() {
                    return; // reader 已丢弃（早退），线程退出
                }
            }
        }));
    }
    (rx, handles)
}

fn decode_member(file: &mut File, m: &Member, is_last: bool) -> Result<MemberLines, String> {
    if file.seek(SeekFrom::Start(m.off)).is_err() {
        return Err(format!("成员 {} 定位失败", m.off));
    }
    let cap = initial_capacity(file, m);
    // initial_capacity 会 seek 到 trailer 读 ISIZE：回到成员起点再 take
    if file.seek(SeekFrom::Start(m.off)).is_err() {
        return Err(format!("成员 {} 定位失败", m.off));
    }
    let limited = (&mut *file).take(m.clen);
    let mut decoder = GzDecoder::new(limited);
    let mut out = Vec::with_capacity(cap);
    match decoder.read_to_end(&mut out) {
        Ok(_) => Ok(split_lines(out)),
        Err(e) => {
            if is_last {
                // 退化末成员容错：截断 trailer 视为干净 EOF，保留已解出的字节；
                // 仅处理解码器明确报错的截断情形。
                Ok(split_lines(out))
            } else {
                Err(format!("gzip 成员 {} 解压失败: {e}", m.off))
            }
        }
    }
}

/// 解压缓冲预容量：用 trailer 的 ISIZE（RFC 1952 §2.3.1，解压后大小 mod 2^32）
/// 一次到位，避免 read_to_end 倍增重分配；ISIZE 缺失/越界（截断成员、单成员
/// 输出 >1GB、压缩比 >8:1 的防御上界）时退回 3×压缩长。
fn initial_capacity(file: &mut File, m: &Member) -> usize {
    let est = 3usize.saturating_mul(m.clen as usize);
    if m.clen < 18 {
        return est; // 连头+trailer 都不完整：无从读 ISIZE
    }
    let mut isz = [0u8; 4];
    if file
        .seek(SeekFrom::Start(m.off + m.clen - 4))
        .and_then(|_| file.read_exact(&mut isz))
        .is_err()
    {
        return est;
    }
    let isize = u32::from_le_bytes(isz) as usize;
    if isize > 0 && isize <= (1 << 30) && isize <= 8usize.saturating_mul(m.clen as usize) {
        isize
    } else {
        est
    }
}

/// 行引用：共享成员缓冲上的字节区间。
#[derive(Clone)]
struct LineRef {
    buf: Arc<Vec<u8>>,
    r: (u32, u32),
}

impl LineRef {
    fn as_bytes(&self) -> &[u8] {
        &self.buf[self.r.0 as usize..self.r.1 as usize]
    }
}

/// 成员级 主线程侧"记录泵"：按成员序接收解压线程切好的行（乱序暂存 pending，
/// 绝不提前 EOF），维护跨成员拼接的行队列；`pop_line` 供 4 行一组建记录。
/// 全部成员消费完后把末尾部分行 flush 为最后一行。
struct LinePump {
    rx: Option<Receiver<(usize, Result<MemberLines, String>)>>,
    handles: Option<Vec<JoinHandle<()>>>,
    n_members: usize,
    next_seq: usize,
    pending: BTreeMap<usize, Result<MemberLines, String>>,
    queue: VecDeque<LineRef>,
    /// 跨成员拼接中的部分行（上一成员末尾无 \n 的段；仅边界处物化，每成员 ≤1 次）。
    partial: Option<Vec<u8>>,
    exhausted: bool,
    err: Option<String>,
}

impl LinePump {
    fn new(
        rx: Receiver<(usize, Result<MemberLines, String>)>,
        handles: Vec<JoinHandle<()>>,
        n_members: usize,
    ) -> Self {
        Self {
            rx: Some(rx),
            handles: Some(handles),
            n_members,
            next_seq: 0,
            pending: BTreeMap::new(),
            queue: VecDeque::new(),
            partial: None,
            exhausted: false,
            err: None,
        }
    }

    /// 队首行；队列空时继续按序摄入成员，直到有行或 EOF（None）。
    fn pop_line(&mut self) -> Result<Option<LineRef>, String> {
        if let Some(e) = &self.err {
            return Err(e.clone());
        }
        if let Some(l) = self.queue.pop_front() {
            return Ok(Some(l));
        }
        while self.next_seq < self.n_members {
            let m = self.recv_next_member()?;
            self.ingest(m);
            if let Some(l) = self.queue.pop_front() {
                return Ok(Some(l));
            }
        }
        if !self.exhausted {
            // 全部成员消费完：末尾部分行（无 \n 结束）作为最后一行，剥 \r
            // 与单线程 read_line_raw 语义一致。
            self.exhausted = true;
            if let Some(mut p) = self.partial.take() {
                if p.last() == Some(&b'\r') {
                    p.pop();
                }
                let len = p.len() as u32;
                self.queue.push_back(LineRef {
                    buf: Arc::new(p),
                    r: (0, len),
                });
            }
            if let Some(l) = self.queue.pop_front() {
                return Ok(Some(l));
            }
        }
        Ok(None)
    }

    /// 按序取下一个成员（乱序到达的暂存 pending；通道关闭 = 解压线程提前退出）。
    fn recv_next_member(&mut self) -> Result<MemberLines, String> {
        loop {
            if let Some(m) = self.pending.remove(&self.next_seq) {
                self.next_seq += 1;
                match &m {
                    Ok(_) => {}
                    Err(e) => self.err = Some(e.clone()),
                }
                return m;
            }
            match self.rx.as_ref().unwrap().recv() {
                Ok((seq, m)) => {
                    self.pending.insert(seq, m);
                }
                Err(_) => {
                    let e = "解压线程提前退出".to_string();
                    self.err = Some(e.clone());
                    return Err(e);
                }
            }
        }
    }

    /// 把成员的行并入队列，处理跨成员部分行拼接：
    /// - 上一成员的部分行并入本成员首行（本成员无行则继续等待）；
    /// - 本成员 tail（无 \n 的尾段）挂起为新的 partial。
    fn ingest(&mut self, m: MemberLines) {
        let mut lines = m.lines;
        if self.partial.is_some() {
            if let Some(&(fs, fe)) = lines.first() {
                let mut p = self.partial.take().unwrap();
                if fs == fe && p.last() == Some(&b'\r') {
                    // 成员边界恰在 \r 与 \n 之间：原行实为 \r\n 结尾，剥 \r
                    p.pop();
                }
                p.extend_from_slice(&m.buf[fs as usize..fe as usize]);
                let len = p.len() as u32;
                self.queue.push_back(LineRef {
                    buf: Arc::new(p),
                    r: (0, len),
                });
                lines.remove(0);
            }
            // lines 为空：partial 保持打开（tail 会并入或等下一成员）
        }
        if let Some((ts, te)) = m.tail {
            let tail = m.buf[ts as usize..te as usize].to_vec();
            match &mut self.partial {
                Some(p) => p.extend_from_slice(&tail),
                None => self.partial = Some(tail),
            }
        }
        for r in lines {
            self.queue.push_back(LineRef {
                buf: m.buf.clone(),
                r,
            });
        }
    }
}

impl Drop for LinePump {
    fn drop(&mut self) {
        // 先关通道让解压线程退出（send 失败即返回），再 join 避免线程泄漏。
        self.rx.take();
        if let Some(handles) = self.handles.take() {
            for h in handles {
                let _ = h.join();
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pair_id_normalization() {
        assert_eq!(normalize_pair_id("read1/1"), "read1");
        assert_eq!(normalize_pair_id("read1/2"), "read1");
        assert_eq!(normalize_pair_id("read1"), "read1");
        assert_eq!(normalize_pair_id(" read1/1 "), "read1");
        // Illumina CASAVA 1.8 描述字段（配对 ID 在空白前）
        assert_eq!(
            normalize_pair_id("INST:1:FC:1:1101 1:N:0:ATCACG"),
            "INST:1:FC:1:1101"
        );
        assert_eq!(
            normalize_pair_id("INST:1:FC:1:1101 2:N:0:ATCACG"),
            "INST:1:FC:1:1101"
        );
    }

    /// 并行测试使用进程 ID 与原子计数器生成唯一临时文件名。
    static TMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

    fn write_tmp(content: &str) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "viroflash_test_{}_{}.fq",
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn next_record_parses_four_line_records() {
        let p = write_tmp("@r1/1 desc\nACGTACGT\n+\nIIIIIIII\n@r2\nNN\n+\n!!\n");
        let mut rd = FastqReader::open(&p).unwrap();
        let a = rd.next_record().unwrap().unwrap();
        assert_eq!(a.id, "r1/1 desc");
        assert_eq!(a.seq, b"ACGTACGT");
        assert_eq!(a.qual, b"IIIIIIII");
        let b = rd.next_record().unwrap().unwrap();
        assert_eq!(b.id, "r2");
        assert_eq!(b.seq, b"NN");
        assert!(rd.next_record().unwrap().is_none());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn next_record_rejects_mismatched_seq_qual() {
        let p = write_tmp("@r1\nACGT\n+\nIII\n");
        let mut rd = FastqReader::open(&p).unwrap();
        let e = rd.next_record().unwrap_err();
        assert!(e.contains("长度不一致"), "{e}");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn pair_span_iter_rejects_unequal_record_counts() {
        let p1 = write_tmp("@a\nA\n+\nI\n@b\nC\n+\nI\n");
        let p2 = write_tmp("@a\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap(); // 非 gz 回退单线程
        assert!(pr.next().unwrap().is_ok()); // 第一条配对
        let e = pr.next().unwrap().unwrap_err();
        assert!(e.contains("记录数不一致"), "{e}");
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn pair_span_iter_rejects_mismatched_ids() {
        let p1 = write_tmp("@a/1\nA\n+\nI\n");
        let p2 = write_tmp("@a/2\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap();
        assert!(pr.next().unwrap().is_ok()); // /1 vs /2 归一化后一致
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
        let p1 = write_tmp("@a\nA\n+\nI\n");
        let p2 = write_tmp("@b\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap();
        let e = pr.next().unwrap().unwrap_err();
        assert!(e.contains("配对 ID 不一致"), "{e}");
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn span_path_matches_owned_path() {
        // 并行泵：next_record_span（零拷贝区间）与 next_record（物化）逐记录一致。
        let rec1 = "@a/1\nACGTACGT\n+\nIIIIIIII\n";
        let rec2 = "@b/1\nTTTTGGGGCCCC\n+\n!!!!!!!!!!!!\n";
        let body = format!("{rec1}{rec2}").into_bytes();
        let mut gz = Vec::new();
        gz.extend_from_slice(&gz_member(&body[..rec1.len()]));
        gz.extend_from_slice(&gz_member(&body[rec1.len()..]));
        let p = write_bin(&gz);
        let mut span = FastqReader::open_parallel(&p, 2).unwrap();
        let mut owned = FastqReader::open_parallel(&p, 2).unwrap();
        let mut n = 0;
        loop {
            match (
                span.next_record_span().unwrap(),
                owned.next_record().unwrap(),
            ) {
                (None, None) => break,
                (Some(x), Some(y)) => {
                    assert_eq!(x.id(), y.id.as_bytes());
                    assert_eq!(x.seq(), y.seq.as_slice());
                    assert_eq!(x.qual(), y.qual.as_slice());
                    n += 1;
                }
                other => panic!("span/owned 记录数不一致: {other:?}"),
            }
        }
        assert_eq!(n, 2);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn pair_ids_eq_matches_normalize_pair_id() {
        let cases = [
            "read1/1",
            "read1/2",
            "read1",
            " read1/1 ",
            "a/1 x",
            "a/2 y",
            "INST:1:FC:1:1101 1:N:0:ATCACG",
            "",
            "/1",
            " ",
            "x/1/2",
            "x/2/1",
        ];
        for &a in &cases {
            for &b in &cases {
                assert_eq!(
                    pair_ids_eq(a.as_bytes(), b.as_bytes()),
                    normalize_pair_id(a) == normalize_pair_id(b),
                    "a={a:?} b={b:?}"
                );
            }
        }
        // 随机 ASCII 串交叉一致性（splitmix64 固定种子）
        fn splitmix(state: &mut u64) -> u64 {
            *state ^= *state << 13;
            *state ^= *state >> 7;
            *state ^= *state << 17;
            *state
        }
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;
        let alphabet: &[u8] = b"ACGTN/12 xyz0123456789\t";
        let rand_bytes = |state: &mut u64| {
            let len = (splitmix(state) % 16) as usize;
            (0..len)
                .map(|_| alphabet[(splitmix(state) as usize) % alphabet.len()])
                .collect::<Vec<u8>>()
        };
        for _ in 0..2000 {
            let s1 = String::from_utf8(rand_bytes(&mut seed)).unwrap();
            let s2 = String::from_utf8(rand_bytes(&mut seed)).unwrap();
            assert_eq!(
                pair_ids_eq(s1.as_bytes(), s2.as_bytes()),
                normalize_pair_id(&s1) == normalize_pair_id(&s2),
                "s1={s1:?} s2={s2:?}"
            );
        }
    }

    // ---- 成员级并行解压测试 ----

    /// 合成一个带 IG extra 子域的 gzip 成员；IG=u32 LE 成员总压缩长。
    fn gz_member(payload: &[u8]) -> Vec<u8> {
        use flate2::write::DeflateEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(payload).unwrap();
        let def = enc.finish().unwrap();
        let mut crc = flate2::Crc::new();
        crc.update(payload);
        let total = 20u32 + def.len() as u32 + 8; // 头(含 FEXTRA) + deflate + trailer
        let mut m = vec![0x1f, 0x8b, 8, 0x04, 0, 0, 0, 0, 0, 0xff];
        m.extend_from_slice(&8u16.to_le_bytes()); // XLEN = 8（单个 IG 子域）
        m.extend_from_slice(b"IG");
        m.extend_from_slice(&4u16.to_le_bytes());
        m.extend_from_slice(&total.to_le_bytes());
        m.extend_from_slice(&def);
        m.extend_from_slice(&crc.sum().to_le_bytes());
        m.extend_from_slice(&(payload.len() as u32).to_le_bytes());
        m
    }

    fn write_bin(content: &[u8]) -> std::path::PathBuf {
        let path = std::env::temp_dir().join(format!(
            "viroflash_test_{}_{}.gz",
            std::process::id(),
            TMP_COUNTER.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::write(&path, content).unwrap();
        path
    }

    #[test]
    fn multi_member_parallel_matches_single_thread() {
        // 3 个 IG 成员（边界切在记录之间），d=2 并行解压与单线程路径逐记录一致。
        let rec1 = "@a/1\nACGTACGT\n+\nIIIIIIII\n";
        let rec2 = "@b/1\nTTTTGGGGCCCC\n+\n!!!!!!!!!!!!\n";
        let rec3 = "@c/1\nNNNNNN\n+\n######\n";
        let body = format!("{rec1}{rec2}{rec3}").into_bytes();
        let mut gz = Vec::new();
        for i in 0..3 {
            let seg = &body[if i == 0 {
                0
            } else {
                rec1.len() + (i - 1) * rec2.len()
            }..if i == 2 {
                body.len()
            } else {
                rec1.len() + i * rec2.len()
            }];
            gz.extend_from_slice(&gz_member(seg));
        }
        let p = write_bin(&gz);
        let mut par = FastqReader::open_parallel(&p, 2).unwrap();
        let mut seq = FastqReader::open(&p).unwrap();
        let mut n = 0;
        loop {
            match (par.next_record().unwrap(), seq.next_record().unwrap()) {
                (None, None) => break,
                (Some(x), Some(y)) => {
                    assert_eq!(x.id, y.id);
                    assert_eq!(x.seq, y.seq);
                    assert_eq!(x.qual, y.qual);
                    n += 1;
                }
                other => panic!("并行/单线程记录数不一致: {other:?}"),
            }
        }
        assert_eq!(n, 3);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn record_spanning_member_boundary() {
        // FASTQ 记录 4 行跨成员边界：成员 0 以半行序列结尾，成员 1 续完。
        let mut gz = gz_member(b"@x/1\nACGTAC");
        gz.extend_from_slice(&gz_member(b"GTACGT\n+\nIIIIIIIIIIII\n"));
        let p = write_bin(&gz);
        let mut rd = FastqReader::open_parallel(&p, 2).unwrap();
        let r = rd.next_record().unwrap().unwrap();
        assert_eq!(r.id, "x/1");
        assert_eq!(r.seq, b"ACGTACGTACGT");
        assert_eq!(r.qual, b"IIIIIIIIIIII");
        assert!(rd.next_record().unwrap().is_none());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn degenerate_last_member_tolerated_as_eof() {
        // 末成员：合法头 + 非法 deflate（0xFF → BTYPE=3）→ 解压报错，按干净 EOF 容错；
        // 前面成员的记录须完整读出且不报错。
        let mut gz = gz_member(b"@a/1\nACGT\n+\nIIII\n");
        let garbage = [0xFFu8; 16];
        let total = 20u32 + garbage.len() as u32; // 头(含 FEXTRA) + 垃圾，无 trailer
        let mut m = vec![0x1f, 0x8b, 8, 0x04, 0, 0, 0, 0, 0, 0xff];
        m.extend_from_slice(&8u16.to_le_bytes());
        m.extend_from_slice(b"IG");
        m.extend_from_slice(&4u16.to_le_bytes());
        m.extend_from_slice(&total.to_le_bytes());
        m.extend_from_slice(&garbage);
        gz.extend_from_slice(&m);
        let p = write_bin(&gz);
        let mut rd = FastqReader::open_parallel(&p, 2).unwrap();
        let r = rd.next_record().unwrap().unwrap();
        assert_eq!(r.seq, b"ACGT");
        assert!(rd.next_record().unwrap().is_none()); // 无错误、干净 EOF
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn ig_missing_falls_back_to_single_thread() {
        // 无 IG 子域的普通多成员 gzip（flate2 默认 FLG=0）→ 索引为空 → 回退单线程，
        // 记录须完整解析（吞吐不降级为错）。
        use flate2::write::GzEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(b"@a/1\nACGT\n+\nIIII\n").unwrap();
        let part1 = enc.finish().unwrap();
        let mut enc = GzEncoder::new(Vec::new(), Compression::default());
        enc.write_all(b"@b/1\nTGCA\n+\n!!!!\n").unwrap();
        let part2 = enc.finish().unwrap();
        let mut gz = part1;
        gz.extend_from_slice(&part2);
        let p = write_bin(&gz);
        assert!(build_member_index(&p).unwrap().is_empty()); // 无 IG → 空索引
        let mut rd = FastqReader::open_parallel(&p, 2).unwrap(); // 回退单线程
        let a = rd.next_record().unwrap().unwrap();
        let b = rd.next_record().unwrap().unwrap();
        assert_eq!(a.id, "a/1");
        assert_eq!(a.seq, b"ACGT");
        assert_eq!(b.id, "b/1");
        assert_eq!(b.seq, b"TGCA");
        assert!(rd.next_record().unwrap().is_none());
        std::fs::remove_file(&p).ok();
    }
}
