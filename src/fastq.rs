//! Streaming FASTQ and FASTQ.gz input with paired-ID validation.
//!
//! The single-threaded path reuses line buffers and never loads the full input. Indexed
//! multi-member gzip files may use bounded parallel decompression: workers decode independent
//! members, while the reader restores member order, joins lines crossing member boundaries, and
//! assembles four-line records. Missing or malformed member indexes fall back to the equivalent
//! single-threaded path. Both paths share record validation and produce byte-equivalent output.

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

#[derive(Debug)]
pub enum SpanRecord {
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

#[derive(Debug)]
pub struct SpanPair {
    pub r1: SpanRecord,
    pub r2: SpanRecord,
}

struct StreamFeed {
    reader: Box<dyn BufRead>,
    line: Vec<u8>,
}

enum Feed {
    Stream(StreamFeed),
    Pump(LinePump),
}

pub struct FastqReader {
    feed: Feed,
}

impl FastqReader {
    pub fn open(path: &Path) -> Result<Self, String> {
        let file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
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

    pub fn open_parallel(path: &Path, d: usize) -> Result<Self, String> {
        if d == 0 || !path.to_string_lossy().ends_with(".gz") {
            return Self::open(path);
        }
        let members = build_member_index(path)?;
        if members.len() < 2 {
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

    pub fn next_record_span(&mut self) -> Result<Option<SpanRecord>, String> {
        match &mut self.feed {
            Feed::Stream(f) => next_record_stream(f).map(|r| r.map(SpanRecord::Owned)),
            Feed::Pump(p) => next_record_pump_span(p),
        }
    }
}

fn read_line_raw(feed: &mut StreamFeed) -> Result<Option<&[u8]>, String> {
    feed.line.clear();
    let n = feed
        .reader
        .read_until(b'\n', &mut feed.line)
        .map_err(|e| format!("Failed to read FASTQ: {e}"))?;
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

fn next_record_stream(feed: &mut StreamFeed) -> Result<Option<FastqRecord>, String> {
    let name_line = match read_line_raw(feed)? {
        None => return Ok(None),
        Some(l) => l.to_vec(),
    };
    let seq = read_line_raw(feed)?
        .ok_or_else(|| "Truncated FASTQ record (missing sequence line)".to_string())?
        .to_vec();
    let plus = read_line_raw(feed)?
        .ok_or_else(|| "Truncated FASTQ record (missing + line)".to_string())?
        .to_vec();
    let qual = read_line_raw(feed)?
        .ok_or_else(|| "Truncated FASTQ record (missing quality line)".to_string())?
        .to_vec();
    assemble_record(name_line, seq, plus, qual).map(Some)
}

fn next_record_pump(pump: &mut LinePump) -> Result<Option<FastqRecord>, String> {
    let name_line = match pump.pop_line()? {
        None => return Ok(None),
        Some(l) => l.as_bytes().to_vec(),
    };
    let seq = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing sequence line)".to_string())?
        .as_bytes()
        .to_vec();
    let plus = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing + line)".to_string())?
        .as_bytes()
        .to_vec();
    let qual = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing quality line)".to_string())?
        .as_bytes()
        .to_vec();
    assemble_record(name_line, seq, plus, qual).map(Some)
}

fn next_record_pump_span(pump: &mut LinePump) -> Result<Option<SpanRecord>, String> {
    let name = match pump.pop_line()? {
        None => return Ok(None),
        Some(l) => l,
    };
    let seq = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing sequence line)".to_string())?;
    let plus = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing + line)".to_string())?;
    let qual = pump
        .pop_line()?
        .ok_or_else(|| "Truncated FASTQ record (missing quality line)".to_string())?;

    if plus.as_bytes().first() != Some(&b'+') {
        return Err(format!(
            "Invalid FASTQ + line: {:?}",
            String::from_utf8_lossy(plus.as_bytes())
        ));
    }
    let name_bytes = name.as_bytes();
    if name_bytes.first() != Some(&b'@') {
        return Err(format!(
            "Invalid FASTQ record: {:?}",
            String::from_utf8_lossy(name_bytes)
        ));
    }

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
            "FASTQ sequence and quality lengths differ: {}",
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

fn assemble_record(
    name_line: Vec<u8>,
    seq: Vec<u8>,
    plus: Vec<u8>,
    qual: Vec<u8>,
) -> Result<FastqRecord, String> {
    if plus.first() != Some(&b'+') {
        return Err(format!(
            "Invalid FASTQ + line: {:?}",
            String::from_utf8_lossy(&plus)
        ));
    }
    if name_line.first() != Some(&b'@') {
        return Err(format!(
            "Invalid FASTQ record: {:?}",
            String::from_utf8_lossy(&name_line)
        ));
    }
    let id =
        String::from_utf8(name_line).map_err(|e| format!("FASTQ ID is not valid UTF-8: {e}"))?;
    let id = id.trim_start_matches('@').trim().to_string();
    if seq.len() != qual.len() {
        return Err(format!("FASTQ sequence and quality lengths differ: {}", id));
    }
    Ok(FastqRecord { id, seq, qual })
}

pub fn normalize_pair_id(id: &str) -> String {
    let token = id.split_ascii_whitespace().next().unwrap_or(id);
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
        .to_string()
}

pub fn pair_ids_eq(a: &[u8], b: &[u8]) -> bool {
    fn norm(id: &[u8]) -> &[u8] {
        let mut s = 0;
        while s < id.len() && id[s].is_ascii_whitespace() {
            s += 1;
        }
        if s == id.len() {
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

pub struct PairSpanIter {
    r1: FastqReader,
    r2: FastqReader,
}

impl PairSpanIter {
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
                Some(Err("R1 and R2 have different record counts".to_string()))
            }
            (Ok(Some(a)), Ok(Some(b))) => {
                if !pair_ids_eq(a.id(), b.id()) {
                    return Some(Err(format!(
                        "Paired IDs do not match: {} vs {}",
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

//

// ---------------------------------------------------------------------------

#[derive(Debug, Clone)]
struct Member {
    off: u64,
    clen: u64,
}

const MAX_XLEN: usize = 64 << 10;
const MAX_NUL_FIELD: usize = 64 << 10;
const MAX_MEMBERS: usize = 16 << 20;

fn build_member_index(path: &Path) -> Result<Vec<Member>, String> {
    let mut file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
    let size = file
        .metadata()
        .map_err(|e| format!("Failed to read metadata for {}: {e}", path.display()))?
        .len();
    let mut members = Vec::new();
    let mut off: u64 = 0;

    loop {
        if members.len() >= MAX_MEMBERS {
            return Ok(Vec::new());
        }

        let Some(hdr) = parse_member_header(&mut file, off)? else {
            return Ok(Vec::new());
        };
        let header_end = hdr.end;
        let Some(ig) = hdr.ig else {
            return Ok(Vec::new());
        };
        if ig == 0 {
            return Ok(Vec::new());
        }

        let next_total = off.checked_add(ig as u64);
        let next_data = header_end.checked_add(ig as u64);
        let next = match (next_total, next_data) {
            (Some(a), Some(_b)) if a == size || next_is_magic(&mut file, a)? => a,
            (_, Some(b)) if b == size || next_is_magic(&mut file, b)? => b,
            _ => return Ok(Vec::new()),
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

struct ParsedHeader {
    end: u64,

    ig: Option<u32>,
}

fn parse_member_header(file: &mut File, off: u64) -> Result<Option<ParsedHeader>, String> {
    let mut hdr = [0u8; 10];
    match file
        .seek(SeekFrom::Start(off))
        .and_then(|_| read_exact_at(file, &mut hdr))
    {
        Ok(()) => {}
        Err(_) => return Ok(None),
    }
    if hdr[0] != 0x1f || hdr[1] != 0x8b || hdr[2] != 8 {
        return Ok(None);
    }
    let flg = hdr[3];
    let mut p = off + 10;
    let mut ig = None;

    if flg & 0x04 != 0 {
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
        let Some(n) = skip_to_nul(file, p)? else {
            return Ok(None);
        };
        p += n;
    }
    if flg & 0x10 != 0 {
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

fn skip_to_nul(file: &mut File, p: u64) -> Result<Option<u64>, String> {
    use std::io::Read;
    file.seek(SeekFrom::Start(p))
        .map_err(|e| format!("Failed to seek within gzip header: {e}"))?;
    let mut buf = [0u8; 256];
    let mut total = 0u64;
    loop {
        let n = file
            .read(&mut buf)
            .map_err(|e| format!("Failed to read gzip header: {e}"))?;
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

struct MemberLines {
    buf: Arc<Vec<u8>>,

    lines: Vec<(u32, u32)>,
    tail: Option<(u32, u32)>,
}

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
                Err(_e) => return,
            };
            loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= members.len() {
                    return;
                }
                let m = &members[i];
                let result = decode_member(&mut file, m, i == members.len() - 1);
                if tx.send((i, result)).is_err() {
                    return;
                }
            }
        }));
    }
    (rx, handles)
}

fn decode_member(file: &mut File, m: &Member, is_last: bool) -> Result<MemberLines, String> {
    if file.seek(SeekFrom::Start(m.off)).is_err() {
        return Err(format!("Failed to seek to gzip member {}", m.off));
    }
    let cap = initial_capacity(file, m);

    if file.seek(SeekFrom::Start(m.off)).is_err() {
        return Err(format!("Failed to seek to gzip member {}", m.off));
    }
    let limited = (&mut *file).take(m.clen);
    let mut decoder = GzDecoder::new(limited);
    let mut out = Vec::with_capacity(cap);
    match decoder.read_to_end(&mut out) {
        Ok(_) => Ok(split_lines(out)),
        Err(e) => {
            if is_last {
                Ok(split_lines(out))
            } else {
                Err(format!("Failed to decompress gzip member {}: {e}", m.off))
            }
        }
    }
}

fn initial_capacity(file: &mut File, m: &Member) -> usize {
    let est = 3usize.saturating_mul(m.clen as usize);
    if m.clen < 18 {
        return est;
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

struct LinePump {
    rx: Option<Receiver<(usize, Result<MemberLines, String>)>>,
    handles: Option<Vec<JoinHandle<()>>>,
    n_members: usize,
    next_seq: usize,
    pending: BTreeMap<usize, Result<MemberLines, String>>,
    queue: VecDeque<LineRef>,

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
                    let e = "Decompression thread exited early".to_string();
                    self.err = Some(e.clone());
                    return Err(e);
                }
            }
        }
    }

    fn ingest(&mut self, m: MemberLines) {
        let mut lines = m.lines;
        if self.partial.is_some() {
            if let Some(&(fs, fe)) = lines.first() {
                let mut p = self.partial.take().unwrap();
                if fs == fe && p.last() == Some(&b'\r') {
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

        assert_eq!(
            normalize_pair_id("INST:1:FC:1:1101 1:N:0:ATCACG"),
            "INST:1:FC:1:1101"
        );
        assert_eq!(
            normalize_pair_id("INST:1:FC:1:1101 2:N:0:ATCACG"),
            "INST:1:FC:1:1101"
        );
    }

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
        assert!(e.contains("lengths differ"), "{e}");
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn pair_span_iter_rejects_unequal_record_counts() {
        let p1 = write_tmp("@a\nA\n+\nI\n@b\nC\n+\nI\n");
        let p2 = write_tmp("@a\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap();
        assert!(pr.next().unwrap().is_ok());
        let e = pr.next().unwrap().unwrap_err();
        assert!(e.contains("different record counts"), "{e}");
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn pair_span_iter_rejects_mismatched_ids() {
        let p1 = write_tmp("@a/1\nA\n+\nI\n");
        let p2 = write_tmp("@a/2\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap();
        assert!(pr.next().unwrap().is_ok());
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
        let p1 = write_tmp("@a\nA\n+\nI\n");
        let p2 = write_tmp("@b\nA\n+\nI\n");
        let mut pr = PairSpanIter::open_parallel(&p1, &p2, 2).unwrap();
        let e = pr.next().unwrap().unwrap_err();
        assert!(e.contains("Paired IDs do not match"), "{e}");
        std::fs::remove_file(&p1).ok();
        std::fs::remove_file(&p2).ok();
    }

    #[test]
    fn span_path_matches_owned_path() {
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
                other => panic!("span/owned different record counts: {other:?}"),
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

    fn gz_member(payload: &[u8]) -> Vec<u8> {
        use flate2::write::DeflateEncoder;
        use flate2::Compression;
        use std::io::Write;
        let mut enc = DeflateEncoder::new(Vec::new(), Compression::default());
        enc.write_all(payload).unwrap();
        let def = enc.finish().unwrap();
        let mut crc = flate2::Crc::new();
        crc.update(payload);
        let total = 20u32 + def.len() as u32 + 8;
        let mut m = vec![0x1f, 0x8b, 8, 0x04, 0, 0, 0, 0, 0, 0xff];
        m.extend_from_slice(&8u16.to_le_bytes());
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
                other => panic!("parallel and single-threaded paths returned different record counts: {other:?}"),
            }
        }
        assert_eq!(n, 3);
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn record_spanning_member_boundary() {
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
        let mut gz = gz_member(b"@a/1\nACGT\n+\nIIII\n");
        let garbage = [0xFFu8; 16];
        let total = 20u32 + garbage.len() as u32;
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
        assert!(rd.next_record().unwrap().is_none());
        std::fs::remove_file(&p).ok();
    }

    #[test]
    fn ig_missing_falls_back_to_single_thread() {
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
        assert!(build_member_index(&p).unwrap().is_empty());
        let mut rd = FastqReader::open_parallel(&p, 2).unwrap();
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
