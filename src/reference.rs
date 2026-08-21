//! 参考序列加载与 minimap2 索引构建。
//! 用户输入 4 类 FASTA：宿主/人类、目标/病毒、诱饵、常见污染。
//! 合并为单一 composite FASTA（contig 重命名为 `{role}_{i}`），并构建 sr 预设索引。

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use minimap2::Aligner;

/// 参考序列类别。污染类在判定时按宿主对待（只作噪声，不产出候选）。
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Role {
    Host,
    Target,
    Decoy,
    Contaminant,
}

impl Role {
    pub fn prefix(self) -> &'static str {
        match self {
            Role::Host => "host",
            Role::Target => "target",
            Role::Decoy => "decoy",
            Role::Contaminant => "contam",
        }
    }
}

#[derive(Debug, Clone)]
pub struct Contig {
    pub name: String,
    pub role: Role,
    pub seq: Vec<u8>, // 大写 ACGTN
    pub gc_frac: f64,
}

impl Contig {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// contig 元数据（索引 manifest 中持久化的最小信息）：下游聚合与统计判定
/// 只依赖 name/role/len/gc，不依赖序列本身，加载索引时无需保留序列。
#[derive(Debug, Clone, PartialEq)]
pub struct ContigMeta {
    pub name: String,
    pub role: Role,
    pub len: u64,
    pub gc_frac: f64,
}

impl From<&Contig> for ContigMeta {
    fn from(c: &Contig) -> Self {
        Self {
            name: c.name.clone(),
            role: c.role,
            len: c.len() as u64,
            gc_frac: c.gc_frac,
        }
    }
}

/// 解析 FASTA（允许序列换行；仅保留 ACGTN，忽略其他字符）。
pub fn parse_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path).map_err(|e| format!("无法打开 {}: {e}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("读取 {} 失败: {e}", path.display()))?;
        if let Some(h) = line.strip_prefix('>') {
            if let Some(hdr) = header.take() {
                records.push((hdr, std::mem::take(&mut seq)));
            }
            header = Some(h.split_whitespace().next().unwrap_or(h).to_string());
        } else if header.is_some() {
            for &c in line.as_bytes() {
                let u = c.to_ascii_uppercase();
                if matches!(u, b'A' | b'C' | b'G' | b'T' | b'N') {
                    seq.push(u);
                }
            }
        }
    }
    if let Some(hdr) = header.take() {
        records.push((hdr, seq));
    }
    Ok(records)
}

pub fn gc_fraction(seq: &[u8]) -> f64 {
    let acgt = seq.iter().filter(|&&c| c != b'N').count();
    if acgt == 0 {
        return 0.0;
    }
    let gc = seq.iter().filter(|&&c| c == b'G' || c == b'C').count();
    gc as f64 / acgt as f64
}

/// 加载 4 类 FASTA，重命名 contig 为 `{role}_{i}`，合并写出 composite FASTA，
/// 并构建 minimap2 sr 索引。返回 `(索引路径, 全部 contig)`。
pub fn build_reference(
    fastas: &[(Role, PathBuf)],
    work_dir: &Path,
    index_threads: usize,
) -> Result<(PathBuf, Vec<Contig>), String> {
    std::fs::create_dir_all(work_dir)
        .map_err(|e| format!("无法创建工作目录 {}: {e}", work_dir.display()))?;
    let composite_path = work_dir.join("composite.fa");
    let mmi_path = work_dir.join("index.mmi");

    let mut contigs = Vec::new();
    for (role, path) in fastas {
        // 每个角色独立编号：host_0、target_0、decoy_0…（避免全局计数导致目标≠target_0）
        for (counter, (orig_header, seq)) in parse_fasta(path)?.into_iter().enumerate() {
            if seq.is_empty() {
                return Err(format!("{} 中存在空序列: {orig_header}", path.display()));
            }
            let name = format!("{}_{}", role.prefix(), counter);
            contigs.push(Contig {
                name: name.clone(),
                role: *role,
                gc_frac: gc_fraction(&seq),
                seq,
            });
        }
    }

    {
        let file = File::create(&composite_path)
            .map_err(|e| format!("无法创建 {}: {e}", composite_path.display()))?;
        let mut writer = BufWriter::new(file);
        for contig in &contigs {
            writeln!(writer, ">{}", contig.name)
                .map_err(|e| format!("写入 composite FASTA 失败: {e}"))?;
            for chunk in contig.seq.chunks(60) {
                writer
                    .write_all(chunk)
                    .map_err(|e| format!("写入 composite FASTA 失败: {e}"))?;
                writer
                    .write_all(b"\n")
                    .map_err(|e| format!("写入 composite FASTA 失败: {e}"))?;
            }
        }
        writer
            .flush()
            .map_err(|e| format!("写入 composite FASTA 失败: {e}"))?;
    }

    let mmi_str = mmi_path
        .to_str()
        .ok_or_else(|| "索引路径非 UTF-8".to_string())?;
    Aligner::builder()
        .sr()
        .with_index_threads(index_threads.max(1))
        .with_index(&composite_path, Some(mmi_str))
        .map_err(|e| format!("minimap2 索引构建失败: {e}"))?;

    Ok((mmi_path, contigs))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_wrapped_fasta() {
        let dir = std::env::temp_dir().join(format!("vf_ref_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("a.fa");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">seq1 desc").unwrap();
            writeln!(f, "ACGTacgt").unwrap();
            writeln!(f, "NNNN").unwrap();
        }
        let records = parse_fasta(&path).unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].0, "seq1");
        assert_eq!(records[0].1, b"ACGTACGTNNNN");
        let _ = std::fs::remove_dir_all(&dir);
    }
}
