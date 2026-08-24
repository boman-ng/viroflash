//! Reference loading and minimap2 index construction.
//! Four FASTA roles—host, target, decoy, and contaminant—are merged into one composite FASTA.
//! Contigs are renamed to `{role}_{i}` before building an `sr` preset index.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use minimap2::Aligner;

/// Competitive alignment requires every role to share one minimap2 shard for MAPQ calculation.
const SINGLE_PART_BATCH_SIZE: u64 = u64::MAX;

pub(crate) fn ensure_single_part_index(part_count: usize) -> Result<(), &'static str> {
    if part_count == 1 {
        Ok(())
    } else {
        Err("The minimap2 index contains multiple shards; rebuild it as a single-shard index for competitive alignment")
    }
}

/// Reference role. Contaminants are background-only and never produce candidates.
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
    pub seq: Vec<u8>, // Uppercase ACGTN
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

/// Minimal contig metadata persisted in the index manifest. Downstream aggregation and testing
/// depend only on name, role, length, and GC content, so loaded indexes need not retain sequences.
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

/// Parse multiline FASTA, retaining only ACGTN characters.
pub fn parse_fasta(path: &Path) -> Result<Vec<(String, Vec<u8>)>, String> {
    let file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
    let reader = BufReader::new(file);
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    let mut header: Option<String> = None;
    let mut seq: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| format!("Failed to read {}: {e}", path.display()))?;
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

/// Load the four FASTA roles, rename contigs to `{role}_{i}`, write a composite FASTA, and build
/// a minimap2 `sr` index. Returns the index path and all contigs.
pub fn build_reference(
    fastas: &[(Role, PathBuf)],
    work_dir: &Path,
    index_threads: usize,
) -> Result<(PathBuf, Vec<Contig>), String> {
    std::fs::create_dir_all(work_dir)
        .map_err(|e| format!("Cannot create work directory {}: {e}", work_dir.display()))?;
    let composite_path = work_dir.join("composite.fa");
    let mmi_path = work_dir.join("index.mmi");

    let mut contigs = Vec::new();
    for (role, path) in fastas {
        // Number each role independently so the first target is always target_0.
        for (counter, (orig_header, seq)) in parse_fasta(path)?.into_iter().enumerate() {
            if seq.is_empty() {
                return Err(format!(
                    "{} contains an empty sequence: {orig_header}",
                    path.display()
                ));
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
            .map_err(|e| format!("Cannot create {}: {e}", composite_path.display()))?;
        let mut writer = BufWriter::new(file);
        for contig in &contigs {
            writeln!(writer, ">{}", contig.name)
                .map_err(|e| format!("Failed to write composite FASTA: {e}"))?;
            for chunk in contig.seq.chunks(60) {
                writer
                    .write_all(chunk)
                    .map_err(|e| format!("Failed to write composite FASTA: {e}"))?;
                writer
                    .write_all(b"\n")
                    .map_err(|e| format!("Failed to write composite FASTA: {e}"))?;
            }
        }
        writer
            .flush()
            .map_err(|e| format!("Failed to write composite FASTA: {e}"))?;
    }

    let mmi_str = mmi_path
        .to_str()
        .ok_or_else(|| "Index path is not valid UTF-8".to_string())?;
    let mut builder = Aligner::builder()
        .sr()
        .with_index_threads(index_threads.max(1));
    builder.idxopt.batch_size = SINGLE_PART_BATCH_SIZE;
    let aligner = builder
        .with_index(&composite_path, Some(mmi_str))
        .map_err(|e| format!("Failed to build minimap2 index: {e}"))?;
    ensure_single_part_index(aligner.idx_parts.len())?;

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
