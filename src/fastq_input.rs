use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::MultiGzDecoder;
use sha2::{Digest, Sha256};

use crate::analysis_profile::hex_sha256;

const FASTQ_READER_BUFFER_BYTES: usize = 1 << 20;

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastqRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Fragment {
    pub ordinal: u64,
    pub id: String,
    pub r1: Vec<u8>,
    pub r2: Option<Vec<u8>>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputCensus {
    pub input_mode: &'static str,
    pub fragments: u64,
    pub input_digest: String,
    pub read_ends_per_fragment: u8,
}

pub struct FragmentReader {
    r1: FastqReader,
    r2: Option<FastqReader>,
    ordinal: u64,
}

impl FragmentReader {
    pub fn open(r1: &Path, r2: Option<&Path>) -> Result<Self, String> {
        Ok(Self {
            r1: FastqReader::open(r1)?,
            r2: r2.map(FastqReader::open).transpose()?,
            ordinal: 0,
        })
    }

    pub fn next_fragment(&mut self) -> Result<Option<Fragment>, String> {
        let left = self.r1.next_record()?;
        let right = self
            .r2
            .as_mut()
            .map(FastqReader::next_record)
            .transpose()?
            .flatten();
        match (left, right, self.r2.is_some()) {
            (None, None, _) => Ok(None),
            (Some(_), None, true) | (None, Some(_), true) => {
                Err("R1 and R2 have different record counts".into())
            }
            (Some(left), Some(right), true) => {
                if normalize_pair_id(&left.id) != normalize_pair_id(&right.id) {
                    return Err(format!(
                        "Paired IDs do not match: {} vs {}",
                        left.id, right.id
                    ));
                }
                let fragment = Fragment {
                    ordinal: self.ordinal,
                    id: normalize_pair_id(&left.id).to_string(),
                    r1: left.sequence,
                    r2: Some(right.sequence),
                };
                self.ordinal += 1;
                Ok(Some(fragment))
            }
            (Some(left), None, false) => {
                let fragment = Fragment {
                    ordinal: self.ordinal,
                    id: normalize_pair_id(&left.id).to_string(),
                    r1: left.sequence,
                    r2: None,
                };
                self.ordinal += 1;
                Ok(Some(fragment))
            }
            (None, Some(_), false) | (Some(_), Some(_), false) => unreachable!(),
        }
    }

    pub fn input_digest(&self) -> String {
        composite_input_digest(
            self.r1.decoded_digest(),
            self.r2.as_ref().map(FastqReader::decoded_digest),
        )
    }
}

struct DigestingReader<R> {
    inner: R,
    hasher: Sha256,
}

impl<R> DigestingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            hasher: Sha256::new(),
        }
    }

    fn digest(&self) -> [u8; 32] {
        self.hasher.clone().finalize().into()
    }
}

impl<R: Read> Read for DigestingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let count = self.inner.read(buffer)?;
        self.hasher.update(&buffer[..count]);
        Ok(count)
    }
}

enum EncodedInput {
    Plain(File),
    Gzip(Box<MultiGzDecoder<File>>),
}

impl Read for EncodedInput {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        match self {
            Self::Plain(reader) => reader.read(buffer),
            Self::Gzip(reader) => reader.read(buffer),
        }
    }
}

struct FastqReader {
    reader: BufReader<DigestingReader<EncodedInput>>,
    path: PathBuf,
    line: Vec<u8>,
}

impl FastqReader {
    fn open(path: &Path) -> Result<Self, String> {
        let file =
            File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?;
        let input = if path.extension().is_some_and(|extension| extension == "gz") {
            EncodedInput::Gzip(Box::new(MultiGzDecoder::new(file)))
        } else {
            EncodedInput::Plain(file)
        };
        Ok(Self {
            reader: BufReader::with_capacity(
                FASTQ_READER_BUFFER_BYTES,
                DigestingReader::new(input),
            ),
            path: path.to_path_buf(),
            line: Vec::new(),
        })
    }

    fn decoded_digest(&self) -> [u8; 32] {
        self.reader.get_ref().digest()
    }

    fn line(&mut self) -> Result<Option<Vec<u8>>, String> {
        self.line.clear();
        let count = self
            .reader
            .read_until(b'\n', &mut self.line)
            .map_err(|error| format!("Failed to read {}: {error}", self.path.display()))?;
        if count == 0 {
            return Ok(None);
        }
        while matches!(self.line.last(), Some(b'\n' | b'\r')) {
            self.line.pop();
        }
        Ok(Some(self.line.clone()))
    }

    fn next_record(&mut self) -> Result<Option<FastqRecord>, String> {
        let name = match self.line()? {
            Some(value) => value,
            None => return Ok(None),
        };
        let sequence = self
            .line()?
            .ok_or_else(|| "Truncated FASTQ record (missing sequence line)".to_string())?;
        let plus = self
            .line()?
            .ok_or_else(|| "Truncated FASTQ record (missing + line)".to_string())?;
        let quality = self
            .line()?
            .ok_or_else(|| "Truncated FASTQ record (missing quality line)".to_string())?;
        if name.first() != Some(&b'@') || plus.first() != Some(&b'+') {
            return Err(format!("Invalid FASTQ record in {}", self.path.display()));
        }
        if sequence.len() != quality.len() {
            return Err(format!(
                "FASTQ sequence and quality lengths differ in {}",
                self.path.display()
            ));
        }
        let id = std::str::from_utf8(&name[1..])
            .map_err(|error| format!("FASTQ ID is not UTF-8: {error}"))?
            .trim();
        if id.is_empty() {
            return Err("FASTQ record ID is empty".into());
        }
        Ok(Some(FastqRecord {
            id: id.to_string(),
            sequence,
        }))
    }
}

fn normalize_pair_id(id: &str) -> &str {
    let token = id.split_ascii_whitespace().next().unwrap_or(id);
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
}

pub fn census_fastq(r1: &Path, r2: Option<&Path>) -> Result<InputCensus, String> {
    let mut reader = FragmentReader::open(r1, r2)?;
    let mut fragments = 0;
    while reader.next_fragment()?.is_some() {
        fragments += 1;
    }
    if fragments == 0 {
        return Err("FASTQ input contains no fragments".into());
    }
    Ok(InputCensus {
        input_mode: if r2.is_some() { "PE" } else { "SE" },
        fragments,
        input_digest: reader.input_digest(),
        read_ends_per_fragment: if r2.is_some() { 2 } else { 1 },
    })
}

fn composite_input_digest(r1: [u8; 32], r2: Option<[u8; 32]>) -> String {
    let mut bytes = b"viroflash-input-v1\0".to_vec();
    let mode = if r2.is_some() {
        b"PE".as_slice()
    } else {
        b"SE".as_slice()
    };
    append_digest_item(&mut bytes, b"mode", Sha256::digest(mode).into());
    append_digest_item(&mut bytes, b"r1", r1);
    if let Some(digest) = r2 {
        append_digest_item(&mut bytes, b"r2", digest);
    }
    hex_sha256(&bytes)
}

fn append_digest_item(bytes: &mut Vec<u8>, tag: &[u8], digest: [u8; 32]) {
    bytes.extend_from_slice(&(tag.len() as u16).to_be_bytes());
    bytes.extend_from_slice(tag);
    bytes.extend_from_slice(&digest);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn paired_census_counts_fragments_not_read_ends() {
        let root = std::env::temp_dir().join(format!("vf-fastq-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        for mate in [1, 2] {
            let mut file = File::create(root.join(format!("r{mate}.fq"))).unwrap();
            writeln!(file, "@a/{mate}\nACGT\n+\nIIII\n@b/{mate}\nTGCA\n+\nIIII").unwrap();
        }
        let census = census_fastq(&root.join("r1.fq"), Some(&root.join("r2.fq"))).unwrap();
        assert_eq!(census.fragments, 2);
        assert_eq!(census.read_ends_per_fragment, 2);
        let _ = std::fs::remove_dir_all(root);
    }
}
