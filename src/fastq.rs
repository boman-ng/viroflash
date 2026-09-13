use std::fs::File;
use std::io::{self, BufRead, Read};
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::mpsc::{sync_channel, Receiver, SyncSender};
use std::thread::JoinHandle;

use flate2::read::MultiGzDecoder;
use sha2::{Digest, Sha256};

use crate::profile::hex_sha256;

const FASTQ_READER_BUFFER_BYTES: usize = 1 << 20;
pub(crate) const FASTQ_BATCH_RECORDS: usize = 1 << 10;
const FASTQ_PREFETCH_BATCHES: usize = 2;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Fragment<'a> {
    pub ordinal: u64,
    pub id: &'a str,
    pub r1: &'a [u8],
    pub r2: Option<&'a [u8]>,
}

#[derive(Debug, Default)]
struct FastqBatch {
    ids: String,
    sequences: Vec<u8>,
    records: Vec<RecordOffsets>,
}

#[derive(Debug)]
struct RecordOffsets {
    id: Range<usize>,
    sequence: Range<usize>,
}

impl FastqBatch {
    fn push(&mut self, id: &str, sequence: &[u8]) {
        let id_start = self.ids.len();
        let sequence_start = self.sequences.len();
        self.ids.push_str(id);
        self.sequences.extend_from_slice(sequence);
        self.records.push(RecordOffsets {
            id: id_start..self.ids.len(),
            sequence: sequence_start..self.sequences.len(),
        });
    }
}

#[derive(Debug)]
pub(crate) struct FragmentBatch {
    ordinals: Vec<u64>,
    r1: FastqBatch,
    r2: Option<FastqBatch>,
}

impl FragmentBatch {
    pub fn len(&self) -> usize {
        self.r1.records.len()
    }

    pub fn fragments(&self) -> impl Iterator<Item = Result<Fragment<'_>, String>> {
        self.r1.records.iter().enumerate().map(|(offset, left)| {
            let left_id = &self.r1.ids[left.id.clone()];
            let right = self.r2.as_ref().map(|batch| {
                let record = &batch.records[offset];
                (
                    &batch.ids[record.id.clone()],
                    &batch.sequences[record.sequence.clone()],
                )
            });
            if let Some((right_id, _)) = right {
                if left_id != right_id {
                    return Err(format!("Paired IDs do not match: {left_id} vs {right_id}"));
                }
            }
            Ok(Fragment {
                ordinal: self.ordinals[offset],
                id: left_id,
                r1: &self.r1.sequences[left.sequence.clone()],
                r2: right.map(|(_, sequence)| sequence),
            })
        })
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputCensus {
    pub input_mode: &'static str,
    pub fragments: u64,
    pub input_digest: String,
    pub read_ends_per_fragment: u8,
}

pub struct FragmentReader {
    r1: PrefetchedFastqReader,
    r2: Option<PrefetchedFastqReader>,
    ordinal: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct FileIdentity {
    length: u64,
    modified: std::time::SystemTime,
    #[cfg(unix)]
    device: u64,
    #[cfg(unix)]
    inode: u64,
    #[cfg(unix)]
    changed: (i64, i64),
}

impl FileIdentity {
    fn from_metadata(metadata: &std::fs::Metadata) -> Result<Self, String> {
        #[cfg(unix)]
        use std::os::unix::fs::MetadataExt;
        Ok(Self {
            length: metadata.len(),
            modified: metadata
                .modified()
                .map_err(|e| format!("Cannot inspect FASTQ modification time: {e}"))?,
            #[cfg(unix)]
            device: metadata.dev(),
            #[cfg(unix)]
            inode: metadata.ino(),
            #[cfg(unix)]
            changed: (metadata.ctime(), metadata.ctime_nsec()),
        })
    }
}

impl FragmentReader {
    pub fn open(r1: &Path, r2: Option<&Path>) -> Result<Self, String> {
        Ok(Self {
            r1: PrefetchedFastqReader::open(r1)?,
            r2: r2.map(PrefetchedFastqReader::open).transpose()?,
            ordinal: 0,
        })
    }

    pub fn next_batch(&mut self) -> Result<Option<FragmentBatch>, String> {
        let left = self.r1.next_batch()?;
        let right = self
            .r2
            .as_mut()
            .map(PrefetchedFastqReader::next_batch)
            .transpose()?
            .flatten();
        match (left, right, self.r2.is_some()) {
            (None, None, _) => Ok(None),
            (Some(_), None, true) | (None, Some(_), true) => {
                Err("R1 and R2 have different record counts".into())
            }
            (Some(r1), r2, _) => {
                if r2
                    .as_ref()
                    .is_some_and(|right| right.records.len() != r1.records.len())
                {
                    return Err("R1 and R2 have different record counts".into());
                }
                let batch = FragmentBatch {
                    ordinals: (self.ordinal..self.ordinal + r1.records.len() as u64).collect(),
                    r1,
                    r2,
                };
                self.ordinal += batch.len() as u64;
                Ok(Some(batch))
            }
            (None, Some(_), false) => unreachable!(),
        }
    }

    pub fn input_digest(&self) -> Result<String, String> {
        Ok(composite_input_digest(
            b"viroflash-input-v1\0",
            self.r1.decoded_digest()?,
            self.r2
                .as_ref()
                .map(PrefetchedFastqReader::decoded_digest)
                .transpose()?,
        ))
    }
}

#[derive(Clone, Copy)]
struct FastqDigests {
    decoded: [u8; 32],
}

enum FastqBatchMessage {
    Records(FastqBatch),
    Complete(FastqDigests),
    Error(String),
}

struct PrefetchedFastqReader {
    receiver: Receiver<FastqBatchMessage>,
    digests: Option<FastqDigests>,
    worker: Option<JoinHandle<()>>,
    path: PathBuf,
    identity: FileIdentity,
}

impl PrefetchedFastqReader {
    fn open(path: &Path) -> Result<Self, String> {
        let file =
            File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?;
        let identity = FileIdentity::from_metadata(
            &file
                .metadata()
                .map_err(|e| format!("Cannot inspect FASTQ: {e}"))?,
        )?;
        let input = if path.extension().is_some_and(|extension| extension == "gz") {
            EncodedInput::Gzip(Box::new(MultiGzDecoder::new(file)))
        } else {
            EncodedInput::Plain(file)
        };
        let reader_path = path.to_path_buf();
        let (sender, receiver) = sync_channel(FASTQ_PREFETCH_BATCHES);
        let worker = std::thread::Builder::new()
            .name("viroflash-fastq-reader".into())
            .spawn(move || {
                std::thread::scope(|scope| {
                    let (decoded_sender, decoded_receiver) = sync_channel(FASTQ_PREFETCH_BATCHES);
                    scope.spawn(move || decode_fastq(input, decoded_sender));
                    let reader = FastqReader {
                        reader: DecodedReader {
                            receiver: decoded_receiver,
                            data: Vec::new(),
                            position: 0,
                            digest: None,
                        },
                        path: reader_path,
                        name: Vec::new(),
                        sequence: Vec::new(),
                        plus: Vec::new(),
                        quality: Vec::new(),
                    };
                    prefetch_fastq(reader, sender);
                });
            })
            .map_err(|error| {
                format!("Cannot start FASTQ reader for {}: {error}", path.display())
            })?;
        Ok(Self {
            receiver,
            digests: None,
            worker: Some(worker),
            path: path.to_path_buf(),
            identity,
        })
    }

    fn next_batch(&mut self) -> Result<Option<FastqBatch>, String> {
        if self.digests.is_some() {
            return Ok(None);
        }
        match self.receiver.recv() {
            Ok(FastqBatchMessage::Records(records)) => Ok(Some(records)),
            Ok(FastqBatchMessage::Complete(digests)) => {
                self.finish_worker()?;
                let metadata = std::fs::metadata(&self.path)
                    .map_err(|e| format!("Cannot inspect FASTQ {}: {e}", self.path.display()))?;
                if FileIdentity::from_metadata(&metadata)? != self.identity {
                    return Err(format!(
                        "FASTQ changed while being read: {}",
                        self.path.display()
                    ));
                }
                self.digests = Some(digests);
                Ok(None)
            }
            Ok(FastqBatchMessage::Error(error)) => {
                self.finish_worker()?;
                Err(error)
            }
            Err(_) => {
                self.finish_worker()?;
                Err(format!(
                    "FASTQ reader closed before end-of-file: {}",
                    self.path.display()
                ))
            }
        }
    }

    fn decoded_digest(&self) -> Result<[u8; 32], String> {
        self.digests.map(|digests| digests.decoded).ok_or_else(|| {
            format!(
                "FASTQ reader has not reached end-of-file: {}",
                self.path.display()
            )
        })
    }

    fn finish_worker(&mut self) -> Result<(), String> {
        let Some(worker) = self.worker.take() else {
            return Ok(());
        };
        worker.join().map_err(|_| {
            format!(
                "FASTQ reader worker panicked while reading {}",
                self.path.display()
            )
        })
    }
}

fn prefetch_fastq(mut reader: FastqReader, sender: SyncSender<FastqBatchMessage>) {
    loop {
        let mut records = FastqBatch {
            records: Vec::with_capacity(FASTQ_BATCH_RECORDS),
            ..FastqBatch::default()
        };
        let mut complete = false;
        while records.records.len() < FASTQ_BATCH_RECORDS {
            match reader.append_record(&mut records) {
                Ok(true) => {}
                Ok(false) => {
                    complete = true;
                    break;
                }
                Err(error) => {
                    if !records.records.is_empty()
                        && sender.send(FastqBatchMessage::Records(records)).is_err()
                    {
                        return;
                    }
                    let _ = sender.send(FastqBatchMessage::Error(error));
                    return;
                }
            }
        }
        if !records.records.is_empty() && sender.send(FastqBatchMessage::Records(records)).is_err()
        {
            return;
        }
        if complete {
            let _ = sender.send(FastqBatchMessage::Complete(FastqDigests {
                decoded: reader.decoded_digest(),
            }));
            return;
        }
    }
}

enum DecodedChunk {
    Data(Vec<u8>),
    Complete([u8; 32]),
}

fn decode_fastq(mut input: EncodedInput, sender: SyncSender<io::Result<DecodedChunk>>) {
    let mut hasher = Sha256::new();
    loop {
        let mut data = Vec::with_capacity(FASTQ_READER_BUFFER_BYTES);
        match input
            .by_ref()
            .take(FASTQ_READER_BUFFER_BYTES as u64)
            .read_to_end(&mut data)
        {
            Ok(0) => {
                let _ = sender.send(Ok(DecodedChunk::Complete(hasher.finalize().into())));
                return;
            }
            Ok(_) => {
                hasher.update(&data);
                if sender.send(Ok(DecodedChunk::Data(data))).is_err() {
                    return;
                }
            }
            Err(error) => {
                let _ = sender.send(Err(error));
                return;
            }
        }
    }
}

struct DecodedReader {
    receiver: Receiver<io::Result<DecodedChunk>>,
    data: Vec<u8>,
    position: usize,
    digest: Option<[u8; 32]>,
}

impl BufRead for DecodedReader {
    fn fill_buf(&mut self) -> io::Result<&[u8]> {
        if self.position == self.data.len() && self.digest.is_none() {
            match self.receiver.recv().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "FASTQ decoder closed before end-of-file",
                )
            })?? {
                DecodedChunk::Data(data) => {
                    self.data = data;
                    self.position = 0;
                }
                DecodedChunk::Complete(digest) => self.digest = Some(digest),
            }
        }
        Ok(&self.data[self.position..])
    }

    fn consume(&mut self, amount: usize) {
        self.position += amount;
    }
}

impl Read for DecodedReader {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let available = self.fill_buf()?;
        let count = available.len().min(buffer.len());
        buffer[..count].copy_from_slice(&available[..count]);
        self.consume(count);
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
    reader: DecodedReader,
    path: PathBuf,
    name: Vec<u8>,
    sequence: Vec<u8>,
    plus: Vec<u8>,
    quality: Vec<u8>,
}

impl FastqReader {
    fn decoded_digest(&self) -> [u8; 32] {
        self.reader
            .digest
            .expect("decoded input reached end-of-file")
    }

    fn append_record(&mut self, batch: &mut FastqBatch) -> Result<bool, String> {
        if !read_fastq_line(&mut self.reader, &self.path, &mut self.name)? {
            return Ok(false);
        }
        if !read_fastq_line(&mut self.reader, &self.path, &mut self.sequence)? {
            return Err("Truncated FASTQ record (missing sequence line)".into());
        }
        if !read_fastq_line(&mut self.reader, &self.path, &mut self.plus)? {
            return Err("Truncated FASTQ record (missing + line)".into());
        }
        if !read_fastq_line(&mut self.reader, &self.path, &mut self.quality)? {
            return Err("Truncated FASTQ record (missing quality line)".into());
        }
        if self.name.first() != Some(&b'@') || self.plus.first() != Some(&b'+') {
            return Err(format!("Invalid FASTQ record in {}", self.path.display()));
        }
        if self.sequence.len() != self.quality.len() {
            return Err(format!(
                "FASTQ sequence and quality lengths differ in {}",
                self.path.display()
            ));
        }
        let id = std::str::from_utf8(&self.name[1..])
            .map_err(|error| format!("FASTQ ID is not UTF-8: {error}"))?
            .trim();
        if id.is_empty() {
            return Err("FASTQ record ID is empty".into());
        }
        // Normalize only at the FASTQ boundary. Candidate IDs are already
        // canonical and may themselves end in /1 or /2.
        batch.push(normalize_pair_id(id), &self.sequence);
        Ok(true)
    }
}

fn read_fastq_line<R: BufRead>(
    reader: &mut R,
    path: &Path,
    line: &mut Vec<u8>,
) -> Result<bool, String> {
    line.clear();
    let count = reader
        .read_until(b'\n', line)
        .map_err(|error| format!("Failed to read {}: {error}", path.display()))?;
    if count == 0 {
        return Ok(false);
    }
    while matches!(line.last(), Some(b'\n' | b'\r')) {
        line.pop();
    }
    Ok(true)
}

fn normalize_pair_id(id: &str) -> &str {
    let token = id.split_ascii_whitespace().next().unwrap_or(id);
    token
        .strip_suffix("/1")
        .or_else(|| token.strip_suffix("/2"))
        .unwrap_or(token)
}

#[cfg(test)]
pub fn census_fastq(r1: &Path, r2: Option<&Path>) -> Result<InputCensus, String> {
    let mut reader = FragmentReader::open(r1, r2)?;
    let mut fragments = 0;
    while let Some(batch) = reader.next_batch()? {
        for fragment in batch.fragments() {
            fragment?;
            fragments += 1;
        }
    }
    if fragments == 0 {
        return Err("FASTQ input contains no fragments".into());
    }
    Ok(InputCensus {
        input_mode: if r2.is_some() { "PE" } else { "SE" },
        fragments,
        input_digest: reader.input_digest()?,
        read_ends_per_fragment: if r2.is_some() { 2 } else { 1 },
    })
}

fn composite_input_digest(domain: &[u8], r1: [u8; 32], r2: Option<[u8; 32]>) -> String {
    let mut bytes = domain.to_vec();
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
    fn batches_preserve_pair_ids_ordinals_sequences_and_raw_digests() {
        let root = std::env::temp_dir().join(format!("vf-fastq-batches-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let count = FASTQ_BATCH_RECORDS * 2 + 3;
        for mate in [1, 2] {
            let mut bytes = Vec::new();
            for ordinal in 0..count {
                let sequence = if ordinal == FASTQ_BATCH_RECORDS {
                    "A".repeat(FASTQ_READER_BUFFER_BYTES + 17)
                } else {
                    "ACGTN".into()
                };
                write!(
                    bytes,
                    "@read-{ordinal}/{mate} comment\r\n{sequence}\r\n+\r\n{}\r\n",
                    "I".repeat(sequence.len())
                )
                .unwrap();
            }
            std::fs::write(root.join(format!("r{mate}.fq")), &bytes).unwrap();
        }
        let r1 = root.join("r1.fq");
        let r2 = root.join("r2.fq");
        let mut reader = FragmentReader::open(&r1, Some(&r2)).unwrap();
        let mut observed = 0;
        while let Some(batch) = reader.next_batch().unwrap() {
            assert!(batch.len() <= FASTQ_BATCH_RECORDS);
            for fragment in batch.fragments() {
                let fragment = fragment.unwrap();
                assert_eq!(fragment.ordinal, observed as u64);
                assert_eq!(fragment.id, format!("read-{observed}"));
                if observed == FASTQ_BATCH_RECORDS {
                    assert_eq!(fragment.r1.len(), FASTQ_READER_BUFFER_BYTES + 17);
                    assert!(fragment.r1.iter().all(|&base| base == b'A'));
                } else {
                    assert_eq!(fragment.r1, b"ACGTN");
                }
                assert_eq!(fragment.r2, Some(fragment.r1));
                observed += 1;
            }
        }
        assert_eq!(observed, count);
        let expected = composite_input_digest(
            b"viroflash-input-v1\0",
            Sha256::digest(std::fs::read(&r1).unwrap()).into(),
            Some(Sha256::digest(std::fs::read(&r2).unwrap()).into()),
        );
        assert_eq!(reader.input_digest().unwrap(), expected);
        let census = census_fastq(&r1, Some(&r2)).unwrap();
        assert_eq!(census.fragments, count as u64);
        assert_eq!(census.input_digest, expected);
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn errors_after_a_full_batch_are_not_lost() {
        let root =
            std::env::temp_dir().join(format!("vf-fastq-batch-errors-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let prefix = "@read/1\nACGT\n+\nIIII\n".repeat(FASTQ_BATCH_RECORDS);
        let r1 = root.join("r1.fq");
        let r2 = root.join("r2.fq");
        for (tail, expected) in [
            ("@tail\nACGT\n+\n", "missing quality"),
            ("@tail\nACGT\n+\nIII\n", "lengths differ"),
            ("@ \nACGT\n+\nIIII\n", "ID is empty"),
            ("tail\nACGT\n+\nIIII\n", "Invalid FASTQ"),
        ] {
            std::fs::write(&r1, format!("{prefix}{tail}")).unwrap();
            assert!(census_fastq(&r1, None).unwrap_err().contains(expected));
        }
        std::fs::write(&r1, format!("{prefix}@left/1\nACGT\n+\nIIII\n")).unwrap();
        std::fs::write(&r2, prefix.replace("/1", "/2")).unwrap();
        assert!(census_fastq(&r1, Some(&r2))
            .unwrap_err()
            .contains("different record counts"));
        std::fs::write(
            &r2,
            format!("{}@right/2\nACGT\n+\nIIII\n", prefix.replace("/1", "/2")),
        )
        .unwrap();
        assert!(census_fastq(&r1, Some(&r2))
            .unwrap_err()
            .contains("Paired IDs do not match"));
        // A parser error must also release a decoder blocked on prefetched data.
        let mut malformed = b"@bad\nACGT\n+\nIII\n".to_vec();
        malformed.resize(FASTQ_READER_BUFFER_BYTES * 5, b'A');
        std::fs::write(&r1, malformed).unwrap();
        assert!(census_fastq(&r1, None)
            .unwrap_err()
            .contains("lengths differ"));
        let mut gzip = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        gzip.write_all(prefix.as_bytes()).unwrap();
        let mut compressed = gzip.finish().unwrap();
        compressed.truncate(compressed.len() - 4);
        let path = root.join("truncated.fq.gz");
        std::fs::write(&path, compressed).unwrap();
        assert!(census_fastq(&path, None).is_err());
        std::fs::remove_dir_all(root).unwrap();
    }

    #[test]
    fn input_changes_during_reading_are_rejected() {
        let root = std::env::temp_dir().join(format!("vf-fastq-change-{}", std::process::id()));
        std::fs::create_dir_all(&root).unwrap();
        let path = root.join("sample.fastq");
        // More batches than the prefetch queue: the reader cannot have completed.
        std::fs::write(
            &path,
            b"@a\nACGT\n+\nIIII\n".repeat(FASTQ_BATCH_RECORDS * 8),
        )
        .unwrap();
        let mut reader = FragmentReader::open(&path, None).unwrap();
        assert!(reader.next_batch().unwrap().is_some());
        let mut file = std::fs::OpenOptions::new()
            .append(true)
            .open(&path)
            .unwrap();
        file.write_all(b"@b\nTGCA\n+\nIIII\n").unwrap();
        file.flush().unwrap();
        let error = loop {
            match reader.next_batch() {
                Ok(Some(_)) => {}
                Ok(None) => panic!("changed input accepted"),
                Err(error) => break error,
            }
        };
        assert!(error.contains("changed while being read"), "{error}");
        let _ = std::fs::remove_dir_all(root);
    }
}
