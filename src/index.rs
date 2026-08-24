//! Build, load, and validate viroflash index directories.
//!
//! Index layout produced by `viroflash index --out <dir>`:
//! ```text
//! <dir>/
//!   ref.mmi        minimap2 `sr` index written through `with_index(fa, Some(out))`
//!   bloom.bin      versioned canonical k-mer Bloom filter for targets and decoys
//!   manifest.json  version, k, role-to-contig metadata, decoy provenance, and checksums
//!   targets.fa     detection-group representatives, the only Target inputs to the MMI
//!   decoys.fa      generated decoys, accompanied by a decoys.tsv metadata report
//! ```
//! Loading reads only the manifest, Bloom filter, and ref.mmi; source sequences are unnecessary.
//! Source FASTA paths and BLAKE3 checksums remain in the manifest for auditability. Construction
//! uses a temporary directory followed by an atomic rename to prevent races and partial indexes.
//!
//! JSON is serialized and parsed directly to preserve the no-serde convention in `report.rs`.
//! The parser intentionally supports only the subset required by this manifest.

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::decoy::{self, DecoyOptions};
use crate::group::{self, DetectionGroup};
use crate::hash;
use crate::prescreen::{self, KmerBloom};
use crate::reference::{self, Contig, ContigMeta, Role};
use crate::report::json_escape;

/// Manifest format version; `bloom.bin` has an independent binary version.
pub const FORMAT_VERSION: u32 = 2;
/// Version 1 remains supported only as a migration window for legacy ungrouped indexes. Closing
/// that window requires raising this value to 2 and removing `singleton_target_groups` and its test.
pub const MIN_SUPPORTED_FORMAT_VERSION: u32 = 1;
pub const MANIFEST_NAME: &str = "manifest.json";
pub const MMI_NAME: &str = "ref.mmi";
pub const BLOOM_NAME: &str = "bloom.bin";
pub const TARGETS_FA_NAME: &str = "targets.fa";
pub const DECOYS_FA_NAME: &str = "decoys.fa";
pub const DECOYS_TSV_NAME: &str = "decoys.tsv";

/// Index construction options shared by `viroflash index` and automatic construction in `run`.
#[derive(Debug, Clone)]
pub struct IndexOptions {
    pub host_fa: PathBuf,
    pub target_fa: PathBuf,
    pub contam_fa: Option<PathBuf>,
    pub decoy_fa: Option<PathBuf>,
    /// Integer ANI percentages used to generate decoys when `decoy_fa` is absent.
    pub decoy_anis: Vec<u8>,
    pub decoy_per_layer: usize,
    pub decoy_seed: u64,
    pub out_dir: PathBuf,
    pub k: usize,
    /// Number of threads used to build the minimap2 index.
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

/// Decoy provenance recorded in the manifest to keep the null distribution auditable.
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

/// Mapping of one target detection group into the composite reference.
/// Representatives and members retain their original target FASTA names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TargetGroupMeta {
    pub contig: String,
    pub representative: String,
    pub members: Vec<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct IndexManifest {
    pub format_version: u32,
    pub k: usize,
    pub created_at_unix: u64,
    pub refs: RefsInfo,
    pub contigs: Vec<ContigMeta>,
    pub target_groups: Vec<TargetGroupMeta>,
    pub bloom: BloomInfo,
}

/// Parsed index shared by the build and load paths, allowing `run` to use identical downstream
/// logic and produce equivalent results through either path.
#[derive(Debug)]
pub struct BuiltIndex {
    pub format_version: u32,
    pub mmi_path: PathBuf,
    pub bloom: KmerBloom,
    pub roles: HashMap<String, Role>,
    pub contigs: Vec<ContigMeta>,
    pub target_groups: Vec<TargetGroupMeta>,
    pub manifest_blake3: String,
}

fn now_unix() -> u64 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_secs())
        .unwrap_or(0)
}

/// Build an index directory. Existing output is rejected rather than overwritten. Construction
/// occurs in a temporary directory that is atomically renamed on success and cleaned on failure.
pub fn build_index(opt: &IndexOptions) -> Result<BuiltIndex, String> {
    validate_index_options(opt)?;
    let label = opt
        .out_dir
        .file_name()
        .map(|v| v.to_string_lossy().into_owned())
        .unwrap_or_else(|| "index".to_string());
    let monitor = crate::perf::PerfMonitor::start("index", label, &opt.out_dir, opt.threads)?;
    let result = build_index_with_monitor_validated(opt, &monitor);
    monitor.complete(result)
}

pub(crate) fn build_index_with_monitor(
    opt: &IndexOptions,
    monitor: &crate::perf::PerfMonitor,
) -> Result<BuiltIndex, String> {
    validate_index_options(opt)?;
    build_index_with_monitor_validated(opt, monitor)
}

fn build_index_with_monitor_validated(
    opt: &IndexOptions,
    monitor: &crate::perf::PerfMonitor,
) -> Result<BuiltIndex, String> {
    monitor.stage("index_prepare");
    if opt.out_dir.exists() {
        return Err(format!(
            "Index directory already exists: {}",
            opt.out_dir.display()
        ));
    }
    let part = PathBuf::from(format!(
        "{}.part.{}",
        opt.out_dir.display(),
        std::process::id()
    ));
    if part.exists() {
        // Remove only a stale temporary directory owned by the same process ID.
        let _ = std::fs::remove_dir_all(&part);
    }
    let result = build_index_into(&part, opt, monitor);
    if result.is_err() {
        let _ = std::fs::remove_dir_all(&part);
        return result;
    }
    if let Err(error) = std::fs::rename(&part, &opt.out_dir) {
        let _ = std::fs::remove_dir_all(&part);
        return Err(format!(
            "Failed to finalize index directory {}: {error}",
            opt.out_dir.display()
        ));
    }
    // The part directory is now final; returned artifact paths must point to its final location.
    let mut built = result?;
    built.mmi_path = opt.out_dir.join(MMI_NAME);
    Ok(built)
}

fn build_index_into(
    part: &Path,
    opt: &IndexOptions,
    monitor: &crate::perf::PerfMonitor,
) -> Result<BuiltIndex, String> {
    std::fs::create_dir_all(part)
        .map_err(|e| format!("Cannot create index directory {}: {e}", part.display()))?;

    // 1. Target groups: retain original sequences for Bloom insertion; write only representatives to the MMI.
    monitor.stage("index_target_groups");
    let target_records = reference::parse_fasta(&opt.target_fa)?;
    let mut original_target_names = HashSet::new();
    for (name, _) in &target_records {
        if !original_target_names.insert(name.as_str()) {
            return Err(format!(
                "Duplicate sequence name in target FASTA; cannot create unique detection-group membership: {name}"
            ));
        }
    }
    let detection_groups =
        group::construct_detection_groups(&target_records, monitor.work_thread_budget())
            .map_err(|e| format!("Failed to build target detection groups: {e}"))?;
    let targets_path = part.join(TARGETS_FA_NAME);
    let target_groups = write_grouped_targets(&targets_path, &target_records, &detection_groups)?;

    // 2. Decoys: use an explicit file or generate them from representatives with fixed parameters.
    //    Generated artifacts remain in the index so its null distribution is reproducible.
    monitor.stage("index_decoy");
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
                target_fa: targets_path.clone(),
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

    // 3. Composite reference and minimap2 index. Build below build/, move the MMI, then discard the composite.
    monitor.stage("index_reference_mmi");
    let build_dir = part.join("build");
    std::fs::create_dir_all(&build_dir)
        .map_err(|e| format!("Cannot create build directory {}: {e}", build_dir.display()))?;
    let mut fastas = vec![
        (Role::Host, opt.host_fa.clone()),
        (Role::Target, targets_path),
        (Role::Decoy, decoy_path.clone()),
    ];
    if let Some(c) = &opt.contam_fa {
        fastas.push((Role::Contaminant, c.clone()));
    }
    let (mmi_built, contigs) =
        reference::build_reference(&fastas, &build_dir, monitor.work_thread_budget().max(1))?;
    let mmi_path = part.join(MMI_NAME);
    std::fs::rename(&mmi_built, &mmi_path)
        .map_err(|e| format!("Failed to move index file: {e}"))?;
    let _ = std::fs::remove_file(build_dir.join("composite.fa"));
    let _ = std::fs::remove_dir(&build_dir);

    // 4. Bloom: insert every original target once plus actual decoys; do not reinsert representatives.
    monitor.stage("index_bloom");
    let original_targets: Vec<Contig> = target_records
        .into_iter()
        .map(|(name, seq)| Contig {
            name,
            role: Role::Target,
            gc_frac: reference::gc_fraction(&seq),
            seq,
        })
        .collect();
    let mut bloom_refs: Vec<&Contig> = original_targets.iter().collect();
    bloom_refs.extend(contigs.iter().filter(|c| c.role == Role::Decoy));
    let bloom = KmerBloom::build(&bloom_refs, opt.k, prescreen::GATE_FPR).ok_or_else(|| {
        format!(
            "The original target and decoy FASTA files contain no valid k-mers with k={} (all sequences are too short or contain non-ACGT characters)",
            opt.k
        )
    })?;
    write_bloom(&part.join(BLOOM_NAME), &bloom)?;

    // 5. Manifest: role-to-contig metadata, detection groups, source checksums, and decoy options.
    monitor.stage("index_manifest");
    let roles: HashMap<String, Role> = contigs.iter().map(|c| (c.name.clone(), c.role)).collect();
    let metas: Vec<ContigMeta> = contigs.iter().map(ContigMeta::from).collect();
    validate_contig_names(&metas)?;
    validate_target_groups(&metas, &target_groups)?;
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
        target_groups: target_groups.clone(),
        bloom: BloomInfo {
            n_inserted: bloom.n_inserted,
            fill_frac: bloom.fill_frac(),
        },
    };
    let text = write_manifest(&manifest);
    let manifest_path = part.join(MANIFEST_NAME);
    std::fs::write(&manifest_path, &text)
        .map_err(|e| format!("Failed to write manifest.json: {e}"))?;

    Ok(BuiltIndex {
        format_version: FORMAT_VERSION,
        mmi_path,
        bloom,
        roles,
        contigs: metas,
        target_groups,
        manifest_blake3: hash::blake3_hex(text.as_bytes()),
    })
}

fn write_grouped_targets(
    path: &Path,
    records: &[(String, Vec<u8>)],
    groups: &[DetectionGroup],
) -> Result<Vec<TargetGroupMeta>, String> {
    let file = File::create(path).map_err(|e| format!("Cannot create {}: {e}", path.display()))?;
    let mut writer = BufWriter::new(file);
    let mut metas = Vec::with_capacity(groups.len());
    for (group_index, detection_group) in groups.iter().enumerate() {
        let (representative, sequence) =
            records.get(detection_group.representative).ok_or_else(|| {
                format!("Representative index is out of bounds for detection group {group_index}")
            })?;
        writeln!(writer, ">{representative}")
            .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
        for chunk in sequence.chunks(60) {
            writer
                .write_all(chunk)
                .and_then(|_| writer.write_all(b"\n"))
                .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
        }

        let members = detection_group
            .members
            .iter()
            .map(|&member_index| {
                records
                    .get(member_index)
                    .map(|(name, _)| name.clone())
                    .ok_or_else(|| format!("Member index is out of bounds for detection group {group_index}: {member_index}"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        metas.push(TargetGroupMeta {
            contig: format!("target_{group_index}"),
            representative: representative.clone(),
            members,
        });
    }
    writer
        .flush()
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
    Ok(metas)
}

fn singleton_target_groups(contigs: &[ContigMeta]) -> Vec<TargetGroupMeta> {
    contigs
        .iter()
        .filter(|contig| contig.role == Role::Target)
        .map(|contig| TargetGroupMeta {
            contig: contig.name.clone(),
            representative: contig.name.clone(),
            members: vec![contig.name.clone()],
        })
        .collect()
}

fn validate_contig_names(contigs: &[ContigMeta]) -> Result<(), String> {
    let mut names = HashSet::new();
    for contig in contigs {
        if !names.insert(contig.name.as_str()) {
            return Err(format!(
                "Duplicate contig name in manifest.json contigs: {}",
                contig.name
            ));
        }
    }
    Ok(())
}

fn validate_target_groups(
    contigs: &[ContigMeta],
    target_groups: &[TargetGroupMeta],
) -> Result<(), String> {
    let target_contigs = contigs
        .iter()
        .filter(|contig| contig.role == Role::Target)
        .map(|contig| contig.name.as_str())
        .collect::<HashSet<_>>();

    let mut grouped_contigs = HashSet::new();
    let mut grouped_members = HashSet::new();
    for target_group in target_groups {
        if !grouped_contigs.insert(target_group.contig.as_str()) {
            return Err(format!(
                "Duplicate contig name in manifest.json target_groups: {}",
                target_group.contig
            ));
        }
        if !target_contigs.contains(target_group.contig.as_str()) {
            return Err(format!(
                "manifest.json target_groups contains unknown Target contig: {}",
                target_group.contig
            ));
        }
        if target_group.members.is_empty() {
            return Err(format!(
                "manifest.json target group {} has no members",
                target_group.contig
            ));
        }
        if !target_group
            .members
            .iter()
            .any(|member| member == &target_group.representative)
        {
            return Err(format!(
                "Representative is not a member of manifest.json target group {}: {}",
                target_group.contig, target_group.representative
            ));
        }
        let mut local_members = HashSet::new();
        for member in &target_group.members {
            if member.is_empty() {
                return Err(format!(
                    "manifest.json target group {} contains an empty member name",
                    target_group.contig
                ));
            }
            if !local_members.insert(member.as_str()) {
                return Err(format!(
                    "manifest.json target group {} contains duplicate member: {member}",
                    target_group.contig
                ));
            }
            if !grouped_members.insert(member.as_str()) {
                return Err(format!(
                    "Member appears in multiple manifest.json target groups: {member}"
                ));
            }
        }
    }

    let missing: Vec<_> = contigs
        .iter()
        .filter(|contig| {
            contig.role == Role::Target && !grouped_contigs.contains(contig.name.as_str())
        })
        .map(|contig| contig.name.as_str())
        .collect();
    if !missing.is_empty() {
        return Err(format!(
            "manifest.json target_groups is missing Target contig: {}",
            missing.join(", ")
        ));
    }
    Ok(())
}

fn file_info(p: &Path) -> Result<FileInfo, String> {
    Ok(FileInfo {
        path: p.display().to_string(),
        blake3: hash::blake3_file_hex(p)?,
    })
}

/// Load an index and validate existence, format version, k consistency, and Bloom integrity.
pub fn load_index(dir: &Path, k: usize) -> Result<BuiltIndex, String> {
    if !dir.is_dir() {
        return Err(format!("Index directory does not exist: {}", dir.display()));
    }
    let manifest_path = dir.join(MANIFEST_NAME);
    let text = std::fs::read_to_string(&manifest_path).map_err(|e| {
        format!(
            "Failed to read {} (not a viroflash index directory or manifest.json is missing): {e}",
            manifest_path.display()
        )
    })?;
    let manifest = parse_manifest(&text)?;
    if manifest.format_version < MIN_SUPPORTED_FORMAT_VERSION {
        return Err(format!(
            "Index format version {} is older than the supported version {}; rebuild the index with the current viroflash version",
            manifest.format_version, MIN_SUPPORTED_FORMAT_VERSION
        ));
    }
    if manifest.format_version > FORMAT_VERSION {
        return Err(format!(
            "Index format version {} is newer than the supported version {}; use a compatible viroflash version",
            manifest.format_version, FORMAT_VERSION
        ));
    }
    if manifest.k != k {
        return Err(format!(
            "Index k={} does not match --k={} (the Bloom filter and index must be rebuilt with the same k)",
            manifest.k, k
        ));
    }
    let mmi_path = dir.join(MMI_NAME);
    if !mmi_path.is_file() {
        return Err(format!("Index is missing {}", MMI_NAME));
    }
    let bloom = read_bloom(&dir.join(BLOOM_NAME))?;
    if bloom.k != manifest.k {
        return Err(format!(
            "{} and {} have mismatched k values ({} / {}); the index is corrupt and must be rebuilt",
            BLOOM_NAME, MANIFEST_NAME, bloom.k, manifest.k
        ));
    }
    let roles: HashMap<String, Role> = manifest
        .contigs
        .iter()
        .map(|c| (c.name.clone(), c.role))
        .collect();
    Ok(BuiltIndex {
        format_version: manifest.format_version,
        mmi_path,
        bloom,
        roles,
        contigs: manifest.contigs,
        target_groups: manifest.target_groups,
        manifest_blake3: hash::blake3_hex(text.as_bytes()),
    })
}

fn validate_index_options(opt: &IndexOptions) -> Result<(), String> {
    if !(1..=prescreen::K_MAX).contains(&opt.k) {
        return Err(format!(
            "--k must be between 1 and {} (the 2-bit encoding limit); got {}",
            prescreen::K_MAX,
            opt.k
        ));
    }
    if opt.threads == 0 {
        return Err("--threads must be greater than 0".into());
    }
    if opt.host_fa.as_os_str().is_empty() || opt.target_fa.as_os_str().is_empty() {
        return Err("Building an index requires --host-fa and --target-fa".into());
    }
    if opt.out_dir.as_os_str().is_empty() {
        return Err("Building an index requires --out".into());
    }
    if opt.decoy_fa.is_none() {
        if opt.decoy_anis.is_empty() {
            return Err("--decoy-ani requires at least one layer".into());
        }
        for &a in &opt.decoy_anis {
            if a == 0 || a >= 100 {
                return Err(format!(
                    "Invalid --decoy-ani layer {a} (must be between 1 and 99)"
                ));
            }
        }
        if opt.decoy_per_layer == 0 {
            return Err("--decoy-per-layer must be greater than 0".into());
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Bloom binary serialization: private, versioned, magic-tagged, and cross-platform deterministic.
// ---------------------------------------------------------------------------

const BLOOM_MAGIC: [u8; 4] = *b"VFB1";
const BLOOM_VERSION: u32 = 1;
const BLOOM_IO_BUFFER_BYTES: usize = 8 * 1024 * 1024;

fn write_bloom(path: &Path, bloom: &KmerBloom) -> Result<(), String> {
    let (words, _mask) = bloom.to_raw();
    let file =
        File::create(path).map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
    let mut writer = BufWriter::new(file);
    writer
        .write_all(&BLOOM_MAGIC)
        .and_then(|_| writer.write_all(&BLOOM_VERSION.to_le_bytes()))
        .and_then(|_| writer.write_all(&(bloom.k as u64).to_le_bytes()))
        .and_then(|_| writer.write_all(&bloom.n_inserted.to_le_bytes()))
        .and_then(|_| writer.write_all(&(words.len() as u64).to_le_bytes()))
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;

    let words_per_chunk = BLOOM_IO_BUFFER_BYTES / 8;
    let mut bytes = Vec::with_capacity(words.len().min(words_per_chunk) * 8);
    for chunk in words.chunks(words_per_chunk) {
        bytes.clear();
        for word in chunk {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        writer
            .write_all(&bytes)
            .map_err(|e| format!("Failed to write {}: {e}", path.display()))?;
    }
    writer
        .flush()
        .map_err(|e| format!("Failed to write {}: {e}", path.display()))
}

fn read_bloom(path: &Path) -> Result<KmerBloom, String> {
    const HEADER_LEN: usize = 4 + 4 + 8 + 8 + 8;
    let mut file = File::open(path)
        .map_err(|e| format!("Index is missing or unreadable {}: {e}", path.display()))?;
    let file_len = file
        .metadata()
        .map_err(|e| format!("Failed to read metadata for {}: {e}", path.display()))?
        .len();
    if file_len < HEADER_LEN as u64 {
        return Err(format!("{} is corrupt (too short)", path.display()));
    }
    let mut header = [0u8; HEADER_LEN];
    file.read_exact(&mut header)
        .map_err(|e| format!("Failed to read {}: {e}", path.display()))?;
    if header[..4] != BLOOM_MAGIC {
        return Err(format!(
            "{} has an invalid magic number (incompatible version or corrupt file)",
            path.display()
        ));
    }
    let version = u32::from_le_bytes(header[4..8].try_into().unwrap());
    if version > BLOOM_VERSION {
        return Err(format!(
            "{} version {version} is newer than supported version {BLOOM_VERSION}; rebuild the index",
            path.display()
        ));
    }
    let k = u64::from_le_bytes(header[8..16].try_into().unwrap()) as usize;
    if !(1..=prescreen::K_MAX).contains(&k) {
        return Err(format!("{} is corrupt (invalid k: {k})", path.display()));
    }
    let n_inserted = u64::from_le_bytes(header[16..24].try_into().unwrap());
    let words_len_u64 = u64::from_le_bytes(header[24..32].try_into().unwrap());
    let words_len = usize::try_from(words_len_u64).map_err(|_| {
        format!(
            "{} is corrupt (bit capacity exceeds the platform limit)",
            path.display()
        )
    })?;
    if words_len == 0 || !words_len.is_power_of_two() {
        return Err(format!(
            "{} is corrupt (bit capacity is not a power of two: {words_len})",
            path.display()
        ));
    }
    let expected_len = words_len_u64
        .checked_mul(8)
        .and_then(|n| n.checked_add(HEADER_LEN as u64))
        .ok_or_else(|| format!("{} is corrupt (length overflow)", path.display()))?;
    if file_len != expected_len {
        return Err(format!(
            "{} is corrupt (length mismatch: expected {}, got {})",
            path.display(),
            expected_len,
            file_len
        ));
    }

    let mut words = Vec::new();
    words.try_reserve_exact(words_len).map_err(|e| {
        format!(
            "Failed to allocate memory while loading {}: {e}",
            path.display()
        )
    })?;
    let chunk_words_capacity = words_len.min(BLOOM_IO_BUFFER_BYTES / 8);
    let mut bytes = vec![0u8; chunk_words_capacity * 8];
    while words.len() < words_len {
        let remaining = words_len - words.len();
        let chunk_words = remaining.min(bytes.len() / 8);
        let chunk_bytes = chunk_words * 8;
        file.read_exact(&mut bytes[..chunk_bytes])
            .map_err(|e| format!("Failed to read {}: {e}", path.display()))?;
        for raw in bytes[..chunk_bytes].chunks_exact(8) {
            words.push(u64::from_le_bytes(raw.try_into().unwrap()));
        }
    }
    Ok(KmerBloom::from_raw(
        words,
        (words_len * 64 - 1) as u64,
        k,
        n_inserted,
    ))
}

// ---------------------------------------------------------------------------
// Dependency-free manifest JSON serialization and parsing.
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
        // Use f64's shortest round-trip representation so parsed GC values remain bit-identical.
        // Six-digit truncation could otherwise change strata at 0.40/0.50/0.60 boundaries.
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
    s.push_str("  \"target_groups\": [\n");
    for (group_index, group) in m.target_groups.iter().enumerate() {
        s.push_str(&format!(
            "    {{\"contig\": \"{}\", \"representative\": \"{}\", \"members\": [",
            json_escape(&group.contig),
            json_escape(&group.representative)
        ));
        for (member_index, member) in group.members.iter().enumerate() {
            if member_index != 0 {
                s.push_str(", ");
            }
            s.push_str(&format!("\"{}\"", json_escape(member)));
        }
        s.push_str(&format!(
            "]}}{}\n",
            if group_index + 1 < m.target_groups.len() {
                ","
            } else {
                ""
            }
        ));
    }
    s.push_str("  ],\n");
    s.push_str(&format!(
        "  \"files\": {{\"mmi\": \"{}\", \"bloom\": \"{}\", \"targets\": \"{}\"}},\n",
        MMI_NAME, BLOOM_NAME, TARGETS_FA_NAME
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

// --- Minimal JSON parser for the subset required by this manifest ---

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
        format!(
            "Failed to parse manifest.json at offset {}: {msg}",
            self.pos
        )
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
            .ok_or_else(|| self.err("Unexpected end of input"))
    }

    fn expect(&mut self, b: u8) -> Result<(), String> {
        if self.peek()? == b {
            self.pos += 1;
            Ok(())
        } else {
            Err(self.err(&format!("Expected '{}'", b as char)))
        }
    }

    fn parse(&mut self) -> Result<Json, String> {
        let v = self.parse_value()?;
        self.ws();
        if self.pos != self.bytes.len() {
            return Err(self.err("Trailing content"));
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
            other => Err(self.err(&format!("Invalid character '{}'", other as char))),
        }
    }

    fn expect_lit(&mut self, lit: &[u8]) -> Result<(), String> {
        for &b in lit {
            if self.peek()? == b {
                self.pos += 1;
            } else {
                return Err(self.err("Invalid literal"));
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
                    return Err(
                        self.err(&format!("Expected ','  or  '}}'; got '{}'", other as char))
                    )
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
                    return Err(self.err(&format!("Expected ','  or  ']'; got '{}'", other as char)))
                }
            }
        }
        Ok(Json::Arr(items))
    }

    fn parse_string(&mut self) -> Result<String, String> {
        self.expect(b'"')?;
        let mut out = String::new();
        loop {
            let byte_pos = self.pos;
            let b = *self
                .bytes
                .get(self.pos)
                .ok_or_else(|| self.err("Unterminated string"))?;
            self.pos += 1;
            match b {
                b'"' => break,
                b'\\' => {
                    let e = *self
                        .bytes
                        .get(self.pos)
                        .ok_or_else(|| self.err("Unterminated escape sequence"))?;
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
                        other => {
                            return Err(
                                self.err(&format!("Invalid escape sequence '\\{}'", other as char))
                            )
                        }
                    }
                }
                0x20..=0x7e => out.push(b as char),
                0x80..=0xff => {
                    let tail = std::str::from_utf8(&self.bytes[byte_pos..])
                        .map_err(|_| self.err("String contains invalid UTF-8"))?;
                    let character = tail
                        .chars()
                        .next()
                        .ok_or_else(|| self.err("Unterminated string"))?;
                    self.pos = byte_pos + character.len_utf8();
                    out.push(character);
                }
                other => {
                    return Err(self.err(&format!(
                        "Invalid control character in string 0x{other:02x}"
                    )))
                }
            }
        }
        Ok(out)
    }

    /// Parse `\uXXXX`, including surrogate pairs.
    fn parse_unicode_escape(&mut self) -> Result<char, String> {
        let hi = self.parse_hex4()?;
        let cp = if (0xD800..=0xDBFF).contains(&hi) {
            // A high surrogate must be followed by \uDC00-\uDFFF.
            if self.bytes.get(self.pos) == Some(&b'\\')
                && self.bytes.get(self.pos + 1) == Some(&b'u')
            {
                self.pos += 2;
                let lo = self.parse_hex4()?;
                if !(0xDC00..=0xDFFF).contains(&lo) {
                    return Err(self.err("Invalid low surrogate"));
                }
                0x10000 + ((hi - 0xD800) << 10) + (lo - 0xDC00)
            } else {
                return Err(self.err("High surrogate is not followed by a low surrogate"));
            }
        } else {
            hi
        };
        char::from_u32(cp).ok_or_else(|| self.err("Invalid Unicode code point"))
    }

    fn parse_hex4(&mut self) -> Result<u32, String> {
        let mut v = 0u32;
        for _ in 0..4 {
            let b = *self
                .bytes
                .get(self.pos)
                .ok_or_else(|| self.err("Truncated \\u escape"))?;
            self.pos += 1;
            v = v * 16
                + match b {
                    b'0'..=b'9' => u32::from(b - b'0'),
                    b'a'..=b'f' => u32::from(b - b'a' + 10),
                    b'A'..=b'F' => u32::from(b - b'A' + 10),
                    _ => return Err(self.err("Invalid hexadecimal digit in \\u escape")),
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
        let text = std::str::from_utf8(&self.bytes[start..self.pos])
            .map_err(|_| self.err("Invalid number"))?;
        text.parse::<f64>()
            .map(Json::Num)
            .map_err(|_| self.err(&format!("Invalid number: {text}")))
    }
}

// --- Extract manifest fields from the parsed representation ---

fn obj_get<'a>(obj: &'a [(String, Json)], key: &str) -> Option<&'a Json> {
    obj.iter().find(|(k, _)| k == key).map(|(_, v)| v)
}

fn as_obj<'a>(v: &'a Json, what: &str) -> Result<&'a [(String, Json)], String> {
    match v {
        Json::Obj(pairs) => Ok(pairs),
        _ => Err(format!("manifest.json field {what} must be an object")),
    }
}

fn as_str<'a>(v: &'a Json, what: &str) -> Result<&'a str, String> {
    match v {
        Json::Str(s) => Ok(s),
        _ => Err(format!("manifest.json field {what} must be a string")),
    }
}

fn as_u64(v: &Json, what: &str) -> Result<u64, String> {
    match v {
        Json::Num(n) if *n >= 0.0 && n.fract() == 0.0 && *n <= u64::MAX as f64 => Ok(*n as u64),
        _ => Err(format!(
            "manifest.json field {what} must be a non-negative integer"
        )),
    }
}

fn as_f64(v: &Json, what: &str) -> Result<f64, String> {
    match v {
        Json::Num(n) => Ok(*n),
        _ => Err(format!("manifest.json field {what} must be a number")),
    }
}

fn role_from_str(s: &str) -> Result<Role, String> {
    match s {
        "host" => Ok(Role::Host),
        "target" => Ok(Role::Target),
        "decoy" => Ok(Role::Decoy),
        "contam" => Ok(Role::Contaminant),
        other => Err(format!("Unknown role in manifest.json: {other}")),
    }
}

fn parse_manifest(text: &str) -> Result<IndexManifest, String> {
    let root = JsonParser::new(text).parse()?;
    let top = as_obj(&root, "root")?;
    let format = as_str(
        obj_get(top, "format").ok_or("manifest.json is missing field format")?,
        "format",
    )?;
    if format != "viroflash.index" {
        return Err(format!(
            "manifest.json format must be \"viroflash.index\"; got \"{format}\""
        ));
    }
    let format_version = as_u64(
        obj_get(top, "version").ok_or("manifest.json is missing field version")?,
        "version",
    )? as u32;
    let k = as_u64(
        obj_get(top, "k").ok_or("manifest.json is missing field k")?,
        "k",
    )? as usize;
    if !(1..=prescreen::K_MAX).contains(&k) {
        return Err(format!("Invalid k in manifest.json: {k}"));
    }
    let created_at_unix = as_u64(
        obj_get(top, "created_at_unix").ok_or("manifest.json is missing field created_at_unix")?,
        "created_at_unix",
    )?;
    let refs_root = as_obj(
        obj_get(top, "refs").ok_or("manifest.json is missing field refs")?,
        "refs",
    )?;
    let parse_file = |key: &str| -> Result<FileInfo, String> {
        let o = as_obj(
            obj_get(refs_root, key)
                .ok_or_else(|| format!("manifest.json refs is missing {key}"))?,
            key,
        )?;
        Ok(FileInfo {
            path: as_str(
                obj_get(o, "path").ok_or_else(|| format!("refs.{key} is missing path"))?,
                "path",
            )?
            .to_string(),
            blake3: as_str(
                obj_get(o, "blake3").ok_or_else(|| format!("refs.{key} is missing blake3"))?,
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
        obj_get(refs_root, "decoy").ok_or("manifest.json refs is missing decoy")?,
        "decoy",
    )?;
    let decoy_source = if let Some(v) = obj_get(decoy, "file") {
        let o = as_obj(v, "decoy.file")?;
        DecoySource::File {
            path: as_str(
                obj_get(o, "path").ok_or("decoy.file is missing path")?,
                "path",
            )?
            .to_string(),
            blake3: as_str(
                obj_get(o, "blake3").ok_or("decoy.file is missing blake3")?,
                "blake3",
            )?
            .to_string(),
        }
    } else if let Some(v) = obj_get(decoy, "generated") {
        let o = as_obj(v, "decoy.generated")?;
        let anis = match obj_get(o, "anis").ok_or("decoy.generated is missing anis")? {
            Json::Arr(items) => items
                .iter()
                .map(|i| as_u64(i, "anis").map(|n| n as u8))
                .collect::<Result<Vec<u8>, _>>()?,
            _ => return Err("decoy.generated.anis must be an array".into()),
        };
        DecoySource::Generated {
            anis,
            per_layer: as_u64(
                obj_get(o, "per_layer").ok_or("decoy.generated is missing per_layer")?,
                "per_layer",
            )? as usize,
            seed: as_u64(
                obj_get(o, "seed").ok_or("decoy.generated is missing seed")?,
                "seed",
            )?,
        }
    } else {
        return Err("manifest.json refs.decoy must be file or generated".into());
    };
    let contigs = match obj_get(top, "contigs").ok_or("manifest.json is missing field contigs")? {
        Json::Arr(items) => items
            .iter()
            .map(|v| {
                let o = as_obj(v, "contigs[]")?;
                Ok(ContigMeta {
                    name: as_str(obj_get(o, "name").ok_or("contig is missing name")?, "name")?
                        .to_string(),
                    role: role_from_str(as_str(
                        obj_get(o, "role").ok_or("contig is missing role")?,
                        "role",
                    )?)?,
                    len: as_u64(obj_get(o, "len").ok_or("contig is missing len")?, "len")?,
                    gc_frac: as_f64(obj_get(o, "gc").ok_or("contig is missing gc")?, "gc")?,
                })
            })
            .collect::<Result<Vec<ContigMeta>, String>>()?,
        _ => return Err("manifest.json contigs must be an array".into()),
    };
    validate_contig_names(&contigs)?;
    let target_groups = match obj_get(top, "target_groups") {
        None if format_version == 1 => singleton_target_groups(&contigs),
        None => {
            return Err(format!(
                "manifest.json format version {format_version} is missing target_groups"
            ));
        }
        Some(Json::Arr(items)) => items
            .iter()
            .map(|value| {
                let group = as_obj(value, "target_groups[]")?;
                let members = match obj_get(group, "members")
                    .ok_or("target_group is missing members")?
                {
                    Json::Arr(items) => items
                        .iter()
                        .map(|item| as_str(item, "target_group.members[]").map(str::to_string))
                        .collect::<Result<Vec<_>, _>>()?,
                    _ => return Err("manifest.json target_group.members must be an array".into()),
                };
                Ok(TargetGroupMeta {
                    contig: as_str(
                        obj_get(group, "contig").ok_or("target_group is missing contig")?,
                        "target_group.contig",
                    )?
                    .to_string(),
                    representative: as_str(
                        obj_get(group, "representative")
                            .ok_or("target_group is missing representative")?,
                        "target_group.representative",
                    )?
                    .to_string(),
                    members,
                })
            })
            .collect::<Result<Vec<_>, String>>()?,
        Some(_) => return Err("manifest.json target_groups must be an array".into()),
    };
    if (2..=FORMAT_VERSION).contains(&format_version) {
        validate_target_groups(&contigs, &target_groups)?;
    }
    let bloom_root = as_obj(
        obj_get(top, "bloom").ok_or("manifest.json is missing field bloom")?,
        "bloom",
    )?;
    let bloom = BloomInfo {
        n_inserted: as_u64(
            obj_get(bloom_root, "n_inserted").ok_or("bloom is missing n_inserted")?,
            "n_inserted",
        )?,
        fill_frac: as_f64(
            obj_get(bloom_root, "fill_frac").ok_or("bloom is missing fill_frac")?,
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
        target_groups,
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

    fn write_records(path: &Path, records: &[(&str, &[u8])]) {
        let mut out = String::new();
        for (name, sequence) in records {
            out.push('>');
            out.push_str(name);
            out.push('\n');
            out.push_str(std::str::from_utf8(sequence).unwrap());
            out.push('\n');
        }
        std::fs::write(path, out).unwrap();
    }

    fn without_target_groups(mut text: String) -> String {
        let groups_start = text.find("  \"target_groups\": [\n").unwrap();
        let files_start = text[groups_start..]
            .find("  \"files\":")
            .map(|offset| groups_start + offset)
            .unwrap();
        text.replace_range(groups_start..files_start, "");
        text
    }

    fn pseudo_dna(len: usize, seed: u64) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = prescreen::splitmix64(state);
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn mutate_substitutions(sequence: &[u8], count: usize) -> Vec<u8> {
        let mut result = sequence.to_vec();
        let mut selected = vec![false; result.len()];
        let mut state = 0x8d26_1c47_b1e9_53f1_u64;
        let mut changed = 0_usize;
        while changed < count.min(result.len()) {
            state = prescreen::splitmix64(state);
            let position = state as usize % result.len();
            if selected[position] {
                continue;
            }
            selected[position] = true;
            result[position] = match result[position] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                _ => b'A',
            };
            changed += 1;
        }
        result
    }

    fn canonical_kmer(sequence: &[u8]) -> Option<u64> {
        let forward = prescreen::encode_kmer(sequence)?;
        let reverse = prescreen::encode_kmer(&prescreen::reverse_complement(sequence))?;
        Some(forward.min(reverse))
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
        assert!(opt.out_dir.join(TARGETS_FA_NAME).is_file());

        let loaded = load_index(&opt.out_dir, opt.k).unwrap();
        assert_eq!(loaded.format_version, FORMAT_VERSION);
        assert_eq!(loaded.contigs, built.contigs);
        assert_eq!(loaded.roles, built.roles);
        assert_eq!(loaded.target_groups, built.target_groups);
        assert_eq!(loaded.manifest_blake3, built.manifest_blake3);
        assert_eq!(loaded.bloom.k, built.bloom.k);
        assert_eq!(loaded.bloom.n_inserted, built.bloom.n_inserted);
        // One contig for each role: host, target, decoy, and contaminant.
        let counts = |c: &[ContigMeta], r: Role| c.iter().filter(|x| x.role == r).count();
        assert_eq!(counts(&loaded.contigs, Role::Host), 1);
        assert_eq!(counts(&loaded.contigs, Role::Target), 1);
        assert_eq!(counts(&loaded.contigs, Role::Decoy), 1);
        assert_eq!(counts(&loaded.contigs, Role::Contaminant), 1);
        let manifest_text = std::fs::read_to_string(opt.out_dir.join(MANIFEST_NAME)).unwrap();
        assert!(manifest_text.contains("\"targets\": \"targets.fa\""));
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn exact_duplicate_targets_fold_and_drive_auto_decoy() {
        let dir = tmp_dir("exact_groups");
        let mut opt = default_index_opts(&dir);
        let sequence = pseudo_dna(3_000, 11);
        write_records(
            &opt.target_fa,
            &[
                ("duplicate-a", sequence.as_slice()),
                ("duplicate-b", sequence.as_slice()),
            ],
        );
        opt.decoy_fa = None;
        opt.decoy_anis = vec![82];
        opt.decoy_per_layer = 1;

        let built = build_index(&opt).unwrap();
        assert_eq!(
            built
                .contigs
                .iter()
                .filter(|contig| contig.role == Role::Target)
                .count(),
            1
        );
        assert_eq!(
            built
                .contigs
                .iter()
                .filter(|contig| contig.role == Role::Decoy)
                .count(),
            1,
            "generated decoys must use the collapsed targets.fa"
        );
        assert_eq!(
            built.target_groups,
            vec![TargetGroupMeta {
                contig: "target_0".into(),
                representative: "duplicate-a".into(),
                members: vec!["duplicate-a".into(), "duplicate-b".into()],
            }]
        );
        assert_eq!(
            reference::parse_fasta(&opt.out_dir.join(TARGETS_FA_NAME)).unwrap(),
            vec![("duplicate-a".into(), sequence)]
        );
        let loaded = load_index(&opt.out_dir, opt.k).unwrap();
        assert_eq!(loaded.target_groups, built.target_groups);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn near_group_metadata_and_folded_member_kmers_survive_bloom() {
        let dir = tmp_dir("near_groups");
        let mut opt = default_index_opts(&dir);
        opt.contam_fa = None;
        let original = pseudo_dna(3_000, 17);
        let variant = mutate_substitutions(&original, 60);
        write_records(
            &opt.target_fa,
            &[
                ("original", original.as_slice()),
                ("variant", variant.as_slice()),
            ],
        );

        let built = build_index(&opt).unwrap();
        assert_eq!(built.target_groups.len(), 1);
        let target_group = &built.target_groups[0];
        assert_eq!(target_group.contig, "target_0");
        assert_eq!(target_group.members, ["original", "variant"]);

        let (representative_sequence, folded_sequence) = match target_group.representative.as_str()
        {
            "original" => (original.as_slice(), variant.as_slice()),
            "variant" => (variant.as_slice(), original.as_slice()),
            other => panic!("unknown representative sequence: {other}"),
        };
        let representative_kmers: HashSet<_> = representative_sequence
            .windows(opt.k)
            .filter_map(canonical_kmer)
            .collect();
        let private_kmer = folded_sequence
            .windows(opt.k)
            .filter_map(canonical_kmer)
            .find(|code| !representative_kmers.contains(code))
            .expect(
                "a collapsed member should have a private k-mer absent from its representative",
            );
        assert!(
            built.bloom.probe(private_kmer),
            "private k-mers from collapsed non-representatives must not produce Bloom false negatives"
        );

        let decoy_kmers: usize = reference::parse_fasta(opt.decoy_fa.as_deref().unwrap())
            .unwrap()
            .iter()
            .map(|(_, sequence)| sequence.len().saturating_sub(opt.k - 1))
            .sum();
        let expected_insertions = 2 * (original.len() - opt.k + 1) + decoy_kmers;
        assert_eq!(
            built.bloom.n_inserted, expected_insertions as u64,
            "the Bloom filter must contain each original target once and each actual decoy, without recounting representatives"
        );

        let representatives = reference::parse_fasta(&opt.out_dir.join(TARGETS_FA_NAME)).unwrap();
        assert_eq!(representatives.len(), 1);
        assert_eq!(representatives[0].0, target_group.representative);
        assert_eq!(representatives[0].1, representative_sequence);
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
                    name: "target_0".into(),
                    role: Role::Target,
                    len: 3_000,
                    gc_frac: 0.51,
                },
                ContigMeta {
                    name: "decoy:tv:ani82:i0".into(),
                    role: Role::Decoy,
                    len: 293116,
                    gc_frac: 0.52,
                },
            ],
            target_groups: vec![TargetGroupMeta {
                contig: "target_0".into(),
                representative: "virus-a".into(),
                members: vec!["virus-a".into(), "virus-b".into()],
            }],
            bloom: BloomInfo {
                n_inserted: 42,
                fill_frac: 0.4123,
            },
        };
        let parsed = parse_manifest(&write_manifest(&m)).unwrap();
        assert_eq!(parsed, m);
        // Round-trip string escapes in paths containing quotes, backslashes, and control characters.
        let mut m2 = m.clone();
        m2.refs.host.path = "a\"b\\c\nd".into();
        m2.target_groups[0].members[1] = "b\"c\\d\ne".into();
        let parsed2 = parse_manifest(&write_manifest(&m2)).unwrap();
        assert_eq!(parsed2, m2);
    }

    #[test]
    fn manifest_v2_rejects_invalid_target_group_contigs() {
        let mut manifest = IndexManifest {
            format_version: FORMAT_VERSION,
            k: 21,
            created_at_unix: 1,
            refs: RefsInfo {
                host: FileInfo {
                    path: "host.fa".into(),
                    blake3: "a".repeat(64),
                },
                target: FileInfo {
                    path: "target.fa".into(),
                    blake3: "b".repeat(64),
                },
                contam: None,
                decoy: DecoySource::File {
                    path: "decoy.fa".into(),
                    blake3: "c".repeat(64),
                },
            },
            contigs: vec![ContigMeta {
                name: "target_0".into(),
                role: Role::Target,
                len: 100,
                gc_frac: 0.5,
            }],
            target_groups: vec![TargetGroupMeta {
                contig: "target_0".into(),
                representative: "virus".into(),
                members: vec!["virus".into()],
            }],
            bloom: BloomInfo {
                n_inserted: 80,
                fill_frac: 0.1,
            },
        };

        let error = parse_manifest(&without_target_groups(write_manifest(&manifest))).unwrap_err();
        assert!(error.contains("missing target_groups"), "error={error}");

        let mut duplicate_non_target = manifest.clone();
        let host = ContigMeta {
            name: "host_0".into(),
            role: Role::Host,
            len: 100,
            gc_frac: 0.5,
        };
        duplicate_non_target.contigs.push(host.clone());
        duplicate_non_target.contigs.push(host);
        let error = parse_manifest(&write_manifest(&duplicate_non_target)).unwrap_err();
        assert!(error.contains("Duplicate contig"), "error={error}");

        manifest
            .target_groups
            .push(manifest.target_groups[0].clone());
        let error = parse_manifest(&write_manifest(&manifest)).unwrap_err();
        assert!(error.contains("Duplicate"), "error={error}");

        manifest.target_groups = vec![TargetGroupMeta {
            contig: "target_1".into(),
            representative: "virus".into(),
            members: vec!["virus".into()],
        }];
        let error = parse_manifest(&write_manifest(&manifest)).unwrap_err();
        assert!(error.contains("unknown"), "error={error}");

        manifest.target_groups.clear();
        let error = parse_manifest(&write_manifest(&manifest)).unwrap_err();
        assert!(error.contains("missing"), "error={error}");

        manifest.target_groups = vec![TargetGroupMeta {
            contig: "target_0".into(),
            representative: "not-a-member".into(),
            members: vec!["virus".into()],
        }];
        let error = parse_manifest(&write_manifest(&manifest)).unwrap_err();
        assert!(
            error.contains("Representative is not a member"),
            "error={error}"
        );

        manifest.target_groups[0] = TargetGroupMeta {
            contig: "target_0".into(),
            representative: "virus".into(),
            members: vec!["virus".into(), "virus".into()],
        };
        let error = parse_manifest(&write_manifest(&manifest)).unwrap_err();
        assert!(error.contains("duplicate member"), "error={error}");
    }

    #[test]
    fn load_v1_manifest_without_groups_synthesizes_singletons() {
        let dir = tmp_dir("v1_fallback");
        let opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        let manifest_path = opt.out_dir.join(MANIFEST_NAME);
        let text = std::fs::read_to_string(&manifest_path).unwrap().replacen(
            &format!("\"version\": {FORMAT_VERSION}"),
            "\"version\": 1",
            1,
        );
        let mut text = without_target_groups(text);
        text = text.replace(", \"targets\": \"targets.fa\"", "");
        std::fs::write(&manifest_path, text).unwrap();

        let loaded = load_index(&opt.out_dir, opt.k).unwrap();
        assert_eq!(loaded.format_version, 1);
        assert_eq!(
            loaded.target_groups,
            vec![TargetGroupMeta {
                contig: "target_0".into(),
                representative: "target_0".into(),
                members: vec!["target_0".into()],
            }]
        );
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_rejects_k_mismatch_and_version() {
        let dir = tmp_dir("mismatch");
        let opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        let err = load_index(&opt.out_dir, opt.k + 1).unwrap_err();
        assert!(err.contains("does not match"), "err={err}");
        // Reject versions outside the supported range after manually changing the version field.
        let mp = opt.out_dir.join(MANIFEST_NAME);
        let original = std::fs::read_to_string(&mp).unwrap();
        let text = original.replacen(
            &format!("\"version\": {FORMAT_VERSION}"),
            "\"version\": 99",
            1,
        );
        std::fs::write(&mp, text).unwrap();
        let err = load_index(&opt.out_dir, opt.k).unwrap_err();
        assert!(err.contains("newer"), "err={err}");

        let text = original.replacen(
            &format!("\"version\": {FORMAT_VERSION}"),
            "\"version\": 0",
            1,
        );
        std::fs::write(&mp, text).unwrap();
        let err = load_index(&opt.out_dir, opt.k).unwrap_err();
        assert!(err.contains("older"), "err={err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn load_rejects_missing_and_corrupt() {
        let dir = tmp_dir("corrupt");
        let opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        // Missing manifest in an empty directory.
        let empty = dir.join("empty");
        std::fs::create_dir_all(&empty).unwrap();
        let err = load_index(&empty, opt.k).unwrap_err();
        assert!(err.contains("manifest.json"), "err={err}");
        // Truncated Bloom file.
        let bp = opt.out_dir.join(BLOOM_NAME);
        let bytes = std::fs::read(&bp).unwrap();
        std::fs::write(&bp, &bytes[..bytes.len() / 2]).unwrap();
        let err = load_index(&opt.out_dir, opt.k).unwrap_err();
        assert!(err.contains("corrupt"), "err={err}");
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn build_rejects_existing_dir_and_bad_params() {
        let dir = tmp_dir("exist");
        let mut opt = default_index_opts(&dir);
        build_index(&opt).unwrap();
        let err = build_index(&opt).unwrap_err();
        assert!(err.contains("already exists"), "err={err}");
        // Failed construction must clean its temporary directory.
        opt.out_dir = dir.join("idx2");
        opt.k = 0;
        let err = build_index(&opt).unwrap_err();
        assert!(err.contains("--k"), "err={err}");
        let leftovers: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().starts_with("idx2"))
            .collect();
        assert!(
            leftovers.is_empty(),
            "failed builds must not leave temporary directories"
        );
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
        // host + target + 2 layers x 2 = 6 contigs; contaminant is omitted.
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
