use std::collections::{HashMap, HashSet};
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use minimap2::Aligner;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

use crate::analysis_profile::{hex_sha256, AnalysisProfile};
use crate::kmer_gate::{BloomSummary, TargetKmerBloom};
use crate::reference_group::{build_reference_groups, parse_fasta, ReferenceGroup};

const INDEX_CONTRACT: &str = "viroflash.reference-index.profile-bound";
const MANIFEST: &str = "manifest.json";
const MMI: &str = "ref.mmi";
const BLOOM: &str = "bloom.bin";
const LEDGER: &str = "reference-groups.tsv";
const COMPOSITE: &str = "reference.fa";
const SINGLE_PART_BATCH_SIZE: u64 = u64::MAX;

#[derive(Debug, Clone)]
pub struct IndexOptions {
    pub host_fa: PathBuf,
    pub target_fa: PathBuf,
    pub out_dir: PathBuf,
    pub threads: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum ReferenceRole {
    Host,
    Target,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct ReferenceContig {
    pub name: String,
    pub role: ReferenceRole,
    pub target_group_ordinal: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct IndexManifest {
    contract_id: String,
    profile_digest: String,
    reference_set_digest: String,
    ledger_digest: String,
    mmi_digest: String,
    bloom_digest: String,
    kmer_length: usize,
    host_fasta_sha256: String,
    target_fasta_sha256: String,
    bloom: BloomSummary,
    contigs: Vec<ReferenceContig>,
}

pub struct ReferenceIndex {
    pub mmi_path: PathBuf,
    pub bloom: TargetKmerBloom,
    pub bloom_summary: BloomSummary,
    pub contigs: HashMap<String, ReferenceContig>,
    pub target_groups: Vec<ReferenceGroup>,
    pub profile_digest: String,
    pub index_digest: String,
}

pub fn build_index(options: &IndexOptions) -> Result<ReferenceIndex, String> {
    if options.threads == 0 {
        return Err("--threads must be greater than zero".into());
    }
    if options.out_dir.exists() {
        return Err(format!(
            "Index directory already exists: {}",
            options.out_dir.display()
        ));
    }
    let temporary = PathBuf::from(format!(
        "{}.part.{}",
        options.out_dir.display(),
        std::process::id()
    ));
    if temporary.exists() {
        return Err(format!(
            "Index staging path already exists: {}",
            temporary.display()
        ));
    }
    std::fs::create_dir_all(&temporary)
        .map_err(|error| format!("Cannot create {}: {error}", temporary.display()))?;
    if let Err(error) = build_into(&temporary, options) {
        let _ = std::fs::remove_dir_all(&temporary);
        return Err(error);
    }
    std::fs::rename(&temporary, &options.out_dir).map_err(|error| {
        format!(
            "Cannot finalize index {}: {error}",
            options.out_dir.display()
        )
    })?;
    load_index(&options.out_dir)
}

fn build_into(directory: &Path, options: &IndexOptions) -> Result<(), String> {
    let profile = AnalysisProfile::FROZEN;
    let host_records = parse_fasta(&options.host_fa)?;
    let target_records = parse_fasta(&options.target_fa)?;
    let groups = build_reference_groups(&target_records);
    let by_id = target_records
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<HashMap<_, _>>();
    let mut contigs = Vec::new();
    let composite = directory.join(COMPOSITE);
    let mut writer = BufWriter::new(
        File::create(&composite)
            .map_err(|error| format!("Cannot create {}: {error}", composite.display()))?,
    );
    for (ordinal, record) in host_records.iter().enumerate() {
        let name = format!("host_{ordinal}");
        write_fasta_record(&mut writer, &name, &record.sequence)?;
        contigs.push(ReferenceContig {
            name,
            role: ReferenceRole::Host,
            target_group_ordinal: None,
        });
    }
    for group in &groups {
        let record = by_id
            .get(group.representative_id.as_str())
            .ok_or_else(|| "ReferenceGroup representative is absent".to_string())?;
        write_fasta_record(&mut writer, &group.contig_name, &record.sequence)?;
        contigs.push(ReferenceContig {
            name: group.contig_name.clone(),
            role: ReferenceRole::Target,
            target_group_ordinal: Some(group.ordinal),
        });
    }
    writer
        .flush()
        .map_err(|error| format!("Failed to write {}: {error}", composite.display()))?;

    let mmi = directory.join(MMI);
    let mmi_text = mmi
        .to_str()
        .ok_or_else(|| "Index path is not valid UTF-8".to_string())?;
    let mut builder = Aligner::builder().sr().with_index_threads(options.threads);
    builder.idxopt.batch_size = SINGLE_PART_BATCH_SIZE;
    let built = builder
        .with_index(&composite, Some(mmi_text))
        .map_err(|error| format!("Failed to build minimap2 index: {error}"))?;
    if built.idx_parts.len() != 1 {
        return Err("minimap2 produced a multi-part index".into());
    }

    let bloom = TargetKmerBloom::build(&target_records, profile.kmer_length);
    bloom.write(&directory.join(BLOOM))?;
    let target_digest = file_digest(&options.target_fa)?;
    let host_digest = file_digest(&options.host_fa)?;
    let ledger_bytes = ledger_bytes(&groups, &target_digest, &profile.digest());
    write_new(&directory.join(LEDGER), &ledger_bytes)?;
    let ledger_digest = hex_sha256(&ledger_bytes);
    let reference_set_digest = reference_set_digest(&host_digest, &target_digest, &ledger_digest)?;
    let manifest = IndexManifest {
        contract_id: INDEX_CONTRACT.into(),
        profile_digest: profile.digest(),
        reference_set_digest,
        ledger_digest,
        mmi_digest: file_digest(&mmi)?,
        bloom_digest: file_digest(&directory.join(BLOOM))?,
        kmer_length: profile.kmer_length,
        host_fasta_sha256: host_digest,
        target_fasta_sha256: target_digest,
        bloom: bloom.summary(),
        contigs,
    };
    let bytes = serde_json::to_vec_pretty(&manifest)
        .map_err(|error| format!("Cannot serialize index manifest: {error}"))?;
    write_new(
        &directory.join(MANIFEST),
        &[bytes.as_slice(), b"\n"].concat(),
    )
}

pub fn load_index(directory: &Path) -> Result<ReferenceIndex, String> {
    let manifest_path = directory.join(MANIFEST);
    let bytes = std::fs::read(&manifest_path).map_err(|error| {
        format!(
            "Cannot read current index manifest {}: {error}",
            manifest_path.display()
        )
    })?;
    let manifest: IndexManifest = serde_json::from_slice(&bytes).map_err(|_| {
        "Index manifest is not the current profile-bound format; rebuild it".to_string()
    })?;
    let profile = AnalysisProfile::FROZEN;
    if manifest.contract_id != INDEX_CONTRACT
        || manifest.profile_digest != profile.digest()
        || manifest.kmer_length != profile.kmer_length
    {
        return Err("Index was not built for the frozen AnalysisProfile; rebuild it".into());
    }
    verify_digest(&directory.join(MMI), &manifest.mmi_digest)?;
    verify_digest(&directory.join(BLOOM), &manifest.bloom_digest)?;
    verify_digest(&directory.join(LEDGER), &manifest.ledger_digest)?;
    let bloom = TargetKmerBloom::read(&directory.join(BLOOM), profile.kmer_length)?;
    if bloom.summary() != manifest.bloom {
        return Err("Bloom summary does not match bloom.bin".into());
    }
    let ledger_bytes = std::fs::read(directory.join(LEDGER))
        .map_err(|error| format!("Cannot read ReferenceGroup ledger: {error}"))?;
    let target_groups = parse_ledger(
        &ledger_bytes,
        &manifest.target_fasta_sha256,
        &manifest.profile_digest,
    )?;
    let contigs = validate_contigs(manifest.contigs, target_groups.len())?;
    let expected_reference_set_digest = reference_set_digest(
        &manifest.host_fasta_sha256,
        &manifest.target_fasta_sha256,
        &manifest.ledger_digest,
    )?;
    if manifest.reference_set_digest != expected_reference_set_digest {
        return Err("Index reference-set digest is inconsistent".into());
    }
    Ok(ReferenceIndex {
        mmi_path: directory.join(MMI),
        bloom,
        bloom_summary: manifest.bloom,
        contigs,
        target_groups,
        profile_digest: manifest.profile_digest,
        index_digest: hex_sha256(&bytes),
    })
}

fn parse_ledger(
    bytes: &[u8],
    target_digest: &str,
    profile_digest: &str,
) -> Result<Vec<ReferenceGroup>, String> {
    const HEADER: &str = "group_ordinal\ttarget_group_id\trepresentative_id\tmember_ordinal\tmember_id\trepresentative_length\ttarget_fasta_sha256\tprofile_digest\n";
    let text = std::str::from_utf8(bytes)
        .map_err(|error| format!("ReferenceGroup ledger is not UTF-8: {error}"))?;
    let rows = text.strip_prefix(HEADER).ok_or_else(|| {
        "ReferenceGroup ledger header does not match the frozen contract".to_string()
    })?;
    if rows.is_empty() || !rows.ends_with('\n') {
        return Err("ReferenceGroup ledger has no canonical data rows".into());
    }
    let mut groups: Vec<ReferenceGroup> = Vec::new();
    let mut members = HashSet::new();
    for row in rows.lines() {
        let fields = row.split('\t').collect::<Vec<_>>();
        if fields.len() != 8 {
            return Err("ReferenceGroup ledger row does not have eight fields".into());
        }
        let group_ordinal = parse_ledger_number(fields[0], "group_ordinal")?;
        let member_ordinal = parse_ledger_number(fields[3], "member_ordinal")?;
        let representative_length = fields[5]
            .parse::<u64>()
            .map_err(|_| "ReferenceGroup ledger representative_length is invalid".to_string())?;
        if representative_length == 0
            || fields[1]
                .strip_prefix("sha256:")
                .is_none_or(|digest| !is_sha256(digest))
            || fields[2].is_empty()
            || fields[4].is_empty()
            || fields[6] != target_digest
            || fields[7] != profile_digest
        {
            return Err("ReferenceGroup ledger row violates the frozen contract".into());
        }
        if !members.insert(fields[4].to_string()) {
            return Err(format!(
                "ReferenceGroup member appears more than once: {}",
                fields[4]
            ));
        }
        if group_ordinal == groups.len() {
            if member_ordinal != 0
                || groups
                    .last()
                    .is_some_and(|previous| previous.target_group_id.as_str() >= fields[1])
            {
                return Err("ReferenceGroup ledger group order is not dense and canonical".into());
            }
            groups.push(ReferenceGroup {
                ordinal: group_ordinal,
                target_group_id: fields[1].to_string(),
                representative_id: fields[2].to_string(),
                member_ids: vec![fields[4].to_string()],
                representative_length,
                contig_name: format!("target_{group_ordinal}"),
            });
        } else if group_ordinal + 1 == groups.len() {
            let group = groups.last_mut().expect("current group exists");
            if group.target_group_id != fields[1]
                || group.representative_id != fields[2]
                || group.representative_length != representative_length
                || member_ordinal != group.member_ids.len()
                || group
                    .member_ids
                    .last()
                    .is_some_and(|previous| previous.as_str() >= fields[4])
            {
                return Err("ReferenceGroup ledger member rows are not canonical".into());
            }
            group.member_ids.push(fields[4].to_string());
        } else {
            return Err("ReferenceGroup ledger group ordinals are not dense".into());
        }
    }
    if groups
        .iter()
        .any(|group| group.representative_id != group.member_ids[0])
    {
        return Err("ReferenceGroup representative is not the first canonical member".into());
    }
    Ok(groups)
}

fn parse_ledger_number(value: &str, field: &str) -> Result<usize, String> {
    value
        .parse()
        .map_err(|_| format!("ReferenceGroup ledger {field} is invalid"))
}

fn is_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
}

fn validate_contigs(
    contigs: Vec<ReferenceContig>,
    group_count: usize,
) -> Result<HashMap<String, ReferenceContig>, String> {
    let mut by_name = HashMap::new();
    let mut target_ordinals = vec![false; group_count];
    let mut host_ordinal = 0;
    for contig in contigs {
        match contig.role {
            ReferenceRole::Host => {
                if contig.target_group_ordinal.is_some()
                    || contig.name != format!("host_{host_ordinal}")
                {
                    return Err("Index host contig mapping is not dense".into());
                }
                host_ordinal += 1;
            }
            ReferenceRole::Target => {
                let ordinal = contig
                    .target_group_ordinal
                    .ok_or_else(|| "Index target contig has no group ordinal".to_string())?;
                if ordinal >= group_count
                    || target_ordinals[ordinal]
                    || contig.name != format!("target_{ordinal}")
                {
                    return Err(
                        "Index target contig/group mapping is not dense and one-to-one".into(),
                    );
                }
                target_ordinals[ordinal] = true;
            }
        }
        if by_name.insert(contig.name.clone(), contig).is_some() {
            return Err("Index contains duplicate contig names".into());
        }
    }
    if host_ordinal == 0 || target_ordinals.iter().any(|mapped| !mapped) {
        return Err("Index does not contain complete HOST and TARGET mappings".into());
    }
    Ok(by_name)
}

fn write_fasta_record(writer: &mut impl Write, name: &str, sequence: &[u8]) -> Result<(), String> {
    writeln!(writer, ">{name}").map_err(|error| error.to_string())?;
    for chunk in sequence.chunks(60) {
        writer
            .write_all(chunk)
            .and_then(|_| writer.write_all(b"\n"))
            .map_err(|error| error.to_string())?;
    }
    Ok(())
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|error| format!("Cannot create {}: {error}", path.display()))?;
    file.write_all(bytes)
        .and_then(|_| file.flush())
        .map_err(|error| format!("Cannot write {}: {error}", path.display()))
}

fn file_digest(path: &Path) -> Result<String, String> {
    let mut file =
        File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?;
    let mut hasher = Sha256::new();
    let mut buffer = [0; 1 << 20];
    loop {
        let count = file
            .read(&mut buffer)
            .map_err(|error| format!("Cannot read {}: {error}", path.display()))?;
        if count == 0 {
            break;
        }
        hasher.update(&buffer[..count]);
    }
    Ok(format!("{:x}", hasher.finalize()))
}

fn verify_digest(path: &Path, expected: &str) -> Result<(), String> {
    if file_digest(path)? == expected {
        Ok(())
    } else {
        Err(format!(
            "Index artifact digest mismatch: {}",
            path.display()
        ))
    }
}

fn ledger_bytes(groups: &[ReferenceGroup], target_digest: &str, profile_digest: &str) -> Vec<u8> {
    let mut ledger = String::from("group_ordinal\ttarget_group_id\trepresentative_id\tmember_ordinal\tmember_id\trepresentative_length\ttarget_fasta_sha256\tprofile_digest\n");
    for group in groups {
        for (member_ordinal, member) in group.member_ids.iter().enumerate() {
            ledger.push_str(&format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
                group.ordinal,
                group.target_group_id,
                group.representative_id,
                member_ordinal,
                member,
                group.representative_length,
                target_digest,
                profile_digest
            ));
        }
    }
    ledger.into_bytes()
}

fn reference_set_digest(host: &str, target: &str, ledger: &str) -> Result<String, String> {
    let mut framed = b"viroflash-reference-set-v1\0".to_vec();
    for (tag, digest) in [
        ("host_fasta", host),
        ("target_fasta", target),
        ("reference_group_ledger", ledger),
    ] {
        framed.extend_from_slice(&(tag.len() as u16).to_be_bytes());
        framed.extend_from_slice(tag.as_bytes());
        let raw = decode_hex(digest)?;
        framed.extend_from_slice(&raw);
    }
    Ok(hex_sha256(&framed))
}

fn decode_hex(value: &str) -> Result<Vec<u8>, String> {
    if value.len() != 64 {
        return Err("Digest is not SHA-256".into());
    }
    (0..value.len())
        .step_by(2)
        .map(|index| {
            u8::from_str_radix(&value[index..index + 2], 16)
                .map_err(|_| "Digest is not hexadecimal".to_string())
        })
        .collect()
}
