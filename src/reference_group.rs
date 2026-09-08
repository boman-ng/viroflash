use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use serde::{Deserialize, Serialize};

use crate::analysis_profile::hex_sha256;

const IUPAC_DNA: &[u8] = b"ACGTMRWSYKVHDBN";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FastaRecord {
    pub id: String,
    pub sequence: Vec<u8>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ReferenceGroup {
    pub ordinal: usize,
    pub target_group_id: String,
    pub representative_id: String,
    pub member_ids: Vec<String>,
    pub representative_length: u64,
    pub contig_name: String,
}

pub fn parse_fasta(path: &Path) -> Result<Vec<FastaRecord>, String> {
    let file =
        File::open(path).map_err(|error| format!("Cannot open {}: {error}", path.display()))?;
    let mut records = Vec::new();
    let mut ids = HashSet::new();
    let mut current_id: Option<String> = None;
    let mut sequence = Vec::new();

    for (line_number, line) in BufReader::new(file).lines().enumerate() {
        let line = line.map_err(|error| format!("Failed to read {}: {error}", path.display()))?;
        if let Some(header) = line.strip_prefix('>') {
            if let Some(id) = current_id.take() {
                push_record(path, id, &mut sequence, &mut records)?;
            }
            let id = header
                .split_ascii_whitespace()
                .next()
                .filter(|value| !value.is_empty())
                .ok_or_else(|| {
                    format!(
                        "{}:{} has an empty FASTA record ID",
                        path.display(),
                        line_number + 1
                    )
                })?
                .to_string();
            if !ids.insert(id.clone()) {
                return Err(format!(
                    "{} contains duplicate FASTA record ID: {id}",
                    path.display()
                ));
            }
            current_id = Some(id);
        } else {
            if current_id.is_none() && !line.trim().is_empty() {
                return Err(format!(
                    "{}:{} has sequence before the first FASTA header",
                    path.display(),
                    line_number + 1
                ));
            }
            for byte in line.bytes().filter(|byte| !byte.is_ascii_whitespace()) {
                let upper = byte.to_ascii_uppercase();
                if !IUPAC_DNA.contains(&upper) {
                    return Err(format!(
                        "{}:{} contains invalid IUPAC DNA symbol {:?}",
                        path.display(),
                        line_number + 1,
                        char::from(byte)
                    ));
                }
                sequence.push(upper);
            }
        }
    }
    if let Some(id) = current_id {
        push_record(path, id, &mut sequence, &mut records)?;
    }
    if records.is_empty() {
        return Err(format!("{} contains no FASTA records", path.display()));
    }
    Ok(records)
}

fn push_record(
    path: &Path,
    id: String,
    sequence: &mut Vec<u8>,
    records: &mut Vec<FastaRecord>,
) -> Result<(), String> {
    if sequence.is_empty() {
        return Err(format!(
            "{} contains an empty FASTA sequence: {id}",
            path.display()
        ));
    }
    records.push(FastaRecord {
        id,
        sequence: std::mem::take(sequence),
    });
    Ok(())
}

pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'M' => b'K',
            b'K' => b'M',
            b'R' => b'Y',
            b'Y' => b'R',
            b'W' => b'W',
            b'S' => b'S',
            b'V' => b'B',
            b'B' => b'V',
            b'H' => b'D',
            b'D' => b'H',
            b'N' => b'N',
            _ => unreachable!("validated IUPAC sequence"),
        })
        .collect()
}

pub fn build_reference_groups(records: &[FastaRecord]) -> Vec<ReferenceGroup> {
    let mut classes: BTreeMap<Vec<u8>, Vec<&FastaRecord>> = BTreeMap::new();
    for record in records {
        let reverse = reverse_complement(&record.sequence);
        let canonical = if record.sequence <= reverse {
            record.sequence.clone()
        } else {
            reverse
        };
        classes.entry(canonical).or_default().push(record);
    }
    let mut grouped = classes
        .into_iter()
        .map(|(canonical, mut members)| {
            members.sort_by(|left, right| left.id.as_bytes().cmp(right.id.as_bytes()));
            ReferenceGroup {
                ordinal: 0,
                target_group_id: format!("sha256:{}", hex_sha256(&canonical)),
                representative_id: members[0].id.clone(),
                member_ids: members.iter().map(|member| member.id.clone()).collect(),
                representative_length: canonical.len() as u64,
                contig_name: String::new(),
            }
        })
        .collect::<Vec<_>>();
    grouped.sort_by(|left, right| left.target_group_id.cmp(&right.target_group_id));
    for (ordinal, group) in grouped.iter_mut().enumerate() {
        group.ordinal = ordinal;
        group.contig_name = format!("target_{ordinal}");
    }
    grouped
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn groups_only_exact_sequences_and_reverse_complements() {
        let records = vec![
            FastaRecord {
                id: "b".into(),
                sequence: b"ACGMRN".to_vec(),
            },
            FastaRecord {
                id: "a".into(),
                sequence: reverse_complement(b"ACGMRN"),
            },
            FastaRecord {
                id: "c".into(),
                sequence: b"ACGMRR".to_vec(),
            },
        ];
        let groups = build_reference_groups(&records);
        assert_eq!(groups.len(), 2);
        assert!(groups.iter().any(|group| group.member_ids == ["a", "b"]));
        assert!(groups.iter().any(|group| group.member_ids == ["c"]));
    }

    #[test]
    fn parser_rejects_invalid_iupac_instead_of_dropping_it() {
        let path =
            std::env::temp_dir().join(format!("viroflash-invalid-{}.fa", std::process::id()));
        let mut file = File::create(&path).unwrap();
        writeln!(file, ">x\nACGTZ").unwrap();
        drop(file);
        let error = parse_fasta(&path).unwrap_err();
        assert!(error.contains("invalid IUPAC"));
        let _ = std::fs::remove_file(path);
    }
}
