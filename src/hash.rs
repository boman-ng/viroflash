//! Thin BLAKE3 checksum wrapper. Index manifests record source FASTA checksums for auditability.
//! The vetted `blake3` crate avoids a bespoke hash implementation and makes checksum overhead
//! negligible even for gigabyte-scale FASTA files.

use std::io::Read;
use std::path::Path;

/// Return the BLAKE3 digest of a byte slice as hexadecimal text.
pub fn blake3_hex(data: &[u8]) -> String {
    blake3::hash(data).to_hex().to_string()
}

/// Stream a file through BLAKE3 with a 1 MiB buffer and return a hexadecimal digest.
pub fn blake3_file_hex(path: &Path) -> Result<String, String> {
    let file = File::open(path).map_err(|e| format!("Cannot open {}: {e}", path.display()))?;
    let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut hasher = blake3::Hasher::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| format!("Failed to read {}: {e}", path.display()))?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
    }
    Ok(hasher.finalize().to_hex().to_string())
}

use std::fs::File;

#[cfg(test)]
mod tests {
    use super::*;

    /// Official BLAKE3 test vector.
    #[test]
    fn official_vectors() {
        assert_eq!(
            blake3_hex(b""),
            "af1349b9f5f9a1a6a0404dea36dcc9499bcb25c9adc112b7cc9a93cae41f3262"
        );
        assert_eq!(
            blake3_hex(b"abc"),
            "6437b3ac38465133ffb63b75273a8db548c558465d79db03fd359c6cd5bd9d85"
        );
    }

    #[test]
    fn incremental_matches_oneshot() {
        // Updating at different offsets must match one-shot hashing.
        let data: Vec<u8> = (0..1000u32).map(|i| (i % 251) as u8).collect();
        for split in [0usize, 1, 63, 64, 65, 127, 128, 500, 999, 1000] {
            let mut h = blake3::Hasher::new();
            h.update(&data[..split]);
            h.update(&data[split..]);
            assert_eq!(
                h.finalize().to_hex().to_string(),
                blake3_hex(&data),
                "split={split}"
            );
        }
    }
}
