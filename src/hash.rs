//! BLAKE3 文件校验和（薄封装）：索引 manifest 中记录来源 FASTA 的校验和，
//! 保证索引可审计。不手写哈希原语，直接使用 vetted `blake3` crate
//! （cargo 自身的校验和也基于 blake3）；GB 级 FASTA 上比标量 SHA-256
//! 快约一个量级，索引构建的校验开销可忽略。

use std::io::Read;
use std::path::Path;

/// 数据块 BLAKE3（hex）。
pub fn blake3_hex(data: &[u8]) -> String {
    blake3::hash(data).to_hex().to_string()
}

/// 文件流式 BLAKE3（hex），1 MiB 缓冲，适配 GB 级 FASTA。
pub fn blake3_file_hex(path: &Path) -> Result<String, String> {
    let file = File::open(path).map_err(|e| format!("无法打开 {}: {e}", path.display()))?;
    let mut reader = std::io::BufReader::with_capacity(1 << 20, file);
    let mut hasher = blake3::Hasher::new();
    let mut buf = vec![0u8; 1 << 20];
    loop {
        let n = reader
            .read(&mut buf)
            .map_err(|e| format!("读取 {} 失败: {e}", path.display()))?;
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

    /// BLAKE3 官方公开测试向量（https://github.com/BLAKE3-team/BLAKE3/blob/master/test_vectors/test_vectors.json）。
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
        // 各 offset 切分的一致性：分块 update ≡ 一次哈希。
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
