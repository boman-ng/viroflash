use std::fs;
use std::io;
use std::path::PathBuf;

fn sequence(seed: u128, length: usize) -> String {
    let mut state = seed;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            state &= u64::MAX as u128;
            b"ACGT"[(state & 3) as usize] as char
        })
        .collect()
}

fn main() -> io::Result<()> {
    let root = PathBuf::from(std::env::args_os().nth(1).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "Expected an output directory")
    })?);
    fs::create_dir_all(&root)?;
    let host = sequence(17, 600);
    let target = sequence(91, 600);
    fs::write(root.join("host.fa"), format!(">host\n{host}\n"))?;
    fs::write(root.join("target.fa"), format!(">target\n{target}\n"))?;
    let mut fastq = String::new();
    for ordinal in 0..20 {
        let start = ordinal % 30;
        let read = &target[start..start + 120];
        fastq.push_str(&format!(
            "@fragment-{ordinal}\n{read}\n+\n{}\n",
            "I".repeat(120)
        ));
    }
    fs::write(root.join("sample.fastq"), fastq)
}
