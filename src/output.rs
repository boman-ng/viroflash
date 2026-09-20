use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;

pub(crate) fn atomic_write(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let temporary = path.with_extension(format!(
        "{}.part.{}",
        path.extension()
            .and_then(|value| value.to_str())
            .unwrap_or("tmp"),
        std::process::id()
    ));
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| format!("Cannot create {}: {error}", temporary.display()))?;
    let result = file
        .write_all(bytes)
        .and_then(|_| file.flush())
        .map_err(|error| format!("Cannot write {}: {error}", temporary.display()))
        .and_then(|_| {
            std::fs::rename(&temporary, path)
                .map_err(|error| format!("Cannot finalize {}: {error}", path.display()))
        });
    if result.is_err() {
        let _ = std::fs::remove_file(temporary);
    }
    result
}
