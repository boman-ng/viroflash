use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;
use std::time::Instant;

use serde::Serialize;
use sysinfo::{get_current_pid, ProcessRefreshKind, ProcessesToUpdate, System};

#[derive(Debug, Clone, Default, Serialize)]
pub struct StageTimes {
    pub pass1_count: u64,
    pub pass2_sample_prescreen_align: u64,
    pub report_write: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct PerformanceReport {
    pub schema_id: &'static str,
    pub status: &'static str,
    pub sample_id: String,
    pub wall_time_ms: u64,
    pub configured_threads: usize,
    pub process_cpu_time_ms: u64,
    pub peak_rss_bytes: u64,
    pub read_bytes: u64,
    pub written_bytes: u64,
    pub pass1_count: u64,
    pub pass2_sample_prescreen_align: u64,
    pub report_write: u64,
    pub input_fragments: u64,
    pub selected_fragments: u64,
    pub prescreen_passed_fragments: u64,
    pub aligned_fragments: u64,
    pub telemetry_status: &'static str,
    pub sample_errors: Vec<String>,
}

pub struct PerformanceMonitor {
    started: Instant,
    cpu_started: u64,
    read_started: u64,
    write_started: u64,
    peak_rss: u64,
    telemetry_available: bool,
}

impl PerformanceMonitor {
    pub fn start() -> Self {
        let metrics = process_metrics();
        let (cpu, rss, read, write) = metrics.unwrap_or_default();
        Self {
            started: Instant::now(),
            cpu_started: cpu,
            read_started: read,
            write_started: write,
            peak_rss: rss,
            telemetry_available: metrics.is_some(),
        }
    }
    fn observe(&mut self) {
        if let Some((_, rss, _, _)) = process_metrics() {
            self.peak_rss = self.peak_rss.max(rss);
        } else {
            self.telemetry_available = false;
        }
    }
    pub fn finish(
        mut self,
        status: &'static str,
        sample_id: String,
        threads: usize,
        stages: StageTimes,
        counts: (u64, u64, u64, u64),
        errors: Vec<String>,
    ) -> PerformanceReport {
        self.observe();
        let metrics = process_metrics();
        let (cpu, _, read, write) = metrics.unwrap_or_default();
        self.telemetry_available &= metrics.is_some();
        let high_water_rss = high_water_rss_bytes();
        self.telemetry_available &= high_water_rss.is_some();
        self.peak_rss = self.peak_rss.max(high_water_rss.unwrap_or_default());
        PerformanceReport {
            schema_id: "viroflash.perf.v1",
            status,
            sample_id,
            wall_time_ms: self.started.elapsed().as_millis() as u64,
            configured_threads: threads,
            process_cpu_time_ms: cpu.saturating_sub(self.cpu_started),
            peak_rss_bytes: self.peak_rss,
            read_bytes: read.saturating_sub(self.read_started),
            written_bytes: write.saturating_sub(self.write_started),
            pass1_count: stages.pass1_count,
            pass2_sample_prescreen_align: stages.pass2_sample_prescreen_align,
            report_write: stages.report_write,
            input_fragments: counts.0,
            selected_fragments: counts.1,
            prescreen_passed_fragments: counts.2,
            aligned_fragments: counts.3,
            telemetry_status: if self.telemetry_available {
                "COMPLETE"
            } else {
                "UNAVAILABLE"
            },
            sample_errors: errors,
        }
    }
}

#[cfg(target_os = "linux")]
fn high_water_rss_bytes() -> Option<u64> {
    parse_linux_high_water_rss(&std::fs::read_to_string("/proc/self/status").ok()?)
}

#[cfg(not(target_os = "linux"))]
fn high_water_rss_bytes() -> Option<u64> {
    None
}

#[cfg(target_os = "linux")]
fn parse_linux_high_water_rss(status: &str) -> Option<u64> {
    let value = status
        .lines()
        .find_map(|line| line.strip_prefix("VmHWM:"))?
        .split_ascii_whitespace()
        .next()?
        .parse::<u64>()
        .ok()?;
    value.checked_mul(1024)
}

fn process_metrics() -> Option<(u64, u64, u64, u64)> {
    let pid = get_current_pid().ok()?;
    let mut system = System::new();
    system.refresh_processes_specifics(
        ProcessesToUpdate::Some(&[pid]),
        true,
        ProcessRefreshKind::nothing()
            .with_cpu()
            .with_memory()
            .with_disk_usage(),
    );
    system.process(pid).map(|process| {
        let disk = process.disk_usage();
        (
            process.accumulated_cpu_time(),
            process.memory(),
            disk.total_read_bytes,
            disk.total_written_bytes,
        )
    })
}

pub fn write_perf_json(path: &Path, report: &PerformanceReport) -> Result<(), String> {
    let mut bytes = serde_json::to_vec_pretty(report)
        .map_err(|error| format!("Cannot serialize perf.json: {error}"))?;
    bytes.push(b'\n');
    let temporary = path.with_extension(format!("json.part.{}", std::process::id()));
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temporary)
        .map_err(|error| format!("Cannot create {}: {error}", temporary.display()))?;
    let result = file
        .write_all(&bytes)
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

pub fn stage_start() -> Instant {
    Instant::now()
}

#[cfg(all(test, target_os = "linux"))]
mod tests {
    use super::*;

    #[test]
    fn parses_linux_resident_set_high_water_mark() {
        assert_eq!(
            parse_linux_high_water_rss("Name:\tviroflash\nVmHWM:\t1234 kB\nVmRSS:\t1000 kB\n"),
            Some(1_263_616)
        );
    }
}
