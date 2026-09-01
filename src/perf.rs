//! In-process performance telemetry with low-frequency targeted refreshes of the viroflash process.
//!
//! `sysinfo` provides portable process CPU, RSS, virtual-memory, I/O, and system-level metrics.
//! PSS, NUMA, cache misses, and memory bandwidth require dedicated external instrumentation.

use std::fs::OpenOptions;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::{mpsc, Mutex};
use std::thread::JoinHandle;
use std::time::{Duration, Instant};

use sysinfo::{
    get_current_pid, Pid, ProcessRefreshKind, ProcessesToUpdate, System,
    MINIMUM_CPU_UPDATE_INTERVAL,
};

use crate::report::json_escape;

const SCHEMA: &str = "viroflash.perf.v1";
const SYSINFO_VERSION: &str = "0.38.4";
const SAMPLE_INTERVAL: Duration = Duration::from_secs(1);

pub fn report_path(out_prefix: &Path) -> PathBuf {
    PathBuf::from(format!("{}.perf.json", out_prefix.display()))
}

#[derive(Debug, Clone, Default)]
struct CpuAggregate {
    sum: f64,
    count: u64,
    peak: f64,
}

impl CpuAggregate {
    fn observe(&mut self, value: f32) {
        let value = f64::from(value);
        self.sum += value;
        self.count += 1;
        self.peak = self.peak.max(value);
    }

    fn mean(&self) -> Option<f64> {
        (self.count > 0).then(|| self.sum / self.count as f64)
    }

    fn peak(&self) -> Option<f64> {
        (self.count > 0).then_some(self.peak)
    }
}

#[derive(Debug, Clone)]
struct RawMetrics {
    at: Instant,
    process_cpu_pct: f32,
    process_cpu_time_ms: u64,
    rss_bytes: u64,
    virtual_bytes: u64,
    read_bytes_total: u64,
    written_bytes_total: u64,
    system_cpu_pct: f32,
    system_memory_total_bytes: u64,
    system_memory_available_bytes: u64,
}

#[derive(Debug, Clone)]
struct StageAggregate {
    name: &'static str,
    started: Instant,
    ended: Option<Instant>,
    process_cpu: CpuAggregate,
    rss_bytes_peak: u64,
}

impl StageAggregate {
    fn wall_time_ms(&self, fallback_end: Instant) -> u64 {
        self.ended
            .unwrap_or(fallback_end)
            .saturating_duration_since(self.started)
            .as_millis() as u64
    }
}

#[derive(Debug)]
struct Collector {
    started: Instant,
    last_cpu_refresh_at: Instant,
    initial: RawMetrics,
    latest: RawMetrics,
    process_cpu: CpuAggregate,
    system_cpu: CpuAggregate,
    rss_bytes_peak: u64,
    virtual_bytes_peak: u64,
    system_memory_available_bytes_min: u64,
    sample_count: u64,
    sample_errors: u64,
    stages: Vec<StageAggregate>,
}

impl Collector {
    fn new(initial: RawMetrics) -> Self {
        let started = initial.at;
        let initial_rss_bytes = initial.rss_bytes;
        Self {
            started,
            last_cpu_refresh_at: started,
            rss_bytes_peak: initial.rss_bytes,
            virtual_bytes_peak: initial.virtual_bytes,
            system_memory_available_bytes_min: initial.system_memory_available_bytes,
            sample_count: 1,
            latest: initial.clone(),
            initial,
            process_cpu: CpuAggregate::default(),
            system_cpu: CpuAggregate::default(),
            sample_errors: 0,
            stages: vec![StageAggregate {
                name: "startup",
                started,
                ended: None,
                process_cpu: CpuAggregate::default(),
                rss_bytes_peak: initial_rss_bytes,
            }],
        }
    }

    fn observe(&mut self, raw: RawMetrics) {
        let cpu_valid = raw.at.saturating_duration_since(self.last_cpu_refresh_at)
            >= MINIMUM_CPU_UPDATE_INTERVAL;
        self.last_cpu_refresh_at = raw.at;
        if cpu_valid {
            self.process_cpu.observe(raw.process_cpu_pct);
            self.system_cpu.observe(raw.system_cpu_pct);
        }
        self.rss_bytes_peak = self.rss_bytes_peak.max(raw.rss_bytes);
        self.virtual_bytes_peak = self.virtual_bytes_peak.max(raw.virtual_bytes);
        self.system_memory_available_bytes_min = self
            .system_memory_available_bytes_min
            .min(raw.system_memory_available_bytes);
        self.sample_count += 1;
        if let Some(stage) = self.stages.last_mut() {
            if cpu_valid {
                stage.process_cpu.observe(raw.process_cpu_pct);
            }
            stage.rss_bytes_peak = stage.rss_bytes_peak.max(raw.rss_bytes);
        }
        self.latest = raw;
    }

    fn sample_error(&mut self) {
        self.sample_errors += 1;
    }

    fn stage(&mut self, name: &'static str, at: Instant) {
        if self.stages.last().is_some_and(|s| s.name == name) {
            return;
        }
        if let Some(stage) = self.stages.last_mut() {
            stage.ended = Some(at.max(stage.started));
        }
        self.stages.push(StageAggregate {
            name,
            started: at,
            ended: None,
            process_cpu: CpuAggregate::default(),
            rss_bytes_peak: self.latest.rss_bytes,
        });
    }

    fn finish(mut self, at: Instant) -> Collected {
        if let Some(stage) = self.stages.last_mut() {
            stage.ended = Some(at.max(stage.started));
        }
        let wall_time_ms = at.saturating_duration_since(self.started).as_millis() as u64;
        let cpu_time_ms = self
            .latest
            .process_cpu_time_ms
            .saturating_sub(self.initial.process_cpu_time_ms);
        let process_cpu_pct_mean =
            (wall_time_ms > 0).then(|| cpu_time_ms as f64 / wall_time_ms as f64 * 100.0);
        Collected {
            wall_time_ms,
            sample_count: self.sample_count,
            valid_cpu_samples: self.process_cpu.count,
            sample_errors: self.sample_errors,
            process_cpu_pct_mean,
            process_cpu_pct_peak: self.process_cpu.peak(),
            process_cpu_time_ms: cpu_time_ms,
            rss_bytes_peak: self.rss_bytes_peak,
            virtual_bytes_peak: self.virtual_bytes_peak,
            read_bytes: self
                .latest
                .read_bytes_total
                .saturating_sub(self.initial.read_bytes_total),
            written_bytes: self
                .latest
                .written_bytes_total
                .saturating_sub(self.initial.written_bytes_total),
            system_cpu_pct_mean: self.system_cpu.mean(),
            system_cpu_pct_peak: self.system_cpu.peak(),
            system_memory_total_bytes: self.latest.system_memory_total_bytes,
            system_memory_available_bytes_min: self.system_memory_available_bytes_min,
            stages: self.stages,
            finished: at,
        }
    }
}

#[derive(Debug)]
struct Collected {
    wall_time_ms: u64,
    sample_count: u64,
    valid_cpu_samples: u64,
    sample_errors: u64,
    process_cpu_pct_mean: Option<f64>,
    process_cpu_pct_peak: Option<f64>,
    process_cpu_time_ms: u64,
    rss_bytes_peak: u64,
    virtual_bytes_peak: u64,
    read_bytes: u64,
    written_bytes: u64,
    system_cpu_pct_mean: Option<f64>,
    system_cpu_pct_peak: Option<f64>,
    system_memory_total_bytes: u64,
    system_memory_available_bytes_min: u64,
    stages: Vec<StageAggregate>,
    finished: Instant,
}

enum Control {
    Stage(&'static str, Instant),
    Stop,
}

enum Backend {
    Periodic {
        tx: mpsc::Sender<Control>,
        handle: JoinHandle<Result<Collected, String>>,
    },
    Boundary(Box<BoundaryBackend>),
}

struct BoundaryBackend {
    state: Mutex<BoundaryState>,
}

struct BoundaryState {
    system: System,
    pid: Pid,
    collector: Collector,
}

pub(crate) struct PerfMonitor {
    command: &'static str,
    label: String,
    configured_threads: usize,
    logical_cpu_count: usize,
    physical_core_count: Option<usize>,
    interval: Duration,
    path: PathBuf,
    backend: Backend,
}

impl PerfMonitor {
    pub(crate) fn start(
        command: &'static str,
        label: String,
        out_prefix: &Path,
        configured_threads: usize,
    ) -> Result<Self, String> {
        let pid = get_current_pid()
            .map_err(|e| format!("Performance monitor cannot obtain the current PID: {e}"))?;
        let path = report_path(out_prefix);
        let interval = SAMPLE_INTERVAL.max(MINIMUM_CPU_UPDATE_INTERVAL);
        let physical_core_count = System::physical_core_count();

        if configured_threads <= 1 {
            let mut system = System::new();
            let initial = refresh(&mut system, pid)?;
            let logical_cpu_count = system.cpus().len();
            return Ok(Self {
                command,
                label,
                configured_threads,
                logical_cpu_count,
                physical_core_count,
                interval,
                path,
                backend: Backend::Boundary(Box::new(BoundaryBackend {
                    state: Mutex::new(BoundaryState {
                        system,
                        pid,
                        collector: Collector::new(initial),
                    }),
                })),
            });
        }

        let (tx, rx) = mpsc::channel();
        let (ready_tx, ready_rx) = mpsc::sync_channel(1);
        let handle = std::thread::Builder::new()
            .name("viroflash-perf".to_string())
            .spawn(move || sampler_loop(pid, interval, rx, ready_tx))
            .map_err(|e| format!("Cannot start performance sampling thread: {e}"))?;
        let logical_cpu_count = match ready_rx.recv() {
            Ok(Ok(n)) => n,
            Ok(Err(e)) => {
                let _ = handle.join();
                return Err(e);
            }
            Err(e) => {
                let _ = handle.join();
                return Err(format!(
                    "Performance sampling thread failed to initialize: {e}"
                ));
            }
        };
        Ok(Self {
            command,
            label,
            configured_threads,
            logical_cpu_count,
            physical_core_count,
            interval,
            path,
            backend: Backend::Periodic { tx, handle },
        })
    }

    /// `--threads` remains the data-pipeline worker budget. The low-frequency telemetry sampler
    /// consumes no mapper or decompressor slot, preserving throughput semantics when enabled.
    pub(crate) fn work_thread_budget(&self) -> usize {
        self.configured_threads
    }

    pub(crate) fn stage(&self, name: &'static str) {
        let at = Instant::now();
        match &self.backend {
            Backend::Periodic { tx, .. } => {
                let _ = tx.send(Control::Stage(name, at));
            }
            Backend::Boundary(boundary) => {
                if let Ok(mut state) = boundary.state.lock() {
                    let pid = state.pid;
                    match refresh(&mut state.system, pid) {
                        Ok(raw) => state.collector.observe(raw),
                        Err(_) => state.collector.sample_error(),
                    }
                    state.collector.stage(name, at);
                }
            }
        }
    }

    pub(crate) fn complete<T>(self, result: Result<T, String>) -> Result<T, String> {
        let success = result.is_ok();
        let perf_result = self.finish(success);
        match (result, perf_result) {
            (Ok(value), Ok(())) => Ok(value),
            (Err(primary), Ok(_)) => Err(primary),
            (Ok(_), Err(perf)) => Err(perf),
            (Err(primary), Err(perf)) => {
                Err(format!("{primary}; Performance report failed: {perf}"))
            }
        }
    }

    fn finish(self, success: bool) -> Result<(), String> {
        let mode = match &self.backend {
            Backend::Periodic { .. } => "periodic",
            Backend::Boundary(_) => "boundary_only",
        };
        let collected = match self.backend {
            Backend::Periodic { tx, handle } => {
                tx.send(Control::Stop)
                    .map_err(|e| format!("Cannot stop performance sampling thread: {e}"))?;
                handle
                    .join()
                    .map_err(|_| "Performance sampling thread exited unexpectedly".to_string())??
            }
            Backend::Boundary(boundary) => {
                let BoundaryBackend { state } = *boundary;
                let BoundaryState {
                    mut system,
                    pid,
                    mut collector,
                } = state.into_inner().map_err(|_| {
                    "Performance boundary-sampling state lock is poisoned".to_string()
                })?;
                match refresh(&mut system, pid) {
                    Ok(raw) => collector.observe(raw),
                    Err(_) => collector.sample_error(),
                }
                collector.finish(Instant::now())
            }
        };
        let metadata = ReportMetadata {
            command: self.command,
            label: &self.label,
            status: if success { "success" } else { "error" },
            configured_threads: self.configured_threads,
            work_thread_budget: self.configured_threads,
            sampling_mode: mode,
            interval_ms: self.interval.as_millis() as u64,
            logical_cpu_count: self.logical_cpu_count,
            physical_core_count: self.physical_core_count,
        };
        write_report(&self.path, &metadata, &collected)?;
        Ok(())
    }
}

fn sampler_loop(
    pid: Pid,
    interval: Duration,
    rx: mpsc::Receiver<Control>,
    ready: mpsc::SyncSender<Result<usize, String>>,
) -> Result<Collected, String> {
    let mut system = System::new();
    let initial = match refresh(&mut system, pid) {
        Ok(raw) => raw,
        Err(e) => {
            let _ = ready.send(Err(e.clone()));
            return Err(e);
        }
    };
    let mut collector = Collector::new(initial);
    if ready.send(Ok(system.cpus().len())).is_err() {
        return Err("Performance monitor initialization receiver is closed".to_string());
    }

    loop {
        match rx.recv_timeout(interval) {
            Ok(Control::Stage(name, at)) => collector.stage(name, at),
            Ok(Control::Stop) => {
                match refresh(&mut system, pid) {
                    Ok(raw) => collector.observe(raw),
                    Err(_) => collector.sample_error(),
                }
                return Ok(collector.finish(Instant::now()));
            }
            Err(mpsc::RecvTimeoutError::Timeout) => match refresh(&mut system, pid) {
                Ok(raw) => collector.observe(raw),
                Err(_) => collector.sample_error(),
            },
            Err(mpsc::RecvTimeoutError::Disconnected) => {
                return Err("Performance monitor control channel closed early".to_string());
            }
        }
    }
}

fn refresh(system: &mut System, pid: Pid) -> Result<RawMetrics, String> {
    system.refresh_memory();
    system.refresh_cpu_usage();
    let pids = [pid];
    let kind = ProcessRefreshKind::nothing()
        .with_cpu()
        .with_memory()
        .with_disk_usage()
        .without_tasks();
    system.refresh_processes_specifics(ProcessesToUpdate::Some(&pids), false, kind);
    let process = system.process(pid).ok_or_else(|| {
        format!(
            "Performance monitor cannot find current process PID {}",
            pid.as_u32()
        )
    })?;
    let disk = process.disk_usage();
    Ok(RawMetrics {
        at: Instant::now(),
        process_cpu_pct: process.cpu_usage(),
        process_cpu_time_ms: process.accumulated_cpu_time(),
        rss_bytes: process.memory(),
        virtual_bytes: process.virtual_memory(),
        read_bytes_total: disk.total_read_bytes,
        written_bytes_total: disk.total_written_bytes,
        system_cpu_pct: system.global_cpu_usage(),
        system_memory_total_bytes: system.total_memory(),
        system_memory_available_bytes: system.available_memory(),
    })
}

struct ReportMetadata<'a> {
    command: &'static str,
    label: &'a str,
    status: &'static str,
    configured_threads: usize,
    work_thread_budget: usize,
    sampling_mode: &'static str,
    interval_ms: u64,
    logical_cpu_count: usize,
    physical_core_count: Option<usize>,
}

fn fmt_opt(value: Option<f64>) -> String {
    value.map_or_else(|| "null".to_string(), |v| format!("{v:.3}"))
}

fn write_report(path: &Path, meta: &ReportMetadata<'_>, c: &Collected) -> Result<(), String> {
    let mut json = String::new();
    json.push_str("{\n");
    json.push_str(&format!("  \"schema\": \"{SCHEMA}\",\n"));
    json.push_str(&format!(
        "  \"command\": \"{}\",\n  \"label\": \"{}\",\n  \"status\": \"{}\",\n",
        meta.command,
        json_escape(meta.label),
        meta.status
    ));
    json.push_str(&format!(
        "  \"run\": {{\"thread_budget\": {}, \"work_thread_budget\": {}, \"sampling_mode\": \"{}\", \"interval_ms\": {}, \"wall_time_ms\": {}, \"sample_count\": {}, \"valid_cpu_samples\": {}}},\n",
        meta.configured_threads,
        meta.work_thread_budget,
        meta.sampling_mode,
        meta.interval_ms,
        c.wall_time_ms,
        c.sample_count,
        c.valid_cpu_samples
    ));
    json.push_str(&format!(
        "  \"process\": {{\"pid\": {}, \"cpu_pct_mean\": {}, \"cpu_pct_peak\": {}, \"cpu_time_ms\": {}, \"rss_bytes_peak\": {}, \"virtual_bytes_peak\": {}, \"read_bytes\": {}, \"written_bytes\": {}}},\n",
        std::process::id(),
        fmt_opt(c.process_cpu_pct_mean),
        fmt_opt(c.process_cpu_pct_peak),
        c.process_cpu_time_ms,
        c.rss_bytes_peak,
        c.virtual_bytes_peak,
        c.read_bytes,
        c.written_bytes
    ));
    json.push_str(&format!(
        "  \"system\": {{\"logical_cpu_count\": {}, \"physical_core_count\": {}, \"cpu_pct_mean\": {}, \"cpu_pct_peak\": {}, \"memory_total_bytes\": {}, \"memory_available_bytes_min\": {}}},\n",
        meta.logical_cpu_count,
        meta.physical_core_count.map_or_else(|| "null".to_string(), |v| v.to_string()),
        fmt_opt(c.system_cpu_pct_mean),
        fmt_opt(c.system_cpu_pct_peak),
        c.system_memory_total_bytes,
        c.system_memory_available_bytes_min
    ));
    json.push_str(&format!(
        "  \"collection\": {{\"sysinfo_version\": \"{}\", \"platform\": \"{}\", \"sample_errors\": {}, \"partial\": {}}},\n",
        SYSINFO_VERSION,
        std::env::consts::OS,
        c.sample_errors,
        c.sample_errors > 0
    ));
    json.push_str("  \"stages\": [\n");
    for (i, stage) in c.stages.iter().enumerate() {
        json.push_str(&format!(
            "    {{\"name\": \"{}\", \"wall_time_ms\": {}, \"cpu_pct_mean\": {}, \"cpu_pct_peak\": {}, \"rss_bytes_peak\": {}}}{}\n",
            stage.name,
            stage.wall_time_ms(c.finished),
            fmt_opt(stage.process_cpu.mean()),
            fmt_opt(stage.process_cpu.peak()),
            stage.rss_bytes_peak,
            if i + 1 < c.stages.len() { "," } else { "" }
        ));
    }
    json.push_str("  ]\n}\n");

    write_report_file(path, json.as_bytes())
}

fn write_report_file(path: &Path, json: &[u8]) -> Result<(), String> {
    let suffix = format!(".part.{}", std::process::id());
    let json_tmp = PathBuf::from(format!("{}{}", path.display(), suffix));
    let result = (|| {
        write_new(&json_tmp, json)?;
        std::fs::rename(&json_tmp, path).map_err(|e| {
            format!(
                "Failed to finalize performance JSON {}: {e}",
                path.display()
            )
        })?;
        Ok(())
    })();
    if result.is_err() {
        let _ = std::fs::remove_file(&json_tmp);
    }
    result
}

fn write_new(path: &Path, bytes: &[u8]) -> Result<(), String> {
    let mut file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(path)
        .map_err(|e| {
            format!(
                "Cannot create temporary performance-report file {}: {e}",
                path.display()
            )
        })?;
    file.write_all(bytes)
        .and_then(|_| file.flush())
        .map_err(|e| format!("Failed to write performance report {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn raw(at: Instant, cpu: f32, cpu_ms: u64, rss: u64) -> RawMetrics {
        RawMetrics {
            at,
            process_cpu_pct: cpu,
            process_cpu_time_ms: cpu_ms,
            rss_bytes: rss,
            virtual_bytes: rss * 2,
            read_bytes_total: cpu_ms * 3,
            written_bytes_total: cpu_ms * 2,
            system_cpu_pct: cpu / 10.0,
            system_memory_total_bytes: 1_000_000,
            system_memory_available_bytes: 900_000 - rss,
        }
    }

    #[test]
    fn collector_tracks_peaks_deltas_and_stages() {
        let start = Instant::now();
        let mut c = Collector::new(raw(start, 0.0, 100, 10));
        c.observe(raw(start + Duration::from_millis(100), 900.0, 200, 20));
        c.observe(raw(start + Duration::from_secs(1), 200.0, 1_100, 50));
        c.stage("map", start + Duration::from_secs(1));
        c.observe(raw(start + Duration::from_secs(2), 300.0, 2_100, 80));
        let got = c.finish(start + Duration::from_secs(2));
        assert_eq!(got.wall_time_ms, 2_000);
        assert_eq!(got.process_cpu_time_ms, 2_000);
        assert_eq!(got.process_cpu_pct_mean, Some(100.0));
        assert_eq!(got.process_cpu_pct_peak, Some(300.0));
        assert_eq!(got.valid_cpu_samples, 2);
        assert_eq!(got.rss_bytes_peak, 80);
        assert_eq!(got.read_bytes, 6_000);
        assert_eq!(got.stages.len(), 2);
        assert_eq!(got.stages[0].wall_time_ms(got.finished), 1_000);
        assert_eq!(got.stages[1].name, "map");
    }

    #[test]
    fn periodic_monitor_stops_immediately_and_writes_reports() {
        let dir = std::env::temp_dir().join(format!("vf_perf_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("sample");
        let monitor = PerfMonitor::start("run", "sample".to_string(), &prefix, 2).unwrap();
        assert_eq!(monitor.work_thread_budget(), 2);
        monitor.stage("test");
        let started = Instant::now();
        monitor.complete(Ok(())).unwrap();
        assert!(started.elapsed() < SAMPLE_INTERVAL);
        let path = report_path(&prefix);
        let json = std::fs::read_to_string(path).unwrap();
        assert!(json.contains("\"schema\": \"viroflash.perf.v1\""));
        assert!(json.contains("\"sampling_mode\": \"periodic\""));
        assert!(json.contains("\"name\": \"test\""));
        std::fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn one_thread_uses_boundary_sampling() {
        let dir = std::env::temp_dir().join(format!("vf_perf_one_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("sample");
        let monitor = PerfMonitor::start("run", "sample".to_string(), &prefix, 1).unwrap();
        assert_eq!(monitor.work_thread_budget(), 1);
        monitor.stage("test");
        monitor.complete(Ok(())).unwrap();
        let json = std::fs::read_to_string(report_path(&prefix)).unwrap();
        assert!(json.contains("\"sampling_mode\": \"boundary_only\""));
        assert!(json.contains("\"sample_count\": 3"));
        assert!(json.contains("\"name\": \"test\""));
        std::fs::remove_dir_all(dir).unwrap();
    }
}
