mod alignment;
mod evidence;
mod fastq;
mod gate;
mod index;
mod integration_evidence;
mod output;
mod pipeline;
mod profile;
mod report;
mod sampling;
mod telemetry;

pub use index::{build_index, IndexOptions};
pub use pipeline::{run_pipeline, RunOptions, RunSummary};
