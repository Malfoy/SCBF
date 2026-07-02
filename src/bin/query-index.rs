use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::PathBuf,
    process::ExitCode,
};

use clap::{Parser, ValueEnum};
use super_counting_bloom::{
    Result, SuperCountingBloomError, SuperCountingBloomIndex, for_each_fastx_record,
};

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Load a SuperCountingBloom index and query FASTA/FASTQ"
)]
struct Cli {
    /// Saved SuperCountingBloom index.
    index: PathBuf,

    /// FASTA/FASTQ queried against the saved index.
    query: PathBuf,

    /// Output CSV path. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Report writes one aggregate row. Fast writes one row per sequence. Slow writes one row per k-mer window.
    #[arg(long, value_enum, default_value_t = Mode::Fast)]
    mode: Mode,
}

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
    Report,
    Fast,
    Slow,
}

fn main() -> ExitCode {
    match run(Cli::parse()) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run(cli: Cli) -> Result<()> {
    let index = SuperCountingBloomIndex::load(&cli.index)?;
    eprintln!(
        "inserted_kmers: {}\tcounter_bits: {}\tthreads: {}",
        index.inserted_kmers(),
        index.config().counter_bits.as_u8(),
        index.config().threads
    );

    let mut writer = open_writer(cli.output.as_ref())?;
    match cli.mode {
        Mode::Report => write_report(&index, &cli.query, &mut writer),
        Mode::Fast => write_fast(&index, &cli.query, &mut writer),
        Mode::Slow => write_slow(&index, &cli.query, &mut writer),
    }
}

fn write_report<W: Write>(
    index: &SuperCountingBloomIndex,
    query: &PathBuf,
    writer: &mut W,
) -> Result<()> {
    let report = index.query_fasta(query)?;
    writeln!(
        writer,
        "records_processed,total_windows,valid_kmers,positive_estimates,estimate_checksum"
    )
    .map_err(stream_error)?;
    writeln!(
        writer,
        "{},{},{},{},{}",
        report.records_processed,
        report.total_windows,
        report.valid_kmers,
        report.positive_estimates,
        report.estimate_checksum
    )
    .map_err(stream_error)
}

fn write_fast<W: Write>(
    index: &SuperCountingBloomIndex,
    query: &PathBuf,
    writer: &mut W,
) -> Result<()> {
    writeln!(
        writer,
        "sequence,length,total_windows,valid_kmers,mean_abundance,median_abundance"
    )
    .map_err(stream_error)?;
    for_each_fastx_record(query, |name, sequence| {
        let summary = index.summarize_sequence(name, sequence);
        writeln!(
            writer,
            "{},{},{},{},{},{}",
            csv_escape(&summary.name),
            summary.sequence_len,
            summary.total_windows,
            summary.valid_kmers,
            format_optional(summary.mean_abundance),
            format_optional(summary.median_abundance)
        )
        .map_err(stream_error)
    })
}

fn write_slow<W: Write>(
    index: &SuperCountingBloomIndex,
    query: &PathBuf,
    writer: &mut W,
) -> Result<()> {
    writeln!(writer, "sequence,kmer_index,abundance").map_err(stream_error)?;
    for_each_fastx_record(query, |name, sequence| {
        let abundances = index.estimate_sequence_abundances(sequence);
        let name = csv_escape(name);
        for (idx, abundance) in abundances.into_iter().enumerate() {
            match abundance {
                Some(value) => writeln!(writer, "{name},{idx},{value}"),
                None => writeln!(writer, "{name},{idx},"),
            }
            .map_err(stream_error)?;
        }
        Ok(())
    })
}

fn open_writer(path: Option<&PathBuf>) -> Result<Box<dyn Write>> {
    match path {
        Some(path) => Ok(Box::new(BufWriter::new(File::create(path).map_err(
            |err| SuperCountingBloomError::Io {
                path: path.display().to_string(),
                message: err.to_string(),
            },
        )?))),
        None => Ok(Box::new(BufWriter::new(io::stdout().lock()))),
    }
}

fn format_optional(value: Option<f64>) -> String {
    value.map(|value| format!("{value:.6}")).unwrap_or_default()
}

fn csv_escape(value: &str) -> String {
    if value.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", value.replace('"', "\"\""))
    } else {
        value.to_string()
    }
}

fn stream_error(err: io::Error) -> SuperCountingBloomError {
    SuperCountingBloomError::Io {
        path: "<stream>".to_string(),
        message: err.to_string(),
    }
}
