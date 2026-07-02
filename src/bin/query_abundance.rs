use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::PathBuf,
    process::ExitCode,
};

use clap::{Parser, ValueEnum};
use super_counting_bloom::{
    CounterBits, Result, SuperCountingBloomConfig, SuperCountingBloomIndex, for_each_fastx_record,
};

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Index reads and query per-sequence k-mer abundances"
)]
struct Cli {
    /// FASTA/FASTQ used to build the abundance index.
    index: PathBuf,

    /// FASTA/FASTQ queried against the index.
    query: PathBuf,

    /// Output CSV path. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Fast mode writes one row per sequence. Slow mode writes one row per k-mer window.
    #[arg(long, value_enum, default_value_t = Mode::Fast)]
    mode: Mode,

    #[arg(short = 'k', long, default_value_t = 31)]
    k: usize,

    #[arg(short = 'm', long, default_value_t = 21)]
    m: usize,

    #[arg(short = 's', long, default_value_t = 25)]
    s: usize,

    #[arg(short = 'H', long, default_value_t = 4)]
    n_hashes: usize,

    #[arg(long, default_value = "16", value_parser = parse_counter_bits)]
    counter_bits: CounterBits,

    #[arg(long, default_value_t = 30)]
    counter_slots_exponent: u8,

    #[arg(long, default_value_t = 9)]
    block_slots_exponent: u8,

    #[arg(short = 't', long)]
    threads: Option<usize>,

    #[arg(long, default_value_t = 4096)]
    queue: usize,
}

#[derive(Clone, Copy, ValueEnum)]
enum Mode {
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
    let cfg = config_from_cli(&cli);
    let build = SuperCountingBloomIndex::build_from_fastx(&cli.index, cfg)?;
    eprintln!(
        "indexed_kmers: {}\tinserted_smers: {}\tcounter_slots: {}\tblocks: {}\tbuild_s: {:.3}",
        build.indexed_kmers,
        build.inserted_smers,
        build.counter_slots,
        build.blocks,
        build.build_seconds
    );

    let mut writer = open_writer(cli.output.as_ref())?;
    match cli.mode {
        Mode::Fast => write_fast(&build.index, &cli.query, &mut writer),
        Mode::Slow => write_slow(&build.index, &cli.query, &mut writer),
    }
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

fn config_from_cli(cli: &Cli) -> SuperCountingBloomConfig {
    SuperCountingBloomConfig {
        k: cli.k,
        m: cli.m,
        s: cli.s,
        n_hashes: cli.n_hashes,
        counter_bits: cli.counter_bits,
        counter_slots_exponent: cli.counter_slots_exponent,
        block_slots_exponent: cli.block_slots_exponent,
        threads: cli.threads.unwrap_or_else(num_cpus::get).max(1),
        queue: cli.queue,
    }
}

fn parse_counter_bits(value: &str) -> std::result::Result<CounterBits, String> {
    match value {
        "8" => Ok(CounterBits::Bits8),
        "16" => Ok(CounterBits::Bits16),
        "32" => Ok(CounterBits::Bits32),
        _ => Err("expected one of: 8, 16, 32".to_string()),
    }
}

fn open_writer(path: Option<&PathBuf>) -> Result<Box<dyn Write>> {
    match path {
        Some(path) => Ok(Box::new(BufWriter::new(File::create(path).map_err(
            |err| super_counting_bloom::SuperCountingBloomError::Io {
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

fn stream_error(err: io::Error) -> super_counting_bloom::SuperCountingBloomError {
    super_counting_bloom::SuperCountingBloomError::Io {
        path: "<stream>".to_string(),
        message: err.to_string(),
    }
}
