use std::{
    fs::File,
    io::{self, BufWriter, Write},
    path::PathBuf,
    process::ExitCode,
};

use clap::Parser;
use super_counting_bloom::{
    CounterBits, Result, SuperCountingBloomConfig, SuperCountingBloomIndex,
};

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Experimental filter-only k-mer spectrum estimator"
)]
struct Cli {
    /// FASTA/FASTQ input used to build the abundance index.
    input: PathBuf,

    /// Output CSV path. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Optional raw counter histogram CSV path.
    #[arg(long)]
    counter_histogram: Option<PathBuf>,

    /// Maximum abundance emitted by the filter-only estimator.
    #[arg(long, default_value_t = 255)]
    max_abundance: u64,

    /// Override the k-mer count scaling. Defaults to positive counters divided by hash probes.
    #[arg(long)]
    kmer_count_hint: Option<u64>,

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
    let build = SuperCountingBloomIndex::build_from_fastx(&cli.input, config_from_cli(&cli))?;
    eprintln!(
        "warning: filter-only spectra are model-based counter-histogram estimates, not a replacement for key-aware spectra"
    );
    eprintln!(
        "indexed_kmers: {}\tinserted_smers: {}\tcounter_slots: {}\tblocks: {}\tbuild_s: {:.3}",
        build.indexed_kmers,
        build.inserted_smers,
        build.counter_slots,
        build.blocks,
        build.build_seconds
    );

    let counter_histogram = build.index.counter_histogram();
    if let Some(path) = cli.counter_histogram.as_ref() {
        write_counter_histogram(path, &counter_histogram)?;
    }

    let kmer_count_hint = cli.kmer_count_hint.unwrap_or_else(|| {
        build
            .index
            .filter_only_count_hint_from_counter_histogram(&counter_histogram, cli.max_abundance)
    });
    let estimate = build
        .index
        .estimate_filter_only_spectrum_from_counter_histogram(
            &counter_histogram,
            cli.max_abundance,
            kmer_count_hint,
        );
    let mut writer = open_writer(cli.output.as_ref())?;
    writeln!(
        writer,
        "estimated_count,estimated_kmers,method,kmer_count_hint"
    )
    .map_err(stream_error)?;
    for &(count, kmers) in &estimate.rows {
        writeln!(
            writer,
            "{count},{kmers},{},{}",
            estimate.method, estimate.kmer_count_hint
        )
        .map_err(stream_error)?;
    }
    Ok(())
}

fn write_counter_histogram(path: &PathBuf, rows: &[(u64, u64)]) -> Result<()> {
    let mut writer = BufWriter::new(File::create(path).map_err(|err| {
        super_counting_bloom::SuperCountingBloomError::Io {
            path: path.display().to_string(),
            message: err.to_string(),
        }
    })?);
    writeln!(writer, "counter_value,counters").map_err(stream_error)?;
    for &(value, counters) in rows {
        writeln!(writer, "{value},{counters}").map_err(stream_error)?;
    }
    Ok(())
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

fn stream_error(err: io::Error) -> super_counting_bloom::SuperCountingBloomError {
    super_counting_bloom::SuperCountingBloomError::Io {
        path: "<stream>".to_string(),
        message: err.to_string(),
    }
}
