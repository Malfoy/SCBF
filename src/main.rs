use std::{path::PathBuf, process::ExitCode};

use clap::Parser;
use super_counting_bloom::{
    CounterBits, SuperCountingBloomConfig, estimate_spectrum_from_fastx_pair,
    stream_estimates_from_fastx_pair, write_spectrum,
};

#[derive(Parser)]
#[command(author, version, about)]
struct Cli {
    /// FASTA/FASTQ input.
    input: PathBuf,

    /// Optional FASTA/FASTQ query input. Defaults to the indexed input.
    #[arg(long)]
    query_input: Option<PathBuf>,

    /// K-mer size.
    #[arg(short = 'k', long, default_value_t = 31)]
    k: usize,

    /// Minimizer length used to group super-k-mers.
    #[arg(short = 'm', long, default_value_t = 21)]
    m: usize,

    /// Fimpera/findere subword length.
    #[arg(short = 's', long, default_value_t = 25)]
    s: usize,

    /// Number of counter probes per s-mer.
    #[arg(short = 'H', long, default_value_t = 4)]
    n_hashes: usize,

    /// Counter width in bits: 8, 16, or 32.
    #[arg(long, default_value = "16", value_parser = parse_counter_bits)]
    counter_bits: CounterBits,

    /// Total number of counters is 2^counter-slots-exponent.
    #[arg(long, default_value_t = 30)]
    counter_slots_exponent: u8,

    /// Counters per SuperCounting Bloom block is 2^block-slots-exponent.
    #[arg(long, default_value_t = 9)]
    block_slots_exponent: u8,

    /// Worker threads used for build and query passes.
    #[arg(short = 't', long)]
    threads: Option<usize>,

    /// Number of parsed chunks buffered between parser and workers.
    #[arg(long, default_value_t = 4096)]
    queue: usize,

    /// Build the index and stream abundance estimates without building a spectrum.
    #[arg(long)]
    stream_query_only: bool,

    /// Write approximate spectrum TSV to this path. Defaults to stdout.
    #[arg(long)]
    spectrum_out: Option<PathBuf>,
}

fn parse_counter_bits(value: &str) -> std::result::Result<CounterBits, String> {
    match value {
        "8" => Ok(CounterBits::Bits8),
        "16" => Ok(CounterBits::Bits16),
        "32" => Ok(CounterBits::Bits32),
        _ => Err("expected one of: 8, 16, 32".to_string()),
    }
}

fn main() -> ExitCode {
    let cli = Cli::parse();
    let cfg = SuperCountingBloomConfig {
        k: cli.k,
        m: cli.m,
        s: cli.s,
        n_hashes: cli.n_hashes,
        counter_bits: cli.counter_bits,
        counter_slots_exponent: cli.counter_slots_exponent,
        block_slots_exponent: cli.block_slots_exponent,
        threads: cli.threads.unwrap_or_else(num_cpus::get).max(1),
        queue: cli.queue,
    };

    let query_input = cli.query_input.as_deref().unwrap_or(&cli.input);
    match run(
        &cli.input,
        query_input,
        cli.spectrum_out.as_deref(),
        cli.stream_query_only,
        &cfg,
    ) {
        Ok(()) => ExitCode::SUCCESS,
        Err(err) => {
            eprintln!("error: {err}");
            ExitCode::FAILURE
        }
    }
}

fn run(
    input: &std::path::Path,
    query_input: &std::path::Path,
    spectrum_out: Option<&std::path::Path>,
    stream_query_only: bool,
    cfg: &SuperCountingBloomConfig,
) -> super_counting_bloom::Result<()> {
    if stream_query_only {
        let report = stream_estimates_from_fastx_pair(input, query_input, cfg)?;
        eprintln!(
            "chunks: {}\tindexed_kmers: {}\tqueried_kmers: {}\tpositive_estimates: {}\tinserted_smers: {}\testimate_checksum: {}\tcounter_bits: {}\tcounter_slots: {}\tblocks: {}\tbuild_s: {:.3}\tquery_s: {:.3}\ttotal_s: {:.3}",
            report.stats.chunks,
            report.stats.indexed_kmers,
            report.stats.queried_kmers,
            report.stats.positive_estimates,
            report.stats.inserted_smers,
            report.stats.estimate_checksum,
            cfg.counter_bits.as_u8(),
            report.stats.counter_slots,
            report.stats.blocks,
            report.stats.build_seconds,
            report.stats.query_seconds,
            report.stats.total_seconds
        );
        return Ok(());
    }

    let report = estimate_spectrum_from_fastx_pair(input, query_input, cfg)?;
    write_spectrum(&report.spectrum, spectrum_out)?;
    eprintln!(
        "chunks: {}\tindexed_kmers: {}\tqueried_kmers: {}\tinserted_smers: {}\tdistinct_keys: {}\tcounter_bits: {}\tcounter_slots: {}\tblocks: {}\tbuild_s: {:.3}\tquery_s: {:.3}\tspectrum_s: {:.3}\ttotal_s: {:.3}",
        report.stats.chunks,
        report.stats.indexed_kmers,
        report.stats.queried_kmers,
        report.stats.inserted_smers,
        report.stats.distinct_keys,
        cfg.counter_bits.as_u8(),
        report.stats.counter_slots,
        report.stats.blocks,
        report.stats.build_seconds,
        report.stats.query_seconds,
        report.stats.spectrum_seconds,
        report.stats.total_seconds
    );
    Ok(())
}
