use std::{path::PathBuf, process::ExitCode};

use clap::Parser;
use super_counting_bloom::{
    CounterBits, Result, SuperCountingBloomBuilder, SuperCountingBloomConfig,
};

#[derive(Parser)]
#[command(author, version, about = "Build and save a SuperCountingBloom index")]
struct Cli {
    /// FASTA/FASTQ input used to build the abundance index.
    input: PathBuf,

    /// Output index path.
    output: PathBuf,

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
    let mut builder = SuperCountingBloomBuilder::new(config_from_cli(&cli))?;
    let report = builder.add_fasta(&cli.input)?;
    let inserted_kmers = builder.inserted_kmers();
    let inserted_smers = builder.inserted_smers();
    let index = builder.into_index()?;
    index.save(&cli.output)?;
    eprintln!(
        "records_processed: {}\trecords_indexed: {}\tinserted_kmers: {}\tinserted_smers: {}\tcounter_bits: {}\tthreads: {}\tadd_s: {:.3}",
        report.records_processed,
        report.records_indexed,
        inserted_kmers,
        inserted_smers,
        index.config().counter_bits.as_u8(),
        index.config().threads,
        report.add_seconds
    );
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
