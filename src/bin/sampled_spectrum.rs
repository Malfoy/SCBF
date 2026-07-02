use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufWriter, Write},
    path::{Path, PathBuf},
    process::ExitCode,
    time::Instant,
};

use clap::Parser;
use super_counting_bloom::{
    CounterBits, Result, SuperCountingBloomConfig, SuperCountingBloomError,
    SuperCountingBloomIndex, for_each_fastx_record,
};

#[derive(Parser)]
#[command(
    author,
    version,
    about = "Estimate an abundance spectrum by querying sampled inserted k-mers"
)]
struct Cli {
    /// FASTA/FASTQ input used to build the abundance index and sample inserted k-mers.
    input: PathBuf,

    /// Output CSV path. Defaults to stdout.
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of inserted k-mer windows sampled after indexing.
    #[arg(long, default_value_t = 1_000_000)]
    sample_size: usize,

    /// Deterministic reservoir-sampling seed.
    #[arg(long, default_value_t = 1)]
    seed: u64,

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

struct SampledKmers {
    kmers: Vec<Vec<u8>>,
    total_inserted_kmers: u64,
}

struct SampledSpectrum {
    rows: Vec<(u64, u64)>,
    queried_kmers: u64,
    invalid_queries: u64,
}

struct ReservoirSampler {
    kmers: Vec<Vec<u8>>,
    capacity: usize,
    seen: u64,
    rng: SplitMix64,
}

#[derive(Clone, Copy)]
struct SplitMix64 {
    state: u64,
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
    if cli.sample_size == 0 {
        return Err(SuperCountingBloomError::InvalidConfig(
            "sample_size must be positive".to_string(),
        ));
    }

    let cfg = config_from_cli(&cli);
    let build = SuperCountingBloomIndex::build_from_fastx(&cli.input, cfg)?;
    eprintln!(
        "indexed_kmers: {}\tinserted_smers: {}\tcounter_slots: {}\tblocks: {}\tbuild_s: {:.3}",
        build.indexed_kmers,
        build.inserted_smers,
        build.counter_slots,
        build.blocks,
        build.build_seconds
    );

    let sample_start = Instant::now();
    let sampled = sample_inserted_kmers(
        &cli.input,
        build.index.config().k,
        cli.sample_size,
        cli.seed,
    )?;
    let sample_seconds = sample_start.elapsed().as_secs_f64();
    if sampled.total_inserted_kmers != build.indexed_kmers {
        eprintln!(
            "warning: sampled valid k-mer windows ({}) differ from indexed k-mers ({})",
            sampled.total_inserted_kmers, build.indexed_kmers
        );
    }

    let query_start = Instant::now();
    let spectrum = query_sampled_kmers(&build.index, &sampled.kmers);
    let query_seconds = query_start.elapsed().as_secs_f64();
    write_scaled_spectrum(cli.output.as_ref(), &sampled, &spectrum)?;

    eprintln!(
        "total_inserted_kmers: {}\tsampled_kmers: {}\tqueried_kmers: {}\tinvalid_queries: {}\tsample_s: {:.3}\tquery_s: {:.3}",
        sampled.total_inserted_kmers,
        sampled.kmers.len(),
        spectrum.queried_kmers,
        spectrum.invalid_queries,
        sample_seconds,
        query_seconds
    );

    Ok(())
}

fn sample_inserted_kmers(
    path: &Path,
    k: usize,
    sample_size: usize,
    seed: u64,
) -> Result<SampledKmers> {
    let mut sampler = ReservoirSampler::new(sample_size, seed);
    for_each_fastx_record(path, |_, sequence| {
        sample_sequence_kmers(sequence, k, &mut sampler);
        Ok(())
    })?;
    Ok(sampler.finish())
}

fn sample_sequence_kmers(sequence: &[u8], k: usize, sampler: &mut ReservoirSampler) {
    if sequence.len() < k {
        return;
    }

    let mut run_start = None;
    for (idx, &base) in sequence.iter().enumerate() {
        if is_acgt(base) {
            if run_start.is_none() {
                run_start = Some(idx);
            }
        } else if let Some(start) = run_start.take() {
            sample_valid_run(sequence, start, idx, k, sampler);
        }
    }
    if let Some(start) = run_start {
        sample_valid_run(sequence, start, sequence.len(), k, sampler);
    }
}

fn sample_valid_run(
    sequence: &[u8],
    start: usize,
    end: usize,
    k: usize,
    sampler: &mut ReservoirSampler,
) {
    if end - start < k {
        return;
    }
    for offset in start..=end - k {
        sampler.consider(&sequence[offset..offset + k]);
    }
}

fn query_sampled_kmers(
    index: &SuperCountingBloomIndex,
    sampled_kmers: &[Vec<u8>],
) -> SampledSpectrum {
    let mut histogram = HashMap::<u64, u64>::new();
    let mut queried_kmers = 0_u64;
    let mut invalid_queries = 0_u64;

    for kmer in sampled_kmers {
        queried_kmers += 1;
        let abundances = index.estimate_sequence_abundances(kmer);
        match abundances.first().copied().flatten() {
            Some(abundance) => *histogram.entry(abundance).or_insert(0) += 1,
            None => invalid_queries += 1,
        }
    }

    let mut rows: Vec<_> = histogram.into_iter().collect();
    rows.sort_unstable_by_key(|&(abundance, _)| abundance);
    SampledSpectrum {
        rows,
        queried_kmers,
        invalid_queries,
    }
}

fn write_scaled_spectrum(
    path: Option<&PathBuf>,
    sampled: &SampledKmers,
    spectrum: &SampledSpectrum,
) -> Result<()> {
    let mut writer = open_writer(path)?;
    writeln!(
        writer,
        "estimated_count,sampled_kmers,linear_scaled_windows,estimated_distinct_kmers"
    )
    .map_err(stream_error)?;

    let sample_count = sampled.kmers.len().max(1) as f64;
    let scale = sampled.total_inserted_kmers as f64 / sample_count;
    for &(abundance, sampled_kmers) in &spectrum.rows {
        let linear_scaled_windows = (sampled_kmers as f64 * scale).round() as u64;
        let estimated_distinct_kmers = if abundance == 0 {
            0
        } else {
            (linear_scaled_windows as f64 / abundance as f64).round() as u64
        };
        writeln!(
            writer,
            "{abundance},{sampled_kmers},{linear_scaled_windows},{estimated_distinct_kmers}"
        )
        .map_err(stream_error)?;
    }
    Ok(())
}

impl ReservoirSampler {
    fn new(capacity: usize, seed: u64) -> Self {
        Self {
            kmers: Vec::with_capacity(capacity),
            capacity,
            seen: 0,
            rng: SplitMix64::new(seed),
        }
    }

    fn consider(&mut self, kmer: &[u8]) {
        self.seen += 1;
        let slot = if self.kmers.len() < self.capacity {
            Some(self.kmers.len())
        } else {
            let replacement = self.rng.gen_below(self.seen);
            (replacement < self.capacity as u64).then_some(replacement as usize)
        };

        if let Some(slot) = slot {
            let retained: Vec<_> = kmer.iter().map(|base| base.to_ascii_uppercase()).collect();
            if slot == self.kmers.len() {
                self.kmers.push(retained);
            } else {
                self.kmers[slot] = retained;
            }
        }
    }

    fn finish(self) -> SampledKmers {
        SampledKmers {
            kmers: self.kmers,
            total_inserted_kmers: self.seen,
        }
    }
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut value = self.state;
        value = (value ^ (value >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        value = (value ^ (value >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        value ^ (value >> 31)
    }

    fn gen_below(&mut self, upper: u64) -> u64 {
        self.next_u64() % upper
    }
}

fn is_acgt(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
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
            |err| SuperCountingBloomError::Io {
                path: path.display().to_string(),
                message: err.to_string(),
            },
        )?))),
        None => Ok(Box::new(BufWriter::new(io::stdout().lock()))),
    }
}

fn stream_error(err: io::Error) -> SuperCountingBloomError {
    SuperCountingBloomError::Io {
        path: "<stream>".to_string(),
        message: err.to_string(),
    }
}
