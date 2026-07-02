use std::{collections::HashMap, path::Path, sync::Arc, thread, time::Instant};

use crossbeam_channel::bounded;
use helicase::{
    Config, FastxParser, HelicaseParser, ParserOptions, dna_format::PackedDNA, input::FromFile,
    parser::Event,
};

use crate::{
    ApproxSpectrum, Counter, CounterBits, FrozenSuperCountingBloom, Result, SuperCountingBloom,
    SuperCountingBloomConfig, SuperCountingBloomError,
};

const PACKED_FASTX_CONFIG: Config = ParserOptions::default()
    .ignore_headers()
    .dna_packed()
    .return_record(false)
    .config();

#[derive(Debug, Clone, Copy, Default)]
pub struct SuperCountingBloomStats {
    pub chunks: u64,
    pub indexed_kmers: u64,
    pub inserted_smers: u64,
    pub queried_kmers: u64,
    pub distinct_keys: usize,
    pub counter_slots: usize,
    pub blocks: usize,
    pub build_seconds: f64,
    pub query_seconds: f64,
    pub spectrum_seconds: f64,
    pub total_seconds: f64,
}

#[derive(Debug)]
pub struct SuperCountingBloomReport {
    pub spectrum: ApproxSpectrum,
    pub stats: SuperCountingBloomStats,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct StreamEstimateStats {
    pub chunks: u64,
    pub indexed_kmers: u64,
    pub inserted_smers: u64,
    pub queried_kmers: u64,
    pub positive_estimates: u64,
    pub estimate_checksum: u64,
    pub counter_slots: usize,
    pub blocks: usize,
    pub build_seconds: f64,
    pub query_seconds: f64,
    pub total_seconds: f64,
}

#[derive(Debug)]
pub struct StreamEstimateReport {
    pub stats: StreamEstimateStats,
}

#[derive(Default)]
pub(crate) struct BuildStats {
    pub(crate) chunks: u64,
    pub(crate) indexed_kmers: u64,
    pub(crate) inserted_smers: u64,
}

struct QueryEstimates<C: Counter> {
    estimates: HashMap<u64, C>,
    queried_kmers: u64,
}

#[derive(Default)]
struct StreamQueryStats {
    queried_kmers: u64,
    positive_estimates: u64,
    estimate_checksum: u64,
}

pub fn estimate_spectrum_from_fastx<P: AsRef<Path>>(
    path: P,
    cfg: &SuperCountingBloomConfig,
) -> Result<SuperCountingBloomReport> {
    estimate_spectrum_from_fastx_pair(path.as_ref(), path.as_ref(), cfg)
}

pub fn estimate_spectrum_from_fastx_pair<I: AsRef<Path>, Q: AsRef<Path>>(
    index_path: I,
    query_path: Q,
    cfg: &SuperCountingBloomConfig,
) -> Result<SuperCountingBloomReport> {
    cfg.validate()?;
    let index_path = index_path.as_ref();
    let query_path = query_path.as_ref();
    match cfg.counter_bits {
        CounterBits::Bits8 => estimate_spectrum_from_fastx_typed::<u8>(index_path, query_path, cfg),
        CounterBits::Bits16 => {
            estimate_spectrum_from_fastx_typed::<u16>(index_path, query_path, cfg)
        }
        CounterBits::Bits32 => {
            estimate_spectrum_from_fastx_typed::<u32>(index_path, query_path, cfg)
        }
    }
}

pub fn stream_estimates_from_fastx<P: AsRef<Path>>(
    path: P,
    cfg: &SuperCountingBloomConfig,
) -> Result<StreamEstimateReport> {
    stream_estimates_from_fastx_pair(path.as_ref(), path.as_ref(), cfg)
}

pub fn stream_estimates_from_fastx_pair<I: AsRef<Path>, Q: AsRef<Path>>(
    index_path: I,
    query_path: Q,
    cfg: &SuperCountingBloomConfig,
) -> Result<StreamEstimateReport> {
    cfg.validate()?;
    let index_path = index_path.as_ref();
    let query_path = query_path.as_ref();
    match cfg.counter_bits {
        CounterBits::Bits8 => stream_estimates_from_fastx_typed::<u8>(index_path, query_path, cfg),
        CounterBits::Bits16 => {
            stream_estimates_from_fastx_typed::<u16>(index_path, query_path, cfg)
        }
        CounterBits::Bits32 => {
            stream_estimates_from_fastx_typed::<u32>(index_path, query_path, cfg)
        }
    }
}

fn estimate_spectrum_from_fastx_typed<C: Counter>(
    index_path: &Path,
    query_path: &Path,
    cfg: &SuperCountingBloomConfig,
) -> Result<SuperCountingBloomReport> {
    let total_start = Instant::now();
    let build_start = Instant::now();
    let (filter, build) = build_frozen_filter_typed::<C>(index_path, cfg)?;
    let build_seconds = build_start.elapsed().as_secs_f64();
    let counter_slots = filter.counter_slots();
    let blocks = filter.blocks();

    let query_start = Instant::now();
    let query = estimate_distinct_kmers(query_path, cfg, Arc::new(filter))?;
    let query_seconds = query_start.elapsed().as_secs_f64();

    let spectrum_start = Instant::now();
    let spectrum = spectrum_from_estimates(&query.estimates);
    let spectrum_seconds = spectrum_start.elapsed().as_secs_f64();
    let total_seconds = total_start.elapsed().as_secs_f64();

    Ok(SuperCountingBloomReport {
        spectrum,
        stats: SuperCountingBloomStats {
            chunks: build.chunks,
            indexed_kmers: build.indexed_kmers,
            inserted_smers: build.inserted_smers,
            queried_kmers: query.queried_kmers,
            distinct_keys: query.estimates.len(),
            counter_slots,
            blocks,
            build_seconds,
            query_seconds,
            spectrum_seconds,
            total_seconds,
        },
    })
}

fn stream_estimates_from_fastx_typed<C: Counter>(
    index_path: &Path,
    query_path: &Path,
    cfg: &SuperCountingBloomConfig,
) -> Result<StreamEstimateReport> {
    let total_start = Instant::now();
    let build_start = Instant::now();
    let (filter, build) = build_frozen_filter_typed::<C>(index_path, cfg)?;
    let build_seconds = build_start.elapsed().as_secs_f64();
    let counter_slots = filter.counter_slots();
    let blocks = filter.blocks();

    let query_start = Instant::now();
    let query = stream_query_kmers(query_path, cfg, Arc::new(filter))?;
    let query_seconds = query_start.elapsed().as_secs_f64();
    let total_seconds = total_start.elapsed().as_secs_f64();

    Ok(StreamEstimateReport {
        stats: StreamEstimateStats {
            chunks: build.chunks,
            indexed_kmers: build.indexed_kmers,
            inserted_smers: build.inserted_smers,
            queried_kmers: query.queried_kmers,
            positive_estimates: query.positive_estimates,
            estimate_checksum: query.estimate_checksum,
            counter_slots,
            blocks,
            build_seconds,
            query_seconds,
            total_seconds,
        },
    })
}

fn build_filter<C: Counter>(
    path: &Path,
    cfg: &SuperCountingBloomConfig,
    filter: Arc<SuperCountingBloom<C>>,
) -> Result<BuildStats> {
    let (sender, receiver) = bounded::<PackedDNA>(cfg.queue);
    let workers: Vec<_> = (0..cfg.threads)
        .map(|_| {
            let receiver = receiver.clone();
            let filter = Arc::clone(&filter);
            let cfg = cfg.clone();
            thread::spawn(move || {
                let mut inserted_smers = 0_u64;
                while let Ok(sequence) = receiver.recv() {
                    inserted_smers += filter.insert_packed_dna(sequence, &cfg);
                }
                inserted_smers
            })
        })
        .collect();

    let mut parser = open_fastx(path)?;
    let mut stats = BuildStats::default();

    while let Some(event) = parser.next() {
        if matches!(event, Event::DnaChunk(_)) && parser.get_dna_len() >= cfg.k {
            let sequence = parser.get_dna_packed_owned();
            stats.chunks += 1;
            stats.indexed_kmers += (sequence.len() + 1 - cfg.k) as u64;
            sender
                .send(sequence)
                .map_err(|err| SuperCountingBloomError::ChannelClosed(err.to_string()))?;
        }
    }
    drop(sender);

    for worker in workers {
        stats.inserted_smers += worker
            .join()
            .map_err(|_| SuperCountingBloomError::WorkerPanic("build"))?;
    }

    Ok(stats)
}

fn estimate_distinct_kmers<C: Counter>(
    path: &Path,
    cfg: &SuperCountingBloomConfig,
    filter: Arc<FrozenSuperCountingBloom<C>>,
) -> Result<QueryEstimates<C>> {
    let (sender, receiver) = bounded::<PackedDNA>(cfg.queue);
    let workers: Vec<_> = (0..cfg.threads)
        .map(|_| {
            let receiver = receiver.clone();
            let filter = Arc::clone(&filter);
            let cfg = cfg.clone();
            thread::spawn(move || {
                let mut estimates = HashMap::new();
                let mut queried_kmers = 0_u64;
                while let Ok(sequence) = receiver.recv() {
                    queried_kmers += filter.estimate_packed_dna(sequence, &cfg, &mut estimates);
                }
                (estimates, queried_kmers)
            })
        })
        .collect();

    let mut parser = open_fastx(path)?;

    while let Some(event) = parser.next() {
        if matches!(event, Event::DnaChunk(_)) && parser.get_dna_len() >= cfg.k {
            let sequence = parser.get_dna_packed_owned();
            sender
                .send(sequence)
                .map_err(|err| SuperCountingBloomError::ChannelClosed(err.to_string()))?;
        }
    }
    drop(sender);

    let mut estimates = HashMap::new();
    let mut queried_kmers = 0_u64;
    for worker in workers {
        let (local_estimates, local_queried) = worker
            .join()
            .map_err(|_| SuperCountingBloomError::WorkerPanic("query"))?;
        queried_kmers += local_queried;
        merge_estimates(&mut estimates, local_estimates);
    }

    Ok(QueryEstimates {
        estimates,
        queried_kmers,
    })
}

fn stream_query_kmers<C: Counter>(
    path: &Path,
    cfg: &SuperCountingBloomConfig,
    filter: Arc<FrozenSuperCountingBloom<C>>,
) -> Result<StreamQueryStats> {
    let (sender, receiver) = bounded::<PackedDNA>(cfg.queue);
    let workers: Vec<_> = (0..cfg.threads)
        .map(|_| {
            let receiver = receiver.clone();
            let filter = Arc::clone(&filter);
            let cfg = cfg.clone();
            thread::spawn(move || {
                let mut stats = StreamQueryStats::default();
                while let Ok(sequence) = receiver.recv() {
                    let abundances = filter.estimate_abundances_packed_dna(sequence, &cfg);
                    stats.queried_kmers += abundances.len() as u64;
                    for estimate in abundances {
                        if estimate.to_u64() > 0 {
                            stats.positive_estimates += 1;
                        }
                        stats.estimate_checksum =
                            stats.estimate_checksum.wrapping_add(estimate.to_u64());
                    }
                }
                stats
            })
        })
        .collect();

    let mut parser = open_fastx(path)?;

    while let Some(event) = parser.next() {
        if matches!(event, Event::DnaChunk(_)) && parser.get_dna_len() >= cfg.k {
            let sequence = parser.get_dna_packed_owned();
            sender
                .send(sequence)
                .map_err(|err| SuperCountingBloomError::ChannelClosed(err.to_string()))?;
        }
    }
    drop(sender);

    let mut stats = StreamQueryStats::default();
    for worker in workers {
        let local = worker
            .join()
            .map_err(|_| SuperCountingBloomError::WorkerPanic("stream query"))?;
        stats.queried_kmers += local.queried_kmers;
        stats.positive_estimates += local.positive_estimates;
        stats.estimate_checksum = stats
            .estimate_checksum
            .wrapping_add(local.estimate_checksum);
    }

    Ok(stats)
}

fn open_fastx(path: &Path) -> Result<FastxParser<'_, PACKED_FASTX_CONFIG>> {
    FastxParser::<PACKED_FASTX_CONFIG>::from_file(path).map_err(|err| {
        SuperCountingBloomError::FastxOpen {
            path: path.display().to_string(),
            message: err.to_string(),
        }
    })
}

pub(crate) fn build_frozen_filter_typed<C: Counter>(
    path: &Path,
    cfg: &SuperCountingBloomConfig,
) -> Result<(FrozenSuperCountingBloom<C>, BuildStats)> {
    let filter = Arc::new(SuperCountingBloom::<C>::new(cfg)?);
    let build = build_filter(path, cfg, Arc::clone(&filter))?;
    let filter = Arc::try_unwrap(filter)
        .map_err(|_| {
            SuperCountingBloomError::InternalState("filter is still shared after build".to_string())
        })?
        .freeze();
    Ok((filter, build))
}

fn merge_estimates<C: Counter>(dst: &mut HashMap<u64, C>, src: HashMap<u64, C>) {
    let reserve = src
        .len()
        .saturating_sub(dst.capacity().saturating_sub(dst.len()));
    dst.reserve(reserve);
    for (key, estimate) in src {
        dst.entry(key).or_insert(estimate);
    }
}

fn spectrum_from_estimates<C: Counter>(estimates: &HashMap<u64, C>) -> ApproxSpectrum {
    let mut spectrum = HashMap::<u64, u64>::new();
    for &count in estimates.values() {
        *spectrum.entry(count.to_u64()).or_insert(0) += 1;
    }
    let mut rows: Vec<_> = spectrum.into_iter().collect();
    rows.sort_unstable_by_key(|&(count, _)| count);
    rows
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!(
            "super_counting_bloom_{name}_{unique}_{}",
            std::process::id()
        ));
        path
    }

    fn small_cfg(bits: CounterBits) -> SuperCountingBloomConfig {
        SuperCountingBloomConfig {
            k: 9,
            m: 5,
            s: 7,
            n_hashes: 3,
            counter_bits: bits,
            counter_slots_exponent: 16,
            block_slots_exponent: 6,
            threads: 2,
            queue: 2,
        }
    }

    #[test]
    fn fastx_spectrum_smoke_for_all_counter_widths() {
        let path = temp_path("reads.fa");
        fs::write(
            &path,
            b">r1\nACGTACGTACGTACGTACGT\n>r2\nACGTACGTACGTACGTACGT\n",
        )
        .unwrap();

        for bits in [CounterBits::Bits8, CounterBits::Bits16, CounterBits::Bits32] {
            let report = estimate_spectrum_from_fastx(&path, &small_cfg(bits)).unwrap();
            assert!(!report.spectrum.is_empty());
            assert_eq!(report.stats.chunks, 2);
            assert_eq!(report.stats.indexed_kmers, 24);
            assert_eq!(report.stats.queried_kmers, 24);
        }

        fs::remove_file(path).unwrap();
    }
}
