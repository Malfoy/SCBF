use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufReader, BufWriter, Read, Write},
    path::Path,
    sync::Arc,
    thread,
    time::Instant,
};

use crossbeam_channel::{Sender, bounded};
use helicase::{
    Config, FastxParser, HelicaseParser, ParserOptions, input::FromFile, parser::Event,
};
use packed_seq::{PackedSeqVec, SeqVec};

use crate::{
    ApproxSpectrum, Counter, CounterBits, FrozenSuperCountingBloom, Result, SuperCountingBloom,
    SuperCountingBloomConfig, SuperCountingBloomError,
    fastx::{BuildStats, build_frozen_filter_typed},
    filter::FrozenFilterLayout,
};

const RECORD_FASTX_CONFIG: Config = ParserOptions::default()
    .dna_string()
    .keep_non_actg()
    .config();
const INDEX_MAGIC: &[u8; 8] = b"SCBIDX\0\0";
const INDEX_FORMAT_VERSION: u32 = 2;

enum TypedIndex {
    Bits8(FrozenSuperCountingBloom<u8>),
    Bits16(FrozenSuperCountingBloom<u16>),
    Bits32(FrozenSuperCountingBloom<u32>),
}

enum MutableTypedIndex {
    Bits8(Arc<SuperCountingBloom<u8>>),
    Bits16(Arc<SuperCountingBloom<u16>>),
    Bits32(Arc<SuperCountingBloom<u32>>),
}

pub struct SuperCountingBloomBuilder {
    cfg: SuperCountingBloomConfig,
    inner: MutableTypedIndex,
    inserted_kmers: u64,
    inserted_smers: u64,
}

pub struct SuperCountingBloomIndex {
    cfg: SuperCountingBloomConfig,
    inner: TypedIndex,
    inserted_kmers: u64,
}

pub struct BuildIndexReport {
    pub index: SuperCountingBloomIndex,
    pub chunks: u64,
    pub indexed_kmers: u64,
    pub inserted_smers: u64,
    pub counter_slots: usize,
    pub blocks: usize,
    pub build_seconds: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct AddReport {
    pub records_processed: u64,
    pub records_indexed: u64,
    pub inserted_kmers: u64,
    pub inserted_smers: u64,
    pub add_seconds: f64,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct QueryFastaReport {
    pub records_processed: u64,
    pub total_windows: u64,
    pub valid_kmers: u64,
    pub positive_estimates: u64,
    pub estimate_checksum: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SequenceAbundanceSummary {
    pub name: String,
    pub sequence_len: usize,
    pub total_windows: usize,
    pub valid_kmers: usize,
    pub mean_abundance: Option<f64>,
    pub median_abundance: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct FilterOnlySpectrumEstimate {
    pub method: &'static str,
    pub kmer_count_hint: u64,
    pub rows: ApproxSpectrum,
}

impl SuperCountingBloomBuilder {
    pub fn new(cfg: SuperCountingBloomConfig) -> Result<Self> {
        cfg.validate()?;
        let inner = match cfg.counter_bits {
            CounterBits::Bits8 => {
                MutableTypedIndex::Bits8(Arc::new(SuperCountingBloom::new(&cfg)?))
            }
            CounterBits::Bits16 => {
                MutableTypedIndex::Bits16(Arc::new(SuperCountingBloom::new(&cfg)?))
            }
            CounterBits::Bits32 => {
                MutableTypedIndex::Bits32(Arc::new(SuperCountingBloom::new(&cfg)?))
            }
        };
        Ok(Self {
            cfg,
            inner,
            inserted_kmers: 0,
            inserted_smers: 0,
        })
    }

    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<u64> {
        let (inserted_kmers, inserted_smers) = match &self.inner {
            MutableTypedIndex::Bits8(filter) => {
                insert_ascii_sequence_typed(filter, &self.cfg, sequence)
            }
            MutableTypedIndex::Bits16(filter) => {
                insert_ascii_sequence_typed(filter, &self.cfg, sequence)
            }
            MutableTypedIndex::Bits32(filter) => {
                insert_ascii_sequence_typed(filter, &self.cfg, sequence)
            }
        };
        self.inserted_kmers = self.inserted_kmers.saturating_add(inserted_kmers);
        self.inserted_smers = self.inserted_smers.saturating_add(inserted_smers);
        Ok(inserted_kmers)
    }

    pub fn add_fasta<P: AsRef<Path>>(&mut self, path: P) -> Result<AddReport> {
        self.cfg.validate()?;
        let start = Instant::now();
        let mut report = match &self.inner {
            MutableTypedIndex::Bits8(filter) => {
                add_fasta_typed(Arc::clone(filter), &self.cfg, path.as_ref())?
            }
            MutableTypedIndex::Bits16(filter) => {
                add_fasta_typed(Arc::clone(filter), &self.cfg, path.as_ref())?
            }
            MutableTypedIndex::Bits32(filter) => {
                add_fasta_typed(Arc::clone(filter), &self.cfg, path.as_ref())?
            }
        };
        report.add_seconds = start.elapsed().as_secs_f64();
        self.inserted_kmers = self.inserted_kmers.saturating_add(report.inserted_kmers);
        self.inserted_smers = self.inserted_smers.saturating_add(report.inserted_smers);
        Ok(report)
    }

    pub fn into_index(self) -> Result<SuperCountingBloomIndex> {
        let inner = match self.inner {
            MutableTypedIndex::Bits8(filter) => {
                TypedIndex::Bits8(unwrap_mutable_filter(filter)?.freeze())
            }
            MutableTypedIndex::Bits16(filter) => {
                TypedIndex::Bits16(unwrap_mutable_filter(filter)?.freeze())
            }
            MutableTypedIndex::Bits32(filter) => {
                TypedIndex::Bits32(unwrap_mutable_filter(filter)?.freeze())
            }
        };
        Ok(SuperCountingBloomIndex {
            cfg: self.cfg,
            inner,
            inserted_kmers: self.inserted_kmers,
        })
    }

    pub fn inserted_kmers(&self) -> u64 {
        self.inserted_kmers
    }

    pub fn inserted_smers(&self) -> u64 {
        self.inserted_smers
    }

    pub fn set_threads(&mut self, threads: usize) -> Result<()> {
        if threads == 0 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "threads must be positive".to_string(),
            ));
        }
        self.cfg.threads = threads;
        Ok(())
    }

    pub fn clear_threads(&mut self) {
        self.cfg.threads = num_cpus::get().max(1);
    }

    pub fn threads(&self) -> usize {
        self.cfg.threads
    }

    pub fn config(&self) -> &SuperCountingBloomConfig {
        &self.cfg
    }
}

impl SuperCountingBloomIndex {
    pub fn build_from_fastx<P: AsRef<Path>>(
        path: P,
        cfg: SuperCountingBloomConfig,
    ) -> Result<BuildIndexReport> {
        cfg.validate()?;
        let path = path.as_ref();
        let start = Instant::now();
        match cfg.counter_bits {
            CounterBits::Bits8 => build_index_typed::<u8>(path, cfg, start),
            CounterBits::Bits16 => build_index_typed::<u16>(path, cfg, start),
            CounterBits::Bits32 => build_index_typed::<u32>(path, cfg, start),
        }
    }

    pub fn config(&self) -> &SuperCountingBloomConfig {
        &self.cfg
    }

    pub fn inserted_kmers(&self) -> u64 {
        self.inserted_kmers
    }

    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<()> {
        let path = path.as_ref();
        let path_name = path.display().to_string();
        let file = File::create(path).map_err(|err| io_path_error(path, err))?;
        let mut writer = BufWriter::new(file);
        write_index(self, &mut writer, &path_name)?;
        writer.flush().map_err(|err| io_name_error(&path_name, err))
    }

    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self> {
        let path = path.as_ref();
        let path_name = path.display().to_string();
        let file = File::open(path).map_err(|err| io_path_error(path, err))?;
        let mut reader = BufReader::new(file);
        read_index(&mut reader, &path_name)
    }

    pub fn estimate_sequence_abundances(&self, sequence: &[u8]) -> Vec<Option<u64>> {
        match &self.inner {
            TypedIndex::Bits8(filter) => {
                estimate_sequence_abundances_typed(filter, &self.cfg, sequence)
            }
            TypedIndex::Bits16(filter) => {
                estimate_sequence_abundances_typed(filter, &self.cfg, sequence)
            }
            TypedIndex::Bits32(filter) => {
                estimate_sequence_abundances_typed(filter, &self.cfg, sequence)
            }
        }
    }

    pub fn summarize_sequence(
        &self,
        name: impl Into<String>,
        sequence: &[u8],
    ) -> SequenceAbundanceSummary {
        let abundances = self.estimate_sequence_abundances(sequence);
        summarize_abundances(name.into(), sequence.len(), abundances)
    }

    pub fn query_fasta<P: AsRef<Path>>(&self, path: P) -> Result<QueryFastaReport> {
        let (sender, receiver) = bounded::<Vec<u8>>(self.cfg.queue);
        let mut worker_reports = Vec::with_capacity(self.cfg.threads);
        let mut parse_result = Ok(());

        thread::scope(|scope| {
            let mut handles = Vec::with_capacity(self.cfg.threads);
            for _ in 0..self.cfg.threads {
                let receiver = receiver.clone();
                handles.push(scope.spawn(move || {
                    let mut report = QueryFastaReport::default();
                    while let Ok(sequence) = receiver.recv() {
                        report.records_processed += 1;
                        let abundances = self.estimate_sequence_abundances(&sequence);
                        report.total_windows += abundances.len() as u64;
                        for abundance in abundances.into_iter().flatten() {
                            report.valid_kmers += 1;
                            if abundance > 0 {
                                report.positive_estimates += 1;
                            }
                            report.estimate_checksum =
                                report.estimate_checksum.wrapping_add(abundance);
                        }
                    }
                    report
                }));
            }
            drop(receiver);

            parse_result = send_fastx_records(path.as_ref(), &sender);
            drop(sender);

            for handle in handles {
                match handle.join() {
                    Ok(report) => worker_reports.push(report),
                    Err(_) => {
                        parse_result = Err(SuperCountingBloomError::WorkerPanic("query_fasta"));
                    }
                }
            }
        });

        parse_result?;
        let mut report = QueryFastaReport::default();
        for worker_report in worker_reports {
            report.records_processed += worker_report.records_processed;
            report.total_windows += worker_report.total_windows;
            report.valid_kmers += worker_report.valid_kmers;
            report.positive_estimates += worker_report.positive_estimates;
            report.estimate_checksum = report
                .estimate_checksum
                .wrapping_add(worker_report.estimate_checksum);
        }
        Ok(report)
    }

    pub fn counter_histogram(&self) -> ApproxSpectrum {
        match &self.inner {
            TypedIndex::Bits8(filter) => filter.counter_histogram(),
            TypedIndex::Bits16(filter) => filter.counter_histogram(),
            TypedIndex::Bits32(filter) => filter.counter_histogram(),
        }
    }

    pub fn estimate_filter_only_spectrum(
        &self,
        max_abundance: u64,
        kmer_count_hint: u64,
    ) -> FilterOnlySpectrumEstimate {
        let counter_histogram = self.counter_histogram();
        self.estimate_filter_only_spectrum_from_counter_histogram(
            &counter_histogram,
            max_abundance,
            kmer_count_hint,
        )
    }

    pub fn estimate_filter_only_spectrum_from_counter_histogram(
        &self,
        counter_histogram: &[(u64, u64)],
        max_abundance: u64,
        kmer_count_hint: u64,
    ) -> FilterOnlySpectrumEstimate {
        let rows = estimate_filter_only_spectrum_from_counter_histogram(
            counter_histogram,
            self.cfg.n_hashes,
            self.cfg.k - self.cfg.s + 1,
            max_abundance,
            kmer_count_hint,
        );
        FilterOnlySpectrumEstimate {
            method: "compound_poisson_counter_deconvolution",
            kmer_count_hint,
            rows,
        }
    }

    pub fn filter_only_count_hint(&self) -> u64 {
        let counter_histogram = self.counter_histogram();
        self.filter_only_count_hint_from_counter_histogram(&counter_histogram, 255)
    }

    pub fn filter_only_count_hint_from_counter_histogram(
        &self,
        counter_histogram: &[(u64, u64)],
        max_abundance: u64,
    ) -> u64 {
        deconvolved_counts(counter_histogram, self.cfg.n_hashes, max_abundance)
            .into_iter()
            .sum::<f64>()
            .round()
            .max(1.0) as u64
    }
}

fn write_index<W: Write>(
    index: &SuperCountingBloomIndex,
    writer: &mut W,
    path_name: &str,
) -> Result<()> {
    write_exact(writer, INDEX_MAGIC, path_name)?;
    write_u32(writer, INDEX_FORMAT_VERSION, path_name)?;
    write_config(writer, &index.cfg, path_name)?;
    write_u64(writer, index.inserted_kmers, path_name)?;
    let layout = match &index.inner {
        TypedIndex::Bits8(filter) => filter.serialized_layout(),
        TypedIndex::Bits16(filter) => filter.serialized_layout(),
        TypedIndex::Bits32(filter) => filter.serialized_layout(),
    };
    write_layout(writer, layout, path_name)?;
    match &index.inner {
        TypedIndex::Bits8(filter) => filter
            .write_counter_data(writer)
            .map_err(|err| io_name_error(path_name, err)),
        TypedIndex::Bits16(filter) => filter
            .write_counter_data(writer)
            .map_err(|err| io_name_error(path_name, err)),
        TypedIndex::Bits32(filter) => filter
            .write_counter_data(writer)
            .map_err(|err| io_name_error(path_name, err)),
    }
}

fn read_index<R: Read>(reader: &mut R, path_name: &str) -> Result<SuperCountingBloomIndex> {
    let mut magic = [0_u8; 8];
    read_exact(reader, &mut magic, path_name)?;
    if &magic != INDEX_MAGIC {
        return Err(SuperCountingBloomError::InvalidIndexFormat(
            "bad SuperCountingBloom index magic".to_string(),
        ));
    }

    let version = read_u32(reader, path_name)?;
    if version != INDEX_FORMAT_VERSION {
        return Err(SuperCountingBloomError::InvalidIndexFormat(format!(
            "unsupported index format version {version}"
        )));
    }

    let cfg = read_config(reader, path_name)?;
    cfg.validate()?;
    let inserted_kmers = read_u64(reader, path_name)?;
    let stored_layout = read_layout(reader, path_name)?;
    let expected_layout = FrozenFilterLayout::from_config(&cfg)?;
    if stored_layout != expected_layout {
        return Err(SuperCountingBloomError::InvalidIndexFormat(format!(
            "stored filter layout does not match config: stored {stored_layout:?}, expected {expected_layout:?}"
        )));
    }

    let inner = match cfg.counter_bits {
        CounterBits::Bits8 => TypedIndex::Bits8(
            FrozenSuperCountingBloom::<u8>::from_counter_data(reader, stored_layout)
                .map_err(|err| io_name_error(path_name, err))?,
        ),
        CounterBits::Bits16 => TypedIndex::Bits16(
            FrozenSuperCountingBloom::<u16>::from_counter_data(reader, stored_layout)
                .map_err(|err| io_name_error(path_name, err))?,
        ),
        CounterBits::Bits32 => TypedIndex::Bits32(
            FrozenSuperCountingBloom::<u32>::from_counter_data(reader, stored_layout)
                .map_err(|err| io_name_error(path_name, err))?,
        ),
    };

    reject_trailing_bytes(reader, path_name)?;
    Ok(SuperCountingBloomIndex {
        cfg,
        inner,
        inserted_kmers,
    })
}

fn write_config<W: Write>(
    writer: &mut W,
    cfg: &SuperCountingBloomConfig,
    path_name: &str,
) -> Result<()> {
    write_usize(writer, cfg.k, path_name)?;
    write_usize(writer, cfg.m, path_name)?;
    write_usize(writer, cfg.s, path_name)?;
    write_usize(writer, cfg.n_hashes, path_name)?;
    write_u8(writer, cfg.counter_bits.as_u8(), path_name)?;
    write_u8(writer, cfg.counter_slots_exponent, path_name)?;
    write_u8(writer, cfg.block_slots_exponent, path_name)?;
    write_usize(writer, cfg.threads, path_name)?;
    write_usize(writer, cfg.queue, path_name)
}

fn read_config<R: Read>(reader: &mut R, path_name: &str) -> Result<SuperCountingBloomConfig> {
    let k = read_usize(reader, path_name, "k")?;
    let m = read_usize(reader, path_name, "m")?;
    let s = read_usize(reader, path_name, "s")?;
    let n_hashes = read_usize(reader, path_name, "n_hashes")?;
    let counter_bits = match read_u8(reader, path_name)? {
        8 => CounterBits::Bits8,
        16 => CounterBits::Bits16,
        32 => CounterBits::Bits32,
        other => {
            return Err(SuperCountingBloomError::InvalidIndexFormat(format!(
                "unsupported counter width {other}"
            )));
        }
    };
    let counter_slots_exponent = read_u8(reader, path_name)?;
    let block_slots_exponent = read_u8(reader, path_name)?;
    let threads = read_usize(reader, path_name, "threads")?;
    let queue = read_usize(reader, path_name, "queue")?;

    Ok(SuperCountingBloomConfig {
        k,
        m,
        s,
        n_hashes,
        counter_bits,
        counter_slots_exponent,
        block_slots_exponent,
        threads,
        queue,
    })
}

fn write_layout<W: Write>(
    writer: &mut W,
    layout: FrozenFilterLayout,
    path_name: &str,
) -> Result<()> {
    write_usize(writer, layout.shard_count, path_name)?;
    write_usize(writer, layout.subblocks, path_name)?;
    write_usize(writer, layout.block_slots, path_name)?;
    write_usize(writer, layout.shard_mask, path_name)?;
    write_u32(writer, layout.shard_shift, path_name)?;
    write_usize(writer, layout.block_slots_mask, path_name)?;
    write_usize(writer, layout.nb_blocks, path_name)?;
    write_usize(writer, layout.block_mask, path_name)?;
    write_usize(writer, layout.n_hashes, path_name)
}

fn read_layout<R: Read>(reader: &mut R, path_name: &str) -> Result<FrozenFilterLayout> {
    Ok(FrozenFilterLayout {
        shard_count: read_usize(reader, path_name, "shard_count")?,
        subblocks: read_usize(reader, path_name, "subblocks")?,
        block_slots: read_usize(reader, path_name, "block_slots")?,
        shard_mask: read_usize(reader, path_name, "shard_mask")?,
        shard_shift: read_u32(reader, path_name)?,
        block_slots_mask: read_usize(reader, path_name, "block_slots_mask")?,
        nb_blocks: read_usize(reader, path_name, "nb_blocks")?,
        block_mask: read_usize(reader, path_name, "block_mask")?,
        n_hashes: read_usize(reader, path_name, "n_hashes")?,
    })
}

fn write_u8<W: Write>(writer: &mut W, value: u8, path_name: &str) -> Result<()> {
    write_exact(writer, &[value], path_name)
}

fn write_u32<W: Write>(writer: &mut W, value: u32, path_name: &str) -> Result<()> {
    write_exact(writer, &value.to_le_bytes(), path_name)
}

fn write_u64<W: Write>(writer: &mut W, value: u64, path_name: &str) -> Result<()> {
    write_exact(writer, &value.to_le_bytes(), path_name)
}

fn write_usize<W: Write>(writer: &mut W, value: usize, path_name: &str) -> Result<()> {
    write_exact(writer, &(value as u64).to_le_bytes(), path_name)
}

fn read_u8<R: Read>(reader: &mut R, path_name: &str) -> Result<u8> {
    let mut bytes = [0_u8; 1];
    read_exact(reader, &mut bytes, path_name)?;
    Ok(bytes[0])
}

fn read_u32<R: Read>(reader: &mut R, path_name: &str) -> Result<u32> {
    let mut bytes = [0_u8; 4];
    read_exact(reader, &mut bytes, path_name)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_u64<R: Read>(reader: &mut R, path_name: &str) -> Result<u64> {
    let mut bytes = [0_u8; 8];
    read_exact(reader, &mut bytes, path_name)?;
    Ok(u64::from_le_bytes(bytes))
}

fn read_usize<R: Read>(reader: &mut R, path_name: &str, field: &str) -> Result<usize> {
    let mut bytes = [0_u8; 8];
    read_exact(reader, &mut bytes, path_name)?;
    usize::try_from(u64::from_le_bytes(bytes)).map_err(|_| {
        SuperCountingBloomError::InvalidIndexFormat(format!("{field} does not fit in usize"))
    })
}

fn write_exact<W: Write>(writer: &mut W, bytes: &[u8], path_name: &str) -> Result<()> {
    writer
        .write_all(bytes)
        .map_err(|err| io_name_error(path_name, err))
}

fn read_exact<R: Read>(reader: &mut R, bytes: &mut [u8], path_name: &str) -> Result<()> {
    reader
        .read_exact(bytes)
        .map_err(|err| io_name_error(path_name, err))
}

fn reject_trailing_bytes<R: Read>(reader: &mut R, path_name: &str) -> Result<()> {
    let mut trailing = [0_u8; 1];
    match reader.read(&mut trailing) {
        Ok(0) => Ok(()),
        Ok(_) => Err(SuperCountingBloomError::InvalidIndexFormat(
            "trailing bytes after serialized index".to_string(),
        )),
        Err(err) => Err(io_name_error(path_name, err)),
    }
}

fn io_path_error(path: &Path, err: io::Error) -> SuperCountingBloomError {
    io_name_error(&path.display().to_string(), err)
}

fn io_name_error(path: &str, err: io::Error) -> SuperCountingBloomError {
    SuperCountingBloomError::Io {
        path: path.to_string(),
        message: err.to_string(),
    }
}

fn add_fasta_typed<C: Counter>(
    filter: Arc<SuperCountingBloom<C>>,
    cfg: &SuperCountingBloomConfig,
    path: &Path,
) -> Result<AddReport> {
    let (sender, receiver) = bounded::<Vec<u8>>(cfg.queue);
    let workers: Vec<_> = (0..cfg.threads)
        .map(|_| {
            let receiver = receiver.clone();
            let filter = Arc::clone(&filter);
            let cfg = cfg.clone();
            thread::spawn(move || {
                let mut report = AddReport::default();
                while let Ok(sequence) = receiver.recv() {
                    report.records_processed += 1;
                    let (inserted_kmers, inserted_smers) =
                        insert_ascii_sequence_typed(filter.as_ref(), &cfg, &sequence);
                    if inserted_kmers > 0 {
                        report.records_indexed += 1;
                    }
                    report.inserted_kmers += inserted_kmers;
                    report.inserted_smers += inserted_smers;
                }
                report
            })
        })
        .collect();
    drop(receiver);

    let send_result = send_fastx_records(path, &sender);
    drop(sender);

    let mut report = AddReport::default();
    for worker in workers {
        let local = worker
            .join()
            .map_err(|_| SuperCountingBloomError::WorkerPanic("add_fasta"))?;
        report.records_processed += local.records_processed;
        report.records_indexed += local.records_indexed;
        report.inserted_kmers += local.inserted_kmers;
        report.inserted_smers += local.inserted_smers;
    }

    send_result?;
    Ok(report)
}

fn send_fastx_records(path: &Path, sender: &Sender<Vec<u8>>) -> Result<()> {
    for_each_fastx_record(path, |_, sequence| {
        sender
            .send(sequence.to_vec())
            .map_err(|err| SuperCountingBloomError::ChannelClosed(err.to_string()))
    })
}

fn insert_ascii_sequence_typed<C: Counter>(
    filter: &SuperCountingBloom<C>,
    cfg: &SuperCountingBloomConfig,
    sequence: &[u8],
) -> (u64, u64) {
    if sequence.len() < cfg.k {
        return (0, 0);
    }

    let mut inserted_kmers = 0_u64;
    let mut inserted_smers = 0_u64;
    let mut run_start = None;
    for (idx, &base) in sequence.iter().enumerate() {
        if is_acgt(base) {
            if run_start.is_none() {
                run_start = Some(idx);
            }
        } else if let Some(start) = run_start.take() {
            let (kmers, smers) = insert_valid_run_typed(filter, cfg, sequence, start, idx);
            inserted_kmers += kmers;
            inserted_smers += smers;
        }
    }
    if let Some(start) = run_start {
        let (kmers, smers) = insert_valid_run_typed(filter, cfg, sequence, start, sequence.len());
        inserted_kmers += kmers;
        inserted_smers += smers;
    }

    (inserted_kmers, inserted_smers)
}

fn insert_valid_run_typed<C: Counter>(
    filter: &SuperCountingBloom<C>,
    cfg: &SuperCountingBloomConfig,
    sequence: &[u8],
    start: usize,
    end: usize,
) -> (u64, u64) {
    if end - start < cfg.k {
        return (0, 0);
    }
    let uppercase: Vec<_> = sequence[start..end]
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect();
    let packed = PackedSeqVec::from_ascii(&uppercase);
    let inserted_kmers = (uppercase.len() + 1 - cfg.k) as u64;
    let inserted_smers = filter.insert_packed_sequence(packed.as_slice(), cfg);
    (inserted_kmers, inserted_smers)
}

fn unwrap_mutable_filter<C: Counter>(
    filter: Arc<SuperCountingBloom<C>>,
) -> Result<SuperCountingBloom<C>> {
    Arc::try_unwrap(filter).map_err(|_| {
        SuperCountingBloomError::InternalState(
            "mutable filter is still shared while freezing".to_string(),
        )
    })
}

fn build_index_typed<C: Counter>(
    path: &Path,
    cfg: SuperCountingBloomConfig,
    start: Instant,
) -> Result<BuildIndexReport>
where
    TypedIndex: From<FrozenSuperCountingBloom<C>>,
{
    let (filter, stats) = build_frozen_filter_typed::<C>(path, &cfg)?;
    let counter_slots = filter.counter_slots();
    let blocks = filter.blocks();
    let index = SuperCountingBloomIndex {
        cfg,
        inner: TypedIndex::from(filter),
        inserted_kmers: stats.indexed_kmers,
    };
    Ok(report_from_build(
        index,
        stats,
        counter_slots,
        blocks,
        start,
    ))
}

fn report_from_build(
    index: SuperCountingBloomIndex,
    stats: BuildStats,
    counter_slots: usize,
    blocks: usize,
    start: Instant,
) -> BuildIndexReport {
    BuildIndexReport {
        index,
        chunks: stats.chunks,
        indexed_kmers: stats.indexed_kmers,
        inserted_smers: stats.inserted_smers,
        counter_slots,
        blocks,
        build_seconds: start.elapsed().as_secs_f64(),
    }
}

impl From<FrozenSuperCountingBloom<u8>> for TypedIndex {
    fn from(filter: FrozenSuperCountingBloom<u8>) -> Self {
        Self::Bits8(filter)
    }
}

impl From<FrozenSuperCountingBloom<u16>> for TypedIndex {
    fn from(filter: FrozenSuperCountingBloom<u16>) -> Self {
        Self::Bits16(filter)
    }
}

impl From<FrozenSuperCountingBloom<u32>> for TypedIndex {
    fn from(filter: FrozenSuperCountingBloom<u32>) -> Self {
        Self::Bits32(filter)
    }
}

fn estimate_sequence_abundances_typed<C: Counter>(
    filter: &FrozenSuperCountingBloom<C>,
    cfg: &SuperCountingBloomConfig,
    sequence: &[u8],
) -> Vec<Option<u64>> {
    if sequence.len() < cfg.k {
        return Vec::new();
    }

    let mut abundances = vec![None; sequence.len() + 1 - cfg.k];
    let mut run_start = None;
    for (idx, &base) in sequence.iter().enumerate() {
        if is_acgt(base) {
            if run_start.is_none() {
                run_start = Some(idx);
            }
        } else if let Some(start) = run_start.take() {
            fill_valid_run(filter, cfg, sequence, start, idx, &mut abundances);
        }
    }
    if let Some(start) = run_start {
        fill_valid_run(
            filter,
            cfg,
            sequence,
            start,
            sequence.len(),
            &mut abundances,
        );
    }

    abundances
}

fn fill_valid_run<C: Counter>(
    filter: &FrozenSuperCountingBloom<C>,
    cfg: &SuperCountingBloomConfig,
    sequence: &[u8],
    start: usize,
    end: usize,
    output: &mut [Option<u64>],
) {
    if end - start < cfg.k {
        return;
    }
    let uppercase: Vec<_> = sequence[start..end]
        .iter()
        .map(|base| base.to_ascii_uppercase())
        .collect();
    let packed = PackedSeqVec::from_ascii(&uppercase);
    let estimates = filter.estimate_abundances_packed_sequence(packed.as_slice(), cfg);
    for (offset, estimate) in estimates.into_iter().enumerate() {
        output[start + offset] = Some(estimate.to_u64());
    }
}

fn summarize_abundances(
    name: String,
    sequence_len: usize,
    abundances: Vec<Option<u64>>,
) -> SequenceAbundanceSummary {
    let total_windows = abundances.len();
    let mut valid: Vec<_> = abundances.into_iter().flatten().collect();
    let valid_kmers = valid.len();
    if valid.is_empty() {
        return SequenceAbundanceSummary {
            name,
            sequence_len,
            total_windows,
            valid_kmers,
            mean_abundance: None,
            median_abundance: None,
        };
    }

    let mean = valid.iter().sum::<u64>() as f64 / valid.len() as f64;
    valid.sort_unstable();
    let mid = valid.len() / 2;
    let median = if valid.len().is_multiple_of(2) {
        (valid[mid - 1] as f64 + valid[mid] as f64) / 2.0
    } else {
        valid[mid] as f64
    };

    SequenceAbundanceSummary {
        name,
        sequence_len,
        total_windows,
        valid_kmers,
        mean_abundance: Some(mean),
        median_abundance: Some(median),
    }
}

pub fn for_each_fastx_record<P, F>(path: P, mut visit: F) -> Result<()>
where
    P: AsRef<Path>,
    F: FnMut(&str, &[u8]) -> Result<()>,
{
    let path = path.as_ref();
    let mut parser = FastxParser::<RECORD_FASTX_CONFIG>::from_file(path).map_err(|err| {
        SuperCountingBloomError::FastxOpen {
            path: path.display().to_string(),
            message: err.to_string(),
        }
    })?;

    let mut record_idx = 0_u64;
    while let Some(event) = parser.next() {
        if matches!(event, Event::Record(_)) {
            record_idx += 1;
            let name = record_name(parser.get_header(), record_idx);
            visit(&name, parser.get_dna_string())?;
        }
    }
    Ok(())
}

fn record_name(header: &[u8], record_idx: u64) -> String {
    let trimmed = header
        .split(|byte| byte.is_ascii_whitespace())
        .next()
        .unwrap_or_default();
    if trimmed.is_empty() {
        format!("record_{record_idx}")
    } else {
        String::from_utf8_lossy(trimmed).into_owned()
    }
}

fn is_acgt(base: u8) -> bool {
    matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
}

fn estimate_filter_only_spectrum_from_counter_histogram(
    counter_histogram: &[(u64, u64)],
    n_hashes: usize,
    _window: usize,
    max_abundance: u64,
    kmer_count_hint: u64,
) -> ApproxSpectrum {
    let deconvolved = deconvolved_counts(counter_histogram, n_hashes, max_abundance);
    let deconvolved_total = deconvolved.iter().sum::<f64>();
    if deconvolved_total <= 0.0 {
        return scaled_positive_counter_histogram(
            counter_histogram,
            max_abundance,
            kmer_count_hint,
        );
    }
    let scale = kmer_count_hint as f64 / deconvolved_total;
    (1..=max_abundance)
        .map(|abundance| {
            let estimated = (deconvolved[(abundance - 1) as usize] * scale).round() as u64;
            (abundance, estimated)
        })
        .collect()
}

fn deconvolved_counts(
    counter_histogram: &[(u64, u64)],
    n_hashes: usize,
    max_abundance: u64,
) -> Vec<f64> {
    let max_abundance = max_abundance as usize;
    if max_abundance == 0 {
        return Vec::new();
    }

    let total_counters = counter_histogram
        .iter()
        .map(|&(_, counters)| counters)
        .sum::<u64>();
    if total_counters == 0 {
        return vec![0.0; max_abundance];
    }

    let mut probabilities = vec![0.0; max_abundance + 1];
    for &(abundance, counters) in counter_histogram {
        let abundance = abundance as usize;
        if abundance <= max_abundance {
            probabilities[abundance] += counters as f64 / total_counters as f64;
        }
    }

    let p0 = probabilities[0];
    if p0 <= 0.0 {
        return Vec::new();
    }

    let mut normalized = vec![0.0; max_abundance + 1];
    normalized[0] = 1.0;
    for abundance in 1..=max_abundance {
        normalized[abundance] = probabilities[abundance] / p0;
    }

    let mut log_coefficients = vec![0.0; max_abundance + 1];
    for abundance in 1..=max_abundance {
        let mut coefficient = abundance as f64 * normalized[abundance];
        for previous in 1..abundance {
            coefficient -=
                previous as f64 * log_coefficients[previous] * normalized[abundance - previous];
        }
        log_coefficients[abundance] = coefficient / abundance as f64;
    }

    let hash_count = n_hashes.max(1) as f64;
    log_coefficients
        .into_iter()
        .skip(1)
        .map(|coefficient| (total_counters as f64 * coefficient / hash_count).max(0.0))
        .collect()
}

fn scaled_positive_counter_histogram(
    counter_histogram: &[(u64, u64)],
    max_abundance: u64,
    kmer_count_hint: u64,
) -> ApproxSpectrum {
    let positive_counters = counter_histogram
        .iter()
        .filter(|&&(abundance, _)| abundance > 0)
        .map(|&(_, counters)| counters)
        .sum::<u64>()
        .max(1) as f64;
    let counter_counts: HashMap<_, _> = counter_histogram.iter().copied().collect();
    (1..=max_abundance)
        .map(|abundance| {
            let counters = counter_counts.get(&abundance).copied().unwrap_or(0);
            let estimated =
                (kmer_count_hint as f64 * counters as f64 / positive_counters).round() as u64;
            (abundance, estimated)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{
        fs,
        path::PathBuf,
        time::{SystemTime, UNIX_EPOCH},
    };

    fn small_cfg() -> SuperCountingBloomConfig {
        SuperCountingBloomConfig {
            k: 5,
            m: 3,
            s: 4,
            n_hashes: 2,
            counter_bits: CounterBits::Bits8,
            counter_slots_exponent: 12,
            block_slots_exponent: 5,
            threads: 1,
            queue: 2,
        }
    }

    #[test]
    fn sequence_abundances_preserve_invalid_windows() {
        let cfg = small_cfg();
        let path = test_fasta("invalid_windows.fa", b">r\nACGTACGTACGT\n");
        let index = SuperCountingBloomIndex::build_from_fastx(&path, cfg)
            .unwrap()
            .index;

        let abundances = index.estimate_sequence_abundances(b"ACGTACNNACGTAC");
        assert_eq!(abundances.len(), 10);
        assert!(abundances[..2].iter().all(Option::is_some));
        assert!(abundances[2..8].iter().all(Option::is_none));
        assert!(abundances[8..].iter().all(Option::is_some));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn fastx_record_iteration_preserves_record_names() {
        let path = test_fasta(
            "records.fa",
            b">contig_1 description\nACGTNN\n>contig_2\nTGCATG\n",
        );
        let mut names = Vec::new();
        let mut sequences = Vec::new();
        for_each_fastx_record(&path, |name, sequence| {
            names.push(name.to_string());
            sequences.push(sequence.to_vec());
            Ok(())
        })
        .unwrap();

        assert_eq!(names, ["contig_1", "contig_2"]);
        assert_eq!(sequences[0], b"ACGTNN");
        assert_eq!(sequences[1], b"TGCATG");
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn sequence_summary_reports_valid_windows() {
        let cfg = small_cfg();
        let path = test_fasta("summary.fa", b">r\nACGTACGTACGT\n");
        let index = SuperCountingBloomIndex::build_from_fastx(&path, cfg)
            .unwrap()
            .index;

        let summary = index.summarize_sequence("query", b"ACGTACGT");
        assert_eq!(summary.name, "query");
        assert_eq!(summary.sequence_len, 8);
        assert_eq!(summary.total_windows, 4);
        assert_eq!(summary.valid_kmers, 4);
        assert!(summary.mean_abundance.unwrap() > 0.0);
        assert!(summary.median_abundance.unwrap() > 0.0);
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn filter_only_spectrum_estimate_is_positive_abundance() {
        let cfg = small_cfg();
        let path = test_fasta("filter_spectrum.fa", b">r\nACGTACGTACGT\n");
        let index = SuperCountingBloomIndex::build_from_fastx(&path, cfg)
            .unwrap()
            .index;

        assert!(index.filter_only_count_hint() > 0);
        let spectrum = index.estimate_filter_only_spectrum(8, 10);
        assert_eq!(spectrum.method, "compound_poisson_counter_deconvolution");
        assert_eq!(spectrum.rows.first().map(|&(count, _)| count), Some(1));
        assert_eq!(spectrum.rows.last().map(|&(count, _)| count), Some(8));
        fs::remove_file(path).unwrap();
    }

    #[test]
    fn save_load_roundtrip_preserves_queries_for_all_counter_widths() {
        for bits in [CounterBits::Bits8, CounterBits::Bits16, CounterBits::Bits32] {
            let mut cfg = small_cfg();
            cfg.counter_bits = bits;
            let reads = test_fasta(
                &format!("roundtrip_reads_{}.fa", bits.as_u8()),
                b">r1\nACGTACGTACGTACGT\n>r2\nTTTTACGTACGTAAAA\n",
            );
            let index_path = temp_path(&format!("roundtrip_{}.scb", bits.as_u8()));
            let index = SuperCountingBloomIndex::build_from_fastx(&reads, cfg)
                .unwrap()
                .index;

            index.save(&index_path).unwrap();
            let loaded = SuperCountingBloomIndex::load(&index_path).unwrap();

            assert_eq!(loaded.config(), index.config());
            assert_eq!(loaded.counter_histogram(), index.counter_histogram());
            assert_eq!(
                loaded.estimate_sequence_abundances(b"ACGTACGTNNACGTACGT"),
                index.estimate_sequence_abundances(b"ACGTACGTNNACGTACGT")
            );

            fs::remove_file(reads).unwrap();
            fs::remove_file(index_path).unwrap();
        }
    }

    #[test]
    fn builder_add_sequence_add_fasta_query_and_save_load() {
        let mut cfg = small_cfg();
        cfg.counter_bits = CounterBits::Bits16;
        cfg.threads = 2;
        cfg.queue = 1;
        let reads = test_fasta(
            "builder_reads.fa",
            b">r1\nACGTACGTACGT\n>r2\nACGTNNACGTAC\n",
        );
        let index_path = temp_path("builder_index.scb");
        let mut builder = SuperCountingBloomBuilder::new(cfg).unwrap();

        assert_eq!(builder.threads(), 2);
        builder.set_threads(1).unwrap();
        assert_eq!(builder.threads(), 1);
        builder.clear_threads();
        assert!(builder.threads() >= 1);
        builder.set_threads(2).unwrap();

        assert_eq!(builder.add_sequence(b"ACGTACNNACGTAC").unwrap(), 4);
        let report = builder.add_fasta(&reads).unwrap();
        assert_eq!(report.records_processed, 2);
        assert_eq!(report.records_indexed, 2);
        assert_eq!(report.inserted_kmers, 10);
        assert_eq!(builder.inserted_kmers(), 14);
        assert!(builder.inserted_smers() > 0);

        let index = builder.into_index().unwrap();
        assert_eq!(index.inserted_kmers(), 14);
        let query = index.query_fasta(&reads).unwrap();
        assert_eq!(query.records_processed, 2);
        assert_eq!(query.total_windows, 16);
        assert_eq!(query.valid_kmers, 10);
        assert!(query.positive_estimates > 0);
        assert!(query.estimate_checksum > 0);

        index.save(&index_path).unwrap();
        let loaded = SuperCountingBloomIndex::load(&index_path).unwrap();
        assert_eq!(loaded.inserted_kmers(), 14);
        assert_eq!(loaded.query_fasta(&reads).unwrap(), query);

        fs::remove_file(reads).unwrap();
        fs::remove_file(index_path).unwrap();
    }

    #[test]
    fn builder_rejects_zero_threads() {
        let cfg = small_cfg();
        let mut builder = SuperCountingBloomBuilder::new(cfg).unwrap();
        match builder.set_threads(0) {
            Err(SuperCountingBloomError::InvalidConfig(_)) => {}
            Err(err) => panic!("expected invalid config, got {err}"),
            Ok(_) => panic!("expected invalid config"),
        }
    }

    #[test]
    fn load_rejects_bad_magic() {
        let path = temp_path("bad_magic.scb");
        fs::write(&path, b"not an index").unwrap();
        match SuperCountingBloomIndex::load(&path) {
            Err(SuperCountingBloomError::InvalidIndexFormat(_)) => {}
            Err(err) => panic!("expected invalid index format, got {err}"),
            Ok(_) => panic!("expected invalid index format"),
        }
        fs::remove_file(path).unwrap();
    }

    fn test_fasta(name: &str, contents: &[u8]) -> PathBuf {
        let path = temp_path(name);
        std::fs::write(&path, contents).unwrap();
        path
    }

    fn temp_path(name: &str) -> PathBuf {
        let mut path = std::env::temp_dir();
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        path.push(format!(
            "super_counting_bloom_api_{name}_{unique}_{}",
            std::process::id(),
        ));
        path
    }
}
