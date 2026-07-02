use std::{
    fs,
    path::PathBuf,
    time::{SystemTime, UNIX_EPOCH},
};

use packed_seq::SeqVec;
use super_counting_bloom::{
    Counter, CounterBits, SuperCountingBloom, SuperCountingBloomConfig,
    estimate_spectrum_from_fastx, estimate_spectrum_from_fastx_pair, stream_estimates_from_fastx,
    stream_estimates_from_fastx_pair, write_spectrum,
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

fn write_reads(name: &str) -> PathBuf {
    let path = temp_path(name);
    fs::write(
        &path,
        b">r1\nACGTACGTACGTACGTACGTACGTACGT\n>r2\nTGCATGCAACGTACGTTGCATGCAACGT\n",
    )
    .unwrap();
    path
}

fn write_negative_reads(name: &str) -> PathBuf {
    let path = temp_path(name);
    fs::write(
        &path,
        b">n1\nTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n>n2\nCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n",
    )
    .unwrap();
    path
}

fn cfg(bits: CounterBits, threads: usize, queue: usize) -> SuperCountingBloomConfig {
    SuperCountingBloomConfig {
        k: 9,
        m: 5,
        s: 7,
        n_hashes: 3,
        counter_bits: bits,
        counter_slots_exponent: 16,
        block_slots_exponent: 6,
        threads,
        queue,
    }
}

fn assert_fastx_spectrum(bits: CounterBits, threads: usize, queue: usize) {
    let path = write_reads("spectrum.fa");
    let cfg = cfg(bits, threads, queue);
    let report = estimate_spectrum_from_fastx(&path, &cfg).unwrap();
    assert_eq!(report.stats.chunks, 2);
    assert_eq!(report.stats.indexed_kmers, 40);
    assert_eq!(report.stats.queried_kmers, 40);
    assert!(!report.spectrum.is_empty());
    fs::remove_file(path).unwrap();
}

fn assert_fastx_stream(bits: CounterBits, threads: usize, queue: usize) {
    let path = write_reads("stream.fa");
    let cfg = cfg(bits, threads, queue);
    let report = stream_estimates_from_fastx(&path, &cfg).unwrap();
    assert_eq!(report.stats.chunks, 2);
    assert_eq!(report.stats.indexed_kmers, 40);
    assert_eq!(report.stats.queried_kmers, 40);
    assert_eq!(report.stats.positive_estimates, 40);
    assert!(report.stats.estimate_checksum > 0);
    fs::remove_file(path).unwrap();
}

fn assert_spectrum_writer(bits: CounterBits, threads: usize, queue: usize) {
    let reads = write_reads("write_spectrum.fa");
    let output = temp_path("spectrum.tsv");
    let cfg = cfg(bits, threads, queue);
    let report = estimate_spectrum_from_fastx(&reads, &cfg).unwrap();
    write_spectrum(&report.spectrum, Some(&output)).unwrap();
    let text = fs::read_to_string(&output).unwrap();
    assert!(text.starts_with("estimated_count\tkmers\n"));
    assert!(text.lines().count() > 1);
    fs::remove_file(reads).unwrap();
    fs::remove_file(output).unwrap();
}

fn assert_paired_query_paths(bits: CounterBits, threads: usize, queue: usize) {
    let index = write_reads("index.fa");
    let positive = write_reads("positive.fa");
    let negative = write_negative_reads("negative.fa");
    let cfg = cfg(bits, threads, queue);
    let positive_report = stream_estimates_from_fastx_pair(&index, &positive, &cfg).unwrap();
    let negative_report = stream_estimates_from_fastx_pair(&index, &negative, &cfg).unwrap();
    assert_eq!(
        positive_report.stats.indexed_kmers,
        negative_report.stats.indexed_kmers
    );
    assert_eq!(positive_report.stats.queried_kmers, 40);
    assert_eq!(negative_report.stats.queried_kmers, 40);
    assert_eq!(positive_report.stats.positive_estimates, 40);
    assert!(positive_report.stats.estimate_checksum >= negative_report.stats.estimate_checksum);
    let spectrum = estimate_spectrum_from_fastx_pair(&index, &positive, &cfg).unwrap();
    assert!(!spectrum.spectrum.is_empty());
    fs::remove_file(index).unwrap();
    fs::remove_file(positive).unwrap();
    fs::remove_file(negative).unwrap();
}

fn assert_saturation<C: Counter>(bits: CounterBits, repeats: usize, floor: u64) {
    let cfg = cfg(bits, 1, 2);
    let sequence = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGTACGTACGT");
    let bloom = SuperCountingBloom::<C>::new(&cfg).unwrap();
    for _ in 0..repeats {
        bloom.insert_packed_sequence(sequence.as_slice(), &cfg);
    }
    let bloom = bloom.freeze();
    let values = bloom.estimate_abundances_packed_sequence(sequence.as_slice(), &cfg);
    assert!(values.iter().all(|&value| value.to_u64() >= floor));
}

macro_rules! pipeline_case {
    ($name:ident, $bits:expr, $threads:expr, $queue:expr) => {
        mod $name {
            use super::*;

            #[test]
            fn fastx_spectrum() {
                assert_fastx_spectrum($bits, $threads, $queue);
            }

            #[test]
            fn fastx_stream() {
                assert_fastx_stream($bits, $threads, $queue);
            }

            #[test]
            fn write_spectrum_output() {
                assert_spectrum_writer($bits, $threads, $queue);
            }

            #[test]
            fn separate_index_and_query_paths() {
                assert_paired_query_paths($bits, $threads, $queue);
            }
        }
    };
}

pipeline_case!(bits8_threads1_queue1, CounterBits::Bits8, 1, 1);
pipeline_case!(bits16_threads1_queue1, CounterBits::Bits16, 1, 1);
pipeline_case!(bits32_threads1_queue1, CounterBits::Bits32, 1, 1);
pipeline_case!(bits8_threads2_queue2, CounterBits::Bits8, 2, 2);
pipeline_case!(bits16_threads2_queue2, CounterBits::Bits16, 2, 2);
pipeline_case!(bits32_threads2_queue2, CounterBits::Bits32, 2, 2);
pipeline_case!(bits8_threads4_queue8, CounterBits::Bits8, 4, 8);
pipeline_case!(bits16_threads4_queue8, CounterBits::Bits16, 4, 8);
pipeline_case!(bits32_threads4_queue8, CounterBits::Bits32, 4, 8);
pipeline_case!(bits8_threads3_queue3, CounterBits::Bits8, 3, 3);
pipeline_case!(bits16_threads3_queue3, CounterBits::Bits16, 3, 3);
pipeline_case!(bits32_threads3_queue3, CounterBits::Bits32, 3, 3);

#[test]
fn u8_saturates() {
    assert_saturation::<u8>(CounterBits::Bits8, 300, u8::MAX as u64);
}

#[test]
fn u16_exceeds_u8_range() {
    assert_saturation::<u16>(CounterBits::Bits16, 300, u8::MAX as u64 + 1);
}

#[test]
fn u32_exceeds_u8_range() {
    assert_saturation::<u32>(CounterBits::Bits32, 300, u8::MAX as u64 + 1);
}

#[test]
fn invalid_zero_k_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.k = 0;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_even_k_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.k = 10;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_large_k_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.k = 33;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_zero_m_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.m = 0;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_m_larger_than_k_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.m = cfg.k + 1;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_zero_s_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.s = 0;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_s_larger_than_k_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.s = cfg.k + 1;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_zero_hash_count_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.n_hashes = 0;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_block_size_equal_to_table_size_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.block_slots_exponent = cfg.counter_slots_exponent;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_zero_threads_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.threads = 0;
    assert!(cfg.validate().is_err());
}

#[test]
fn invalid_zero_queue_is_rejected() {
    let mut cfg = cfg(CounterBits::Bits8, 1, 1);
    cfg.queue = 0;
    assert!(cfg.validate().is_err());
}
