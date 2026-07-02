use std::collections::HashMap;

use packed_seq::SeqVec;
use super_counting_bloom::{
    Counter, CounterBits, FrozenSuperCountingBloom, SuperCountingBloom, SuperCountingBloomConfig,
};

fn cfg(k: usize, m: usize, s: usize, bits: CounterBits) -> SuperCountingBloomConfig {
    SuperCountingBloomConfig {
        k,
        m,
        s,
        n_hashes: 3,
        counter_bits: bits,
        counter_slots_exponent: 16,
        block_slots_exponent: 6,
        threads: 2,
        queue: 4,
    }
}

fn packed(sequence: &[u8]) -> packed_seq::PackedSeqVec {
    packed_seq::PackedSeqVec::from_ascii(sequence)
}

fn freeze_bloom<C: Counter>(
    cfg: &SuperCountingBloomConfig,
    sequence: &packed_seq::PackedSeqVec,
    repeats: usize,
) -> FrozenSuperCountingBloom<C> {
    let filter = SuperCountingBloom::<C>::new(cfg).unwrap();
    for _ in 0..repeats {
        let inserted = filter.insert_packed_sequence(sequence.as_slice(), cfg);
        assert!(inserted > 0);
    }
    filter.freeze()
}

fn assert_positive<C: Counter>(cfg: &SuperCountingBloomConfig, sequence: &[u8]) {
    let sequence = packed(sequence);
    let bloom = freeze_bloom::<C>(cfg, &sequence, 1);
    let estimates = bloom.estimate_abundances_packed_sequence(sequence.as_slice(), cfg);
    assert_eq!(estimates.len(), sequence.len() + 1 - cfg.k);
    assert!(estimates.iter().all(|&estimate| estimate.to_u64() > 0));
}

fn assert_vector_length<C: Counter>(cfg: &SuperCountingBloomConfig, sequence: &[u8]) {
    let sequence = packed(sequence);
    let bloom = freeze_bloom::<C>(cfg, &sequence, 2);
    assert_eq!(
        bloom
            .estimate_abundances_packed_sequence(sequence.as_slice(), cfg)
            .len(),
        sequence.len() + 1 - cfg.k
    );
}

fn assert_stream_matches_vector<C: Counter>(cfg: &SuperCountingBloomConfig, sequence: &[u8]) {
    let sequence = packed(sequence);
    let bloom = freeze_bloom::<C>(cfg, &sequence, 3);
    let vector = bloom.estimate_abundances_packed_sequence(sequence.as_slice(), cfg);
    let checksum = vector
        .iter()
        .fold(0_u64, |acc, estimate| acc.wrapping_add(estimate.to_u64()));
    assert_eq!(
        bloom.stream_estimate_packed_sequence(sequence.as_slice(), cfg),
        (vector.len() as u64, checksum)
    );
}

fn assert_map_query<C: Counter>(cfg: &SuperCountingBloomConfig, sequence: &[u8]) {
    let sequence = packed(sequence);
    let bloom = freeze_bloom::<C>(cfg, &sequence, 2);
    let mut estimates = HashMap::new();
    let queried = bloom.estimate_packed_sequence(sequence.as_slice(), cfg, &mut estimates);
    assert_eq!(queried, (sequence.len() + 1 - cfg.k) as u64);
    assert!(!estimates.is_empty());
    assert!(estimates.values().all(|&estimate| estimate.to_u64() > 0));
}

fn assert_repeats_raise<C: Counter>(cfg: &SuperCountingBloomConfig, sequence: &[u8]) {
    let sequence = packed(sequence);
    let once = freeze_bloom::<C>(cfg, &sequence, 1);
    let more = freeze_bloom::<C>(cfg, &sequence, 4);
    let once = once.estimate_abundances_packed_sequence(sequence.as_slice(), cfg);
    let more = more.estimate_abundances_packed_sequence(sequence.as_slice(), cfg);
    assert_eq!(once.len(), more.len());
    assert!(
        more.iter()
            .zip(once.iter())
            .all(|(&more, &once)| more.to_u64() >= once.to_u64())
    );
}

macro_rules! backend_case {
    ($name:ident, $counter:ty, $bits:expr, $k:expr, $m:expr, $s:expr, $sequence:expr) => {
        mod $name {
            use super::*;

            #[test]
            fn positive_estimates() {
                assert_positive::<$counter>(&cfg($k, $m, $s, $bits), $sequence);
            }

            #[test]
            fn abundance_vector_length() {
                assert_vector_length::<$counter>(&cfg($k, $m, $s, $bits), $sequence);
            }

            #[test]
            fn stream_matches_materialized_vector() {
                assert_stream_matches_vector::<$counter>(&cfg($k, $m, $s, $bits), $sequence);
            }

            #[test]
            fn map_query_consistency() {
                assert_map_query::<$counter>(&cfg($k, $m, $s, $bits), $sequence);
            }

            #[test]
            fn repeated_inserts_raise_estimates() {
                assert_repeats_raise::<$counter>(&cfg($k, $m, $s, $bits), $sequence);
            }
        }
    };
}

backend_case!(
    u8_k9_s7_repeat,
    u8,
    CounterBits::Bits8,
    9,
    5,
    7,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u16_k9_s7_repeat,
    u16,
    CounterBits::Bits16,
    9,
    5,
    7,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u32_k9_s7_repeat,
    u32,
    CounterBits::Bits32,
    9,
    5,
    7,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u8_k9_s9_repeat,
    u8,
    CounterBits::Bits8,
    9,
    5,
    9,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u16_k9_s9_repeat,
    u16,
    CounterBits::Bits16,
    9,
    5,
    9,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u32_k9_s9_repeat,
    u32,
    CounterBits::Bits32,
    9,
    5,
    9,
    b"ACGTACGTACGTACGTACGTACGTACGT"
);
backend_case!(
    u8_k11_s8_mixed,
    u8,
    CounterBits::Bits8,
    11,
    5,
    8,
    b"ACGTTGCATGTCAGTACGATCGTACGTTAGC"
);
backend_case!(
    u16_k11_s8_mixed,
    u16,
    CounterBits::Bits16,
    11,
    5,
    8,
    b"ACGTTGCATGTCAGTACGATCGTACGTTAGC"
);
backend_case!(
    u32_k11_s8_mixed,
    u32,
    CounterBits::Bits32,
    11,
    5,
    8,
    b"ACGTTGCATGTCAGTACGATCGTACGTTAGC"
);
backend_case!(
    u8_k13_s9_balanced,
    u8,
    CounterBits::Bits8,
    13,
    7,
    9,
    b"TGCATGCAACGTACGTTGCATGCAACGTACGT"
);
backend_case!(
    u16_k13_s9_balanced,
    u16,
    CounterBits::Bits16,
    13,
    7,
    9,
    b"TGCATGCAACGTACGTTGCATGCAACGTACGT"
);
backend_case!(
    u32_k13_s9_balanced,
    u32,
    CounterBits::Bits32,
    13,
    7,
    9,
    b"TGCATGCAACGTACGTTGCATGCAACGTACGT"
);
