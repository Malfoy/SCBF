use std::{
    collections::{HashMap, VecDeque},
    io::{self, Read, Write},
    ops::Range,
    sync::Mutex,
};

use helicase::dna_format::PackedDNA;
use packed_seq::Seq;

use crate::{
    Counter, Result, SuperCountingBloomConfig, SuperCountingBloomError, counter::SerializedCounter,
};

pub type KmerKey = u64;

struct CounterShard<C: Counter> {
    counters: Vec<C>,
    subblocks: usize,
    block_slots: usize,
}

impl<C: Counter> CounterShard<C> {
    fn new(subblocks: usize, block_slots: usize) -> Self {
        Self {
            counters: vec![C::default(); subblocks * block_slots],
            subblocks,
            block_slots,
        }
    }

    #[inline(always)]
    fn index(&self, subblock: usize, address: usize) -> usize {
        debug_assert!(subblock < self.subblocks);
        debug_assert!(address < self.block_slots);
        subblock * self.block_slots + address
    }

    #[inline(always)]
    fn increment(&mut self, subblock: usize, address: usize) {
        let idx = self.index(subblock, address);
        self.counters[idx].saturating_increment();
    }

    #[inline(always)]
    fn get(&self, subblock: usize, address: usize) -> C {
        self.counters[self.index(subblock, address)]
    }

    fn add_counter_histogram(&self, histogram: &mut HashMap<u64, u64>) {
        for &counter in &self.counters {
            *histogram.entry(counter.to_u64()).or_insert(0) += 1;
        }
    }
}

pub struct SuperCountingBloom<C: Counter> {
    shards: Vec<Mutex<CounterShard<C>>>,
    shard_mask: usize,
    shard_shift: u32,
    block_slots_mask: usize,
    nb_blocks: usize,
    block_mask: usize,
    n_hashes: usize,
}

pub struct FrozenSuperCountingBloom<C: Counter> {
    shards: Vec<CounterShard<C>>,
    shard_mask: usize,
    shard_shift: u32,
    block_slots_mask: usize,
    nb_blocks: usize,
    block_mask: usize,
    n_hashes: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct FrozenFilterLayout {
    pub shard_count: usize,
    pub subblocks: usize,
    pub block_slots: usize,
    pub shard_mask: usize,
    pub shard_shift: u32,
    pub block_slots_mask: usize,
    pub nb_blocks: usize,
    pub block_mask: usize,
    pub n_hashes: usize,
}

impl FrozenFilterLayout {
    pub(crate) fn from_config(cfg: &SuperCountingBloomConfig) -> Result<Self> {
        cfg.validate()?;
        let counter_slots = cfg.counter_slots()?;
        let block_slots = cfg.block_slots()?;
        let nb_blocks = counter_slots / block_slots;
        if !nb_blocks.is_power_of_two() {
            return Err(SuperCountingBloomError::InvalidConfig(
                "number of blocks must be a power of two".to_string(),
            ));
        }

        let shard_count = nb_blocks.min(1024);
        let subblocks = nb_blocks / shard_count;
        Ok(Self {
            shard_count,
            subblocks,
            block_slots,
            shard_mask: shard_count - 1,
            shard_shift: shard_count.trailing_zeros(),
            block_slots_mask: block_slots - 1,
            nb_blocks,
            block_mask: nb_blocks - 1,
            n_hashes: cfg.n_hashes,
        })
    }
}

impl<C: Counter> SuperCountingBloom<C> {
    pub fn new(cfg: &SuperCountingBloomConfig) -> Result<Self> {
        cfg.validate()?;
        let counter_slots = cfg.counter_slots()?;
        let block_slots = cfg.block_slots()?;
        let nb_blocks = counter_slots / block_slots;
        if !nb_blocks.is_power_of_two() {
            return Err(SuperCountingBloomError::InvalidConfig(
                "number of blocks must be a power of two".to_string(),
            ));
        }

        let shard_count = nb_blocks.min(1024);
        let subblocks = nb_blocks / shard_count;
        let shards = (0..shard_count)
            .map(|_| Mutex::new(CounterShard::new(subblocks, block_slots)))
            .collect();

        Ok(Self {
            shards,
            shard_mask: shard_count - 1,
            shard_shift: shard_count.trailing_zeros(),
            block_slots_mask: block_slots - 1,
            nb_blocks,
            block_mask: nb_blocks - 1,
            n_hashes: cfg.n_hashes,
        })
    }

    pub fn freeze(self) -> FrozenSuperCountingBloom<C> {
        let shards = self
            .shards
            .into_iter()
            .map(|shard| shard.into_inner().expect("counter shard mutex poisoned"))
            .collect();

        FrozenSuperCountingBloom {
            shards,
            shard_mask: self.shard_mask,
            shard_shift: self.shard_shift,
            block_slots_mask: self.block_slots_mask,
            nb_blocks: self.nb_blocks,
            block_mask: self.block_mask,
            n_hashes: self.n_hashes,
        }
    }

    pub fn counter_slots(&self) -> usize {
        self.shards
            .iter()
            .map(|shard| {
                shard
                    .lock()
                    .expect("counter shard mutex poisoned")
                    .counters
                    .len()
            })
            .sum()
    }

    pub fn blocks(&self) -> usize {
        self.nb_blocks
    }

    pub fn counter_histogram(&self) -> Vec<(u64, u64)> {
        let mut histogram = HashMap::<u64, u64>::new();
        for shard in &self.shards {
            shard
                .lock()
                .expect("counter shard mutex poisoned")
                .add_counter_histogram(&mut histogram);
        }
        let mut rows: Vec<_> = histogram.into_iter().collect();
        rows.sort_unstable_by_key(|&(count, _)| count);
        rows
    }

    pub fn insert_packed_dna(&self, sequence: PackedDNA, cfg: &SuperCountingBloomConfig) -> u64 {
        self.insert_packed_sequence(sequence.as_packed_seq(), cfg)
    }

    pub fn insert_packed_sequence(
        &self,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
    ) -> u64 {
        if sequence.len() < cfg.k {
            return 0;
        }

        let mut inserted_smers = 0_u64;
        let (starts, minimizers) = super_kmers(sequence, cfg.k, cfg.m);
        for i in 0..starts.len() {
            let start_kmer = starts[i] as usize;
            let end_kmer = if i + 1 < starts.len() {
                starts[i + 1] as usize
            } else {
                sequence.len() + 1 - cfg.k
            };
            let block = xorshift_u64(minimizers[i]) as usize & self.block_mask;
            let shard_idx = block & self.shard_mask;
            let subblock = block >> self.shard_shift;
            let mut shard = self.shards[shard_idx]
                .lock()
                .expect("counter shard mutex poisoned");

            let end_smer = end_kmer + (cfg.k - cfg.s);
            for pos in start_kmer..end_smer {
                let smer = canonical_kmer(sequence, cfg.s, pos);
                let mut hash = xorshift_u64(smer);
                for _ in 0..self.n_hashes {
                    shard.increment(subblock, hash as usize & self.block_slots_mask);
                    hash = xorshift_u64(hash);
                }
                inserted_smers += 1;
            }
        }

        inserted_smers
    }
}

impl<C: Counter> FrozenSuperCountingBloom<C> {
    pub(crate) fn serialized_layout(&self) -> FrozenFilterLayout {
        let first_shard = self
            .shards
            .first()
            .expect("frozen filter always has at least one shard");
        FrozenFilterLayout {
            shard_count: self.shards.len(),
            subblocks: first_shard.subblocks,
            block_slots: first_shard.block_slots,
            shard_mask: self.shard_mask,
            shard_shift: self.shard_shift,
            block_slots_mask: self.block_slots_mask,
            nb_blocks: self.nb_blocks,
            block_mask: self.block_mask,
            n_hashes: self.n_hashes,
        }
    }

    pub fn counter_slots(&self) -> usize {
        self.shards.iter().map(|shard| shard.counters.len()).sum()
    }

    pub fn blocks(&self) -> usize {
        self.nb_blocks
    }

    pub fn counter_histogram(&self) -> Vec<(u64, u64)> {
        let mut histogram = HashMap::<u64, u64>::new();
        for shard in &self.shards {
            shard.add_counter_histogram(&mut histogram);
        }
        let mut rows: Vec<_> = histogram.into_iter().collect();
        rows.sort_unstable_by_key(|&(count, _)| count);
        rows
    }

    pub fn estimate_packed_dna(
        &self,
        sequence: PackedDNA,
        cfg: &SuperCountingBloomConfig,
        estimates: &mut HashMap<KmerKey, C>,
    ) -> u64 {
        self.estimate_packed_sequence(sequence.as_packed_seq(), cfg, estimates)
    }

    pub fn estimate_abundances_packed_dna(
        &self,
        sequence: PackedDNA,
        cfg: &SuperCountingBloomConfig,
    ) -> Vec<C> {
        self.estimate_abundances_packed_sequence(sequence.as_packed_seq(), cfg)
    }

    pub fn stream_estimate_packed_dna(
        &self,
        sequence: PackedDNA,
        cfg: &SuperCountingBloomConfig,
    ) -> (u64, u64) {
        self.stream_estimate_packed_sequence(sequence.as_packed_seq(), cfg)
    }

    pub fn estimate_packed_sequence(
        &self,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
        estimates: &mut HashMap<KmerKey, C>,
    ) -> u64 {
        if sequence.len() < cfg.k {
            return 0;
        }

        let mut queried = 0_u64;
        let (starts, minimizers) = super_kmers(sequence, cfg.k, cfg.m);
        for i in 0..starts.len() {
            let start_kmer = starts[i] as usize;
            let end_kmer = if i + 1 < starts.len() {
                starts[i + 1] as usize
            } else {
                sequence.len() + 1 - cfg.k
            };
            let block = xorshift_u64(minimizers[i]) as usize & self.block_mask;
            let shard_idx = block & self.shard_mask;
            let subblock = block >> self.shard_shift;
            let shard = &self.shards[shard_idx];
            queried += self.estimate_super_kmer_in_shard(
                shard,
                subblock,
                sequence,
                cfg,
                start_kmer..end_kmer,
                estimates,
            );
        }
        queried
    }

    pub fn estimate_abundances_packed_sequence(
        &self,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
    ) -> Vec<C> {
        if sequence.len() < cfg.k {
            return Vec::new();
        }

        let mut abundances = Vec::with_capacity(sequence.len() + 1 - cfg.k);
        let (starts, minimizers) = super_kmers(sequence, cfg.k, cfg.m);
        for i in 0..starts.len() {
            let start_kmer = starts[i] as usize;
            let end_kmer = if i + 1 < starts.len() {
                starts[i + 1] as usize
            } else {
                sequence.len() + 1 - cfg.k
            };
            let block = xorshift_u64(minimizers[i]) as usize & self.block_mask;
            let shard_idx = block & self.shard_mask;
            let subblock = block >> self.shard_shift;
            let shard = &self.shards[shard_idx];
            self.estimate_abundances_super_kmer_in_shard(
                shard,
                subblock,
                sequence,
                cfg,
                start_kmer..end_kmer,
                &mut abundances,
            );
        }
        abundances
    }

    pub fn stream_estimate_packed_sequence(
        &self,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
    ) -> (u64, u64) {
        if sequence.len() < cfg.k {
            return (0, 0);
        }

        let mut queried = 0_u64;
        let mut checksum = 0_u64;
        let (starts, minimizers) = super_kmers(sequence, cfg.k, cfg.m);
        for i in 0..starts.len() {
            let start_kmer = starts[i] as usize;
            let end_kmer = if i + 1 < starts.len() {
                starts[i + 1] as usize
            } else {
                sequence.len() + 1 - cfg.k
            };
            let block = xorshift_u64(minimizers[i]) as usize & self.block_mask;
            let shard_idx = block & self.shard_mask;
            let subblock = block >> self.shard_shift;
            let shard = &self.shards[shard_idx];
            let (local_queried, local_checksum) = self.stream_estimate_super_kmer_in_shard(
                shard,
                subblock,
                sequence,
                cfg,
                start_kmer..end_kmer,
            );
            queried += local_queried;
            checksum = checksum.wrapping_add(local_checksum);
        }
        (queried, checksum)
    }

    fn estimate_super_kmer_in_shard(
        &self,
        shard: &CounterShard<C>,
        subblock: usize,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
        kmer_range: Range<usize>,
        estimates: &mut HashMap<KmerKey, C>,
    ) -> u64 {
        let start_kmer = kmer_range.start;
        let end_kmer = kmer_range.end;
        let window = cfg.k - cfg.s + 1;
        let kmer_count = end_kmer - start_kmer;
        let smer_count = kmer_count + window - 1;
        let mut smer_estimates = Vec::with_capacity(smer_count);

        for offset in 0..smer_count {
            let smer = canonical_kmer(sequence, cfg.s, start_kmer + offset);
            smer_estimates.push(self.estimate_smer_in_shard(shard, subblock, smer));
        }

        let mut deque = VecDeque::<usize>::new();
        let mut queried = 0_u64;

        for offset in 0..smer_estimates.len() {
            while let Some(&back) = deque.back() {
                if smer_estimates[back] >= smer_estimates[offset] {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(offset);

            while let Some(&front) = deque.front() {
                if front + window <= offset {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if offset + 1 >= window {
                let kmer_offset = offset + 1 - window;
                let kmer_pos = start_kmer + kmer_offset;
                let key = canonical_kmer(sequence, cfg.k, kmer_pos);
                let estimate = smer_estimates[*deque.front().expect("non-empty minimum deque")];
                estimates.entry(key).or_insert(estimate);
                queried += 1;
            }
        }

        queried
    }

    fn estimate_abundances_super_kmer_in_shard(
        &self,
        shard: &CounterShard<C>,
        subblock: usize,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
        kmer_range: Range<usize>,
        abundances: &mut Vec<C>,
    ) {
        let start_kmer = kmer_range.start;
        let end_kmer = kmer_range.end;
        let window = cfg.k - cfg.s + 1;
        let kmer_count = end_kmer - start_kmer;
        let smer_count = kmer_count + window - 1;
        let mut smer_estimates = Vec::with_capacity(smer_count);

        for offset in 0..smer_count {
            let smer = canonical_kmer(sequence, cfg.s, start_kmer + offset);
            smer_estimates.push(self.estimate_smer_in_shard(shard, subblock, smer));
        }

        let mut deque = VecDeque::<usize>::new();

        for offset in 0..smer_estimates.len() {
            while let Some(&back) = deque.back() {
                if smer_estimates[back] >= smer_estimates[offset] {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(offset);

            while let Some(&front) = deque.front() {
                if front + window <= offset {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if offset + 1 >= window {
                let estimate = smer_estimates[*deque.front().expect("non-empty minimum deque")];
                abundances.push(estimate);
            }
        }
    }

    fn stream_estimate_super_kmer_in_shard(
        &self,
        shard: &CounterShard<C>,
        subblock: usize,
        sequence: packed_seq::PackedSeq<'_>,
        cfg: &SuperCountingBloomConfig,
        kmer_range: Range<usize>,
    ) -> (u64, u64) {
        let start_kmer = kmer_range.start;
        let end_kmer = kmer_range.end;
        let window = cfg.k - cfg.s + 1;
        let kmer_count = end_kmer - start_kmer;
        let smer_count = kmer_count + window - 1;
        let mut smer_estimates = Vec::with_capacity(smer_count);

        for offset in 0..smer_count {
            let smer = canonical_kmer(sequence, cfg.s, start_kmer + offset);
            smer_estimates.push(self.estimate_smer_in_shard(shard, subblock, smer));
        }

        let mut deque = VecDeque::<usize>::new();
        let mut queried = 0_u64;
        let mut checksum = 0_u64;

        for offset in 0..smer_estimates.len() {
            while let Some(&back) = deque.back() {
                if smer_estimates[back] >= smer_estimates[offset] {
                    deque.pop_back();
                } else {
                    break;
                }
            }
            deque.push_back(offset);

            while let Some(&front) = deque.front() {
                if front + window <= offset {
                    deque.pop_front();
                } else {
                    break;
                }
            }

            if offset + 1 >= window {
                let estimate = smer_estimates[*deque.front().expect("non-empty minimum deque")];
                checksum = checksum.wrapping_add(estimate.to_u64());
                queried += 1;
            }
        }

        (queried, checksum)
    }

    #[inline(always)]
    fn estimate_smer_in_shard(&self, shard: &CounterShard<C>, subblock: usize, smer: u64) -> C {
        let mut hash = xorshift_u64(smer);
        let mut estimate = C::max_value();
        for _ in 0..self.n_hashes {
            estimate = estimate.min(shard.get(subblock, hash as usize & self.block_slots_mask));
            hash = xorshift_u64(hash);
        }
        estimate
    }
}

impl<C: SerializedCounter> FrozenSuperCountingBloom<C> {
    pub(crate) fn write_counter_data<W: Write>(&self, writer: &mut W) -> io::Result<()> {
        for shard in &self.shards {
            C::write_slice(writer, &shard.counters)?;
        }
        Ok(())
    }

    pub(crate) fn from_counter_data<R: Read>(
        reader: &mut R,
        layout: FrozenFilterLayout,
    ) -> io::Result<Self> {
        let shard_len = layout.subblocks * layout.block_slots;
        let mut shards = Vec::with_capacity(layout.shard_count);
        for _ in 0..layout.shard_count {
            shards.push(CounterShard {
                counters: C::read_vec(reader, shard_len)?,
                subblocks: layout.subblocks,
                block_slots: layout.block_slots,
            });
        }

        Ok(Self {
            shards,
            shard_mask: layout.shard_mask,
            shard_shift: layout.shard_shift,
            block_slots_mask: layout.block_slots_mask,
            nb_blocks: layout.nb_blocks,
            block_mask: layout.block_mask,
            n_hashes: layout.n_hashes,
        })
    }
}

pub(crate) fn super_kmers(
    sequence: packed_seq::PackedSeq<'_>,
    k: usize,
    m: usize,
) -> (Vec<u32>, Vec<u64>) {
    let mut super_starts = Vec::new();
    let mut minimizer_positions = Vec::new();
    let output = simd_minimizers::canonical_minimizers(m, k - m + 1)
        .super_kmers(&mut super_starts)
        .run(sequence, &mut minimizer_positions);
    let minimizers = output.values_u64().collect();
    (super_starts, minimizers)
}

#[inline(always)]
pub fn canonical_kmer(sequence: packed_seq::PackedSeq<'_>, len: usize, pos: usize) -> u64 {
    sequence
        .read_kmer(len, pos)
        .min(sequence.read_revcomp_kmer(len, pos))
}

#[inline(always)]
pub fn xorshift_u64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

#[cfg(test)]
mod tests {
    use super::*;
    use packed_seq::SeqVec;

    fn test_cfg() -> SuperCountingBloomConfig {
        SuperCountingBloomConfig {
            k: 9,
            m: 5,
            s: 7,
            n_hashes: 3,
            counter_bits: crate::CounterBits::Bits8,
            counter_slots_exponent: 16,
            block_slots_exponent: 6,
            threads: 1,
            queue: 4,
        }
    }

    fn estimates_after_repeated_inserts<C: Counter>(
        cfg: &SuperCountingBloomConfig,
        seq: &packed_seq::PackedSeqVec,
        repeats: usize,
    ) -> HashMap<KmerKey, C> {
        let filter = SuperCountingBloom::<C>::new(cfg).unwrap();
        for _ in 0..repeats {
            filter.insert_packed_sequence(seq.as_slice(), cfg);
        }
        let filter = filter.freeze();

        let mut estimates = HashMap::new();
        filter.estimate_packed_sequence(seq.as_slice(), cfg, &mut estimates);
        estimates
    }

    fn assert_inserted_sequence_has_positive_estimates<C: Counter>() {
        let cfg = test_cfg();
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGT");
        let estimates = estimates_after_repeated_inserts::<C>(&cfg, &seq, 1);
        assert!(!estimates.is_empty());
        assert!(estimates.values().all(|&count| count.to_u64() > 0));
    }

    #[test]
    fn inserted_sequence_has_positive_estimates_for_all_counter_widths() {
        assert_inserted_sequence_has_positive_estimates::<u8>();
        assert_inserted_sequence_has_positive_estimates::<u16>();
        assert_inserted_sequence_has_positive_estimates::<u32>();
    }

    #[test]
    fn abundance_vector_has_one_entry_per_kmer() {
        let cfg = test_cfg();
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGT");
        let filter = SuperCountingBloom::<u16>::new(&cfg).unwrap();
        filter.insert_packed_sequence(seq.as_slice(), &cfg);
        let filter = filter.freeze();

        let abundances = filter.estimate_abundances_packed_sequence(seq.as_slice(), &cfg);
        assert_eq!(abundances.len(), seq.len() + 1 - cfg.k);
        assert!(abundances.iter().all(|&count| count > 0));
    }

    #[test]
    fn repeated_sequence_increases_estimate_for_all_counter_widths() {
        let cfg = test_cfg();
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGT");
        for estimates in [
            estimates_after_repeated_inserts::<u8>(&cfg, &seq, 2)
                .values()
                .map(|&x| x.to_u64())
                .collect::<Vec<_>>(),
            estimates_after_repeated_inserts::<u16>(&cfg, &seq, 2)
                .values()
                .map(|&x| x.to_u64())
                .collect::<Vec<_>>(),
            estimates_after_repeated_inserts::<u32>(&cfg, &seq, 2)
                .values()
                .map(|&x| x.to_u64())
                .collect::<Vec<_>>(),
        ] {
            assert!(estimates.iter().all(|&count| count >= 2));
        }
    }

    #[test]
    fn counters_saturate_at_u8_max() {
        let cfg = test_cfg();
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGT");
        let estimates = estimates_after_repeated_inserts::<u8>(&cfg, &seq, 300);
        assert!(estimates.values().all(|&count| count == u8::MAX));
    }

    #[test]
    fn wider_counters_can_exceed_u8_max() {
        let cfg = test_cfg();
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTACGTACGTACGTACGT");
        let estimates16 = estimates_after_repeated_inserts::<u16>(&cfg, &seq, 300);
        let estimates32 = estimates_after_repeated_inserts::<u32>(&cfg, &seq, 300);
        assert!(
            estimates16
                .values()
                .all(|&count| count.to_u64() > u64::from(u8::MAX))
        );
        assert!(
            estimates32
                .values()
                .all(|&count| count.to_u64() > u64::from(u8::MAX))
        );
    }

    #[test]
    fn canonical_key_matches_reverse_complement() {
        let seq = packed_seq::PackedSeqVec::from_ascii(b"ACGTAC");
        let rc = packed_seq::PackedSeqVec::from_ascii(b"GTACGT");
        assert_eq!(
            canonical_kmer(seq.as_slice(), 6, 0),
            canonical_kmer(rc.as_slice(), 6, 0)
        );
    }
}
