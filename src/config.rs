use crate::{Result, SuperCountingBloomError};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CounterBits {
    Bits8,
    Bits16,
    Bits32,
}

impl CounterBits {
    pub const fn as_u8(self) -> u8 {
        match self {
            Self::Bits8 => 8,
            Self::Bits16 => 16,
            Self::Bits32 => 32,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SuperCountingBloomConfig {
    pub k: usize,
    pub m: usize,
    pub s: usize,
    pub n_hashes: usize,
    pub counter_bits: CounterBits,
    pub counter_slots_exponent: u8,
    pub block_slots_exponent: u8,
    pub threads: usize,
    pub queue: usize,
}

impl Default for SuperCountingBloomConfig {
    fn default() -> Self {
        Self {
            k: 31,
            m: 21,
            s: 25,
            n_hashes: 4,
            counter_bits: CounterBits::Bits16,
            counter_slots_exponent: 30,
            block_slots_exponent: 9,
            threads: num_cpus::get().max(1),
            queue: 4096,
        }
    }
}

impl SuperCountingBloomConfig {
    pub fn validate(&self) -> Result<()> {
        if self.k == 0 || self.k > 32 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "k must be in 1..=32".to_string(),
            ));
        }
        if self.m == 0 || self.m > self.k {
            return Err(SuperCountingBloomError::InvalidConfig(
                "m must be in 1..=k".to_string(),
            ));
        }
        if self.s == 0 || self.s > self.k || self.s > 32 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "s must be in 1..=k and <= 32".to_string(),
            ));
        }
        if self.k.is_multiple_of(2) {
            return Err(SuperCountingBloomError::InvalidConfig(
                "canonical SIMD minimizers require odd k".to_string(),
            ));
        }
        if self.n_hashes == 0 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "n_hashes must be positive".to_string(),
            ));
        }
        if self.block_slots_exponent >= self.counter_slots_exponent {
            return Err(SuperCountingBloomError::InvalidConfig(
                "block_slots_exponent must be smaller than counter_slots_exponent".to_string(),
            ));
        }
        if self.threads == 0 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "threads must be positive".to_string(),
            ));
        }
        if self.queue == 0 {
            return Err(SuperCountingBloomError::InvalidConfig(
                "queue must be positive".to_string(),
            ));
        }
        Ok(())
    }

    pub(crate) fn counter_slots(&self) -> Result<usize> {
        checked_pow2(self.counter_slots_exponent, "counter_slots_exponent")
    }

    pub(crate) fn block_slots(&self) -> Result<usize> {
        checked_pow2(self.block_slots_exponent, "block_slots_exponent")
    }
}

pub(crate) fn checked_pow2(exponent: u8, name: &str) -> Result<usize> {
    1_usize
        .checked_shl(u32::from(exponent))
        .ok_or_else(|| SuperCountingBloomError::InvalidConfig(format!("{name} is too large")))
}
