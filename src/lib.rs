//! SuperBloom-style counting Bloom estimator for DNA k-mer abundance spectra.
//!
//! The high-level entry point is [`estimate_spectrum_from_fastx`]. Lower-level
//! users can build a typed [`SuperCountingBloom`] directly, freeze it, and query
//! packed sequences through [`FrozenSuperCountingBloom`].

mod api;
mod config;
mod counter;
mod error;
mod fastx;
mod filter;
mod spectrum;

pub use api::{
    AddReport, BuildIndexReport, FilterOnlySpectrumEstimate, QueryFastaReport,
    SequenceAbundanceSummary, SuperCountingBloomBuilder, SuperCountingBloomIndex,
    for_each_fastx_record,
};
pub use config::{CounterBits, SuperCountingBloomConfig};
pub use counter::Counter;
pub use error::SuperCountingBloomError;
pub use fastx::{
    StreamEstimateReport, StreamEstimateStats, SuperCountingBloomReport, SuperCountingBloomStats,
    estimate_spectrum_from_fastx, estimate_spectrum_from_fastx_pair, stream_estimates_from_fastx,
    stream_estimates_from_fastx_pair,
};
pub use filter::{FrozenSuperCountingBloom, SuperCountingBloom, canonical_kmer, xorshift_u64};
pub use spectrum::{ApproxSpectrum, write_spectrum};

pub type Result<T> = std::result::Result<T, SuperCountingBloomError>;
