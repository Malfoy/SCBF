# SuperCountingBloom

![tests](https://github.com/EtienneC-K/SuperCountingBloom/workflows/tests/badge.svg)

A Rust implementation of a SuperCounting Bloom filter for streaming DNA
`k`-mer abundance indexing, abundance querying, and approximate `k`-mer spectra.


## What This Crate Provides

- A **library API** to:
  - build a SuperCountingBloom index from explicit parameters,
  - insert DNA from memory (`add_sequence`) or from FASTA/FASTQ (`add_fasta`),
  - query abundance from memory (`estimate_sequence_abundances`) or FASTA/FASTQ (`query_fasta`),
  - compute approximate `k`-mer abundance spectra,
  - serialize/deserialize frozen indexes (`save` / `load`),
  - choose 8-, 16-, or 32-bit saturating counters.
- CLI binaries for spectrum estimation, saved-index build/query workflows, and compact API examples.

## Build, Test, Run

### Prerequisites

- Rust toolchain
- Standard build tools for crates in `Cargo.toml`

### Build

```bash
cargo build -r
```

### Run tests

```bash
cargo test
```

## Fully Commented Example (Library Usage)

```rust
use super_counting_bloom::{
    CounterBits, SuperCountingBloomBuilder, SuperCountingBloomConfig, SuperCountingBloomIndex,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // 1) Start from explicit parameters.
    // Defaults are: k=31, m=21, s=25, h=4, 16-bit counters,
    // 2^30 counter slots, 2^9 counter slots per block.
    let config = SuperCountingBloomConfig {
        k: 31,
        m: 21,
        s: 25,
        n_hashes: 4,
        counter_bits: CounterBits::Bits16,
        counter_slots_exponent: 30,
        block_slots_exponent: 9,
        threads: 16,
        queue: 4096,
    };

    // 2) Build a mutable index.
    let mut builder = SuperCountingBloomBuilder::new(config)?;

    // 3) Configure threads used by parallel FASTA/FASTQ indexing.
    builder.set_threads(16)?;

    // 4) Insert one in-memory DNA sequence.
    let added = builder.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGTA")?;
    println!("added {added} k-mers from memory");

    // 5) Insert FASTA/FASTQ from disk.
    let add_report = builder.add_fasta("reads.fastq")?;
    println!("indexed {} records", add_report.records_indexed);

    // 6) Freeze into a query index.
    let index = builder.into_index()?;

    // 7) Query from memory.
    let abundances = index.estimate_sequence_abundances(
        b"ACGTACGTACGTACGTACGTACGTACGTACGTA",
    );
    let valid = abundances.iter().flatten().count();
    println!("memory query valid windows: {valid}/{}", abundances.len());

    // 8) Query an entire FASTA/FASTQ.
    let q = index.query_fasta("contigs.fa")?;
    println!("fasta valid k-mers: {}", q.valid_kmers);

    // 9) Save and load.
    index.save("/tmp/demo.scb")?;
    let loaded = SuperCountingBloomIndex::load("/tmp/demo.scb")?;
    println!("loaded inserted k-mers: {}", loaded.inserted_kmers());

    Ok(())
}
```

## Parameter Suggestions

All geometry is explicit in `SuperCountingBloomConfig`, and threading is part of the config or controlled at runtime on `SuperCountingBloomBuilder`.

### Defaults

`SuperCountingBloomConfig::default()` uses:

- `k = 31`
- `m = 21`
- `s = 25` (`k - 6`)
- `n_hashes = 4`
- `counter_bits = CounterBits::Bits16`
- `counter_slots_exponent = 30`
- `block_slots_exponent = 9`
- `threads = num_cpus::get()`
- `queue = 4096`

### Meaning and Influence

- `k` (k-mer length, default `31`)
  - Higher `k`: usually more specific abundance estimates, but more sensitivity to sequencing errors/variants.
  - Lower `k`: more tolerant matching, but higher chance of ambiguous signal.
  - Current compact key path supports `1..=32`; canonical SIMD minimizer tie-breaking requires odd `k`.

- `m` (minimizer length, default `21`)
  - Higher `m`: generally more minimizer changes and less super-k-mer grouping.
  - Lower `m`: larger super-k-mer groups and stronger locality, but potentially more block pressure.

- `s` (fimpera/findere-like subword length, default `25`)
  - Rule of thumb: `s = k - 6` is a strong baseline.
  - Higher `s` (closer to `k`): behavior gets closer to direct k-mer evidence.
  - Lower `s`: more overlapping subword evidence per k-mer, but more counter probes.

- `n_hashes` (default `4`)
  - Higher values reduce counter-collision overestimation up to a point.
  - Too high harms speed because every `s`-mer touches more counters.

- `counter_bits` (default `16`)
  - `8`: maximum count 255; 
  - `16`: maximum count 65,535; 
  - `32`: maximum count 4,294,967,295;

- `counter_slots_exponent` (total counters = `2^counter_slots_exponent`, default `30`)
  - Main memory/accuracy knob.
  - Higher values use more RAM and usually reduce overestimation.
  - Approximate counter memory is `2^counter_slots_exponent * counter_bits / 8` bytes.
  - With the default `2^30` counters, memory is about 1 GiB with 8-bit counters, 2 GiB with 16-bit counters, and 4 GiB with 32-bit counters.

- `block_slots_exponent` (counters per local block = `2^block_slots_exponent`, default `9`)
  - Controls locality granularity.
  - Smaller blocks improve locality but can increase block pressure.
  - Larger blocks reduce per-block pressure but can reduce cache efficiency.

- `threads`
  - Controls parallelism for FASTA/FASTQ insertion and aggregate query operations.
  - `SuperCountingBloomBuilder::set_threads(n)` changes the runtime value before `add_fasta`.

- `queue`
  - Bounded parser-to-worker queue length.
  - Larger values can smooth parsing/worker imbalance at the cost of buffered memory.

## What Happens Under the Hood

1. **Blocked layout for locality**  
   The filter is divided into local counter blocks so accesses for related queries stay localized.

2. **Minimizer-based super-k-mer grouping**  
   Consecutive k-mers often share a minimizer. SuperCountingBloom maps grouped k-mers to the same counter block, amortizing random accesses over sequence streaming.

3. **fimpera/findere-style abundance evidence through `s`-mers**  
   A queried `k`-mer is decomposed into overlapping `s`-mers. The `k`-mer abundance estimate is the minimum over the corresponding `s`-mer counter estimates.

4. **Counter-width specialization**  
   The crate monomorphizes 8-, 16-, and 32-bit saturating counter paths.

## Serialization

- Binary format with magic header: `SCBIDX`.
- `save` serializes config + inserted k-mer count + deterministic layout metadata + frozen counter shards.
- `load` validates the header, version, counter width, and layout before reconstructing the query-ready index.

## Performance Notes

- Querying full sequences benefits most from streaming locality and super-k-mer reuse.
- Isolated random k-mer queries behave closer to blocked counting Bloom behavior.
- Larger `counter_slots_exponent` values lower overestimation but increase memory.
- More hash probes can reduce overestimation but increase memory traffic.
- 8-bit counters are fast and compact, but saturation should be checked on high-depth data.

## Citation

Citation information for SuperCountingBloom will be added when the manuscript is available.
