# Benchmark Scripts

This folder provides Python scripts to benchmark the `super-counting-bloom`
binary by sweeping one parameter at a time while keeping the other parameters
fixed.

Each sweep script:

- takes a FASTA/FASTQ input;
- runs multiple benchmark points;
- writes a TSV file;
- writes a PNG plot when `matplotlib` is available.

## Prerequisites

- Python 3
- `matplotlib` for plots
- Rust toolchain

## Available scripts

- `sweep_counter_bits.py`
- `sweep_memory.py`
- `run_all_sweeps.py`

## Example: one sweep

```bash
python benchmark/sweep_memory.py \
  --input reads.fastq \
  --binary target/release/super-counting-bloom
```

Outputs:

- `benchmark/results/sweep_memory.tsv`
- `benchmark/results/sweep_memory.png`

## Example: run all sweeps

```bash
python benchmark/run_all_sweeps.py \
  --input reads.fastq \
  --binary target/release/super-counting-bloom
```
