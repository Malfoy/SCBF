#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--binary", default="target/release/super-counting-bloom")
    parser.add_argument("--out-dir", default="benchmark/results")
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    args = parser.parse_args()

    common = [
        "--input",
        args.input,
        "--binary",
        args.binary,
        "--out-dir",
        args.out_dir,
        "--threads",
        str(args.threads),
    ]
    for script in ["sweep_counter_bits.py", "sweep_memory.py"]:
        subprocess.run(["python", f"benchmark/{script}", *common], check=True)


if __name__ == "__main__":
    main()
