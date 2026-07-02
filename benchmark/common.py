#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import time
from pathlib import Path


METRIC_RE = re.compile(r"([A-Za-z0-9_]+):\s+([^\t]+)")


def base_parser(description: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("--input", required=True, help="FASTA/FASTQ input path")
    parser.add_argument("--binary", default="target/release/super-counting-bloom")
    parser.add_argument("--out-dir", default="benchmark/results")
    parser.add_argument("--k", type=int, default=31)
    parser.add_argument("--m", type=int, default=21)
    parser.add_argument("--s", type=int, default=27)
    parser.add_argument("--n-hashes", type=int, default=4)
    parser.add_argument("--counter-bits", type=int, default=16)
    parser.add_argument("--counter-slots-exponent", type=int, default=30)
    parser.add_argument("--block-slots-exponent", type=int, default=9)
    parser.add_argument("--threads", type=int, default=os.cpu_count() or 1)
    return parser


def run_point(args: argparse.Namespace, **overrides: object) -> dict[str, object]:
    params = {
        "k": args.k,
        "m": args.m,
        "s": args.s,
        "n_hashes": args.n_hashes,
        "counter_bits": args.counter_bits,
        "counter_slots_exponent": args.counter_slots_exponent,
        "block_slots_exponent": args.block_slots_exponent,
        "threads": args.threads,
    }
    params.update(overrides)

    command = [
        args.binary,
        args.input,
        "--k",
        str(params["k"]),
        "--m",
        str(params["m"]),
        "--s",
        str(params["s"]),
        "--n-hashes",
        str(params["n_hashes"]),
        "--counter-bits",
        str(params["counter_bits"]),
        "--counter-slots-exponent",
        str(params["counter_slots_exponent"]),
        "--block-slots-exponent",
        str(params["block_slots_exponent"]),
        "--threads",
        str(params["threads"]),
        "--stream-query-only",
    ]

    start = time.perf_counter()
    completed = subprocess.run(command, check=False, text=True, capture_output=True)
    elapsed = time.perf_counter() - start

    row: dict[str, object] = dict(params)
    row["elapsed_seconds"] = elapsed
    row["returncode"] = completed.returncode
    row["stderr"] = completed.stderr.strip()

    for match in METRIC_RE.finditer(completed.stderr):
        key = match.group(1)
        value = match.group(2)
        try:
            row[key] = float(value) if "." in value else int(value)
        except ValueError:
            row[key] = value

    return row


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_lines(path: Path, rows: list[dict[str, object]], x: str, y: str, label: str) -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        return

    groups: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        groups.setdefault(str(row.get(label, "all")), []).append(row)

    fig, ax = plt.subplots(figsize=(7, 4))
    for name, group in groups.items():
        group = sorted(group, key=lambda row: row[x])
        ax.plot([row[x] for row in group], [row[y] for row in group], marker="o", label=name)
    ax.set_xlabel(x.replace("_", " "))
    ax.set_ylabel(y.replace("_", " "))
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
