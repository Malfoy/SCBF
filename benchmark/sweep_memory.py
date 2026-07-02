#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from common import base_parser, plot_lines, run_point, write_tsv


def main() -> None:
    parser = base_parser("Sweep counter table size")
    parser.add_argument("--values", default="26,27,28,29,30")
    args = parser.parse_args()
    rows = [
        run_point(args, counter_slots_exponent=int(value))
        for value in args.values.split(",")
        if value.strip()
    ]
    out_dir = Path(args.out_dir)
    write_tsv(out_dir / "sweep_memory.tsv", rows)
    plot_lines(out_dir / "sweep_memory.png", rows, "counter_slots_exponent", "total_s", "counter_bits")


if __name__ == "__main__":
    main()
