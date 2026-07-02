#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from common import base_parser, plot_lines, run_point, write_tsv


def main() -> None:
    parser = base_parser("Sweep counter widths")
    parser.add_argument("--values", default="8,16,32")
    args = parser.parse_args()
    rows = [
        run_point(args, counter_bits=int(value))
        for value in args.values.split(",")
        if value.strip()
    ]
    out_dir = Path(args.out_dir)
    write_tsv(out_dir / "sweep_counter_bits.tsv", rows)
    plot_lines(out_dir / "sweep_counter_bits.png", rows, "counter_bits", "total_s", "structure")


if __name__ == "__main__":
    main()
