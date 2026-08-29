#!/usr/bin/env python3
"""Plot a histogram of a numeric distance column from a subtyping CSV."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--column", default="FirstChoiceDistance")
    parser.add_argument("--bin-width", type=float, default=1.0)
    parser.add_argument("--output-png", type=Path)
    parser.add_argument(
        "--skip-empty",
        action="store_true",
        help="Print a message and skip output when the selected column has no numeric values.",
    )
    return parser.parse_args()


def default_output(input_csv: Path, column: str) -> Path:
    stem = input_csv.stem.removesuffix("_Genotype_Distances")
    return input_csv.parent / "figures" / f"{stem}_{column}_Histogram.png"


def main() -> None:
    args = parse_args()
    if args.bin_width <= 0:
        raise ValueError("--bin-width must be positive")

    values: list[float] = []
    with args.input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or args.column not in reader.fieldnames:
            raise ValueError(f"{args.input_csv} does not contain a {args.column!r} column")
        for row in reader:
            value = row.get(args.column, "").strip()
            if value:
                values.append(float(value))
    if not values:
        if args.skip_empty:
            print(f"Skipped: no numeric values found in {args.column!r} for {args.input_csv}")
            return
        raise ValueError(f"No numeric values found in {args.column!r}")

    lower = math.floor(min(values) / args.bin_width) * args.bin_width
    upper = math.ceil(max(values) / args.bin_width) * args.bin_width
    bin_count = max(1, math.ceil((upper - lower) / args.bin_width))
    bins = [lower + index * args.bin_width for index in range(bin_count + 1)]
    if len(bins) == 1:
        bins.append(lower + args.bin_width)

    output_png = args.output_png or default_output(args.input_csv, args.column)
    output_png.parent.mkdir(parents=True, exist_ok=True)
    mean = statistics.mean(values)
    median = statistics.median(values)

    figure, axis = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    axis.hist(values, bins=bins, color="#277da1", edgecolor="white", linewidth=0.7)
    axis.axvline(mean, color="#d1495b", linewidth=2, label=f"Mean: {mean:.2f}")
    axis.axvline(median, color="#f4a261", linewidth=2, linestyle="--", label=f"Median: {median:.2f}")
    axis.set_title(f"{args.input_csv.stem}: {args.column} distribution")
    axis.set_xlabel(f"{args.column} (100 − BLAST % identity)")
    axis.set_ylabel("Number of accessions")
    axis.set_xlim(lower, upper)
    axis.grid(axis="y", alpha=0.25)
    axis.legend()
    axis.text(
        0.99,
        0.97,
        f"n = {len(values):,}\nrange = {min(values):.2f}–{max(values):.2f}\nbin width = {args.bin_width:g}",
        transform=axis.transAxes,
        ha="right",
        va="top",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "edgecolor": "#bbbbbb"},
    )
    figure.savefig(output_png, dpi=200)
    plt.close(figure)

    print(f"Values plotted: {len(values):,}")
    print(f"Mean: {mean:.3f}")
    print(f"Median: {median:.3f}")
    print(f"Output: {output_png}")


if __name__ == "__main__":
    main()
