#!/usr/bin/env python3
"""Rename gene CSV files in lexicographic filename order to numbered parts."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument("--prefix", required=True)
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not directory.is_dir():
        raise SystemExit(f"Directory not found: {directory}")
    sources = sorted(
        (path for path in directory.iterdir() if path.is_file() and path.suffix == ".csv"),
        key=lambda path: path.name,
    )
    if not sources:
        raise SystemExit(f"No CSV files found in {directory}")
    targets = [
        directory / f"{args.prefix}_part_{number:04d}.csv"
        for number in range(1, len(sources) + 1)
    ]
    if len(set(targets)) != len(targets):
        raise SystemExit("Target filenames are not unique")
    temporary = [directory / f".{source.name}.renaming" for source in sources]
    for source, temporary_path in zip(sources, temporary, strict=True):
        source.rename(temporary_path)
    for temporary_path, target in zip(temporary, targets, strict=True):
        temporary_path.rename(target)
        print(target)


if __name__ == "__main__":
    main()
