#!/usr/bin/env python3
"""Retain coverage rows whose AA overlap contains one or more RAS positions."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


OVERLAP_RE = re.compile(r"^\s*(\d+)\s*-\s*(\d+)\s*$")


def parse_positions(value: str) -> tuple[int, ...]:
    try:
        positions = tuple(sorted({int(item.strip()) for item in value.split(",") if item.strip()}))
    except ValueError as error:
        raise argparse.ArgumentTypeError("--ras-positions must be comma-separated integers") from error
    if not positions:
        raise argparse.ArgumentTypeError("--ras-positions must not be empty")
    return positions


def positions_in_overlap(value: str, ras_positions: tuple[int, ...]) -> list[int]:
    """Return RAS positions inside an inclusive ReferenceOverlapAA range."""
    match = OVERLAP_RE.fullmatch(value or "")
    if not match:
        return []
    start, end = sorted((int(match.group(1)), int(match.group(2))))
    return [position for position in ras_positions if start <= position <= end]


def accession_key(value: str) -> str:
    """Normalize an accession, ignoring a version suffix."""
    return value.strip().split(".", 1)[0].upper()


def staged_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an Accession column")
        return {
            accession_key(row.get("Accession", ""))
            for row in reader
            if accession_key(row.get("Accession", ""))
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, required=True)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--ras-positions", type=parse_positions, required=True)
    parser.add_argument(
        "--staged-accessions-csv",
        type=Path,
        help="Optional CSV with an Accession column; retain only its accessions.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected_accessions = (
        staged_accessions(args.staged_accessions_csv) if args.staged_accessions_csv else None
    )
    if selected_accessions is not None:
        print(f"Staged accessions: {len(selected_accessions):,}")
    with args.input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {args.input_csv}")
        required = {"Accession", "ReferenceOverlapAA"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{args.input_csv} is missing columns: {', '.join(sorted(missing))}")
        fieldnames = [*reader.fieldnames, "CoveredRASPositions"]
        kept_rows: list[dict[str, str]] = []
        for row in reader:
            positions = positions_in_overlap(row.get("ReferenceOverlapAA", ""), args.ras_positions)
            if (
                positions
                and (
                    selected_accessions is None
                    or accession_key(row.get("Accession", "")) in selected_accessions
                )
            ):
                row["CoveredRASPositions"] = ";".join(map(str, positions))
                kept_rows.append(row)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept_rows)

    accessions = {accession_key(row["Accession"]) for row in kept_rows if row["Accession"].strip()}
    print(f"Kept accessions: {len(accessions):,}")
    print(f"Kept rows: {len(kept_rows):,}")
    print(f"Output: {args.output_csv}")


if __name__ == "__main__":
    main()
