#!/usr/bin/env python3
"""Retain NS3 coverage rows whose AA overlap contains an NS3 RAS position."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


NS3_RAS_POSITIONS = (36, 41, 43, 54, 55, 56, 80, 122, 155, 156, 158, 166, 168, 170, 175)
OVERLAP_RE = re.compile(r"^\s*(\d+)\s*-\s*(\d+)\s*$")
REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT = REPO_ROOT / "HCVData" / "nonComet-PerGene" / "NS3_AllSeq_NonComet_Coverage.csv"
DEFAULT_PROFILE_ACCESSIONS = (
    REPO_ROOT
    / "outputs"
    / "comet-NS3-all-ras"
    / "16_build-complete-profiles"
    / "NS3_Profile_Accessions_QC_Pass.csv"
)
DEFAULT_OUTPUT = (
    REPO_ROOT
    / "HCVData"
    / "subtyping-comparison-all-ras"
    / "01-filter-ns3-ras-overlap"
    / "NS3_AllSeq_NonComet_Coverage_RAS_Overlap.csv"
)


def ras_positions_in_overlap(value: str) -> list[int]:
    """Return NS3 RAS positions inside an inclusive ReferenceOverlapAA range."""
    match = OVERLAP_RE.fullmatch(value or "")
    if not match:
        return []
    start, end = sorted((int(match.group(1)), int(match.group(2))))
    return [position for position in NS3_RAS_POSITIONS if start <= position <= end]


def accession_key(value: str) -> str:
    return value.strip().split(".", 1)[0].upper()


def profile_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an accession column")
        return {
            accession_key(row.get("accession", ""))
            for row in reader
            if accession_key(row.get("accession", ""))
        }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--profile-accessions-csv", type=Path, default=DEFAULT_PROFILE_ACCESSIONS)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    selected_accessions = profile_accessions(args.profile_accessions_csv)
    print(f"QC-passed profile accessions: {len(selected_accessions):,}")
    with args.input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {args.input_csv}")
        required = {"Accession", "ReferenceOverlapAA"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{args.input_csv} is missing columns: {', '.join(sorted(missing))}")
        fieldnames = [*reader.fieldnames, "CoveredNS3RASPositions"]
        kept_rows: list[dict[str, str]] = []
        for row in reader:
            positions = ras_positions_in_overlap(row.get("ReferenceOverlapAA", ""))
            if (
                positions
                and accession_key(row.get("Accession", "")) in selected_accessions
            ):
                row["CoveredNS3RASPositions"] = ";".join(map(str, positions))
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
