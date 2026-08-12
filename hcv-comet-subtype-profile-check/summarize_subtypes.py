#!/usr/bin/env python3
"""Summarize Comet subtype assignments and profile inclusion by gene.

For NS3, NS5A, and NS5B, this script counts subtypes in the corresponding
Comet assignment CSV and in the accessions selected for profile building. It
writes one CSV per gene to this script's folder by default.
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


GENES = ("NS3", "NS5A", "NS5B")
COMMON_SUBTYPES = {"1a", "1b", "2a", "2b", "2c", "3a", "4a", "5a", "6a"}
REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignment-dir", type=Path, default=REPO_ROOT / "Comet Subtyping")
    parser.add_argument("--profile-dir", type=Path, default=REPO_ROOT / "outputs/comet")
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument("--gene", choices=GENES, action="append", help="Gene to summarize; repeat as needed.")
    return parser.parse_args()


def normalized_subtype(value: str | None) -> str:
    return (value or "").strip() or "Unassigned"


def subtype_counts(path: Path, accession_column: str) -> Counter[str]:
    with path.open(newline="", encoding="utf-8-sig") as input_file:
        reader = csv.DictReader(input_file)
        required = {accession_column, "subtype"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain columns: {', '.join(sorted(required))}")
        return Counter(normalized_subtype(row.get("subtype")) for row in reader)


def summarize_gene(gene: str, assignment_dir: Path, profile_dir: Path, output_dir: Path) -> Path:
    assignment_path = assignment_dir / f"{gene}.csv"
    profile_path = profile_dir / f"{gene}_Profile_Accessions.csv"
    if not assignment_path.exists():
        raise FileNotFoundError(f"Comet assignment file not found: {assignment_path}")
    if not profile_path.exists():
        raise FileNotFoundError(f"Profile-accession file not found: {profile_path}")

    assignment_counts = subtype_counts(assignment_path, "name")
    profile_counts = subtype_counts(profile_path, "accession")
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / f"{gene}_subtype_profile_summary.csv"
    with output_path.open("w", newline="", encoding="utf-8") as output_file:
        writer = csv.DictWriter(
            output_file,
            fieldnames=["Subtype", "NumAccs", "common subtype", "Number in profile"],
        )
        writer.writeheader()
        for subtype in sorted(assignment_counts.keys() | profile_counts.keys()):
            writer.writerow(
                {
                    "Subtype": subtype,
                    "NumAccs": assignment_counts[subtype],
                    "common subtype": "Yes" if subtype in COMMON_SUBTYPES else "No",
                    "Number in profile": profile_counts[subtype],
                }
            )
    return output_path


def main() -> None:
    args = parse_args()
    for gene in args.gene or GENES:
        output_path = summarize_gene(gene, args.assignment_dir, args.profile_dir, args.output_dir)
        print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
