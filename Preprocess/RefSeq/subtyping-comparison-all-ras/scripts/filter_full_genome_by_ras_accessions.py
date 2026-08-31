#!/usr/bin/env python3
"""Filter full-genome coverage files to per-gene RAS-overlap accessions."""

from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STEP_OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras"
FULL_GENOME_DIR = REPO_ROOT / "HCVData" / "nonComet-Full-genome"
OUTPUT_DIR = STEP_OUTPUT_DIR / "04-filter-full-genome-by-ras-accessions"
RAS_ACCESSION_FILES = {
    "NS3": STEP_OUTPUT_DIR
    / "01-filter-ns3-ras-overlap"
    / "NS3_AllSeq_NonComet_Coverage_RAS_Overlap.csv",
    "NS5A": STEP_OUTPUT_DIR
    / "02-filter-ns5a-ras-overlap"
    / "NS5A_AllSeq_NonComet_Coverage_RAS_Overlap.csv",
    "NS5B": STEP_OUTPUT_DIR
    / "03-filter-ns5b-ras-overlap"
    / "NS5B_AllSeq_NonComet_Coverage_RAS_Overlap.csv",
}


def accession_key(value: str) -> str:
    return value.strip().split(".", 1)[0].upper()


def accession_set(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an Accession column")
        return {accession_key(row.get("Accession", "")) for row in reader if accession_key(row.get("Accession", ""))}


def filter_coverage(gene: str, selected_accessions: set[str]) -> tuple[Path, int]:
    input_csv = FULL_GENOME_DIR / f"{gene}_AllSeq_NonComet_Coverage.csv"
    output_csv = OUTPUT_DIR / f"{gene}_AllSeq_NonComet_Coverage_RAS_Overlap.csv"
    with input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise ValueError(f"{input_csv} must contain an Accession column")
        rows = [row for row in reader if accession_key(row.get("Accession", "")) in selected_accessions]
        fieldnames = reader.fieldnames
    with output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return output_csv, len({accession_key(row["Accession"]) for row in rows})


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for gene, ras_accession_csv in RAS_ACCESSION_FILES.items():
        selected_accessions = accession_set(ras_accession_csv)
        output_csv, kept_count = filter_coverage(gene, selected_accessions)
        print(f"{gene} kept accessions: {kept_count:,}")
        print(f"{gene} output: {output_csv}")


if __name__ == "__main__":
    main()
