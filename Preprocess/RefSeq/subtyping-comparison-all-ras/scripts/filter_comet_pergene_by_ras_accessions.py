#!/usr/bin/env python3
"""Filter COMET per-gene calls to per-gene RAS-overlap accessions."""

from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STEP_OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras"
COMET_DIR = REPO_ROOT / "HCVData" / "Comet-PerGene"
OUTPUT_DIR = STEP_OUTPUT_DIR / "05-filter-comet-pergene-by-ras-accessions"
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
    """Normalize a plain accession or the accession prefix of a COMET name."""
    return value.strip().split("|", 1)[0].split(".", 1)[0].upper()


def accession_set(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an Accession column")
        return {accession_key(row.get("Accession", "")) for row in reader if accession_key(row.get("Accession", ""))}


def filter_comet(gene: str, selected_accessions: set[str]) -> tuple[Path, int]:
    input_csv = COMET_DIR / f"{gene}.csv"
    output_csv = OUTPUT_DIR / f"{gene}_RAS_Overlap.csv"
    with input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "name" not in reader.fieldnames:
            raise ValueError(f"{input_csv} must contain a name column")
        rows = [row for row in reader if accession_key(row.get("name", "")) in selected_accessions]
        fieldnames = reader.fieldnames
    with output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return output_csv, len({accession_key(row["name"]) for row in rows})


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for gene, ras_accession_csv in RAS_ACCESSION_FILES.items():
        selected_accessions = accession_set(ras_accession_csv)
        output_csv, kept_count = filter_comet(gene, selected_accessions)
        print(f"{gene} kept accessions: {kept_count:,}")
        print(f"{gene} output: {output_csv}")


if __name__ == "__main__":
    main()
