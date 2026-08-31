#!/usr/bin/env python3
"""Filter full-genome COMET calls to each gene's RAS-overlap accessions."""

from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STEP_OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras"
COMET_CSV = REPO_ROOT / "HCVData" / "Comet-Full-genome" / "all_comet_subtype.csv"
OUTPUT_DIR = STEP_OUTPUT_DIR / "06-filter-comet-full-genome-by-ras-accessions"
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
    return value.strip().split("|", 1)[0].split(".", 1)[0].upper()


def accession_set(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "Accession" not in reader.fieldnames:
            raise ValueError(f"{path} must contain an Accession column")
        return {accession_key(row.get("Accession", "")) for row in reader if accession_key(row.get("Accession", ""))}


def main() -> None:
    with COMET_CSV.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "name" not in reader.fieldnames:
            raise ValueError(f"{COMET_CSV} must contain a name column")
        comet_rows = list(reader)
        fieldnames = reader.fieldnames

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for gene, ras_accession_csv in RAS_ACCESSION_FILES.items():
        selected_accessions = accession_set(ras_accession_csv)
        rows = [row for row in comet_rows if accession_key(row.get("name", "")) in selected_accessions]
        output_csv = OUTPUT_DIR / f"{gene}_Comet_FullGenome_RAS_Overlap.csv"
        with output_csv.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.DictWriter(destination, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        kept_count = len({accession_key(row["name"]) for row in rows})
        print(f"{gene} kept accessions: {kept_count:,}")
        print(f"{gene} output: {output_csv}")


if __name__ == "__main__":
    main()
