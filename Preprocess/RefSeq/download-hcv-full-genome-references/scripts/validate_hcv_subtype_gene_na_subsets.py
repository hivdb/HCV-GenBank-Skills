#!/usr/bin/env python3
"""Verify subtype per-gene NA references are continuous subsets of full genomes."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
GENES = ("NS3", "NS5A_NTD", "NS5B")
DEFAULT_REFERENCE_DIR = REPO_ROOT / "HCVData" / "Subtype-Ref"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence).upper()))
            header, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line).upper())
    if header:
        records.append((header, "".join(sequence).upper()))
    return records


def header_field(header: str, field: str) -> str:
    prefix = f"{field}="
    return next((item[len(prefix):] for item in header.split("|") if item.startswith(prefix)), "")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument(
        "--report-csv", type=Path,
        default=DEFAULT_REFERENCE_DIR / "HCV_Subtype_Gene_NA_Subset_Check.csv",
    )
    args = parser.parse_args()

    full_genomes = {
        header_field(header, "accession"): sequence
        for header, sequence in read_fasta(args.reference_dir / "HCV_Subtype_FullGenome_Refs.fasta")
    }
    fields = ["Gene", "Subtype", "Accession", "GeneNALength", "FullGenomeStartNA", "FullGenomeEndNA", "Status"]
    report_rows: list[dict[str, str]] = []
    failed = 0
    for gene in GENES:
        gene_path = args.reference_dir / f"HCV_Subtype_Refs_{gene}_NA.fasta"
        for header, sequence in read_fasta(gene_path):
            accession = header_field(header, "accession")
            subtype = header_field(header, "subtype")
            genome = full_genomes.get(accession, "")
            start = genome.find(sequence) if genome else -1
            status = "PASS" if start >= 0 else "FAIL"
            if status == "FAIL":
                failed += 1
            report_rows.append({
                "Gene": gene,
                "Subtype": subtype,
                "Accession": accession,
                "GeneNALength": str(len(sequence)),
                "FullGenomeStartNA": str(start + 1) if start >= 0 else "",
                "FullGenomeEndNA": str(start + len(sequence)) if start >= 0 else "",
                "Status": status,
            })
    args.report_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.report_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"report_csv={args.report_csv.resolve().relative_to(REPO_ROOT)}")
    print(f"checked_rows={len(report_rows)}")
    print(f"failed_rows={failed}")
    if failed:
        raise SystemExit("Per-gene NA subset validation failed")


if __name__ == "__main__":
    main()
