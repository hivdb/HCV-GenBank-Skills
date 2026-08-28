#!/usr/bin/env python3
"""Split HCV genotype per-gene NA references and validate them against full genomes."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
GENES = ("NS3", "NS5A", "NS5B")
DEFAULT_REFERENCE_DIR = REPO_ROOT / "HCVData" / "Genotype-Ref"


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


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def full_header_field(header: str, field: str) -> str:
    prefix = f"{field}="
    return next((item[len(prefix):] for item in header.split("|") if item.startswith(prefix)), "")


def parse_gene_header(header: str) -> tuple[str, str, str]:
    match = re.fullmatch(r"HCV([1-8])(NS3|NS5A|NS5B)\s*\|\s*(\S+)", header)
    if not match:
        raise ValueError(f"Invalid genotype per-gene FASTA header: {header}")
    return match.groups()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument(
        "--per-gene-fasta", type=Path,
        default=DEFAULT_REFERENCE_DIR / "HCV_GT_PerGene_NA_RefSeqs.fasta",
    )
    parser.add_argument(
        "--full-genome-fasta", type=Path,
        default=DEFAULT_REFERENCE_DIR / "HCV_GT_FullGenome_RefSeqs.fasta",
    )
    parser.add_argument(
        "--report-csv", type=Path,
        default=DEFAULT_REFERENCE_DIR / "HCV_GT_Gene_NA_Subset_Check.csv",
    )
    args = parser.parse_args()

    full_genomes = {
        full_header_field(header, "accession"): sequence
        for header, sequence in read_fasta(args.full_genome_fasta)
    }
    records_by_gene: dict[str, list[tuple[str, str]]] = defaultdict(list)
    report_rows: list[dict[str, str]] = []
    failures = 0
    for header, sequence in read_fasta(args.per_gene_fasta):
        genotype, gene, accession = parse_gene_header(header)
        records_by_gene[gene].append((header, sequence))
        genome = full_genomes.get(accession, "")
        start = genome.find(sequence) if genome else -1
        status = "PASS" if start >= 0 else "FAIL"
        if status == "FAIL":
            failures += 1
        report_rows.append({
            "Gene": gene,
            "Genotype": genotype,
            "Accession": accession,
            "GeneNALength": str(len(sequence)),
            "FullGenomeStartNA": str(start + 1) if start >= 0 else "",
            "FullGenomeEndNA": str(start + len(sequence)) if start >= 0 else "",
            "Status": status,
        })
    for gene in GENES:
        records = records_by_gene[gene]
        if len(records) != 8:
            raise ValueError(f"Expected eight {gene} records, found {len(records)}")
        write_fasta(args.reference_dir / f"HCV_GT_Refs_{gene}_NA.fasta", records)
    fields = ["Gene", "Genotype", "Accession", "GeneNALength", "FullGenomeStartNA", "FullGenomeEndNA", "Status"]
    args.report_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.report_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"report_csv={args.report_csv.resolve().relative_to(REPO_ROOT)}")
    for gene in GENES:
        print(f"{gene}_output={args.reference_dir.resolve().relative_to(REPO_ROOT) / f'HCV_GT_Refs_{gene}_NA.fasta'}")
    print(f"checked_rows={len(report_rows)}")
    print(f"failed_rows={failures}")
    if failures:
        raise SystemExit("Genotype per-gene NA subset validation failed")


if __name__ == "__main__":
    main()
