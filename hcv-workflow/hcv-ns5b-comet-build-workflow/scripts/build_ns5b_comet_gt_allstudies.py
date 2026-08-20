#!/usr/bin/env python3
"""Create the NS5B genotype workbook directly from Comet assignments."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from openpyxl import Workbook, load_workbook


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", required=True)
    parser.add_argument("--comet-genotype-csv", required=True)
    parser.add_argument("--priority-assignments-csv", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def fasta_accessions(path: Path) -> list[str]:
    return [line[1:].strip().split(maxsplit=1)[0] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]


def load_comet(path: Path) -> dict[str, str]:
    assignments: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("accession") or "").strip()
            genotype = (row.get("genotype") or "").strip()
            if accession and genotype:
                assignments[accession] = genotype
                assignments.setdefault(accession.split(".", 1)[0], genotype)
    return assignments


def load_priority_genotypes(path: Path) -> dict[str, tuple[str, str, str, str]]:
    assignments: dict[str, tuple[str, str, str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = str(row.get("Accession") or "").strip()
            genotype = str(row.get("ClosestGenotype") or "").strip().lower()
            subtype = str(row.get("ClosestSubtype") or "").strip().lower()
            if accession:
                assignments[accession.split(".", 1)[0]] = ("", "", accession, genotype)
    return assignments


def main() -> int:
    args = parse_args()
    assignments = load_comet(Path(args.comet_genotype_csv))
    noncomet_priority_genotypes = load_priority_genotypes(Path(args.priority_assignments_csv))
    rows: list[tuple[str, str, str, str, str]] = []
    seen_accessions: set[str] = set()
    override_count = 0
    for fasta_path in sorted(Path(args.fasta_dir).glob("*.fasta")):
        refid, _, refname = fasta_path.stem.partition("_")
        for accession in fasta_accessions(fasta_path):
            accession_key = accession.split(".", 1)[0]
            priority_assignment = noncomet_priority_genotypes.get(accession_key)
            genotype = priority_assignment[3] if priority_assignment else assignments.get(accession) or assignments.get(accession_key)
            if genotype:
                source = "Non-Comet priority genotype override" if priority_assignment else "Comet"
                rows.append((refid, refname, accession, genotype, source))
                seen_accessions.add(accession_key)
                if priority_assignment:
                    override_count += 1

    addition_count = 0
    for accession_key, (refid, refname, accession, genotype) in noncomet_priority_genotypes.items():
        if accession_key not in seen_accessions:
            rows.append((refid, refname, accession, genotype, "Non-Comet priority genotype addition"))
            addition_count += 1

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "NS5B_GT_AllStudies.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "NS5B_GT_AllStudies"
    sheet.append(["RefID", "RefName", "GenBankAccession", "BestGT", "BestGTAssignmentSource"])
    for row in rows:
        sheet.append(row)
    workbook.save(output_path)
    print(json.dumps({"combined_xlsx": str(output_path.resolve()), "master_row_count": len(rows), "noncomet_priority_genotype_override_count": override_count, "noncomet_priority_genotype_addition_count": addition_count}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
