#!/usr/bin/env python3
"""Create the NS5B genotype workbook directly from Comet assignments."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from openpyxl import Workbook


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", required=True)
    parser.add_argument("--comet-genotype-csv", required=True)
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


def main() -> int:
    args = parse_args()
    assignments = load_comet(Path(args.comet_genotype_csv))
    rows: list[tuple[str, str, str, str]] = []
    for fasta_path in sorted(Path(args.fasta_dir).glob("*.fasta")):
        refid, _, refname = fasta_path.stem.partition("_")
        for accession in fasta_accessions(fasta_path):
            genotype = assignments.get(accession) or assignments.get(accession.split(".", 1)[0])
            if genotype:
                rows.append((refid, refname, accession, genotype))

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "NS5B_GT_AllStudies.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "NS5B_GT_AllStudies"
    sheet.append(["RefID", "RefName", "GenBankAccession", "BestGT", "BestGTAssignmentSource"])
    for row in rows:
        sheet.append([*row, "Comet"])
    workbook.save(output_path)
    print(json.dumps({"combined_xlsx": str(output_path.resolve()), "master_row_count": len(rows), "comet_best_gt_count": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
