#!/usr/bin/env python3
"""Select non-COMET assignments that take priority over COMET calls."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--local-assignments-csv", required=True)
    parser.add_argument("--fasta-dir", required=True, help="RefID-prefixed FASTA pool used to resolve priority accessions.")
    parser.add_argument("--comet-unassigned-csv", required=True, help="Current-workflow COMET missing/unassigned report.")
    parser.add_argument("--output-csv", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with Path(args.comet_unassigned_csv).open(encoding="utf-8-sig", newline="") as source:
        comet_unassigned_accessions = {
            (row.get("accession") or "").strip().split(".", 1)[0]
            for row in csv.DictReader(source)
            if (row.get("column_name") or "").strip().lower() == "unassigned"
        }
    refid_by_accession: dict[str, tuple[str, str]] = {}
    for fasta_path in sorted(Path(args.fasta_dir).glob("*.fasta")):
        refid, _, refname = fasta_path.stem.partition("_")
        for line in fasta_path.read_text(encoding="utf-8").splitlines():
            if line.startswith(">"):
                accession = line[1:].strip().split(maxsplit=1)[0].split(".", 1)[0]
                refid_by_accession.setdefault(accession, (refid, refname))
    with Path(args.local_assignments_csv).open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise RuntimeError(f"No header found in {args.local_assignments_csv}")
        rows = []
        for row in reader:
            accession = (row.get("accession") or "").strip()
            genotype = (row.get("genotype") or "").strip().removeprefix("GT")
            subtype = (row.get("subtype") or "").strip().lower()
            if not accession or not genotype.isdigit() or not subtype:
                continue
            accession_key = accession.split(".", 1)[0]
            if subtype == "1d" or genotype in {"7", "8"} or subtype.startswith(("7", "8")) or accession_key in comet_unassigned_accessions:
                rows.append({"Accession": accession, "ClosestGenotype": genotype, "ClosestSubtype": subtype})
        fieldnames = ["Accession", "ClosestGenotype", "ClosestSubtype", "RefID", "RefName"]
    for row in rows:
        refid, refname = refid_by_accession.get((row.get("Accession") or "").strip().split(".", 1)[0], ("", ""))
        row["RefID"], row["RefName"] = refid, refname
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output_csv": str(output.resolve()), "priority_assignment_count": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
