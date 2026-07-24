#!/usr/bin/env python3
"""Prepare authoritative NS5B Comet assignments, optionally filtering FASTA records."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


SUBTYPE_RE = re.compile(r"^([1-8])([a-z][a-z0-9]*)$", re.IGNORECASE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comet-csv", required=True)
    parser.add_argument("--fasta-dir", required=True)
    parser.add_argument("--genotype-output-csv", required=True)
    parser.add_argument("--subtype-output-csv", required=True)
    parser.add_argument("--not-found-output-csv", required=True)
    parser.add_argument(
        "--not-found-fasta-output",
        required=True,
        help="FASTA file containing raw records missing from Comet or marked unassigned.",
    )
    parser.add_argument(
        "--remove-unassigned",
        action="store_true",
        help="Remove missing and unassigned records after reporting their counts.",
    )
    return parser.parse_args()


def accession_from_header(header: str) -> str:
    return header[1:].strip().split(maxsplit=1)[0]


def read_fasta(path: Path) -> list[list[str]]:
    records: list[list[str]] = []
    current: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(">"):
            if current:
                records.append(current)
            current = [line]
        else:
            current.append(line)
    if current:
        records.append(current)
    return records


def load_comet(path: Path) -> dict[str, tuple[str, str] | None]:
    calls: dict[str, tuple[str, str] | None] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("name") or "").strip()
            subtype = (row.get("subtype") or "").strip().lower()
            if not accession or (row.get("virus") or "").strip().upper() != "HCV":
                continue
            match = SUBTYPE_RE.fullmatch(subtype)
            call = (match.group(1), subtype) if match else None
            calls[accession] = call
            calls.setdefault(accession.split(".", 1)[0], call)
    return calls


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession", "genotype", "subtype", "column_name"])
        writer.writeheader()
        writer.writerows(rows)


def write_fasta(path: Path, records: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(line for record in records for line in record), encoding="utf-8")


def main() -> int:
    args = parse_args()
    fasta_dir = Path(args.fasta_dir)
    calls = load_comet(Path(args.comet_csv))
    genotype_rows: list[dict[str, str]] = []
    subtype_rows: list[dict[str, str]] = []
    missing_rows: list[dict[str, str]] = []
    missing_or_unassigned_records: list[list[str]] = []
    kept_records_by_path: dict[Path, list[list[str]]] = {}
    total_accessions = assigned = unassigned = missing = 0

    for fasta_path in sorted(fasta_dir.glob("*.fasta")):
        kept: list[list[str]] = []
        for record in read_fasta(fasta_path):
            accession = accession_from_header(record[0])
            total_accessions += 1
            call = calls.get(accession) or calls.get(accession.split(".", 1)[0])
            if call is None:
                status = "unassigned" if accession in calls or accession.split(".", 1)[0] in calls else "missing"
                unassigned += status == "unassigned"
                missing += status == "missing"
                missing_rows.append({"accession": accession, "genotype": "", "subtype": "", "column_name": status})
                missing_or_unassigned_records.append(record)
                continue
            genotype, subtype = call
            assigned += 1
            kept.append(record)
            genotype_rows.append({"accession": accession, "genotype": genotype, "subtype": "", "column_name": "Comet NS5B"})
            subtype_rows.append({"accession": accession, "genotype": genotype, "subtype": subtype, "column_name": "Comet NS5B"})
        kept_records_by_path[fasta_path] = kept

    print(f"staged_fasta_accession_count={total_accessions}")
    print(f"comet_assigned_accession_count={assigned}")
    print(f"comet_missing_accession_count={missing}")
    print(f"comet_unassigned_accession_count={unassigned}")
    if args.remove_unassigned:
        print("removing_missing_or_unassigned_accessions=true")
        for fasta_path, kept in kept_records_by_path.items():
            fasta_path.write_text("".join(line for record in kept for line in record), encoding="utf-8")
    else:
        print("removing_missing_or_unassigned_accessions=false")
    write_csv(Path(args.genotype_output_csv), genotype_rows)
    write_csv(Path(args.subtype_output_csv), subtype_rows)
    write_csv(Path(args.not_found_output_csv), missing_rows)
    write_fasta(Path(args.not_found_fasta_output), missing_or_unassigned_records)
    print(f"comet_genotype_assignments={Path(args.genotype_output_csv)}")
    print(f"comet_subtype_assignments={Path(args.subtype_output_csv)}")
    print(f"comet_not_found_or_unassigned={Path(args.not_found_output_csv)}")
    print(f"comet_not_found_or_unassigned_fasta={Path(args.not_found_fasta_output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
