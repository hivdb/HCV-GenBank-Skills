#!/usr/bin/env python3
"""Add aggregated non-unassigned COMET subtypes to each Ref.csv row."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def accession_key(value: str | None) -> str:
    return (value or "").strip().split(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-csv", type=Path, default=Path("Ref.csv"))
    parser.add_argument("--accessions-csv", type=Path, default=Path("Accessions.csv"))
    parser.add_argument("--comet-csv", type=Path, default=Path("all_comet_subtype.csv"))
    parser.add_argument("--output-csv", type=Path, default=Path("Ref_with_CometSubtypes.csv"))
    args = parser.parse_args()

    with args.comet_csv.open(newline="", encoding="utf-8-sig") as handle:
        comet_subtypes = {
            accession_key(row.get("name")): (row.get("subtype") or "").strip()
            for row in csv.DictReader(handle)
            if accession_key(row.get("name"))
        }

    subtypes_by_refid: dict[str, set[str]] = defaultdict(set)
    with args.accessions_csv.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            refid = (row.get("RefID") or "").strip()
            subtype = comet_subtypes.get(accession_key(row.get("Accession")), "")
            if refid and subtype and "unassigned" not in subtype.casefold():
                subtypes_by_refid[refid].add(subtype)

    with args.ref_csv.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "RefID" not in reader.fieldnames:
            raise ValueError(f"{args.ref_csv} must contain a RefID column")
        fields = [*reader.fieldnames, "CometSubtypes"]
        with args.output_csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=fields)
            writer.writeheader()
            for row in reader:
                row["CometSubtypes"] = "; ".join(sorted(subtypes_by_refid.get((row.get("RefID") or "").strip(), set())))
                writer.writerow(row)

    print(args.output_csv)


if __name__ == "__main__":
    main()
