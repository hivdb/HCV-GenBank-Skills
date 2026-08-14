#!/usr/bin/env python3
"""Add COMET calls and selected non-COMET priority subtypes to each Ref.csv row."""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


PRIORITY_SUBTYPES = {"1d", "7a", "7b", "8a"}
COVERAGE_COLUMNS = {
    "NS3": "IncludeNS3Pos36_175",
    "NS5A": "IncludeNS5APos26_93",
    "NS5B": "IncludeNS5BPos150_321",
}
REPO_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = REPO_ROOT / "HCVData" / "HCV-all-seq-subtype"


def accession_key(value: str | None) -> str:
    return (value or "").strip().split(".", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-csv", type=Path, default=DATA_DIR / "Ref.csv")
    parser.add_argument("--accessions-csv", type=Path, default=DATA_DIR / "Accessions.csv")
    parser.add_argument("--comet-csv", type=Path, default=DATA_DIR / "all_comet_subtype.csv")
    parser.add_argument(
        "--coverage-csv",
        type=Path,
        action="append",
        help="Coverage CSV to supply priority non-COMET subtypes; repeat as needed.",
    )
    parser.add_argument("--output-csv", type=Path, default=DATA_DIR / "Ref_with_CometSubtypes.csv")
    args = parser.parse_args()

    with args.comet_csv.open(newline="", encoding="utf-8-sig") as handle:
        comet_subtypes = {
            accession_key(row.get("name")): (row.get("subtype") or "").strip()
            for row in csv.DictReader(handle)
            if accession_key(row.get("name"))
        }

    subtypes_by_refid: dict[str, set[str]] = defaultdict(set)
    refid_by_accession: dict[str, str] = {}
    with args.accessions_csv.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            refid = (row.get("RefID") or "").strip()
            accession = accession_key(row.get("Accession"))
            if accession and refid:
                refid_by_accession[accession] = refid
            subtype = comet_subtypes.get(accession, "")
            if refid and subtype and "unassigned" not in subtype.casefold():
                subtypes_by_refid[refid].add(subtype)

    coverage_paths = args.coverage_csv or [
        DATA_DIR / "NS3_AllSeq_NonComet_Coverage.csv",
        DATA_DIR / "NS5A_AllSeq_NonComet_Coverage.csv",
        DATA_DIR / "NS5B_AllSeq_NonComet_Coverage.csv",
    ]
    priority_added = 0
    coverage_by_refid: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    for coverage_path in coverage_paths:
        gene = next((value for value in COVERAGE_COLUMNS if coverage_path.name.startswith(f"{value}_")), None)
        if gene is None:
            raise ValueError(f"Cannot infer gene from coverage filename: {coverage_path}")
        with coverage_path.open(newline="", encoding="utf-8-sig") as handle:
            for row in csv.DictReader(handle):
                accession = accession_key(row.get("Accession"))
                subtype = (row.get("ClosestSubtype") or "").strip()
                refid = refid_by_accession.get(accession, "")
                if refid and subtype.casefold() in PRIORITY_SUBTYPES:
                    before = len(subtypes_by_refid[refid])
                    subtypes_by_refid[refid].add(subtype)
                    priority_added += len(subtypes_by_refid[refid]) - before
                if refid and (row.get("ReferenceOverlapAA") or "").strip():
                    coverage_by_refid[refid][gene].add(accession)

    with args.ref_csv.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames or "RefID" not in reader.fieldnames:
            raise ValueError(f"{args.ref_csv} must contain a RefID column")
        fields = [*reader.fieldnames, "CometSubtypes", *COVERAGE_COLUMNS.values()]
        with args.output_csv.open("w", newline="", encoding="utf-8") as destination:
            writer = csv.DictWriter(destination, fieldnames=fields)
            writer.writeheader()
            for row in reader:
                row["CometSubtypes"] = "; ".join(sorted(subtypes_by_refid.get((row.get("RefID") or "").strip(), set())))
                refid = (row.get("RefID") or "").strip()
                for gene, column in COVERAGE_COLUMNS.items():
                    row[column] = len(coverage_by_refid[refid][gene])
                writer.writerow(row)

    print(f"{args.output_csv} ({priority_added} RefID/subtype priority additions)")


if __name__ == "__main__":
    main()
