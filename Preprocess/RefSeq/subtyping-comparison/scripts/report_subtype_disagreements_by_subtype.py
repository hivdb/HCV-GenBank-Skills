#!/usr/bin/env python3
"""Summarize non-unanimous subtype calls by the per-gene COMET subtype."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

from export_csv_directory_to_excel import write_csv_directory_workbook

REPO_ROOT = Path(__file__).resolve().parents[4]
STEP11_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison" / "11-report-subtype-agreement"
OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison" / "12-report-subtype-disagreements"
GENES = ("NS3", "NS5A", "NS5B")
DISAGREEMENT_STATUSES = ("AllFiveDisagree", "AvailableCallsDisagree")
COMET_SUBTYPE_COLUMN = "CometPerGeneSubtype"
OTHER_SUBTYPE_COLUMNS = (
    "PerGeneSubtype",
    "FullGenomeSubtype",
    "CometFullGenomeSubtype",
    "GenBankSubtype",
)


def subtype_key(value: str) -> str:
    """Normalize a subtype value used in the report."""
    return value.strip().upper() or "(BLANK)"


def build_report(gene: str) -> None:
    input_csv = STEP11_DIR / f"{gene}_Subtype_Agreement_By_Accession.csv"
    output_csv = OUTPUT_DIR / f"{gene}_Subtype_Disagreement_Summary.csv"
    counts: dict[str, Counter[str]] = {}
    other_subtypes: dict[str, set[str]] = {}

    with input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        required = {COMET_SUBTYPE_COLUMN, *OTHER_SUBTYPE_COLUMNS, "AgreementStatus"}
        if not reader.fieldnames or (missing := required - set(reader.fieldnames)):
            raise ValueError(f"{input_csv} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            status = row.get("AgreementStatus", "").strip()
            if status not in DISAGREEMENT_STATUSES:
                continue
            comet_subtype = subtype_key(row.get(COMET_SUBTYPE_COLUMN, ""))
            counts.setdefault(comet_subtype, Counter())[status] += 1
            other_subtypes.setdefault(comet_subtype, set()).update(
                subtype_key(row.get(column, ""))
                for column in OTHER_SUBTYPE_COLUMNS
                if row.get(column, "").strip()
            )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    for report in ("MoreThanHalfAgree", "MoreThanHalfDisagree"):
        split_csv = OUTPUT_DIR / f"{gene}_Subtype_{report}.csv"
        if split_csv.exists():
            split_csv.unlink()
    with output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=(
                "CometSubtype",
                *DISAGREEMENT_STATUSES,
                "TotalDisagree",
                "OtherMethodSubtypes",
            ),
            lineterminator="\n",
        )
        writer.writeheader()
        for subtype in sorted(counts):
            writer.writerow(
                {
                    "CometSubtype": subtype,
                    "AllFiveDisagree": counts[subtype]["AllFiveDisagree"],
                    "AvailableCallsDisagree": counts[subtype]["AvailableCallsDisagree"],
                    "TotalDisagree": sum(counts[subtype].values()),
                    "OtherMethodSubtypes": ";".join(sorted(other_subtypes[subtype])),
                }
            )
    print(f"{gene} disagreement COMET subtypes: {len(counts):,}")
    print(f"{gene} output: {output_csv}")


def main() -> None:
    for gene in GENES:
        build_report(gene)
    output_xlsx = OUTPUT_DIR / "Subtype_Disagreement_Reports.xlsx"
    write_csv_directory_workbook(OUTPUT_DIR, output_xlsx)
    print(f"Excel output: {output_xlsx}")


if __name__ == "__main__":
    main()
