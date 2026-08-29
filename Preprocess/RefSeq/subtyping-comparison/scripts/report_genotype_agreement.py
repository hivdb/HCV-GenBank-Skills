#!/usr/bin/env python3
"""Report agreement among five genotype sources in Step 8 merged results."""

from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STEP8_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison" / "08-merge-ns3-subtyping-sources"
OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison" / "10-report-genotype-agreement"
GENES = ("NS3", "NS5A", "NS5B")
METHOD_COLUMNS = (
    "PerGeneGenotype",
    "FullGenomeGenotype",
    "CometPerGeneGenotype",
    "CometFullGenomeGenotype",
    "GenBankGenotype",
)


def genotype(value: str) -> str:
    return value.strip().upper()


def agreement_status(present_count: int, distinct_count: int) -> str:
    if present_count < 2:
        return "InsufficientCalls"
    if present_count == len(METHOD_COLUMNS):
        return "AllFiveAgree" if distinct_count == 1 else "AllFiveDisagree"
    return "AvailableCallsAgree" if distinct_count == 1 else "AvailableCallsDisagree"


def build_reports(gene: str) -> None:
    input_csv = STEP8_DIR / f"{gene}_Subtyping_Sources_Merged.csv"
    detail_csv = OUTPUT_DIR / f"{gene}_Genotype_Agreement_By_Accession.csv"
    summary_csv = OUTPUT_DIR / f"{gene}_Genotype_Agreement_Summary.csv"

    with input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {input_csv}")
        required = {"Accession", *METHOD_COLUMNS}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{input_csv} is missing columns: {', '.join(sorted(missing))}")

        detail_rows: list[dict[str, str | int]] = []
        statuses: Counter[str] = Counter()
        for row in reader:
            calls = {column: genotype(row.get(column, "")) for column in METHOD_COLUMNS}
            present_calls = [value for value in calls.values() if value]
            unique_calls = sorted(set(present_calls))
            status = agreement_status(len(present_calls), len(unique_calls))
            statuses[status] += 1
            detail_rows.append(
                {
                    "Accession": row["Accession"].strip(),
                    **calls,
                    "PresentMethodCount": len(present_calls),
                    "DistinctGenotypeCount": len(unique_calls),
                    "GenotypesObserved": ";".join(unique_calls),
                    "AgreementStatus": status,
                }
            )

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    detail_fields = (
        "Accession",
        *METHOD_COLUMNS,
        "PresentMethodCount",
        "DistinctGenotypeCount",
        "GenotypesObserved",
        "AgreementStatus",
    )
    with detail_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=detail_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(detail_rows)

    total = len(detail_rows)
    summary_rows = [
        {
            "Gene": gene,
            "AgreementStatus": status,
            "AccessionCount": statuses[status],
            "PercentOfAccessions": f"{(statuses[status] / total * 100) if total else 0:.2f}",
        }
        for status in (
            "AllFiveAgree",
            "AllFiveDisagree",
            "AvailableCallsAgree",
            "AvailableCallsDisagree",
            "InsufficientCalls",
        )
    ]
    with summary_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=("Gene", "AgreementStatus", "AccessionCount", "PercentOfAccessions"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"{gene} accessions: {total:,}")
    for row in summary_rows:
        print(f"{gene} {row['AgreementStatus']}: {row['AccessionCount']:,}")
    print(f"Detail output: {detail_csv}")
    print(f"Summary output: {summary_csv}")


def main() -> None:
    for gene in GENES:
        build_reports(gene)


if __name__ == "__main__":
    main()
