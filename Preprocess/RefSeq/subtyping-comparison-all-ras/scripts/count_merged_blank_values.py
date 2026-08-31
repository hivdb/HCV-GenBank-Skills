#!/usr/bin/env python3
"""Count blank values in each column of the Step 8 merged subtype results."""

from __future__ import annotations

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
STEP8_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras" / "08-merge-ns3-subtyping-sources"
OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras" / "09-count-merged-blank-values"
GENES = ("NS3", "NS5A", "NS5B")


def count_blanks(input_csv: Path, output_csv: Path, gene: str) -> None:
    with input_csv.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        if not reader.fieldnames:
            raise ValueError(f"No header found in {input_csv}")
        fieldnames = reader.fieldnames
        blank_counts = {column: 0 for column in fieldnames}
        total_rows = 0
        for row in reader:
            total_rows += 1
            for column in fieldnames:
                if not row.get(column, "").strip():
                    blank_counts[column] += 1

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(
            destination,
            fieldnames=("Gene", "Column", "BlankCount", "TotalRows"),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(
            {
                "Gene": gene,
                "Column": column,
                "BlankCount": blank_counts[column],
                "TotalRows": total_rows,
            }
            for column in fieldnames
        )

    print(f"{gene} rows: {total_rows:,}")
    for column in fieldnames:
        print(f"{gene} {column} blank values: {blank_counts[column]:,}")
    print(f"Output: {output_csv}")


def main() -> None:
    for gene in GENES:
        count_blanks(
            STEP8_DIR / f"{gene}_Subtyping_Sources_Merged.csv",
            OUTPUT_DIR / f"{gene}_Blank_Value_Counts.csv",
            gene,
        )


if __name__ == "__main__":
    main()
