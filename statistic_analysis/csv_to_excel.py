#!/usr/bin/env python3
"""Convert a CSV file to an Excel workbook using XlsxWriter."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import xlsxwriter


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument(
        "output_xlsx",
        type=Path,
        nargs="?",
        help="Optional destination; defaults to <gene>_RASs_BySequence.xlsx.",
    )
    parser.add_argument("--sheet-name", default="raw-data")
    return parser.parse_args()


def default_output_path(input_csv: Path) -> Path:
    suffix = "_Profile_Accession_AA_Calls"
    stem = input_csv.stem
    gene = stem.removesuffix(suffix)
    return input_csv.with_name(f"{gene}_RASs_BySequence.xlsx")


def main() -> int:
    args = parse_args()
    output_xlsx = args.output_xlsx or default_output_path(args.input_csv)
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with args.input_csv.open(encoding="utf-8", newline="") as input_handle:
        reader = csv.reader(input_handle)
        with xlsxwriter.Workbook(output_xlsx) as workbook:
            worksheet = workbook.add_worksheet(args.sheet_name)
            header_format = workbook.add_format({"bold": True})
            for row_number, row in enumerate(reader):
                worksheet.write_row(
                    row_number, 0, row, header_format if row_number == 0 else None
                )
                row_count += 1
            worksheet.freeze_panes(1, 0)
            worksheet.autofilter(0, 0, max(row_count - 1, 0), max(len(row) - 1, 0))
    print(f"Wrote {row_count - 1} data rows to {output_xlsx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
