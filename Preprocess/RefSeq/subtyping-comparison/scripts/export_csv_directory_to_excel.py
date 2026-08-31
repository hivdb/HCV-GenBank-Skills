#!/usr/bin/env python3
"""Export every CSV in a report directory to a worksheet in one Excel file."""

from __future__ import annotations

import csv
from pathlib import Path

import xlsxwriter


def sheet_name(csv_path: Path, existing_names: set[str]) -> str:
    """Create a unique Excel-compatible worksheet name from a CSV filename."""
    base = csv_path.stem.replace("_Agreement_By_Accession", " Detail").replace(
        "_Agreement_Summary", " Summary"
    ).replace("_Subtype_Disagreement_Summary", " Summary").replace("_", " ")
    candidate = base[:31]
    suffix = 2
    while candidate in existing_names:
        suffix_text = f" ({suffix})"
        candidate = f"{base[:31 - len(suffix_text)]}{suffix_text}"
        suffix += 1
    return candidate


def write_csv_directory_workbook(csv_directory: Path, output_xlsx: Path) -> None:
    """Write sorted CSV files in *csv_directory* to separate workbook sheets."""
    csv_files = sorted(csv_directory.glob("*.csv"))
    if not csv_files:
        raise ValueError(f"No CSV files found in {csv_directory}")

    workbook = xlsxwriter.Workbook(output_xlsx)
    header_format = workbook.add_format({"bold": True, "bg_color": "#D9EAF7"})
    used_names: set[str] = set()
    try:
        for csv_path in csv_files:
            name = sheet_name(csv_path, used_names)
            used_names.add(name)
            worksheet = workbook.add_worksheet(name)
            worksheet.freeze_panes(1, 0)
            with csv_path.open(encoding="utf-8-sig", newline="") as source:
                reader = csv.reader(source)
                header = next(reader, None)
                if header is None:
                    continue
                worksheet.write_row(0, 0, header, header_format)
                worksheet.autofilter(0, 0, 0, len(header) - 1)
                for row_number, row in enumerate(reader, start=1):
                    worksheet.write_row(row_number, 0, row)
                for column, value in enumerate(header):
                    worksheet.set_column(column, column, min(max(len(value) + 2, 12), 28))
    finally:
        workbook.close()

