#!/usr/bin/env python3
"""Merge selected active COMET NA distance matrices into one paper-ready workbook."""

from __future__ import annotations

import argparse
from copy import copy
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.styles import PatternFill


REPO_ROOT = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_DIR = REPO_ROOT / "outputs/hcv-na-distance-matrix-merge"
SOURCES = (
    (
        "NS3 One RAS",
        REPO_ROOT
        / "outputs/comet-NS3-one-ras/30_build-paired-distance-matrices/NS3_GT_NA_Distance_RAS.xlsx",
        REPO_ROOT
        / "outputs/comet-NS3-one-ras/30_build-paired-distance-matrices/NS3_Subtype_NA_Distance_RAS.xlsx",
        "DDEBF7",
    ),
    (
        "NS5A One RAS",
        REPO_ROOT
        / "outputs/comet-NS5A-one-ras/29_build-paired-distance-matrices/NS5A_GT_NA_Distance_RAS.xlsx",
        REPO_ROOT
        / "outputs/comet-NS5A-one-ras/29_build-paired-distance-matrices/NS5A_Subtype_NA_Distance_RAS.xlsx",
        "E2F0D9",
    ),
    (
        "NS5B Five or more",
        REPO_ROOT
        / "outputs/comet-NS5B-position-282-four-ras/29_build-paired-distance-matrices/NS5B_GT_NA_Distance_RAS.xlsx",
        REPO_ROOT
        / "outputs/comet-NS5B-position-282-four-ras/29_build-paired-distance-matrices/NS5B_Subtype_NA_Distance_RAS.xlsx",
        "FCE4D6",
    ),
)


def copy_vertical_block(
    source,
    destination,
    start_row: int,
    start_column: int,
    title: str,
    background_fill: PatternFill,
) -> int:
    """Copy a titled worksheet block and return its ending row."""
    title_cell = destination.cell(start_row, start_column, title)
    title_font = copy(source["A1"].font)
    title_font.bold = True
    title_cell.font = title_font
    title_cell.fill = copy(background_fill)
    for row in source.iter_rows():
        for cell in row:
            target = destination.cell(
                start_row + cell.row,
                start_column + cell.column - 1,
                cell.value,
            )
            if cell.has_style:
                target._style = copy(cell._style)
            target.fill = copy(background_fill)
            if cell.number_format:
                target.number_format = cell.number_format
    for column, dimension in source.column_dimensions.items():
        target_column = destination.cell(
            1, start_column + source[column][0].column - 1
        ).column_letter
        current_width = destination.column_dimensions[target_column].width
        if dimension.width and (
            current_width is None or dimension.width > current_width
        ):
            destination.column_dimensions[target_column].width = dimension.width
    return start_row + source.max_row


def relevant_subtype_sheets(workbook) -> list[str]:
    return [
        name
        for name in workbook.sheetnames
        if name.startswith("GT") and not name.endswith("_counts")
    ]


def merge_distance_matrices(output_path: Path) -> None:
    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "GT_NA_Distance"
    column = 1
    for gene_title, gt_path, subtype_path, fill_color in SOURCES:
        gt_book = load_workbook(gt_path, data_only=False)
        subtype_book = load_workbook(subtype_path, data_only=False)
        try:
            gt_sheet = gt_book["distance_matrix"]
            subtype_sheets = [
                subtype_book[name] for name in relevant_subtype_sheets(subtype_book)
            ]
            background_fill = PatternFill("solid", fgColor=fill_color)
            gene_width = max(
                [gt_sheet.max_column]
                + [subtype_sheet.max_column for subtype_sheet in subtype_sheets]
            )

            gene_cell = worksheet.cell(1, column, gene_title)
            gene_font = copy(gt_sheet["A1"].font)
            gene_font.bold = True
            gene_cell.font = gene_font
            for gene_column in range(column, column + gene_width):
                worksheet.cell(1, gene_column).fill = copy(background_fill)

            row = (
                copy_vertical_block(
                    gt_sheet,
                    worksheet,
                    2,
                    column,
                    "Genotype",
                    background_fill,
                )
                + 2
            )
            for subtype_sheet in subtype_sheets:
                row = (
                    copy_vertical_block(
                        subtype_sheet,
                        worksheet,
                        row,
                        column,
                        f"{subtype_sheet.title} Subtypes",
                        background_fill,
                    )
                    + 2
                )

            column += gene_width + 1
        finally:
            gt_book.close()
            subtype_book.close()
    workbook.save(output_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "NA_Distance_Matrics.xlsx"
    merge_distance_matrices(output)
    print(f"output={output.resolve()}")


if __name__ == "__main__":
    main()
