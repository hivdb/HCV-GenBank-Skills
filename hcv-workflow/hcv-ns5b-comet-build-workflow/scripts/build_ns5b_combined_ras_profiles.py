#!/usr/bin/env python3
"""Combine NS5B genotype and subtype RAS profile workbooks by genotype."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path

from openpyxl import Workbook, load_workbook
from openpyxl.cell.rich_text import CellRichText, TextBlock
from openpyxl.cell.text import InlineFont
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter


VARIANT_RE = re.compile(r"([A-Z*])(\d+(?:\.\d+)?)")
TOTAL_SEQUENCE_RE = re.compile(r"\(\s*(\d+)\s*,")
MIN_TOTAL_SEQUENCES = 10
FULL_PROFILE_SUMMARY_MIN_TOTAL_SEQUENCES = 10


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Combine NS5B GT and subtype RAS profiles, retaining subtype amino-acid "
            "variants strictly above a frequency threshold."
        )
    )
    parser.add_argument("--gt-ras-profile-workbook", required=True)
    parser.add_argument("--subtype-ras-profile-workbook", required=True)
    parser.add_argument("--output-xlsx", default="outputs/NS5B_Combined_RAS_Profiles.xlsx")
    parser.add_argument("--subtype-frequency-threshold", type=float, default=1.0)
    return parser.parse_args()


def read_profile_workbook(path: Path) -> tuple[list[object], list[list[object]]]:
    workbook = load_workbook(path, read_only=True, data_only=True)
    worksheet = workbook.active
    rows = list(worksheet.iter_rows(values_only=True))
    workbook.close()
    if not rows:
        raise RuntimeError(f"Profile workbook is empty: {path}")
    return list(rows[0]), [list(row) for row in rows[1:] if row and row[0]]


def genotype_from_label(label: object) -> str | None:
    match = re.match(r"GT(\d+)(?:\s|_)", str(label))
    return match.group(1) if match else None


def has_minimum_total_sequences(label: object) -> bool:
    if genotype_from_label(label) in {"7", "8"}:
        return True
    match = TOTAL_SEQUENCE_RE.search(str(label))
    return match is not None and int(match.group(1)) >= MIN_TOTAL_SEQUENCES


def total_sequences(label: object) -> int | None:
    match = TOTAL_SEQUENCE_RE.search(str(label))
    return int(match.group(1)) if match is not None else None


def variants_to_rich_text(
    value: object, threshold: float | None = None, excluded_amino_acids: set[str] | None = None
) -> CellRichText | str:
    """Render amino-acid frequencies with superscript frequencies, as in RAS profiles."""
    variants = [
        (amino_acid, frequency)
        for amino_acid, frequency in VARIANT_RE.findall(str(value or ""))
        if (threshold is None or float(frequency) >= threshold)
        and amino_acid not in (excluded_amino_acids or set())
    ]
    if not variants:
        return ""
    parts: list[str | TextBlock] = []
    for amino_acid, frequency in variants:
        parts.extend((amino_acid, TextBlock(InlineFont(vertAlign="superscript"), frequency)))
    return CellRichText(*parts)


def most_frequent_variant(value: object) -> tuple[str, str] | None:
    variants = VARIANT_RE.findall(str(value or ""))
    if not variants:
        return None
    return max(variants, key=lambda variant: float(variant[1]))


def most_frequent_variant_to_rich_text(value: object) -> CellRichText | str:
    variant = most_frequent_variant(value)
    if variant is None:
        return ""
    amino_acid, frequency = variant
    return CellRichText(amino_acid, TextBlock(InlineFont(vertAlign="superscript"), frequency))


def mean_diff(value_rows: list[object], gt_variants: list[tuple[str, str] | None], threshold: float) -> float:
    """Sum displayed non-consensus amino-acid percentages and express as a decimal."""
    displayed_percent = sum(
        float(frequency)
        for value, gt_variant in zip(value_rows, gt_variants)
        for amino_acid, frequency in VARIANT_RE.findall(str(value or ""))
        if float(frequency) >= threshold and (gt_variant is None or amino_acid != gt_variant[0])
    )
    return displayed_percent / 100.0


def write_combined_workbook(
    output_path: Path,
    positions: list[object],
    gt_rows: list[list[object]],
    subtype_rows: list[list[object]],
    threshold: float,
) -> tuple[int, int]:
    gt_by_number = {
        genotype: row
        for row in gt_rows
        if (genotype := genotype_from_label(row[0])) is not None
        and has_minimum_total_sequences(row[0])
    }
    subtypes_by_gt: dict[str, list[list[object]]] = defaultdict(list)
    for row in subtype_rows:
        genotype = genotype_from_label(row[0])
        if genotype is not None and has_minimum_total_sequences(row[0]):
            subtypes_by_gt[genotype].append(row)

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Combined_RAS_Profile"
    header_fill = PatternFill(fill_type="solid", fgColor="F2F2F2")
    gt_fill = PatternFill(fill_type="solid", fgColor="D9EAF7")
    bold = Font(bold=True)

    worksheet.append([*positions, "MeanDiff"])
    for cell in worksheet[1]:
        cell.fill = header_fill
        cell.font = bold

    output_rows = 0
    for genotype in sorted(gt_by_number, key=int):
        gt_row = gt_by_number[genotype]
        worksheet.append(
            [gt_row[0]] + [most_frequent_variant_to_rich_text(value) for value in gt_row[1:]] + [None]
        )
        output_rows += 1
        for cell in worksheet[worksheet.max_row]:
            cell.fill = gt_fill
            cell.font = bold

        gt_amino_acids = [most_frequent_variant(value) for value in gt_row[1:]]
        for subtype_row in subtypes_by_gt.get(genotype, []):
            subtype_values = subtype_row[1:]
            worksheet.append(
                [subtype_row[0]]
                + [
                    variants_to_rich_text(
                        value,
                        threshold,
                        {gt_variant[0]} if gt_variant is not None else None,
                    )
                    for value, gt_variant in zip(subtype_values, gt_amino_acids)
                ]
                + [mean_diff(subtype_values, gt_amino_acids, threshold)]
            )
            output_rows += 1
            worksheet.cell(worksheet.max_row, worksheet.max_column).number_format = "0.0"
        worksheet.append([])

    for row in worksheet.iter_rows():
        for cell in row:
            cell.alignment = Alignment(horizontal="center")
    worksheet.column_dimensions["A"].width = 24
    for column in range(2, worksheet.max_column + 1):
        worksheet.column_dimensions[get_column_letter(column)].width = 12

    output_path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output_path)
    return len(gt_by_number), output_rows


def main() -> int:
    args = parse_args()
    gt_path = Path(args.gt_ras_profile_workbook).expanduser()
    subtype_path = Path(args.subtype_ras_profile_workbook).expanduser()
    output_path = Path(args.output_xlsx).expanduser()
    if not gt_path.is_file() or not subtype_path.is_file():
        raise RuntimeError("Both GT and subtype RAS profile workbooks are required")

    positions, gt_rows = read_profile_workbook(gt_path)
    subtype_positions, subtype_rows = read_profile_workbook(subtype_path)
    if positions != subtype_positions:
        raise RuntimeError("GT and subtype RAS profile workbooks have different position rows")
    genotype_count, output_rows = write_combined_workbook(
        output_path, positions, gt_rows, subtype_rows, args.subtype_frequency_threshold
    )
    full_profile_subtype_count = len(subtype_rows)
    full_profile_subtype_at_least_10_count = sum(
        (count := total_sequences(row[0])) is not None
        and count >= FULL_PROFILE_SUMMARY_MIN_TOTAL_SEQUENCES
        for row in subtype_rows
    )
    print(
        json.dumps(
            {
                "output_xlsx": str(output_path.resolve()),
                "genotype_count": genotype_count,
                "profile_row_count": output_rows,
                "full_profile_subtype_count": full_profile_subtype_count,
                "full_profile_subtype_at_least_10_sequence_count": full_profile_subtype_at_least_10_count,
                "subtype_frequency_threshold": args.subtype_frequency_threshold,
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
