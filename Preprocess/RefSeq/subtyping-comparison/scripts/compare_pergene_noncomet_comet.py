#!/usr/bin/env python3
"""Compare COMET, per-gene BLAST, and full-genome BLAST subtype calls."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path

import xlsxwriter


GENES = ("NS3", "NS5A", "NS5B")
COMPARISONS = {
    "Comet_vs_PerGene": ("Comet", "PerGene"),
    "Comet_vs_FullGenome": ("Comet", "FullGenome"),
    "PerGene_vs_FullGenome": ("PerGene", "FullGenome"),
    "Comet_vs_PerGene_vs_FullGenome": ("Comet", "PerGene", "FullGenome"),
}
REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_PERGENE_DIR = REPO_ROOT / "HCVData" / "nonComet-PerGene"
DEFAULT_FULLGENOME_DIR = REPO_ROOT / "HCVData" / "nonComet-Full-genome"
DEFAULT_COMET_CSV = REPO_ROOT / "HCVData" / "Comet-Full-genome" / "all_comet_subtype.csv"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "HCVData" / "subtyping-comparison"


def accession_key(value: object) -> str:
    return str(value or "").strip().split(".", 1)[0].upper()


def normalized(value: object) -> str:
    return str(value or "").strip().lower().removeprefix("gt")


def read_comet_calls(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"name", "subtype"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        calls: dict[str, str] = {}
        for row in reader:
            accession = accession_key(row.get("name"))
            subtype = normalized(row.get("subtype"))
            if not accession:
                continue
            if accession in calls and calls[accession] != subtype:
                raise ValueError(f"{path} has conflicting COMET subtypes for {accession}")
            calls[accession] = subtype
    return calls


def read_blast_calls(path: Path) -> dict[str, str]:
    """Return subtype calls only for rows with a non-blank genotype assignment."""
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Accession", "ClosestGenotype", "ClosestSubtype"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        calls: dict[str, str] = {}
        for row in reader:
            if not normalized(row.get("ClosestGenotype")):
                continue
            accession = accession_key(row.get("Accession"))
            subtype = normalized(row.get("ClosestSubtype"))
            if accession in calls and calls[accession] != subtype:
                raise ValueError(f"{path} has conflicting subtype calls for {accession}")
            if accession:
                calls[accession] = subtype
    return calls


def compare_calls(
    gene: str,
    condition: str,
    methods: tuple[str, ...],
    calls_by_method: dict[str, dict[str, str]],
) -> tuple[dict[str, int | str], list[dict[str, str]]]:
    # Per-gene BLAST defines the accession universe for every comparison. This
    # keeps NS3, NS5A, and NS5B comparisons scoped to the corresponding
    # PerGene file instead of all FullGenome or COMET accessions.
    accessions = set(calls_by_method["PerGene"])
    for method in methods:
        accessions.intersection_update(calls_by_method[method])
    same_count = 0
    differences: list[dict[str, str]] = []
    for accession in sorted(accessions):
        selected_calls = {method: calls_by_method[method][accession] for method in methods}
        if len(set(selected_calls.values())) == 1:
            same_count += 1
            continue
        differences.append(
            {
                "Gene": gene,
                "Accession": accession,
                "CometSubtype": calls_by_method["Comet"].get(accession, ""),
                "BlastPerGeneSubtype": calls_by_method["PerGene"].get(accession, ""),
                "BlastFullGenomeSubtype": calls_by_method["FullGenome"].get(accession, ""),
                "CompareCondition": condition,
            }
        )
    return (
        {
            "Gene": gene,
            "SameCount": same_count,
            "DifferentCount": len(differences),
            "CompareCondition": condition,
        },
        differences,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pergene-dir", type=Path, default=DEFAULT_PERGENE_DIR)
    parser.add_argument("--fullgenome-dir", type=Path, default=DEFAULT_FULLGENOME_DIR)
    parser.add_argument("--comet-csv", type=Path, default=DEFAULT_COMET_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def write_csv(path: Path, fields: list[str], rows: list[dict[str, int | str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_count_workbook(path: Path, fields: list[str], rows: list[dict[str, int | str]]) -> None:
    """Write subtype-difference counts with one comparison condition per sheet."""
    with xlsxwriter.Workbook(str(path)) as workbook:
        header_format = workbook.add_format({"bold": True, "bg_color": "#5B9BD5"})
        for condition in COMPARISONS:
            worksheet = workbook.add_worksheet(condition)
            worksheet.write_row(0, 0, fields, header_format)
            condition_rows = [
                row for row in rows if row["CompareCondition"] == condition
            ]
            for row_number, row in enumerate(condition_rows, start=1):
                worksheet.write_row(row_number, 0, [row[field] for field in fields])
            worksheet.freeze_panes(1, 0)
            for column_number, field in enumerate(fields):
                values = [field, *(str(row[field]) for row in condition_rows)]
                worksheet.set_column(column_number, column_number, min(max(map(len, values)) + 2, 40))
            if condition_rows:
                worksheet.add_table(
                    0,
                    0,
                    len(condition_rows),
                    len(fields) - 1,
                    {
                        "columns": [{"header": field} for field in fields],
                        "style": "Table Style Medium 2",
                    },
                )
            else:
                worksheet.autofilter(0, 0, 0, len(fields) - 1)


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comet_calls = read_comet_calls(args.comet_csv)
    summary_rows: list[dict[str, int | str]] = []
    difference_rows: list[dict[str, str]] = []
    for gene in GENES:
        calls_by_method = {
            "Comet": comet_calls,
            "PerGene": read_blast_calls(
                args.pergene_dir / f"{gene}_AllSeq_NonComet_Coverage.csv"
            ),
            "FullGenome": read_blast_calls(
                args.fullgenome_dir / f"{gene}_AllSeq_NonComet_Coverage.csv"
            ),
        }
        for condition, methods in COMPARISONS.items():
            summary, differences = compare_calls(gene, condition, methods, calls_by_method)
            summary_rows.append(summary)
            difference_rows.extend(differences)

    write_csv(
        args.output_dir / "Subtyping_Comparison_summary.csv",
        ["Gene", "SameCount", "DifferentCount", "CompareCondition"],
        summary_rows,
    )
    write_csv(
        args.output_dir / "Subtyping_Comparison_differences.csv",
        [
            "Gene",
            "Accession",
            "CometSubtype",
            "BlastPerGeneSubtype",
            "BlastFullGenomeSubtype",
            "CompareCondition",
        ],
        difference_rows,
    )
    subtype_counts = Counter(
        (
            row["Gene"],
            row["CometSubtype"],
            row["BlastPerGeneSubtype"],
            row["BlastFullGenomeSubtype"],
            row["CompareCondition"],
        )
        for row in difference_rows
    )
    count_rows = [
        {
            "Gene": gene,
            "CometSubtype": comet_subtype,
            "BlastPerGeneSubtype": pergene_subtype,
            "BlastFullGenomeSubtype": fullgenome_subtype,
            "DifferentCount": count,
            "CompareCondition": condition,
        }
        for (gene, comet_subtype, pergene_subtype, fullgenome_subtype, condition), count in sorted(
            subtype_counts.items(), key=lambda item: (item[0][4], item[0][0], -item[1], item[0])
        )
    ]
    count_fields = [
        "Gene",
        "CometSubtype",
        "BlastPerGeneSubtype",
        "BlastFullGenomeSubtype",
        "DifferentCount",
        "CompareCondition",
    ]
    write_csv(
        args.output_dir / "Subtyping_Comparison_difference_subtype_counts.csv",
        count_fields,
        count_rows,
    )
    write_count_workbook(
        args.output_dir / "Subtyping_Comparison_difference_subtype_counts.xlsx",
        count_fields,
        count_rows,
    )
    for row in summary_rows:
        print(
            f"{row['Gene']} {row['CompareCondition']}: "
            f"same={row['SameCount']}, different={row['DifferentCount']}"
        )


if __name__ == "__main__":
    main()
