#!/usr/bin/env python3
"""Filter per-gene non-COMET calls and summarize agreement with COMET."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from pathlib import Path


GENES = ("NS3", "NS5A", "NS5B")
REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_NONCOMET_DIR = REPO_ROOT / "HCVData" / "nonComet-PerGene"
DEFAULT_COMET_CSV = REPO_ROOT / "HCVData" / "HCV-all-seq-subtype" / "all_comet_subtype.csv"
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


def compare_gene(
    gene: str, noncomet_path: Path, comet_calls: dict[str, str]
) -> tuple[dict[str, int | str], list[dict[str, str]]]:
    with noncomet_path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Accession", "ClosestGenotype", "ClosestSubtype"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{noncomet_path} is missing columns: {', '.join(sorted(missing))}"
            )
        rows_with_genotype: list[dict[str, str]] = []
        for row in reader:
            if normalized(row.get("ClosestGenotype")):
                rows_with_genotype.append(row)

    same_count = 0
    differences: list[dict[str, str]] = []
    for row in rows_with_genotype:
        comet_subtype = comet_calls.get(accession_key(row.get("Accession")), "")
        noncomet_subtype = normalized(row.get("ClosestSubtype"))
        if comet_subtype == noncomet_subtype:
            same_count += 1
        else:
            differences.append(
                {
                    "Gene": gene,
                    "Accession": str(row.get("Accession") or "").strip(),
                    "CometSubtype": comet_subtype,
                    "BlastPerGeneSubtype": noncomet_subtype,
                    "CompareCondition": "Comet_vs_PerGene",
                }
            )
    return (
        {
            "Gene": gene,
            "SameCount": same_count,
            "DifferentCount": len(differences),
            "CompareCondition": "Comet_vs_PerGene",
        },
        differences,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--noncomet-dir", type=Path, default=DEFAULT_NONCOMET_DIR)
    parser.add_argument("--comet-csv", type=Path, default=DEFAULT_COMET_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comet_calls = read_comet_calls(args.comet_csv)
    summary_rows = []
    difference_rows = []
    for gene in GENES:
        summary, differences = compare_gene(
            gene,
            args.noncomet_dir / f"{gene}_AllSeq_NonComet_Coverage.csv",
            comet_calls,
        )
        summary_rows.append(summary)
        difference_rows.extend(differences)
    report_path = args.output_dir / "Subtyping_Comparison_summary.csv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["Gene", "SameCount", "DifferentCount", "CompareCondition"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    differences_path = args.output_dir / "Subtyping_Comparison_differences.csv"
    with differences_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Gene",
                "Accession",
                "CometSubtype",
                "BlastPerGeneSubtype",
                "CompareCondition",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(difference_rows)
    subtype_counts = Counter(
        (
            row["Gene"],
            row["CometSubtype"],
            row["BlastPerGeneSubtype"],
        )
        for row in difference_rows
    )
    subtype_counts_path = args.output_dir / "Subtyping_Comparison_difference_subtype_counts.csv"
    with subtype_counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "Gene",
                "CometSubtype",
                "BlastPerGeneSubtype",
                "DifferentCount",
                "CompareCondition",
            ],
            lineterminator="\n",
        )
        writer.writeheader()
        for (gene, comet_subtype, blast_subtype), count in sorted(
            subtype_counts.items(), key=lambda item: (-item[1], item[0])
        ):
            writer.writerow(
                {
                    "Gene": gene,
                    "CometSubtype": comet_subtype,
                    "BlastPerGeneSubtype": blast_subtype,
                    "DifferentCount": count,
                    "CompareCondition": "Comet_vs_PerGene",
                }
            )
    for row in summary_rows:
        print(
            f"{row['Gene']}: same={row['SameCount']}, "
            f"different={row['DifferentCount']}"
        )
    print(report_path)
    print(differences_path)
    print(subtype_counts_path)


if __name__ == "__main__":
    main()
