#!/usr/bin/env python3
"""Compare all-sequence non-COMET coverage calls with folder BLAST assignments."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
GENES = ("NS3", "NS5A", "NS5B")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--coverage-dir", type=Path, default=REPO_ROOT / "HCVData/HCV-all-seq-subtype"
    )
    parser.add_argument(
        "--assignments-dir",
        type=Path,
        default=REPO_ROOT / "archived-skills/outputs/folder_assignments",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT
        / "archived-skills/outputs/noncomet-coverage-assignment-comparison",
    )
    return parser.parse_args()


def accession_key(value: object) -> str:
    return str(value or "").strip().split(".", 1)[0].upper()


def normalized_call(value: object) -> str:
    return str(value or "").strip().lower().removeprefix("gt")


def read_calls(
    path: Path, accession_column: str, genotype_column: str, subtype_column: str
) -> tuple[dict[str, dict[str, str]], int]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {accession_column, genotype_column, subtype_column}
        missing = sorted(required - set(reader.fieldnames or []))
        if missing:
            raise RuntimeError(
                f"{path} is missing required columns: {', '.join(missing)}"
            )
        rows: dict[str, dict[str, str]] = {}
        input_rows = 0
        for row in reader:
            input_rows += 1
            key = accession_key(row.get(accession_column))
            if not key:
                continue
            call = {
                "Accession": str(row.get(accession_column) or "").strip(),
                "Genotype": str(row.get(genotype_column) or "").strip(),
                "Subtype": str(row.get(subtype_column) or "").strip(),
            }
            if key in rows and rows[key] != call:
                raise RuntimeError(f"{path} has conflicting rows for accession {key}")
            rows[key] = call
    return rows, input_rows


def write_csv(path: Path, fields: list[str], rows: list[dict[str, str | int]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_rows: list[dict[str, str | int]] = []

    for gene in GENES:
        coverage_path = args.coverage_dir / f"{gene}_AllSeq_NonComet_Coverage.csv"
        assignment_path = args.assignments_dir / f"{gene}_assignments.csv"
        if not coverage_path.is_file() or not assignment_path.is_file():
            missing = [
                str(path)
                for path in (coverage_path, assignment_path)
                if not path.is_file()
            ]
            raise SystemExit(f"Missing required input file(s): {', '.join(missing)}")

        coverage, coverage_rows = read_calls(
            coverage_path, "Accession", "ClosestGenotype", "ClosestSubtype"
        )
        assignments, assignment_rows = read_calls(
            assignment_path, "accession", "genotype", "subtype"
        )
        shared = sorted(set(coverage) & set(assignments))
        only_in_coverage = sorted(set(coverage) - set(assignments))
        only_in_assignments = sorted(set(assignments) - set(coverage))
        differences: list[dict[str, str | int]] = []
        genotype_difference_count = 0
        subtype_difference_count = 0

        for key in shared:
            coverage_call = coverage[key]
            assignment_call = assignments[key]
            genotype_differs = normalized_call(
                coverage_call["Genotype"]
            ) != normalized_call(assignment_call["Genotype"])
            subtype_differs = normalized_call(
                coverage_call["Subtype"]
            ) != normalized_call(assignment_call["Subtype"])
            if not genotype_differs and not subtype_differs:
                continue
            genotype_difference_count += int(genotype_differs)
            subtype_difference_count += int(subtype_differs)
            differences.append(
                {
                    "Accession": key,
                    "FromDbAccession": coverage_call["Accession"],
                    "FromFastaAccession": assignment_call["Accession"],
                    "FromDbGenotype": coverage_call["Genotype"],
                    "FromFastaGenotype": assignment_call["Genotype"],
                    "FromDbSubtype": coverage_call["Subtype"],
                    "FromFastaSubtype": assignment_call["Subtype"],
                    "DifferenceType": "BOTH"
                    if genotype_differs and subtype_differs
                    else "GENOTYPE"
                    if genotype_differs
                    else "SUBTYPE",
                }
            )

        write_csv(
            args.output_dir / f"{gene}_genotype_subtype_differences.csv",
            [
                "Accession",
                "FromDbAccession",
                "FromFastaAccession",
                "FromDbGenotype",
                "FromFastaGenotype",
                "FromDbSubtype",
                "FromFastaSubtype",
                "DifferenceType",
            ],
            differences,
        )
        write_csv(
            args.output_dir / f"{gene}_only_in_from_db.csv",
            ["Accession", "FromDbAccession", "FromDbGenotype", "FromDbSubtype"],
            [
                {
                    "Accession": key,
                    "FromDbAccession": coverage[key]["Accession"],
                    "FromDbGenotype": coverage[key]["Genotype"],
                    "FromDbSubtype": coverage[key]["Subtype"],
                }
                for key in only_in_coverage
            ],
        )
        write_csv(
            args.output_dir / f"{gene}_only_in_from_fasta.csv",
            [
                "Accession",
                "FromFastaAccession",
                "FromFastaGenotype",
                "FromFastaSubtype",
            ],
            [
                {
                    "Accession": key,
                    "FromFastaAccession": assignments[key]["Accession"],
                    "FromFastaGenotype": assignments[key]["Genotype"],
                    "FromFastaSubtype": assignments[key]["Subtype"],
                }
                for key in only_in_assignments
            ],
        )
        summary_rows.append(
            {
                "Gene": gene,
                "FromDbInputRows": coverage_rows,
                "FromDbUniqueAccessions": len(coverage),
                "FromFastaInputRows": assignment_rows,
                "FromFastaUniqueAccessions": len(assignments),
                "SharedAccessions": len(shared),
                "OnlyInFromDb": len(only_in_coverage),
                "OnlyInFromFasta": len(only_in_assignments),
                "GenotypeDifferences": genotype_difference_count,
                "SubtypeDifferences": subtype_difference_count,
                "AnyGenotypeOrSubtypeDifference": len(differences),
            }
        )

    summary_fields = [
        "Gene",
        "FromDbInputRows",
        "FromDbUniqueAccessions",
        "FromFastaInputRows",
        "FromFastaUniqueAccessions",
        "SharedAccessions",
        "OnlyInFromDb",
        "OnlyInFromFasta",
        "GenotypeDifferences",
        "SubtypeDifferences",
        "AnyGenotypeOrSubtypeDifference",
    ]
    write_csv(
        args.output_dir / "accession_count_summary.csv", summary_fields, summary_rows
    )
    for row in summary_rows:
        print(
            f"{row['Gene']}: from_db={row['FromDbUniqueAccessions']} from_fasta={row['FromFastaUniqueAccessions']} "
            f"shared={row['SharedAccessions']} only_from_db={row['OnlyInFromDb']} "
            f"only_from_fasta={row['OnlyInFromFasta']} differences={row['AnyGenotypeOrSubtypeDifference']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
