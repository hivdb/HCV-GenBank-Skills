#!/usr/bin/env python3
"""Compare Comet and local-alignment profile-input assignments."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--comet-csv", required=True)
    parser.add_argument("--local-csv", required=True)
    parser.add_argument("--summary-csv", required=True)
    parser.add_argument("--differences-csv", required=True)
    return parser.parse_args()


def load_rows(path: Path) -> dict[str, dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"accession", "genotype", "subtype"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise RuntimeError(f"Columns missing from {path}: {', '.join(sorted(missing))}")
        return {
            row["accession"].strip(): {
                "genotype": row["genotype"].strip(),
                "subtype": row["subtype"].strip(),
            }
            for row in reader
            if row.get("accession", "").strip()
        }


def write_summary(path: Path, metrics: list[tuple[str, str | int | float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerows(metrics)


def main() -> int:
    args = parse_args()
    comet_path = Path(args.comet_csv)
    local_path = Path(args.local_csv)
    summary_path = Path(args.summary_csv)
    differences_path = Path(args.differences_csv)

    if not local_path.is_file():
        metrics: list[tuple[str, str | int | float]] = [
            ("status", "local_profile_input_csv_missing"),
            ("local_profile_input_csv", str(local_path)),
        ]
        write_summary(summary_path, metrics)
        differences_path.parent.mkdir(parents=True, exist_ok=True)
        with differences_path.open("w", encoding="utf-8", newline="") as handle:
            csv.writer(handle).writerow(["accession", "comet_genotype", "local_genotype", "comet_subtype", "local_subtype"])
        print("profile_assignment_comparison_status=local_profile_input_csv_missing")
        print(f"profile_assignment_comparison_summary={summary_path}")
        return 0

    comet = load_rows(comet_path)
    local = load_rows(local_path)
    shared = sorted(set(comet).intersection(local))
    same_gt = diff_gt = same_subtype = diff_subtype = 0
    differences: list[list[str]] = []
    for accession in shared:
        comet_row = comet[accession]
        local_row = local[accession]
        genotype_same = comet_row["genotype"] == local_row["genotype"]
        subtype_same = comet_row["subtype"] == local_row["subtype"]
        same_gt += genotype_same
        diff_gt += not genotype_same
        same_subtype += subtype_same
        diff_subtype += not subtype_same
        if not genotype_same or not subtype_same:
            differences.append([accession, comet_row["genotype"], local_row["genotype"], comet_row["subtype"], local_row["subtype"]])

    shared_count = len(shared)
    genotype_diff_pct = 100.0 * diff_gt / shared_count if shared_count else 0.0
    subtype_diff_pct = 100.0 * diff_subtype / shared_count if shared_count else 0.0
    metrics = [
        ("status", "ok"),
        ("comet_profile_accession_count", len(comet)),
        ("local_profile_accession_count", len(local)),
        ("shared_accession_count", shared_count),
        ("comet_only_accession_count", len(set(comet).difference(local))),
        ("local_only_accession_count", len(set(local).difference(comet))),
        ("same_genotype_count", same_gt),
        ("different_genotype_count", diff_gt),
        ("different_genotype_percent", f"{genotype_diff_pct:.2f}"),
        ("same_subtype_count", same_subtype),
        ("different_subtype_count", diff_subtype),
        ("different_subtype_percent", f"{subtype_diff_pct:.2f}"),
    ]
    write_summary(summary_path, metrics)
    differences_path.parent.mkdir(parents=True, exist_ok=True)
    with differences_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "comet_genotype", "local_genotype", "comet_subtype", "local_subtype"])
        writer.writerows(differences)

    for metric, value in metrics[6:]:
        print(f"profile_assignment_{metric}={value}")
    print(f"profile_assignment_comparison_summary={summary_path}")
    print(f"profile_assignment_differences={differences_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
