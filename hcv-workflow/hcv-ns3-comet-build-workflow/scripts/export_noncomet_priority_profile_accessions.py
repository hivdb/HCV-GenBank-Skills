#!/usr/bin/env python3
"""Report profile-used accessions whose subtype is prioritized from non-Comet."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def accession_key(value: object) -> str:
    return str(value or "").strip().split(".", 1)[0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile-accessions-csv", required=True)
    parser.add_argument("--comet-subtype-csv", required=True)
    parser.add_argument("--noncomet-coverage-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with Path(args.profile_accessions_csv).open(encoding="utf-8", newline="") as handle:
        profile_accessions = {
            accession_key(row["accession"]): row["accession"].strip()
            for row in csv.DictReader(handle)
            if accession_key(row.get("accession"))
        }
    with Path(args.comet_subtype_csv).open(encoding="utf-8", newline="") as handle:
        comet_subtypes = {
            accession_key(row.get("name") or row.get("accession")): str(
                row.get("subtype") or ""
            )
            .strip()
            .lower()
            for row in csv.DictReader(handle)
        }

    priority_subtypes: dict[str, str] = {}
    with Path(args.noncomet_coverage_csv).open(
        encoding="utf-8-sig", newline=""
    ) as handle:
        for row in csv.DictReader(handle):
            accession = accession_key(row.get("Accession"))
            genotype = str(row.get("ClosestGenotype") or "").strip().lower()
            subtype = str(row.get("ClosestSubtype") or "").strip().lower()
            if accession and genotype in {"7", "8"}:
                priority_subtypes[accession] = subtype

    rows = [
        (profile_accessions[key], comet_subtypes.get(key, ""), subtype)
        for key, subtype in priority_subtypes.items()
        if key in profile_accessions
    ]
    rows.sort(key=lambda row: row[0])
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["AccessionID", "CometSubtype", "NonCometSubtype"])
        writer.writerows(rows)
    print(
        json.dumps(
            {
                "output_csv": str(output.resolve()),
                "profile_used_noncomet_priority_accession_count": len(rows),
            }
        )
    )


if __name__ == "__main__":
    main()
