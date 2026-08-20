#!/usr/bin/env python3
"""Select non-COMET assignments that take priority over COMET calls."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--noncomet-coverage-csv", required=True)
    parser.add_argument("--output-csv", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    with Path(args.noncomet_coverage_csv).open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        fieldnames = reader.fieldnames
        if not fieldnames:
            raise RuntimeError(f"No header found in {args.noncomet_coverage_csv}")
        rows = [
            row for row in reader
            if (row.get("Accession") or "").strip()
            and (
                (row.get("ClosestSubtype") or "").strip().lower() == "1d"
                or (row.get("ClosestGenotype") or "").strip().lower() in {"7", "8"}
                or (row.get("ClosestSubtype") or "").strip().lower().startswith(("7", "8"))
            )
        ]
    output = Path(args.output_csv)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps({"output_csv": str(output.resolve()), "priority_assignment_count": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
