#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


RULES_CSV = Path(__file__).with_name("refid_metadata_filters.csv")
REPO_ROOT = Path(__file__).resolve().parents[3]
ACCESSION_LIST_DIR = REPO_ROOT / "HCVData" / "Ref-selection" / "NS5_Ref_filter" / "NS5A"


@dataclass(frozen=True)
class FilterRule:
    refid: str
    field: str
    operator: str
    value: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split included_accessions_metadata.csv into per-RefID CSV files for "
            "NS5A-specific filters."
        )
    )
    parser.add_argument("--input-csv", required=True, help="Path to included_accessions_metadata.csv")
    parser.add_argument("--output-dir", required=True, help="Directory for per-RefID CSV files")
    parser.add_argument(
        "--rules-csv",
        default=RULES_CSV,
        help="CSV containing RefID, Field, Operator, and Value columns",
    )
    parser.add_argument(
        "--accession-list-dir",
        default=ACCESSION_LIST_DIR,
        help="Directory containing accession-list CSVs referenced by in_accession_list rules",
    )
    return parser.parse_args()


def load_rows(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "RefID" not in fieldnames:
            raise RuntimeError(f"Column 'RefID' was not found in {path}")
        return fieldnames, list(reader)


def load_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        accession_column = "Accession" if "Accession" in fieldnames else "IsolateID"
        if accession_column not in fieldnames:
            raise RuntimeError(f"Column 'Accession' or 'IsolateID' was not found in {path}")
        return {
            accession
            for row in reader
            if (accession := (row.get(accession_column) or "").strip())
        }


def load_rules(path: Path) -> list[FilterRule]:
    if not path.is_file():
        raise RuntimeError(f"Rules CSV was not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"RefID", "Field", "Operator", "Value"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise RuntimeError(f"Rules CSV must contain {', '.join(sorted(expected))}: {path}")
        rules = [
            FilterRule(
                refid=(row["RefID"] or "").strip(),
                field=(row["Field"] or "").strip(),
                operator=(row["Operator"] or "").strip(),
                value=(row["Value"] or "").strip(),
            )
            for row in reader
        ]
    if not rules or any(not all((rule.refid, rule.field, rule.operator, rule.value)) for rule in rules):
        raise RuntimeError(f"Rules CSV contains an empty required value: {path}")
    supported = {"contains", "not_contains", "equals", "in_accession_list"}
    invalid = sorted({rule.operator for rule in rules} - supported)
    if invalid:
        raise RuntimeError(f"Unsupported rule operator(s): {', '.join(invalid)}")
    if len({rule.refid for rule in rules}) != len(rules):
        raise RuntimeError(f"Rules CSV must contain one rule per RefID: {path}")
    return rules


def rule_description(rule: FilterRule) -> str:
    if rule.operator == "in_accession_list":
        return f"{rule.field} in {rule.value}"
    if rule.operator == "not_contains":
        return f"{rule.field} does not contain {rule.value}"
    if rule.operator == "equals":
        return f"{rule.field} equals {rule.value}"
    return f"{rule.field} contains {rule.value}"


def row_is_kept(row: dict[str, str], rule: FilterRule, accession_list: set[str] | None) -> bool:
    cell = (row.get(rule.field) or "").strip()
    if rule.operator == "in_accession_list":
        return cell in (accession_list or set())
    if rule.operator == "contains":
        return rule.value.casefold() in cell.casefold()
    if rule.operator == "not_contains":
        return rule.value.casefold() not in cell.casefold()
    return cell.casefold() == rule.value.casefold()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_csv = Path(args.input_csv).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    rules_csv = Path(args.rules_csv).expanduser()
    accession_list_dir = Path(args.accession_list_dir).expanduser()

    if not input_csv.is_file():
        raise RuntimeError(f"Input CSV was not found: {input_csv}")
    output_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = load_rows(input_csv)
    rules = load_rules(rules_csv)
    unknown_fields = sorted({rule.field for rule in rules} - set(fieldnames))
    if unknown_fields:
        raise RuntimeError(f"Rule field(s) not found in {input_csv}: {', '.join(unknown_fields)}")
    input_accessions = {
        (row.get("Accession") or "").strip()
        for row in rows
        if (row.get("Accession") or "").strip()
    }
    accession_filters = {
        rule.refid: load_accessions(accession_list_dir / rule.value)
        for rule in rules
        if rule.operator == "in_accession_list"
    }
    rows_by_refid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        refid = (row.get("RefID") or "").strip()
        if refid:
            rows_by_refid[refid].append(row)

    summary_rows: list[dict[str, str | int]] = []
    output_accessions: set[str] = set()
    for rule in sorted(rules, key=lambda item: (int(item.refid), item.refid)):
        ref_rows = rows_by_refid[rule.refid]
        kept_rows = [row for row in ref_rows if row_is_kept(row, rule, accession_filters.get(rule.refid))]
        output_accessions.update(
            (row.get("Accession") or "").strip()
            for row in kept_rows
            if (row.get("Accession") or "").strip()
        )
        output_csv = output_dir / f"RefID_{rule.refid}_metadata.csv"
        write_csv(output_csv, fieldnames, kept_rows)
        summary_rows.append(
            {
                "RefID": rule.refid,
                "Filter": rule_description(rule),
                "TotalRows": len(ref_rows),
                "KeptRows": len(kept_rows),
                "RemovedRows": len(ref_rows) - len(kept_rows),
            }
        )

    print(f"refid_count={len(summary_rows)}")
    print(f"input_accession_count={len(input_accessions)}")
    print(f"output_accession_count={len(output_accessions)}")
    print(f"input_row_count={len(rows)}")
    print(f"filtered_refids={','.join(rule.refid for rule in rules)}")
    print(f"output_dir={output_dir.resolve()}")
    for row in summary_rows:
        print(
            "filter_result="
            f"RefID:{row['RefID']},"
            f"Filter:{row['Filter']},"
            f"TotalRows:{row['TotalRows']},"
            f"KeptRows:{row['KeptRows']},"
            f"RemovedRows:{row['RemovedRows']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
