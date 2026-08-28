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
    parser.add_argument(
        "--input-csv", required=True, help="Path to included_accessions_metadata.csv"
    )
    parser.add_argument(
        "--output-dir", required=True, help="Directory for per-RefID CSV files"
    )
    parser.add_argument(
        "--source-fasta-dir",
        help="Directory containing RefID FASTAs before COMET filtering",
    )
    parser.add_argument(
        "--comet-fasta-dir",
        help="Directory containing RefID FASTAs after COMET filtering",
    )
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
            raise RuntimeError(
                f"Column 'Accession' or 'IsolateID' was not found in {path}"
            )
        return {
            accession
            for row in reader
            if (accession := (row.get(accession_column) or "").strip())
        }


def fasta_accession_count(fasta_dir: Path | None, refid: str) -> int | None:
    if fasta_dir is None:
        return None
    fasta_paths = list(fasta_dir.glob(f"{refid}_*.fasta"))
    if not fasta_paths:
        return None
    if len(fasta_paths) != 1:
        raise RuntimeError(
            f"Expected one RefID {refid} FASTA under {fasta_dir}, found {len(fasta_paths)}"
        )
    with fasta_paths[0].open(encoding="utf-8") as handle:
        return sum(1 for line in handle if line.startswith(">"))


def load_rules(path: Path) -> list[FilterRule]:
    if not path.is_file():
        raise RuntimeError(f"Rules CSV was not found: {path}")
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        expected = {"RefID", "Field", "Operator", "Value"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise RuntimeError(
                f"Rules CSV must contain {', '.join(sorted(expected))}: {path}"
            )
        rules = [
            FilterRule(
                refid=(row["RefID"] or "").strip(),
                field=(row["Field"] or "").strip(),
                operator=(row["Operator"] or "").strip(),
                value=(row["Value"] or "").strip(),
            )
            for row in reader
        ]
    if not rules or any(
        not all((rule.refid, rule.field, rule.operator, rule.value)) for rule in rules
    ):
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


def row_is_kept(
    row: dict[str, str], rule: FilterRule, accession_list: set[str] | None
) -> bool:
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
    source_fasta_dir = (
        Path(args.source_fasta_dir).expanduser() if args.source_fasta_dir else None
    )
    comet_fasta_dir = (
        Path(args.comet_fasta_dir).expanduser() if args.comet_fasta_dir else None
    )

    if not input_csv.is_file():
        raise RuntimeError(f"Input CSV was not found: {input_csv}")
    for label, fasta_dir in (("source", source_fasta_dir), ("COMET", comet_fasta_dir)):
        if fasta_dir is not None and not fasta_dir.is_dir():
            raise RuntimeError(f"{label} FASTA directory was not found: {fasta_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)

    fieldnames, rows = load_rows(input_csv)
    rules = load_rules(rules_csv)
    unknown_fields = sorted({rule.field for rule in rules} - set(fieldnames))
    if unknown_fields:
        raise RuntimeError(
            f"Rule field(s) not found in {input_csv}: {', '.join(unknown_fields)}"
        )
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

    summary_rows: list[dict[str, str | int | None]] = []
    output_accessions: set[str] = set()
    for rule in sorted(rules, key=lambda item: (int(item.refid), item.refid)):
        ref_rows = rows_by_refid[rule.refid]
        source_fasta_count = fasta_accession_count(source_fasta_dir, rule.refid)
        comet_fasta_count = fasta_accession_count(comet_fasta_dir, rule.refid)
        total_rows = (
            source_fasta_count if source_fasta_count is not None else len(ref_rows)
        )
        comet_excluded_count = (
            source_fasta_count - comet_fasta_count
            if source_fasta_count is not None and comet_fasta_count is not None
            else None
        )
        listed_accession_count: int | None = None
        listed_present_in_input_count: int | None = None
        listed_absent_from_input_count: int | None = None
        if rule.refid in accession_filters:
            listed_accessions = accession_filters[rule.refid]
            input_refid_accessions = {
                (row.get("Accession") or "").strip()
                for row in ref_rows
                if (row.get("Accession") or "").strip()
            }
            listed_accession_count = len(listed_accessions)
            listed_present_in_input_count = len(
                listed_accessions & input_refid_accessions
            )
            listed_absent_from_input_count = len(
                listed_accessions - input_refid_accessions
            )
        kept_rows = [
            row
            for row in ref_rows
            if row_is_kept(row, rule, accession_filters.get(rule.refid))
        ]
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
                "TotalRows": total_rows,
                "CometRetainedRows": comet_fasta_count,
                "CometExcludedRows": comet_excluded_count,
                "KeptRows": len(kept_rows),
                "RemovedRows": len(ref_rows) - len(kept_rows),
                "ListedAccessions": listed_accession_count,
                "ListedPresentInInput": listed_present_in_input_count,
                "ListedAbsentFromInput": listed_absent_from_input_count,
            }
        )

    print(f"refid_count={len(summary_rows)}")
    print(f"input_accession_count={len(input_accessions)}")
    print(f"output_accession_count={len(output_accessions)}")
    print(f"input_row_count={len(rows)}")
    print(
        f"filtered_refids={','.join(rule.refid for rule in sorted(rules, key=lambda rule: int(rule.refid)))}"
    )
    print(f"filtered_refid_count={len(summary_rows)}")
    print(f"output_dir={output_dir.resolve()}")
    for row in summary_rows:
        print("filter_result:")
        print(f"  RefID: {row['RefID']}")
        print(f"  Source FASTA rows: {row['TotalRows']}")
        if row["CometRetainedRows"] is not None:
            print(f"  After COMET: {row['CometRetainedRows']}")
            print(f"  Removed by COMET: {row['CometExcludedRows']}")
        print(f"  RefID filter: {row['Filter']}")
        if row["ListedAccessions"] is not None:
            print(f"  Listed accessions: {row['ListedAccessions']}")
            print(
                f"  Listed accessions present after COMET: {row['ListedPresentInInput']}"
            )
            print(
                f"  Listed accessions absent after COMET: {row['ListedAbsentFromInput']}"
            )
        print(f"  After RefID filter: {row['KeptRows']}")
        print(f"  Removed by RefID filter: {row['RemovedRows']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
