#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
ACCESSION_LIST_DIR = REPO_ROOT / "HCVData" / "Ref-selection" / "NS5_Ref_filter" / "NS3"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Split included_accessions_metadata.csv into one CSV per RefID, applying "
            "RefID-specific row filters and writing a filter summary."
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
        "--accession-list-dir",
        default=ACCESSION_LIST_DIR,
        help="Directory containing manual accession-list CSVs",
    )
    return parser.parse_args()


def sanitize_filename(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return text.strip("._-") or "unknown"


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


def refid_filter_description(refid: str) -> str:
    if refid == "30":
        return "source_isolate contains Day1"
    if refid == "142":
        return "source_isolate contains baseline"
    if refid == "192":
        return "source_isolate contains day1"
    if refid == "346":
        return "source_isolate contains baseline/D0"
    if refid == "600":
        return "source_isolate does not contain failure"
    if refid == "661":
        return "source_isolation_source == plasma"
    if refid == "884":
        return "source_isolate contains Pre-TH"
    if refid == "943":
        return "source_isolate contains Day 1"
    if refid == "1356":
        return "source_isolate does not contain IC"
    if refid == "2008":
        return "source_isolate does not contain chimpanzee"
    if refid == "2110":
        return "source_isolate contains T0"
    if refid == "2116":
        return "source_collection_date before 2011"
    if refid == "2150":
        return "source_isolate contains b"
    if refid == "2168":
        return "source_isolate contains pre"
    if refid == "2178":
        return "source_isolation_source == plasma"
    raise KeyError(refid)


def text_contains(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() in (row.get(column) or "").casefold()


def text_does_not_contain(row: dict[str, str], column: str, needle: str) -> bool:
    return needle.casefold() not in (row.get(column) or "").casefold()


def text_equals(row: dict[str, str], column: str, value: str) -> bool:
    return (row.get(column) or "").strip().casefold() == value.casefold()


def first_year(value: str) -> int | None:
    match = re.search(r"\b(19|20)\d{2}\b", value)
    return int(match.group(0)) if match else None


def row_is_kept(refid: str, row: dict[str, str]) -> bool:
    if refid == "30":
        return text_contains(row, "source_isolate", "Day1")
    if refid == "142":
        return text_contains(row, "source_isolate", "baseline")
    if refid == "192":
        return text_contains(row, "source_isolate", "day1")
    if refid == "346":
        return text_contains(row, "source_isolate", "baseline/D0")
    if refid == "600":
        return text_does_not_contain(row, "source_isolate", "failure")
    if refid == "661":
        return text_equals(row, "source_isolation_source", "plasma")
    if refid == "884":
        return text_contains(row, "source_isolate", "Pre-TH")
    if refid == "943":
        return text_contains(row, "source_isolate", "Day 1")
    if refid == "1356":
        return text_does_not_contain(row, "source_isolate", "IC")
    if refid == "2008":
        return text_does_not_contain(row, "source_isolate", "chimpanzee")
    if refid == "2110":
        return text_contains(row, "source_isolate", "T0")
    if refid == "2116":
        year = first_year(row.get("source_collection_date") or "")
        return year is not None and year < 2011
    if refid == "2150":
        return text_contains(row, "source_isolate", "b")
    if refid == "2168":
        return text_contains(row, "source_isolate", "pre")
    if refid == "2178":
        return text_equals(row, "source_isolation_source", "plasma")
    raise KeyError(refid)


def filtered_refids() -> set[str]:
    return {
        "30",
        "142",
        "192",
        "346",
        "600",
        "661",
        "884",
        "943",
        "1356",
        "2008",
        "2110",
        "2116",
        "2150",
        "2168",
        "2178",
    }


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_csv = Path(args.input_csv).expanduser()
    output_dir = Path(args.output_dir).expanduser()
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
    input_accessions = {
        (row.get("Accession") or "").strip()
        for row in rows
        if (row.get("Accession") or "").strip()
    }
    accession_filters = {
        "85": load_accessions(accession_list_dir / "85.csv"),
    }
    rows_by_refid: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        refid = (row.get("RefID") or "").strip()
        if refid:
            rows_by_refid[refid].append(row)

    summary_rows: list[dict[str, str | int | None]] = []
    output_accessions: set[str] = set()
    for refid in sorted(
        filtered_refids(),
        key=lambda value: (int(value) if value.isdigit() else 10**12, value),
    ):
        ref_rows = rows_by_refid[refid]
        source_fasta_count = fasta_accession_count(source_fasta_dir, refid)
        comet_fasta_count = fasta_accession_count(comet_fasta_dir, refid)
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
        if refid in accession_filters:
            listed_accessions = accession_filters[refid]
            input_refid_accessions = {
                (row.get("Accession") or "").strip()
                for row in ref_rows
                if (row.get("Accession") or "").strip()
            }
            kept_rows = [
                row
                for row in ref_rows
                if (row.get("Accession") or "").strip() in listed_accessions
            ]
            listed_accession_count = len(listed_accessions)
            listed_present_in_input_count = len(
                listed_accessions & input_refid_accessions
            )
            listed_absent_from_input_count = len(
                listed_accessions - input_refid_accessions
            )
        else:
            kept_rows = [row for row in ref_rows if row_is_kept(refid, row)]
        output_accessions.update(
            (row.get("Accession") or "").strip()
            for row in kept_rows
            if (row.get("Accession") or "").strip()
        )
        output_csv = output_dir / f"RefID_{sanitize_filename(refid)}_metadata.csv"
        write_csv(output_csv, fieldnames, kept_rows)
        summary_rows.append(
            {
                "RefID": refid,
                "Filter": refid_filter_description(refid),
                "TotalRows": total_rows,
                "CometRetainedRows": comet_fasta_count,
                "CometExcludedRows": comet_excluded_count,
                "KeptRows": len(kept_rows),
                "RemovedRows": len(ref_rows) - len(kept_rows),
                "ListedAccessions": listed_accession_count,
                "ListedPresentInInput": listed_present_in_input_count,
                "ListedAbsentFromInput": listed_absent_from_input_count,
                "OutputCSV": str(output_csv.resolve()),
            }
        )

    print(f"refid_count={len(summary_rows)}")
    print(f"input_accession_count={len(input_accessions)}")
    print(f"output_accession_count={len(output_accessions)}")
    print(f"input_row_count={len(rows)}")
    print(f"filtered_refids={','.join(sorted(filtered_refids(), key=int))}")
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
