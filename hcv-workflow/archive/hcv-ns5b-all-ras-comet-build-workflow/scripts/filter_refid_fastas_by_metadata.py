#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Filter copied RefID FASTA files in place using per-RefID metadata CSV "
            "files produced by split_refid_metadata_csv.py."
        )
    )
    parser.add_argument(
        "--metadata-dir",
        required=True,
        help="Directory containing RefID_*_metadata.csv files",
    )
    parser.add_argument(
        "--fasta-dir",
        required=True,
        help="Directory containing copied RefID FASTA files",
    )
    parser.add_argument(
        "--kept-accessions-output",
        required=True,
        help="CSV manifest of accessions retained after RefID metadata filtering.",
    )
    return parser.parse_args()


def refid_from_metadata_path(path: Path) -> str | None:
    match = re.fullmatch(r"RefID_(.+)_metadata\.csv", path.name)
    return match.group(1) if match else None


def refid_from_fasta_path(path: Path) -> str | None:
    match = re.match(r"([^_]+)_", path.name)
    return match.group(1) if match else None


def refid_sort_key(refid: str) -> tuple[int, int | str]:
    return (0, int(refid)) if refid.isdigit() else (1, refid)


def load_metadata_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "Accession" not in fieldnames:
            raise RuntimeError(f"Column 'Accession' was not found in {path}")
        return {
            accession
            for row in reader
            if (accession := (row.get("Accession") or "").strip())
        }


def read_fasta(path: Path) -> list[tuple[str, list[str]]]:
    records: list[tuple[str, list[str]]] = []
    header: str | None = None
    sequence_lines: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    records.append((header, sequence_lines))
                header = line.rstrip("\n")
                sequence_lines = []
            else:
                sequence_lines.append(line.rstrip("\n"))
    if header is not None:
        records.append((header, sequence_lines))
    return records


def header_accession(header: str) -> str:
    return header[1:].strip().split()[0]


def collect_fasta_accessions(fasta_dir: Path) -> set[str]:
    return {
        header_accession(header)
        for fasta_path in fasta_dir.glob("*.fasta")
        for header, _sequence_lines in read_fasta(fasta_path)
    }


def write_fasta(path: Path, records: list[tuple[str, list[str]]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, sequence_lines in records:
            handle.write(f"{header}\n")
            for sequence_line in sequence_lines:
                handle.write(f"{sequence_line}\n")


def write_accession_manifest(path: Path, accessions: set[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["accession"])
        writer.writeheader()
        writer.writerows({"accession": accession} for accession in sorted(accessions))


def main() -> int:
    args = parse_args()
    metadata_dir = Path(args.metadata_dir).expanduser()
    fasta_dir = Path(args.fasta_dir).expanduser()

    if not metadata_dir.is_dir():
        raise RuntimeError(f"Metadata directory was not found: {metadata_dir}")
    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory was not found: {fasta_dir}")

    metadata_paths = sorted(metadata_dir.glob("RefID_*_metadata.csv"))
    input_accessions = collect_fasta_accessions(fasta_dir)
    fasta_by_refid = {
        refid: path
        for path in sorted(fasta_dir.glob("*.fasta"))
        if (refid := refid_from_fasta_path(path))
    }

    filtered_fasta_refids = 0
    missing_fasta_count = 0
    per_refid_counts: list[tuple[str, int, int]] = []

    for metadata_path in metadata_paths:
        refid = refid_from_metadata_path(metadata_path)
        if not refid:
            continue
        fasta_path = fasta_by_refid.get(refid)
        if fasta_path is None:
            missing_fasta_count += 1
            continue

        allowed_accessions = load_metadata_accessions(metadata_path)
        records = read_fasta(fasta_path)
        kept_records = [
            record
            for record in records
            if header_accession(record[0]) in allowed_accessions
        ]
        write_fasta(fasta_path, kept_records)

        filtered_fasta_refids += 1
        per_refid_counts.append((refid, len(records), len(kept_records)))

    output_accessions = collect_fasta_accessions(fasta_dir)
    kept_accessions_output = Path(args.kept_accessions_output).expanduser()
    write_accession_manifest(kept_accessions_output, output_accessions)
    print(f"staged_fasta_accessions_before_filter={len(input_accessions)}")
    print(f"staged_fasta_accessions_after_filter={len(output_accessions)}")
    print(
        f"staged_fasta_accessions_removed={len(input_accessions - output_accessions)}"
    )
    print(f"refid_fasta_files_filtered={filtered_fasta_refids}")
    print(f"filter_rules_without_matching_fasta={missing_fasta_count}")
    sorted_per_refid_counts = sorted(
        per_refid_counts, key=lambda item: refid_sort_key(item[0])
    )
    print(
        f"filtered_refids={','.join(refid for refid, _before_count, _after_count in sorted_per_refid_counts)}"
    )
    for refid, before_count, after_count in sorted_per_refid_counts:
        print(
            f"refid_filter_result=RefID:{refid},"
            f"BeforeRows:{before_count},AfterRows:{after_count},"
            f"RemovedRows:{before_count - after_count}"
        )
    print(f"kept_accessions_manifest={kept_accessions_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
