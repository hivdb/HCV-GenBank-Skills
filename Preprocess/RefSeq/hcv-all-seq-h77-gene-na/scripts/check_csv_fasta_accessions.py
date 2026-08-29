#!/usr/bin/env python3
"""Find CSV and FASTA files with identical accession coverage."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


def accession(value: str) -> str:
    return value.strip().split("|", 1)[0].split()[0]


def csv_accessions(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError("CSV has no header")
        column = next(
            (name for name in reader.fieldnames if name.lower() in {"accession", "name"}),
            None,
        )
        if column is None:
            raise ValueError("CSV needs an Accession or name column")
        return {accession(row[column]) for row in reader if row.get(column, "").strip()}


def fasta_accessions(path: Path) -> set[str]:
    return {
        accession(line[1:])
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, required=True)
    parser.add_argument(
        "--rename-csv-to-fasta",
        action="store_true",
        help="Rename each uniquely matched CSV to its matching FASTA basename.",
    )
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not directory.is_dir():
        raise SystemExit(f"Directory not found: {directory}")
    csv_files = sorted(directory.glob("*.csv"))
    fasta_files = sorted(directory.glob("*.fasta"))
    if not csv_files:
        raise SystemExit(f"No CSV files found in {directory}")
    csv_sets: dict[Path, frozenset[str]] = {}
    errors = 0
    for csv_file in csv_files:
        try:
            csv_sets[csv_file] = frozenset(csv_accessions(csv_file))
        except ValueError as exc:
            print(f"ERROR {csv_file.name}: {exc}")
            errors += 1
    fasta_sets = {
        fasta_file: frozenset(fasta_accessions(fasta_file))
        for fasta_file in fasta_files
    }
    fastas_by_accessions: dict[frozenset[str], list[Path]] = {}
    for fasta_file, accessions in fasta_sets.items():
        fastas_by_accessions.setdefault(accessions, []).append(fasta_file)
    matched_fastas: set[Path] = set()
    matched_pairs: list[tuple[Path, Path]] = []
    for csv_file, accessions in csv_sets.items():
        matches = fastas_by_accessions.get(accessions, [])
        if len(matches) == 1:
            fasta_file = matches[0]
            matched_fastas.add(fasta_file)
            matched_pairs.append((csv_file, fasta_file))
            print(
                f"MATCH {csv_file.name} <-> {fasta_file.name}: "
                f"{len(accessions)} accessions"
            )
        elif not matches:
            print(f"UNMATCHED CSV {csv_file.name}")
            errors += 1
        else:
            names = ", ".join(path.name for path in matches)
            print(f"AMBIGUOUS CSV {csv_file.name}: {names}")
            errors += 1
    for fasta_file in fasta_files:
        if fasta_file not in matched_fastas:
            print(f"UNMATCHED FASTA {fasta_file.name}")
            errors += 1
    if errors:
        raise SystemExit(f"Unmatched or ambiguous files: {errors}")
    if args.rename_csv_to_fasta:
        temporary_paths = [
            csv_file.with_name(f".{csv_file.name}.matching-rename")
            for csv_file, _ in matched_pairs
        ]
        if any(path.exists() for path in temporary_paths):
            raise SystemExit("Temporary rename file already exists; no files were renamed")
        for (csv_file, _), temporary_path in zip(
            matched_pairs, temporary_paths, strict=True
        ):
            csv_file.rename(temporary_path)
        for (_, fasta_file), temporary_path in zip(
            matched_pairs, temporary_paths, strict=True
        ):
            target = fasta_file.with_suffix(".csv")
            temporary_path.rename(target)
            print(f"RENAMED {target.name}")


if __name__ == "__main__":
    main()
