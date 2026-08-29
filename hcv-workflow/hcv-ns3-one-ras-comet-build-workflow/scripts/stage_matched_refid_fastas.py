#!/usr/bin/env python3
"""Copy matched RefID FASTA files into a pipeline staging directory."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]


def display_path(path: Path) -> str:
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--matched-files",
        required=True,
        help="Newline-delimited matched FASTA file list",
    )
    parser.add_argument(
        "--output-dir",
        required=True,
        help="Empty staging directory for copied FASTA files",
    )
    return parser.parse_args()


def accessions_in_fasta(path: Path) -> set[str]:
    accessions: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                accession = line[1:].strip().split(maxsplit=1)[0]
                if accession:
                    accessions.add(accession)
    return accessions


def write_accessions(path: Path, accessions: set[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["Accession"])
        writer.writerows((accession,) for accession in sorted(accessions))


def main() -> int:
    args = parse_args()
    matched_files_path = Path(args.matched_files).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve()
    if not matched_files_path.is_file():
        raise RuntimeError(
            f"Matched FASTA file list was not found: {matched_files_path}"
        )
    if not output_dir.is_dir():
        raise RuntimeError(f"Staging directory was not found: {output_dir}")

    source_paths = [
        Path(line).expanduser().resolve()
        for line in matched_files_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    for source_path in source_paths:
        if not source_path.is_file():
            raise RuntimeError(f"Matched FASTA file was not found: {source_path}")

    source_accessions: set[str] = set()
    staged_accessions: set[str] = set()
    for source_path in source_paths:
        source_accessions.update(accessions_in_fasta(source_path))
        destination = output_dir / source_path.name
        shutil.copy2(source_path, destination)
        staged_accessions.update(accessions_in_fasta(destination))

    accession_csv = output_dir.parent / "staged_accessions.csv"
    write_accessions(accession_csv, staged_accessions)

    print(f"staged_file_count={len(source_paths)}")
    print(f"staged_fasta_accessions_input={len(source_accessions)}")
    print(f"staged_fasta_accessions_output={len(staged_accessions)}")
    print(f"staged_accessions_csv={display_path(accession_csv)}")
    print(f"staged_fasta_dir={display_path(output_dir)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
