#!/usr/bin/env python3
"""Align each subtype per-gene nucleotide FASTA to its subtype-1a sequence."""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
GENES = ("NS3", "NS5A_NTD", "NS5B")
DEFAULT_INPUT_DIR = REPO_ROOT / "HCVData" / "Subtype-Ref"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "HCVData" / "Subtype-Ref"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line).upper())
    if header:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def is_subtype_1a(header: str) -> bool:
    return "subtype=1a" in header.split("|")


def align_gene(gene: str, input_dir: Path, output_dir: Path, temp_dir: Path) -> tuple[int, int, str]:
    input_path = input_dir / f"HCV_Subtype_Refs_{gene}_NA.fasta"
    records = read_fasta(input_path)
    anchors = [record for record in records if is_subtype_1a(record[0])]
    if not anchors:
        raise ValueError(f"{input_path} has no subtype=1a record")
    ordered_records = [anchors[0], *(record for record in records if record != anchors[0])]
    ordered_input = temp_dir / f"{gene}_with_1a_first.fasta"
    write_fasta(ordered_input, ordered_records)
    output_path = output_dir / f"HCV_Subtype_Refs_{gene}_NA_Aligned_to_1a.fasta"
    with output_path.open("w", encoding="utf-8") as output_handle:
        subprocess.run(
            ["mafft", "--auto", "--inputorder", str(ordered_input)],
            check=True,
            stdout=output_handle,
            stderr=subprocess.DEVNULL,
        )
    aligned = read_fasta(output_path)
    if len(aligned) != len(records) or not is_subtype_1a(aligned[0][0]):
        raise RuntimeError(f"{output_path} does not preserve the subtype-1a anchor and all input records")
    lengths = {len(sequence) for _, sequence in aligned}
    if len(lengths) != 1:
        raise RuntimeError(f"{output_path} is not a multiple sequence alignment")
    anchor_accession = next(
        item.split("=", 1)[1] for item in anchors[0][0].split("|") if item.startswith("accession=")
    )
    return len(aligned), lengths.pop(), anchor_accession


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--gene", choices=GENES, action="append", help="Align only this gene; repeat as needed.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    def display_path(path: Path) -> Path:
        try:
            return path.resolve().relative_to(REPO_ROOT)
        except ValueError:
            return path.resolve()

    with tempfile.TemporaryDirectory(prefix="hcv_subtype_1a_alignment_") as temp:
        temp_dir = Path(temp)
        for gene in args.gene or GENES:
            records, alignment_length, anchor_accession = align_gene(gene, args.input_dir, args.output_dir, temp_dir)
            print(f"{gene}_aligned_records={records}")
            print(f"{gene}_alignment_length={alignment_length}")
            print(f"{gene}_anchor_accession={anchor_accession}")
            print(f"{gene}_output={display_path(args.output_dir / f'HCV_Subtype_Refs_{gene}_NA_Aligned_to_1a.fasta')}")


if __name__ == "__main__":
    main()
