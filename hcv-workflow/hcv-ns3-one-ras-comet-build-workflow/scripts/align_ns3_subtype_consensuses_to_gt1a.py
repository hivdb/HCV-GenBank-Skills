#!/usr/bin/env python3
"""Align NS3 subtype consensus sequences to the fixed GT1_1a coordinate system."""

from __future__ import annotations

import argparse
import subprocess
import tempfile
from pathlib import Path


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    sequence: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(sequence)))
            header = line[1:].strip()
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        records.append((header, "".join(sequence)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n{sequence}\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--output-fasta", type=Path, required=True)
    parser.add_argument("--reference-header", default="GT1_1a")
    parser.add_argument("--mafft-bin", default="mafft")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    records = read_fasta(args.input_fasta)
    reference_records = [record for record in records if record[0] == args.reference_header]
    if len(reference_records) != 1:
        raise RuntimeError(
            f"Expected exactly one reference record named {args.reference_header!r}, found {len(reference_records)}"
        )
    other_records = [record for record in records if record[0] != args.reference_header]
    if not other_records:
        raise RuntimeError("At least one non-reference subtype consensus is required")

    with tempfile.TemporaryDirectory(prefix="ns3_gt1a_alignment_") as directory:
        temp_dir = Path(directory)
        reference_path = temp_dir / "reference.fasta"
        additions_path = temp_dir / "other_consensuses.fasta"
        write_fasta(reference_path, reference_records)
        write_fasta(additions_path, other_records)
        result = subprocess.run(
            [args.mafft_bin, "--quiet", "--add", str(additions_path), "--keeplength", str(reference_path)],
            check=True,
            text=True,
            capture_output=True,
        )

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    args.output_fasta.write_text(result.stdout, encoding="utf-8")
    aligned_records = read_fasta(args.output_fasta)
    lengths = {len(sequence) for _, sequence in aligned_records}
    if len(aligned_records) != len(records) or len(lengths) != 1:
        raise RuntimeError("MAFFT output does not contain one equally sized alignment record per input consensus")
    aligned_reference = aligned_records[0] if aligned_records else ("", "")
    if aligned_reference[0] != args.reference_header or "-" in aligned_reference[1]:
        raise RuntimeError("MAFFT did not preserve the requested gap-free reference sequence as the first record")

    print(f"output_fasta={args.output_fasta}")
    print(f"reference_header={args.reference_header}")
    print(f"record_count={len(aligned_records)}")
    print(f"aligned_length={lengths.pop()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
