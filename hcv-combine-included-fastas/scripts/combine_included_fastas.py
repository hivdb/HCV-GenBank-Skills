#!/usr/bin/env python3
"""Combine staged HCV workflow FASTA files into one FASTA per gene."""

from __future__ import annotations

import argparse
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SKILL_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Combine NS3, NS5A, and NS5B included_refid_fastas directories."
    )
    parser.add_argument(
        "--ns3-dir",
        type=Path,
        default=REPO_ROOT / "temp/hcv-ns3-build-workflow/run_ns3_pipeline/included_refid_fastas",
    )
    parser.add_argument(
        "--ns5a-dir",
        type=Path,
        default=REPO_ROOT / "temp/hcv-ns5a-build-workflow/run_ns5a_pipeline/included_refid_fastas",
    )
    parser.add_argument(
        "--ns5b-dir",
        type=Path,
        default=REPO_ROOT / "temp/hcv-ns5b-build-workflow/run_ns5b_pipeline/included_refid_fastas",
    )
    parser.add_argument("--output-dir", type=Path, default=SKILL_ROOT / "assets")
    return parser.parse_args()


def combine(gene: str, source_dir: Path, output_dir: Path) -> tuple[Path, int, int]:
    if not source_dir.is_dir():
        raise RuntimeError(f"{gene} input directory was not found: {source_dir}")

    input_paths = sorted(source_dir.glob("*.fasta"))
    if not input_paths:
        raise RuntimeError(f"No .fasta files were found for {gene}: {source_dir}")

    output_path = output_dir / f"{gene}.fasta"
    record_count = 0
    with output_path.open("w", encoding="utf-8", newline="\n") as output_handle:
        for input_path in input_paths:
            text = input_path.read_text(encoding="utf-8")
            if text and not text.endswith("\n"):
                text += "\n"
            output_handle.write(text)
            record_count += sum(1 for line in text.splitlines() if line.startswith(">"))
    return output_path, len(input_paths), record_count


def main() -> int:
    args = parse_args()
    output_dir = args.output_dir.expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)

    for gene, source_dir in (("NS3", args.ns3_dir), ("NS5A", args.ns5a_dir), ("NS5B", args.ns5b_dir)):
        output_path, input_count, record_count = combine(gene, source_dir.expanduser(), output_dir)
        print(f"{gene}: input_files={input_count} records={record_count} output={output_path.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
