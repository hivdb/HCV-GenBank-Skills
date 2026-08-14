#!/usr/bin/env python3
"""Pairwise-align subtype AA references to the first FASTA record.

This intentionally does *not* make a multiple sequence alignment.  For each
input record after the first, the output contains one independently aligned
reference/query FASTA pair.  The reference record is repeated once for every
query so inserted gaps remain specific to that pairwise comparison.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from Bio import Align


DEFAULT_GENES = ("NS3", "NS5A_NTD", "NS5B")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=Path("HCVData/Reference_seqs"),
        help="Directory containing HCV_Subtype_Refs_<GENE>_AA.fasta files.",
    )
    parser.add_argument(
        "--gene",
        choices=DEFAULT_GENES,
        action="append",
        help="Gene to process. Repeat to select several; default processes all genes.",
    )
    parser.add_argument(
        "--output-suffix",
        default="_FirstSeq_Pairwise_Aligned.fasta",
        help="Suffix appended after HCV_Subtype_Refs_<GENE>_AA in output names.",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    if len(records) < 2:
        raise RuntimeError(f"Expected at least two records in {path}")
    return records


def aligned_strings(alignment: Align.Alignment, target: str, query: str) -> tuple[str, str]:
    """Return Biopython's gapped target and query rows.

    ``Alignment.coordinates`` can combine an indel and matching residues in a
    single segment, so indexing the alignment is safer than expanding those
    coordinates ourselves.
    """
    return str(alignment[0]), str(alignment[1])


def write_record(handle, header: str, sequence: str) -> None:
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 70):
        handle.write(sequence[start : start + 70] + "\n")


def build_aligner() -> Align.PairwiseAligner:
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner


def process_gene(reference_dir: Path, gene: str, suffix: str) -> tuple[Path, int]:
    input_path = reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta"
    output_path = reference_dir / f"HCV_Subtype_Refs_{gene}_AA{suffix}"
    records = read_fasta(input_path)
    reference_header, reference_sequence = records[0]
    aligner = build_aligner()

    with output_path.open("w", encoding="utf-8") as handle:
        for pair_index, (query_header, query_sequence) in enumerate(records[1:], start=1):
            alignment = aligner.align(reference_sequence, query_sequence)[0]
            aligned_reference, aligned_query = aligned_strings(
                alignment, reference_sequence, query_sequence
            )
            if len(aligned_reference) != len(aligned_query):
                raise RuntimeError(f"Unequal pairwise alignment lengths for {query_header}")
            write_record(
                handle,
                f"{reference_header}|pairwise_target_for={pair_index}",
                aligned_reference,
            )
            write_record(handle, f"{query_header}|pairwise_query={pair_index}", aligned_query)
    return output_path, len(records) - 1


def main() -> int:
    args = parse_args()
    reference_dir = args.reference_dir.expanduser()
    genes = args.gene or DEFAULT_GENES
    for gene in genes:
        output_path, pair_count = process_gene(reference_dir, gene, args.output_suffix)
        print(f"{gene}: {pair_count} pairwise alignments -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
