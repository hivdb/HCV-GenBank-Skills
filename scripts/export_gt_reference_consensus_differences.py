#!/usr/bin/env python3
"""Export aligned amino acids between GT references and Comet GT consensuses."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from Bio import Align


GENES = ("NS3", "NS5A_NTD", "NS5B")
VALID_AAS = set("ACDEFGHIKLMNPQRSTVWY*")


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    header: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(chunks).upper()
            header = line[1:].strip()
            chunks = []
        else:
            chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def load_references(path: Path) -> dict[tuple[str, str], str]:
    references: dict[tuple[str, str], str] = {}
    for header, sequence in read_fasta(path).items():
        match = re.fullmatch(r"genotype ([1-8]) \| (NS3|NS5A_NTD|NS5B)", header)
        if match:
            references[(match.group(1), match.group(2))] = sequence
    missing = [f"GT{gt} {gene}" for gt in "12345678" for gene in GENES if (gt, gene) not in references]
    if missing:
        raise ValueError(f"Missing reference sequences: {', '.join(missing)}")
    return references


def alignment_strings(reference: str, consensus: str, coordinates: object) -> tuple[str, str]:
    reference_chunks: list[str] = []
    consensus_chunks: list[str] = []
    ref_points = coordinates[0]  # type: ignore[index]
    consensus_points = coordinates[1]  # type: ignore[index]
    for index in range(len(ref_points) - 1):
        ref_start, ref_end = int(ref_points[index]), int(ref_points[index + 1])
        con_start, con_end = int(consensus_points[index]), int(consensus_points[index + 1])
        ref_span, con_span = ref_end - ref_start, con_end - con_start
        if ref_span and con_span:
            reference_chunks.append(reference[ref_start:ref_end])
            consensus_chunks.append(consensus[con_start:con_end])
        elif ref_span:
            reference_chunks.append(reference[ref_start:ref_end])
            consensus_chunks.append("-" * ref_span)
        elif con_span:
            reference_chunks.append("-" * con_span)
            consensus_chunks.append(consensus[con_start:con_end])
    return "".join(reference_chunks), "".join(consensus_chunks)


def covered_alignment_pairs(reference: str, consensus: str) -> dict[int, tuple[str, str]]:
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(reference, consensus)[0]
    aligned_reference, aligned_consensus = alignment_strings(reference, consensus, alignment.coordinates)
    pairs: dict[int, tuple[str, str]] = {}
    reference_position = 0
    for reference_aa, consensus_aa in zip(aligned_reference, aligned_consensus):
        if reference_aa != "-":
            reference_position += 1
        if reference_aa in VALID_AAS and consensus_aa in VALID_AAS:
            pairs[reference_position] = (reference_aa, consensus_aa)
    return pairs


def header_fields(header: str) -> dict[str, str]:
    return dict(field.split("=", 1) for field in header.split("|") if "=" in field)


def write_subtype_alignment(
    gene: str, reference_path: Path, consensus_path: Path, output_path: Path
) -> tuple[int, int]:
    consensuses = read_fasta(consensus_path)
    rows: list[tuple[str, str, str, dict[int, tuple[str, str]]]] = []
    skipped = 0
    for header, reference in read_fasta(reference_path).items():
        fields = header_fields(header)
        genotype = fields.get("genotype", "")
        subtype = fields.get("subtype", "")
        accession = fields.get("accession", "")
        consensus = consensuses.get(f"GT{genotype}_{subtype}")
        if not genotype or not subtype or consensus is None:
            skipped += 1
            continue
        rows.append((genotype, subtype, accession, covered_alignment_pairs(reference, consensus)))
    rows.sort(key=lambda row: (int(row[0]), row[1], row[2]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        for index, (genotype, subtype, accession, pairs) in enumerate(rows):
            positions = sorted(pairs)
            writer.writerow(["Gene", "Genotype", "Subtype", "ReferenceAccession", "Sequence"] + positions)
            writer.writerow(
                [gene, f"GT{genotype}", subtype, accession, "Reference"]
                + [pairs[position][0] for position in positions]
            )
            writer.writerow(
                [gene, f"GT{genotype}", subtype, accession, "Consensus"]
                + [
                    consensus_aa if consensus_aa != reference_aa else ""
                    for reference_aa, consensus_aa in (pairs[position] for position in positions)
                ]
            )
            if index != len(rows) - 1:
                writer.writerow([])
    return len(rows), skipped


def write_differences(gene: str, refs: dict[tuple[str, str], str], consensus_path: Path, output_path: Path) -> int:
    consensuses = read_fasta(consensus_path)
    sequences: dict[str, dict[int, tuple[str, str]]] = {}
    for genotype in "12345678":
        reference = refs[(genotype, gene)]
        consensus = consensuses.get(f"GT{genotype}")
        if consensus is None:
            raise ValueError(f"{consensus_path} has no GT{genotype} record")
        pairs = covered_alignment_pairs(reference, consensus)
        sequences[genotype] = pairs
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        for genotype in "12345678":
            pairs = sequences[genotype]
            # Include every position covered by the aligned consensus, not just mismatches.
            positions = sorted(pairs)
            writer.writerow(["Gene", "Genotype", "Sequence"] + positions)
            values = [pairs[position][0] for position in positions]
            writer.writerow([gene, f"GT{genotype}", "Reference"] + values)
            values = [
                consensus_aa if consensus_aa != reference_aa else ""
                for reference_aa, consensus_aa in (pairs[position] for position in positions)
            ]
            writer.writerow([gene, f"GT{genotype}", "Consensus"] + values)
            if genotype != "8":
                writer.writerow([])
    return sum(len(pairs) for pairs in sequences.values())


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-fasta", type=Path, required=True)
    parser.add_argument("--consensus-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--subtype-reference-dir",
        type=Path,
        help="Optional directory containing HCV_Subtype_Refs_<gene>_AA.fasta files.",
    )
    args = parser.parse_args()

    refs = load_references(args.reference_fasta)
    for gene in GENES:
        consensus_gene = "NS5A" if gene == "NS5A_NTD" else gene
        output = args.output_dir / f"HCV_GT_Ref_vs_Comet_GT_Consensus_Differences_{gene}.csv"
        count = write_differences(
            gene,
            refs,
            args.consensus_dir / f"{consensus_gene}_GT_Consensus.fasta",
            output,
        )
        print(f"{output}: {count} aligned reference-consensus amino-acid pairs")
    if args.subtype_reference_dir:
        for gene in GENES:
            consensus_gene = "NS5A" if gene == "NS5A_NTD" else gene
            reference_path = args.subtype_reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta"
            output = args.output_dir / f"HCV_Subtype_Ref_vs_Comet_Subtype_Consensus_Aligned_{gene}.csv"
            matched, skipped = write_subtype_alignment(
                gene,
                reference_path,
                args.consensus_dir / f"{consensus_gene}_Subtype_Consensus.fasta",
                output,
            )
            print(f"{output}: {matched} matched subtype references; {skipped} without a Comet consensus")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
