#!/usr/bin/env python3
"""Apply ICTV subtype names and add missing HCV subtype AA references from GenBank."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from Bio import Align, SeqIO
from Bio.Seq import Seq


ADDITIONS = {
    "4a": ("4", "Y11604"),
    "5b": ("5", "PQ899568"),
    "6b": ("6", "D84262"),
    "7a": ("7", "EF108306"),
}
GENES = ("NS3", "NS5A_NTD", "NS5B")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks)))
            header, chunks = line[1:], []
        elif line:
            chunks.append(line)
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def reference_aa_by_gene(gt_reference_fasta: Path) -> dict[tuple[str, str], str]:
    refs: dict[tuple[str, str], str] = {}
    for header, nucleotide_sequence in read_fasta(gt_reference_fasta):
        match = re.match(r"HCV([1-8])(NS3|NS5A|NS5B)(?:\s|$)", header)
        if not match:
            continue
        genotype, gene = match.groups()
        aa = str(Seq(nucleotide_sequence).translate()).rstrip("*")
        refs[(genotype, "NS5A_NTD" if gene == "NS5A" else gene)] = (
            aa[:213] if gene == "NS5A" else aa
        )
    return refs


def polyprotein_translation(genbank_path: Path) -> str:
    record = SeqIO.read(genbank_path, "genbank")
    for feature in record.features:
        if feature.type == "CDS" and "translation" in feature.qualifiers:
            return str(feature.qualifiers["translation"][0]).replace(" ", "").upper()
    raise ValueError(f"No translated CDS found in {genbank_path}")


def extract_gene(reference_aa: str, polyprotein_aa: str) -> str:
    aligner = Align.PairwiseAligner(mode="local")
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -5.0
    aligner.extend_gap_score = -1.0
    alignment = aligner.align(reference_aa, polyprotein_aa)[0]
    start, end = int(alignment.coordinates[1][0]), int(alignment.coordinates[1][-1])
    sequence = polyprotein_aa[start:end]
    if len(sequence) < len(reference_aa) * 0.8:
        raise ValueError(
            f"Low-coverage extraction: {len(sequence)} AA for a {len(reference_aa)} AA reference"
        )
    return sequence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, required=True)
    parser.add_argument("--gt-reference-fasta", type=Path, required=True)
    parser.add_argument("--genbank-dir", type=Path, required=True)
    args = parser.parse_args()

    gene_refs = reference_aa_by_gene(args.gt_reference_fasta)
    additions: dict[str, dict[str, str]] = {}
    for subtype, (genotype, accession) in ADDITIONS.items():
        polyprotein = polyprotein_translation(args.genbank_dir / f"{accession}.gb")
        additions[subtype] = {
            gene: extract_gene(gene_refs[(genotype, gene)], polyprotein)
            for gene in GENES
        }

    for gene in GENES:
        path = args.reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta"
        records = []
        for header, sequence in read_fasta(path):
            header = (
                header.replace("subtype=4w3", "subtype=4w")
                .replace("subtype=6xa5", "subtype=6xa")
                .replace("genotypeName=Genotype4w3", "genotypeName=Genotype4w")
                .replace("genotypeName=Genotype6xa5", "genotypeName=Genotype6xa")
            )
            records.append((header, sequence))
        existing = {
            header.split("|accession=")[1].split("|")[0]
            for header, _ in records
            if "|accession=" in header
        }
        for subtype, (genotype, accession) in ADDITIONS.items():
            if accession not in existing:
                header = (
                    f"gene={gene}|genotype={genotype}|subtype={subtype}|accession={accession}"
                    f"|genotypeName=Genotype{subtype}|source=ICTV_Jan2026"
                )
                records.append((header, additions[subtype][gene]))
        records.sort(
            key=lambda item: (
                int(re.search(r"genotype=(\d+)", item[0]).group(1)),
                item[0],
            )
        )
        write_fasta(path, records)
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
