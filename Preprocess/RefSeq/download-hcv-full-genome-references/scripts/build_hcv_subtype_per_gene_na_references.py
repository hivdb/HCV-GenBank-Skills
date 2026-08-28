#!/usr/bin/env python3
"""Extract subtype per-gene NA references by BLASTX against per-gene AA references."""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
GENES = ("NS3", "NS5A_NTD", "NS5B")
DEFAULT_AA_REFERENCE_DIR = REPO_ROOT / "HCVData" / "Reference_seqs"
DEFAULT_FULL_GENOME_FASTA = REPO_ROOT / "HCVData" / "Subtype-Ref" / "HCV_Subtype_FullGenome_Refs.fasta"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "HCVData" / "Subtype-Ref"
MIN_AA_COVERAGE = 0.80


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence).upper()))
            header, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line).upper())
    if header:
        records.append((header, "".join(sequence).upper()))
    return records


def fasta_fields(header: str) -> dict[str, str]:
    return {key: value for key, value in (item.split("=", 1) for item in header.split("|") if "=" in item)}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"))[::-1]


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")


def write_aa_database_fasta(path: Path, records: list[tuple[str, str]]) -> dict[str, tuple[str, str]]:
    by_accession: dict[str, tuple[str, str]] = {}
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            accession = fasta_fields(header).get("accession", "")
            if not accession:
                raise ValueError(f"AA reference header has no accession: {header}")
            if accession in by_accession:
                raise ValueError(f"Duplicate AA reference accession: {accession}")
            by_accession[accession] = (header, sequence)
            handle.write(f">{accession}\n{sequence}\n")
    return by_accession


def build_gene_records(
    gene: str, aa_reference_dir: Path, full_genome_fasta: Path, output_dir: Path, work_dir: Path
) -> tuple[int, int]:
    aa_path = aa_reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta"
    aa_by_accession = write_aa_database_fasta(work_dir / f"{gene}_aa.fasta", read_fasta(aa_path))
    genomes = {
        fasta_fields(header).get("accession", ""): sequence
        for header, sequence in read_fasta(full_genome_fasta)
    }
    query_fasta = work_dir / f"{gene}_genomes.fasta"
    write_fasta(query_fasta, [(accession, sequence) for accession, sequence in genomes.items() if accession])
    database = work_dir / f"{gene}_aa_db"
    subprocess.run(["makeblastdb", "-in", str(work_dir / f"{gene}_aa.fasta"), "-dbtype", "prot", "-out", str(database)], check=True, stdout=subprocess.DEVNULL)
    result_path = work_dir / f"{gene}_blastx.tsv"
    subprocess.run(
        [
            "blastx", "-query", str(query_fasta), "-db", str(database), "-evalue", "1e-6",
            "-max_hsps", "1", "-max_target_seqs", "1000", "-outfmt",
            "6 qseqid sseqid qstart qend pident length bitscore qframe", "-out", str(result_path),
        ],
        check=True,
    )
    best: dict[str, tuple[float, int, int, int, int]] = {}
    for line in result_path.read_text(encoding="utf-8").splitlines():
        qid, sid, qstart, qend, _pident, length, bitscore, qframe = line.split("\t")
        if qid != sid:
            continue
        score = (float(bitscore), int(length))
        if qid not in best or score > best[qid][:2]:
            best[qid] = (*score, int(qstart), int(qend), int(qframe))

    output_records: list[tuple[str, str]] = []
    missing: list[str] = []
    for accession, (header, aa_sequence) in aa_by_accession.items():
        if accession not in genomes or accession not in best:
            missing.append(accession)
            continue
        _bitscore, aligned_aa, qstart, qend, qframe = best[accession]
        coverage = aligned_aa / len(aa_sequence)
        if coverage < MIN_AA_COVERAGE:
            missing.append(accession)
            continue
        start, end = min(qstart, qend) - 1, max(qstart, qend)
        nucleotide = genomes[accession][start:end]
        if qframe < 0:
            nucleotide = reverse_complement(nucleotide)
        if len(nucleotide) % 3:
            missing.append(accession)
            continue
        output_records.append(
            (f"{header}|fullGenomeNA={start + 1}-{end}|alignmentAACoverage={coverage:.3f}", nucleotide)
        )
    if missing:
        raise RuntimeError(f"{gene}: unable to extract {len(missing)} accession(s): {', '.join(missing)}")
    output_path = output_dir / f"HCV_Subtype_Refs_{gene}_NA.fasta"
    write_fasta(output_path, output_records)
    return len(output_records), sum(len(sequence) for _, sequence in output_records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aa-reference-dir", type=Path, default=DEFAULT_AA_REFERENCE_DIR)
    parser.add_argument("--full-genome-fasta", type=Path, default=DEFAULT_FULL_GENOME_FASTA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    def display_path(path: Path) -> Path:
        try:
            return path.resolve().relative_to(REPO_ROOT)
        except ValueError:
            return path.resolve()

    with tempfile.TemporaryDirectory(prefix="hcv_subtype_gene_na_") as temp:
        work_dir = Path(temp)
        for gene in GENES:
            count, nucleotides = build_gene_records(
                gene, args.aa_reference_dir, args.full_genome_fasta, args.output_dir, work_dir
            )
            print(f"{gene}_reference_count={count}")
            print(f"{gene}_nucleotide_count={nucleotides}")
            print(f"{gene}_output={display_path(args.output_dir / f'HCV_Subtype_Refs_{gene}_NA.fasta')}")


if __name__ == "__main__":
    main()
