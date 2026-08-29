#!/usr/bin/env python3
"""Extract H77 NS3, NS5A, and NS5B nucleotide matches from an HCV FASTA."""

from __future__ import annotations

import argparse
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


GENES = ("NS3", "NS5A", "NS5B")
REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT = REPO_ROOT / "HCVData" / "HCV-all-seq-subtype" / "all.fasta"
DEFAULT_REFERENCE = REPO_ROOT / "HCVData" / "HCV-Ref-H77-Genotype1.fasta"
DEFAULT_OUTPUT = REPO_ROOT / "outputs" / "hcv-all-seq-h77-gene-na"
COMPLEMENT = str.maketrans("ACGTUNWSMKRYBDHVacgtunwsmkrybdhv", "TGCAANWSKMYRVHDBtgcaanwskmyrvhdb")


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if accession is not None:
                records[accession] = "".join(sequence).upper()
            accession = line[1:].split()[0].split(".")[0]
            if accession in records:
                raise ValueError(f"Duplicate accession in FASTA: {accession}")
            sequence = []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line))
    if accession is not None:
        records[accession] = "".join(sequence).upper()
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def best_hits(
    input_fasta: Path, reference_fasta: Path, work_dir: Path, threads: int,
    min_aligned_nt: int,
) -> dict[str, dict[str, list[str]]]:
    database = work_dir / "h77_genes"
    output = work_dir / "hits.tsv"
    run(
        [
            "makeblastdb", "-in", str(reference_fasta), "-dbtype", "prot",
            "-out", str(database), "-parse_seqids",
        ]
    )
    run(
        [
            "blastx", "-query", str(input_fasta), "-db", str(database),
            "-seg", "no", "-num_threads", str(threads),
            "-evalue", "1e-6", "-max_hsps", "1", "-max_target_seqs", "5",
            "-outfmt", "6 qseqid sseqid qstart qend sstart send bitscore length pident",
            "-out", str(output),
        ]
    )
    best: dict[str, dict[str, tuple[tuple[float, int], list[str]]]] = defaultdict(dict)
    for line in output.read_text(encoding="utf-8").splitlines():
        hit = line.split("\t")
        if len(hit) != 9 or hit[1] not in GENES or int(hit[7]) * 3 < min_aligned_nt:
            continue
        score = (float(hit[6]), int(hit[7]))
        previous = best[hit[0]].get(hit[1])
        if previous is None or score > previous[0]:
            best[hit[0]][hit[1]] = (score, hit)
    return {
        accession: {gene: item[1] for gene, item in hits.items()}
        for accession, hits in best.items()
    }


def write_gene_fastas(
    records: dict[str, str], hits: dict[str, dict[str, list[str]]], output_dir: Path
) -> dict[str, int]:
    output_dir.mkdir(parents=True, exist_ok=True)
    counts = {gene: 0 for gene in GENES}
    handles = {
        gene: (output_dir / f"HCV_AllSeq_H77_{gene}_NA.fasta").open(
            "w", encoding="utf-8"
        )
        for gene in GENES
    }
    try:
        for accession, sequence in records.items():
            for gene, hit in hits.get(accession, {}).items():
                qstart, qend = int(hit[2]), int(hit[3])
                start, end = sorted((qstart, qend))
                matched = sequence[start - 1 : end]
                if qstart > qend:
                    matched = reverse_complement(matched)
                handles[gene].write(
                    f">{accession}|gene={gene}|query={qstart}-{qend}|"
                    f"h77_aa={hit[4]}-{hit[5]}|extracted_nt={len(matched)}|"
                    f"aligned_aa={hit[7]}|pident={hit[8]}\n"
                    f"{matched}\n"
                )
                counts[gene] += 1
    finally:
        for handle in handles.values():
            handle.close()
    return counts


def split_gene_fastas(output_dir: Path, records_per_file: int) -> dict[str, int]:
    """Split each extracted gene FASTA into fixed-size record batches."""
    batch_counts: dict[str, int] = {}
    for gene in GENES:
        source = output_dir / f"HCV_AllSeq_H77_{gene}_NA.fasta"
        if not source.is_file():
            raise FileNotFoundError(f"Gene FASTA not found: {source}")
        gene_dir = output_dir / gene
        gene_dir.mkdir(parents=True, exist_ok=True)
        for stale_file in gene_dir.glob(f"HCV_AllSeq_H77_{gene}_NA_part_*.fasta"):
            stale_file.unlink()
        part = 0
        records_in_part = records_per_file
        handle = None
        try:
            for line in source.open(encoding="utf-8"):
                if line.startswith(">"):
                    if records_in_part == records_per_file:
                        if handle is not None:
                            handle.close()
                        part += 1
                        batch_path = (
                            gene_dir / f"HCV_AllSeq_H77_{gene}_NA_part_{part:04d}.fasta"
                        )
                        handle = batch_path.open("w", encoding="utf-8")
                        records_in_part = 0
                    records_in_part += 1
                if handle is not None:
                    handle.write(line)
        finally:
            if handle is not None:
                handle.close()
        batch_counts[gene] = part
    return batch_counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--h77-reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--min-aligned-nt", type=int, default=100)
    parser.add_argument("--records-per-file", type=int, default=1000)
    parser.add_argument(
        "--split-only",
        action="store_true",
        help="Split existing extracted gene FASTAs without rerunning BLASTX.",
    )
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1")
    if args.min_aligned_nt < 1:
        raise SystemExit("--min-aligned-nt must be at least 1")
    if args.records_per_file < 1:
        raise SystemExit("--records-per-file must be at least 1")
    input_fasta = args.input_fasta.resolve()
    reference_fasta = args.h77_reference.resolve()
    output_dir = args.output_dir.resolve()
    if args.split_only:
        batch_counts = split_gene_fastas(output_dir, args.records_per_file)
        for gene in GENES:
            print(f"{gene}_batch_files={batch_counts[gene]}")
            print(output_dir / gene)
        return
    if not input_fasta.is_file() or not reference_fasta.is_file():
        raise SystemExit("Input FASTA and H77 reference FASTA must both exist")
    records = parse_fasta(input_fasta)
    with tempfile.TemporaryDirectory(prefix="hcv_h77_gene_na_") as temp:
        hits = best_hits(
            input_fasta, reference_fasta, Path(temp), args.threads, args.min_aligned_nt
        )
    counts = write_gene_fastas(records, hits, output_dir)
    batch_counts = split_gene_fastas(output_dir, args.records_per_file)
    for gene in GENES:
        print(f"{gene}_matched_accessions={counts[gene]}")
        print(output_dir / f"HCV_AllSeq_H77_{gene}_NA.fasta")
        print(f"{gene}_batch_files={batch_counts[gene]}")
        print(output_dir / gene)


if __name__ == "__main__":
    main()
