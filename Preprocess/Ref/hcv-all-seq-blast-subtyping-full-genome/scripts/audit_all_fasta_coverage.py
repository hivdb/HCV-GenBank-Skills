#!/usr/bin/env python3
"""Write accession-level non-COMET assignment and target-range coverage tables."""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
import time
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path


GENES = {"NS3": (36, 175), "NS5A": (26, 93), "NS5B": (150, 321)}
# H77 full-genome nucleotide coordinates (one-based) for the first codon of
# each gene.  Coverage is evaluated in the reference coordinate system.
GENE_START_NT = {"NS3": 3420, "NS5A": 6258, "NS5B": 7602}
REPO_ROOT = Path(__file__).resolve().parents[4]
GT_REFERENCES = REPO_ROOT / "HCVData" / "Genotype-Ref" / "HCV_GT_FullGenome_RefSeqs.fasta"
SUBTYPE_REFERENCES = (
    REPO_ROOT / "HCVData" / "Subtype-Ref" / "HCV_Subtype_FullGenome_Refs.fasta"
)


def fasta_accessions(path: Path) -> list[str]:
    accessions: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            accessions.append(line[1:].split()[0].split(".")[0])
    if len(accessions) != len(set(accessions)):
        raise ValueError("Input FASTA contains duplicate accession headers.")
    return accessions


def fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = None
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if accession is not None:
                records[accession] = "".join(sequence).upper()
            accession = line[1:].split()[0].split(".")[0]
            sequence = []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line))
    if accession is not None:
        records[accession] = "".join(sequence).upper()
    return records


def fasta_header_records(path: Path) -> dict[str, str]:
    """Return sequences keyed by their complete FASTA header."""
    records: dict[str, str] = {}
    header = None
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records[header] = "".join(sequence).upper()
            header = line[1:].strip()
            sequence = []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line))
    if header is not None:
        records[header] = "".join(sequence).upper()
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for identifier, sequence in records.items():
            handle.write(f">{identifier}\n{sequence}\n")


def blast(
    query: Path, database: Path, output: Path, threads: int, label: str
) -> list[list[str]]:
    outfmt = "6 qseqid sseqid length mismatch gaps pident evalue bitscore"
    run_with_spinner(
        [
            "blastn",
            "-query",
            str(query),
            "-db",
            str(database),
            "-dust",
            "no",
            "-task",
            "blastn",
            "-num_threads",
            str(threads),
            "-evalue",
            "1e-6",
            "-max_hsps",
            "1",
            "-max_target_seqs",
            "1000",
            "-outfmt",
            outfmt,
            "-out",
            str(output),
        ],
        label,
    )
    return [
        line.split("\t")
        for line in output.read_text(encoding="utf-8").splitlines()
        if line
    ]


def make_database(input_fasta: Path, database: Path) -> None:
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(input_fasta),
            "-dbtype",
            "nucl",
            "-out",
            str(database),
            "-parse_seqids",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


def assignment_calls(
    input_fasta: Path, work_dir: Path, threads: int, min_aligned_nt: int
) -> dict[str, dict[str, dict[str, str]]]:
    """Assign full-genome genotype and subtype calls once, then reuse by gene."""
    work_dir.mkdir(parents=True, exist_ok=True)
    records = fasta_records(input_fasta)
    genotype_db = work_dir / "genotype_refs"
    make_database(GT_REFERENCES, genotype_db)
    genotype_hits = blast(
        input_fasta,
        genotype_db,
        work_dir / "genotype.tsv",
        threads,
        "assigning genotypes",
    )
    best_genotype: dict[str, tuple[tuple[float, int], str, list[str]]] = {}
    for hit in genotype_hits:
        query, subject = hit[0], hit[1]
        match = re.search(r"genotype=([1-8])(?:\||$)", subject)
        if not match or int(hit[2]) < min_aligned_nt:
            continue
        genotype, score = match.group(1), (float(hit[7]), int(hit[2]))
        if query not in best_genotype or score > best_genotype[query][0]:
            best_genotype[query] = (score, genotype, hit)

    subtype_by_genotype: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for header, sequence in fasta_header_records(SUBTYPE_REFERENCES).items():
        match = re.search(r"subtype=([1-8][A-Za-z0-9]*)", header)
        if match and sequence:
            subtype_by_genotype[match.group(1)[0]].append((match.group(1), sequence))

    query_ids_by_genotype: dict[str, list[str]] = defaultdict(list)
    for accession, (_, genotype, _) in best_genotype.items():
        query_ids_by_genotype[genotype].append(accession)

    calls: dict[str, dict[str, str]] = {}
    for genotype, accessions in query_ids_by_genotype.items():
        subtype_records = {
            f"S{index}": sequence
            for index, (_, sequence) in enumerate(subtype_by_genotype[genotype])
        }
        subtype_labels = {
            f"S{index}": subtype
            for index, (subtype, _) in enumerate(subtype_by_genotype[genotype])
        }
        subtype_fasta = work_dir / f"{genotype}_subtypes.fasta"
        subtype_db = work_dir / f"{genotype}_subtypes"
        write_fasta(subtype_fasta, subtype_records)
        make_database(subtype_fasta, subtype_db)
        query_fasta = work_dir / f"{genotype}_queries.fasta"
        write_fasta(query_fasta, {accession: records[accession] for accession in accessions})
        subtype_hits = blast(
            query_fasta,
            subtype_db,
            work_dir / f"{genotype}_subtypes.tsv",
            threads,
            f"assigning genotype {genotype} subtypes",
        )
        best_subtype: dict[str, tuple[tuple[float, int], list[str]]] = {}
        for hit in subtype_hits:
            score = (float(hit[7]), int(hit[2]))
            if int(hit[2]) >= min_aligned_nt and (
                hit[0] not in best_subtype or score > best_subtype[hit[0]][0]
            ):
                best_subtype[hit[0]] = (score, hit)
        for accession in accessions:
            _, assigned_genotype, genotype_hit = best_genotype[accession]
            call = {
                "genotype": assigned_genotype,
                "genotype_pident": genotype_hit[5],
                "subtype": "",
                "subtype_pident": "",
            }
            if accession in best_subtype:
                subtype_hit = best_subtype[accession][1]
                call.update(
                    subtype=subtype_labels[subtype_hit[1]],
                    subtype_pident=subtype_hit[5],
                )
            calls[accession] = call
    return {gene: dict(calls) for gene in GENES}


def run_with_spinner(command: list[str], label: str) -> None:
    """Show activity while BLAST runs as a batch without per-record callbacks."""
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL)
    frames = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
    frame = 0
    while process.poll() is None:
        print(
            f"\r{frames[frame % len(frames)]} {label}",
            end="",
            file=sys.stderr,
            flush=True,
        )
        frame += 1
        time.sleep(0.2)
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command)
    print(f"\r✓ {label} complete", file=sys.stderr, flush=True)


def show_accession_progress(completed: int, total: int) -> None:
    width = 30
    filled = int(width * completed / total) if total else width
    bar = "█" * filled + " " * (width - filled)
    print(
        f"\r[{bar}] {completed:,}/{total:,} accessions",
        end="",
        file=sys.stderr,
        flush=True,
    )


def coverage_hits(
    input_fasta: Path, work_dir: Path, threads: int
) -> dict[str, dict[str, tuple[int, int, int, int]]]:
    work_dir.mkdir(parents=True, exist_ok=True)
    database = work_dir / "gt_refs"
    subprocess.run(
        [
            "makeblastdb",
            "-in",
            str(GT_REFERENCES),
            "-dbtype",
            "nucl",
            "-out",
            str(database),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
    )
    result = work_dir / "coverage.tsv"
    outfmt = "6 qseqid sseqid qstart qend sstart send bitscore length"
    run_with_spinner(
        [
            "blastn",
            "-query",
            str(input_fasta),
            "-db",
            str(database),
            "-dust",
            "no",
            "-task",
            "blastn",
            "-num_threads",
            str(threads),
            "-evalue",
            "1e-6",
            "-max_hsps",
            "1",
            "-max_target_seqs",
            "1000",
            "-outfmt",
            outfmt,
            "-out",
            str(result),
        ],
        "mapping target-gene coverage",
    )
    best: dict[str, dict[str, tuple[float, int, int, int, int]]] = {
        gene: {} for gene in GENES
    }
    for line in result.read_text(encoding="utf-8").splitlines():
        qid, sid, qstart, qend, sstart, send, bitscore, length = line.split("\t")
        if not re.search(r"genotype=[1-8](?:\||$)", sid):
            continue
        score = (float(bitscore), int(length))
        for gene in GENES:
            if qid not in best[gene] or score > best[gene][qid][:2]:
                best[gene][qid] = (
                    *score,
                    int(qstart),
                    int(qend),
                    int(sstart),
                    int(send),
                )
    return {
        gene: {accession: values[2:] for accession, values in hits.items()}
        for gene, hits in best.items()
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "HCVData" / "nonComet-Full-genome",
    )
    parser.add_argument("--min-aligned-nt", type=int, default=100)
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="BLAST worker threads to use (default: 4).",
    )
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1")
    input_fasta = args.input_fasta.resolve()
    accessions = fasta_accessions(input_fasta)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Auditing {len(accessions):,} FASTA records", file=sys.stderr)
    with tempfile.TemporaryDirectory(prefix="hcv_allseq_coverage_") as temp:
        temp_dir = Path(temp)
        assignments_by_gene = assignment_calls(
            input_fasta, temp_dir / "assignments", args.threads, args.min_aligned_nt
        )
        hits_by_gene = coverage_hits(input_fasta, temp_dir / "coverage", args.threads)
        fields = [
            "Accession",
            "ClosestGenotype",
            "ClosestGenotypePident",
            "ClosestSubtype",
            "ClosestSubtypePident",
            "ReferenceOverlapAA",
            "FullyCover",
        ]
        output_paths = {
            gene: args.output_dir / f"{gene}_AllSeq_NonComet_Coverage.csv"
            for gene in GENES
        }
        with ExitStack() as stack:
            writers = {}
            for gene, output in output_paths.items():
                handle = stack.enter_context(
                    output.open("w", newline="", encoding="utf-8")
                )
                writers[gene] = csv.DictWriter(
                    handle, fieldnames=fields, lineterminator="\n"
                )
                writers[gene].writeheader()
            show_accession_progress(0, len(accessions))
            for completed, accession in enumerate(accessions, start=1):
                for gene, (start_aa, end_aa) in GENES.items():
                    assignment = assignments_by_gene[gene].get(accession, {})
                    hits = hits_by_gene[gene]
                    gene_start_nt = GENE_START_NT[gene]
                    target_start_nt = gene_start_nt + (start_aa - 1) * 3
                    target_end_nt = gene_start_nt + end_aa * 3 - 1
                    reference_overlap = fully_cover = ""
                    if accession in hits:
                        qstart, qend, sstart, send = hits[accession]
                        reference_start, reference_end = sorted((sstart, send))
                        overlap_start, overlap_end = (
                            max(reference_start, target_start_nt),
                            min(reference_end, target_end_nt),
                        )
                        if overlap_start <= overlap_end:
                            aa_start = (overlap_start - 1) // 3 + 1
                            aa_end = (overlap_end - 1) // 3 + 1
                            reference_overlap = f"{aa_start}-{aa_end}"
                            if aa_start == start_aa and aa_end == end_aa:
                                fully_cover = "Yes"
                    writers[gene].writerow(
                        {
                            "Accession": accession,
                            "ClosestGenotype": assignment.get("genotype", ""),
                            "ClosestGenotypePident": assignment.get(
                                "genotype_pident", ""
                            ),
                            "ClosestSubtype": assignment.get("subtype", ""),
                            "ClosestSubtypePident": assignment.get(
                                "subtype_pident", ""
                            ),
                            "ReferenceOverlapAA": reference_overlap,
                            "FullyCover": fully_cover,
                        }
                    )
                show_accession_progress(completed, len(accessions))
            print(file=sys.stderr)
        for output in output_paths.values():
            print(output)


if __name__ == "__main__":
    main()
