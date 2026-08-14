#!/usr/bin/env python3
"""Write accession-level non-COMET assignment and target-range coverage tables."""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
import time
from contextlib import ExitStack
from pathlib import Path


GENES = {"NS3": (36, 175), "NS5A": (26, 93), "NS5B": (150, 321)}
REPO_ROOT = Path(__file__).resolve().parents[2]
ASSIGNER = REPO_ROOT / "hcv-folder-genotype-subtype-assignment" / "scripts" / "assign_folder_genotype_subtype.py"
GT_REFERENCES = REPO_ROOT / "HCV_GT_RefSeqs.fasta"
SUBTYPE_REFERENCES = REPO_ROOT / "HCV_Subtype_Refs_By_Genome_NA.json"


def fasta_accessions(path: Path) -> list[str]:
    accessions: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            accessions.append(line[1:].split()[0].split(".")[0])
    if len(accessions) != len(set(accessions)):
        raise ValueError("Input FASTA contains duplicate accession headers.")
    return accessions


def assignments(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {row["accession"]: row for row in csv.DictReader(handle)}


def run_with_spinner(command: list[str], label: str) -> None:
    """Show activity while BLAST runs as a batch without per-record callbacks."""
    process = subprocess.Popen(command, stdout=subprocess.DEVNULL)
    frames = "⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏"
    frame = 0
    while process.poll() is None:
        print(f"\r{frames[frame % len(frames)]} {label}", end="", file=sys.stderr, flush=True)
        frame += 1
        time.sleep(0.2)
    if process.returncode:
        raise subprocess.CalledProcessError(process.returncode, command)
    print(f"\r✓ {label} complete", file=sys.stderr, flush=True)


def show_accession_progress(completed: int, total: int) -> None:
    width = 30
    filled = int(width * completed / total) if total else width
    bar = "█" * filled + " " * (width - filled)
    print(f"\r[{bar}] {completed:,}/{total:,} accessions", end="", file=sys.stderr, flush=True)


def coverage_hits(input_fasta: Path, work_dir: Path, threads: int) -> dict[str, dict[str, tuple[int, int, int, int]]]:
    work_dir.mkdir(parents=True, exist_ok=True)
    database = work_dir / "gt_refs"
    subprocess.run(["makeblastdb", "-in", str(GT_REFERENCES), "-dbtype", "nucl", "-out", str(database)], check=True, stdout=subprocess.DEVNULL)
    result = work_dir / "coverage.tsv"
    outfmt = "6 qseqid sseqid qstart qend sstart send bitscore length"
    run_with_spinner([
        "blastn", "-query", str(input_fasta), "-db", str(database), "-dust", "no", "-task", "blastn", "-num_threads", str(threads),
        "-evalue", "1e-6", "-max_hsps", "1", "-max_target_seqs", "1000", "-outfmt", outfmt,
        "-out", str(result),
    ], "mapping target-gene coverage")
    best: dict[str, dict[str, tuple[float, int, int, int, int]]] = {gene: {} for gene in GENES}
    for line in result.read_text(encoding="utf-8").splitlines():
        qid, sid, qstart, qend, sstart, send, bitscore, length = line.split("\t")
        gene = next((candidate for candidate in GENES if any(sid.split()[0].startswith(f"HCV{genotype}{candidate}") for genotype in range(1, 9))), None)
        if gene is None:
            continue
        score = (float(bitscore), int(length))
        if qid not in best[gene] or score > best[gene][qid][:2]:
            best[gene][qid] = (*score, int(qstart), int(qend), int(sstart), int(send))
    return {gene: {accession: values[2:] for accession, values in hits.items()} for gene, hits in best.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=REPO_ROOT / "outputs" / "local_alignment")
    parser.add_argument("--min-aligned-nt", type=int, default=200)
    parser.add_argument("--threads", type=int, default=4, help="BLAST worker threads to use (default: 4).")
    args = parser.parse_args()
    if args.threads < 1:
        raise SystemExit("--threads must be at least 1")
    input_fasta = args.input_fasta.resolve()
    accessions = fasta_accessions(input_fasta)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Auditing {len(accessions):,} FASTA records", file=sys.stderr)
    with tempfile.TemporaryDirectory(prefix="hcv_allseq_coverage_") as temp:
        temp_dir = Path(temp)
        assignment_dir = temp_dir / "assignments"
        run_with_spinner([
            sys.executable, str(ASSIGNER), "--fasta-dir", str(input_fasta.parent), "--output-dir", str(assignment_dir),
            "--gt-reference-fasta", str(GT_REFERENCES), "--subtype-json", str(SUBTYPE_REFERENCES),
            "--min-aligned-nt", str(args.min_aligned_nt), "--threads", str(args.threads),
        ], "assigning genotype and subtype")
        hits_by_gene = coverage_hits(input_fasta, temp_dir / "coverage", args.threads)
        fields = ["Accession", "ClosestGenotype", "ClosestSubtype", "ReferenceOverlapAA"]
        assignments_by_gene = {gene: assignments(assignment_dir / f"{gene}_assignments.csv") for gene in GENES}
        output_paths = {gene: args.output_dir / f"{gene}_AllSeq_NonComet_Coverage.csv" for gene in GENES}
        with ExitStack() as stack:
            writers = {}
            for gene, output in output_paths.items():
                handle = stack.enter_context(output.open("w", newline="", encoding="utf-8"))
                writers[gene] = csv.DictWriter(handle, fieldnames=fields)
                writers[gene].writeheader()
            show_accession_progress(0, len(accessions))
            for completed, accession in enumerate(accessions, start=1):
                for gene, (start_aa, end_aa) in GENES.items():
                    assignment = assignments_by_gene[gene].get(accession, {})
                    hits = hits_by_gene[gene]
                    target_start_nt, target_end_nt = (start_aa - 1) * 3 + 1, end_aa * 3
                    reference_overlap = query_overlap = ""
                    if accession in hits:
                        qstart, qend, sstart, send = hits[accession]
                        reference_start, reference_end = sorted((sstart, send))
                        overlap_start, overlap_end = max(reference_start, target_start_nt), min(reference_end, target_end_nt)
                        if overlap_start <= overlap_end:
                            aa_start = (overlap_start - 1) // 3 + 1
                            aa_end = (overlap_end - 1) // 3 + 1
                            reference_overlap = f"{aa_start}-{aa_end}"
                    writers[gene].writerow({
                        "Accession": accession, "ClosestGenotype": assignment.get("genotype", ""),
                        "ClosestSubtype": assignment.get("subtype", ""), "ReferenceOverlapAA": reference_overlap,
                    })
                show_accession_progress(completed, len(accessions))
            print(file=sys.stderr)
        for output in output_paths.values():
            print(output)


if __name__ == "__main__":
    main()
