#!/usr/bin/env python3
"""Create an NS3 Comet genotype workbook with nucleotide-distance diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from openpyxl import Workbook


BLAST_OUTFMT = "6 qseqid sseqid length mismatch gaps pident evalue bitscore"
REFERENCE_GTS = tuple(str(index) for index in range(1, 9))
FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta-dir", required=True)
    parser.add_argument("--comet-genotype-csv", required=True)
    parser.add_argument("--reference-fasta", required=True, help="NS3 nucleotide references for GT1-GT8")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--temp-dir", default="temp/hcv-ns3-comet-build-workflow/build_ns3_comet_gt_allstudies")
    parser.add_argument("--min-aligned-nt", type=int, default=200)
    return parser.parse_args()


def fasta_accessions(path: Path) -> list[str]:
    return [line[1:].strip().split(maxsplit=1)[0] for line in path.read_text(encoding="utf-8").splitlines() if line.startswith(">")]


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks).upper()))
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(re.sub(r"\s+", "", line))
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def load_ns3_references(path: Path) -> dict[str, str]:
    refs: dict[str, str] = {}
    for header, sequence in parse_fasta(path):
        match = re.match(r"HCV([1-8])NS3(?:\s|$)", header)
        if match:
            refs[match.group(1)] = sequence
    missing = [gt for gt in REFERENCE_GTS if gt not in refs]
    if missing:
        raise RuntimeError(f"Missing NS3 nucleotide references for GTs: {', '.join(missing)}")
    return refs


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n{sequence}\n")


def build_reference_db(job_dir: Path, references: dict[str, str]) -> Path:
    fasta_path = job_dir / "ns3_gt_references.fasta"
    write_fasta(fasta_path, [(f"GT{gt}", references[gt]) for gt in REFERENCE_GTS])
    db_prefix = job_dir / "ns3_gt_references_db"
    subprocess.run(
        ["makeblastdb", "-in", str(fasta_path), "-dbtype", "nucl", "-out", str(db_prefix), "-parse_seqids"],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    return db_prefix


def nucleotide_distances(fasta_path: Path, db_prefix: Path, min_aligned_nt: int) -> dict[str, dict[str, Any]]:
    blast_path = fasta_path.with_suffix(fasta_path.suffix + ".ns3_gt_distance.tsv")
    subprocess.run(
        [
            "blastn", "-query", str(fasta_path), "-db", str(db_prefix), "-dust", "no", "-task", "blastn",
            "-evalue", "1e-6", "-max_hsps", "1", "-max_target_seqs", "8", "-outfmt", BLAST_OUTFMT,
            "-out", str(blast_path),
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        text=True,
    )
    best: dict[tuple[str, str], dict[str, Any]] = {}
    for line in blast_path.read_text(encoding="utf-8").splitlines():
        values = line.split("\t")
        if len(values) != 8:
            continue
        query, subject = values[0], values[1]
        gt_match = re.fullmatch(r"GT([1-8])", subject)
        if gt_match is None:
            continue
        length, mismatch, gaps, bitscore = int(values[2]), int(values[3]), int(values[4]), float(values[7])
        if length < min_aligned_nt:
            continue
        hit = {"distance": (mismatch + gaps) / length, "aligned_nt": length, "bitscore": bitscore}
        key = (query, gt_match.group(1))
        current = best.get(key)
        if current is None or (hit["distance"], -hit["aligned_nt"], -hit["bitscore"]) < (
            current["distance"], -current["aligned_nt"], -current["bitscore"]
        ):
            best[key] = hit
    blast_path.unlink(missing_ok=True)
    return {query: {gt: best[(query, gt)] for gt in REFERENCE_GTS if (query, gt) in best} for query, _ in parse_fasta(fasta_path)}


def load_comet(path: Path) -> dict[str, str]:
    assignments: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            accession = (row.get("accession") or "").strip()
            genotype = (row.get("genotype") or "").strip()
            if accession and genotype:
                assignments[accession] = genotype
                assignments.setdefault(accession.split(".", 1)[0], genotype)
    return assignments


def main() -> int:
    args = parse_args()
    assignments = load_comet(Path(args.comet_genotype_csv))
    if args.min_aligned_nt < 1:
        raise SystemExit("--min-aligned-nt must be at least 1")
    fasta_dir = Path(args.fasta_dir)
    if not fasta_dir.is_dir():
        raise RuntimeError(f"FASTA directory does not exist: {fasta_dir}")
    references = load_ns3_references(Path(args.reference_fasta))
    temp_dir = Path(args.temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)
    job_dir = Path(tempfile.mkdtemp(prefix="ns3_comet_gt_na_distance_", dir=temp_dir))
    try:
        db_prefix = build_reference_db(job_dir, references)
        rows: list[dict[str, Any]] = []
        for fasta_path in sorted(path for path in fasta_dir.iterdir() if path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS):
            refid, _, refname = fasta_path.stem.partition("_")
            distances_by_accession = nucleotide_distances(fasta_path, db_prefix, args.min_aligned_nt)
            for accession in fasta_accessions(fasta_path):
                genotype = assignments.get(accession) or assignments.get(accession.split(".", 1)[0])
                if genotype:
                    distances = distances_by_accession.get(accession, {})
                    comet_gt = re.match(r"([1-8])", genotype)
                    comet_hit = distances.get(comet_gt.group(1)) if comet_gt else None
                    rows.append({
                        "RefID": refid,
                        "RefName": refname,
                        "GenBankAccession": accession,
                        "BestGT": genotype,
                        "BestGTAssignmentSource": "Comet",
                        "BestGTDistance": comet_hit["distance"] if comet_hit else "",
                        "AlignedNT": comet_hit["aligned_nt"] if comet_hit else "",
                        **{f"GT{gt}_Distance": distances.get(gt, {}).get("distance", "") for gt in REFERENCE_GTS},
                        **{f"GT{gt}_AlignedNT": distances.get(gt, {}).get("aligned_nt", "") for gt in REFERENCE_GTS},
                    })
    finally:
        shutil.rmtree(job_dir, ignore_errors=True)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "NS3_GT_AllStudies.xlsx"
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "NS3_GT_AllStudies"
    fields = [
        "RefID", "RefName", "GenBankAccession", *[f"GT{gt}_Distance" for gt in REFERENCE_GTS],
        "BestGT", "BestGTAssignmentSource", "BestGTDistance", "AlignedNT",
        *[f"GT{gt}_AlignedNT" for gt in REFERENCE_GTS],
    ]
    sheet.append(fields)
    for row in rows:
        sheet.append([row[field] for field in fields])
    workbook.save(output_path)
    print(json.dumps({
        "combined_xlsx": str(output_path.resolve()),
        "master_row_count": len(rows),
        "comet_best_gt_count": len(rows),
        "distance_definition": "(BLAST nucleotide mismatches + gaps) / aligned nucleotides",
        "min_aligned_nt": args.min_aligned_nt,
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
