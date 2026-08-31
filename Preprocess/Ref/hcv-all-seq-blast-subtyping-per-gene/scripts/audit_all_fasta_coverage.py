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
import xlsxwriter
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path


GENES = {"NS3": (36, 175), "NS5A": (26, 93), "NS5B": (150, 321)}
REPO_ROOT = Path(__file__).resolve().parents[4]
GT_REFERENCES = {
    "NS3": REPO_ROOT / "HCVData" / "Genotype-Ref" / "HCV_GT_Refs_NS3_NA.fasta",
    "NS5A": REPO_ROOT / "HCVData" / "Genotype-Ref" / "HCV_GT_Refs_NS5A_NA.fasta",
    "NS5B": REPO_ROOT / "HCVData" / "Genotype-Ref" / "HCV_GT_Refs_NS5B_NA.fasta",
}
SUBTYPE_REFERENCES = {
    "NS3": REPO_ROOT / "HCVData" / "Subtype-Ref" / "HCV_Subtype_Refs_NS3_NA.fasta",
    "NS5A": REPO_ROOT
    / "HCVData"
    / "Subtype-Ref"
    / "HCV_Subtype_Refs_NS5A_NTD_NA.fasta",
    "NS5B": REPO_ROOT / "HCVData" / "Subtype-Ref" / "HCV_Subtype_Refs_NS5B_NA.fasta",
}


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


def distance(hit: list[str]) -> str:
    """Return BLAST percent-identity distance, rounded for CSV output."""
    return f"{100 - float(hit[5]):.3f}"


def ranked_hits(
    hits: list[list[str]], label_for_subject, min_aligned_nt: int
) -> dict[str, list[tuple[str, list[str]]]]:
    """Keep the strongest hit for each query and reference label."""
    best: dict[str, dict[str, tuple[tuple[float, int], list[str]]]] = defaultdict(dict)
    for hit in hits:
        label = label_for_subject(hit[1])
        if not label or int(hit[2]) < min_aligned_nt:
            continue
        score = (float(hit[7]), int(hit[2]))
        current = best[hit[0]].get(label)
        if current is None or score > current[0]:
            best[hit[0]][label] = (score, hit)
    return {
        accession: [
            (label, item[1])
            for label, item in sorted(values.items(), key=lambda item: item[1][0], reverse=True)
        ]
        for accession, values in best.items()
    }


def write_choice_report(
    path: Path, labels: list[str], label_kind: str, calls: dict[str, dict[str, str]],
    distance_key: str, assigned_genotype: bool = False,
) -> None:
    distance_fields = [f"{label_kind}{label}Distance" for label in labels]
    fields = ["Accession", "ClosestSubtypeAlignedNT"]
    if assigned_genotype:
        fields.append("AssignedGenotype")
    fields += distance_fields + [
        f"FirstChoice{label_kind}",
        "FirstChoiceDistance",
        f"SecondChoice{label_kind}",
        "SecondChoiceDistance",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for accession, call in calls.items():
            distances = call.get(distance_key, {})
            row = {
                "Accession": accession,
                "ClosestSubtypeAlignedNT": call.get("subtype_aligned_nt", ""),
            }
            if assigned_genotype:
                row["AssignedGenotype"] = call.get("genotype", "")
            row.update(
                {f"{label_kind}{label}Distance": distances.get(label, "") for label in labels}
            )
            choices = call.get(f"{distance_key}_choices", [])
            for number, choice in enumerate(choices[:2], start=1):
                row[f"{'First' if number == 1 else 'Second'}Choice{label_kind}"] = choice[0]
                row[f"{'First' if number == 1 else 'Second'}ChoiceDistance"] = choice[1]
            writer.writerow(row)


def write_choice_sheet(
    workbook, sheet_name: str, labels: list[str], label_kind: str,
    calls: dict[str, dict[str, str]], distance_key: str, assigned_genotype: bool,
) -> None:
    distance_fields = [f"{label_kind}{label}Distance" for label in labels]
    fields = ["Accession", "ClosestSubtypeAlignedNT"]
    if assigned_genotype:
        fields.append("AssignedGenotype")
    fields += distance_fields + [
        f"FirstChoice{label_kind}", "FirstChoiceDistance",
        f"SecondChoice{label_kind}", "SecondChoiceDistance",
    ]
    worksheet = workbook.add_worksheet(sheet_name)
    header_format = workbook.add_format({"bold": True, "bg_color": "#5B9BD5"})
    worksheet.write_row(0, 0, fields, header_format)
    rows = []
    for accession, call in calls.items():
        distances = call.get(distance_key, {})
        row = {
            "Accession": accession,
            "ClosestSubtypeAlignedNT": call.get("subtype_aligned_nt", ""),
        }
        if assigned_genotype:
            row["AssignedGenotype"] = call.get("genotype", "")
        row.update({f"{label_kind}{label}Distance": distances.get(label, "") for label in labels})
        for number, choice in enumerate(call.get(f"{distance_key}_choices", [])[:2], start=1):
            prefix = "First" if number == 1 else "Second"
            row[f"{prefix}Choice{label_kind}"] = choice[0]
            row[f"{prefix}ChoiceDistance"] = choice[1]
        rows.append(row)
    for row_number, row in enumerate(rows, start=1):
        worksheet.write_row(row_number, 0, [row.get(field, "") for field in fields])
    worksheet.freeze_panes(1, 0)
    worksheet.autofilter(0, 0, max(len(rows), 1), len(fields) - 1)
    for column, field in enumerate(fields):
        width = max([len(field), *[len(str(row.get(field, ""))) for row in rows]]) + 2
        worksheet.set_column(column, column, min(width, 40))


def write_distance_workbook(
    path: Path, genotype_labels: list[str], subtype_labels_by_genotype: dict[str, list[str]],
    calls: dict[str, dict[str, str]],
) -> None:
    with xlsxwriter.Workbook(str(path)) as workbook:
        write_choice_sheet(
            workbook, "Genotype", genotype_labels, "Genotype", calls,
            "genotype_distances", False,
        )
        for genotype in genotype_labels:
            calls_for_genotype = {
                accession: call for accession, call in calls.items()
                if call.get("genotype") == genotype
            }
            write_choice_sheet(
                workbook, f"Subtype_{genotype}",
                subtype_labels_by_genotype.get(genotype, []), "Subtype",
                calls_for_genotype, "subtype_distances", True,
            )


def assignment_calls(
    input_fasta: Path, work_dir: Path, threads: int, min_aligned_nt: int
) -> tuple[
    dict[str, dict[str, dict[str, str]]],
    dict[str, list[str]],
    dict[str, dict[str, list[str]]],
]:
    """Assign each gene using its genotype and subtype nucleotide references."""
    work_dir.mkdir(parents=True, exist_ok=True)
    records = fasta_records(input_fasta)
    results: dict[str, dict[str, dict[str, str]]] = {gene: {} for gene in GENES}
    genotype_labels_by_gene: dict[str, list[str]] = {}
    subtype_labels_by_gene: dict[str, dict[str, list[str]]] = {}
    for gene in GENES:
        genotype_db = work_dir / f"{gene}_genotype_refs"
        make_database(GT_REFERENCES[gene], genotype_db)
        genotype_hits = blast(
            input_fasta,
            genotype_db,
            work_dir / f"{gene}_genotype.tsv",
            threads,
            f"assigning {gene} genotypes",
        )
        genotype_for_subject = lambda subject: (
            match.group(1)
            if (match := re.search(rf"HCV([1-8]){gene}", subject))
            else ""
        )
        genotype_labels = sorted(
            {
                label
                for header in fasta_header_records(GT_REFERENCES[gene])
                if (label := genotype_for_subject(header))
            }
        )
        genotype_labels_by_gene[gene] = genotype_labels
        genotype_rankings = ranked_hits(
            genotype_hits, genotype_for_subject, min_aligned_nt
        )
        best_genotype = {
            accession: ranking[0]
            for accession, ranking in genotype_rankings.items()
            if ranking
        }

        subtype_by_genotype: dict[str, list[tuple[str, str]]] = defaultdict(list)
        for header, sequence in fasta_header_records(SUBTYPE_REFERENCES[gene]).items():
            match = re.search(r"subtype=([1-8][A-Za-z0-9]*)", header)
            if match and sequence:
                subtype_by_genotype[match.group(1)[0]].append((match.group(1), sequence))
        query_ids_by_genotype: dict[str, list[str]] = defaultdict(list)
        subtype_labels_by_gene[gene] = {
            genotype: sorted(subtype for subtype, _ in subtypes)
            for genotype, subtypes in subtype_by_genotype.items()
        }
        for accession, (genotype, _) in best_genotype.items():
            query_ids_by_genotype[genotype].append(accession)
        for genotype, accessions in query_ids_by_genotype.items():
            subtype_records = {
                f"S{index}": sequence
                for index, (_, sequence) in enumerate(subtype_by_genotype[genotype])
            }
            subtype_labels = {
                f"S{index}": subtype
                for index, (subtype, _) in enumerate(subtype_by_genotype[genotype])
            }
            subtype_fasta, subtype_db = (
                work_dir / f"{gene}_{genotype}_subtypes.fasta",
                work_dir / f"{gene}_{genotype}_subtypes",
            )
            write_fasta(subtype_fasta, subtype_records)
            make_database(subtype_fasta, subtype_db)
            query_fasta = work_dir / f"{gene}_{genotype}_queries.fasta"
            write_fasta(
                query_fasta, {accession: records[accession] for accession in accessions}
            )
            subtype_hits = blast(
                query_fasta,
                subtype_db,
                work_dir / f"{gene}_{genotype}_subtypes.tsv",
                threads,
                f"assigning {gene} genotype {genotype} subtypes",
            )
            subtype_rankings = ranked_hits(
                subtype_hits, lambda subject: subtype_labels[subject], min_aligned_nt
            )
            for accession in accessions:
                assigned_genotype, genotype_hit = best_genotype[accession]
                genotype_choices = [
                    (label, distance(hit)) for label, hit in genotype_rankings[accession]
                ]
                call = {
                    "genotype": assigned_genotype,
                    "genotype_pident": genotype_hit[5],
                    "genotype_aligned_nt": genotype_hit[2],
                    "subtype": "",
                    "subtype_pident": "",
                    "subtype_aligned_nt": "",
                    "genotype_distances": dict(genotype_choices),
                    "genotype_distances_choices": genotype_choices,
                    "subtype_distances": {},
                    "subtype_distances_choices": [],
                }
                subtype_choices = [
                    (label, distance(hit))
                    for label, hit in subtype_rankings.get(accession, [])
                ]
                call["subtype_distances"] = dict(subtype_choices)
                call["subtype_distances_choices"] = subtype_choices
                if subtype_choices:
                    subtype_hit = subtype_rankings[accession][0][1]
                    call.update(
                        subtype=subtype_choices[0][0],
                        subtype_pident=subtype_hit[5],
                        subtype_aligned_nt=subtype_hit[2],
                    )
                results[gene][accession] = call
    return results, genotype_labels_by_gene, subtype_labels_by_gene


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
    best: dict[str, dict[str, tuple[float, int, int, int, int]]] = {
        gene: {} for gene in GENES
    }
    for gene in GENES:
        database = work_dir / f"{gene}_genotype_refs"
        make_database(GT_REFERENCES[gene], database)
        result = work_dir / f"{gene}_coverage.tsv"
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
            f"mapping {gene} coverage",
        )
        for line in result.read_text(encoding="utf-8").splitlines():
            qid, sid, qstart, qend, sstart, send, bitscore, length = line.split("\t")
            if not re.search(rf"HCV[1-8]{gene}", sid):
                continue
            score = (float(bitscore), int(length))
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
        default=REPO_ROOT / "HCVData" / "nonComet-PerGene",
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
        assignments_by_gene, genotype_labels_by_gene, subtype_labels_by_gene = assignment_calls(
            input_fasta, temp_dir / "assignments", args.threads, args.min_aligned_nt
        )
        hits_by_gene = coverage_hits(input_fasta, temp_dir / "coverage", args.threads)
        fields = [
            "Accession",
            "ClosestGenotype",
            "ClosestGenotypePident",
            "ClosestGenotypeAlignedNT",
            "ClosestSubtype",
            "ClosestSubtypePident",
            "ClosestSubtypeAlignedNT",
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
                    target_start_nt = (start_aa - 1) * 3 + 1
                    target_end_nt = end_aa * 3
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
                            "ClosestGenotypeAlignedNT": assignment.get(
                                "genotype_aligned_nt", ""
                            ),
                            "ClosestSubtype": assignment.get("subtype", ""),
                            "ClosestSubtypePident": assignment.get(
                                "subtype_pident", ""
                            ),
                            "ClosestSubtypeAlignedNT": assignment.get(
                                "subtype_aligned_nt", ""
                            ),
                            "ReferenceOverlapAA": reference_overlap,
                            "FullyCover": fully_cover,
                        }
                    )
                show_accession_progress(completed, len(accessions))
            print(file=sys.stderr)
        for output in output_paths.values():
            print(output)
        for gene, calls in assignments_by_gene.items():
            genotype_report = args.output_dir / f"{gene}_Genotype_Distances.csv"
            write_choice_report(
                genotype_report,
                genotype_labels_by_gene[gene],
                "Genotype",
                calls,
                "genotype_distances",
            )
            print(genotype_report)
            for genotype in genotype_labels_by_gene[gene]:
                subtype_report = (
                    args.output_dir
                    / f"{gene}_Subtype_Distances_Genotype_{genotype}.csv"
                )
                calls_for_genotype = {
                    accession: call
                    for accession, call in calls.items()
                    if call.get("genotype") == genotype
                }
                write_choice_report(
                    subtype_report,
                    subtype_labels_by_gene[gene].get(genotype, []),
                    "Subtype",
                    calls_for_genotype,
                    "subtype_distances",
                    assigned_genotype=True,
                )
                print(subtype_report)
            workbook = args.output_dir / f"{gene}_Subtyping_Distances.xlsx"
            write_distance_workbook(
                workbook,
                genotype_labels_by_gene[gene],
                subtype_labels_by_gene[gene],
                calls,
            )
            print(workbook)


if __name__ == "__main__":
    main()
