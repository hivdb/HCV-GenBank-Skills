#!/usr/bin/env python3
"""Rebuild HCV subtype AA reference FASTAs from annotated GenBank proteins.

For each accession listed in the existing subtype-reference FASTAs, this script
uses an annotated gene CDS translation when available. Otherwise it extracts
the gene from the annotated polyprotein translation by local AA alignment to
the matching genotype AA reference. GenBank files are cached locally.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

from Bio import Align, SeqIO


EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
GENES = ("NS3", "NS5A_NTD", "NS5B")
MIN_COVERAGE = 0.8


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-dir", type=Path, default=Path("HCVData/Reference_seqs"))
    parser.add_argument("--hcv-fasta", type=Path, default=Path("HCVData/HCV-Ref-H77-Genotype1.fasta"))
    parser.add_argument(
        "--gt-reference-fasta",
        type=Path,
        default=Path("HCVData/Reference_seqs/HCV_GT_Refs_NS3_NS5A_NTD_NS5B_AA.fasta"),
        help="Per-gene genotype AA references used to extract genes from polyprotein CDSs.",
    )
    parser.add_argument("--genbank-dir", type=Path, default=Path("HCVData/Reference_seqs/genbank_records"))
    parser.add_argument("--email", default="", help="Optional NCBI contact email")
    parser.add_argument("--tool", default="hcv-subtype-reference-rebuild")
    parser.add_argument("--refresh-genbank", action="store_true")
    parser.add_argument("--gene", choices=GENES, action="append", help="Process only this gene; repeat as needed.")
    parser.add_argument("--report-xlsx", type=Path, default=Path("outputs/reference_seqs/HCV_Subtype_Reference_Rebuild_QC.xlsx"))
    return parser.parse_args()


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header: str | None = None
    chunks: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(chunks).upper()))
            header, chunks = line[1:].strip(), []
        elif line.strip():
            chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks).upper()))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def header_fields(header: str) -> dict[str, str]:
    return dict(item.split("=", 1) for item in header.split("|") if "=" in item)


def fetch_genbank(accession: str, email: str, tool: str) -> str:
    query = {"db": "nuccore", "id": accession, "rettype": "gbwithparts", "retmode": "text", "tool": tool}
    if email:
        query["email"] = email
    try:
        with urlopen(f"{EFETCH_URL}?{urlencode(query)}") as response:
            return response.read().decode("utf-8")
    except (HTTPError, URLError) as exc:
        raise RuntimeError(f"Unable to download GenBank record {accession}: {exc}") from exc


def get_record(accession: str, cache_dir: Path, email: str, tool: str, refresh: bool):
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / f"{accession}.gb"
    if refresh or not path.exists():
        path.write_text(fetch_genbank(accession, email, tool), encoding="utf-8")
    return SeqIO.read(path, "genbank")


def hcv_gene_references(path: Path) -> dict[str, str]:
    refs = {record.id.upper(): str(record.seq).upper() for record in SeqIO.parse(path, "fasta")}
    missing = [gene for gene in ("NS3", "NS5A", "NS5B") if gene not in refs]
    if missing:
        raise RuntimeError(f"Missing H77 genotype-1 FASTA gene references: {', '.join(missing)}")
    refs["NS5A_NTD"] = refs["NS5A"][:213]
    return refs


def genotype_gene_references(path: Path) -> dict[tuple[str, str], str]:
    references: dict[tuple[str, str], str] = {}
    for header, sequence in read_fasta(path):
        match = re.fullmatch(r"genotype ([1-8]) \| (NS3|NS5A_NTD|NS5B)", header)
        if match:
            references[(match.group(1), match.group(2))] = sequence
    missing = [f"GT{genotype} {gene}" for genotype in "12345678" for gene in GENES if (genotype, gene) not in references]
    if missing:
        raise RuntimeError(f"Missing genotype AA references: {', '.join(missing)}")
    return references


def feature_translation(feature) -> str:
    translation = feature.qualifiers.get("translation", [])
    return str(translation[0]).replace(" ", "").upper() if translation else ""


def feature_labels(feature) -> str:
    values: list[str] = []
    for key in ("gene", "product", "note", "region_name"):
        values.extend(str(value) for value in feature.qualifiers.get(key, []))
    return " ".join(values).upper()


def annotated_gene_translation(record, gene: str) -> str | None:
    aliases = ("NS5A",) if gene == "NS5A_NTD" else (gene,)
    for feature in record.features:
        if feature.type != "CDS":
            continue
        translation = feature_translation(feature)
        labels = feature_labels(feature)
        # A polyprotein product lists every gene (including NS5A) in its
        # description. It is not an individual gene CDS and must be aligned
        # as a whole below, rather than truncated from its N terminus.
        if translation and "POLYPROTEIN" not in labels and any(alias in labels for alias in aliases):
            return translation[:213] if gene == "NS5A_NTD" else translation
    return None


def polyprotein_translation(record) -> str:
    translations = [feature_translation(feature) for feature in record.features if feature.type == "CDS"]
    translations = [translation for translation in translations if translation]
    if not translations:
        raise RuntimeError(f"{record.id} has no annotated CDS translation")
    return max(translations, key=len)


def aligned_query_span(reference: str, query: str) -> tuple[str, float, float]:
    aligner = Align.PairwiseAligner(mode="local")
    aligner.match_score, aligner.mismatch_score = 2.0, -1.0
    aligner.open_gap_score, aligner.extend_gap_score = -5.0, -1.0
    alignment = aligner.align(reference, query)[0]
    coordinates = alignment.coordinates
    ref_start, ref_end = int(coordinates[0][0]), int(coordinates[0][-1])
    query_start, query_end = int(coordinates[1][0]), int(coordinates[1][-1])
    coverage = (ref_end - ref_start) / len(reference)
    identity = sum(a == b for a, b in zip(str(alignment[0]), str(alignment[1])) if a != "-" and b != "-")
    comparable = sum(a != "-" and b != "-" for a, b in zip(str(alignment[0]), str(alignment[1])))
    return query[query_start:query_end], coverage, identity / comparable if comparable else 0.0


def extract_gene(record, gene: str, reference: str) -> tuple[str, str, float, float]:
    annotated = annotated_gene_translation(record, gene)
    if annotated:
        sequence, coverage, identity = aligned_query_span(reference, annotated)
        return sequence, "annotated_gene_CDS", coverage, identity
    sequence, coverage, identity = aligned_query_span(reference, polyprotein_translation(record))
    return sequence, "annotated_polyprotein_alignment", coverage, identity


def write_report(path: Path, rows: list[dict[str, object]]) -> None:
    from openpyxl import Workbook

    workbook = Workbook()
    worksheet = workbook.active
    worksheet.title = "Rebuild_QC"
    fields = ["gene", "genotype", "subtype", "accession", "source", "coverage", "identity", "aa_length", "stop_codon_count", "status", "message"]
    worksheet.append(fields)
    for row in rows:
        worksheet.append([row.get(field, "") for field in fields])
    path.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(path)


def main() -> int:
    args = parse_args()
    genes = tuple(args.gene or GENES)
    genotype_references = genotype_gene_references(args.gt_reference_fasta)
    input_records = {gene: read_fasta(args.reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta") for gene in genes}
    accession_cache: dict[str, object] = {}
    rebuilt: dict[str, list[tuple[str, str]]] = {gene: [] for gene in genes}
    report_rows: list[dict[str, object]] = []
    total = sum(len(records) for records in input_records.values())
    completed = 0
    for gene in genes:
        for header, old_sequence in input_records[gene]:
            completed += 1
            fields = header_fields(header)
            accession = fields.get("accession", "")
            print(f"[{completed}/{total}] {gene} {accession}", file=sys.stderr, flush=True)
            row: dict[str, object] = {"gene": gene, "genotype": fields.get("genotype", ""), "subtype": fields.get("subtype", ""), "accession": accession}
            try:
                if accession not in accession_cache:
                    accession_cache[accession] = get_record(accession, args.genbank_dir, args.email, args.tool, args.refresh_genbank)
                reference = genotype_references[(fields["genotype"], gene)]
                sequence, source, coverage, identity = extract_gene(accession_cache[accession], gene, reference)
                if coverage < MIN_COVERAGE:
                    raise RuntimeError(f"low alignment coverage ({coverage:.1%})")
                rebuilt[gene].append((header, sequence))
                row.update(source=source, coverage=coverage, identity=identity, aa_length=len(sequence), stop_codon_count=sequence.count("*"), status="rebuilt", message="")
            except Exception as exc:
                rebuilt[gene].append((header, old_sequence))
                row.update(source="", coverage="", identity="", aa_length=len(old_sequence), stop_codon_count=old_sequence.count("*"), status="retained_original", message=str(exc))
            report_rows.append(row)
    for gene in genes:
        write_fasta(args.reference_dir / f"HCV_Subtype_Refs_{gene}_AA.fasta", rebuilt[gene])
    write_report(args.report_xlsx, report_rows)
    print(f"Rebuilt FASTAs in {args.reference_dir}; QC report: {args.report_xlsx}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
