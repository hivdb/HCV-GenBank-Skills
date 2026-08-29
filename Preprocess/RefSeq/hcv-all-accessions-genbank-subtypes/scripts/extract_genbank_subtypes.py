#!/usr/bin/env python3
"""Extract genotype and subtype annotations for FASTA accessions from GenBank."""

from __future__ import annotations

import argparse
import csv
import io
import re
from collections.abc import Iterable
from pathlib import Path

from Bio import SeqIO


ACCESSION_RE = re.compile(r"^ACCESSION\s+([^\s]+)")
GENOTYPE_RE = re.compile(r"\bgenotype\s*[:=]\s*([0-9]+(?:[a-z][a-z0-9-]*)?)\b", re.I)
SUBTYPE_RE = re.compile(r"\bsubtype\s*[:=]\s*([0-9]+[a-z][a-z0-9-]*)\b", re.I)
GENOTYPE_ONLY_RE = re.compile(r"^([0-9]+)$")
GENOTYPE_SUBTYPE_RE = re.compile(r"^([0-9]+)([a-z][a-z0-9-]*)$", re.I)


def normalize_accession(accession: str) -> str:
    """Return the accession without an optional dot-version suffix."""
    return accession.split(".", 1)[0].strip()


def fasta_accessions(path: Path) -> list[str]:
    seen: set[str] = set()
    accessions: list[str] = []
    for record in SeqIO.parse(path, "fasta"):
        accession = normalize_accession(record.id)
        if accession and accession not in seen:
            seen.add(accession)
            accessions.append(accession)
    return accessions


def selected_genbank_records(paths: Iterable[Path], wanted: set[str]):
    """Yield wanted flatfile records, avoiding full parsing of unrelated records."""
    for file_number, path in enumerate(paths, start=1):
        print(f"Scanning file {file_number}: {path.name}", flush=True)
        record_lines: list[str] = []
        selected = False
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("LOCUS"):
                    record_lines = [line]
                    selected = False
                    continue
                if not record_lines:
                    continue
                record_lines.append(line)
                match = ACCESSION_RE.match(line)
                if match:
                    selected = normalize_accession(match.group(1)) in wanted
                if line == "//\n" or line.rstrip("\r\n") == "//":
                    if selected:
                        yield path, "".join(record_lines)
                    record_lines = []
                    selected = False


def qualifier_values(record) -> list[str]:
    values: list[str] = []
    for feature in record.features:
        if feature.type == "source":
            for entries in feature.qualifiers.values():
                values.extend(str(value) for value in entries)
    return values


def split_annotation(value: str) -> tuple[str, str]:
    value = value.strip().lower()
    if match := GENOTYPE_SUBTYPE_RE.fullmatch(value):
        return match.group(1), value
    if GENOTYPE_ONLY_RE.fullmatch(value):
        return value, ""
    return "", ""


def genotype_subtype(record) -> tuple[str, str]:
    genotype = ""
    subtype = ""
    for value in qualifier_values(record):
        if not genotype:
            match = GENOTYPE_RE.search(value)
            if match:
                genotype, inferred_subtype = split_annotation(match.group(1))
                subtype = subtype or inferred_subtype
        if not subtype:
            match = SUBTYPE_RE.search(value)
            if match:
                parsed_genotype, subtype = split_annotation(match.group(1))
                genotype = genotype or parsed_genotype
    return genotype, subtype


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genbank-dir", type=Path)
    parser.add_argument("--output-csv", type=Path, required=True)
    parser.add_argument("--file-start", type=int, default=0, help="Zero-based archive file offset.")
    parser.add_argument("--file-stop", type=int, help="Exclusive archive file offset.")
    parser.add_argument("--matches-csv", type=Path, help="Write matched metadata only, for batch processing.")
    parser.add_argument("--metadata-input", type=Path, action="append", default=[], help="Previously written --matches-csv file; may be repeated.")
    args = parser.parse_args()

    accessions = fasta_accessions(args.fasta)
    wanted = set(accessions)
    print(f"Loaded {len(accessions):,} unique FASTA accessions.", flush=True)
    metadata: dict[str, tuple[str, str]] = {}
    parsed_records = 0
    if args.metadata_input:
        for input_csv in args.metadata_input:
            with input_csv.open(newline="", encoding="utf-8") as handle:
                for row in csv.DictReader(handle):
                    metadata.setdefault(row["accession"], (row["genotype"], row["subtype"]))
    else:
        if args.genbank_dir is None:
            parser.error("--genbank-dir is required unless --metadata-input is supplied")
        archive_files = sorted(args.genbank_dir.rglob("*.seq"))[args.file_start : args.file_stop]
        for path, flatfile in selected_genbank_records(archive_files, wanted):
            try:
                record = SeqIO.read(io.StringIO(flatfile), "genbank")
            except Exception as error:  # Keep processing a large, imperfect archive.
                print(f"Warning: could not parse record in {path}: {error}")
                continue
            parsed_records += 1
            accession = normalize_accession(record.id)
            if accession in wanted and accession not in metadata:
                metadata[accession] = genotype_subtype(record)

    if args.matches_csv:
        args.matches_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.matches_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["accession", "genotype", "subtype"])
            writer.writerows((accession, *values) for accession, values in metadata.items())
        print(f"Matched records written: {len(metadata):,} -> {args.matches_csv}")
        return

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["accession", "genotype", "subtype"])
        writer.writerows((accession, *metadata.get(accession, ("", ""))) for accession in accessions)

    annotated = sum(bool(genotype or subtype) for genotype, subtype in metadata.values())
    print(f"FASTA accessions: {len(accessions):,}")
    print(f"Matched GenBank records: {len(metadata):,} (parsed {parsed_records:,})")
    print(f"Records with genotype or subtype: {annotated:,}")
    print(f"CSV rows written: {len(accessions):,} -> {args.output_csv}")


if __name__ == "__main__":
    main()
