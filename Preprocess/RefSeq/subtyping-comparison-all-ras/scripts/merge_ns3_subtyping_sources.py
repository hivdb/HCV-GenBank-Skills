#!/usr/bin/env python3
"""Outer-merge NS3 RAS-overlap genotype/subtype calls from five sources."""

from __future__ import annotations

import csv
import re
import argparse
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
OUTPUT_ROOT = REPO_ROOT / "HCVData" / "subtyping-comparison-all-ras"
GENES = ("NS3", "NS5A", "NS5B")
RAS_STEP_DIRS = {
    "NS3": "01-filter-ns3-ras-overlap",
    "NS5A": "02-filter-ns5a-ras-overlap",
    "NS5B": "03-filter-ns5b-ras-overlap",
}
MERGE_OUTPUT_DIR = OUTPUT_ROOT / "08-merge-ns3-subtyping-sources"
GENOTYPE_RE = re.compile(r"^(\d+)")


def accession_key(value: str) -> str:
    return value.strip().split("|", 1)[0].split(".", 1)[0].upper()


def comet_genotype(subtype: str) -> str:
    match = GENOTYPE_RE.match(subtype.strip())
    return match.group(1) if match else ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene", choices=GENES, default="NS3")
    return parser.parse_args()


def input_paths(gene: str) -> dict[str, Path]:
    return {
        "PerGene": OUTPUT_ROOT
        / RAS_STEP_DIRS[gene]
        / f"{gene}_AllSeq_NonComet_Coverage_RAS_Overlap.csv",
        "FullSeq": OUTPUT_ROOT
        / "04-filter-full-genome-by-ras-accessions"
        / f"{gene}_AllSeq_NonComet_Coverage_RAS_Overlap.csv",
        "CometPerGene": OUTPUT_ROOT
        / "05-filter-comet-pergene-by-ras-accessions"
        / f"{gene}_RAS_Overlap.csv",
        "CometFullSeq": OUTPUT_ROOT
        / "06-filter-comet-full-genome-by-ras-accessions"
        / f"{gene}_Comet_FullGenome_RAS_Overlap.csv",
        "GenBank": OUTPUT_ROOT
        / "07-filter-genbank-subtypes-by-ras-accessions"
        / f"{gene}_Genbank_Subtypes_RAS_Overlap.csv",
    }


def add_record(
    records: dict[str, tuple[str, str, str]],
    accession: str,
    genotype: str,
    subtype: str,
    subtype_aligned_nt: str,
    source: Path,
) -> None:
    record = (genotype, subtype, subtype_aligned_nt)
    if accession in records and records[accession] != record:
        raise ValueError(f"Conflicting calls for {accession} in {source}")
    records[accession] = record


def read_blast(path: Path) -> dict[str, tuple[str, str, str]]:
    records: dict[str, tuple[str, str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        required = {"Accession", "ClosestGenotype", "ClosestSubtype"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain {', '.join(sorted(required))}")
        for row in reader:
            accession = accession_key(row.get("Accession", ""))
            if accession:
                add_record(
                    records,
                    accession,
                    row.get("ClosestGenotype", "").strip(),
                    row.get("ClosestSubtype", "").strip(),
                    row.get("ClosestSubtypeAlignedNT", "").strip(),
                    path,
                )
    return records


def read_comet(path: Path) -> dict[str, tuple[str, str, str]]:
    records: dict[str, tuple[str, str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        required = {"name", "subtype"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain {', '.join(sorted(required))}")
        for row in reader:
            accession = accession_key(row.get("name", ""))
            subtype = row.get("subtype", "").strip()
            if accession:
                add_record(records, accession, comet_genotype(subtype), subtype, "", path)
    return records


def read_genbank(path: Path) -> dict[str, tuple[str, str, str]]:
    records: dict[str, tuple[str, str, str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as source:
        reader = csv.DictReader(source)
        required = {"accession", "genotype", "subtype"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain {', '.join(sorted(required))}")
        for row in reader:
            accession = accession_key(row.get("accession", ""))
            if accession:
                add_record(
                    records,
                    accession,
                    row.get("genotype", "").strip(),
                    row.get("subtype", "").strip(),
                    "",
                    path,
                )
    return records


def main() -> None:
    args = parse_args()
    inputs = input_paths(args.gene)
    calls = {
        "PerGene": read_blast(inputs["PerGene"]),
        "FullSeq": read_blast(inputs["FullSeq"]),
        "CometPerGene": read_comet(inputs["CometPerGene"]),
        "CometFullSeq": read_comet(inputs["CometFullSeq"]),
        "GenBank": read_genbank(inputs["GenBank"]),
    }
    accession_order = list(calls["PerGene"])
    known_accessions = set(accession_order)
    for source_calls in calls.values():
        accession_order.extend(sorted(set(source_calls) - known_accessions))
        known_accessions.update(source_calls)

    fields = ["Accession"]
    for source in calls:
        fields.extend((f"{source}Genotype", f"{source}Subtype"))
        if source in {"PerGene", "FullSeq"}:
            fields.append(f"{source}ClosestSubtypeAlignedNT")
    output_csv = MERGE_OUTPUT_DIR / f"{args.gene}_Subtyping_Sources_Merged.csv"
    missing_csv = MERGE_OUTPUT_DIR / f"{args.gene}_Missing_Accessions_By_Source.csv"
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    missing_rows = [
        {
            "Accession": accession,
            "MissingSources": ";".join(source for source, source_calls in calls.items() if accession not in source_calls),
        }
        for accession in accession_order
        if any(accession not in source_calls for source_calls in calls.values())
    ]
    with missing_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=["Accession", "MissingSources"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(missing_rows)
    with output_csv.open("w", encoding="utf-8", newline="") as destination:
        writer = csv.DictWriter(destination, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for accession in accession_order:
            row = {"Accession": accession}
            for source, source_calls in calls.items():
                genotype, subtype, subtype_aligned_nt = source_calls.get(
                    accession, ("", "", "")
                )
                row[f"{source}Genotype"] = genotype
                row[f"{source}Subtype"] = subtype
                if source in {"PerGene", "FullSeq"}:
                    row[f"{source}ClosestSubtypeAlignedNT"] = subtype_aligned_nt
            writer.writerow(row)

    print(f"Merged {args.gene} accessions: {len(accession_order):,}")
    for source, source_calls in calls.items():
        print(f"{source} missing accessions: {len(known_accessions - set(source_calls)):,}")
    print(f"Accessions missing from any source: {len(missing_rows):,}")
    print(f"Missing-source report: {missing_csv}")
    print(f"Output: {output_csv}")


if __name__ == "__main__":
    main()
