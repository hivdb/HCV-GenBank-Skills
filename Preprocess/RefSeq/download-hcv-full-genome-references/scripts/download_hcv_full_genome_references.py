#!/usr/bin/env python3
"""Download complete HCV genotype and subtype reference genomes from NCBI."""

from __future__ import annotations

import argparse
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_GT_GENE_FASTA = REPO_ROOT / "HCVData" / "HCV_GT_RefSeqs.fasta"
DEFAULT_SUBTYPE_JSON = REPO_ROOT / "HCVData" / "HCV_Subtype_Refs_By_Genome_NA.json"
DEFAULT_GENOTYPE_OUTPUT_DIR = REPO_ROOT / "HCVData" / "Genotype-Ref"
DEFAULT_SUBTYPE_OUTPUT_DIR = REPO_ROOT / "HCVData" / "Subtype-Ref"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def fasta_records(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence)))
            header, sequence = line[1:].strip(), []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line).upper())
    if header:
        records.append((header, "".join(sequence)))
    return records


def genotype_accessions(path: Path) -> list[tuple[str, str]]:
    by_genotype: dict[str, str] = {}
    for header, _ in fasta_records(path):
        match = re.match(r"HCV([1-8])NS(?:3|5A|5B)\s*\|\s*(\S+)", header)
        if not match:
            continue
        genotype, accession = match.groups()
        existing = by_genotype.setdefault(genotype, accession)
        if existing != accession:
            raise ValueError(f"Conflicting genotype {genotype} reference accessions: {existing}, {accession}")
    missing = sorted(set("12345678") - set(by_genotype))
    if missing:
        raise ValueError(f"{path} has no genotype references for: {', '.join(missing)}")
    return [(genotype, by_genotype[genotype]) for genotype in sorted(by_genotype, key=int)]


def subtype_accessions(path: Path) -> list[tuple[str, str]]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    result: list[tuple[str, str]] = []
    seen: set[str] = set()
    for row in rows:
        accession = str(row.get("accession", "")).strip()
        subtype = re.sub(r"^Genotype", "", str(row.get("genotypeName", "")).strip())
        if not accession or not subtype or accession in seen:
            continue
        result.append((subtype, accession))
        seen.add(accession)
    if not result:
        raise ValueError(f"{path} has no subtype reference accessions")
    return result


def cached_or_download(accession: str, cache_dir: Path, force: bool) -> tuple[str, str]:
    cache_path = cache_dir / f"{accession}.fasta"
    if cache_path.is_file() and not force:
        records = fasta_records(cache_path)
        if len(records) == 1 and records[0][1]:
            return records[0]
    url = f"{NCBI_EFETCH}?db=nuccore&id={accession}&rettype=fasta&retmode=text"
    for attempt in range(3):
        try:
            with urlopen(url, timeout=60) as response:
                payload = response.read().decode("utf-8")
            cache_path.write_text(payload, encoding="utf-8")
            records = fasta_records(cache_path)
            if len(records) != 1 or not records[0][1]:
                raise ValueError(f"NCBI returned no FASTA sequence for {accession}")
            return records[0]
        except (HTTPError, URLError, TimeoutError) as error:
            if attempt == 2:
                raise RuntimeError(f"Unable to download {accession}: {error}") from error
            time.sleep(2 ** attempt)
    raise AssertionError("unreachable")


def write_combined_fasta(
    path: Path, label_name: str, references: list[tuple[str, str]], cache_dir: Path, force: bool
) -> int:
    records: list[tuple[str, str]] = []
    for label, accession in references:
        _, sequence = cached_or_download(accession, cache_dir, force)
        records.append((f"{label_name}={label}|accession={accession}", sequence))
    with path.open("w", encoding="utf-8") as handle:
        for header, sequence in records:
            handle.write(f">{header}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start:start + 80] + "\n")
    return len(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gt-gene-fasta", type=Path, default=DEFAULT_GT_GENE_FASTA)
    parser.add_argument("--subtype-json", type=Path, default=DEFAULT_SUBTYPE_JSON)
    parser.add_argument("--genotype-output-dir", type=Path, default=DEFAULT_GENOTYPE_OUTPUT_DIR)
    parser.add_argument("--subtype-output-dir", type=Path, default=DEFAULT_SUBTYPE_OUTPUT_DIR)
    parser.add_argument("--force", action="store_true", help="Redownload references even when cached.")
    args = parser.parse_args()

    genotype_output_dir = args.genotype_output_dir.resolve()
    subtype_output_dir = args.subtype_output_dir.resolve()
    genotype_cache_dir = genotype_output_dir / "genbank_records"
    subtype_cache_dir = subtype_output_dir / "genbank_records"
    genotype_output_dir.mkdir(parents=True, exist_ok=True)
    subtype_output_dir.mkdir(parents=True, exist_ok=True)
    genotype_cache_dir.mkdir(parents=True, exist_ok=True)
    subtype_cache_dir.mkdir(parents=True, exist_ok=True)
    gt_count = write_combined_fasta(
        genotype_output_dir / "HCV_GT_FullGenome_RefSeqs.fasta",
        "genotype",
        genotype_accessions(args.gt_gene_fasta),
        genotype_cache_dir,
        args.force,
    )
    subtype_count = write_combined_fasta(
        subtype_output_dir / "HCV_Subtype_FullGenome_Refs.fasta",
        "subtype",
        subtype_accessions(args.subtype_json),
        subtype_cache_dir,
        args.force,
    )
    def display_path(path: Path) -> Path:
        try:
            return path.relative_to(REPO_ROOT)
        except ValueError:
            return path

    print(f"genotype_reference_count={gt_count}")
    print(f"subtype_reference_count={subtype_count}")
    print(f"genotype_output={display_path(genotype_output_dir / 'HCV_GT_FullGenome_RefSeqs.fasta')}")
    print(f"subtype_output={display_path(subtype_output_dir / 'HCV_Subtype_FullGenome_Refs.fasta')}")


if __name__ == "__main__":
    main()
