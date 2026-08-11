#!/usr/bin/env python3
"""Update the Comet subtype-genome JSON to the ICTV January 2026 subtype set."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from Bio import SeqIO

ADDITIONS = {
    "4a": ("Y11604", "ICTV26"),
    "5b": ("PQ899568", "ICTV26"),
    "6b": ("D84262", "ICTV26"),
    "7a": ("EF108306", "ICTV26"),
}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=Path, required=True)
    parser.add_argument("--genbank-dir", type=Path, required=True)
    args = parser.parse_args()
    rows = json.loads(args.json.read_text(encoding="utf-8"))
    for row in rows:
        row["genotypeName"] = {"Genotype4w3": "Genotype4w", "Genotype6xa5": "Genotype6xa"}.get(
            row.get("genotypeName"), row.get("genotypeName")
        )
    existing = {str(row.get("accession")) for row in rows}
    for subtype, (accession, author_year) in ADDITIONS.items():
        if accession in existing:
            continue
        record = SeqIO.read(args.genbank_dir / f"{accession}.gb", "genbank")
        sequence = str(record.seq).upper()[3419:]
        rows.append(
            {
                "accession": accession,
                "authorYear": author_year,
                "country": "UNK",
                "firstNA": 3420,
                "lastNA": len(record.seq),
                "genotypeName": f"Genotype{subtype}",
                "year": 2026,
                "sequence": sequence,
            }
        )
    rows.sort(key=lambda row: (str(row["genotypeName"]), str(row["accession"])))
    args.json.write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
