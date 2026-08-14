#!/usr/bin/env python3
"""Filter combined HCV FASTAs to accessions used by profile-build workbooks."""

from __future__ import annotations

import argparse
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def cell_value(cell: ET.Element, shared_strings: list[str]) -> str:
    cell_type = cell.attrib.get("t")
    if cell_type == "inlineStr":
        return "".join(node.text or "" for node in cell.findall(".//main:t", NS))
    value = cell.find("main:v", NS)
    text = "" if value is None else (value.text or "")
    return shared_strings[int(text)] if cell_type == "s" and text else text


def profile_accessions(workbook_path: Path) -> set[str]:
    with zipfile.ZipFile(workbook_path) as archive:
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(node.text or "" for node in item.findall(".//main:t", NS))
                for item in root.findall("main:si", NS)
            ]
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        target_by_id = {
            item.attrib["Id"]: item.attrib["Target"].lstrip("/")
            for item in relationships
        }
        sheet = next(iter(workbook.find("main:sheets", NS)))
        sheet_path = target_by_id[sheet.attrib[f"{{{NS['rel']}}}id"]]
        sheet_root = ET.fromstring(archive.read(sheet_path))

    rows: list[dict[str, str]] = []
    for row in sheet_root.findall(".//main:sheetData/main:row", NS):
        rows.append(
            {
                cell.attrib["r"].rstrip("0123456789"): cell_value(cell, shared_strings)
                for cell in row.findall("main:c", NS)
            }
        )
    columns = {name: column for column, name in rows[0].items()}
    accession_column = columns["AccessionID"]
    aa_column = columns["AASequence"]
    return {
        row[accession_column]
        for row in rows[1:]
        if row.get(accession_column) and row.get(aa_column)
    }


def filter_fasta(fasta_path: Path, accessions: set[str]) -> tuple[int, int]:
    kept: list[str] = []
    record: list[str] = []
    original_count = 0
    kept_count = 0

    def flush() -> None:
        nonlocal kept_count
        if not record:
            return
        accession = record[0][1:].split(maxsplit=1)[0]
        if accession in accessions:
            kept.extend(record)
            kept_count += 1

    for line in fasta_path.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(">"):
            flush()
            record = [line]
            original_count += 1
        else:
            record.append(line)
    flush()
    fasta_path.write_text("".join(kept), encoding="utf-8")
    return original_count, kept_count


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--profile-workbook", type=Path, required=True)
    args = parser.parse_args()
    total, kept = filter_fasta(args.fasta, profile_accessions(args.profile_workbook))
    print(f"{args.fasta.name}: retained={kept} removed={total - kept}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
