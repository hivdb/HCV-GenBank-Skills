#!/usr/bin/env python3
"""Find likely duplicate RefID records from a workbook sheet and Accessions CSV.

The candidate criteria mirror HCV-blasthit's search_same_ref.py: same
first-author surname, representative worksheet year within one year, shared
six-character accession prefix, and a PMID, title, or journal match. Records
with the same non-empty PMID are also direct candidates on that evidence alone.
"""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from difflib import SequenceMatcher
from pathlib import Path
from typing import Iterable
from xml.etree import ElementTree
from zipfile import ZIP_DEFLATED, ZipFile
from xml.sax.saxutils import escape


REPO_ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = REPO_ROOT / "HCVData" / "HCV-all-seq-subtype"
DEFAULT_WORKBOOK = REPO_ROOT / "HCVData" / "HCV_BlastHists_202604_data.xlsx"
DEFAULT_REF_CSV = DATA_DIR / "Ref.csv"
OUTPUT_DIR = REPO_ROOT / "outputs" / "Ref_same"
DEFAULT_REF_WITH_REFNAMES_OUTPUT = OUTPUT_DIR / "01_ref_with_refname" / "Ref_with_RefName.csv"
DEFAULT_MERGED_ORIGINAL_PMID_DIR = OUTPUT_DIR / "02_original_pmid_merge"
DEFAULT_FOUND_PMID_OUTPUT = OUTPUT_DIR / "03_found_pmid" / "Found_PMID.csv"
DEFAULT_ORIGINAL_SHEETS_REPLACED_PMID_DIR = OUTPUT_DIR / "04_original_sheets_pmid_replaced"
DEFAULT_REFNAME_SUMMARY = OUTPUT_DIR / "05_refname_counts" / "RefName_duplicate_counts.csv"
DEFAULT_OUTPUT = OUTPUT_DIR / "06_same_ref_candidates" / "same_ref_candidates.csv"
DEFAULT_GROUPS_OUTPUT = OUTPUT_DIR / "07_same_ref_candidate_groups" / "same_ref_candidate_groups.csv"
DEFAULT_GROUP_REFNAMES_OUTPUT = (
    OUTPUT_DIR / "08_same_ref_candidate_group_refnames" / "same_ref_candidate_groups_refnames.csv"
)
DEFAULT_GROUPKEY_REF_WITH_REFNAME_OUTPUT = (
    OUTPUT_DIR / "09_ref_with_refname_groupkey_deduplication" / "Ref_with_RefName_groupkey_deduplicated.csv"
)
DEFAULT_ORIGINAL_SHEETS_GROUPKEY_DEDUPLICATION_DIR = (
    OUTPUT_DIR / "10_original_sheets_groupkey_deduplication"
)
DEFAULT_DEDUPLICATED_REFNAME_SUMMARY = (
    OUTPUT_DIR / "11_ref_with_refname_refname_counts" / "RefName_duplicate_counts.csv"
)
DEFAULT_DUPLICATE_REFNAME_ROWS_DIR = OUTPUT_DIR / "12_duplicate_refname_rows"
PMID_SHEETS = ("Original", "Original_NS5A", "Original_NS3", "Original_NS5B")
SPREADSHEET_NS = {"x": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def normalize_text(value: str | None) -> str:
    """Lowercase text after removing punctuation and collapsing whitespace."""
    cleaned = re.sub(r"\s+", " ", (value or "").strip().lower())
    return re.sub(r"[^a-z0-9 ]", "", cleaned)


def normalize_pmid(value: str | None) -> str:
    """Return a PMID for exact matching, retaining all stored characters."""
    return (value or "").strip()


def author_key(author: str | None) -> str:
    normalized = normalize_text(author)
    if not normalized:
        return ""
    first_author = re.split(r"\bet al\b|;|,| and ", normalized, maxsplit=1)[0].strip()
    tokens = first_author.split()
    return tokens[-1] if tokens else ""


def parse_year(value: str | None) -> int | None:
    match = re.match(r"^(\d{4})", str(value or "").strip())
    return int(match.group(1)) if match else None


def refname_year(ref_name: str) -> str:
    """Extract the four-digit publication year embedded in a RefName."""
    match = re.search(r"\b(?:19|20)\d{2}\b", ref_name)
    return match.group(0) if match else ""


def accession_prefix(accession: str, length: int = 6) -> str:
    token = re.sub(r"[^A-Za-z0-9]", "", accession.strip().upper())
    return token[:length] if len(token) >= length else token


def similar(left: str, right: str, threshold: float) -> tuple[bool, float]:
    if not left or not right:
        return False, 0.0
    ratio = SequenceMatcher(None, left, right).ratio()
    return ratio >= threshold, ratio


def journal_parts(value: str | None) -> tuple[str, int | None]:
    text = (value or "").strip()
    if not text:
        return "", None
    match = re.search(r"\b(19|20)\d{2}\b", text)
    year = int(match.group(0)) if match else None
    remaining = re.sub(r"\b(19|20)\d{2}\b", " ", text)
    remaining = re.sub(
        r"\b(?:jan|january|feb|february|mar|march|apr|april|may|jun|june|jul|july|aug|august|sep|sept|september|oct|october|nov|november|dec|december)\b",
        " ",
        remaining,
        flags=re.IGNORECASE,
    )
    remaining = re.sub(r"\b\d{1,2}\b", " ", remaining)
    remaining = re.sub(r"[();:,./-]+", " ", remaining)
    return normalize_text(remaining), year


@dataclass
class RefRecord:
    ref_id: int
    author: str
    title: str
    journal: str
    pmid: str
    ref_name: str
    ref_name_key: str
    author_key: str
    title_norm: str
    journal_norm: str
    journal_year: int | None
    years: set[int] = field(default_factory=set)
    accessions: set[str] = field(default_factory=set)
    prefixes6: set[str] = field(default_factory=set)

    @property
    def rep_year(self) -> int | None:
        return min(self.years) if self.years else None


@dataclass(frozen=True)
class Candidate:
    left_id: int
    right_id: int
    reason: str
    shared_prefixes: str
    pmid_match: bool
    title_similarity: float
    journal_similarity: float


def require_columns(path: Path, columns: Iterable[str], fieldnames: list[str] | None) -> None:
    missing = sorted(set(columns) - set(fieldnames or []))
    if missing:
        raise ValueError(f"{path} is missing required column(s): {', '.join(missing)}")


def column_index(cell_reference: str) -> int:
    """Convert an Excel cell reference such as AB12 to its zero-based column index."""
    letters = re.match(r"[A-Z]+", cell_reference)
    if letters is None:
        raise ValueError(f"Invalid Excel cell reference: {cell_reference!r}")
    index = 0
    for letter in letters.group(0):
        index = index * 26 + ord(letter) - ord("A") + 1
    return index - 1


def read_shared_strings(archive: ZipFile) -> list[str]:
    if "xl/sharedStrings.xml" not in archive.namelist():
        return []
    root = ElementTree.fromstring(archive.read("xl/sharedStrings.xml"))
    return ["".join(node.itertext()) for node in root.findall("x:si", SPREADSHEET_NS)]


def worksheet_path(archive: ZipFile, sheet_name: str) -> str:
    workbook = ElementTree.fromstring(archive.read("xl/workbook.xml"))
    sheet = next(
        (node for node in workbook.findall("x:sheets/x:sheet", SPREADSHEET_NS)
         if node.get("name", "").casefold() == sheet_name.casefold()),
        None,
    )
    if sheet is None:
        available = ", ".join(node.get("name", "") for node in workbook.findall("x:sheets/x:sheet", SPREADSHEET_NS))
        raise ValueError(f"Workbook has no sheet named {sheet_name!r}; available sheets: {available}")
    relation_id = sheet.get("{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id")
    relationships = ElementTree.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
    relation = next((node for node in relationships if node.get("Id") == relation_id), None)
    if relation is None:
        raise ValueError(f"Workbook relationship missing for sheet {sheet_name!r}")
    return "xl/" + relation.get("Target", "").lstrip("/").removeprefix("xl/")


def read_worksheet_rows(workbook_path: Path, sheet_name: str) -> Iterable[dict[str, str]]:
    """Yield non-empty worksheet rows using only the Python standard library."""
    with ZipFile(workbook_path) as archive:
        strings = read_shared_strings(archive)
        root = ElementTree.fromstring(archive.read(worksheet_path(archive, sheet_name)))
        rows = root.findall(".//x:sheetData/x:row", SPREADSHEET_NS)
        if not rows:
            return
        values_by_row: list[list[str]] = []
        for row in rows:
            values: list[str] = []
            for cell in row.findall("x:c", SPREADSHEET_NS):
                index = column_index(cell.get("r", ""))
                values.extend([""] * (index + 1 - len(values)))
                value = cell.findtext("x:v", default="", namespaces=SPREADSHEET_NS)
                if cell.get("t") == "s" and value:
                    value = strings[int(value)]
                elif cell.get("t") == "inlineStr":
                    value = "".join(cell.itertext())
                values[index] = value
            values_by_row.append(values)
        headers = values_by_row[0]
        for values in values_by_row[1:]:
            yield {header: values[index] if index < len(values) else "" for index, header in enumerate(headers) if header}


def read_ref_records(input_csv: Path) -> dict[int, RefRecord]:
    records: dict[int, RefRecord] = {}
    required = ("RefID", "RefName", "Year", "PMID", "Author", "Title", "Journal")
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(input_csv, required, reader.fieldnames)
        for row in reader:
            raw_ref_id = (row.get("RefID") or "").strip()
            if not raw_ref_id:
                continue
            try:
                ref_id = int(raw_ref_id)
            except ValueError as error:
                raise ValueError(f"{input_csv} has non-integer RefID: {raw_ref_id!r}") from error
            if ref_id in records:
                raise ValueError(f"{input_csv} contains duplicate RefID {ref_id}")
            journal = (row.get("Journal") or "").strip()
            journal_norm, journal_year = journal_parts(journal)
            author = (row.get("Author") or "").strip()
            title = (row.get("Title") or "").strip()
            ref_name = (row.get("RefName") or "").strip()
            records[ref_id] = RefRecord(
                ref_id=ref_id,
                author=author,
                title=title,
                journal=journal,
                pmid=normalize_pmid(row.get("PMID")),
                ref_name=ref_name,
                ref_name_key=normalize_text(ref_name),
                author_key=author_key(author),
                title_norm=normalize_text(title),
                journal_norm=journal_norm,
                journal_year=journal_year,
            )
            year = parse_year(row.get("Year"))
            if year is not None:
                records[ref_id].years.add(year)
    return records


def add_accession_data(records: dict[int, RefRecord], accessions_csv: Path) -> None:
    with accessions_csv.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        require_columns(accessions_csv, ("RefID", "Accession"), reader.fieldnames)
        for row in reader:
            raw_ref_id = (row.get("RefID") or "").strip()
            if not raw_ref_id:
                continue
            try:
                record = records.get(int(raw_ref_id))
            except ValueError as error:
                raise ValueError(f"{accessions_csv} has non-integer RefID: {raw_ref_id!r}") from error
            if record is None:
                continue
            accession = (row.get("Accession") or "").strip().upper()
            if accession:
                record.accessions.add(accession)
                prefix = accession_prefix(accession)
                if prefix:
                    record.prefixes6.add(prefix)


def journal_match(left: RefRecord, right: RefRecord) -> tuple[bool, float]:
    matches, ratio = similar(left.journal_norm, right.journal_norm, 0.90)
    if not matches:
        return False, ratio
    if left.journal_year is not None and right.journal_year is not None:
        return abs(left.journal_year - right.journal_year) <= 1, ratio
    return True, ratio


def compare(left: RefRecord, right: RefRecord) -> Candidate | None:
    if left.pmid and right.pmid and left.pmid != right.pmid:
        return None
    pmid_match = bool(left.pmid and left.pmid == right.pmid)
    title_match, title_ratio = similar(left.title_norm, right.title_norm, 0.92)
    journal_is_match, journal_ratio = journal_match(left, right)
    prefixes = sorted(left.prefixes6 & right.prefixes6)
    legacy_match = (
        bool(left.author_key and left.author_key == right.author_key)
        and left.rep_year is not None
        and right.rep_year is not None
        and abs(left.rep_year - right.rep_year) <= 1
        and bool(prefixes)
        and (pmid_match or title_match or journal_is_match)
    )
    same_refname = bool(left.ref_name_key and left.ref_name_key == right.ref_name_key)
    same_refname_title = same_refname and title_match
    same_refname_journal = (
        same_refname and bool(left.journal_norm) and left.journal_norm == right.journal_norm
    )
    # A shared non-empty PMID identifies the same publication even if the
    # remaining bibliographic metadata differ or are absent.
    if not (pmid_match or legacy_match or same_refname_title or same_refname_journal):
        return None
    reason = []
    if pmid_match:
        reason.append(f"same_pmid={left.pmid}")
    if legacy_match:
        reason.extend([
            f"author={left.author_key}",
            f"year={left.rep_year}/{right.rep_year}",
            f"prefix6={';'.join(prefixes[:5])}",
        ])
        if title_match:
            reason.append("title_similar")
        if journal_is_match:
            reason.append("journal_similar")
    if same_refname_title:
        reason.append("same_refname_title_similar")
    if same_refname_journal:
        reason.append("same_refname_journal_same")
    return Candidate(
        left.ref_id, right.ref_id, "|".join(reason), ";".join(prefixes), pmid_match,
        title_ratio, journal_ratio,
    )


def filter_candidates_by_pmid(
    records: dict[int, RefRecord], candidates: list[Candidate]
) -> list[Candidate]:
    """Keep candidate edges only when their connected group has at most one PMID."""
    parent: dict[int, int] = {}
    component_pmids: dict[int, set[str]] = {}

    def find(ref_id: int) -> int:
        parent.setdefault(ref_id, ref_id)
        component_pmids.setdefault(ref_id, {records[ref_id].pmid} if records[ref_id].pmid else set())
        if parent[ref_id] != ref_id:
            parent[ref_id] = find(parent[ref_id])
        return parent[ref_id]

    accepted: list[Candidate] = []
    for candidate in candidates:
        left_root = find(candidate.left_id)
        right_root = find(candidate.right_id)
        if left_root == right_root:
            accepted.append(candidate)
            continue
        if len(component_pmids[left_root] | component_pmids[right_root]) > 1:
            continue
        parent[right_root] = left_root
        component_pmids[left_root].update(component_pmids[right_root])
        accepted.append(candidate)
    return accepted


def find_candidates(records: dict[int, RefRecord]) -> list[Candidate]:
    blocks: dict[str, list[RefRecord]] = {}
    for record in records.values():
        if record.author_key:
            blocks.setdefault(f"author:{record.author_key}", []).append(record)
        if record.ref_name_key:
            blocks.setdefault(f"refname:{record.ref_name_key}", []).append(record)
        if record.pmid:
            blocks.setdefault(f"pmid:{record.pmid}", []).append(record)
    pair_ids: set[tuple[int, int]] = set()
    for block in blocks.values():
        block.sort(key=lambda record: record.ref_id)
        for index, left in enumerate(block):
            pair_ids.update((left.ref_id, right.ref_id) for right in block[index + 1:])
    candidates: list[Candidate] = []
    for left_id, right_id in sorted(pair_ids):
        candidate = compare(records[left_id], records[right_id])
        if candidate:
            candidates.append(candidate)
    return filter_candidates_by_pmid(records, candidates)


def component_keys(candidates: list[Candidate]) -> dict[int, int]:
    adjacency: dict[int, set[int]] = {}
    for candidate in candidates:
        adjacency.setdefault(candidate.left_id, set()).add(candidate.right_id)
        adjacency.setdefault(candidate.right_id, set()).add(candidate.left_id)
    keys: dict[int, int] = {}
    for ref_id in sorted(adjacency):
        if ref_id in keys:
            continue
        stack, component = [ref_id], []
        while stack:
            current = stack.pop()
            if current in keys:
                continue
            component.append(current)
            keys[current] = -1
            stack.extend(adjacency[current] - keys.keys())
        key = min(component)
        for member in component:
            keys[member] = key
    return keys


def write_candidates(records: dict[int, RefRecord], candidates: list[Candidate], output_csv: Path) -> None:
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    keys = component_keys(candidates)
    fields = [
        "GroupKeyRefID", "RefID1", "RefID2", "AuthorKey", "Year1", "Year2",
        "SharedPrefix6", "MatchingPMID", "TitleSimilarity", "JournalSimilarity",
        "RefName1", "RefName2", "Author1", "Author2", "PMID1", "PMID2", "Title1", "Title2",
        "Journal1", "Journal2", "AccessionCount1", "AccessionCount2", "Reason",
    ]
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for candidate in candidates:
            left, right = records[candidate.left_id], records[candidate.right_id]
            writer.writerow({
                "GroupKeyRefID": keys[candidate.left_id], "RefID1": left.ref_id,
                "RefID2": right.ref_id, "AuthorKey": left.author_key,
                "Year1": left.rep_year or "", "Year2": right.rep_year or "",
                "SharedPrefix6": candidate.shared_prefixes,
                "MatchingPMID": "Yes" if candidate.pmid_match else "",
                "TitleSimilarity": f"{candidate.title_similarity:.3f}",
                "JournalSimilarity": f"{candidate.journal_similarity:.3f}",
                "RefName1": left.ref_name, "RefName2": right.ref_name,
                "Author1": left.author, "Author2": right.author, "PMID1": left.pmid,
                "PMID2": right.pmid, "Title1": left.title, "Title2": right.title,
                "Journal1": left.journal, "Journal2": right.journal,
                "AccessionCount1": len(left.accessions), "AccessionCount2": len(right.accessions),
                "Reason": candidate.reason,
            })


def write_candidate_groups(candidates: list[Candidate], output_csv: Path) -> int:
    """Write connected candidate groups with their complete RefID membership."""
    groups: dict[int, list[int]] = {}
    for ref_id, group_key in component_keys(candidates).items():
        groups.setdefault(group_key, []).append(ref_id)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["GroupKeyRefID", "RefIDCount", "RefIDs"])
        writer.writerows(
            [group_key, len(ref_ids), "; ".join(str(ref_id) for ref_id in sorted(ref_ids))]
            for group_key, ref_ids in sorted(groups.items())
        )
    return len(groups)


def read_refname_by_ref_id(input_csv: Path) -> dict[int, str]:
    """Return non-empty RefNames from an enriched Ref CSV, indexed by RefID."""
    refnames: dict[int, str] = {}
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(input_csv, ("RefID", "RefName"), reader.fieldnames)
        for row in reader:
            raw_ref_id = (row.get("RefID") or "").strip()
            ref_name = (row.get("RefName") or "").strip()
            if raw_ref_id and ref_name:
                refnames[int(raw_ref_id)] = ref_name
    return refnames


def write_candidate_group_refnames(
    groups_csv: Path, ref_with_refnames_csv: Path, output_csv: Path
) -> int:
    """Add a deduplicated RefNameList and its count to every candidate-group row."""
    refnames = read_refname_by_ref_id(ref_with_refnames_csv)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with groups_csv.open(newline="", encoding="utf-8") as input_handle, output_csv.open(
        "w", newline="", encoding="utf-8"
    ) as output_handle:
        reader = csv.DictReader(input_handle)
        require_columns(groups_csv, ("GroupKeyRefID", "RefIDs"), reader.fieldnames)
        writer = csv.DictWriter(
            output_handle,
            fieldnames=["GroupKeyRefID", "RefIDs", "RefNameList", "RefNameListCount"],
            lineterminator="\n",
        )
        writer.writeheader()
        count = 0
        for row in reader:
            ref_ids = [int(value.strip()) for value in (row["RefIDs"] or "").split(";") if value.strip()]
            refname_list = list(dict.fromkeys(refnames[ref_id] for ref_id in ref_ids if ref_id in refnames))
            writer.writerow({
                "GroupKeyRefID": row["GroupKeyRefID"],
                "RefIDs": row["RefIDs"],
                "RefNameList": "; ".join(refname_list),
                "RefNameListCount": len(refname_list),
            })
            count += 1
    return count


def read_groupkey_by_ref_id(groups_csv: Path) -> dict[int, int]:
    """Map every RefID in a candidate group to that group's key RefID."""
    group_keys: dict[int, int] = {}
    with groups_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(groups_csv, ("GroupKeyRefID", "RefIDs"), reader.fieldnames)
        for row in reader:
            group_key = int((row["GroupKeyRefID"] or "").strip())
            for value in (row["RefIDs"] or "").split(";"):
                if not value.strip():
                    continue
                ref_id = int(value.strip())
                existing = group_keys.setdefault(ref_id, group_key)
                if existing != group_key:
                    raise ValueError(f"RefID {ref_id} belongs to multiple candidate groups")
    return group_keys


def write_groupkey_deduplicated_original_rows(
    groups_csv: Path, source_csv: Path, output_csv: Path
) -> tuple[int, int]:
    """Remove non-key group members while leaving retained enriched-CSV rows unchanged."""
    group_keys = read_groupkey_by_ref_id(groups_csv)
    with source_csv.open(newline="", encoding="utf-8") as input_handle:
        reader = csv.DictReader(input_handle)
        require_columns(source_csv, ("RefID",), reader.fieldnames)
        fieldnames = reader.fieldnames or []
        source_rows = list(reader)

    source_ref_ids = {int((row["RefID"] or "").strip()) for row in source_rows}
    missing_group_keys = sorted(set(group_keys.values()) - source_ref_ids)
    if missing_group_keys:
        raise ValueError(f"{source_csv} has no row for group key(s): {missing_group_keys}")

    non_key_group_members = {
        ref_id for ref_id, group_key in group_keys.items() if ref_id != group_key
    }
    retained_rows = [
        row for row in source_rows if int((row["RefID"] or "").strip()) not in non_key_group_members
    ]

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(retained_rows)
    return len(source_rows), len(retained_rows)


def write_ref_with_refnames(
    ref_csv: Path, workbook_path: Path, sheet_name: str, output_csv: Path
) -> int:
    """Append Original-sheet RefNames and their extracted years to Ref.csv."""
    refname_by_id = {
        int((row["RefID"] or "").strip()): (row.get("RefName") or "").strip()
        for row in read_worksheet_rows(workbook_path, sheet_name)
        if (row.get("RefID") or "").strip()
    }
    with ref_csv.open(newline="", encoding="utf-8-sig") as input_handle:
        reader = csv.DictReader(input_handle)
        require_columns(ref_csv, ("RefID", "MedlineID"), reader.fieldnames)
        fieldnames = reader.fieldnames or []
        if "PMID" in fieldnames:
            raise ValueError(f"{ref_csv} cannot contain both MedlineID and PMID")
        output_fields = ["PMID" if field == "MedlineID" else field for field in fieldnames]
        if "RefName" not in output_fields:
            output_fields.append("RefName")
        if "Year" not in output_fields:
            output_fields.append("Year")
        rows = list(reader)

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=output_fields, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            ref_id = int((row["RefID"] or "").strip())
            if ref_id not in refname_by_id:
                raise ValueError(f"{workbook_path} has no RefName row for RefID {ref_id}")
            output_row = dict(row)
            output_row["PMID"] = output_row.pop("MedlineID")
            output_row["RefName"] = refname_by_id[ref_id]
            output_row["Year"] = refname_year(output_row["RefName"])
            writer.writerow(output_row)
    return len(rows)


def write_refname_summary(input_csv: Path, output_csv: Path) -> int:
    """Write the RefName summary CSV."""
    groups: dict[str, list[str]] = {}
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(input_csv, ("RefID", "RefName"), reader.fieldnames)
        for row in reader:
            ref_name = (row.get("RefName") or "").strip()
            ref_id = (row.get("RefID") or "").strip()
            if ref_name and ref_id:
                groups.setdefault(ref_name, []).append(ref_id)
    refname_groups = sorted(
        ((name, sorted(ref_ids, key=int)) for name, ref_ids in groups.items()),
        key=lambda item: (-len(item[1]), item[0].casefold()),
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["RefName", "RefIDCount", "RefIDs"])
        writer.writerows(
            [name, len(ref_ids), "; ".join(ref_ids)]
            for name, ref_ids in refname_groups
        )
    return len(refname_groups)


def write_deduplicated_refname_summary(input_csv: Path, output_csv: Path) -> int:
    """Write the step-01 RefName summary for retained step-05 rows."""
    groups: dict[str, list[str]] = {}
    with input_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(input_csv, ("RefID", "RefName"), reader.fieldnames)
        for row in reader:
            ref_name = (row["RefName"] or "").strip()
            ref_id = (row["RefID"] or "").strip()
            if ref_name and ref_id:
                groups.setdefault(ref_name, []).append(ref_id)
    refname_groups = sorted(
        ((name, sorted(ref_ids, key=int)) for name, ref_ids in groups.items()),
        key=lambda item: (-len(item[1]), item[0].casefold()),
    )
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["RefName", "RefIDCount", "RefIDs"])
        writer.writerows(
            [name, len(ref_ids), "; ".join(ref_ids)]
            for name, ref_ids in refname_groups
        )
    return len(refname_groups)


def refname_filename(ref_name: str, ref_id_count: int) -> str:
    """Return a filesystem-safe filename that identifies a duplicate RefName group."""
    safe_name = re.sub(r"[^A-Za-z0-9._-]+", "_", ref_name).strip("._") or "unnamed"
    return f"{safe_name}__RefIDCount_{ref_id_count}.csv"


def write_duplicate_refname_rows(
    summary_csv: Path, original_csv: Path, output_dir: Path
) -> int:
    """Write one original-row CSV for each step-06 RefName group with two or more RefIDs."""
    with original_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(original_csv, ("RefID",), reader.fieldnames)
        fieldnames = reader.fieldnames or []
        original_rows = list(reader)

    output_dir.mkdir(parents=True, exist_ok=True)
    for stale_output in output_dir.rglob("*.csv"):
        stale_output.unlink()
    for stale_dir in sorted((path for path in output_dir.iterdir() if path.is_dir()), reverse=True):
        stale_dir.rmdir()

    used_filenames: set[str] = set()
    count = 0
    with summary_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(summary_csv, ("RefName", "RefIDCount", "RefIDs"), reader.fieldnames)
        for group in reader:
            ref_id_count = int((group["RefIDCount"] or "").strip())
            if ref_id_count < 2:
                continue
            ref_ids = {int(value.strip()) for value in (group["RefIDs"] or "").split(";") if value.strip()}
            if len(ref_ids) != ref_id_count:
                raise ValueError(f"{summary_csv} has an invalid RefIDCount for {group['RefName']!r}")
            filename = refname_filename(group["RefName"], ref_id_count)
            if filename in used_filenames:
                raise ValueError(f"Duplicate output filename for RefName group: {filename}")
            used_filenames.add(filename)
            rows = [row for row in original_rows if int((row["RefID"] or "").strip()) in ref_ids]
            if len(rows) != ref_id_count:
                raise ValueError(f"{original_csv} is missing RefID(s) for {group['RefName']!r}")
            count_dir = output_dir / f"RefIDCount_{ref_id_count}"
            count_dir.mkdir(exist_ok=True)
            with (count_dir / filename).open("w", newline="", encoding="utf-8") as output_handle:
                writer = csv.DictWriter(output_handle, fieldnames=fieldnames, lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
            count += 1
    return count


def write_merged_original_pmids(workbook_path: Path, output_dir: Path) -> tuple[int, int, int]:
    """Merge PMIDs from all Original sheets and write the requested Original CSV copies."""
    pmids_by_ref_id: dict[int, list[str]] = {}
    original_rows: list[dict[str, str]] | None = None
    fieldnames: list[str] | None = None
    for sheet_name in PMID_SHEETS:
        rows = list(read_worksheet_rows(workbook_path, sheet_name))
        if not rows:
            raise ValueError(f"{workbook_path} sheet {sheet_name!r} is empty")
        require_columns(workbook_path, ("RefID", "PMID"), list(rows[0]))
        if sheet_name == "Original":
            original_rows = rows
            fieldnames = list(rows[0])
        for row in rows:
            raw_ref_id = (row["RefID"] or "").strip()
            if not raw_ref_id:
                continue
            ref_id = int(raw_ref_id)
            pmid = (row["PMID"] or "").strip()
            if pmid and pmid not in pmids_by_ref_id.setdefault(ref_id, []):
                pmids_by_ref_id[ref_id].append(pmid)

    if original_rows is None or fieldnames is None:
        raise AssertionError("Original sheet was not loaded")
    merged_rows = []
    for row in original_rows:
        merged_row = dict(row)
        ref_id = int((row["RefID"] or "").strip())
        merged_row["PMID"] = "; ".join(pmids_by_ref_id.get(ref_id, []))
        merged_rows.append(merged_row)

    output_dir.mkdir(parents=True, exist_ok=True)

    def write_rows(filename: str, rows: list[dict[str, str]]) -> None:
        with (output_dir / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    with_pmid = [row for row in merged_rows if row["PMID"]]
    numeric_pmid = [row for row in with_pmid if row["PMID"].isdigit()]
    non_numeric_pmid = [row for row in with_pmid if not row["PMID"].isdigit()]
    write_rows("Original_merged_PMID.csv", merged_rows)
    write_rows("Original_with_PMID.csv", with_pmid)
    write_rows("Original_with_numeric_PMID.csv", numeric_pmid)
    write_rows("Original_with_non_numeric_PMID.csv", non_numeric_pmid)
    return len(merged_rows), len(numeric_pmid), len(non_numeric_pmid)


def write_found_pmids(step1_csv: Path, merged_pmid_csv: Path, output_csv: Path) -> int:
    """Fill blank step-01 PMIDs from the merged step-02 Original PMID values."""
    with merged_pmid_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(merged_pmid_csv, ("RefID", "PMID"), reader.fieldnames)
        pmids_by_ref_id = {
            int((row["RefID"] or "").strip()): (row["PMID"] or "").strip()
            for row in reader if (row.get("RefID") or "").strip() and (row.get("PMID") or "").strip()
        }
    with step1_csv.open(newline="", encoding="utf-8") as input_handle:
        reader = csv.DictReader(input_handle)
        require_columns(step1_csv, ("RefID", "PMID"), reader.fieldnames)
        fieldnames = reader.fieldnames or []
        source_rows = list(reader)

    replacements = 0
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in source_rows:
            output_row = dict(row)
            if not (output_row["PMID"] or "").strip():
                found_pmid = pmids_by_ref_id.get(int((output_row["RefID"] or "").strip()), "")
                if found_pmid:
                    output_row["PMID"] = found_pmid
                    replacements += 1
            writer.writerow(output_row)
    return replacements


def write_found_pmid_report(
    workbook_path: Path, original_ref_csv: Path, found_pmid_csv: Path, output_csv: Path
) -> dict[str, int]:
    """Count per-sheet RefIDs whose blank step-01 PMID was filled in step 03."""
    def read_pmids(path: Path) -> dict[int, str]:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            require_columns(path, ("RefID", "PMID"), reader.fieldnames)
            return {
                int((row["RefID"] or "").strip()): (row["PMID"] or "").strip()
                for row in reader
                if (row.get("RefID") or "").strip()
            }

    original_pmids_by_ref_id = read_pmids(original_ref_csv)
    found_pmids_by_ref_id = read_pmids(found_pmid_csv)

    found_by_sheet: dict[str, int] = {}
    for sheet_name in PMID_SHEETS:
        rows = list(read_worksheet_rows(workbook_path, sheet_name))
        if not rows:
            raise ValueError(f"{workbook_path} sheet {sheet_name!r} is empty")
        require_columns(workbook_path, ("RefID",), list(rows[0]))
        found = 0
        for row in rows:
            raw_ref_id = (row["RefID"] or "").strip()
            if not raw_ref_id:
                continue
            ref_id = int(raw_ref_id)
            if not original_pmids_by_ref_id.get(ref_id, "") and found_pmids_by_ref_id.get(ref_id, ""):
                found += 1
        found_by_sheet[sheet_name] = found

    output_csv.parent.mkdir(parents=True, exist_ok=True)
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["Sheet", "FoundPMIDCount"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {"Sheet": sheet_name, "FoundPMIDCount": count}
            for sheet_name, count in found_by_sheet.items()
        )
    return found_by_sheet


def write_original_sheets_with_replaced_pmids(
    workbook_path: Path, merged_pmid_csv: Path, output_dir: Path
) -> dict[str, int]:
    """Replace each original-sheet PMID from step 02's RefID-keyed PMID values."""
    with merged_pmid_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(merged_pmid_csv, ("RefID", "PMID"), reader.fieldnames)
        pmid_by_ref_id = {
            int((row["RefID"] or "").strip()): (row["PMID"] or "").strip()
            for row in reader
            if (row.get("RefID") or "").strip()
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    updates_by_sheet: dict[str, int] = {}
    for sheet_name in PMID_SHEETS:
        rows = list(read_worksheet_rows(workbook_path, sheet_name))
        if not rows:
            raise ValueError(f"{workbook_path} sheet {sheet_name!r} is empty")
        fieldnames = list(rows[0])
        require_columns(workbook_path, ("RefID", "PMID"), fieldnames)
        updates = 0
        output_rows = []
        for row in rows:
            output_row = dict(row)
            raw_ref_id = (output_row["RefID"] or "").strip()
            if raw_ref_id:
                replacement = pmid_by_ref_id.get(int(raw_ref_id))
                if replacement is not None:
                    original_pmid = (output_row["PMID"] or "").strip()
                    if original_pmid != replacement:
                        updates += 1
                    output_row["PMID"] = replacement
            output_rows.append(output_row)
        output_csv = output_dir / f"{sheet_name}_PMID_replaced.csv"
        with output_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(output_rows)
        updates_by_sheet[sheet_name] = updates
    return updates_by_sheet


def excel_column_name(index: int) -> str:
    """Return the one-based Excel column name for a zero-based index."""
    name = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(ord("A") + remainder) + name
    return name


def write_csv_tables_xlsx(
    tables: dict[str, tuple[list[str], list[dict[str, str]]]], output_xlsx: Path
) -> None:
    """Write CSV-shaped tables as an Excel workbook using inline-string cells."""
    output_xlsx.parent.mkdir(parents=True, exist_ok=True)
    content_types = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
        '<Default Extension="xml" ContentType="application/xml"/>',
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    ]
    workbook_sheets = []
    workbook_relationships = []
    with ZipFile(output_xlsx, "w", ZIP_DEFLATED) as archive:
        for index, (sheet_name, (fieldnames, rows)) in enumerate(tables.items(), start=1):
            content_types.append(
                f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
                'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
            )
            workbook_sheets.append(
                f'<sheet name="{escape(sheet_name)}" sheetId="{index}" r:id="rId{index}"/>'
            )
            workbook_relationships.append(
                f'<Relationship Id="rId{index}" '
                'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
                f'Target="worksheets/sheet{index}.xml"/>'
            )
            all_rows = [dict(zip(fieldnames, fieldnames)), *rows]
            xml_rows = []
            for row_index, row in enumerate(all_rows, start=1):
                cells = []
                for column_index, fieldname in enumerate(fieldnames):
                    value = escape(str(row.get(fieldname, "")))
                    cell_reference = f"{excel_column_name(column_index)}{row_index}"
                    cells.append(
                        f'<c r="{cell_reference}" t="inlineStr"><is><t xml:space="preserve">{value}</t></is></c>'
                    )
                xml_rows.append(f'<row r="{row_index}">{"".join(cells)}</row>')
            worksheet = (
                '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
                f'<sheetData>{"".join(xml_rows)}</sheetData></worksheet>'
            )
            archive.writestr(f"xl/worksheets/sheet{index}.xml", worksheet)
        content_types.append('</Types>')
        archive.writestr("[Content_Types].xml", "".join(content_types))
        archive.writestr(
            "_rels/.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
            'Target="xl/workbook.xml"/></Relationships>',
        )
        archive.writestr(
            "xl/workbook.xml",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f'<sheets>{"".join(workbook_sheets)}</sheets></workbook>',
        )
        archive.writestr(
            "xl/_rels/workbook.xml.rels",
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{"".join(workbook_relationships)}</Relationships>',
        )


def write_groupkey_deduplicated_original_sheets(
    groups_csv: Path, input_dir: Path, output_dir: Path, found_pmid_report_csv: Path
) -> dict[str, tuple[int, int]]:
    """Apply group-key deduplication to every step-04 Original-sheet CSV."""
    group_keys = read_groupkey_by_ref_id(groups_csv)
    output_dir.mkdir(parents=True, exist_ok=True)
    tables: dict[str, tuple[list[str], list[dict[str, str]]]] = {}
    row_counts: dict[str, tuple[int, int]] = {}
    for sheet_name in PMID_SHEETS:
        input_csv = input_dir / f"{sheet_name}_PMID_replaced.csv"
        with input_csv.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            require_columns(input_csv, ("RefID",), reader.fieldnames)
            fieldnames = reader.fieldnames or []
            source_rows = list(reader)
        retained_rows = []
        for row in source_rows:
            raw_ref_id = (row["RefID"] or "").strip()
            if not raw_ref_id:
                retained_rows.append(row)
                continue
            ref_id = int(raw_ref_id)
            if group_keys.get(ref_id, ref_id) == ref_id:
                retained_rows.append(row)
        output_csv = output_dir / f"{sheet_name}_groupkey_deduplicated.csv"
        with output_csv.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(retained_rows)
        tables[sheet_name] = (fieldnames, retained_rows)
        row_counts[sheet_name] = (len(source_rows), len(retained_rows))
    tables["Summary"] = (
        ["Sheet", "RowsBefore", "RowsAfter", "RowsRemoved"],
        [
            {
                "Sheet": sheet_name,
                "RowsBefore": str(before),
                "RowsAfter": str(after),
                "RowsRemoved": str(before - after),
            }
            for sheet_name, (before, after) in row_counts.items()
        ],
    )
    with found_pmid_report_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        require_columns(found_pmid_report_csv, ("Sheet", "FoundPMIDCount"), reader.fieldnames)
        tables["Found_PMID_report"] = (reader.fieldnames or [], list(reader))
    write_csv_tables_xlsx(tables, output_dir / "Original_sheets_groupkey_deduplicated.xlsx")
    return row_counts


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, default=DEFAULT_WORKBOOK)
    parser.add_argument("--sheet", default="Original")
    parser.add_argument("--accessions-csv", type=Path, default=DATA_DIR / "Accessions.csv")
    parser.add_argument("--ref-csv", type=Path, default=DEFAULT_REF_CSV)
    parser.add_argument(
        "--ref-with-refnames-output-csv", type=Path, default=DEFAULT_REF_WITH_REFNAMES_OUTPUT,
        help="CSV of Ref.csv with RefNames (default: outputs/Ref_same/01_ref_with_refname/Ref_with_RefName.csv).",
    )
    parser.add_argument(
        "--refname-summary-csv", type=Path, default=DEFAULT_REFNAME_SUMMARY,
        help="CSV for RefName counts (default: outputs/Ref_same/05_refname_counts/RefName_duplicate_counts.csv).",
    )
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument(
        "--groups-output-csv", type=Path, default=DEFAULT_GROUPS_OUTPUT,
        help="CSV for grouped candidate RefIDs (default: outputs/Ref_same/07_same_ref_candidate_groups/same_ref_candidate_groups.csv).",
    )
    parser.add_argument(
        "--group-refnames-output-csv", type=Path, default=DEFAULT_GROUP_REFNAMES_OUTPUT,
        help="CSV for candidate groups with RefNames (default: outputs/Ref_same/08_same_ref_candidate_group_refnames/same_ref_candidate_groups_refnames.csv).",
    )
    parser.add_argument(
        "--groupkey-ref-with-refname-output-csv", type=Path, default=DEFAULT_GROUPKEY_REF_WITH_REFNAME_OUTPUT,
        help="CSV with Ref_with_RefName rows deduplicated to group keys (default: outputs/Ref_same/09_ref_with_refname_groupkey_deduplication/Ref_with_RefName_groupkey_deduplicated.csv).",
    )
    parser.add_argument(
        "--deduplicated-refname-summary-csv", type=Path, default=DEFAULT_DEDUPLICATED_REFNAME_SUMMARY,
        help="CSV of RefName counts from step 09 (default: outputs/Ref_same/11_ref_with_refname_refname_counts/RefName_duplicate_counts.csv).",
    )
    parser.add_argument(
        "--duplicate-refname-rows-dir", type=Path, default=DEFAULT_DUPLICATE_REFNAME_ROWS_DIR,
        help="Directory for per-RefName duplicate-row CSVs (default: outputs/Ref_same/12_duplicate_refname_rows).",
    )
    parser.add_argument(
        "--merged-original-pmid-output-dir", type=Path, default=DEFAULT_MERGED_ORIGINAL_PMID_DIR,
        help="Directory for merged Original PMID CSVs (default: outputs/Ref_same/02_original_pmid_merge).",
    )
    parser.add_argument(
        "--found-pmid-output-csv", type=Path, default=DEFAULT_FOUND_PMID_OUTPUT,
        help="CSV of step-01 rows with filled PMIDs (default: outputs/Ref_same/03_found_pmid/Found_PMID.csv).",
    )
    parser.add_argument(
        "--original-sheets-replaced-pmid-output-dir", type=Path,
        default=DEFAULT_ORIGINAL_SHEETS_REPLACED_PMID_DIR,
        help="Directory for step-02 PMID-replaced copies of all Original sheets (default: outputs/Ref_same/04_original_sheets_pmid_replaced).",
    )
    parser.add_argument(
        "--original-sheets-groupkey-deduplication-output-dir", type=Path,
        default=DEFAULT_ORIGINAL_SHEETS_GROUPKEY_DEDUPLICATION_DIR,
        help="Directory for group-key-deduplicated step-04 sheet CSVs and Excel workbook (default: outputs/Ref_same/10_original_sheets_groupkey_deduplication).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    ref_with_refnames_rows = write_ref_with_refnames(
        args.ref_csv, args.workbook, args.sheet, args.ref_with_refnames_output_csv
    )
    merged_original_rows, numeric_pmid_rows, non_numeric_pmid_rows = write_merged_original_pmids(
        args.workbook, args.merged_original_pmid_output_dir
    )
    found_pmid_replacements = write_found_pmids(
        args.ref_with_refnames_output_csv,
        args.merged_original_pmid_output_dir / "Original_merged_PMID.csv",
        args.found_pmid_output_csv,
    )
    found_pmid_rows_by_sheet = write_found_pmid_report(
        args.workbook,
        args.ref_with_refnames_output_csv,
        args.found_pmid_output_csv,
        args.found_pmid_output_csv.parent / "Found_PMID_report.csv",
    )
    updated_pmid_rows_by_sheet = write_original_sheets_with_replaced_pmids(
        args.workbook,
        args.merged_original_pmid_output_dir / "Original_merged_PMID.csv",
        args.original_sheets_replaced_pmid_output_dir,
    )
    records = read_ref_records(args.found_pmid_output_csv)
    add_accession_data(records, args.accessions_csv)
    candidates = find_candidates(records)
    refname_groups = write_refname_summary(
        args.found_pmid_output_csv, args.refname_summary_csv
    )
    write_candidates(records, candidates, args.output_csv)
    candidate_groups = write_candidate_groups(candidates, args.groups_output_csv)
    candidate_group_refnames = write_candidate_group_refnames(
        args.groups_output_csv, args.found_pmid_output_csv, args.group_refnames_output_csv
    )
    original_rows_before, original_rows_after = write_groupkey_deduplicated_original_rows(
        args.group_refnames_output_csv, args.found_pmid_output_csv,
        args.groupkey_ref_with_refname_output_csv
    )
    original_sheets_row_counts = write_groupkey_deduplicated_original_sheets(
        args.group_refnames_output_csv,
        args.original_sheets_replaced_pmid_output_dir,
        args.original_sheets_groupkey_deduplication_output_dir,
        args.found_pmid_output_csv.parent / "Found_PMID_report.csv",
    )
    deduplicated_refname_groups = write_deduplicated_refname_summary(
        args.groupkey_ref_with_refname_output_csv, args.deduplicated_refname_summary_csv
    )
    duplicate_refname_files = write_duplicate_refname_rows(
        args.deduplicated_refname_summary_csv, args.groupkey_ref_with_refname_output_csv,
        args.duplicate_refname_rows_dir,
    )
    def print_step(step: int, name: str) -> None:
        print(f"\n=== Step {step:02d}: {name} ===")

    def display_path(path: Path) -> Path:
        """Prefer a repository-relative path in console output."""
        try:
            return path.resolve().relative_to(REPO_ROOT)
        except ValueError:
            return path.resolve()

    print_step(1, "Ref with RefName")
    print(f"output_csv={display_path(args.ref_with_refnames_output_csv)}")
    print(f"rows={ref_with_refnames_rows}")

    print_step(2, "Original PMID merge")
    print(f"output_dir={display_path(args.merged_original_pmid_output_dir)}")
    print(f"merged_original_rows={merged_original_rows}")
    print(f"numeric_pmid_rows={numeric_pmid_rows}")
    print(f"non_numeric_pmid_rows={non_numeric_pmid_rows}")

    print_step(3, "Found PMID")
    print(f"output_csv={display_path(args.found_pmid_output_csv)}")
    print(f"report_csv={display_path(args.found_pmid_output_csv.parent / 'Found_PMID_report.csv')}")
    print(f"found_pmid_replacements={found_pmid_replacements}")
    for sheet_name, found in found_pmid_rows_by_sheet.items():
        print(f"found_pmid_rows_{sheet_name}={found}")

    print_step(4, "Original sheets PMID replaced")
    print(f"output_dir={display_path(args.original_sheets_replaced_pmid_output_dir)}")
    for sheet_name, updates in updated_pmid_rows_by_sheet.items():
        print(f"pmid_rows_updated_{sheet_name}={updates}")

    print_step(5, "RefName counts")
    print(f"output_csv={display_path(args.refname_summary_csv)}")
    print(f"refname_groups={refname_groups}")

    print_step(6, "Same-reference candidates")
    print(f"output_csv={display_path(args.output_csv)}")
    print(f"refs_total={len(records)}")
    print(f"candidate_pairs={len(candidates)}")

    print_step(7, "Same-reference candidate groups")
    print(f"output_csv={display_path(args.groups_output_csv)}")
    print(f"candidate_groups={candidate_groups}")

    print_step(8, "Candidate groups with RefNames")
    print(f"output_csv={display_path(args.group_refnames_output_csv)}")
    print(f"candidate_group_refnames={candidate_group_refnames}")

    print_step(9, "Ref-with-RefName group-key deduplication")
    print(f"output_csv={display_path(args.groupkey_ref_with_refname_output_csv)}")
    print(f"original_rows_before={original_rows_before}")
    print(f"original_rows_after={original_rows_after}")

    print_step(10, "Original sheets group-key deduplication")
    print(f"output_dir={display_path(args.original_sheets_groupkey_deduplication_output_dir)}")
    for sheet_name, (before, after) in original_sheets_row_counts.items():
        print(f"groupkey_dedup_{sheet_name}_rows_before={before}")
        print(f"groupkey_dedup_{sheet_name}_rows_after={after}")

    print_step(11, "Deduplicated RefName counts")
    print(f"output_csv={display_path(args.deduplicated_refname_summary_csv)}")
    print(f"deduplicated_refname_groups={deduplicated_refname_groups}")

    print_step(12, "Duplicate RefName rows")
    print(f"output_dir={display_path(args.duplicate_refname_rows_dir)}")
    print(f"duplicate_refname_files={duplicate_refname_files}")


if __name__ == "__main__":
    main()
