---
name: search-same-ref
description: Find likely duplicate HCV RefID records from Ref.csv and its accession mapping, with auditable matching evidence. Use when reviewing or consolidating duplicate reference records; do not use for accession-sequence duplicate detection.
---

# Search same HCV references

Run the bundled script from the repository root. It reads the `Original` sheet
from `HCV_BlastHists_202604_data.xlsx` plus the corresponding `Accessions.csv`.
Results are written under `outputs/Ref_same/` in this order, with each
step's order and name as its subfolder:

1. `01_ref_with_refname/Ref_with_RefName.csv` — `HCVData/HCV-all-seq-subtype/Ref.csv`
   with `MedlineID` renamed to `PMID` and `RefName` appended by matching its
   RefID to the `Original` worksheet. `Year` is extracted from `RefName`.
2. `02_original_pmid_merge/` — merges non-empty PMIDs by RefID from `Original`,
   `Original_NS5A`, `Original_NS3`, and `Original_NS5B` into an Original-sheet
   CSV, plus copies containing all non-empty PMIDs, numeric-only PMIDs, and
   non-numeric PMIDs.
3. `03_found_pmid/Found_PMID.csv` — copies step 01 and fills blank `PMID`
   values from step 02 by RefID. `Found_PMID_report.csv` stores the found-PMID
   count for each Original worksheet: its RefIDs with blank step-01 PMID values
   that have a non-empty PMID in `Found_PMID.csv`. The script prints those
   counts and `found_pmid_replacements`.
4. `04_original_sheets_pmid_replaced/` — loads `Original`, `Original_NS5A`,
   `Original_NS3`, and `Original_NS5B` from the workbook and writes one CSV per
   sheet after replacing its `PMID` by `RefID` from step 02's
   `Original_merged_PMID.csv`. The script prints the changed-row count for each
   sheet.
5. `05_refname_counts/RefName_duplicate_counts.csv` — one row per non-empty
   `RefName`, sorted by descending record count, with `RefName`, `RefIDCount`,
   and `RefIDs` columns.
6. `06_same_ref_candidates/same_ref_candidates.csv` — directly matched RefID
   candidate pairs.
7. `07_same_ref_candidate_groups/same_ref_candidate_groups.csv` — each
   candidate `GroupKeyRefID`, `RefIDCount`, and complete semicolon-separated
   `RefIDs` list.
8. `08_same_ref_candidate_group_refnames/same_ref_candidate_groups_refnames.csv`
   — step 07 with a deduplicated `RefNameList` and `RefNameListCount`, looked
   up from step 03 for every grouped RefID.
9. `09_ref_with_refname_groupkey_deduplication/Ref_with_RefName_groupkey_deduplicated.csv`
   — copies step 03 and removes every non-key group-member row from step 08.
   All columns and values in retained rows are unchanged. The
   script prints `original_rows_before` and `original_rows_after`.
10. `10_original_sheets_groupkey_deduplication/` — applies the step-09 group-key
   deduplication to every step-04 sheet CSV. It writes one unchanged-column CSV
   per sheet plus `Original_sheets_groupkey_deduplicated.xlsx`, containing those
   four CSV tables as worksheets and a `Summary` worksheet with the before,
   after, and removed row counts. It also copies step 03's `Found_PMID_report.csv`
   to a separate `Found_PMID_report` worksheet. The script prints each sheet's row
   count before and after deduplication.
11. `11_ref_with_refname_refname_counts/RefName_duplicate_counts.csv` — applies step
   05's `RefName`, `RefIDCount`, and `RefIDs` summary to the retained step-09
   rows.
12. `12_duplicate_refname_rows/` — one CSV for every step-11 RefName group with
   `RefIDCount >= 2`. Each CSV contains the unchanged matching step-09 rows and
   is stored in `RefIDCount_<count>/` as `<RefName>__RefIDCount_<count>.csv`,
   with unsafe filename characters replaced by underscores.

```bash
python3 Preprocess/Ref/search_same_ref/search_same_ref.py
```

The default match rules reproduce the legacy database script: the first-author
surname and representative worksheet year must agree (year difference at most
one), the records must share a six-character accession prefix, and they must
also match on PMID, highly similar title, or sufficiently similar journal.

It also adds RefName-based candidate paths: an exact non-empty RefName with a
title similarity of at least 0.92, or an exact non-empty RefName with the same
normalized journal. These paths do not require the legacy author, year, or
accession gates.

Any two RefIDs with the same exact non-empty stored PMID value are direct
candidates, regardless of their other metadata. Any pair with different
non-empty PMID values is excluded. The different-PMID rule is also enforced
while building connected groups, so blank-PMID records cannot transitively join
records with different known PMIDs.

Each output row is a direct candidate pair, rather than an inferred statement
that every pair in a transitive group matches. `GroupKeyRefID` identifies the
connected candidate group. `Reason` records all evidence for that pair;
`MatchingPMID`, `TitleSimilarity`, and `JournalSimilarity` expose the specific
metadata result. Review candidates before merging records.

Use `--workbook`, `--sheet`, `--accessions-csv`, or `--output-csv` to analyze
another export. The worksheet must provide `RefID`, `RefName`, `Year`, `PMID`,
`Title`, `Author`, and `Journal`; the accession CSV must provide `RefID` and
`Accession`.
