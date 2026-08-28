---
name: hcv-original-refid-pmid-audit
description: Combine RefID/PMID pairs from the Original_NS3, Original_NS5A, and Original_NS5B workbook sheets, then identify records whose matching Ref.csv entry has a blank PMID (MedlineID). Use when auditing missing reference PMIDs in HCV source data.
---

# HCV Original RefID PMID Audit

Run the audit from the repository root:

```bash
uv run python hcv-workflow/hcv-original-refid-pmid-audit/scripts/audit_original_refid_pmids.py
```

The script reads these sources by default:

- `HCVData/HCV_BlastHists_202604_data.xlsx`, sheets `Original_NS3`,
  `Original_NS5A`, and `Original_NS5B`
- `HCVData/HCV-all-seq-subtype/Ref.csv`

It writes CSV outputs to `outputs/hcv-original-refid-pmid-audit/`:

- `combined_refid_pmid_pairs.csv`: every source row with a nonblank `RefID` and
  `PMID`, retaining its source sheet and Excel row number.
- `found_refids_with_blank_ref_pmid.csv`: combined rows whose `RefID` exists in
  `Ref.csv` and whose `MedlineID` (the Ref.csv PMID field) is blank.
- `NS3/`, `NS5A/`, and `NS5B/`: each contains its own
  `found_refids_with_blank_ref_pmid.csv` and one-row `report.csv`.
- `PMID_found_summary_by_sheet.csv`: the three gene reports together.

Each report has `Total rows` in column C, then `Final with PMID`, `Found PMID`,
and `PMID from genbank`. `Found PMID` counts matching Ref.csv rows with an
empty PMID whose Original-sheet row has a populated PMID. `PMID from genbank`
is `Final with PMID − Found PMID`. The terminal prints the same per-gene
summary. To use alternate inputs or an output directory, pass `--workbook`,
`--ref-csv`, or `--output-dir`.
