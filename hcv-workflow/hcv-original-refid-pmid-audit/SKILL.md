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

The terminal output reports the number of found rows. To use alternate inputs or
an output directory, pass `--workbook`, `--ref-csv`, or `--output-dir`.
