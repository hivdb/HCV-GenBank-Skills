---
name: hcv-profile-subtype-accession-summary
description: Aggregate NS3, NS5A, and NS5B profile-accession CSVs into per-genotype/subtype accession counts. Use when updating or reviewing the cross-gene profile count summary.
---

# HCV Profile Subtype Accession Summary

Run the summary script from the repository root after the gene profile workflows have produced their subtype RAS profile workbooks.

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py
```

The script reads the NS3, NS5A, and NS5B `*_Subtype_RAS_Profiles.xlsx` files. It uses the per-subtype RAS-coverage count shown in each profile label, matching the combined profile count. It writes both a cross-gene genotype/subtype comparison and a gene-by-gene subtype list as CSV and Excel files under `outputs/hcv-profile-subtype-accession-summary/`; both omit counts below 10. The workbooks color displayed counts blue; CSV does not support formatting.
