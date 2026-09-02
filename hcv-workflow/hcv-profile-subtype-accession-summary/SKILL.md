---
name: hcv-profile-subtype-accession-summary
description: Aggregate the active One RAS NS3, NS5A, and NS5B workflows into per-genotype/subtype accession counts. Use when updating or reviewing cross-gene profile count summaries.
---

# HCV Profile Subtype Accession Summary

Run the summary script from the repository root after the gene profile workflows have produced their subtype RAS profile workbooks.

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py
```

The script reads the NS3, NS5A, and NS5B `*_Subtype_RAS_Profiles.xlsx` files. It uses the per-subtype RAS-coverage count shown in each profile label, matching the combined profile count. It writes both a cross-gene genotype/subtype comparison and a gene-by-gene subtype list as CSV and Excel files. Every subtype in an input profile is retained. The workbooks color counts of 10 or more blue; CSV does not support formatting.

## Default summary

Running with no arguments generates the active One RAS summary under `outputs/hcv-profile-subtype-accession-summary/one-ras/`. It retains every subtype, including those with fewer than 10 accessions. It uses the stage-23 `*_Subtype_RAS_Profiles.xlsx` workbooks from these active workflows.

The default summary also copies `HCV_Profile_Subtype_Accession_Counts_one-ras.xlsx` as `Table1_Gene_Subtype_Counts.xlsx` in the `one-ras/` folder.

| Summary | NS3 workflow | NS5A workflow | NS5B workflow | Output subfolder |
| --- | --- | --- | --- | --- |
| One RAS | `hcv-ns3-one-ras-comet-build-workflow` | `hcv-ns5a-one-ras-comet-build-workflow` | `hcv-ns5b-position-282-include-or-short-comet-build-workflow` | `one-ras/` |

Run it to generate the default One RAS summary:

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py \
  --output-dir outputs/hcv-profile-subtype-accession-summary
```

For one custom combination, pass all three profile inputs: `--ns3-ras-profile`, `--ns5a-ras-profile`, and `--ns5b-ras-profile`.
