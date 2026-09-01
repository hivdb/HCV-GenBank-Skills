---
name: hcv-profile-subtype-accession-summary
description: Aggregate three predefined NS3, NS5A, and NS5B workflow combinations into per-genotype/subtype accession counts. Use when updating or reviewing cross-gene profile count summaries.
---

# HCV Profile Subtype Accession Summary

Run the summary script from the repository root after the gene profile workflows have produced their subtype RAS profile workbooks.

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py
```

The script reads the NS3, NS5A, and NS5B `*_Subtype_RAS_Profiles.xlsx` files. It uses the per-subtype RAS-coverage count shown in each profile label, matching the combined profile count. It writes both a cross-gene genotype/subtype comparison and a gene-by-gene subtype list as CSV and Excel files. Counts below 10 are omitted except that every genotype 7 and genotype 8 subtype present in an input profile is retained. The workbooks color counts of 10 or more blue; CSV does not support formatting.

## Summary versions

Running with no arguments generates these three version subfolders under `outputs/hcv-profile-subtype-accession-summary/`. Each also contains an `all-subtypes/` subfolder that retains every subtype, including counts below 10. Each uses the stage-23 `*_Subtype_RAS_Profiles.xlsx` workbooks from the listed workflows. The two generated Excel workbook filenames include the version condition.

| Version | NS3 workflow | NS5A workflow | NS5B workflow | Summary output subfolder |
| --- | --- | --- | --- | --- |
| All RAS | `hcv-ns3-all-ras-comet-build-workflow` | `hcv-ns5a-all-ras-comet-build-workflow` | `hcv-ns5b-all-ras-comet-build-workflow` | `all-ras/` |
| One RAS | `hcv-ns3-one-ras-comet-build-workflow` | `hcv-ns5a-one-ras-comet-build-workflow` | `hcv-ns5b-position-282-comet-build-workflow` | `one-ras/` |
| Position 282 plus five RAS positions | `hcv-ns3-comet-build-workflow` | `hcv-ns5a-comet-build-workflow` | `hcv-ns5b-position-282-five-ras-comet-build-workflow` | `position-282-five-ras/` |

Run it to generate all three versions:

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py \
  --output-dir outputs/hcv-profile-subtype-accession-summary
```

For one custom combination, pass all three profile inputs: `--ns3-ras-profile`, `--ns5a-ras-profile`, and `--ns5b-ras-profile`.
