---
name: hcv-profile-subtype-accession-summary
description: Aggregate NS3, NS5A, and one selected NS5B profile variant into per-genotype/subtype accession counts. Use when updating or reviewing cross-gene profile count summaries.
---

# HCV Profile Subtype Accession Summary

Run the summary script from the repository root after the gene profile workflows have produced their subtype RAS profile workbooks. NS3 and NS5A inputs are shared; generate a separate summary for each NS5B profile variant.

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py
```

The script reads the NS3, NS5A, and NS5B `*_Subtype_RAS_Profiles.xlsx` files. It uses the per-subtype RAS-coverage count shown in each profile label, matching the combined profile count. It writes both a cross-gene genotype/subtype comparison and a gene-by-gene subtype list as CSV and Excel files. Counts below 10 are omitted except that every genotype 7 and genotype 8 subtype present in an input profile is retained. The workbooks color counts of 10 or more blue; CSV does not support formatting.

## NS5B variants

Use the stage-23 `NS5B_Subtype_RAS_Profiles.xlsx` input and save each result in its matching subfolder under `outputs/hcv-profile-subtype-accession-summary/`:

| Variant | NS5B input | Summary output subfolder |
| --- | --- | --- |
| All RAS | `outputs/comet-NS5B-all-ras/23_build-subtype-ras-profile/NS5B_Subtype_RAS_Profiles.xlsx` | `all-ras/` |
| Position 282 | `outputs/comet-NS5B-position-282/23_build-subtype-ras-profile/NS5B_Subtype_RAS_Profiles.xlsx` | `position-282/` |
| Position 282 plus four RAS positions | `outputs/comet-NS5B-position-282-four-ras/23_build-subtype-ras-profile/NS5B_Subtype_RAS_Profiles.xlsx` | `position-282-four-ras/` |

Running the script with no arguments generates all three variants in these subfolders:

```bash
.venv/bin/python hcv-workflow/hcv-profile-subtype-accession-summary/scripts/build_profile_subtype_accession_summary.py \
  --output-dir outputs/hcv-profile-subtype-accession-summary
```

Pass `--ns5b-ras-profile` only to generate one custom NS5B summary directly in `--output-dir`.
