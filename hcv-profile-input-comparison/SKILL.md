---
name: hcv-profile-input-comparison
description: Generate Comet profile-input accession CSVs and compare genotype/subtype assignments with local-alignment results for HCV NS3, NS5A, and NS5B. Use when profile-input CSVs or cross-workflow assignment comparison reports are needed.
---

# HCV Profile Input Comparison

Generate the three Comet profile-input CSVs from their retained amino-acid source workbooks, then compare them with the matching local-alignment CSVs.

## Run

Run after the normal and Comet workflows finish:

```bash
.venv/bin/python hcv-profile-input-comparison/scripts/build_and_compare_profile_inputs.py --repo-root .
```

The script reads `outputs/comet/<GENE>_Profile_Input_Source.xlsx`, writes the Comet `Profile_Input_Accessions.csv`, and writes comparison and accession-level differences CSVs under `outputs/comet/`.

Rows without usable amino-acid data, or with genotype/subtype beginning with `unassign`, are excluded. The local input is read from `outputs/local_alignment/<GENE>_Profile_Input_Accessions.csv`.
