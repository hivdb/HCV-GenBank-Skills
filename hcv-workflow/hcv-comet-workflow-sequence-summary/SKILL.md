---
name: hcv-comet-workflow-sequence-summary
description: Generate an Excel summary of NS3, NS5A, and NS5B COMET workflow sequence-inclusion rules and current accession counts. Use when comparing COMET build workflow eligibility across HCV genes.
---

# HCV COMET Workflow Sequence Summary

Generate the summary with:

```bash
uv run python hcv-workflow/hcv-comet-workflow-sequence-summary/scripts/build_comet_workflow_sequence_summary.py
```

The script reads each workflow's `15_report-profile-input-counts/profile_input_counts.json` count file and writes `outputs/HCV_Comet_Workflow_Sequence_Inclusion_Summary.xlsx`. Some NS5B count files include a descriptive line before the JSON payload; the script handles that format.

The workbook contains these columns: workflow name, gene, sequence-inclusion method for the combined profile, and included accession count. Keep all nine COMET workflow variants in the report: three each for NS3, NS5A, and NS5B.

If a required count file is missing, stop and report the missing workflow output; do not substitute a stale or inferred count.
