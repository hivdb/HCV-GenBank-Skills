---
name: hcv-for-paper-collection-workflow
description: Collect the requested HCV Comet publication workbooks into the repository for_paper folder with standardized table and figure filenames.
---

# HCV For-Paper Collection Workflow

Run this workflow after the active Comet workflows have produced their profile outputs. It copies the configured source workbooks into the repository-root `for_paper/` directory and applies the publication filenames.

```bash
.venv/bin/python hcv-workflow/hcv-for-paper-collection-workflow/scripts/collect_for_paper.py
```

The workflow is reproducible: source files remain in their original numbered workflow directories, while the collected copies are written to `for_paper/`.

## Collected outputs

- `Table S1 - References.xlsx`
- `Table S2 - NS3.xlsx`
- `Table S3 - NS5A.xlsx`
- `Table S4 - NS5B.xlsx`
- `Figure S1 - NS3 Profile - all subtypes.xlsx`
- `Figure S3 - NS5A Profile - all subtypes.xlsx`
- `Figure S5 - NS5B Profile - all subtypes.xlsx`

Use `--repo-root` to run against a different checkout and `--output-dir` to override the destination folder.
