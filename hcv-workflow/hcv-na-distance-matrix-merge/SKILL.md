---
name: hcv-na-distance-matrix-merge
description: Merge the active NS3 One RAS, NS5A One RAS, and NS5B position-282 four-RAS NA distance matrices into side-by-side Excel workbooks. Use for preparing cross-gene NA distance-matrix tables.
---

# HCV NA Distance Matrix Merge

Run from the repository root after Step 29 has completed for the three source workflows:

```bash
.venv/bin/python hcv-workflow/hcv-na-distance-matrix-merge/scripts/merge_na_distance_matrices.py
```

The script reads these Step-29 workbooks:

- NS3 One RAS: `outputs/comet-NS3-one-ras/.../NS3_GT_NA_Distance_RAS.xlsx` and `NS3_Subtype_NA_Distance_RAS.xlsx`
- NS5A One RAS: `outputs/comet-NS5A-one-ras/.../NS5A_GT_NA_Distance_RAS.xlsx` and `NS5A_Subtype_NA_Distance_RAS.xlsx`
- NS5B position-282 four-RAS: `outputs/comet-NS5B-position-282-four-ras/.../NS5B_GT_NA_Distance_RAS.xlsx` and `NS5B_Subtype_NA_Distance_RAS.xlsx`

It writes to `outputs/hcv-na-distance-matrix-merge/`:

- `GT_NA_Distance_Merged.xlsx`: all three genotype `distance_matrix` sheets arranged left to right on one worksheet.
- `Subtype_NA_Distance_Merged.xlsx`: one worksheet per genotype; the corresponding NS3, NS5A, and NS5B subtype matrices are arranged left to right.

Only distance-matrix worksheets are included. Sequence-count, metadata, and exclusion worksheets are intentionally excluded.
