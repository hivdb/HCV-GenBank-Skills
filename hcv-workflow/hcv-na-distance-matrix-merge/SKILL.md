---
name: hcv-na-distance-matrix-merge
description: Merge the active NS3 One RAS, NS5A One RAS, and NS5B position-282 four-RAS genotype and subtype NA distance matrices into one gene-grouped Excel workbook. Use for preparing cross-gene NA distance-matrix tables.
---

# HCV NA Distance Matrix Merge

Run from the repository root after the paired-distance-matrix step has completed for the three source workflows:

```bash
.venv/bin/python hcv-workflow/hcv-na-distance-matrix-merge/scripts/merge_na_distance_matrices.py
```

The script reads these paired-distance-matrix workbooks:

- NS3 One RAS (Step 30): `outputs/comet-NS3-one-ras/.../NS3_GT_NA_Distance_RAS.xlsx` and `NS3_Subtype_NA_Distance_RAS.xlsx`
- NS5A One RAS: `outputs/comet-NS5A-one-ras/.../NS5A_GT_NA_Distance_RAS.xlsx` and `NS5A_Subtype_NA_Distance_RAS.xlsx`
- NS5B position-282 four-RAS: `outputs/comet-NS5B-position-282-four-ras/.../NS5B_GT_NA_Distance_RAS.xlsx` and `NS5B_Subtype_NA_Distance_RAS.xlsx`

It writes one final workbook to `outputs/hcv-na-distance-matrix-merge/`:

- `NA_Distance_Matrics.xlsx`: one worksheet with NS3, NS5A, and NS5B arranged left to right. Within each gene block, the genotype `distance_matrix` is first and that gene's genotype-specific subtype matrices are stacked underneath it. The displayed NS5B heading is `NS5B Five or more`.

Use a distinct light background for every copied cell in each gene block: light blue (`DDEBF7`) for NS3, light green (`E2F0D9`) for NS5A, and light orange (`FCE4D6`) for NS5B. Keep spacer rows and columns unfilled.

Only distance-matrix worksheets are included. Sequence-count, metadata, and exclusion worksheets are intentionally excluded.
