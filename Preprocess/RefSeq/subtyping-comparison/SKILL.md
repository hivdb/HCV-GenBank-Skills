---
name: subtyping-comparison
description: Filter assigned HCV per-gene non-COMET BLAST calls and summarize subtype agreement with COMET. Use when comparing NS3, NS5A, and NS5B per-gene calls against all_comet_subtype.csv.
---

# HCV per-gene non-COMET versus COMET comparison

Run from the repository root:

```bash
.venv/bin/python Preprocess/RefSeq/subtyping-comparison/scripts/compare_pergene_noncomet_comet.py
```

The workflow reads the three coverage CSVs from `HCVData/nonComet-PerGene/`
and COMET calls from `HCVData/HCV-all-seq-subtype/all_comet_subtype.csv`.

For each gene, it removes rows with a blank `ClosestGenotype`, joins COMET by
accession, and compares `ClosestSubtype` with the COMET subtype. A blank subtype
is counted as different from a non-blank subtype.

Outputs are written to `HCVData/subtyping-comparison/`:

- `Subtyping_Comparison_summary.csv`, containing `Gene`, `SameCount`, and
  `DifferentCount`;
- `Subtyping_Comparison_differences.csv`, containing `Gene`, `Accession`,
  `CometSubtype`, and `BlastPerGeneSubtype` for every difference.
- `Subtyping_Comparison_difference_subtype_counts.csv`, grouping differences
  by gene and the COMET/BLAST subtype pair with `DifferentCount`.

Every output includes a final `CompareCondition` column set to
`Comet_vs_PerGene`.

Use `--noncomet-dir`, `--comet-csv`, or `--output-dir` to override the default
paths.
