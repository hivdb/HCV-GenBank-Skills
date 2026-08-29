---
name: subtyping-comparison
description: Compare HCV COMET, per-gene BLAST, and full-genome BLAST subtype calls. Use when evaluating pairwise or three-way subtype agreement for NS3, NS5A, and NS5B.
---

# HCV COMET, per-gene, and full-genome subtyping comparison

Run from the repository root:

```bash
.venv/bin/python Preprocess/RefSeq/subtyping-comparison/scripts/compare_pergene_noncomet_comet.py
```

The workflow reads the three coverage CSVs from both
`HCVData/nonComet-PerGene/` and `HCVData/nonComet-Full-genome/`, plus COMET
calls from `HCVData/Comet-Full-genome/all_comet_subtype.csv`.

For each BLAST method and gene, it excludes rows with a blank `ClosestGenotype`.
The assigned accessions in each PerGene CSV define the comparison universe for
that gene. It then joins calls by accession and creates four comparison
conditions:

- `Comet_vs_PerGene`
- `Comet_vs_FullGenome`
- `PerGene_vs_FullGenome`
- `Comet_vs_PerGene_vs_FullGenome`

A blank subtype is counted as different from a non-blank subtype.

Outputs are written to `HCVData/subtyping-comparison/`:

- `Subtyping_Comparison_summary.csv`, containing `Gene`, `SameCount`,
  `DifferentCount`, and `CompareCondition`;
- `Subtyping_Comparison_differences.csv`, containing `Gene`, `Accession`,
  `CometSubtype`, `BlastPerGeneSubtype`, `BlastFullGenomeSubtype`, and
  `CompareCondition` for every difference.
- `Subtyping_Comparison_difference_subtype_counts.csv`, grouping differences
  by gene, method subtype calls, and comparison condition with `DifferentCount`.
- `Subtyping_Comparison_difference_subtype_counts.xlsx`, containing the same
  subtype-difference counts with one comparison condition per worksheet.

Every output includes a final `CompareCondition` column set to
`Comet_vs_PerGene`.

Use `--noncomet-dir`, `--comet-csv`, or `--output-dir` to override the default
paths.
