---
name: hcv-folder-genotype-subtype-assignment
description: Assign HCV genotype then genotype-matched subtype for every sequence in a folder of FASTA files. Detect NS3, NS5A, and NS5B sequence regions and write one accession-level assignment CSV per detected gene. Use when a user provides a FASTA folder and needs sequence-based HCV genotype/subtype assignments without study-workbook inputs.
---

# HCV Folder Genotype/Subtype Assignment

Run `assign_folder_genotype_subtype/assign_folder_genotype_subtype.py` when the user supplies the input folder.

```bash
uv run python hcv-workflow/hcv-folder-genotype-subtype-assignment/assign_folder_genotype_subtype/assign_folder_genotype_subtype.py \
  --fasta-dir /path/to/fasta_folder
```

The script searches `.fa`, `.fasta`, and `.fna` files recursively. It first assigns genotype by nucleotide alignment to the bundled gene/genotype references, then assigns the nearest subtype only among references for that genotype. It emits `NS3_assignments.csv`, `NS5A_assignments.csv`, `NS5B_assignments.csv`, and `assignment_summary.json`. Preserve every result row; do not infer an assignment for an absent gene or a sequence failing the minimum aligned-nucleotide threshold.

Defaults use `HCV_GT_RefSeqs.fasta`, `HCV_Subtype_Refs_By_Genome_NA.json`, and 200 aligned nt. Override them only when the user requests different references or threshold.
The default output directory is `outputs/folder_assignments`; pass `--output-dir` to override it.

## Related coverage tables

For the documented generation chain of the all-sequence non-COMET coverage tables, read [2026-08-27-noncomet-coverage-generation.md](references/2026-08-27-noncomet-coverage-generation.md).
