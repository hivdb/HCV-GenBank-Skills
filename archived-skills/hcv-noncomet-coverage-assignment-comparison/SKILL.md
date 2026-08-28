---
name: hcv-noncomet-coverage-assignment-comparison
description: Compare NS3, NS5A, and NS5B all-sequence non-COMET coverage calls with folder BLAST assignment CSVs. Use when accession counts or genotype/subtype disagreements need review.
---

# HCV non-COMET coverage versus folder-assignment comparison

Compare the three paired source files with the bundled script. It normalizes accession versions before matching and writes a summary plus one disagreement CSV per gene.

```bash
uv run python archived-skills/hcv-noncomet-coverage-assignment-comparison/scripts/compare_coverage_and_assignments.py
```

Default input pairs are:

- `HCVData/HCV-all-seq-subtype/NS3_AllSeq_NonComet_Coverage.csv` and `archived-skills/outputs/folder_assignments/NS3_assignments.csv`
- `HCVData/HCV-all-seq-subtype/NS5A_AllSeq_NonComet_Coverage.csv` and `archived-skills/outputs/folder_assignments/NS5A_assignments.csv`
- `HCVData/HCV-all-seq-subtype/NS5B_AllSeq_NonComet_Coverage.csv` and `archived-skills/outputs/folder_assignments/NS5B_assignments.csv`

Outputs default to `archived-skills/outputs/noncomet-coverage-assignment-comparison/`:

- `accession_count_summary.csv`: total and unique accession counts, overlap, and accession-set differences for each gene; the coverage-table side is labeled `FromDb` and the folder-FASTA side is labeled `FromFasta`.
- `{NS3,NS5A,NS5B}_genotype_subtype_differences.csv`: only accessions present in both sources whose genotype and/or subtype differs.
- `{NS3,NS5A,NS5B}_only_in_from_db.csv`: accessions absent from the corresponding folder assignment CSV.
- `{NS3,NS5A,NS5B}_only_in_from_fasta.csv`: accessions absent from the corresponding coverage CSV.

Blank genotype or subtype values are treated as values during comparison, so an assigned call versus a blank call is reported as a difference.
