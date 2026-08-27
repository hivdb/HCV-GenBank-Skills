# 2026-08-27: folder genotype/subtype assignment outputs

## Files in this directory

- `NS3_assignments.csv`
- `NS5A_assignments.csv`
- `NS5B_assignments.csv`
- `assignment_summary.json`

## How they are generated

Run the folder assignment script from the repository root:

```bash
.venv/bin/python hcv-workflow/hcv-folder-genotype-subtype-assignment/scripts/assign_folder_genotype_subtype.py \
  --fasta-dir /path/to/input_fasta_folder \
  --output-dir outputs/folder_assignments
```

The input directory is searched recursively for `.fa`, `.fasta`, and `.fna` files. Each record is evaluated independently for NS3, NS5A, and NS5B.

For each detected gene, the script uses `blastn` in two stages:

1. Assign genotype against `HCVData/Reference_seqs/HCV_GT_RefSeqs.fasta`.
2. Assign the closest subtype only among the subtype references for that genotype in `HCVData/Reference_seqs/HCV_Subtype_Refs_By_Genome_NA.json`.

The default minimum comparable length is 200 nt. A blank subtype means no subtype assignment passed the workflow criteria; it should not be inferred from genotype alone.

## CSV columns

Each gene-specific CSV contains:

| Column | Meaning |
| --- | --- |
| `accession` | Accession parsed from the FASTA record. |
| `header` | Original FASTA header. |
| `source_fasta` | Input FASTA file containing the record. |
| `gene` | Detected gene: NS3, NS5A, or NS5B. |
| `genotype` | First-stage BLAST genotype call. |
| `genotype_aligned_nt` / `genotype_pident` | Genotype alignment length and percent identity. |
| `subtype` | Second-stage, genotype-matched BLAST subtype call. |
| `subtype_reference_accession` | Best matching subtype reference accession. |
| `subtype_aligned_nt` / `subtype_pident` | Subtype alignment length and percent identity. |
