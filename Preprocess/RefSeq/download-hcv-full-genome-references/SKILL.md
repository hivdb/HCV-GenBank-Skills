---
name: download-hcv-full-genome-references
description: Download full HCV genotype and subtype nucleotide reference genomes and create combined FASTA files under HCVData/Genotype-Ref and HCVData/Subtype-Ref. Use when these complete reference sequences need to be refreshed or inspected; do not use for gene-only reference extraction.
---

# Download HCV full-genome references

Run the bundled downloader from the repository root:

```bash
python3 Preprocess/RefSeq/download-hcv-full-genome-references/scripts/download_hcv_full_genome_references.py
```

It derives the eight genotype reference accessions from
`HCVData/HCV_GT_RefSeqs.fasta` and subtype reference accessions/subtype labels
from `HCVData/HCV_Subtype_Refs_By_Genome_NA.json`. It downloads their complete
nucleotide records from NCBI and writes:

- `HCVData/Genotype-Ref/HCV_GT_FullGenome_RefSeqs.fasta`
- `HCVData/Subtype-Ref/HCV_Subtype_FullGenome_Refs.fasta`

Individual records are cached under each matching output directory's
`genbank_records/` folder; reruns reuse valid cache files. Use `--force` to
redownload every reference, or `--genotype-output-dir` and
`--subtype-output-dir` to choose alternate locations. The downloader verifies
that each FASTA record has a sequence and reports the reference counts and
output paths.

## Build subtype per-gene nucleotide references

After the full subtype-genome FASTA exists, use the per-gene subtype AA FASTAs
to identify the matching gene sequence in each subtype genome and write
per-gene NA FASTAs:

```bash
python3 Preprocess/RefSeq/download-hcv-full-genome-references/scripts/build_hcv_subtype_per_gene_na_references.py
```

This uses BLASTX to align each complete subtype genome to its own NS3,
NS5A_NTD, or NS5B AA reference. It writes:

- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS3_NA.fasta`
- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5A_NTD_NA.fasta`
- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5B_NA.fasta`

Each output header retains the corresponding subtype AA reference metadata and
adds the full-genome nucleotide interval and BLASTX alignment coverage.
`makeblastdb` and `blastx` must be available on `PATH`.

## Align subtype gene NA references to subtype 1a

After building the per-gene NA files, create one 1a-anchored multiple alignment
per gene:

```bash
python3 Preprocess/RefSeq/download-hcv-full-genome-references/scripts/align_hcv_subtype_gene_na_to_1a.py
```

It puts the first subtype-1a reference in each input file first, aligns every
subtype gene NA sequence with MAFFT, and writes:

- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS3_NA_Aligned_to_1a.fasta`
- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5A_NTD_NA_Aligned_to_1a.fasta`
- `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5B_NA_Aligned_to_1a.fasta`

`mafft` must be available on `PATH`.
Use `--gene NS3`, `--gene NS5A_NTD`, or `--gene NS5B` to run selected genes.

## Validate per-gene NA sequences against full subtype genomes

Check that every unaligned per-gene NA reference is a continuous exact subset
of the matching full subtype-genome NA sequence:

```bash
python3 Preprocess/RefSeq/download-hcv-full-genome-references/scripts/validate_hcv_subtype_gene_na_subsets.py
```

It writes `HCVData/Subtype-Ref/HCV_Subtype_Gene_NA_Subset_Check.csv` with one
row per gene/accession, including the one-based full-genome interval and a
`Status` value. The command fails if any sequence is absent or not continuous.

## Split and validate genotype per-gene NA references

Split the combined genotype per-gene NA FASTA by gene and confirm each sequence
is a continuous exact subset of its matching full genotype genome:

```bash
python3 Preprocess/RefSeq/download-hcv-full-genome-references/scripts/split_and_validate_hcv_genotype_gene_na_references.py
```

It writes `HCV_GT_Refs_NS3_NA.fasta`, `HCV_GT_Refs_NS5A_NA.fasta`, and
`HCV_GT_Refs_NS5B_NA.fasta` plus
`HCV_GT_Gene_NA_Subset_Check.csv` under `HCVData/Genotype-Ref/`.
