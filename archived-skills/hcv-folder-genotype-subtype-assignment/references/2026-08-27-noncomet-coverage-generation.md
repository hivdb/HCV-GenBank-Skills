# 2026-08-27: all-sequence non-COMET coverage tables

## Purpose

These files provide non-COMET genotype/subtype calls and gene-position coverage for the complete HCV sequence collection:

- `HCVData/nonComet-Full-genome/NS3_AllSeq_NonComet_Coverage.csv`
- `HCVData/nonComet-Full-genome/NS5A_AllSeq_NonComet_Coverage.csv`
- `HCVData/nonComet-Full-genome/NS5B_AllSeq_NonComet_Coverage.csv`

They are subsequently used by the NS3, NS5A, and NS5B COMET workflows as priority non-COMET subtype annotations. In particular, the workflows use genotype 7/8 and subtype 1d records from these tables to supplement or override COMET calls.

## Generation chain

```text
HCVData/HCV-all-seq-subtype/all.fasta
  -> hcv-folder-genotype-subtype-assignment (BLAST genotype then genotype-matched subtype)
  -> accession-level NS3/NS5A/NS5B assignments in a temporary directory
  -> coverage mapping against the gene reference targets
  -> *_AllSeq_NonComet_Coverage.csv
```

Run the coverage-audit script from the repository root:

```bash
.venv/bin/python Preprocess/Ref/hcv-all-seq-blast-subtyping-full-genome/scripts/audit_all_fasta_coverage.py \
  --input-fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --output-dir HCVData/nonComet-Full-genome \
  --threads 4
```

The audit script invokes this folder assignment script internally:

```text
archived-skills/hcv-folder-genotype-subtype-assignment/scripts/assign_folder_genotype_subtype.py
```

It runs `blastn` in two stages: genotype against the gene/genotype references, then subtype against only reference sequences from the assigned genotype.

## Output columns

Each coverage CSV contains:

| Column | Meaning |
| --- | --- |
| `Accession` | Sequence accession from `all.fasta`. |
| `ClosestGenotype` | BLAST genotype assignment. |
| `ClosestSubtype` | BLAST subtype assignment within that genotype. |
| `ReferenceOverlapAA` | Overlapping amino-acid interval in reference coordinates. |
| `FullyCover` | `Yes` only when the entire target interval is covered. |

Target intervals are NS3 36–175, NS5A 26–93, and NS5B 150–321.

`all.fasta` currently contains the nucleotide sequences corresponding to the `NASeq` field in `HCVData/HCV-all-seq-subtype/Accessions.csv`. It is a local ignored file; no generator for that FASTA export is currently tracked in this repository.
