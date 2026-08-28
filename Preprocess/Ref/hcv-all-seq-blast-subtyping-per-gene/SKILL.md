---
name: hcv-all-seq-blast-subtyping-per-gene
description: Assign non-COMET HCV genotype and closest within-genotype subtype against separate NS3, NS5A, and NS5B nucleotide references, then report target-range coverage. Use when gene-specific reference FASTAs are required.
---

# HCV all-sequence BLAST subtyping with per-gene references

Run the bundled script from the repository root. It performs genotype and
genotype-matched subtype assignment directly from the input FASTA, using
separate reference FASTAs for each gene, then writes one CSV per gene in the
requested output directory.

- Genotype references: `HCVData/Genotype-Ref/HCV_GT_Refs_NS3_NA.fasta`, `HCVData/Genotype-Ref/HCV_GT_Refs_NS5A_NA.fasta`, and `HCVData/Genotype-Ref/HCV_GT_Refs_NS5B_NA.fasta`
- Subtype references: `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS3_NA.fasta`, `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5A_NTD_NA.fasta`, and `HCVData/Subtype-Ref/HCV_Subtype_Refs_NS5B_NA.fasta`

```bash
.venv/bin/python Preprocess/Ref/hcv-all-seq-blast-subtyping-per-gene/scripts/audit_all_fasta_coverage.py \
  --input-fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --output-dir HCVData/nonComet-PerGene \
  --threads 4
```

The output files are `NS3_AllSeq_NonComet_Coverage.csv`, `NS5A_AllSeq_NonComet_Coverage.csv`, and `NS5B_AllSeq_NonComet_Coverage.csv`.

Each table has seven columns: `Accession`, `ClosestGenotype`, `ClosestGenotypePident`, `ClosestSubtype`, `ClosestSubtypePident`, `ReferenceOverlapAA`, and `FullyCover`. The percent-identity columns are the BLAST percent identities for the best gene-specific genotype and genotype-matched subtype hits. Blank assignment fields mean the sequence did not meet the assignment threshold.

`ReferenceOverlapAA` is blank when the best genotype alignment does not overlap the requested target range. When it overlaps, it reports the overlapping reference amino-acid interval, including partial overlap.

`FullyCover` is `Yes` only when the best genotype alignment spans the whole requested target range; it is blank for partial or absent overlap.

The script displays a live stage-level progress bar. Its BLAST searches use four workers by default; change `--threads` only when the available CPU capacity requires it.

Keep the default targets unless explicitly requested otherwise:

- NS3: amino-acid positions 36-175
- NS5A: amino-acid positions 26-93
- NS5B: amino-acid positions 150-321
