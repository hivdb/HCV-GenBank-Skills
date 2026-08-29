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

The output files are `NS3_AllSeq_NonComet_Coverage.csv`, `NS5A_AllSeq_NonComet_Coverage.csv`, and `NS5B_AllSeq_NonComet_Coverage.csv`. Each gene also has a `*_Genotype_Distances.csv` file, one `*_Subtype_Distances_Genotype_<genotype>.csv` file per genotype, and a `*_Subtyping_Distances.xlsx` workbook.

Each table has seven columns: `Accession`, `ClosestGenotype`, `ClosestGenotypePident`, `ClosestSubtype`, `ClosestSubtypePident`, `ReferenceOverlapAA`, and `FullyCover`. The percent-identity columns are the BLAST percent identities for the best gene-specific genotype and genotype-matched subtype hits. Blank assignment fields mean the sequence did not meet the assignment threshold.

Each genotype-distance CSV has one row per accession, a distance column for every reference genotype, and the first and second choices with their distances. Each genotype-specific subtype-distance CSV does the same only for subtype references and accessions of that genotype. Distance is BLAST nucleotide percent-identity distance (`100 - pident`); blank values did not meet the minimum aligned-length threshold.

Each gene's workbook combines its distance CSVs: `Genotype` contains the genotype-distance table, and every `Subtype_<genotype>` sheet contains that genotype's subtype-distance table.

## Plot genotype-distance groups

Use the histogram script to inspect the distributions of `FirstChoiceDistance`
and `SecondChoiceDistance`. It uses one-distance-unit bins, which makes the
value groups directly comparable, and writes PNG files to a `figures`
subfolder beside the input CSV by default.
Use `--skip-empty` when processing every distance CSV so that files without a
second-choice value are reported and skipped instead of stopping the batch.

```bash
uv run python Preprocess/Ref/hcv-all-seq-blast-subtyping-per-gene/scripts/plot_distance_histogram.py \
  --input-csv HCVData/nonComet-PerGene/NS3_Genotype_Distances.csv
uv run python Preprocess/Ref/hcv-all-seq-blast-subtyping-per-gene/scripts/plot_distance_histogram.py \
  --input-csv HCVData/nonComet-PerGene/NS3_Genotype_Distances.csv \
  --column SecondChoiceDistance
```

For the NS3 input above, the image is
`HCVData/nonComet-PerGene/figures/NS3_FirstChoiceDistance_Histogram.png` and
`HCVData/nonComet-PerGene/figures/NS3_SecondChoiceDistance_Histogram.png`.

`ReferenceOverlapAA` is blank when the best genotype alignment does not overlap the requested target range. When it overlaps, it reports the overlapping reference amino-acid interval, including partial overlap.

`FullyCover` is `Yes` only when the best genotype alignment spans the whole requested target range; it is blank for partial or absent overlap.

The script displays a live stage-level progress bar. Its BLAST searches use four workers by default; change `--threads` only when the available CPU capacity requires it.

The default `--min-aligned-nt` is 100. Use the option to set a different
minimum aligned nucleotide length for genotype and subtype calls.

Keep the default targets unless explicitly requested otherwise:

- NS3: amino-acid positions 36-175
- NS5A: amino-acid positions 26-93
- NS5B: amino-acid positions 150-321
