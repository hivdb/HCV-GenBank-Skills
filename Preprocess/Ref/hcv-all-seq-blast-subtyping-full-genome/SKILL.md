---
name: hcv-all-seq-blast-subtyping-full-genome
description: Assign non-COMET HCV genotype and closest within-genotype subtype for every record in an all.fasta-style nucleotide FASTA, then report overlap with NS3 positions 36-175, NS5A positions 26-93, and NS5B positions 150-321. Use when an accession-level coverage table is needed from an HCV sequence collection without relying on COMET calls.
---

# HCV all-sequence BLAST subtyping and full-genome coverage

Run the bundled script from the repository root. It performs genotype and
genotype-matched subtype assignment directly from the input FASTA, using the
full-genome genotype and subtype references below, then writes one CSV per gene
in the requested output directory.

Genotype and subtype calls are made once from the full-genome references, then
reused in all three gene-specific output tables. The coverage result is also
mapped once and evaluated against each gene target range.

- Genotype references: `HCVData/Genotype-Ref/HCV_GT_FullGenome_RefSeqs.fasta`
- Subtype references: `HCVData/Subtype-Ref/HCV_Subtype_FullGenome_Refs.fasta`

```bash
.venv/bin/python Preprocess/Ref/hcv-all-seq-blast-subtyping-full-genome/scripts/audit_all_fasta_coverage.py \
  --input-fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --output-dir HCVData/nonComet-Full-genome \
  --threads 4
```

The output files are `NS3_AllSeq_NonComet_Coverage.csv`, `NS5A_AllSeq_NonComet_Coverage.csv`, and `NS5B_AllSeq_NonComet_Coverage.csv`, plus `Genotype_Distances.csv`, one `Subtype_Distances_Genotype_<genotype>.csv` file per genotype, and `Subtyping_Distances.xlsx`.

Each table has seven columns: `Accession`, `ClosestGenotype`, `ClosestGenotypePident`, `ClosestSubtype`, `ClosestSubtypePident`, `ReferenceOverlapAA`, and `FullyCover`. The percent-identity columns are the BLAST percent identities for the best full-genome genotype and genotype-matched subtype hits. Blank assignment fields mean the sequence did not meet the assignment threshold.

`Genotype_Distances.csv` has one row per accession, a distance column for every genotype reference, and the first and second genotype choices with their distances. Each genotype-specific subtype CSV has the same information for only that genotype's subtype references and assigned accessions. Distance is BLAST nucleotide percent-identity distance (`100 - pident`); blank values did not meet the minimum aligned-length threshold.

`Subtyping_Distances.xlsx` combines these distance tables: its `Genotype` sheet contains `Genotype_Distances.csv`, and each `Subtype_<genotype>` sheet contains the corresponding genotype-specific subtype table.

## Plot genotype-distance groups

Use the shared histogram script for `FirstChoiceDistance` and
`SecondChoiceDistance`. It writes figures to a `figures` subfolder beside the
input CSV. Use `--skip-empty` when processing all distance CSVs because some
subtype tables can have no second-choice values.

```bash
uv run python Preprocess/Ref/hcv-all-seq-blast-subtyping-per-gene/scripts/plot_distance_histogram.py \
  --input-csv HCVData/nonComet-Full-genome/Genotype_Distances.csv
uv run python Preprocess/Ref/hcv-all-seq-blast-subtyping-per-gene/scripts/plot_distance_histogram.py \
  --input-csv HCVData/nonComet-Full-genome/Genotype_Distances.csv \
  --column SecondChoiceDistance
```

`ReferenceOverlapAA` is blank when the best genotype alignment does not overlap the requested target range. When it overlaps, it reports the overlapping reference amino-acid interval, including partial overlap.

`FullyCover` is `Yes` only when the best genotype alignment spans the whole requested target range; it is blank for partial or absent overlap.

The script displays a live stage-level progress bar. Its BLAST searches use four workers by default; change `--threads` only when the available CPU capacity requires it.

The default `--min-aligned-nt` is 100. Use the option to set a different
minimum aligned nucleotide length for genotype and subtype calls.

Keep the default targets unless explicitly requested otherwise:

- NS3: amino-acid positions 36-175
- NS5A: amino-acid positions 26-93
- NS5B: amino-acid positions 150-321
