---
name: hcv-all-accessions-genbank-subtypes
description: Extract HCV genotype and subtype metadata for accessions in a FASTA file from a local GenBank flatfile archive, producing a three-column CSV.
---

# HCV All Accessions GenBank Subtypes

Use this skill when an HCV FASTA accession list must be matched to local GenBank
flatfiles and represented as `accession`, `genotype`, and `subtype` in a CSV.

## Run

Use the bundled Biopython script. It preserves FASTA accession order, removes a
version suffix before matching GenBank records, and leaves metadata blank when
the record or a usable genotype/subtype annotation is absent.

```bash
uv run python Preprocess/RefSeq/hcv-all-accessions-genbank-subtypes/scripts/extract_genbank_subtypes.py \
  --fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --genbank-dir /path/to/GenBankFiles \
  --output-csv HCVData/Genbank-subtype/all_accessions_genotypes_subtypes.csv
```

The archive is scanned as text to avoid constructing Biopython objects for
unrelated records. Every selected record is parsed with Biopython. Genotype and
subtype are collected from source-feature qualifiers (including notes); a value
such as `4n` is represented as genotype `4` and subtype `4n`.
