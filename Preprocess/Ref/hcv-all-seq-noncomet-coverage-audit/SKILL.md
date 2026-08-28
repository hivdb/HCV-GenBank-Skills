---
name: hcv-all-seq-noncomet-coverage-audit
description: Assign non-COMET HCV genotype and closest within-genotype subtype for every record in an all.fasta-style nucleotide FASTA, then report overlap with NS3 positions 36-175, NS5A positions 26-93, and NS5B positions 150-321. Use when an accession-level coverage table is needed from an HCV sequence collection without relying on COMET calls.
---

# HCV all-sequence non-COMET coverage audit

Run the bundled script from the repository root. It performs genotype and genotype-matched subtype assignment directly from the input FASTA, then writes one CSV per gene in the requested output directory.

```bash
.venv/bin/python Preprocess/Ref/hcv-all-seq-noncomet-coverage-audit/scripts/audit_all_fasta_coverage.py \
  --input-fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --output-dir outputs/local_alignment \
  --threads 4
```

The output files are `NS3_AllSeq_NonComet_Coverage.csv`, `NS5A_AllSeq_NonComet_Coverage.csv`, and `NS5B_AllSeq_NonComet_Coverage.csv`.

Each table has seven columns: `Accession`, `ClosestGenotype`, `ClosestGenotypePident`, `ClosestSubtype`, `ClosestSubtypePident`, `ReferenceOverlapAA`, and `FullyCover`. The percent-identity columns are the BLAST percent identities for the best genotype and genotype-matched subtype hits. Blank assignment fields mean the sequence did not meet the assignment threshold for that gene.

`ReferenceOverlapAA` is blank when the best genotype alignment does not overlap the requested target range. When it overlaps, it reports the overlapping reference amino-acid interval, including partial overlap.

`FullyCover` is `Yes` only when the best genotype alignment spans the whole requested target range; it is blank for partial or absent overlap.

The script displays a live stage-level progress bar. Its BLAST searches use four workers by default; change `--threads` only when the available CPU capacity requires it.

Keep the default targets unless explicitly requested otherwise:

- NS3: amino-acid positions 36-175
- NS5A: amino-acid positions 26-93
- NS5B: amino-acid positions 150-321
