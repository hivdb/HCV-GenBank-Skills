# All-RAS subtyping comparison outputs

This directory contains outputs from the `subtyping-comparison-all-ras` workflow.
It is isolated from `HCVData/subtyping-comparison`.

The workflow begins with each gene's COMET all-RAS
`*_Profile_Accessions_QC_Pass.csv`, then retains only per-gene BLAST
records that cover at least one relevant RAS position before comparing subtype
calls across per-gene BLAST, full-genome BLAST, per-gene COMET, full-genome
COMET, and GenBank.
