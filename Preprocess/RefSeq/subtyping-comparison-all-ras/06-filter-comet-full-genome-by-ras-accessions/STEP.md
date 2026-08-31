# Step 6: Filter full-genome COMET calls by RAS-overlap accessions

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/filter_comet_full_genome_by_ras_accessions.py
```

This step filters `HCVData/Comet-Full-genome/all_comet_subtype.csv` separately
with the NS3, NS5A, and NS5B accession lists retained by steps 1–3. It writes
one filtered full-genome COMET CSV per gene to
`HCVData/subtyping-comparison-all-ras/06-filter-comet-full-genome-by-ras-accessions/`
and prints the retained unique accession count for each gene.
