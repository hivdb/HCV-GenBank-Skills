# Step 4: Filter full-genome coverage by RAS-overlap accessions

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/filter_full_genome_by_ras_accessions.py
```

For NS3, NS5A, and NS5B, this step uses accessions retained by steps 1–3 to
filter the matching `HCVData/nonComet-Full-genome/<gene>_AllSeq_NonComet_Coverage.csv`
file. It writes one filtered full-genome coverage CSV per gene to
`HCVData/subtyping-comparison-all-ras/04-filter-full-genome-by-ras-accessions/` and
prints the retained accession count for each gene.
