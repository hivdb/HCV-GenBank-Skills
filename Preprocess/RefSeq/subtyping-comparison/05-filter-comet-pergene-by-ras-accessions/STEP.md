# Step 5: Filter COMET per-gene calls by RAS-overlap accessions

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/filter_comet_pergene_by_ras_accessions.py
```

For NS3, NS5A, and NS5B, this step extracts the accession from the COMET
`name` value and retains calls whose accession was kept by steps 1–3. It writes
one filtered COMET CSV per gene to
`HCVData/subtyping-comparison/05-filter-comet-pergene-by-ras-accessions/` and
prints the retained unique accession count for each gene.
