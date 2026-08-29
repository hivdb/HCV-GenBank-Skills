# Step 1: Filter NS3 coverage by RAS overlap

Run the filter from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/filter_ns3_ras_overlap.py
```

It retains rows from `HCVData/nonComet-PerGene/NS3_AllSeq_NonComet_Coverage.csv`
whose `ReferenceOverlapAA` interval covers at least one NS3 RAS position. It
also requires the accession to be listed in
`outputs/comet-NS3-one-ras/03_stage-refid-fastas/staged_accessions.csv`.
The output is written to
`HCVData/subtyping-comparison/01-filter-ns3-ras-overlap/` and the script prints
the number of retained unique accessions.
