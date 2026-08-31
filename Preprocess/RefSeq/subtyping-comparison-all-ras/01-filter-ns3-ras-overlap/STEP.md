# Step 1: Filter NS3 coverage by RAS overlap

Run the filter from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/filter_ns3_ras_overlap.py
```

It retains rows from `HCVData/nonComet-PerGene/NS3_AllSeq_NonComet_Coverage.csv`
whose `ReferenceOverlapAA` interval covers at least one NS3 RAS position. It
also requires the accession to be listed in
`outputs/comet-NS3-all-ras/16_build-complete-profiles/NS3_Profile_Accessions_QC_Pass.csv`.
The output is written to
`HCVData/subtyping-comparison-all-ras/01-filter-ns3-ras-overlap/` and the script prints
the number of retained unique accessions.
