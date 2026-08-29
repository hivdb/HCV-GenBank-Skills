# Step 3: Filter NS5B coverage by RAS overlap

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/filter_ras_overlap.py \
  --input-csv HCVData/nonComet-PerGene/NS5B_AllSeq_NonComet_Coverage.csv \
  --staged-accessions-csv outputs/comet-NS5B-position-282/03_stage-refid-fastas/staged_accessions.csv \
  --output-csv HCVData/subtyping-comparison/03-filter-ns5b-ras-overlap/NS5B_AllSeq_NonComet_Coverage_RAS_Overlap.csv \
  --ras-positions 150,159,206,282,316,320,321
```

The step retains rows whose `ReferenceOverlapAA` covers at least one NS5B RAS
position and whose accession is in the NS5B staged-accession list, then prints
the staged-list and retained unique accession counts.
