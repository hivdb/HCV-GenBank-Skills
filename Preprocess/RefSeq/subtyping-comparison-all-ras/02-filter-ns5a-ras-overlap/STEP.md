# Step 2: Filter NS5A coverage by RAS overlap

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/filter_ras_overlap.py \
  --input-csv HCVData/nonComet-PerGene/NS5A_AllSeq_NonComet_Coverage.csv \
  --profile-accessions-csv outputs/comet-NS5A-all-ras/16_build-complete-profiles/NS5A_Profile_Accessions_QC_Pass.csv \
  --output-csv HCVData/subtyping-comparison-all-ras/02-filter-ns5a-ras-overlap/NS5A_AllSeq_NonComet_Coverage_RAS_Overlap.csv \
  --ras-positions 24,26,28,29,30,31,32,38,58,62,92,93
```

The step retains rows whose `ReferenceOverlapAA` covers at least one NS5A RAS
position and whose accession is in the NS5A QC-passed profile accession list,
then prints the selected-list and retained unique accession counts.
