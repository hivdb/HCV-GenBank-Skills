# Step 3: Filter NS5B coverage by RAS overlap

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/filter_ras_overlap.py \
  --input-csv HCVData/nonComet-PerGene/NS5B_AllSeq_NonComet_Coverage.csv \
  --profile-accessions-csv outputs/comet-NS5B-all-ras/16_build-complete-profiles/NS5B_Profile_Accessions_QC_Pass.csv \
  --output-csv HCVData/subtyping-comparison-all-ras/03-filter-ns5b-ras-overlap/NS5B_AllSeq_NonComet_Coverage_RAS_Overlap.csv \
  --ras-positions 150,159,206,282,316,320,321
```

The step retains rows whose `ReferenceOverlapAA` covers at least one NS5B RAS
position and whose accession is in the NS5B QC-passed profile accession list,
then prints the selected-list and retained unique accession counts.
