# Step 12: Report subtype disagreements by subtype

Run from the repository root:

```bash
python3 Preprocess/RefSeq/subtyping-comparison-all-ras/scripts/report_subtype_disagreements_by_subtype.py
```

For NS3, NS5A, and NS5B, this step groups non-unanimous subtype calls by
`CometFullSeqSubtype`. Each output CSV has `CometSubtype`,
`AllFiveDisagree`, `AvailableCallsDisagree`, `TotalDisagree`, and
`OtherMethodSubtypes` as its final column. `OtherMethodSubtypes` is the
deduplicated list of calls from per-gene BLAST, full-sequence BLAST, full-genome
COMET, and GenBank for the disagreement accessions in that COMET subtype.
The three output CSVs are merged into `Subtype_Disagreement_Reports.xlsx`,
with one worksheet per gene.
