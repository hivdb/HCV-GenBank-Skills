# Step 11: Report subtype agreement

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/report_subtype_agreement.py
```

For each Step 8 merged result, this step compares the subtype calls from
per-gene BLAST, full-sequence BLAST, per-gene COMET, full-sequence COMET, and
GenBank. It writes accession-level and summary reports per gene. The
accession-level report includes `PresentMethodCount` and `BlankMethods`; the
summary splits non-unanimous calls into `MoreThanHalfAgree` (one call has
majority support) and `MoreThanHalfDisagree` (no call has majority support),
with plain-language agreement-status descriptions.
All six report CSVs are also merged into `Subtype_Agreement_Reports.xlsx`,
with one worksheet per CSV.

Each gene's Step 11 accession list is also used to filter
`HCVData/nonComet-Full-genome/Subtyping_Distances.xlsx`. The matching rows are
saved separately as `<GENE>_LocalFullSeq_Subtyping_Distances.xlsx` in this
output directory.
For each gene, additional filtered workbooks are created for COMET full-sequence
subtypes `1L`, `3B`, `3H`, and `4R`.
