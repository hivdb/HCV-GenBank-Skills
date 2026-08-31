# Step 11: Report subtype agreement

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/report_subtype_agreement.py
```

For each Step 8 merged result, this step compares the subtype calls from
per-gene BLAST, full-genome BLAST, per-gene COMET, full-genome COMET, and
GenBank. It writes accession-level and summary reports per gene. The
accession-level report includes `PresentMethodCount` and `BlankMethods`; the
summary splits non-unanimous calls into `MoreThanHalfAgree` (one call has
majority support) and `MoreThanHalfDisagree` (no call has majority support),
with plain-language agreement-status descriptions.
All six report CSVs are also merged into `Subtype_Agreement_Reports.xlsx`,
with one worksheet per CSV.
