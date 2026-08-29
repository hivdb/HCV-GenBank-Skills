# Step 10: Report genotype agreement

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/report_genotype_agreement.py
```

For each Step 8 merged result, this step compares the genotype calls from
per-gene BLAST, full-genome BLAST, per-gene COMET, full-genome COMET, and
GenBank. It writes an accession-level report and a summary report per gene.
Calls are considered in agreement when all nonblank genotype values for that
accession are identical. The report separately identifies complete five-method
agreement, disagreement, and accessions with fewer than two available calls.
