# Step 7: Filter GenBank subtypes by RAS-overlap accessions

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/filter_genbank_subtypes_by_ras_accessions.py
```

This step filters `HCVData/Genbank-subtype/all_accessions_genotypes_subtypes.csv`
separately with the NS3, NS5A, and NS5B accession lists retained by steps 1–3.
It writes one CSV per gene to
`HCVData/subtyping-comparison/07-filter-genbank-subtypes-by-ras-accessions/`
and prints the retained unique accession count for each gene.
