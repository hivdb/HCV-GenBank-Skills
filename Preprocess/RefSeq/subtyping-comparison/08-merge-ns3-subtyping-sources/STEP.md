# Step 8: Merge subtyping sources

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/merge_ns3_subtyping_sources.py
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/merge_ns3_subtyping_sources.py --gene NS5A
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/merge_ns3_subtyping_sources.py --gene NS5B
```

This step outer-merges the RAS-overlap outputs from per-gene BLAST,
full-genome BLAST, per-gene COMET, full-genome COMET, and GenBank. The output
has one accession per row and paired genotype/subtype columns for each source.
COMET genotypes are derived from the numeric prefix of the COMET subtype.
The script also prints, for each source, the number of merged accessions absent
from that source file. Before writing the merged table, it writes
`<gene>_Missing_Accessions_By_Source.csv` for accessions absent from one or more
sources.
