# Step 9: Count blank values in merged subtyping results

Run from the repository root:

```bash
uv run python Preprocess/RefSeq/subtyping-comparison/scripts/count_merged_blank_values.py
```

For each NS3, NS5A, and NS5B Step 8 merged result, this step writes a separate
`<gene>_Blank_Value_Counts.csv` file with the blank-value count and total row
count for every column. A value is blank when it is empty or contains only
whitespace.
