# filter_included_ns5b_ref_selection

Filters `HCVData/HCV_BlastHists_202604_data.xlsx`, sheet `Original_NS5B`, to rows whose `Status` value begins with `include`, ignoring case.

It writes `HCVData/Ref-selection/IncludedNS5BRefs_StatusInclude.xlsx` with the filtered rows in the `Original_NS5B_Included` sheet.

```bash
.venv/bin/python Preprocess/Ref/filter_included_ns5b_ref_selection/filter_included_ns5b_ref_selection.py
```
