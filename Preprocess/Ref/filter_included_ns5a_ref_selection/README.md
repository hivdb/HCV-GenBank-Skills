# filter_included_ns5a_ref_selection

Filters `HCVData/HCV_BlastHists_202604_data.xlsx`, sheet `Original_NS5A`, to rows whose `Status` value begins with `include`, ignoring case.

It writes `HCVData/Ref-selection/IncludedNS5ARefs_StatusInclude.xlsx` with the filtered rows in the `Original_NS5A_Included` sheet.

```bash
.venv/bin/python Preprocess/Ref/filter_included_ns5a_ref_selection/filter_included_ns5a_ref_selection.py
```
