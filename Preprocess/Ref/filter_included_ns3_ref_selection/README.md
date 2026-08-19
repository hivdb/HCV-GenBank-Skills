# filter_included_ns3_ref_selection

Filters `HCVData/HCV_BlastHists_202604_data.xlsx`, sheet `Original_NS3`, to rows whose `Status` value begins with `include`, ignoring case.

It writes `HCVData/Ref-selection/IncludedNS3Refs_StatusInclude.xlsx` with the filtered rows in the `Original_NS3_Included` sheet.

```bash
.venv/bin/python Preprocess/Ref/filter_included_ns3_ref_selection/filter_included_ns3_ref_selection.py
```
