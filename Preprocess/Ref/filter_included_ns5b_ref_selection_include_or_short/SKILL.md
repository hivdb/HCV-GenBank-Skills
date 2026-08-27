---
name: filter-included-ns5b-ref-selection-include-or-short
description: Build an NS5B reference-selection workbook that retains Status values beginning with include (case-insensitive) or Short (case-sensitive).
---

# NS5B include-or-Short reference selection

Use `filter_included_ns5b_ref_selection_include_or_short.py` to create a new NS5B selection workbook from the `Original_NS5B` sheet.

- Keep statuses starting with `include`, ignoring case.
- Also keep statuses starting exactly with `Short`; do not treat `short` or `SHORT` as matches.
- The default output is `HCVData/Ref-selection/IncludedNS5BRefs_StatusIncludeOrShort.xlsx`; do not overwrite the existing include-only workbook.

Run the script with the project environment. Use `--input-xlsx` or `--output-xlsx` only when a different source or destination is explicitly required.
