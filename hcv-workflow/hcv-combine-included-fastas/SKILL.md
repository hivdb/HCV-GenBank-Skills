---
name: hcv-combine-included-fastas
description: Combine the staged included_refid_fastas files from the HCV NS3, NS5A, and NS5B workflows into one gene-named FASTA per workflow. Use when a user needs a combined NS3.fasta, NS5A.fasta, or NS5B.fasta from the current filtered study FASTA folders.
---

# Combine HCV Included FASTAs

Run the bundled script to concatenate the current filtered study FASTAs. It writes `NS3.fasta`, `NS5A.fasta`, and `NS5B.fasta` in `assets/` by default.

```bash
python hcv-workflow/hcv-combine-included-fastas/combine_included_fastas/combine_included_fastas.py
```

The defaults read these folders under the repository root:

- `outputs/comet-NS3/temp/run_ns3_pipeline/included_refid_fastas`
- `outputs/comet-NS5A/temp/run_ns5a_pipeline/included_refid_fastas`
- `outputs/comet-NS5B-all-ras/08_filter-refid-fastas/included_refid_fastas`

The script accepts `--ns3-dir`, `--ns5a-dir`, `--ns5b-dir`, and `--output-dir` when a run uses different locations. It concatenates sorted `*.fasta` files without changing FASTA headers or sequences, and reports per-gene input-file and record counts.
