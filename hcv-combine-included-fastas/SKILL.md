---
name: hcv-combine-included-fastas
description: Combine the staged included_refid_fastas files from the HCV NS3, NS5A, and NS5B workflows into one gene-named FASTA per workflow. Use when a user needs a combined NS3.fasta, NS5A.fasta, or NS5B.fasta from the current filtered study FASTA folders.
---

# Combine HCV Included FASTAs

Run the bundled script to concatenate the current filtered study FASTAs. It writes `NS3.fasta`, `NS5A.fasta`, and `NS5B.fasta` in `assets/` by default.

```bash
python hcv-combine-included-fastas/scripts/combine_included_fastas.py
```

The defaults read these folders under the repository root:

- `temp/hcv-ns3-build-workflow/run_ns3_pipeline/included_refid_fastas`
- `temp/hcv-ns5a-build-workflow/run_ns5a_pipeline/included_refid_fastas`
- `temp/hcv-ns5b-build-workflow/run_ns5b_pipeline/included_refid_fastas`

The script accepts `--ns3-dir`, `--ns5a-dir`, `--ns5b-dir`, and `--output-dir` when a run uses different locations. It concatenates sorted `*.fasta` files without changing FASTA headers or sequences, and reports per-gene input-file and record counts.
