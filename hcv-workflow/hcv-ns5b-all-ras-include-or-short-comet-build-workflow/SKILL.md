---
name: hcv-ns5b-all-ras-include-or-short-comet-build-workflow
description: Use this skill for HCV NS5B all-RAS COMET builds using the Include-or-Short reference-selection workbook.
---

# HCV NS5B All-RAS Include-or-Short COMET Build Workflow

Use this skill for the full NS5B high-throughput build workflow. The first step reads the configured Excel worksheet, discovers matching RefID FASTA files, and stages those files for downstream NS5B build steps.

## Script Order

1. `find_refid_fastas/find_refid_fastas.py`
2. copy matched FASTA files to `included_refid_fastas/`
3. `filter_accessions_metadata_by_fasta/filter_accessions_metadata_by_fasta.py`
4. `split_refid_metadata_csv/split_refid_metadata_csv.py`
5. `filter_refid_fastas_by_metadata/filter_refid_fastas_by_metadata.py`
6. `build_ns5b_gt_allstudies/build_ns5b_gt_allstudies.py`
7. `build_ns5b_subtype_allstudies_wseqs/build_ns5b_subtype_allstudies_wseqs.py`
8. `build_ns5b_subtype_with_gt_aa/build_ns5b_subtype_with_gt_aa.py`
9. `validate_ns5b_profile_alignment/validate_ns5b_profile_alignment.py`
10. `build_ns5b_completeprofiles_tabspergt/build_ns5b_completeprofiles_tabspergt.py`
11. `export_ns5b_consensus_fasta/export_ns5b_consensus_fasta.py`
12. `align_ns5b_subtype_consensuses_to_gt1a/align_ns5b_subtype_consensuses_to_gt1a.py`
13. `export_gt_reference_consensus_differences/export_gt_reference_consensus_differences.py`
14. `build_ns5b_gt_ras_profiles/build_ns5b_gt_ras_profiles.py`
15. `build_ns5b_subtype_ras_profiles/build_ns5b_subtype_ras_profiles.py`
16. Build the combined RAS profile and update its coverage labels.
    - `build_ns5b_combined_ras_profiles/build_ns5b_combined_ras_profiles.py`
    - `replace_comet_profile_coverage_range_with_mean_diff/replace_comet_profile_coverage_range_with_mean_diff.py`
17. Build a separate combined RAS profile using only accessions with at least 90% non-`X` coverage across NS5B positions 150–321.

Use the Python orchestrator for complete runs or selected, resumable stages:

```bash
python hcv-workflow/hcv-ns5b-all-ras-include-or-short-comet-build-workflow/scripts/run_ns5b_pipeline.py --list-steps
python hcv-workflow/hcv-ns5b-all-ras-include-or-short-comet-build-workflow/scripts/run_ns5b_pipeline.py
python hcv-workflow/hcv-ns5b-all-ras-include-or-short-comet-build-workflow/scripts/run_ns5b_pipeline.py --step discover-refid-fastas
```

`--list-steps` prints all named stages. `--step <name>` runs only that stage; repeat the flag to run selected stages in order. Selected stages expect their prerequisite outputs to exist. The legacy shell wrapper remains available for compatibility but is no longer the preferred entry point.

Configuration stays in the repository base folder. The wrapper loads:

1. `.env`
2. `pipeline.local.toml`
3. built-in fallbacks

Explicit environment variables provided by the caller take precedence over `pipeline.local.toml`.
The TOML loader is bundled at `load_pipeline_defaults/load_pipeline_defaults.py` and is called with the explicit root config path.
The variant inherits settings from `[common]` and `[ns5b_comet]`. Its `[ns5b_all_ras_include_or_short]` configuration selects `HCVData/Ref-selection/IncludedNS5BRefs_StatusIncludeOrShort.xlsx` and its isolated output directory.
Each stage writes its outputs under a numbered directory in `outputs/comet-NS5B-all-ras-include-or-short/`.
Step 11 caches RefID FASTAs and parallelizes accession preparation/AA extraction and BLASTX with four workers by default; set `NS5B_WORKERS` or pass `--workers` to tune this.

## Inputs

- Excel workbook and configured worksheet containing `RefID`
- FASTA pool directory containing RefID-prefixed FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- `HCVData/HCV-all-seq-subtype/NS5B_AllSeq_NonComet_Coverage.csv` for mandatory non-COMET subtype 1d and genotype 7/8 overrides and additions
- `Accessions_metadata.csv` for filtering metadata to accessions present in included FASTA files

The discovery step keeps every row with a non-empty `RefID`. It does not require or filter on patient-count columns, `NS5BCount`, or `Notes`.
After discovery, the runner copies all matched RefID FASTA files into `03_stage-refid-fastas/included_refid_fastas/`. Step 4 copies them into `04_prepare-comet-assignments/included_refid_fastas/` and removes records missing from COMET or marked unassigned. Steps 6–8 apply metadata filtering to that COMET-filtered copy; downstream workbooks use the filtered copy under `08_filter-refid-fastas/`.
After every stage following staging, the runner prints the current number of unique GT7 and GT8 subtypes and their accession counts. Before the subtype workbook exists, these counts are calculated from the same COMET-plus-priority assignments it will use.
The metadata filtering step writes `included_accessions_metadata.csv` and reports any FASTA accessions missing from `Accessions_metadata.csv` in `missing_accessions_from_metadata.txt`; both files live in the parent folder of `included_refid_fastas/`.
The per-RefID metadata split step writes CSVs only for RefIDs that have explicit filters under `refid_metadata/`. Current filters: `30` source_isolate contains `day1`; `142` accession is listed in `HCVData/Ref-selection/NS5_Ref_filter/NS5B/142.csv`; `346` source_isolate contains `baseline`; `891` source_isolate contains a token from `Ha01` through `Ha97`; and `943` source_isolate contains `day 1`.
The per-RefID FASTA filtering step reads `refid_metadata/RefID_<RefID>_metadata.csv`, keeps only matching `Accession` records in the corresponding copied FASTA file under `included_refid_fastas/`, and prints per-RefID and total before/after record counts.
The priority-assignment stage uses `HCVData/Reference_seqs/HCV_Subtype_Refs_AA_Accession_Subtype.csv` as the highest-priority accession/subtype source, then selects non-COMET calls for retained accessions called `1d` and for genotype 7 or 8 accessions. The following genotype and subtype steps consume that single selection to override COMET calls or add accessions absent from COMET. Amino-acid extraction reads the configured FASTA pool for added accessions.

Complete-profile construction retains an accession only when it has callable amino-acid coverage at every NS5B RAS position: 150, 159, 206, 282, 316, 320, and 321. Missing, `X`, stop (`*`), or non-standard calls at any of these positions exclude that accession from the profile and its downstream RAS and distance reports.
The `build-90pct-range-coverage-combined-profile` step is separate from the standard combined profile. It applies the same RAS requirements plus at least 90% non-`X` coverage across positions 150–321, and writes its own profile-accession list and intermediate profiles.

## Outputs

The final GT7/GT8 local-assignment comparison step reads `outputs/folder_assignments/NS5B_assignments.csv` and writes the workflow/local subtype comparison workbook and CSV in its own numbered output directory.

The workflow writes NS5B outputs under `outputs/comet-NS5B-all-ras/`, including:

- `NS5B_GT_AllStudies.xlsx`
- `NS5B_matched_fasta_files.txt`
- discovery `filtered_rows.xlsx` under `outputs/comet-NS5B-all-ras/temp/.../find_refid_fastas/...`
- copied included RefID FASTA files under `03_stage-refid-fastas/`, the COMET-filtered copy under `04_prepare-comet-assignments/`, and the metadata-filtered copy and `kept_accessions.csv` under `08_filter-refid-fastas/`
- `NS5B_NonComet_Priority_Assignments.csv` under `05_select-noncomet-priority-assignments/`
- `included_accessions_metadata.csv`
- `missing_accessions_from_metadata.txt`
- `refid_metadata/RefID_<RefID>_metadata.csv`
- filtered copied RefID FASTA files in `included_refid_fastas/` for RefIDs with metadata filters
- `NS5B_Subtype_AllStudies_WSeqs.xlsx`
- `NS5B_Subtype_With_GT_AA.xlsx`
- `NS5B_Profile_Input_Alignment_QC.xlsx` (profile input with per-accession alignment QC columns)
- `NS5B_Profile_Accessions_QC_Pass.csv` (accessions retained for complete-profile construction)
- `NS5B_QC_Passed_Genotype_Mutation_Burden_Summary.csv` (per-genotype mutation burden among QC-passed input rows)
- `NS5B_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS5B_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS5B_Subtype_CompleteProfiles_Merged.xlsx` (one merged subtype table with `Subtype`, `NS5BPosition`, `NumSeqsIncludingPosition`, `AminoAcid`, `CountWithAA`, and `PctWithAA`)
- `NS5B_GT_Consensus.fasta`
- `NS5B_Subtype_Consensus.fasta`
- `NS5B_Subtype_Consensus_Aligned_to_GT1_1a.fasta` (all subtype consensuses aligned to GT1_1a coordinates)
- `NS5B_GT_RAS_Profiles.xlsx`
- `NS5B_Subtype_RAS_Profiles.xlsx`
- `NS5B_Subtype_RAS_Profiles_Explicit_AA.xlsx` (all reportable subtype amino acids at RAS positions; used by the ICTV publication step)
- `NS5B_Combined_RAS_Profiles_90Pct_Range_Coverage.xlsx` under `26_build-90pct-range-coverage-combined-profile/` (combined profile restricted to at least 90% coverage of positions 150–321)
- paired AA/NA RAS and position-range distance workbooks under `outputs/`

## Operating Rules

- Keep NS5B scripts together in this skill folder.
- Use `scripts/run_ns5b_pipeline.py` for complete runs or `--step <name>` for a specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Keep temporary outputs under `outputs/comet-NS5B-all-ras/temp/` so they do not mix with other skills.
- Preserve the order above because later reports consume earlier workbooks.
