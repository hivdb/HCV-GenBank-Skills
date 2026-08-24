---
name: hcv-ns5a-comet-build-workflow
description: Use this skill when the user wants to run or inspect the HCV NS5A build scripts that discover RefID FASTA files, create genotype/subtype study workbooks, source-feature summaries, complete profile workbooks, and genotype/subtype RAS profile reports.
---

# HCV NS5A Build Workflow

Use this skill for the full NS5A high-throughput build workflow. The first step reads the configured Excel worksheet, discovers matching RefID FASTA files, and stages those files for downstream NS5A build steps.

## Workflow Chart

See `NS5A_workflow.svg` in this skill folder.
See `scripts/comet_ns5a_workflow_steps.csv` for the named Python-runner steps and their function scripts.

## Workflow Steps and Python Calls

Run the workflow in this order. For an individual step, invoke its listed script with `"$PYTHON_BIN" scripts/<script>.py ...`; use `scripts/run_ns5a_pipeline.py` as the source of truth for its arguments and paths. Steps with multiple listed calls require every call, in order.

1. Load configured defaults.
   - Python: `load_pipeline_defaults.py ns5a_comet <pipeline.local.toml> <repo-root>`
2. Prepare the output and temporary directories.
   - Python: `prepare_ns5a_pipeline_workdirs.py`
3. Discover RefID FASTA files.
   - Python: `find_refid_fastas.py`
4. Stage the matched FASTA files.
   - Python: `stage_matched_refid_fastas.py`
5. Create COMET assignment files and remove missing or unassigned records.
   - Python: `prepare_comet_ns5a_assignments.py`
6. Select priority non-COMET assignments for later genotype and subtype workbooks.
   - Python: `select_noncomet_priority_assignments.py`
7. Filter master accession metadata to the COMET-filtered FASTA accessions.
   - Python: `filter_accessions_metadata_by_fasta.py`
8. Apply RefID-specific metadata rules.
   - Python: `split_refid_metadata_csv.py`
9. Filter COMET-filtered FASTA records using the RefID metadata files.
   - Python: `filter_refid_fastas_by_metadata.py`
10. Build the COMET genotype study workbook.
   - Python: `build_ns5a_comet_gt_allstudies.py`
   - Python: `add_gt_counts_sheet.py`
11. Build the COMET subtype study workbook.
   - Python: `build_ns5a_comet_subtype_allstudies.py`
12. Extract genotype-position amino-acid sequences from the selected FASTA pool.
   - Python: `build_ns5a_subtype_with_gt_aa.py`
13. Validate profile alignment coordinates.
   - Python: `validate_ns5a_profile_alignment.py`
14. Summarize QC-passed genotype mutation burden and calculate profile-input counts.
   - Python: `build_qc_passed_genotype_mutation_burden_summary.py`
   - Python: `build_ns5a_completeprofiles_tabspergt.py --report-only`
15. Build complete genotype and subtype profile workbooks, then identify priority non-COMET profile accessions.
   - Python: `build_ns5a_completeprofiles_tabspergt.py`
   - Python: `export_noncomet_priority_profile_accessions.py`
16. Export genotype and subtype consensus FASTA files.
   - Python: `export_ns5a_consensus_fasta.py`
17. Align subtype consensus sequences to the fixed GT1_1a coordinate system.
   - Python: `align_ns5a_subtype_consensuses_to_gt1a.py`
18. Compare COMET consensus sequences with genotype and subtype references.
   - Python: `export_gt_reference_consensus_differences.py`
   - Python: `build_ns5a_subtype_consensus_reference_distance.py`
19. Build genotype and subtype RAS profiles.
   - Python: `build_ns5a_gt_ras_profiles.py`
   - Python: `build_ns5a_subtype_ras_profiles.py`
20. Build the combined RAS profile, update its coverage labels, and create its COMET coverage and sequence-audit reports.
   - Python: `build_ns5a_combined_ras_profiles.py`
   - Python: `replace_comet_profile_coverage_range_with_mean_diff.py`
   - Python: `build_comet_subtype_ras_coverage_report.py`
   - Python: `build_comet_workflow_sequence_audit.py`
21. Summarize subtype RAS differences from genotype consensus.
   - Python: `build_ns5a_subtype_ras_consensus_difference_summary.py`
22. Build the genotype amino-acid consensus distance matrix.
   - Python: `build_ns5a_gt_aa_distance_matrix.py`
23. Build subtype amino-acid distance matrices.
   - Python: `build_ns5a_subtype_aa_distance_matrices.py`
24. Build paired genotype and subtype AA/NA distance matrices for RAS positions and positions 24–93.
   - Python: `build_ns5a_aa_distance_matrices.py` (run once for each position set)
   - Python: `build_ns5a_na_distance_matrices.py` (run once for each position set)
25. Build genotype and subtype RAS entropy reports.
   - Python: `build_ns5a_ras_entropy.py`
26. Create an annotated combined RAS profile with `MeanDiff` and `PositionDiff` values.
   - Python: `add_combined_profile_nonconsensus_row.py`
27. Publish the shared ICTV reference/consensus comparison report.
   - Python: `add_subtype_consensus_mutation_summaries.py`

Run the Python orchestrator for full workflows or a named stage:

```bash
python hcv-workflow/hcv-ns5a-comet-build-workflow/scripts/run_ns5a_pipeline.py --list-steps
python hcv-workflow/hcv-ns5a-comet-build-workflow/scripts/run_ns5a_pipeline.py
python hcv-workflow/hcv-ns5a-comet-build-workflow/scripts/run_ns5a_pipeline.py --step discover-refid-fastas
```

`--list-steps` prints all available named stages. `--step <name>` runs only that named stage; repeat `--step` to run selected stages in the order supplied. Each selected stage expects its prior outputs to already exist.

Configuration stays in the repository base folder. The Python runner loads:

1. `.env`
2. `pipeline.local.toml`
3. built-in fallbacks

Explicit environment variables provided by the caller take precedence over `pipeline.local.toml`.
The TOML loader is bundled at `load_pipeline_defaults/load_pipeline_defaults.py` and is called with the explicit root config path.
Set `sheet_name` in the `[ns5a]` section of `pipeline.local.toml` to choose the input worksheet for discovery and genotype assignment.
Each step writes its files under `outputs/comet-NS5A/<order>_<step-name>/` (for example, `09_build-genotype-workbook/`). The runner prints each step as a Markdown heading with its order, name, and short explanation, followed by input, excluded, and final-included accession counts.

## Inputs

- Excel workbook and configured worksheet containing `RefID`, `RefName`, and patient-count fields
- FASTA pool directory containing RefID-prefixed FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- `outputs/folder_assignments/NS5A_assignments.csv` for local-assignment fallback, including subtype 1d, genotype 7/8, and COMET-unassigned accessions
- `Accessions_metadata.csv` for filtering metadata to accessions present in included FASTA files
- optional GenBank directory for source-feature extraction if the commented source-feature steps are re-enabled

The discovery step keeps every row with a `RefID` from the configured `IncludedNS5ARefs_StatusInclude.xlsx` selection workbook; it does not apply a `Num Pts` filter.
After discovery, the runner copies matched FASTA files into `03_stage-refid-fastas/included_refid_fastas/`. Step 4 copies them into `04_prepare-comet-assignments/included_refid_fastas/` and removes records missing from COMET or marked unassigned. Steps 6–8 apply metadata filtering to that COMET-filtered copy; downstream workbooks use the filtered copy under `08_filter-refid-fastas/`.
After every stage following staging, the runner prints the current number of unique GT7 and GT8 subtypes and their accession counts. Before the subtype workbook exists, these counts are calculated from the same COMET-plus-priority assignments it will use.
The metadata filtering step writes its metadata CSVs in `04_filter-accession-metadata/`; the per-RefID rules and CSVs are in `05_split-refid-metadata/refid_metadata/`.
The priority-assignment stage selects local calls for retained accessions called `1d`, for genotype 7 or 8 accessions, and whenever COMET marks the subtype unassigned. The following genotype and subtype steps consume that single selection to override COMET calls or add accessions absent from COMET. Amino-acid extraction reads the configured FASTA pool for added accessions.

## Outputs

The workflow writes NS5A outputs under `outputs/comet-NS5A/`, in numbered step folders, including:

- `NS5A_GT_AllStudies.xlsx`
- `NS5A_matched_fasta_files.txt`
- discovery `filtered_rows.xlsx` under `02_discover-refid-fastas/`
- copied included RefID FASTA files under `03_stage-refid-fastas/`, the COMET-filtered copy under `04_prepare-comet-assignments/`, and the metadata-filtered copy and `kept_accessions.csv` under `08_filter-refid-fastas/`
- `NS5A_NonComet_Priority_Assignments.csv` under `05_select-noncomet-priority-assignments/`
- `included_accessions_metadata.csv`
- `missing_accessions_from_metadata.txt`
- `refid_metadata/RefID_<RefID>_metadata.csv`
- filtered copied RefID FASTA files in `included_refid_fastas/` for RefIDs with metadata filters
- source-feature CSV/XLSX outputs only if the commented source-feature steps are re-enabled
- `NS5A_Subtype_AllStudies_WSeqs.xlsx`
- `NS5A_Subtype_With_GT_AA.xlsx`
- `NS5A_Profile_Input_Alignment_QC.xlsx` (profile input with per-accession alignment QC columns)
- `NS5A_QC_Passed_Genotype_Mutation_Burden_Summary.csv` (per-genotype mutation burden among QC-passed input rows)
- `NS5A_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS5A_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS5A_Subtype_CompleteProfiles_Merged.xlsx` (one merged subtype table with `Subtype`, `NS5APosition`, `NumSeqsIncludingPosition`, `AminoAcid`, `CountWithAA`, and `PctWithAA`)
- `NS5A_GT_Consensus.fasta`
- `NS5A_Subtype_Consensus.fasta`
- `NS5A_Subtype_Consensus_Aligned_to_GT1_1a.fasta` (all subtype consensuses aligned to GT1_1a coordinates)
- `NS5A_GT_RAS_Profiles.xlsx`
- `NS5A_Subtype_RAS_Profiles.xlsx`
- `NS5A_Subtype_RAS_Profiles_Explicit_AA.xlsx` (all reportable subtype amino acids at RAS positions; used by the ICTV publication step)
- `NS5A_Combined_RAS_Profiles.xlsx` (base combined profile)
- `NS5A_Combined_RAS_Profiles_Annotated.xlsx` (derived profile with `MeanDiff`, `PositionDiff`, and updated coverage labels)
- paired AA/NA RAS and position-range distance workbooks under `28_build-paired-distance-matrices/`

## Operating Rules

- Keep NS5A scripts together in this skill folder.
- Use `scripts/run_ns5a_pipeline.py` for complete runs or `--step <name>` for one specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Keep every generated artifact in its numbered step folder under `outputs/comet-NS5A/`.
- Preserve the order above because later reports consume earlier workbooks.
