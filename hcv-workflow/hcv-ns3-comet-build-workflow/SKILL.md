---
name: hcv-ns3-comet-build-workflow
description: Use this skill when the user wants to run or inspect the HCV NS3 COMET build scripts that discover RefID FASTA files, create genotype/subtype study workbooks, complete profile workbooks, and genotype/subtype RAS profile reports.
---

# HCV NS3 Build Workflow

Use this skill for the full NS3 high-throughput build workflow. The first step reads the configured Excel worksheet, discovers matching RefID FASTA files, and stages those files for downstream NS3 build steps.

## Script Order

1. `find_refid_fastas/find_refid_fastas.py`
2. copy matched FASTA files to `included_refid_fastas/`
3. `filter_accessions_metadata_by_fasta/filter_accessions_metadata_by_fasta.py`
4. `split_refid_metadata_csv/split_refid_metadata_csv.py`
5. `filter_refid_fastas_by_metadata/filter_refid_fastas_by_metadata.py`
6. `build_ns3_gt_allstudies/build_ns3_gt_allstudies.py`
7. `build_ns3_subtype_allstudies_wseqs/build_ns3_subtype_allstudies_wseqs.py`
8. `build_ns3_subtype_with_gt_aa/build_ns3_subtype_with_gt_aa.py`
9. `validate_ns3_profile_alignment/validate_ns3_profile_alignment.py`
10. `build_ns3_completeprofiles_tabspergt/build_ns3_completeprofiles_tabspergt.py`
11. `export_ns3_consensus_fasta/export_ns3_consensus_fasta.py`
12. `align_ns3_subtype_consensuses_to_gt1a/align_ns3_subtype_consensuses_to_gt1a.py`
13. `export_gt_reference_consensus_differences/export_gt_reference_consensus_differences.py`
14. `build_ns3_gt_ras_profiles/build_ns3_gt_ras_profiles.py`
15. `build_ns3_subtype_ras_profiles/build_ns3_subtype_ras_profiles.py`
16. Build the combined RAS profile and update its coverage labels.
    - `build_ns3_combined_ras_profiles/build_ns3_combined_ras_profiles.py`
    - `replace_comet_profile_coverage_range_with_mean_diff/replace_comet_profile_coverage_range_with_mean_diff.py`
17. `add_combined_profile_nonconsensus_row/add_combined_profile_nonconsensus_row.py`
18. `build_ns3_subtype_ras_consensus_difference_summary/build_ns3_subtype_ras_consensus_difference_summary.py`
19. `build_ns3_subtype_profile_coverage_report/build_ns3_subtype_profile_coverage_report.py` (GT5 subtype 5a coverage audit)
20. `build_ns3_gt7_gt8_step_audit.py` (GT7/GT8 retention and exclusion audit compared with the previous workflow step)

Use the Python orchestrator for complete runs or selected, resumable stages:

```bash
python hcv-workflow/hcv-ns3-comet-build-workflow/scripts/run_ns3_pipeline.py --list-steps
python hcv-workflow/hcv-ns3-comet-build-workflow/scripts/run_ns3_pipeline.py
python hcv-workflow/hcv-ns3-comet-build-workflow/scripts/run_ns3_pipeline.py --step discover-refid-fastas
```

`--list-steps` prints all named stages. `--step <name>` runs only that stage; repeat the flag to run selected stages in order. Selected stages expect their prerequisite outputs to exist.

Configuration stays in the repository base folder. The wrapper loads:

1. `.env`
2. `pipeline.local.toml`
3. built-in fallbacks

Explicit environment variables provided by the caller take precedence over `pipeline.local.toml`.
The TOML loader is bundled at `load_pipeline_defaults/load_pipeline_defaults.py` and is called with the explicit root config path.
Set `sheet_name` in the `[ns3]` section of `pipeline.local.toml` to choose the input worksheet for discovery and genotype assignment.
Each stage writes its outputs under a numbered directory in `outputs/comet-NS3/`.

## Inputs

- Excel workbook and configured worksheet containing `RefID`
- FASTA pool directory containing RefID-prefixed FASTA files
- `HCV_GT_RefSeqs.fasta`
- `HCV_Subtype_Refs_By_Genome_NA.json`
- `HCV_GT_Refs_By_Gene_AA.json`
- `HCVData/HCV-all-seq-subtype/NS3_AllSeq_NonComet_Coverage.csv` for mandatory non-COMET subtype 1d and genotype 7/8 overrides and additions
- `Accessions_metadata.csv` for filtering metadata to accessions present in included FASTA files

The discovery step keeps every row with a non-empty `RefID`. It does not require or filter on `NumPatients`, `Num Pts`, `NS3Count`, or `Notes`.
After discovery, the runner copies all matched RefID FASTA files into `03_stage-refid-fastas/included_refid_fastas/`. Step 4 copies them into `04_prepare-comet-assignments/included_refid_fastas/` and removes records missing from COMET or marked unassigned. Steps 6–8 apply metadata filtering to that COMET-filtered copy; downstream workbooks use the filtered copy under `08_filter-refid-fastas/`.
After every stage following staging, the runner prints the current number of unique GT7 and GT8 subtypes and their accession counts. Before the subtype workbook exists, these counts are calculated from the same COMET-plus-priority assignments it will use.
The metadata filtering step writes `included_accessions_metadata.csv` and reports any FASTA accessions missing from `Accessions_metadata.csv` in `missing_accessions_from_metadata.txt`; both files live in the parent folder of `included_refid_fastas/`.
The per-RefID metadata split step writes CSVs only for RefIDs that have explicit filters under `refid_metadata/` and prints filter, kept row count, and total row count. Current filters: `30` source_isolate contains `Day1`; `85` accession is listed in `85.csv`; `142` source_isolate contains `baseline`; `192` source_isolate contains `day 1`; `346` source_isolate contains `baseline/D0`; `499` source_isolate contains `HCC`; `600` source_isolate does not contain `failure`; `661` source_isolation_source equals `plasma`; `884` source_isolate contains `Pre-TH`; `943` source_isolate contains `Day 1`; `1356` source_isolate does not contain `IC`; `2008` source_isolate does not contain `chimpanzee`; `2110` source_isolate contains `T0`; `2116` source_collection_date is before 2011; `2138` source_isolate contains `Week 0`; `2150` source_isolate contains `b`; `2168` source_isolate contains `pre`; `2178` source_isolation_source equals `plasma`. The manual accession list is a durable input in `HCVData/Ref-selection/NS5_Ref_filter/NS3/`.
The per-RefID FASTA filtering step reads `refid_metadata/RefID_<RefID>_metadata.csv`, keeps only matching `Accession` records in the corresponding copied FASTA file under `included_refid_fastas/`, and prints per-RefID and total before/after record counts.
The priority-assignment stage uses `HCVData/Reference_seqs/HCV_Subtype_Refs_AA_Accession_Subtype.csv` as the highest-priority accession/subtype source, then selects non-COMET calls for retained accessions called `1d` and for genotype 7 or 8 accessions. The following genotype and subtype steps consume that single selection to override COMET calls or add accessions absent from COMET. Amino-acid extraction reads the configured FASTA pool for added accessions.
The pre-profile alignment QC is an eligibility gate: only rows marked `PASS` in `NS3_Profile_Input_Alignment_QC.xlsx` are used for genotype/subtype profiles and their profile-accession list. The gate excludes coordinate-span failures, missing AA/coordinates, and sequences with at least 30% genotype-reference disagreement across at least 150 comparable AA positions. All excluded rows remain in the QC workbook and CSV, grouped by status and reason.

## Outputs

The workflow writes NS3 outputs under `outputs/comet-NS3/`, including:

- `NS3_GT_AllStudies.xlsx` (Comet genotype calls plus per-GT NS3 nucleotide distances and aligned-nucleotide counts)
- `NS3_matched_fasta_files.txt`
- discovery `filtered_rows.xlsx` under `outputs/comet-NS3/temp/.../find_refid_fastas/...`
- copied included RefID FASTA files under `03_stage-refid-fastas/`, the COMET-filtered copy under `04_prepare-comet-assignments/`, and the metadata-filtered copy and `kept_accessions.csv` under `08_filter-refid-fastas/`
- `NS3_NonComet_Priority_Assignments.csv` under `05_select-noncomet-priority-assignments/`
- `included_accessions_metadata.csv`
- `missing_accessions_from_metadata.txt`
- `refid_metadata/RefID_<RefID>_metadata.csv`
- filtered copied RefID FASTA files in `included_refid_fastas/` for RefIDs with metadata filters
- `NS3_Subtype_AllStudies_WSeqs.xlsx`
- `NS3_Subtype_With_GT_AA.xlsx`
- `NS3_Profile_Input_Alignment_QC.xlsx` (profile input with per-accession alignment QC columns and a `Flagged_Accessions` sheet)
- `NS3_QC_Passed_Genotype_Mutation_Burden_Summary.csv` (per-genotype mutation burden among QC-passed input rows)
- `NS3_Profile_Alignment_QC_Flagged_Accessions.csv` (flagged accessions and RAS positions requiring review)
- `NS3_Profile_Accessions_QC_Pass.csv` (accessions retained for complete-profile construction)
- `NS3_GT_NA_Distance_RAS.xlsx` and `NS3_Subtype_NA_Distance_RAS.xlsx` (full pairwise nucleotide-distance matrices over NS3 RAS codons only; the subtype workbook has one matrix per genotype and requires at least 10 sequences per subtype)
- `NS3_GT_AA_Distance_RAS.xlsx` and `NS3_Subtype_AA_Distance_RAS.xlsx` (equivalent full pairwise amino-acid distance matrices over NS3 RAS positions only)
- `NS3_GT_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_Subtype_CompleteProfiles_TabsPerGT.xlsx`
- `NS3_Subtype_CompleteProfiles_Merged.xlsx` (one merged subtype table with `Subtype`, `NS3Position`, `NumSeqsIncludingPosition`, `AminoAcid`, `CountWithAA`, and `PctWithAA`)
- `NS3_GT_Consensus.fasta`
- `NS3_Subtype_Consensus.fasta`
- `NS3_Subtype_Consensus_Aligned_to_GT1_1a.fasta` (all subtype consensuses aligned to GT1_1a coordinates)
- `NS3_GT_RAS_Profiles.xlsx`
- `NS3_Subtype_RAS_Profiles.xlsx`
- `NS3_Subtype_RAS_Profiles_Explicit_AA.xlsx` (all reportable subtype amino acids at RAS positions; used by the ICTV publication step)
- `NS3_Combined_RAS_Profiles.xlsx`
- `NS3_Subtype_RAS_Consensus_Difference_Summary.xlsx` (per-subtype mean and median RAS AA differences from genotype consensus)
- `NS3_GT5_5a_Profile_Coverage.xlsx`, `NS3_GT5_5a_Profile_Position_Coverage.csv`, and `NS3_GT5_5a_Profile_Position_Coverage.png` (accession-level, per-position, and charted subtype 5a coverage across full NS3 positions 1-631; the chart includes ambiguous and stop calls)
- `NS3_GT7_GT8_Step_Sequence_Audit.csv`, `NS3_GT7_GT8_Step_Sequence_Audit_Accessions.csv`, and `NS3_GT7_GT8_Step_Sequence_Audit_Summary.xlsx` (per-step GT7/GT8 kept counts, exclusions compared with the preceding step, accession-level evidence, and an Excel key-changes summary)

## Operating Rules

- Keep NS3 scripts together in this skill folder.
- Use `scripts/run_ns3_pipeline.py` for complete runs or `--step <name>` for a specific build step.
- Keep `.env` and `pipeline.local.toml` in the repository root; do not copy them into this skill folder.
- Keep temporary outputs under `outputs/comet-NS3/temp/` so they do not mix with other skills.
- Preserve the order above because later reports consume earlier workbooks.
