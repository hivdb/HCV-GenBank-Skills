# Top-Level Directory Map

| Folder | Purpose |
|---|---|
| `.git` | Git repository history and metadata. |
| `.venv` | Local Python environment. |
| `HCVData` | HCV input data: accession metadata, sequence and genotype/subtype references in `Reference_seqs`, COMET CSVs, reference-selection workbooks, BlastHist workbooks, and the all-sequence subtype/coverage dataset. Generated reference comparison and mutation reports are in `outputs/reference_seqs`. |
| `outputs` | Generated pipeline reports, including local-alignment and COMET results, plus the `temp` subfolder for temporary and intermediate files. |
| `notes` | Workflow and distance-analysis notes. |
| `add_combined_profile_nonconsensus_row` | Adds normal-amino-acid non-consensus fractions to combined profiles. |
| `add_gt_counts_sheet` | Adds genotype count information to profile workbooks. |
| `add_subtype_consensus_mutation_summaries` | Adds subtype reference-to-consensus mutation summaries. |
| `align_subtype_refs_to_first` | Pairwise-aligns subtype reference FASTAs to the first record. |
| `build_comet_subtype_ras_coverage_report` | Builds COMET subtype RAS coverage reports. |
| `build_comet_workflow_sequence_audit` | Builds COMET workflow sequence and filtering audits. |
| `build_qc_passed_genotype_mutation_burden_summary` | Summarizes QC-passed genotype mutation burden. |
| `check_subtype_ref_edges` | Checks subtype reference edge cases and consistency. |
| `collect_genbank_by_fasta` | Collects GenBank records matching FASTA accessions. |
| `detect_accession_hcv_genes` | Detects NS3, NS5A, and NS5B gene coverage. |
| `export_check_refid_filter_summary` | Exports RefID `Check` filtering summaries. |
| `export_gt_reference_consensus_differences` | Exports reference-to-COMET-consensus comparisons. |
| `export_noncomet_priority_profile_accessions` | Exports non-COMET priority profile accessions. |
| `export_ras_consensus_mutations_csv` | Exports RAS consensus mutations as CSV. |
| `export_ref_selection_accessions_missing_from_comet` | Exports selected-reference accessions absent from COMET CSVs. |
| `extract_gt_refs_aa_to_fasta` | Extracts genotype amino-acid references to FASTA. |
| `extract_subtype_consensus_boundaries` | Extracts subtype consensus coordinate boundaries. |
| `rebuild_subtype_reference_aas_from_genbank` | Rebuilds subtype amino-acid references from GenBank proteins. |
| `replace_comet_profile_coverage_range_with_mean_diff` | Replaces combined-profile coverage ranges with non-X coverage. |
| `update_ictv_subtype_genome_json` | Updates the ICTV subtype-genome JSON. |
| `update_ictv_subtype_reference_fastas` | Updates ICTV subtype names and reference FASTAs. |
| `hcv-accessions-metadata-csv` | Creates RefID/accession metadata tables from FASTA and local GenBank data. |
| `hcv-all-seq-noncomet-coverage-audit` | Non-COMET genotype/subtype assignment and target-position coverage audit for `all.fasta`. |
| `hcv-combine-included-fastas` | Combines selected reference FASTAs into gene-level files. |
| `hcv-comet-subtype-profile-check` | Compares COMET subtype calls with profile-builder inclusion by gene. |
| `hcv-excel-refid-fasta-discovery` | Intended RefID/FASTA discovery workflow; currently appears empty or incomplete. |
| `hcv-folder-genotype-subtype-assignment` | Assigns genotype first, then closest subtype, to FASTAs in a folder. |
| `hcv-gene-genotype-subtype-ref-alignment` | Prepares and aligns genotype and subtype references for NS3, NS5A, and NS5B. |
| `hcv-metadata-subtype-consensus-workflow` | Builds subtype profiles and consensus sequences from metadata and RefID FASTAs. |
| `hcv-profile-input-comparison` | Compares profile inputs and local versus COMET subtype assignments. |
| `hcv-ns3-build-workflow` | Local NS3 reference/profile workflow. |
| `hcv-ns3-comet-build-workflow` | COMET-based NS3 reference/profile workflow. |
| `hcv-ns5a-build-workflow` | Local NS5A reference/profile workflow. |
| `hcv-ns5a-comet-build-workflow` | COMET-based NS5A reference/profile workflow. |
| `hcv-ns5b-build-workflow` | Local NS5B reference/profile workflow. |
| `hcv-ns5b-comet-build-workflow` | COMET-based NS5B reference/profile workflow. |

The main active areas are `HCVData`, the six per-gene workflow folders, and `outputs`.
