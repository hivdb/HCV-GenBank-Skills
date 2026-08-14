# Script Directory Map

| Folder | Purpose |
|---|---|
| `HCVData/HCV-all-seq-subtype` | Adds COMET subtype and coverage information to reference-accession tables. |
| `Preprocess/align_subtype_refs_to_first` | Pairwise-aligns subtype reference FASTAs to the first record. |
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

Only directories that contain scripts are listed. The `hcv-workflow` directory is intentionally excluded.
