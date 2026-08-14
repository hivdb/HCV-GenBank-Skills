# Script Directory Map

| Folder | Purpose |
|---|---|
| `HCVData` | Contains the all-sequence subtype/coverage table utility. |
| `Preprocess` | Contains subtype-reference alignment/validation and GenBank-record collection utilities. |
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

Only top-level directories that contain scripts are listed. The `hcv-workflow` directory is intentionally excluded.
