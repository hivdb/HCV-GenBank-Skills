# Statistic-analysis inputs

This directory contains copies of the per-accession amino-acid call exports
used for downstream statistical analysis. Each `*_Profile_Accession_AA_Calls.csv`
file has one usable amino-acid call per accession and selected RAS position,
with its genotype, subtype, genotype consensus call, and observed amino acid.

Each matching `*_Profile_Accession_AA_Calls_Accession_Totals.csv` file contains
the overall number of unique accessions in that export. It does not contain a
genotype or subtype breakdown.

## Source workbooks and workflow outputs

| Copied files | Source workflow output |
| --- | --- |
| `NS3_Profile_Accession_AA_Calls*.csv` | `outputs/comet-NS3-one-ras/report/` (generated in step 17, `merge-subtype-complete-profiles`) |
| `NS5A_Profile_Accession_AA_Calls*.csv` | `outputs/comet-NS5A-one-ras/report/` (generated in step 17, `merge-subtype-complete-profiles`) |
| `NS5B_Profile_Accession_AA_Calls*.csv` | `outputs/comet-NS5B-position-282-four-ras/report/` (generated in step 17, `merge-subtype-complete-profiles`) |

The source reports are copies of the corresponding step-17 outputs. The
exporter is `hcv-workflow/export_profile_accession_aa_calls.py`.
