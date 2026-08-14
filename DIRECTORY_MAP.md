# Top-Level Directory Map

| Folder | Purpose |
|---|---|
| `.git` | Git repository history and metadata. |
| `.venv` | Local Python environment. |
| `HCVData` | HCV input data: accession metadata, sequence and genotype/subtype references, COMET CSVs, reference-selection workbooks, BlastHist workbooks, and the all-sequence subtype/coverage dataset. |
| `Reference_seqs` | Genotype/subtype reference FASTAs, consensus comparisons, and mutation reports. |
| `input_data` | Input staging folder; currently essentially empty. |
| `outputs` | Generated pipeline reports, including local-alignment and COMET results. |
| `temp` | Temporary and intermediate generated files. |
| `notes` | Workflow and distance-analysis notes. |
| `scripts` | Shared scripts for consensus, audits, profile comparison, and reporting. |
| `genbank-accession-list-metadata` | Skill/workflow for turning GenBank accession lists into metadata and cohort summaries. |
| `genbank-gene-split-alignment` | Extracts GenBank sequences, aligns them, and splits results by gene. |
| `genbank-reference-alignment` | Aligns one or more accessions against reference FASTAs. |
| `genbank-single-accession-extractor` | Downloads or reads one accession and extracts sequence plus source metadata. |
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

The main active areas are `HCVData`, `Reference_seqs`, the six per-gene workflow folders, and `outputs`.
