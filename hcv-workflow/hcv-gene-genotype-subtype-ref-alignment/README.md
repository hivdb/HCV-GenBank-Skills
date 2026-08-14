# HCV Gene / Genotype / Subtype Reference Alignment

## Important: do not use this workflow to create production AA references

`build_hcv_gene_subtype_refs/build_hcv_gene_subtype_refs.py` locates per-gene sequences by
aligning subtype nucleotide genomes to genotype nucleotide references.  That
NA-alignment approach is not reliable enough for the production subtype
amino-acid references.

To rebuild `HCVData/Reference_seqs/HCV_Subtype_Refs_NS3_AA.fasta`,
`HCV_Subtype_Refs_NS5A_NTD_AA.fasta`, and
`HCV_Subtype_Refs_NS5B_AA.fasta`, use
`Preprocess/rebuild_subtype_reference_aas_from_genbank/rebuild_subtype_reference_aas_from_genbank.py` from the repository
root.  It uses each accession's GenBank protein annotation and then extracts
the gene AA sequence by alignment to `HCV-Ref-H77-Genotype1.fasta`.

Keep this folder's NA-alignment script only for legacy NA extraction,
diagnostics, and alignment-quality reporting.
