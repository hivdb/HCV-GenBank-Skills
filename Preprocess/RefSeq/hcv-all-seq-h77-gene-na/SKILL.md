---
name: hcv-all-seq-h77-gene-na
description: Extract H77-matched NS3, NS5A, and NS5B nucleotide segments for every accession in an HCV FASTA. Use when gene FASTAs are needed from all.fasta-style sequences using the H77 genotype-1 gene references.
---

# Extract H77-matched HCV gene nucleotide sequences

Run from the repository root. The script BLASTs every sequence in
`HCVData/HCV-all-seq-subtype/all.fasta` against the NS3, NS5A, and NS5B
amino-acid references in `HCVData/HCV-Ref-H77-Genotype1.fasta` with BLASTX.
For the best match of each gene per accession, it writes the continuous matched
query nucleotide segment in H77 orientation.

```bash
.venv/bin/python Preprocess/RefSeq/hcv-all-seq-h77-gene-na/scripts/extract_h77_gene_na.py \
  --input-fasta HCVData/HCV-all-seq-subtype/all.fasta \
  --h77-reference HCVData/HCV-Ref-H77-Genotype1.fasta \
  --output-dir outputs/hcv-all-seq-h77-gene-na \
  --threads 4
```

It writes:

- `outputs/hcv-all-seq-h77-gene-na/HCV_AllSeq_H77_NS3_NA.fasta`
- `outputs/hcv-all-seq-h77-gene-na/HCV_AllSeq_H77_NS5A_NA.fasta`
- `outputs/hcv-all-seq-h77-gene-na/HCV_AllSeq_H77_NS5B_NA.fasta`

It then splits each gene FASTA into 1,000-record batches in its gene folder:

- `outputs/hcv-all-seq-h77-gene-na/NS3/HCV_AllSeq_H77_NS3_NA_part_0001.fasta`
- `outputs/hcv-all-seq-h77-gene-na/NS5A/HCV_AllSeq_H77_NS5A_NA_part_0001.fasta`
- `outputs/hcv-all-seq-h77-gene-na/NS5B/HCV_AllSeq_H77_NS5B_NA_part_0001.fasta`

Use `--records-per-file` to change the batch size. To recreate batches from
the existing gene FASTAs without rerunning BLASTX, use `--split-only`.

## Rename gene CSV files

To rename all CSV files directly in a gene folder, sort them by filename and
assign sequential part numbers, run:

```bash
.venv/bin/python Preprocess/RefSeq/hcv-all-seq-h77-gene-na/scripts/rename_gene_csv_parts.py \
  --directory outputs/hcv-all-seq-h77-gene-na/NS3 \
  --prefix HCV_AllSeq_H77_NS3_NA
```

This renames the files to `HCV_AllSeq_H77_NS3_NA_part_0001.csv`,
`HCV_AllSeq_H77_NS3_NA_part_0002.csv`, and so on. It does not rename FASTA
files in the same folder.

## Find matching CSV and FASTA files

To find the correct CSV–FASTA pairs in a gene folder, run:

```bash
.venv/bin/python Preprocess/RefSeq/hcv-all-seq-h77-gene-na/scripts/check_csv_fasta_accessions.py \
  --directory outputs/hcv-all-seq-h77-gene-na/NS3
```

The check extracts the accession before `|` from both CSV `name` or
`Accession` values and FASTA headers, then compares every CSV accession set to
every FASTA accession set. It prints each full `MATCH` regardless of filename.
It prints `UNMATCHED CSV` or `UNMATCHED FASTA` only for files with no complete
accession-set match, and exits unsuccessfully if any unmatched or ambiguous
files remain.

Add `--rename-csv-to-fasta` to rename each uniquely matched CSV to the
matching FASTA basename while preserving the `.csv` extension. The rename is
not performed if any CSV or FASTA remains unmatched or ambiguous.

Each output header contains the accession, query coordinates, H77 amino-acid
coordinates, extracted nucleotide length, aligned amino-acid length, and
percent identity. The default 100-nt minimum alignment length can be changed
with `--min-aligned-nt`.
