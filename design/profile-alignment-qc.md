# Profile-alignment QC

All NS3, NS5A, and NS5B COMET build workflows apply this gate after amino-acid extraction and before complete-profile construction. Only records with `AlignmentQCStatus = PASS` contribute to genotype, subtype, and RAS profiles.

## Reference comparison

Each extracted amino-acid sequence is compared position-by-position with the genotype reference selected by its `ClosestGT` assignment in `HCV_GT_Refs_By_Gene_AA.json`. The comparison is gene-specific (`NS3`, `NS5A_NTD`, or `NS5B`); it does not use a subtype reference or subtype consensus.

Only standard amino acids in both the accession and reference are comparable. `X`, `*`, and non-standard amino acids do not count toward the comparable-position or mutation counts.

## QC order

Records are evaluated in this order:

1. `NO_AA_SEQUENCE` when `AASequence` is empty.
2. `MISSING_AA_COORDINATES` when `StartAAPosition` or `EndAAPosition` is missing.
3. `FLAGGED` with `coordinate_span_length_mismatch` when `EndAAPosition - StartAAPosition + 1` differs from the extracted sequence length.
4. `HIGH_DIVERGENCE` when the coordinate check passes, the minimum number of comparable amino acids is met, and genotype-reference disagreement is at least the gene threshold.
5. `PASS` otherwise.

The shared high-divergence threshold is **30%** for all three genes:

| Gene | Minimum comparable amino acids | Divergence threshold |
| --- | ---: | ---: |
| NS3 | 150 | ≥30% |
| NS5A | 150 | ≥30% |
| NS5B | 250 | ≥30% |

Non-pass records stay in the profile-input QC workbook and the flagged-accessions CSV, but are excluded from profile construction.

## NS5B post-QC eligibility

After a record passes alignment QC, NS5B workflows apply their own profile-coverage rule:

- The all-RAS workflow requires callable standard amino acids at positions 150, 159, 206, 282, 316, 320, and 321.
- The position-282-plus-four-RAS workflow requires a callable amino acid at 282 and at least four callable calls among the other six positions.
- The separate 90%-range report applies the corresponding RAS rule plus at least 90% non-`X` coverage across positions 150–321.
