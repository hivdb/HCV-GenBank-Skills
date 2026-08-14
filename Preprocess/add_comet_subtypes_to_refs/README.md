# add_comet_subtypes_to_refs

Adds COMET subtype and gene-coverage summaries to the reference-selection table.

The script reads the `HCVData/HCV-all-seq-subtype` dataset by default:

- `Ref.csv`: rows to enrich; must contain `RefID`.
- `Accessions.csv`: maps accessions to RefIDs.
- `all_comet_subtype.csv`: supplies COMET subtype calls.
- `NS3_AllSeq_NonComet_Coverage.csv`, `NS5A_AllSeq_NonComet_Coverage.csv`, and `NS5B_AllSeq_NonComet_Coverage.csv`: supply priority non-COMET subtypes and per-gene coverage.

It writes `Ref_with_CometSubtypes.csv` in the same dataset folder. The output adds:

- `CometSubtypes`: all non-`unassigned` COMET calls for accessions belonging to the RefID, plus priority non-COMET subtypes `1d`, `7a`, `7b`, and `8a`.
- `IncludeNS3Pos36_175`, `IncludeNS5APos26_93`, and `IncludeNS5BPos150_321`: numbers of accessions with a nonblank `ReferenceOverlapAA` in the corresponding coverage file.

Run from any directory with:

```bash
.venv/bin/python Preprocess/add_comet_subtypes_to_refs/add_comet_subtypes_to_refs.py
```

Use the `--ref-csv`, `--accessions-csv`, `--comet-csv`, `--coverage-csv`, and `--output-csv` options to override the default dataset paths.
