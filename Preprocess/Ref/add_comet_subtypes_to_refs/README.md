# add_comet_subtypes_to_refs

Adds COMET subtype and gene-coverage summaries to the reference-selection table.

The script reads the `HCVData/HCV-all-seq-subtype` dataset by default, with
non-COMET coverage files read from `HCVData/nonComet-Full-genome`:

- `Ref.csv`: rows to enrich; must contain `RefID`.
- `Accessions.csv`: maps accessions to RefIDs.
- `HCVData/Comet-Full-genome/all_comet_subtype.csv`: supplies COMET subtype calls.
- `HCVData/nonComet-Full-genome/NS3_AllSeq_NonComet_Coverage.csv`, `HCVData/nonComet-Full-genome/NS5A_AllSeq_NonComet_Coverage.csv`, and `HCVData/nonComet-Full-genome/NS5B_AllSeq_NonComet_Coverage.csv`: supply priority non-COMET subtypes and per-gene coverage. Only rows with `FullyCover=Yes` contribute to the coverage counts.

It writes `Ref_with_CometSubtypes.csv` in the same dataset folder. The output adds:

- `CometSubtypes`: all non-`unassigned` COMET calls for accessions belonging to the RefID, plus priority non-COMET subtypes `1d`, `7a`, `7b`, and `8a`.

It also updates `HCV_BlastHists_202604_data_Aug19.xlsx`. In each `Original_NS3`, `Original_NS5A`, and `Original_NS5B` sheet, non-empty `RefID` rows receive `CometSubtypes` and the gene-specific `Include...` coverage count. Blank subtype values are written as `(unassigned)`.
- `IncludeNS3Pos36_175` and `IncludeNS5APos26_93`: numbers of accessions with `FullyCover=Yes` in the corresponding coverage file. `Includes S282 + all positions` counts NS5B accessions whose `ReferenceOverlapAA` covers all seven RAS positions: 150, 159, 206, 282, 316, 320, and 321.
- In `Original_NS5B`, three additional counts are calculated from `ReferenceOverlapAA`: `Includes S282`, `Includes S282 + 4 other RAS positions`, and `Includes S282 + 5 other RAS positions`. The latter two require S282 plus at least four or five, respectively, of 150, 159, 206, 316, 320, and 321.

Run from any directory with:

```bash
.venv/bin/python Preprocess/Ref/add_comet_subtypes_to_refs/add_comet_subtypes_to_refs.py
```

Use the `--ref-csv`, `--accessions-csv`, `--comet-csv`, `--coverage-csv`, and `--output-csv` options to override the default dataset paths.
