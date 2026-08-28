---
name: hcv-sierra-subtyping
description: Submit HCV nucleotide FASTA records to an HCV-compatible Sierra GraphQL service in concurrent batches and save resumable subtype calls. Use when an HCV Sierra deployment and its endpoint URLs are available.
---

# HCV Sierra subtyping

This workflow starts with `HCVData/HCV-all-seq-subtype/all.fasta` and writes
`HCVData/sierra-subtyping/Sierra_HCV_Subtypes.csv`.

It adapts the batching, concurrent endpoint use, HTTP retry, subtype-display
parsing, and incremental rerun behavior from
`HIVDB3-blasthit/backend/scripts/populate_sierra_subtype.py`. It does not use
that script's MySQL database logic.

An **HCV-compatible** Sierra GraphQL endpoint is required. Do not use an
HIV-only Sierra deployment: it cannot provide HCV subtype calls.

## Local Sierra container

```bash
docker pull hivdb/sierra:latest
docker run -d --name sierra -p 8111:8080 hivdb/sierra:latest dev
```

The local GraphQL endpoint is:

```text
http://localhost:8111/sierra/rest/graphql
```

View logs or stop the container with:

```bash
docker logs -f sierra
docker stop sierra
```

Run from the repository root:

```bash
.venv/bin/python Preprocess/RefSeq/hcv-sierra-subtyping/scripts/subtype_hcv_fasta_with_sierra.py \
  --endpoints https://your-hcv-sierra.example/rest/graphql
```

For a small connectivity test, add `--limit 100`. The service accepts at most
120 sequences per request; the default batch size is 100. A rerun sends only
accessions absent from the existing output CSV. Sierra `Unknown` calls are
saved as blank `SierraSubtype` values and are not retried automatically.
