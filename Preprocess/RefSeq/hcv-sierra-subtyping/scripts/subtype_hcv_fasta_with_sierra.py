#!/usr/bin/env python3
"""Call an HCV-compatible Sierra GraphQL service for each FASTA sequence."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import re
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT_FASTA = REPO_ROOT / "HCVData" / "HCV-all-seq-subtype" / "all.fasta"
DEFAULT_OUTPUT_CSV = REPO_ROOT / "HCVData" / "sierra-subtyping" / "Sierra_HCV_Subtypes.csv"
SIERRA_BATCH_CAP = 120
DEFAULT_BATCH_SIZE = 100
DEFAULT_WORKERS = 8
HTTP_TIMEOUT = 180
HTTP_RETRIES = 2
DISPLAY_RE = re.compile(r"^(?P<sub>[^()]+?)\s*\(.+\)$")
QUERY = (
    "query Q($sequences:[UnalignedSequenceInput]!){"
    "sequenceAnalysis(sequences:$sequences){"
    "inputSequence{header} bestMatchingSubtype{display}}}"
)


def parse_subtype(display: str | None) -> str:
    """Return Sierra's subtype label without its trailing distance."""
    if not display:
        return ""
    match = DISPLAY_RE.match(display)
    subtype = match.group("sub").strip() if match else display.strip()
    return "" if subtype.lower() == "unknown" else subtype


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    header = ""
    sequence: list[str] = []
    seen: set[str] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(sequence).upper()))
            header = line[1:].split()[0].split(".", 1)[0]
            if not header or header in seen:
                raise ValueError(f"{path} has a missing or duplicate accession: {header!r}")
            seen.add(header)
            sequence = []
        elif line.strip():
            sequence.append(re.sub(r"\s+", "", line))
    if header:
        records.append((header, "".join(sequence).upper()))
    return records


def read_existing(path: Path) -> dict[str, str]:
    if not path.is_file():
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Accession", "SierraSubtype"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        return {
            str(row["Accession"]).strip(): str(row.get("SierraSubtype") or "").strip()
            for row in reader
            if str(row.get("Accession") or "").strip()
        }


def call_sierra(url: str, batch: list[tuple[str, str]]) -> dict[str, str]:
    sequences = [{"header": accession, "sequence": sequence} for accession, sequence in batch]
    payload = json.dumps({"variables": {"sequences": sequences}, "query": QUERY}).encode()
    last_error: Exception | None = None
    for attempt in range(HTTP_RETRIES + 1):
        request = urllib.request.Request(
            url, data=payload, headers={"Content-Type": "application/json"}
        )
        try:
            response = urllib.request.urlopen(request, timeout=HTTP_TIMEOUT)
            document = json.loads(response.read())
            break
        except urllib.error.HTTPError as error:
            body = error.read().decode("utf-8", "replace")[:400]
            if 400 <= error.code < 500:
                raise RuntimeError(f"{url} -> HTTP {error.code}: {body}") from error
            last_error = RuntimeError(f"{url} -> HTTP {error.code}: {body}")
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            last_error = error
        time.sleep(1.5 * (attempt + 1))
    else:
        raise last_error or RuntimeError(f"{url} is unreachable")

    if "errors" in document:
        raise RuntimeError(f"{url} -> GraphQL errors: {document['errors']}")
    calls: dict[str, str] = {}
    for entry in document["data"]["sequenceAnalysis"]:
        accession = str(entry["inputSequence"]["header"])
        match = entry.get("bestMatchingSubtype")
        calls[accession] = parse_subtype(match.get("display") if match else None)
    return calls


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-fasta", type=Path, default=DEFAULT_INPUT_FASTA)
    parser.add_argument("--output-csv", type=Path, default=DEFAULT_OUTPUT_CSV)
    parser.add_argument(
        "--endpoints",
        required=True,
        help="comma-separated HCV-compatible Sierra GraphQL endpoint URLs",
    )
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    parser.add_argument("--limit", type=int, help="maximum number of pending sequences")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if not 1 <= args.batch_size <= SIERRA_BATCH_CAP:
        raise SystemExit(f"--batch-size must be between 1 and {SIERRA_BATCH_CAP}")
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    endpoints = [value.strip() for value in args.endpoints.split(",") if value.strip()]
    if not endpoints:
        raise SystemExit("--endpoints must contain at least one URL")

    existing = read_existing(args.output_csv)
    records = read_fasta(args.input_fasta)
    pending = [(accession, sequence) for accession, sequence in records if accession not in existing]
    if args.limit is not None:
        pending = pending[: args.limit]
    print(f"input={len(records):,} existing={len(existing):,} pending={len(pending):,}")
    if not pending:
        return 0

    batches = [pending[index : index + args.batch_size] for index in range(0, len(pending), args.batch_size)]
    work = list(zip(itertools.cycle(endpoints), batches))
    failures = 0
    start = time.monotonic()
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(call_sierra, endpoint, batch): batch for endpoint, batch in work
        }
        for completed, future in enumerate(as_completed(futures), start=1):
            batch = futures[future]
            try:
                calls = future.result()
            except Exception as error:
                failures += 1
                print(f"! batch failed ({batch[0][0]}): {error}", file=sys.stderr)
                continue
            for accession, _ in batch:
                existing[accession] = calls.get(accession, "")
            if completed % 10 == 0 or completed == len(futures):
                print(f"completed_batches={completed}/{len(futures)}")

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["Accession", "SierraSubtype"], lineterminator="\n"
        )
        writer.writeheader()
        for accession, _ in records:
            if accession in existing:
                writer.writerow({"Accession": accession, "SierraSubtype": existing[accession]})
    elapsed = time.monotonic() - start
    print(f"saved={args.output_csv} failures={failures} elapsed_seconds={elapsed:.1f}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
