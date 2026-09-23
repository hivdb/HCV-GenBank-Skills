#!/usr/bin/env python3
"""Stratified amino-acid covariation analysis for the wide HCV profiles.

For each RAS--RAS and RAS--non-RAS position pair, this program calculates
within-subtype mutual information (MI), combines it with both sequence and
subtype-balanced weights, and obtains permutation P values by shuffling one
position within subtype.  P values are Benjamini--Hochberg adjusted separately
for each scan family and weighting scheme.

The default command analyzes every ``*_Wide_*.csv`` file beside this script.
It creates results for the complete data and for an exact-profile-deduplicated
sensitivity data set.  Incomplete profiles are deliberately retained by the
exact-deduplication filter because their true complete profile is unknown.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import time
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, NamedTuple


RAS_POSITIONS = {
    "NS3": (36, 41, 43, 54, 55, 56, 80, 122, 155, 156, 158, 166, 168, 170, 175),
    "NS5A": (24, 26, 28, 29, 30, 31, 32, 38, 58, 62, 92, 93),
    "NS5B": (150, 159, 206, 282, 316, 320, 321),
}


class SubtypePair(NamedTuple):
    subtype: str
    x: tuple[str, ...]
    y: tuple[str, ...]

    @property
    def n(self) -> int:
        return len(self.x)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "inputs",
        nargs="*",
        type=Path,
        help="Wide-profile CSV files. Defaults to every *_Wide_*.csv in this directory.",
    )
    parser.add_argument("--permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=20260923)
    parser.add_argument("--min-n", type=int, default=20,
                        help="Minimum paired observations in a subtype (default: 20).")
    parser.add_argument("--min-minor", type=int, default=3,
                        help="Minimum observations outside the modal AA at each position (default: 3).")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Q-value cutoff for amino-acid-combination details (default: 0.05).")
    parser.add_argument("--min-combination-count", type=int, default=3)
    parser.add_argument("--include-all-position-scan", action="store_true",
                        help="Also analyze every pair of positions (a separate testing family).")
    parser.add_argument(
        "--near-duplicate-differences", type=int,
        metavar="K",
        help="Optional complete-linkage sensitivity analysis: retain one complete profile per subtype cluster with at most K differences. This can be slow.",
    )
    return parser.parse_args()


def gene_for(path: Path) -> str:
    name = path.name.upper()
    matches = [gene for gene in RAS_POSITIONS if gene in name]
    if len(matches) != 1:
        raise ValueError(f"Cannot identify one of {sorted(RAS_POSITIONS)} from {path.name}")
    return matches[0]


def read_profile(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"{path} has no header")
        positions = [name for name in reader.fieldnames if name and name.strip().isdigit()]
        required = {"Genotype", "Subtype"}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"{path} lacks required columns: {', '.join(sorted(missing))}")
        if not positions:
            raise ValueError(f"{path} has no numeric position columns")
        rows: list[dict[str, str]] = []
        for line_number, row in enumerate(reader, start=2):
            cleaned = {key: (value or "").strip() for key, value in row.items() if key is not None}
            if not cleaned["Subtype"]:
                raise ValueError(f"{path}:{line_number} has an empty Subtype")
            rows.append(cleaned)
    return rows, sorted(positions, key=int)


def exact_deduplicate(rows: list[dict[str, str]], positions: list[str]) -> list[dict[str, str]]:
    """Keep one fully observed profile per subtype; keep all incomplete rows."""
    seen: set[tuple[str, tuple[str, ...]]] = set()
    retained: list[dict[str, str]] = []
    for row in rows:
        profile = tuple(row[position] for position in positions)
        if not all(profile):
            retained.append(row)
            continue
        key = (row["Subtype"], profile)
        if key not in seen:
            seen.add(key)
            retained.append(row)
    return retained


def profile_distance(left: tuple[str, ...], right: tuple[str, ...], maximum: int) -> int:
    differences = 0
    for a, b in zip(left, right):
        if a != b:
            differences += 1
            if differences > maximum:
                return differences
    return differences


def complete_linkage_deduplicate(
    rows: list[dict[str, str]], positions: list[str], maximum_differences: int
) -> list[dict[str, str]]:
    """Retain one row per complete-linkage cluster, separately within subtype.

    This intentionally straightforward implementation checks every relevant
    profile pair. It is optional because large subtypes can make it expensive.
    """
    if maximum_differences < 0:
        raise ValueError("--near-duplicate-differences must be non-negative")
    complete: dict[str, list[dict[str, str]]] = defaultdict(list)
    incomplete: list[dict[str, str]] = []
    for row in rows:
        if all(row[position] for position in positions):
            complete[row["Subtype"]].append(row)
        else:
            incomplete.append(row)

    retained = list(incomplete)
    for subtype in sorted(complete):
        candidates = complete[subtype]
        # Agglomerative complete linkage: clusters merge only when every
        # cross-cluster profile pair is within the prespecified threshold.
        clusters: list[list[int]] = [[index] for index in range(len(candidates))]
        profiles = [tuple(row[position] for position in positions) for row in candidates]
        while True:
            merge: tuple[int, int] | None = None
            for left in range(len(clusters)):
                for right in range(left + 1, len(clusters)):
                    if all(
                        profile_distance(profiles[i], profiles[j], maximum_differences)
                        <= maximum_differences
                        for i in clusters[left] for j in clusters[right]
                    ):
                        merge = (left, right)
                        break
                if merge:
                    break
            if merge is None:
                break
            left, right = merge
            clusters[left].extend(clusters[right])
            del clusters[right]
        retained.extend(candidates[cluster[0]] for cluster in clusters)
    return retained


def mutual_information(x: Iterable[str], y: Iterable[str]) -> float:
    pairs = Counter(zip(x, y))
    x_counts = Counter()
    y_counts = Counter()
    for (x_value, y_value), count in pairs.items():
        x_counts[x_value] += count
        y_counts[y_value] += count
    if len(x_counts) < 2 or len(y_counts) < 2:
        return 0.0
    total = sum(pairs.values())
    return sum(
        (count / total) * math.log2((count * total) / (x_counts[x_value] * y_counts[y_value]))
        for (x_value, y_value), count in pairs.items()
    )


def eligible_subtypes(
    rows: list[dict[str, str]], pos1: str, pos2: str, min_n: int, min_minor: int
) -> list[SubtypePair]:
    grouped: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for row in rows:
        x, y = row[pos1], row[pos2]
        if x and y:
            grouped[row["Subtype"]].append((x, y))
    eligible: list[SubtypePair] = []
    for subtype in sorted(grouped):
        pairs = grouped[subtype]
        if len(pairs) < min_n:
            continue
        x, y = zip(*pairs)
        if len(x) - max(Counter(x).values()) < min_minor:
            continue
        if len(y) - max(Counter(y).values()) < min_minor:
            continue
        eligible.append(SubtypePair(subtype, x, y))
    return eligible


def weighted_mi(values: list[float], sizes: list[int], weighting: str) -> float:
    if weighting == "sequence_weighted":
        return sum(value * size for value, size in zip(values, sizes)) / sum(sizes)
    if weighting == "subtype_balanced":
        return sum(values) / len(values)
    raise ValueError(f"Unknown weighting: {weighting}")


def analyze_pair(
    rows: list[dict[str, str]], pos1: str, pos2: str, permutations: int,
    min_n: int, min_minor: int, rng: random.Random,
) -> list[dict[str, object]]:
    subtypes = eligible_subtypes(rows, pos1, pos2, min_n, min_minor)
    base = {
        "pos1": int(pos1), "pos2": int(pos2), "n_subtypes": len(subtypes),
        "n_sequences": sum(item.n for item in subtypes),
    }
    if not subtypes:
        return [{**base, "weighting": weighting, "observed_mi_bits": "", "mean_permuted_mi_bits": "",
                 "adjusted_mi_bits": "", "p_value": ""}
                for weighting in ("sequence_weighted", "subtype_balanced")]

    sizes = [item.n for item in subtypes]
    observed_values = [mutual_information(item.x, item.y) for item in subtypes]
    observed = {
        weighting: weighted_mi(observed_values, sizes, weighting)
        for weighting in ("sequence_weighted", "subtype_balanced")
    }
    null_values = {weighting: [] for weighting in observed}
    for _ in range(permutations):
        permuted_mi = []
        for item in subtypes:
            shuffled = list(item.y)
            rng.shuffle(shuffled)
            permuted_mi.append(mutual_information(item.x, shuffled))
        for weighting in observed:
            null_values[weighting].append(weighted_mi(permuted_mi, sizes, weighting))

    result = []
    for weighting, observed_mi in observed.items():
        null = null_values[weighting]
        mean_null = sum(null) / permutations
        result.append({
            **base,
            "weighting": weighting,
            "observed_mi_bits": observed_mi,
            "mean_permuted_mi_bits": mean_null,
            "adjusted_mi_bits": observed_mi - mean_null,
            "p_value": (1 + sum(value >= observed_mi for value in null)) / (permutations + 1),
        })
    return result


def benjamini_hochberg(rows: list[dict[str, object]]) -> None:
    """Add q_value, separately for each weighting scheme, ignoring ineligible pairs."""
    for weighting in ("sequence_weighted", "subtype_balanced"):
        valid = [row for row in rows if row["weighting"] == weighting and row["p_value"] != ""]
        valid.sort(key=lambda row: float(row["p_value"]))
        total = len(valid)
        running = 1.0
        for rank in range(total, 0, -1):
            row = valid[rank - 1]
            running = min(running, float(row["p_value"]) * total / rank)
            row["q_value"] = running
        for row in rows:
            if row["weighting"] == weighting and "q_value" not in row:
                row["q_value"] = ""


def subtype_mi_rows(
    rows: list[dict[str, str]], selected_pairs: set[tuple[int, int]], min_n: int, min_minor: int
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for pos1, pos2 in sorted(selected_pairs):
        for item in eligible_subtypes(rows, str(pos1), str(pos2), min_n, min_minor):
            output.append({"pos1": pos1, "pos2": pos2, "subtype": item.subtype,
                           "n_sequences": item.n,
                           "mi_bits": mutual_information(item.x, item.y)})
    return output


def combination_rows(
    rows: list[dict[str, str]], selected_pairs: set[tuple[int, int]], min_n: int,
    min_minor: int, minimum_count: int,
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for pos1, pos2 in sorted(selected_pairs):
        for item in eligible_subtypes(rows, str(pos1), str(pos2), min_n, min_minor):
            pairs = Counter(zip(item.x, item.y))
            x_counts, y_counts = Counter(item.x), Counter(item.y)
            for (aa1, aa2), observed in sorted(pairs.items()):
                if observed < minimum_count:
                    continue
                expected = x_counts[aa1] * y_counts[aa2] / item.n
                output.append({
                    "pos1": pos1, "aa1": aa1, "pos2": pos2, "aa2": aa2,
                    "subtype": item.subtype, "n_sequences": item.n,
                    "observed_count": observed, "expected_count_independent": expected,
                    "lift_observed_over_expected": observed / expected if expected else "",
                    "log2_lift": math.log2(observed / expected) if expected else "",
                })
    return output


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            formatted = {
                key: (f"{value:.10g}" if isinstance(value, float) else value)
                for key, value in row.items()
            }
            writer.writerow(formatted)


def scan(
    rows: list[dict[str, str]], pairs: list[tuple[int, int]], args: argparse.Namespace,
    rng: random.Random,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    started = time.monotonic()
    for index, (pos1, pos2) in enumerate(pairs, start=1):
        result.extend(analyze_pair(rows, str(pos1), str(pos2), args.permutations,
                                   args.min_n, args.min_minor, rng))
        if index % 100 == 0 or index == len(pairs):
            elapsed_seconds = time.monotonic() - started
            percent = 100 * index / len(pairs)
            print(
                f"    progress: {index}/{len(pairs)} pairs ({percent:.1f}%), "
                f"elapsed {elapsed_seconds / 60:.1f} min",
                flush=True,
            )
    benjamini_hochberg(result)
    return sorted(result, key=lambda row: (
        row["weighting"], float(row["q_value"]) if row["q_value"] != "" else float("inf"),
        -float(row["adjusted_mi_bits"]) if row["adjusted_mi_bits"] != "" else 0,
    ))


def write_scan_outputs(
    output_prefix: Path, scan_name: str, rows: list[dict[str, object]], data_rows: list[dict[str, str]],
    args: argparse.Namespace,
) -> tuple[Path, int, int]:
    result_path = output_prefix.with_name(f"{output_prefix.name}_{scan_name}.csv")
    fields = ["pos1", "pos2", "weighting", "n_subtypes", "n_sequences", "observed_mi_bits",
              "mean_permuted_mi_bits", "adjusted_mi_bits", "p_value", "q_value"]
    write_csv(result_path, rows, fields)
    selected = {
        (int(row["pos1"]), int(row["pos2"])) for row in rows
        if row["q_value"] != "" and float(row["q_value"]) <= args.alpha
    }
    details_prefix = output_prefix.with_name(f"{output_prefix.name}_{scan_name}_significant")
    write_csv(
        details_prefix.with_name(f"{details_prefix.name}_subtype_mi.csv"),
        subtype_mi_rows(data_rows, selected, args.min_n, args.min_minor),
        ["pos1", "pos2", "subtype", "n_sequences", "mi_bits"],
    )
    write_csv(
        details_prefix.with_name(f"{details_prefix.name}_aa_combinations.csv"),
        combination_rows(data_rows, selected, args.min_n, args.min_minor, args.min_combination_count),
        ["pos1", "aa1", "pos2", "aa2", "subtype", "n_sequences", "observed_count",
         "expected_count_independent", "lift_observed_over_expected", "log2_lift"],
    )
    eligible_pairs = len({
        (int(row["pos1"]), int(row["pos2"])) for row in rows
        if row["p_value"] != ""
    })
    return result_path, eligible_pairs, len(selected)


def analyze_input(path: Path, args: argparse.Namespace, seed_offset: int) -> None:
    gene = gene_for(path)
    rows, position_names = read_profile(path)
    positions = [int(position) for position in position_names]
    ras = RAS_POSITIONS[gene]
    absent = sorted(set(ras) - set(positions))
    if absent:
        raise ValueError(f"{path} is missing {gene} RAS positions: {absent}")
    non_ras = [position for position in positions if position not in ras]
    families = {
        "ras_ras": list(itertools.combinations(ras, 2)),
        "ras_nonras": [(ras_position, other) for ras_position in ras for other in non_ras],
    }
    if args.include_all_position_scan:
        families["all_positions"] = list(itertools.combinations(positions, 2))

    data_sets: list[tuple[str, list[dict[str, str]]]] = [
        ("all_sequences", rows),
        ("exact_deduplicated", exact_deduplicate(rows, position_names)),
    ]
    if args.near_duplicate_differences is not None:
        data_sets.append((
            f"complete_linkage_{args.near_duplicate_differences}diff",
            complete_linkage_deduplicate(rows, position_names, args.near_duplicate_differences),
        ))

    complete_profiles = sum(all(row[position] for position in position_names) for row in rows)
    print(
        f"\nStarting {path.name}\n"
        f"  gene: {gene}; rows: {len(rows)}; subtypes: {len({row['Subtype'] for row in rows})}; "
        f"complete profiles: {complete_profiles}\n"
        f"  RAS positions: {','.join(map(str, ras))}",
        flush=True,
    )
    summary: list[dict[str, object]] = []
    for data_index, (label, data_rows) in enumerate(data_sets):
        print(f"  data set: {label} ({len(data_rows)} rows)", flush=True)
        summary.append({"dataset": label, "n_rows": len(data_rows), "n_subtypes": len({r['Subtype'] for r in data_rows})})
        prefix = path.with_name(f"{path.stem}_covariation_{label}")
        for family_index, (scan_name, pairs) in enumerate(families.items()):
            print(f"  scan: {scan_name} ({len(pairs)} pairs)", flush=True)
            rng = random.Random(args.seed + seed_offset + data_index * 10_000 + family_index)
            results = scan(data_rows, pairs, args, rng)
            result_path, eligible_pairs, significant_pairs = write_scan_outputs(
                prefix, scan_name, results, data_rows, args
            )
            print(
                f"    complete: {eligible_pairs} eligible pairs; "
                f"{significant_pairs} pairs with q <= {args.alpha:g}; wrote {result_path.name}",
                flush=True,
            )
    write_csv(path.with_name(f"{path.stem}_covariation_dataset_summary.csv"), summary,
              ["dataset", "n_rows", "n_subtypes"])


def main() -> None:
    args = parse_args()
    started = time.monotonic()
    if args.permutations < 1:
        raise SystemExit("--permutations must be at least 1")
    if args.min_n < 1 or args.min_minor < 1:
        raise SystemExit("--min-n and --min-minor must be positive")
    script_dir = Path(__file__).resolve().parent
    # Result filenames retain the input stem, so exclude them explicitly on a
    # rerun rather than accidentally treating a previous result as input.
    inputs = args.inputs or sorted(
        path for path in script_dir.glob("*_Wide_*.csv")
        if "_covariation_" not in path.stem
    )
    if not inputs:
        raise SystemExit("No input CSV files found")
    print(
        f"Covariation analysis: {len(inputs)} input file(s); {args.permutations} permutations; "
        f"min_n={args.min_n}; min_minor={args.min_minor}; seed={args.seed}",
        flush=True,
    )
    for index, path in enumerate(inputs):
        analyze_input(path.resolve(), args, index * 100_000)
    elapsed_seconds = time.monotonic() - started
    hours, remainder = divmod(int(round(elapsed_seconds)), 3600)
    minutes, seconds = divmod(remainder, 60)
    print(
        f"\nCovariation analysis complete: {len(inputs)} input file(s) finished in "
        f"{hours:d}h {minutes:02d}m {seconds:02d}s.",
        flush=True,
    )


if __name__ == "__main__":
    main()
