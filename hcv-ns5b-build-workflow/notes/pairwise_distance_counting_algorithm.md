# Pairwise NA and AA distance calculation

The distance workbooks report the **mean pairwise mismatch proportion** across the selected positions. They do not report a median and do not compare consensus sequences.

## Inputs included in a comparison

Only accessions in the final profile-accession CSV are considered. A sequence must cover every selected position and have an unambiguous call at each one:

- amino-acid (AA) comparisons accept the 20 standard amino acids;
- nucleotide (NA) comparisons accept `A`, `C`, `G`, and `T`.

For NA comparisons, each selected amino-acid RAS position contributes its three codon bases. For AA comparisons, each selected RAS position contributes one amino-acid call.

## Direct pair-by-pair definition

For two groups, compare every sequence in group A with every sequence in group B at every selected position.

```text
mean distance = total mismatching calls / total compared calls
```

The diagonal of a distance matrix compares every unique pair within the same group. For a group of `n` sequences, that is `n × (n - 1) / 2` pairs.

## Faster count-based calculation

The implementation gives the identical mean without explicitly iterating through all sequence pairs.

At each selected position, it counts each observed call. For example:

```text
Group A: A = 80, G = 20
Group B: A = 70, G = 30
```

There are `80 × 30` individual `A`-versus-`G` comparisons and `20 × 70` individual `G`-versus-`A` comparisons. Both categories are mismatches. The matching categories are `80 × 70` (`A` versus `A`) and `20 × 30` (`G` versus `G`).

Equivalently:

```text
all cross-group pairs = size(A) × size(B)
matching pairs         = sum(call_count_A × call_count_B) for each call
mismatching pairs      = all cross-group pairs - matching pairs
```

For within-group (diagonal) cells, the same calculation uses combinations rather than cross-products:

```text
all pairs      = n × (n - 1) / 2
matching pairs = sum(call_count × (call_count - 1) / 2) for each call
mismatching    = all pairs - matching pairs
```

The mismatching and comparable pair counts are added over all selected positions, then divided once at the end.

## Why the two methods are equal

Every individual sequence comparison at one position belongs to exactly one call-pair category. For example, every `A` call in group A compared with every `G` call in group B is one of the `count(A in A) × count(G in B)` comparisons. Count products therefore count every direct comparison exactly once; they only group equivalent comparisons together.

The approach preserves the **overall mean** across all positions and pairs. It does not preserve the distance of each named sequence pair, which is not needed for these group-level matrix cells.

### Two-position example

| Sequence | RAS 1 | RAS 2 |
| --- | --- | --- |
| A1 | A | T |
| A2 | G | T |
| B1 | A | C |
| B2 | G | C |

Direct comparisons:

| Pair | Differences | Compared positions | Pair distance |
| --- | ---: | ---: | ---: |
| A1 vs B1 | 1 | 2 | 50% |
| A1 vs B2 | 2 | 2 | 100% |
| A2 vs B1 | 2 | 2 | 100% |
| A2 vs B2 | 1 | 2 | 50% |

The direct mean is `(1 + 2 + 2 + 1) / (4 pairs × 2 positions) = 6 / 8 = 75%`.

Count-based calculation:

- At RAS 1, two of four cross-group comparisons differ.
- At RAS 2, all four cross-group comparisons differ.
- Total: six differences out of eight comparisons, or `6 / 8 = 75%`.

The result is therefore the same as direct pair-by-pair comparison.
