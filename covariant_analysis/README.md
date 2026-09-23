# Covariation analysis

This directory contains three wide amino-acid profile inputs and
`analyze_covariation.py`, a dependency-free Python program for identifying
within-subtype amino-acid associations. It analyzes NS3, NS5A, and NS5B
separately, rather than combining genes.

## Basic analysis logic

Imagine each virus sequence is a row of letter boxes. Each box is a position
in a viral protein, and its letter is an amino acid. We ask a simple question:
**when one box has a certain letter, does another box tend to have a certain
letter too?** If so, the two positions may be covarying.

### A small example

Suppose we look only at 100 viruses from subtype 1a and compare positions 30
and 93:

| Position 30 | Position 93 | Number of viruses |
|---|---|---:|
| R | H | 15 |
| R | Y | 5 |
| Q | H | 5 |
| Q | Y | 75 |

There are 20 viruses with `R` at position 30 and 20 with `H` at position 93.
If the two positions had nothing to do with each other, we would expect about
4 viruses to have both `R` and `H` (`20% × 20% × 100`). We see 15 instead.
That suggests `30R` and `93H` travel together more often than chance would
predict.

Mutual information (MI) turns all four counts in the table into one number.
It gets larger when knowing the letter at one position helps us guess the
letter at the other position. MI does not require choosing one letter as
“normal,” and it can use more than two possible letters.

### How the fairness check works

Different subtypes naturally have different amino-acid letters. To avoid
mistaking these subtype differences for covariation, the program compares
positions only within the same subtype.

Then it performs a fairness check called a **shuffle**. It keeps the same 20
`H` letters and 80 `Y` letters at position 93, but randomly gives them to the
100 subtype-1a viruses again. Position 30 is not changed. This keeps the same
number of each letter but breaks the original position-30/position-93 pairs.

The program calculates MI after this shuffle many times:

- **199 shuffles** are useful for a quick practice run.
- **999 shuffles** are used for final results because they give a more precise
  answer.

If the real MI is larger than nearly all shuffled MI values, the two positions
are covarying more than expected by chance. The P value says how often a
shuffle produced MI at least as large as the real MI. The q value adjusts this
answer because the program tests many position pairs.

### How results from subtypes are combined

The program gives two summaries:

- **Sequence-weighted:** subtypes with more sequences, such as 1a, 1b, and
  3a, contribute more to the answer.
- **Subtype-balanced:** every adequately sampled subtype gets one equal vote.

Looking at both summaries helps us see whether a result is broad across
subtypes or mostly comes from one large subtype. For important pairs, the
program also lists the actual amino-acid combinations that caused the result.

Finally, the analysis is repeated after removing identical profiles within a
subtype. This checks whether repeated or very similar viruses are making a
result look stronger than it really is. The analysis can show that two amino
acids are associated, but it cannot prove that they directly interact; they
may simply have been inherited together in a virus family.

## Step-by-step workflow

1. Check the input files.

   The default run uses every original `*_Wide_*.csv` file in this directory.
   Each input must contain `Genotype`, `Subtype`, and numeric amino-acid
   position columns. Blank amino-acid calls are treated as missing. A subtype
   must be present for every row.

2. Confirm the complete position range used for each gene.

   - NS3: positions 36--175
   - NS5A: positions 24--93
   - NS5B: positions 150--321

   The program uses every position in these ranges to identify exact or
   near-duplicate profiles and for the RAS--non-RAS and optional all-position
   scans. The focused RAS--RAS scan uses the gene-specific RAS subset within
   its complete range; the script stops with an error if an expected RAS column
   is absent.

3. Run an exploratory analysis first.

   Use 199 permutations to assess run time and inspect provisional results.
   These P values are exploratory and should not be used as final results.

   ```sh
   python3 analyze_covariation.py --permutations 199
   ```

   While running, the program prints the current input file, gene, number of
   rows and subtypes, current data set (`all_sequences` or
   `exact_deduplicated`), scan family, percent of pairs completed, elapsed
   time, number of eligible pairs, number of pairs with `q <= 0.05`, and the
   result filename it wrote. After every input file is complete, it also
   prints the total elapsed run time.

### Command options

Run `python3 analyze_covariation.py --help` to see these options in the
terminal. The options below can be combined; for example,
`--permutations 199 --min-n 50` uses both settings.

| Option | Default | What it does | When to change it |
|---|---:|---|---|
| `inputs` | All original `*_Wide_*.csv` files in this folder | Optional file names placed after the command; analyzes only those files. | Use when rerunning one gene, for example `python3 analyze_covariation.py NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93.csv`. |
| `--permutations N` | 999 | Number of shuffled datasets used to calculate P values. | Use 199 for a quicker exploratory run; use 999 or more for final results. |
| `--seed N` | 20260923 | Starting number for the random shuffles. The same seed gives reproducible results when all other settings are unchanged. | Change only to repeat the analysis with a different random permutation sample. Record the seed if changed. |
| `--min-n N` | 20 | Minimum number of sequences in a subtype with nonmissing calls at both tested positions. | Raise it, for example to 50, to focus on better-sampled subtypes; lower it only with caution because small groups give unstable estimates. |
| `--min-minor N` | 3 | Each position must have at least this many calls outside its most common amino acid within a subtype. | Raise it to exclude rarer variation; lower it only for exploratory work. A position with no variation cannot covary. |
| `--alpha N` | 0.05 | Q-value cutoff for writing the detailed subtype-MI and amino-acid-combination files. It does not change the main scan results. | Use a smaller value, such as 0.01, for a stricter detailed-results list. |
| `--min-combination-count N` | 3 | Minimum observed count for an amino-acid combination to appear in the detailed combination file. | Raise it to hide very rare combinations; lower it to inspect rarer combinations cautiously. |
| `--near-duplicate-differences K` | Off | Adds a sensitivity data set with one representative per complete-linkage cluster whose members differ at no more than `K` positions. | Try `K = 1` after exact deduplication. This can be slow for large subtypes. |
| `--include-all-position-scan` | Off | Adds a scan of every pair of positions in the full gene range. | Use only after the focused RAS scans, because it performs many more tests and takes longer. |

4. Review the primary RAS--RAS results.

   For each gene and data set, open the file ending in
   `_covariation_<dataset>_ras_ras.csv`. There are two rows for every position
   pair:

   - `sequence_weighted`: large subtypes contribute more to the overall MI.
   - `subtype_balanced`: each eligible subtype has equal weight.

   Prioritize pairs with a small `q_value`, a positive `adjusted_mi_bits`, and
   adequate `n_sequences` and `n_subtypes`. `adjusted_mi_bits` is observed MI
   minus mean permuted MI; it is the primary effect-size column.

5. Rerun the final analysis with 999 permutations.

   Use the same seed and thresholds for the final result. The default command
   already uses 999 permutations, at least 20 paired observations per subtype,
   and at least 3 calls outside the modal amino acid at each position.

   ```sh
   python3 analyze_covariation.py
   ```

   The program performs RAS--RAS and RAS--non-RAS scans as separate
   Benjamini--Hochberg testing families, independently for each weighting
   method. RAS--non-RAS scanning is substantially larger and may take time.

6. Inspect the amino acids responsible for significant pairs.

   For any scan with q values at or below 0.05, inspect the accompanying:

   - `*_significant_subtype_mi.csv`: MI and sample size for each contributing
     subtype.
   - `*_significant_aa_combinations.csv`: subtype-specific amino-acid-pair
     counts, independence-expected counts, observed/expected lift, and
     log2-lift.

   This determines whether a result occurs in several subtypes or is driven by
   one subtype, and identifies combinations such as `30R` with `93H`. Do not
   describe a position-pair result as one universal amino-acid combination
   without this inspection.

7. Compare the exact-deduplication sensitivity analysis.

   The program produces both `all_sequences` and `exact_deduplicated` results.
   Exact deduplication retains one copy of each fully observed profile within a
   subtype, while retaining all incomplete profiles. Compare effect sizes,
   q-values, and the number of retained records in
   `*_covariation_dataset_summary.csv`. A signal that disappears after a large
   reduction in one subtype should be interpreted cautiously.

8. Optionally perform the near-duplicate sensitivity analysis.

   This retains one representative from each complete-linkage cluster within
   subtype. Every pair of complete profiles in a cluster differs at no more
   than the specified number of positions. Start with one difference:

   ```sh
   python3 analyze_covariation.py --near-duplicate-differences 1
   ```

   This step can be slow for large subtypes. Specify the threshold before
   examining its association results.

9. Optionally conduct the all-position scan only after reviewing the focused
   scans.

   This examines every position pair and is a separate, much larger
   multiple-testing family.

   ```sh
   python3 analyze_covariation.py --include-all-position-scan
   ```

10. Report results as associations, with appropriate limitations.

    Report the position pair, weighting method, adjusted MI, P value, q value,
    contributing subtype sample sizes, leading amino-acid combinations, and
    results before and after deduplication. Stratification prevents differences
    in subtype amino-acid frequencies from creating an association, but it
    does not remove shared ancestry or study/geographic clustering within a
    subtype. Therefore, these results identify amino-acid associations, not
    direct biological interactions.

## Output naming

For an input named `NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93.csv`, the
primary outputs are named like:

```text
NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93_covariation_all_sequences_ras_ras.csv
NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93_covariation_exact_deduplicated_ras_ras.csv
NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93_covariation_all_sequences_ras_nonras.csv
NS5A_Profile_Accession_AA_Calls_Wide_Pos24_93_covariation_dataset_summary.csv
```

The script ignores prior files with `_covariation_` in their names when it
automatically finds input CSV files, so it is safe to rerun without treating
an earlier result as an input.

## Presentation-ready primary report

After the covariation analysis has finished, create one concise table for the
primary RAS--RAS results across all genes:

```sh
python3 build_primary_covariation_report.py
```

This writes `covariation_primary_ras_ras_report.csv`. By default, it retains a
position pair when its q value is at most 0.05 under either
`sequence_weighted` or `subtype_balanced` analysis, and puts both summaries on
one row. The full scan files remain the complete record of every tested pair.

It also writes `covariation_primary_ras_ras_report.xlsx`, with separate `NS3`,
`NS5A`, and `NS5B` sheets in one workbook, plus a `Read Me` sheet explaining
the columns and weighting methods. The workbook has no color styling; it uses
frozen headers, filters, readable column widths, and the same displayed
precision as the CSV.
The report generator requires the `openpyxl` package. In this project, run it
with the included environment:

```sh
./.venv/bin/python build_primary_covariation_report.py
```

Values smaller than 1 are formatted with two significant digits without
scientific notation. For example, `0.01234` is reported as `0.012`, `0.5` as
`0.50`, and `0.001` as `0.0010`.

To include every pair, including non-significant results:

```sh
python3 build_primary_covariation_report.py --include-nonsignificant
```
