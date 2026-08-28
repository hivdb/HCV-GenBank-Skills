# Genotype and subtype amino-acid predictability

## Purpose

For every resistance-associated amino-acid position, determine whether knowing a
sequence's genotype improves prediction of its amino acid, and whether knowing
its subtype provides additional predictive information after genotype is known.

The primary unit of analysis is one amino-acid position. Positions must not be
combined into a single haplotype for this primary analysis, because each has a
different amino-acid distribution and coverage. A haplotype analysis may be
reported separately as a secondary analysis.

## Input population

Use the merged complete-profile workbook produced after profile construction:

- NS3: `outputs/comet-NS3/17_merge-subtype-complete-profiles/NS3_Subtype_CompleteProfiles_Merged.xlsx`
- NS5A: `outputs/comet-NS5A/17_merge-subtype-complete-profiles/NS5A_Subtype_CompleteProfiles_Merged.xlsx`
- NS5B: `outputs/comet-NS5B/17_merge-subtype-complete-profiles/NS5B_Subtype_CompleteProfiles_Merged.xlsx`

For each row it supplies:

- `Subtype`: the final subtype assignment;
- gene-specific position column (`NS3Position`, `NS5APosition`, or
  `NS5BPosition`): the amino-acid position;
- `NumSeqsIncludingPosition`: final-profile sequences that cover that position
  for the subtype;
- `AminoAcid` and `CountWithAA`: the observed amino-acid distribution.

Derive genotype from subtype (for example, `1a` belongs to genotype `1`). Only
final-profile accessions represented by this workbook are included. At each
position, exclude sequences without coverage at that position.

The primary estimand is a random accession with a non-ambiguous amino acid at
the position: subtypes are weighted by their analyzed amino-acid counts. A
secondary, clearly labelled analysis may give each adequately sampled subtype
equal weight; that instead describes a random subtype.

## Calculation for one position

Let $A$ be the amino acid at the current position, $G$ genotype, and $S$
subtype. For every amino acid $a$, convert counts to probabilities:

$$
P(A=a)=\frac{\operatorname{CountWithAA}(a)}{\operatorname{NumSeqsIncludingPosition}}.
$$

Calculate the Shannon entropy (in bits):

$$
H(A)=-\sum_a P(A=a)\log_2P(A=a).
$$

This is the uncertainty in the amino acid before any genotype or subtype is
known. It is zero when every covered sequence has the same amino acid, and
increases as the position becomes more diverse.

Next calculate the entropy separately within each genotype, then take an
analyzed-sequence-weighted average:

$$
H(A\mid G)=\sum_g P(G=g)H(A\mid G=g).
$$

This is the uncertainty remaining after genotype is known. The information
provided by genotype is:

$$
I(A;G)=H(A)-H(A\mid G).
$$

Then calculate entropy separately within every subtype and take an
analyzed-sequence-weighted average:

$$
H(A\mid G,S)=\sum_{g,s}P(G=g,S=s)H(A\mid G=g,S=s).
$$

This is the uncertainty remaining after both genotype and subtype are known.
The extra information supplied by subtype once genotype is known is:

$$
I(A;S\mid G)=H(A\mid G)-H(A\mid G,S).
$$

Report the percentage of genotype-conditioned uncertainty resolved by subtype:

$$
100\times\frac{H(A\mid G)-H(A\mid G,S)}{H(A\mid G)}.
$$

If $H(A\mid G)=0$, report this percentage as not applicable: genotype already
fully predicts the amino acid at that position.

## Interpretation

- $I(A;G)$ is the amino-acid uncertainty removed by knowing genotype.
- $I(A;S\mid G)$ is the additional uncertainty removed by knowing subtype
  after genotype is already available. This is the central result.
- A subtype percentage of 0% means subtype adds no prediction beyond genotype.
- A subtype percentage of 100% means subtype resolves all uncertainty that
  remained after knowing genotype.

Report, for every position: covered sequence count, eligible subtype count,
$H(A)$, $H(A\mid G)$, $H(A\mid G,S)$, $I(A;G)$, $I(A;S\mid G)$, and
the percentage resolved by subtype.

## Safeguards and secondary analyses

Include every subtype present in the final subtype profile that has at least one
covered, non-ambiguous amino-acid observation at the current position. Report
each subtype's coverage so readers can see when estimates for very small groups
are unstable. A minimum-coverage sensitivity analysis may be reported
separately, but it must not replace the all-subtypes primary analysis.

Use bootstrap confidence intervals for information measures when uncertainty is
needed. Formal tests across multiple positions require false-discovery-rate
control. True cross-validated log loss needs accession-level data; the merged
workbook has aggregated counts and therefore supports the entropy analysis but
not accession-level cross-validation.
