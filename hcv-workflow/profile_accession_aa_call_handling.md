# Accession-level amino-acid call export

`*_Profile_Accession_AA_Calls.csv` is emitted by the
`merge-subtype-complete-profiles` step.  It contains one row per usable
amino-acid call for each accession and position that was included in complete
profile construction.  The `GT_Consensus` call is selected from that
accession's genotype profile: highest `PctWithAA` wins, and exact ties are
resolved by lexical one-letter amino-acid order.  A position without a
genotype-profile call is reported as `X` for `GT_Consensus`.

Only the 20 standard one-letter amino acids and `*` are emitted in `AA`.
`*` is retained as an observed stop codon; it is not treated as a mixture.
`X` (unknown translation) and any other non-standard symbol are omitted from
the accession-level file, matching the complete-profile rule that unknown
calls do not contribute an amino-acid count.  The step summary records the
number of omitted symbols.

The pipeline's translated input uses one amino-acid character per codon and
does not encode amino-acid mixtures.  Consequently, symbols such as `B`,
`Z`, `J`, or textual forms such as `A/S` have no defined mixture expansion;
they are treated as non-standard calls and omitted rather than being split or
silently assigned to a standard amino acid.

## Ambiguous DNA codons and amino-acid mixtures

An ambiguous nucleotide codon can encode several amino acids.  For example,
`CNN` expands to the codons for `H`, `L`, `P`, `Q`, and `R`.  This is a known
multi-amino-acid mixture, not missing sequence data.

At present, the COMET translation code converts every codon containing a
non-`ACGT` nucleotide to `X`.  Thus `CNN`, `AAN`, and other ambiguous codons
are indistinguishable in the AA profile from an unresolved amino-acid call.
The current RAS ambiguity audit reports these as `XCount`; it does not yet
expand the nucleotide codon into its possible amino-acid set.

When mixture-aware profile support is added, it must derive the possible
amino acids from the corresponding `NASequence` codon and retain the set at
that one position (for example, `CNN` as `H/L/P/Q/R`).  It must not insert
multiple characters into `AASequence`, because that would shift every
following amino-acid position.  A mixed call should be counted once for
position coverage and retained as a mixture in downstream reports; the
combined profile must not select one component as the accession's sole call.
