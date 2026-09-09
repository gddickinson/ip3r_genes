# Appendix C — the sensitivity tables

Four results in this thesis rest on a threshold that somebody chose, and in
each case the count was recomputed across the whole range of that threshold
and committed. This appendix says where those tables are, what each axis is a
knob on, and what moving it costs.

The principle behind all four is the same, and it is worth stating once. **A
threshold-dependent result must say which settings manufacture it, and it must
walk the threshold across its own measured uncertainty before any invented
value.** A robustness claim asserted in prose is not checkable; a grid is.

## C.1 The loss count

`results/loss_counts/sensitivity_matrix.tsv` — 384 rows.

Four axes. **Coding**: the family-level presence character against the
paralogue-resolved one. **Evidence**: an ordered eight-rung ladder, each rung
named after what it refuses, whose first three rungs are the reconstruction
calibration's own measured gap edges and midpoint rather than round numbers.
**The cyclostome rule** on and off. **The contiguity bar** on and off. Every
cell is computed under three branch-length schemes.

The count is zero at the operating point on both codings. The family-level
coding manufactures a loss in 2 of 32 settings and the paralogue-resolved one
in 18, up to 45 loss edges.

`manufactured_losses.tsv` (144 rows) says the same thing at cell resolution:
every genome × cell that ever reads absent under any setting, the state and
rule that held it at the operating point, and the loosest setting that
releases it. That table is what makes the sensitivity claim checkable one cell
at a time rather than one count at a time.

**The branch-length axis is carried and reported as invariant.** Parsimony
counts edges and cannot read a length, so no scheme changes any count — and a
matrix that quietly dropped that axis would be indistinguishable from one that
had tested it and found nothing.

## C.2 Copy number

`results/duplication/copy_sensitivity.tsv` — 56 rows.

Every copy count in the duplication chapter recomputed at seven coverage bars.
A duplication claim that survives only one bar is a claim about the bar, and
the teleost result is the one this table exists to protect: it holds across
the range.

## C.3 Gene architecture

`results/gene_architecture/architecture_sensitivity.tsv` — 24 rows.

Exon counts, spans and intron statistics recomputed at six coverage bars. The
result the table protects is that the exon count is conserved and the span is
not, which is a comparison between two quantities measured on the same loci
and therefore robust to the bar by construction — but the bar decides *which*
loci, and the table says so.

## C.4 Annotation states

`results/annotation_audit/state_sensitivity.tsv` — 21 rows.

Every locus state re-counted across the whole range of the completeness bar.
The audit inherits that bar from the sweep rather than deriving a second one,
so this table is what validates the inherited value: it shows where each state
boundary sits in the distribution and what a different bar would have done to
the headline.

## C.5 Two thresholds with no sensitivity table, and why

**The reconstruction bar** has none, because the evidence ladder in C.1 *is*
its sensitivity analysis: its first three rungs are the two edges of the gap
the calibration measured and the midpoint between them, and moving across the
whole gap changes no cell in 927.

**The identity floor for what counts as a locus** outside the vertebrates has
none, because the calibration that produced it reported that neither identity
nor coverage separates its two populations cleanly, and a sensitivity grid
around a threshold whose calibration refused to separate would be a grid
around nothing. What was done instead was to add a second, independent axis —
every recorded cluster scored against both family profiles — and to report
both statistics.
