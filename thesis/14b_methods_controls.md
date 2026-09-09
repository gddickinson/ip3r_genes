## 14.9 The constructed negative controls are worth taking as a body of work

The paper reports these in one Methods paragraph and a count. There are more
of them than there are paragraphs in the paper's Results, and the ones that
mattered mattered a great deal.

There are **437 named constructed checks across 19 suites.** The inventory in
Appendix B is derived from the test modules themselves rather than typed, so a
control that is deleted disappears from the appendix and one that is added
appears in it. The unit counted is a named check, meaning one assertion
reported by name as the suite runs, and where a suite groups its checks into
named test functions instead, that is what is counted and the table says so,
because the two units are not the same size and adding them together silently
would be worse than either.

Every suite runs **before anything is written**, and a failure refuses the
build.

### What the controls are checks on

Almost none of them checks that a routine returns the right answer. They check
that it **refuses**, and that it **can act**.

That emphasis is not stylistic. In this project almost every rule returns a
plausible number when it is wrong. An interval routine that scores against
gene spans instead of coding blocks credits an annotation with sequence it
never called. A minus-strand intron read without reverse-complementing both
ends is non-canonical, and getting that wrong makes every minus-strand gene in
the sweep look broken, which is a systematic signal that would discredit
exactly the alignments a chapter is defending. A codon model will fit a
frame-shifted alignment and return a perfectly plausible rate. A duplication
detector that ignores the family call measures paralogy at a specificity of
0.16 and a sensitivity of 0.997, which looks like an excellent instrument.

**A rule that fires on nothing is indistinguishable from a rule that cannot
fire.** Several results in this thesis are zeros: no locus encodes the same
part of the protein twice, no block pair is a frameshift pair, no cell reaches
the absence state, the merge fires on none of 2,146 loci, and the chimera
screen rejects none of 57 candidates. Each of those zeros has a constructed
control that builds the case the rule is supposed to catch and requires it
caught, which is the only thing that makes the zero a measurement.

### Six controls changed something

**The chimera screen's own control failed on its first run, and the test was
wrong.** It asked the domain-envelope rule to reject a truncation, which the
envelope cannot and should not do, because a protein truncated to 40 % of its
length aligns 100 % of itself to the profile in one clean segment. Truncation
is the length band's job and fusion is the envelope's.

**A reassignment rule's positive and negative halves both failed**, and the
rule was wrong rather than the test. It climbed the tree from a mislabelled
tip until some group dominated, which hands a tip to whichever clade is
largest three nodes up.

**A ranked-table order-invariance check caught real hash-seeded
nondeterminism** in a committed table.

**A within-protein control was being selected by a name-prefix test** that
swept 225 residues of the element the ligand question is about into the
control set.

**A self-test damaged the artefact it tested.** One suite called the real
routine on constructed rounds and wrote its two synthetic rows over a
committed table, after which the report read two runs where there are seven,
and nothing failed. Two rules follow: any routine a control exercises takes a
write flag and the control passes it false, and the suite checks the checksum
of every committed table before and after it runs.

**A mutation test passed on broken code**, twice, in one suite. A
duplicate-column case has to be built where the two columns hold the same
residues, or an earlier guard catches it and the one being tested is never
exercised. And a byte-identity check on a saved figure passes vacuously
whenever both saves land in the same second, so the property is now tested
directly.

### Mutation testing is the check on the checks

Where a suite's own report records it, the suite was **mutation-tested**,
meaning the rule it guards was deliberately broken and the suite required to
catch it. Across the tasks that recorded this, every deliberate breakage was
caught by the check responsible, giving a total of dozens of mutations across
the analysis chapters, with each one named in its task's committed statistics
rather than summarised.

That is the check on the checks. A suite that has never been shown to fail is
a suite whose passing means nothing.

## 14.10 One rule cannot be automated: look at the figure

One decision in this project cannot be automated, and it is stated as an
instruction rather than a rule: **look at the figure.**

Every main and supplementary figure was opened and read against its own
legend, because errors of the kind "the legend says four and the figure draws
three" are invisible to every table check in the project. A figure is
generated from a committed table, so its data cannot drift. Its description
can, and nothing downstream notices.

**Twenty-three figures were inspected and 26 findings recorded**, comprising
16 legends corrected and 10 figures corrected. They are recorded as data, one
row per finding, with what the legend said, what the figure shows, the
committed table the correction was re-derived from, and which of the two was
changed.

The findings are worth characterising because they are all of one kind. One
key described two colours where the figure draws four. One legend said
"lineages grouped by kingdom" where three of the groups are not kingdoms.
Another named the class with the most non-blue area where the figure shows a
different class ahead on the fraction it plots. One figure carried a bold
panel letter with no second panel to go with it. And one legend called four
labels the extended paralogue clades where the figure carries five nodes. None
of these is a data error, and none would have been caught by re-reading a
table.

**The mechanical half of that inspection is now automated and runs on every
build.** Every figure must have a legend, every legend a figure, and each
multi-panel figure's legend letters must match its panel files. A later pass
added the other mechanical half, requiring every figure to be cited by a
sentence of the body and the citations to run in ascending order, because
before it two figures were cited by nothing at all and the methods figure was
numbered last and first mentioned in the third results section.

## 14.11 Two joins behind the supplementary figures are checked in code

A supplementary figure exists so a reader can check a join, so both joins were
checked first, in code, as hard failures.

**The trimming column map** requires that every trimmed column be the input
column the map names, for every sequence, across 1,797 columns and 134
sequences, walked exhaustively. The trimming program writes no column map of
its own. One was recovered from a reporting flag, and it is **not trusted**,
because an off-by-one would have no other symptom: every column would still
map to a column, and every residue claim downstream would be renumbered with
nothing to show for it.

**The residue at the column** requires that every variant's reference amino
acid be the residue its own paralogue's per-residue table holds there, and
every aligned partner the residue the other paralogue's table holds, across
1,780 variants and 2,699 partners.

Both passed, which is why a residue-level figure in this thesis prints
letters.

**A structure may carry a position from another numbering only if it earns
it.** A structure carries a human variant position only when every residue it
shares with the human table carries the same amino acid, and two of four
candidates fail, so one paralogue's pathogenic positions are simply not drawn
on coordinates that are not theirs.

## 14.12 How this document is built, in eight stages

Running `python scripts/s25_assemble.py` runs eight stages in dependency order
and exits non-zero on any failure.

The **guards** stage breaks every build guard on purpose and checks that each
fires, and §14.13 and Appendix E describe it.

The **assign** stage commits the chapter grouping and fails on a results
directory that is neither assigned to a chapter nor explicitly excluded, which
is the no-leftovers rule enforced rather than asserted, and which means a
future task's results cannot be silently left out of this thesis.

The **controls** stage derives the negative-control inventory from the test
modules themselves.

The **refs** stage audits every reference added for this document, and §14.14
is that audit.

The **figures** stage copies every placed figure from the results directory
that committed it, in both formats, and fails on a missing file, a slug used
twice, or a committed figure that the thesis neither places nor explicitly
excludes.

The **claims** stage re-verifies every load-bearing number against its
committed table through the manuscript's own engine, and additionally requires
each declared number to **appear in the chapter that declares it**, which is
what stops a ledger being padded with checks the text never makes.

The **stitch** stage concatenates the chapters, numbers the figures per
chapter in order of first placement, resolves the citation keys, and fails on
a missing chapter, a placed figure the map does not declare, a declared figure
no chapter places, a figure placed in a chapter other than its own, or a
figure reference that resolves to nothing.

The **pdf** stage typesets the result and reads the typesetter's own log for
dropped glyphs.

## 14.13 Every guard is broken on purpose, on every build

The requirement was that each guard be broken on purpose once. Doing that by
hand proves it for one afternoon, so it is instead a build stage that runs
first, and it carries two properties that make it worth the file.

**Each case declares a fragment the guard's own message must contain**, so a
case that starts failing for a different reason is a failure rather than a
pass. The suite's own first run had four cases returning non-zero for the
wrong reason: two path bugs in the harness, one guard reading a module-level
constant the sandbox could not redirect, and one mutation applied to a list
the driver had already copied.

**It runs against a sandboxed copy of the thesis directory and checks the
SHA-256 of every committed file before and after**, because this project has
already had a self-test overwrite the committed table it was testing.

All 15 fire, and `thesis/guard_check.tsv` records what each said.

## 14.14 The reference audit caught nine citations that named the wrong paper

The literature review carries 137 references, audited when it was built. This
document needed 58 more, almost all of them methods and tools, and the rule
for them is the review's rule unchanged: **a reference added and not audited
is worse than no reference, because it launders an assumption into a
citation.**

So nothing bibliographic is typed. Each new reference is declared by an
identifier alone, meaning a PubMed identifier or a digital object identifier
where the work predates indexing, together with a distinctive phrase that must
appear in the title that identifier resolves to, and a note of what the thesis
rests on it for. The build resolves each identifier against a live record,
admits it **only if the title carries the declared phrase**, and then writes
the bibliographic row from the fetched record.

**That check caught nine of 76.** Nine identifiers written from memory
resolved perfectly well, to entirely different papers. A gene-duplication
inference algorithm's identifier returned a Bayesian phylogenetics program. A
reconciliation method's returned a paper on real-time gene expression. A
morphological-character likelihood model's returned a paper on species names
in phylogenetic nomenclature. A vertebrate ancestral-genome reconstruction's
returned a paper on structured RNAs.

Every one of those would have entered a bibliography looking entirely normal,
with the right journal era, a plausible author and a number nobody checks. **A
bibliography assembled by hand has no way to notice**, and that is the whole
argument for the rule.

The audit is committed as a table with one row per new reference, giving the
source database, the identifier queried, the phrase required, the title
returned, and the verdict. All 58 now read *verified*.
