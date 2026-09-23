# Paper rules

*These rules were written and committed before any results directory was
assigned to a paper. The assignment they produced is `paper_assignment.tsv`,
and the build (`scripts/s26_assemble.py`) enforces the rules that can be
checked mechanically.*

The single manuscript reports seven bodies of results in 17,807 words, which
leaves each one a few hundred. S26 reports the same work as a series of papers.
The risk runs the other way from the thesis's. The thesis risked becoming the
paper with padding added. A series risks becoming one paper cut into slices,
where each slice leans on the others for its controls, its scope or its point.
These rules are meant to catch that.

The thesis (S25) ran first and committed `thesis/chapter_assignment.tsv`. This
series starts from that table and records every place where it departs from
it, with the reason (`s26_assign.DEPARTURES`). A departure that is not declared
fails the build.

## The six rules

**P1: a paper answers one question.** The question can be written as one
sentence with no *and* in it. Each paper declares its question in its
configuration module, and the build rejects a question that contains the word
"and". A paper may make several claims, as long as each one is a step towards
answering its question.

**P2: a paper measures its own controls.** Every control that one of a paper's
central claims depends on is measured inside that paper, not cited from a
sibling. Each paper lists its controls with the committed table that holds each
one. The build fails if a control's table sits in a results entry that belongs
to a different paper. A paper whose negative control lives in another paper is
really a section of that paper.

**P3: each paper declares its own scope.** The paper that uses a denominator
states it in its own methods. Two papers may share a scope, but neither may
take it over by reference. Each paper declares its scope as claims in the
ledger: the genome count, the proteome count, the orthologue count. Those
claims must appear in the paper's own text, not in a sibling's.

**P4: a paper has enough figures to show its result.** It has four to seven
main figures, all copied from committed results directories. Enforced: a
grouping with fewer than four is a section, and one that needs more than seven
is two papers. Extended Data and Supplementary figures are not limited by this
rule.

**P5: each result is primary in exactly one paper.** A result is claimed in
exactly one paper. Other papers may cite that paper, and the citations are what
make the set a series rather than a set of slices. Enforced at two levels:

- *Entries.* Every entry under `results/` is either primary in exactly one
  paper or explicitly excluded with a reason. A new results directory that
  nobody has assigned fails the build.
- *Figures.* A paper may place only figures from results entries that are
  primary in that paper. To show a sibling's result, it cites the sibling.

**P5s: split a directory only when a control forces it.** Assignment is by
directory unless one directory's tables answer the central questions of two
papers. The split is allowed only when P2 forces it: the file group moved is a
control that another paper's central claim cannot stand without. It is then
done by file glob, and the build requires every file in the directory to match
exactly one group.

**P6: each paper survives on its own.** For each paper, write down what it
claims if none of the others is ever published. A paper that cannot answer that
is a slice and must be merged with a neighbour. The answer is committed as
prose in each paper's configuration and in `papers/README.md`. Its mechanical
half is enforced: the answer names the ledger claims it rests on, and each of
those claims must have its source in a results entry that is primary in that
paper. If a standalone claim rests on a sibling's table, the paper does not
stand alone.

## The citation graph between papers

Papers cite one another in their text as `{paper:<id>}`, and the build draws
the dependency graph from those citations rather than from a declaration. Two
rules follow from P6:

- **The graph has no cycles.** Two papers that each need the other are one
  paper.
- **A paper cites only papers that come earlier in submission order**, so every
  paper can be read when it appears. The submission order is declared in
  `s26_lib.SERIES`, and the build checks it against the graph.

## One ledger

Every load-bearing number in the series is declared once, in one ledger
(`scripts/s26_claims*.py`), with a `papers` column naming every paper that
states it. When a number is quoted in two papers, it is checked once. It cannot
disagree between them, because there is only one row. Numbers the manuscript or
thesis ledgers already declare are carried by identifier, not restated. The
same number declared twice under two identifiers is a build error. Every paper
listed against a claim must actually state its value (the thesis's D73 rule),
so the paper column cannot grow a paper that never quotes the number.

## What these rules produced

See `paper_assignment.tsv` (one row per results entry, or per file group for a
split directory, with the paper it is primary in, the rule that placed it, and
the rule that kept it out of the paper it came closest to), and
`papers/README.md`.
