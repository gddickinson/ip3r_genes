# Appendix D. The reference audit

The audit table is `thesis/reference_audit.tsv`, with one row per reference
added for this document, giving the source database, the identifier queried,
the phrase the title had to carry, the title returned, the year, the journal,
the first author, the verdict, and what the thesis rests on the reference for.

**76 references were added and all 76 verified.** Nine required a
corrected identifier before they could be.

## D.1 What the rule for a new reference requires

The literature review's 137 references were audited when it was built, and
this document inherits them as a frozen baseline, committed as
`thesis/reference_baseline.txt` before any new reference was added. Every key
cited by this thesis that is not in that baseline is a new reference, and the
build fails on a new reference with no verified audit row.

The rule for a new reference is the review's rule unchanged: it is **audited
on entry**. What that means in practice is that nothing bibliographic is
typed. Each new reference is declared by an identifier alone, together with a
distinctive phrase that must appear in the title that identifier resolves to.
The build resolves the identifier against a live record, using PubMed for
anything indexed there and Crossref otherwise, admits it only if the title
carries the phrase, and then writes the bibliographic row from the fetched
record.

So a mistyped identifier fails the audit instead of quietly putting a
different paper in the bibliography, and no author, title, year or journal can
be wrong, because none of them is entered by hand.

Responses are archived, so the audit re-runs offline.

## D.2 What the audit caught

**Nine of the 76 identifiers resolved to entirely different papers.**

A gene-duplication inference algorithm's identifier returned a Bayesian
phylogenetics program. A reconciliation method's returned a paper on
statistical challenges in real-time gene expression. A morphological-character
likelihood model's returned a paper on species names in phylogenetic
nomenclature. A vertebrate ancestral-genome reconstruction's returned a paper
on structured RNAs in the human genome. A zebrafish gene-cluster paper's
returned a paper on a poly(ADP-ribose) polymerase at telomeres. Two timetree
papers, a gene-loss review and a Dollo-parsimony paper resolved elsewhere or
not at all.

One of the nine was a different failure and is worth separating: the
identifier was right and the phrase was wrong. That one was corrected by
adjusting the phrase, because the record the identifier returned was the
intended work.

Every one of the other eight would have entered a bibliography looking
entirely normal, with the right journal era, a plausible author list, and a
number nobody checks. That is the argument for the rule in one sentence: **a
bibliography assembled by hand has no way to notice.**

## D.3 What the new references are

By audit class, there are 12 tools this project actually ran, whose versions
are recorded in the toolchain manifest; 21 methods, models or statistical
procedures implemented or applied; 10 databases, cited at the release the
project used; and 33 substantive biological or evolutionary claims.

The tool and method references are the bulk of the addition and the reason it
was needed at all. The manuscript cites 29 references and none of them is a
tool, because a paper can name its software in a Methods section and be
understood. A chapter that explains why an alignment was run single-threaded,
or why a branch-site likelihood ratio is tested against a mixture, is citing a
method and must say whose.

## D.4 What is not audited here

The 137 inherited references are not re-audited. They were audited once,
against their primary sources, in the session that built the literature
review, and 19 atomic claims resting on them were scored, giving 12 verified,
3 qualified, 2 struck and 1 downgraded to an open question. Chapter 2 reports
that audit and what it changed.

Re-auditing them would be a second session's work for no new information. What
this appendix adds is that the new references have been through a check the
inherited ones were not: their identifiers have been resolved against a live
record and required to return the paper the declaration names.
