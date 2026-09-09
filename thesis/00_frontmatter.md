# A genome-scale census of the inositol 1,4,5-trisphosphate receptor family

*Range, origin, retention, constraint and record quality across 503 genomes
and 7,691 reference proteomes.*

**Author:** Claude (Opus 5, Anthropic)

**Correspondence:** george.dickinson@gmail.com

The analyses this document reports were carried out over 35 working sessions
in the repository it is built from, and the document was written from their
committed outputs. The correspondent above directed the project and is
the point of contact for it, since the author has no address of its own.

---

## Abstract

The inositol 1,4,5-trisphosphate receptor (IP₃R, gene family *ITPR*) is the
endoplasmic reticulum's ligand-gated calcium-release channel. In vertebrates
it exists as three paralogues, ITPR1, ITPR2 and ITPR3, whose origin,
distribution and retention had never been measured against a declared search
space. This thesis reports that measurement, and it gives equal space to the
instruments that made it possible.

The family was enumerated across 7,691 reference proteomes and 503 genomes, of
which 309 are vertebrate and 194 are not. The ryanodine receptors were
searched alongside at every stage. That control is not decoration: RYR1
carries all four of the Pfam signatures that define an IP₃ receptor, is twice
its length, and is returned by every query this project ran. Separating the
two families is therefore treated as a positive test at every stage rather
than as an assumption, and where an instrument could not make that test the
result is reported as undecided rather than as an absence.

Three results follow. First, the family is ancestrally eukaryotic and has been
lost repeatedly and independently. It is absent from all 384 land-plant and
all 1,353 Dikarya reference proteomes, from Apicomplexa, Microsporidia,
Rhodophyta, diatoms and four other fungal phyla, and from every archaeal and
bacterial proteome swept. All 35 clade-level absences hold at assembly level
in genomes where a measured positive control recovers a comparably long,
deeply conserved gene. Second, the three vertebrate paralogues are not lost at
all. Across 927 genome by paralogue cells, none reaches the absence state,
Dollo parsimony places no loss anywhere on the tree, and what is reported
instead is the set of analytical settings that would manufacture one. Third,
the family's public record is substantially worse than the family itself. Most
full-length *ITPR* protein records carry no usable gene symbol, and three
quarters of the genes demonstrated here cannot be reached by any
protein-database search.

The duplications that produced the three paralogues are dated and their
asymmetry is explained. ITPR2 and ITPR3 are sisters. The first split sits on
the vertebrate stem and the second on the gnathostome stem. The surviving
two-round paralogon links run through ITPR1, which is also the only paralogue
that was doubled and retained after the teleost genome duplication, and the
only one held under roughly twice the purifying selection of its sisters. The
three paralogues share most of their intron positions in every genome that
carries them, an enrichment of nearly two orders of magnitude over a null
drawn from the alignment itself, while the ryanodine receptors, which carry
every diagnostic domain of the family, share one. Constraint mapped onto the
cryo-EM channel identifies the gate and the selectivity filter as the least
changeable elements, places the ten measured IP₃ contacts outside the Pfam
domain named after the ligand, and finds the constrained unit at the ligand
site to be a pocket rather than the contacts themselves.

Three things here are usable rather than only informative. The archive audit
produces a list of 297 coordinate-level corrections with the evidence a
curator would need to act on each. The per-residue constraint map is scored as
a variant classifier and identifies which of four conservation layers a
clinical resource for this family should quote. And five of the results are
about method rather than about receptors, including a measured false-negative
rate for genome-scale orthologue sweeps and a measurement of what a protein
bait panel's breadth is actually worth.

## How this document is organised, and why it is long

The same work is also reported as a paper of 17,807 words. That paper is the
right length for its claims and the wrong length for its instruments, and the
difference between the two is the reason this document exists.

Almost every number in the paper rests on a threshold, and almost every
threshold in this project was measured rather than chosen. Examples include
the intron length above which a gene reads as a fragment, the alignment margin
that separates two families of the same architecture, the coverage bar a
reassembled pseudogene has to clear, and the identity floor below which a
locus is not a locus. Each of those measurements has a result of its own, and
several of them changed the answer. In a paper each is a Methods sentence.
Here each sits in the chapter that first depends on it, together with what it
was measured against and what guessing would have cost.

The same applies to what did not work. A thesis is the only format in which a
measured dead end is worth its page, and this one carries several: an
exhaustive substitution-model scan that was started, timed and abandoned at 11
of up to 1,232 models; an attribution threshold inherited from a sister
project and overturned because it came from a family half as similar to itself
as this one; a duplication detector that reached a specificity of 0.16 before
it was scoped correctly; and a synteny-caller threshold that optimising the
obvious quantity would have set at the loosest value on offer. None of these
is filed in a chapter of failures. Each sits in the chapter that owns the
instrument it was measured on, because that is where it explains why a number
can be believed.

Chapters 1 and 2 set out the receptor and audit what the literature could be
trusted to say about it before any of this was measured. Chapters 3 to 5 build
the search and report the family's range. Chapters 6 to 9 cover its history:
the tree, the duplications, the gene's own architecture, and a retention
result that required more care to state than any positive finding in the
project. Chapters 10 to 12 cover the protein: selection, constraint, and the
one part of it the ryanodine receptors do not share. Chapter 13 audits the
archive. Chapter 14 is the methods, written as argument rather than as
procedure, and it carries the project's decision log and its constructed
negative controls. Chapter 15 discusses what the whole says and what it does
not, and Chapter 16 describes how the project itself was carried out.

## What is checked in this document, and by what

Nothing here is typed twice. Every figure is copied from the results directory
that committed it and is never re-plotted. Every load-bearing number is
declared with the committed table it comes from and the operation that
recovers it, and the build re-reads those tables on every run. Every citation
is a stable key resolved against one reference table, and the bibliography is
rendered from that table rather than written. The claims ledger additionally
requires each declared number to appear in the chapter that declares it, which
is what stops a ledger being padded with checks the text never makes.

Running `python scripts/s25_assemble.py` rebuilds the whole document and exits
non-zero on a missing figure, a missing chapter, a cited key with no reference
row, a reference added without an audit, a glyph the document font cannot set,
or a failed claim. Every one of those guards is broken on purpose on every
build to confirm that it fires. Chapter 14 explains how.

## Where the data, the code and the outputs are

Every number in this document is recovered from a committed table, and every
one of those tables is in the repository this document is built from, under
`results/`, one directory per analysis task. The analysis code is under
`scripts/`, one file per task, and the search application it was built on is
under `src/`. The exact version of every external program the project ran is
recorded in `results/toolchain_manifest.txt`, and the versions of the public
databases queried are recorded in each task's own statistics file.

Bulk inputs are deliberately not in the repository, because 503 genome
assemblies and 7,691 reference proteomes are not a reviewable artefact. What
is committed instead is the manifest that declares exactly which ones were
used, with the rule that admitted each, so the download is reproducible from a
table rather than from a description. `manuscript/deposit_notes.md` lists each
excluded bulk class together with the command that regenerates it.

## How this project was carried out

This project was carried out by Claude Code, which is Anthropic's Claude
running as an agent in a terminal with access to the file system, the shell
and the network. Claude Code wrote every line of the analysis software, chose
and calibrated every threshold, ran every search, generated every table and
figure, recorded every methodological decision, and wrote this document and
the paper that accompanies it.

The human collaborator set the goal, provided the machine and the storage,
ran the small number of commands that need an interactive login, and read and
corrected the drafts. He did not write the code, select the genomes, choose
the statistical tests or draft the text.

That division of labour is unusual enough to be worth describing rather than
leaving for a reader to infer, and it is also relevant to how the work should
be judged, because an agent that can run for hours unattended can produce a
great deal of plausible output. **Nearly every methodological rule in this
thesis exists because an autonomous worker needs to be stopped from
convincing itself.** The measured thresholds, the positive controls carried
through every stage, the negative controls that must fire, the reports
rendered from tables rather than written, and the claims ledger that re-reads
every number are all answers to the same question: how does anyone, including
the agent, tell whether what it just produced is true?

Chapter 16 describes the whole arrangement in detail, covering how the work
was organised into tasks, how the literature was surveyed and audited, how the
decisions about what to download and what to measure were made, how the
analyses ran, and how the results, the figures and these documents were
produced.
