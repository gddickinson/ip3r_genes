# The thesis package

`python scripts/s25_assemble.py` builds everything here and exits non-zero on
a failed guard. The same work is reported as a paper in `manuscript/`, and
this directory holds its long form: about 64,000 words across 16 chapters and
5 appendices, with 103 figures and 80 references.

## What is written by hand, and what is generated

**Hand-written:** the 27 numbered chapter source files, and
`chapter_rules.md`. Nothing else.

**Generated, on every build:**

| file | stage | what it is |
|---|---|---|
| `thesis.md` | `stitch` | the assembled document, do not edit |
| `itpr_family_thesis.pdf` | `pdf` | the typeset document, 173 pages |
| `figures/` | `figures` | every placed figure, copied from the results tree, png + pdf |
| `figure_manifest.tsv` | `figures` | chapter, number, slug, source, drawn width, SHA-256 |
| `chapter_assignment.tsv` | `assign` | every results directory → chapter, with the rule |
| `control_inventory.tsv` | `controls` | the negative controls, parsed out of the test modules |
| `production_stats.tsv` | `production` | how the project was made, measured from the repository |
| `reference_audit.tsv` | `refs` | every reference added for the thesis, resolved live |
| `claims_check.tsv` | `claims` | every load-bearing number re-verified |
| `guard_check.tsv` | `guards` | every build guard broken on purpose, and what it said |

`reference_baseline.txt` is frozen: the reference keys that existed before this
document was written. Anything cited that is not in it must carry a verified
audit row.

## The stages

    guards     break every build guard on purpose; check each fires
    assign     commit the chapter grouping (rule T6)
    controls   derive the control inventory from the test modules
    production measure how the project itself was made
    refs      audit every new reference against a live record
    figures   copy every placed figure from the results tree
    claims    re-verify every number, and that its chapter states it
    stitch    concatenate, number the figures, resolve the citations
    pdf       typeset, and read the log for dropped glyphs

    python scripts/s25_assemble.py --only claims
    python scripts/s25_assemble.py --from stitch
    python scripts/s25_assemble.py --refs-offline    # re-audit from the archive

## What fails the build

The build fails on any of the following. A results directory that is neither
assigned to a chapter nor explicitly excluded. A missing chapter file. A
declared figure with no file, a slug declared twice, a committed figure
neither placed nor excluded, a figure placed in the wrong chapter, or a
reference to a figure that nothing places. A cited key with no reference row,
a new reference with no verified audit row, or a reference audited and never
cited. A number that no committed table produces, or a number declared in the
ledger that its own chapter never states. And a glyph the document font cannot
set.

All of them are exercised on every build by `scripts/s25_test_guards.py`,
which runs in a sandboxed copy and verifies that it altered no committed file.

## Reading it

Chapters 1–2 set out the receptor and audit the literature the project started
from. Chapters 3–5 build the search and report the family's range. Chapters
6–9 are its history. Chapters 10–12 are the protein. Chapter 13 is the
archive. Chapter 14 is the methods written as argument, and carries the
project's recorded decisions. Chapter 15 discusses the whole, and Chapter 16
describes how the project was carried out with Claude Code.

Appendix A is the correction list, B the control inventory, C the sensitivity
tables, D the reference audit, and E the chapter/paper assignment with the
guard log.
