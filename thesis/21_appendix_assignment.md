# Appendix E. How the results were assigned to chapters, and the build guards

## E.1 What the assignment table records

The assignment is committed as `thesis/chapter_assignment.tsv`, with one row
per entry under `results/`, giving the chapter it is primary in, the rule that
placed it, what it is, and the chapter it might plausibly have gone to instead
with the reason it did not.

There are **37 results directories assigned across 13 chapters, with one
excluded.**

The rules are written out in `thesis/chapter_rules.md`, which was committed
before any chapter prose was written, and four of the seven are enforced by
the build rather than asserted.

**Rule T4 says a results directory is primary in exactly one chapter.** It is
enforced by construction, because the assignment is a mapping.

**Rule T6 says a grouping that only holds leftovers is an appendix.** It is
enforced in the strong direction: every entry under `results/` is either
assigned to a chapter with a rule and a reason, or named in an exclusion list
saying why it is not a result. A directory that exists and is unassigned fails
the build, which also means a future task's results cannot be silently left
out of this thesis. One entry is excluded, namely the dashboard's progress
panel, which is session state.

**Rule T7 says every figure is copied from a committed results directory and
never re-plotted.** It is enforced, in that a missing file, a slug used twice,
or a committed figure the thesis neither places nor explicitly excludes all
fail the build. **The thesis places 103 figures**, every one copied from the
results directory that committed it in both formats, with the source and the
SHA-256 of each recorded in `thesis/figure_manifest.tsv`. Eight figures are
explicitly excluded, all of them the search application's own bundle plots
from smoke-test runs, which are drawn from whatever one search returned rather
than from a committed analysis table.

**Rule T3 allows four to twelve figures per results chapter.** It is checked
by inspection against the figure manifest rather than by the build. The
introduction, the methods chapter and the discussion are exempt as exposition,
and the methods chapter carries one figure.

## E.2 The paper series should expect one departure from this assignment

This task ran before S26, so this table is the one S26 starts from. There is
one place where a paper series will have to choose differently, and it is
recorded here as a disagreement rather than resolved.

**The methods results, meaning the measurement of what the search was worth,
sit entirely in Chapter 4** with the search they measure. S26's own starting
proposal splits them across two papers: a retention paper taking the
false-negative rate, and an archive paper taking the contribution half.

Both are right about their own format, and the reason they differ is not a
matter of taste.

**The thesis's reason.** The false-negative rate is a property of the search,
in that 15.2 % of cells whose gene is independently known present were not
found by the ledger. A document that explained the search in one chapter and
then withheld the number saying what that search was worth until five chapters
later would be hiding the instrument's error bar from the chapter that builds
the instrument.

**The series' reason.** A retention claim without its own false-negative rate
is a paper whose central control lives somewhere else, which is exactly what
S26's second rule forbids, since every control a paper's claims rest on must
be measured inside that paper.

The rule that forces the choice is the one both documents share, namely that a
result is primary in exactly one place, and it is shared deliberately so that
a disagreement between the two groupings is visible rather than invisible.

**What a reader holding both documents should take from it.** The thesis and
the series are not the same object cut two ways. A chapter can be a stage in a
longer construction and a paper cannot. Where the two groupings differ, the
difference is a fact about the two formats, and this appendix is where it is
recorded rather than smoothed over.

## E.3 Every guard is broken on purpose, and this is what each one said

Every guard in the build is broken on purpose on **every** build, by
`scripts/s25_test_guards.py`, which runs as the first stage and writes
`thesis/guard_check.tsv`. **All 18 fire with the message they are supposed
to**, and the suite verifies that it altered no committed file while doing it,
because a self-test that damages what it tests is a failure mode this project
has already had once.

Each case declares a fragment the guard's own message must contain, so a case
that starts failing for a different reason is a failure rather than a pass.

| guard | how it was broken | what it said |
|---|---|---|
| the assignment (T6) | a results directory removed from the assignment | `results/<name> is neither assigned to a chapter nor listed in ASSIGNMENT_EXCLUDED (rule T6)` |
| the assignment (existence) | an assignment pointing at a directory that does not exist | `results/<name> is assigned but does not exist` |
| a missing chapter | a chapter file moved aside | `missing chapter file thesis/<name>` |
| a missing figure | a declared figure's source file moved aside | `<slug>: missing <path>.png` |
| a duplicated slug | one figure declared twice | `slug '<slug>' is declared more than once` |
| an unexcluded figure | a committed figure removed from the exclusion list | `<path>.png is committed under a chapter's results tree but is neither placed in the thesis nor listed in UNPLACED` |
| an unplaced figure | a declared figure's placement deleted from its chapter | `<slug> is declared for chapter N but no chapter places it` |
| a borrowed figure | a figure placed in a chapter other than its own | `thesis/<file> places <slug>, which is declared for chapter N` |
| a dangling figure reference | a reference to a figure that is never placed | `{fig:<slug>} is referenced but the figure is never placed` |
| an undeclared reference | a citation key with no row in the reference table | `cited keys with no row in references.tsv: <key>` |
| an unaudited reference | a new key cited with its audit row removed | `cited keys added after the frozen baseline with no verified audit row: <key>` |
| a decorative reference | a key's single citation deleted | `references added for the thesis and never cited: <key>` |
| a failed claim | a claim's expected value altered | `[FAIL] <id> <claim>: expected X, found Y` |
| a padded claim | a claim declared for a chapter whose text does not state its value | `chapter N does not state 'X', so either the number is missing from the text or the claim is padding` |
| a dropped glyph | a character the document font cannot set | `N distinct missing characters, so the document font lacks a glyph the text uses` |
| an em-dash | an em-dash reintroduced into a chapter source | `N em-dash(es) in a chapter source` |
| a legend repeating itself | a sentence duplicated inside one figure legend | `legend <slug> says the same thing twice` |
| a legend repeating its neighbour | a legend given the text of the paragraph beside it | `legend <slug> restates the paragraph after it` |

Four of these are the ones this thesis added over the manuscript's build: the
padded-claim guard and the three prose guards at the foot of the table. The
padded-claim guard closes the ledger in the direction the manuscript's cannot,
because it fails not only when the document states a number no table produces,
but when the ledger declares a check the document never makes. The three prose
guards exist because the two editorial passes on this document each introduced
a defect that every other guard passed over: an em-dash carrying a clause a
reader needed, and a figure legend saying the same thing twice.
