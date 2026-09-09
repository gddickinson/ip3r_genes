# 16. How this project was carried out with Claude Code

## 16.1 What Claude Code did, and what it did not

This project was carried out by Claude Code, which is Anthropic's Claude
running as an agent in a terminal with access to the file system, the shell
and the network. It reads and writes files, runs commands, inspects their
output, and decides what to do next.

**Claude Code did the following.** It designed the project plan and wrote it
down as a task ledger with dependencies. It wrote all 373 Python files under
`scripts/`, which run to over a hundred thousand lines, together with the
9,243 lines of the search application under `src/`. It audited the
literature the project started from, claim by claim, against primary
sources. It decided which genomes and proteomes to download, wrote the
manifests that declare them, and fetched 624 GB of assemblies and 60 GB of
proteomes. It chose every threshold in the work and, in most cases, measured
it rather than choosing it. It ran every search, alignment, tree, selection
test and structural comparison. It produced all 444 committed result tables
and all 111 committed figures. It recorded 83 numbered methodological
decisions as it made them. And it wrote the literature review, the
manuscript, this thesis and every one of the 28 rendered task reports that
stand behind them, which together run to 109,250 words.

Every figure in that list is measured rather than recalled. The build
re-derives them from the repository on every run and commits them as
`thesis/production_stats.tsv`, so they move as the project does and this
chapter cannot describe a repository that no longer exists.

**The human collaborator did the following.** He stated the goal, which was a
publication-grade genome-scale census of one gene family. He provided the
machine, the external storage and the network access. He ran the small number
of commands that require an interactive login or elevated privileges. He read
the drafts and corrected them, including two corrections to this document that
are recorded in the decisions log. And at the start of each working session he
said, in effect, continue.

He did not write the code, choose the genomes, select the statistical tests,
set the thresholds or draft the text.

**What follows is a description of that arrangement rather than a defence of
it.** An agent that runs for hours unattended can produce a great deal of
plausible output, and the interesting question is not whether it can but how
anybody, including the agent, can tell whether the output is true. Almost
every methodological rule in this thesis is an answer to that question, and
§16.9 sets out what the answer looks like in practice.

## 16.2 The work was organised as a ledger of one-task sessions

The organising artefact is a single file, `PUBLICATION_ROADMAP.md`, which
Claude Code wrote in the first session and has amended in every session since.
It contains four things: a statement of the goal and how the claim is scoped;
a session protocol; a task ledger; and a decisions log.

**The session protocol is a fixed procedure that starts and ends every working
session.** At the start it opens a progress dashboard, pulls the repository,
verifies that the external storage is attached and stops if it is not, then
picks exactly one task, which is the single row already in progress, or
otherwise the topmost pending row whose dependencies are all complete. It
announces that task and works only that task until its completion criteria
pass. At the end it updates the ledger and the results column, appends to the
session log, appends a plain-language entry to a findings document, refreshes
the repository's README, commits and pushes.

That protocol exists because of a specific failure mode. An agent with a long
list of things to do will do several of them partially, and a project of this
size dies of half-finished tasks rather than of hard ones. **Working one task
to its stated completion criteria, and splitting a task in the ledger when it
proves too large, is the rule that kept it finishable.** Several tasks were
split this way, and each split is visible in the ledger as two rows with the
reason recorded.

**The ledger is 36 tasks, of which 34 are complete.** Each row carries the
task, its dependencies, its status with a date, and a results column holding
the load-bearing numbers and the paths they came from. The two incomplete rows
are the human-gated deposit, which needs a person to create a public archive
record, and the paper series, which is the next task.

Alongside the ledger, a briefs document holds a detailed specification for
every row: the goal, the steps, the completion criteria, the outputs, and, for
the later rows, what would make that task a failure. Claude Code wrote those
briefs in advance of doing the work, which matters because a brief written
after the fact describes what happened rather than what was intended.

**35 sessions were logged across 9 working days**, the first on 2026-08-18 and
the most recent on 2026-09-09. The session log is a running technical record
of what ran, what resulted and what is next. It is written for the next
session rather than for a reader, and it is the mechanism by which an agent
with no memory between sessions resumes work that is already in progress.

## 16.3 The literature was surveyed and then audited claim by claim

The first substantive task was not a search. It was to write down what the
project was taking as known, and then to check it.

Claude Code assembled a background document in which every statement carries a
tag: `[db]` if it can be re-derived from a database with a stated query,
`[lit]` if it rests on a published source, and `[open]` if it is a question
the project answers. It then extracted nineteen atomic claims from the
`[lit]`-tagged statements and checked each against its primary sources,
and re-queried every `[db]` number against the live database it came from.

Chapter 2 reports the outcome. Twelve claims were verified, three were
qualified, two were corrected and one was demoted to an open question. The
demoted one mattered most: the background asserted that the IP₃ and ryanodine
receptor families had independently expanded to three paralogues each, no
primary source establishes it, and leaving it in place would have made three
later chapters into a confirmation of something already assumed.

**Claude Code then wrote a full literature review**, running to 137 references
and twelve generated figures, as a separate document. That review is not
decoration. It is where the project's claim about what is already known lives,
and it is built by a script from numbered source files with citations resolved
against one reference table, so its text and its bibliography cannot drift
apart.

Two habits from that stage run through everything after it. **A number is
quoted with the query that produces it**, and **a claim is tagged by how far
it can be trusted**, so that a later chapter can tell the difference between a
fact it inherited and a fact it measured.

## 16.4 What to download was decided by rule and committed as a manifest

The project needed 503 genome assemblies and 7,691 reference proteomes. Which
ones was not a judgement call, because a judgement call cannot be audited.

**The denominator is declared before the search runs, and it is derived from
the data rather than hand-listed.** For the vertebrate scope, Claude Code
wrote two rules: one best assembly per vertebrate order, ranked by a stated
function, and a margin rule with four positive tests that reads the census to
find species whose protein records are missing, fragmentary or absent. The
union is 309 genomes. For the non-vertebrate scope it wrote five rules under
one principle, which is to sample most finely where the negative claim is, and
one of those rules re-derives its own target list from a committed presence
table. That rule found three absences the earlier summary had never named,
including a metazoan clade with no record at all.

**Both manifests are committed tables with a reason column on every row**, so
the scope rebuilds itself if the census changes, and a reader can see why any
given genome is in the study.

The downloads themselves were run by a script Claude Code wrote for the
purpose, which verifies every archive against its published checksums,
retries download, extraction and verification together because a truncated
archive surfaces as a checksum failure rather than a download error, and can
operate in a fetch, search and purge cycle so that a 624 GB scope can be
processed at one genome of peak storage.

## 16.5 What to analyse was decided by a menu written before the analyses ran

Between building the search and running the analyses, Claude Code wrote a
document listing what the harvested data could answer. Each entry names its
inputs, its method, its deliverable and the caveat that would invalidate it,
and the entries became ledger rows.

Writing that menu before doing the work is the difference between a study and
a fishing expedition. **An analysis chosen after seeing the data is chosen
partly because of what the data shows**, and the menu is the record that these
were not.

The menu is also where the project's scope discipline is visible. Several
entries were written and then answered as underpowered rather than dropped,
because a measurement that did not reach its question is a result about the
evidence and an omission is indistinguishable from a section nobody ran.

## 16.6 The analyses ran unattended, and each one is a driver with stages

Every analysis task is a set of small modules and one driver script with
ordered, resumable stages. The pattern is the same throughout: a driver with
`--only` and `--from` flags, stages that are individually re-runnable, a
negative-control suite that runs before anything is written and refuses the
build if it fails, and a report and figures rendered at the end from the
tables the stages committed.

**Three properties made unattended operation possible.**

Stages are resumable, so a run that is interrupted resumes at the stage that
failed rather than from the beginning. The genome sweep is resumable per
genome, and the selection suite is resumable per job.

**A failed stage does not stop the stages after it that do not depend on it**,
and the report marks an unfinished section as unfinished. Waking up to a
partial analysis that says which parts are partial is worth more than waking
up to nothing.

And expensive work is cached at the level that makes a re-run cheap, with the
rule that **a cache may hold what a parser found and never what a rule
decided**, so that changing a rule invalidates the verdicts without
invalidating the parse.

The largest single run was the selection suite, which ran for 22.7 hours
unattended. The genome sweep processed 309 vertebrate assemblies and then 194
non-vertebrate ones. **15.2 compute-hours are recorded across the eight search
channels**, measured on clean recomputes rather than on cached re-runs,
because a cached run makes a budget look free.

During long runs the drivers write a small progress file that a dashboard
reads, so the human collaborator could see what was happening without
interrupting it.

## 16.7 Thresholds were measured rather than chosen, which is where the agent's judgement went

The part of this project that most needed judgement was not which analysis to
run. It was where to put each threshold, and the answer Claude Code arrived at
early and applied throughout is that **a threshold should be measured against
a population the threshold has not already filtered, and reported with its
separation**.

Chapter 14 sets out the rules in full. What is worth saying here is that this
is where the agent's own reasoning is most visible, and it is visible because
it was written down at the time. **87 numbered decisions** are recorded in the
roadmap, each with the incident that produced it. Several changed an answer.
An attribution margin inherited from a sister project was measured and
overturned, because it came from a family whose paralogues are half as similar
to each other as this one's. An intron-length parameter was measured across
eleven species because it is not a performance setting but a threshold on
whether a gene reads as a fragment. A calibration refused to separate its two
populations and sent a chapter back to find the contamination in its decoy.

**The decisions log is also where the agent recorded being wrong.** A kill
criterion written for the project's central hazard turned out to have a
sensitivity of zero when it was finally scored. A duplication detector reached
a specificity of 0.16 before it was scoped correctly. A self-test overwrote
the committed table it was testing. Each is in the log with the fix and the
rule that generalises it.
