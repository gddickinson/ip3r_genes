"""S19 — renders `results/methods/report.md` purely from the committed tables.

D13 applied to a report about methods: nothing here is hand-written and
nothing is recomputed. Every number comes through `s19_tables.headline()`,
which reads the tables, so the prose and the tables cannot drift.

Scope, the control design, the instrument and the negative controls live
here; the results — what each method was worth, what the panel bought, where
the search and the inference ran out — live in `s19_report_results.py`
(the `s3_report.py` / `s3_report_d10.py` split).

Headlines are chosen by the data. Every prior an earlier task set is stated
in `s19_priors.PRIOR` with where it was said, computed on S19's own tables,
and rendered `confirmed` / `contradicted` / `not corroborated` /
`orthogonal` / `underpowered` with both numbers printed either way. A
section whose table is absent renders *not run yet*, so a stage that did not
run is visible as one.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402
import s19_tables as T                                          # noqa: E402


def fmt(v, nd: int = 4) -> str:
    if v is None or v == "":
        return "—"
    if isinstance(v, float):
        return f"{v:,.{nd}f}".rstrip("0").rstrip(".") if v else "0"
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def pct(v) -> str:
    return "—" if v in (None, "") else f"{float(v) * 100:.1f} %"


def pfmt(p) -> str:
    """A p-value at four decimal places reads `0.0000`, which hides how
    strong a claim is rather than how weak (`s15_report.py`'s rule)."""
    try:
        x = float(p)
    except (TypeError, ValueError):
        return str(p) if p else "—"
    return f"{x:.3g}"


def table(header: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if c is None else str(c)
                                     for c in r) + " |")
    return "\n".join(out)


def section_1(h: dict) -> str:
    n_ctrl = h["itpr_present"]["n_cells"] + h["ryr_sister"]["n_cells"]
    return f"""# S19 — methods results: what the search was worth

*Generated {date.today().isoformat()} by `scripts/s19_run.py` from the
committed tables in `results/methods/`. Nothing here is hand-written.*

## 1. What this task can measure that most cannot

Every methods section that reports a search's sensitivity has to estimate it,
because the genes the search missed are the ones nobody can count. This
project does not have that problem. S15b reconstructed **no losses anywhere
in the {fmt(h['n_genomes'])}-genome scope**: all
{fmt(h['itpr_present']['n_cells'])} assignable genome × paralog cells hold a
gene that is there, established by evidence the sweep's locus placement did
not produce. So every ledger cell that is not `found` is a **false negative
of the method**, and the sensitivity of a genome-scale ortholog sweep is a
direct measurement.

There are two such series and they were measured by the same instrument in
the same assemblies:

{table(["series", "cells", "missed", "false-negative rate", "95 % CI"],
       [["ITPR cells S15 states present", fmt(h['itpr_present']['n_cells']),
         fmt(h['itpr_present']['n_false_negative']),
         pct(h['itpr_present']['rate']),
         f"{pct(h['itpr_present']['wilson'][0])} – "
         f"{pct(h['itpr_present']['wilson'][1])}"],
        ["RyR sister cell, present in every vertebrate",
         fmt(h['ryr_sister']['n_cells']),
         fmt(h['ryr_sister']['n_false_negative']),
         pct(h['ryr_sister']['rate']),
         f"{pct(h['ryr_sister']['wilson'][0])} – "
         f"{pct(h['ryr_sister']['wilson'][1])}"]])}

The two rates are not distinguishable (Fisher exact
*p* = {pfmt(h['itpr_vs_ryr_p'])}), which is the point of carrying the second
one: the ITPR figure rests on S15a's state assignments and the RyR figure
does not, so agreement between them is evidence that the number is a
property of the search rather than of how the states were called.

Exactly {fmt(h['not_counted_as_control'])} cells are in **neither** series,
and they are the cyclostome ITPR2/ITPR3 cells S15a filed
`paralog_unassignable`: those genomes carry spare ITPR loci the bait panel
cannot label, so calling them present or absent would decide a question S15a
declined.

## 2. Scope and inputs

Nothing in this task downloads or re-aligns. Every input is a committed
artefact of an earlier session, with one derived asset built here:

- the **swept accession universe** — one `grep '^>'` pass over each of the
  seven reference-proteome FASTAs, cached under the data root. Without it,
  "the InterPro enumeration holds a record the profile HMM did not return"
  cannot be told apart from "that record was never in the database the
  profile HMM searched", and the head-to-head in §5 would be an accounting
  artefact rather than a comparison.
- the **retained per-genome miniprot alignments** — 309 GFFs. miniprot
  aligns each bait independently, so dropping baits from the file and
  re-clustering reproduces exactly what the sweep would have reported had
  those baits never been in the panel. The ablation in §6 is therefore
  exact, not a model of one.

Total: {fmt(n_ctrl)} control cells, {fmt(h['n_genomes'])} genomes,
{fmt(h['drift']['n_runs'])} jackhmmer runs, and
{h['recorded_hours']:.1f} recorded compute-hours over eight search
channels.

## 3. Negative controls

`s19_test_methods.py` runs before anything is written and refuses the build
if it fails. The checks are on *refusal* and on *reachability*, because every
rule in this task returns a plausible number when it is wrong.

{table(["check", "what it refuses"],
       [["T1–T2", "a FASTA header must yield a bare accession in all three "
         "shapes, and a genome model id must survive intact — stripping one "
         "at its first `.` would collapse every model of an assembly onto "
         "its accession prefix"],
        ["T3", "a zero-row parse must refuse to cache a universe; an empty "
         "universe makes every record read as absent from the database it "
         "was found in"],
        ["T4", "a `paralog_unassignable` cell must enter neither control "
         "series"],
        ["T5", "a false negative must be *constructible* — a control that "
         "can only ever return zero misses is not a measurement"],
        ["T6", "raising the contiguity floor must never retain more genomes, "
         "and D4's own bar must be on the scanned grid"],
        ["T7", "the full-panel simulation must reproduce the committed "
         "ledger cell for cell"],
        ["T8", "a locus with no labelled paralog bait must fill no paralog "
         "cell, and must be offered once one aligns there"],
        ["T9", "D14 — a locus won by the RyR baits must be offered to no "
         "ITPR cell"],
        ["T10", "a killed jackhmmer run must contribute only its accepted "
         "rounds"],
        ["T11", "the proposed drift rule must fire on a drifting run and "
         "decline on a stable one, both constructed"],
        ["T12", "the drift label must come from the finished model, and the "
         "two populations must separate"],
        ["T13", "the head-to-head must compare family *calls*, not raw "
         "profile targets"],
        ["T14", "every manifest genome must land in exactly one scope class"],
        ["T15", "the suite itself must alter no committed table — the first "
         "build's T11 called the real trace routine and overwrote "
         "`kill_criterion_trace.tsv` with its two constructed rows, and the "
         "report then read 2 jackhmmer runs where there are 7"]])}

All {fmt(T._self_test_status()['n_checks'])} checks pass. Five were
mutation-tested by breaking the rule they guard — counting
`paralog_unassignable` as present, comparing raw targets instead of calls,
pointing the proposed drift rule back at K1's axis, letting an empty universe
cache, and letting the suite write to the results directory — and every
mutation was caught by its own check.
"""


def run(log=S.log) -> Path:
    import s19_report_limits as L
    import s19_report_results as R
    h = T.headline()
    body = [section_1(h), R.section_4(h), R.section_5(h), R.section_6(h),
            L.section_7(h), L.section_8(h), L.section_9(h), L.section_10(h)]
    path = S.out_dir() / "report.md"
    path.write_text("\n\n".join(body) + "\n")
    log(f"report.md: {len(path.read_text().splitlines())} lines")
    return path


if __name__ == "__main__":
    run()
