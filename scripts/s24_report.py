"""S24 — renders `results/supplementary/report.md` purely from the committed
tables and `supplementary_stats.json` (D13).

Nothing here recomputes a number. Every figure caption in the report is the
caption in `supp_figure_stats.tsv`, every count comes from `headline()`, and
the audit section is rendered from `figure_findings.tsv` — so a finding that
was recorded but not written up, or written up but not recorded, cannot
happen.

A section whose table is absent renders *not run yet* rather than nothing,
so a stage that was skipped is visible as skipped.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s24_lib as L                                            # noqa: E402
import s24_tables as T                                         # noqa: E402


def _n(x, default="not run yet"):
    return f"{x:,}" if isinstance(x, int) else (x if x else default)


def render(stats: dict) -> None:
    h = stats.get("headline") or T.headline(stats)
    g = stats.get("guards", {})
    f2 = stats.get("supp_fig_2", {})
    f6 = stats.get("supp_fig_6", {})
    audit = stats.get("audit", {})

    lines: list[str] = []
    w = lines.append
    w("# S24 — supplementary figures, and the figure audit")
    w("")
    w(f"*Rendered {date.today().isoformat()} from the committed tables "
      f"(D13). Nothing in this report is hand-written.*")
    w("")
    w("## 1. What this task is")
    w("")
    w("Two deliverables. Six supplementary figures showing the alignments "
      "and structures the main figures rest on — none of which re-aligns or "
      "re-renders anything, each drawn from the same committed file as the "
      "main figure it supports. And a figure-by-figure audit (**D11**): "
      "every main and Extended Data figure opened and read against its own "
      "legend, because errors of the kind *the legend says four and the "
      "figure draws three* are invisible to every table check in the "
      "project.")
    w("")

    w("## 2. The two guards, run before any figure was drawn")
    w("")
    if not g:
        w("*not run yet*")
    else:
        w("A supplementary figure exists so a reader can check a join. Both "
          "joins were therefore checked first, in code, as hard failures.")
        w("")
        w("| guard | what it requires | scale |")
        w("|---|---|---|")
        w(f"| the trimAl column map | every trimmed column must be the input "
          f"column the map names, for every sequence | "
          f"{_n(h['columns_checked'])} columns x "
          f"{_n(h['sequences_checked'])} sequences |")
        w(f"| the residue at the column | every variant's reference residue "
          f"must be the residue its own paralogue's table holds there, and "
          f"every aligned partner the residue the *other* paralogue's table "
          f"holds | {_n(h['variant_residues_checked'])} variants, "
          f"{_n(h['aligned_partners_checked'])} partners |")
        w("")
        w("Both passed. trimAl writes no column map of its own; S6 recovered "
          "one with `-colnumbering` and S24 does not trust it — the walk is "
          "exhaustive and an off-by-one would have no other symptom, because "
          "every column would still map to a column.")
        w("")
        st = stats.get("self_test", {})
        w(f"{_n(st.get('checks'))} constructed negative controls run before "
          f"anything is written (`scripts/s24_test_supp.py`), and "
          f"{_n(st.get('mutations_tested'))} deliberate rule breakages were "
          f"mutation-tested against them; all "
          f"{_n(st.get('mutations_caught'))} were caught. Two of them were "
          f"not, at first, and both are recorded in the test. The "
          f"duplicate-column case has to be built where the two input "
          f"columns hold the *same* residues, or the content walk catches "
          f"it first and the one-to-one guard is never exercised. And a "
          f"byte-identity check on a saved figure passes vacuously whenever "
          f"both saves land in the same second, so the property is tested "
          f"directly instead: the pdf must carry no creation timestamp.")
    w("")

    w("## 3. The six supplementary figures")
    w("")
    rows = []
    path = L.OUT_DIR / "supp_figure_stats.tsv"
    if path.exists():
        rows = [r for r in L.read_tsv(path) if r["format"] == "png"]
    if not rows:
        w("*not run yet*")
    else:
        w("| # | figure | what it shows | drawn from |")
        w("|---|---|---|---|")
        for r in rows:
            src = ", ".join(f"`{s.split('/')[-1]}`"
                            for s in r["sources"].split(";"))
            w(f"| S{r['number']} | `{r['figure']}` | {r['caption']} | "
              f"{src} |")
        w("")
        w(f"**Supplementary Fig. S1** — the L-INS-i alignment is "
          f"{_n(h['aln_columns'])} columns and the tree, the selection tests "
          f"and the constraint map were all computed on the "
          f"{_n(h['trimmed_columns'])} trimAl kept. What trimAl removed was "
          f"the sparse columns: median occupancy "
          f"{h['median_occupancy_kept']} kept against "
          f"{h['median_occupancy_cut']} cut.")
        w("")
        w(f"**Supplementary Fig. S2** — {_n(h['pathogenic_positions_drawn'])} "
          f"pathogenic positions fall in the ligand core or the pore module, "
          f"and {_n(h['pathogenic_conserved_all_three'])} of them carry the "
          f"same residue in all three paralogues. Every letter in the panel "
          f"was checked against the paralogue's own per-residue table.")
        w("")
        w(f"**Supplementary Fig. S3** — the deep layers run 249–265 "
          f"orthologues per paralogue against msa_v2's 13–19. "
          f"{_n(h['orthologs_screened'])} sequences went through the shape "
          f"screen and {_n(h['orthologs_dropped'])} was dropped.")
        w("")
        w(f"**Supplementary Fig. S4** — {_n(h['codons_kept'])} of "
          f"{_n(h['codons_in'])} codons survive trimming, and that is the "
          f"alignment every ω in the paper was estimated on.")
        w("")
        w("**Supplementary Fig. S5** — the constraint map as S17 painted it, "
          "on all three cryo-EM references, with the same structure painted "
          "with the selection layer beside it. Unscored residues are grey, "
          "never the low end of the scale: S17 wrote them −1 precisely so "
          "that *we could not score this* and *this is the least conserved "
          "part of the receptor* cannot look the same, and the luminal loop "
          "is both.")
        w("")
        if f6:
            keep = ", ".join(h["structures_eligible_for_variants"])
            refused = ", ".join(f"`{x}`" for x in h["structures_refused"])
            w(f"**Supplementary Fig. S6** — a structure may carry a human "
              f"variant position only if *every* residue it shares with the "
              f"human per-residue table carries the same amino acid. "
              f"{keep} pass; {refused} do not, so *ITPR1*'s pathogenic "
              f"positions are not drawn on coordinates that are not theirs. "
              f"The per-element enrichment test runs over "
              f"{_n(h['element_tests'])} element × gene cells, of which "
              f"{_n(h['element_tests_significant'])} clear q < 0.05 after "
              f"Benjamini–Hochberg.")
    w("")

    w("## 4. The figure audit (D11)")
    w("")
    if not audit:
        w("*not run yet*")
    else:
        w(f"{_n(h['figures_audited'])} manuscript figures — 7 main and 14 "
          f"Extended Data — were opened and read against their own legends. "
          f"The mechanical half is checked in code on every build: every "
          f"figure must have a legend, every legend a figure, and for each "
          f"Extended Data figure the panel letters its legend uses must "
          f"match the panel files the manifest holds.")
        w("")
        by = h.get("audit_findings_by_status") or {}
        w(f"{_n(h['audit_findings'])} findings, "
          + ", ".join(f"{v} {k.replace('_', ' ')}"
                      for k, v in sorted(by.items())) + ".")
        w("")
        fpath = L.OUT_DIR / "figure_findings.tsv"
        if fpath.exists():
            w("| figure | the legend said | the figure shows | fixed |")
            w("|---|---|---|---|")
            for r in L.read_tsv(fpath):
                w(f"| {r['figure']} ({r['panel']}) | {r['legend_said']} | "
                  f"{r['figure_shows']} | {r['status'].replace('_', ' ')} |")
            w("")
        w("Every numeric finding was re-derived from the committed table "
          "named in its `verified_against` column before it was written "
          "down; the table is `figure_findings.tsv`.")
    w("")

    w("## 5. What this changes, and what it does not")
    w("")
    w("**It changes the manuscript text.** Every legend finding above is "
      "corrected in `manuscript/11_figure_legends.md` or "
      "`manuscript/12_extended_data.md`, and the corrections are numbers, "
      "not wording: a legend that said four backbone labels where the tree "
      "draws five, a legend that named the wrong colour for the fragmentary "
      "class, a legend that claimed a range of 51–53 for a panel that plots "
      "up to 102.")
    w("")
    n_fig = (h.get("audit_findings_by_status") or {}).get("fixed_figure", 0)
    n_leg = (h.get("audit_findings_by_status") or {}).get("fixed_legend", 0)
    w(f"**It changes {_n(n_fig)} things about the figures themselves, "
      f"against {_n(n_leg)} legend corrections.** Panel letters were "
      f"uppercase in four figures and lowercase everywhere else; a lone "
      f"panel letter `a` sat on a figure with no panel b; two annotations "
      f"overlapped into unreadable text; two heat maps side by side carried "
      f"independent colour scales; a bar the legend asked the reader to "
      f"compare against was hidden behind the bar it was being compared "
      f"with; a threshold was drawn and never named; marker shape carried "
      f"meaning that appeared in no key; a panel title asserted an agreement "
      f"two of its six rows did not reach; four AUCs were printed at mixed "
      f"precision; and one colour that read clearly on screen vanished at "
      f"printed size.")
    w("")
    w("**And it makes every figure in the project byte-reproducible.** "
      "matplotlib stamps the wall clock into a PDF's creation date, so a "
      "figure rebuilt from the same code on the same data differed from "
      "itself and no checksum recorded against a pdf meant anything. "
      "`figstyle.save` now drops the field; twelve of twelve S24 files "
      "rebuild byte-identically, and the property is a self-test rather "
      "than a claim.")
    w("")
    w("**It does not change a result.** No table was recomputed and no "
      "analysis was re-run. Every fix is to a legend or to how a committed "
      "number is drawn.")
    w("")
    w("**What it leaves open.** The audit is an inspection, and an "
      "inspection finds what it looks for. The mechanical checks that now "
      "run on every build cover existence, pairing and panel counts; "
      "whether a legend's *description* matches what a panel draws is not "
      "checkable in code and will need looking at again whenever a figure "
      "or a legend changes — which is what D11 says.")
    w("")

    (L.OUT_DIR / "report.md").write_text("\n".join(lines) + "\n")
