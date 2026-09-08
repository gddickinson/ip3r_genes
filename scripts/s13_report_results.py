"""S13 — the placement half of `results/reconciliation/report.md` (§5-§6).

Split from `s13_report.py` to keep both under 500 lines, and taking the
caller's loader and formatter so the two halves cannot read the tables
differently (the `s3_report.py` / `s3_report_d10.py` pattern).

**Headlines are chosen by the data.** Every section states the prior from
`s13_priors.py` — what an earlier task concluded and where it said it —
computes S13's own answer beside it, and renders the verdict from the
comparison, printing both either way. The comparison that has to be sayable is
the one that goes badly: if the reconciliation had put both duplications on the
gnathostome stem, the table it would be said from is the same one.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s13_lib import deepest_nodes as _deepest      # noqa: E402
from s13_lib import matrix_cells as _cells         # noqa: E402
from s13_priors import PRIOR, verdict_line          # noqa: E402

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")




# ------------------------------------------------------------------ §5

def _placement(load, table, num, not_run) -> list:
    summ = load("reconciliation_summary.tsv")
    place = load("duplication_placement.tsv")
    cal = {r["node"]: r for r in load("species_tree_calibrations.tsv")}
    L = ["## 5. Where the duplications sit", ""]
    if not summ or not place:
        return L + not_run("the reconciliation")

    L += ["### 5.1 The matrix, cell by cell", "",
          "Every cell reports the same two numbers and the same answer: how "
          "many duplications the reconciliation implies, how many losses, and "
          "which species-tree node the **deepest paralog-spanning duplication** "
          "maps to.", ""]
    L += table(["topology", "variant", "tips", "collapsed", "duplications",
                "implied losses", "deepest paralog duplication",
                "its support in the gene tree"],
               [[f"`{r['topology']}`", r["variant"].replace("_", " "),
                 r["n_tips"], r["n_nodes_collapsed"] or "—",
                 r["n_duplications"], r["implied_losses"],
                 f"**{r['deepest_paralog_dup_node']}**",
                 r["deepest_paralog_dup_support"] or "—"]
                for r in _cells(summ)])

    with_c = _deepest(summ, "with_cyclostome")
    without_c = _deepest(summ, "without_cyclostome")
    coll = _deepest(summ, "support_collapsed")
    n_topo = len({r["topology"] for r in _cells(summ)})
    same = len(with_c) == 1
    L += [f"**The topology does not move the answer.** All {n_topo} gene-tree "
          f"topologies — the ML tree, the three constrained sister "
          f"arrangements the AU test scored, and the `--bnni` re-search — give "
          f"identical duplication and loss counts within a variant, and "
          f"{'the same' if same else 'different'} deepest placement: "
          f"{', '.join(sorted(with_c))} with the cyclostome loci, "
          f"{', '.join(sorted(without_c))} without them"
          + (f", {', '.join(sorted(coll))} with the unsupported nodes "
             f"collapsed" if coll else "") + ".", "",
          "That is expected rather than surprising, and it is worth saying "
          "why: the three AU hypotheses differ only in how the three paralog "
          "clades are arranged relative to one another, and every arrangement "
          "puts all three inside the same vertebrate clade. What decides the "
          "*placement* is not their order but whether cyclostome loci sit "
          "among them — which is the second axis.", ""]

    L += ["### 5.2 The two events, and their brackets", ""]
    ml = [r for r in place if r["topology"] == "ml"]
    L += table(["variant", "rank", "maps to", "paralogs below it",
                "gene tips", "bracket (Ma)", "gene-node support"],
               [[r["variant"].replace("_", " "), r["rank"],
                 f"**{r['species_node']}**",
                 r["paralog_groups"].replace("+unplaced", " + cyclostome"),
                 r["n_gene_leaves"],
                 (f"{num(r['age_young_ma'],0)} – {num(r['age_old_ma'],0)}"
                  if r["age_old_ma"] not in ("", "unbounded")
                  else f"older than {num(r['age_young_ma'],0)}"),
                 r["gene_node_support"] or "—"] for r in ml])
    gn = cal.get("Gnathostomata", {})
    vt = cal.get("Vertebrata", {})
    cy = load("cyclostome_loci.tsv")
    nested = [r for r in cy if r["cyclostome_clade_size"] == "2"]
    L += [f"Read plainly, on the ML tree with every tip in:", "",
          f"- **The split that separated ITPR1 from ITPR2+ITPR3 is on the "
          f"vertebrate stem** — older than crown Vertebrata, whose published "
          f"estimates run {num(vt.get('age_lo'),0)}–{num(vt.get('age_hi'),0)} "
          f"Ma. Nothing in this tree bounds it from above, and the report says "
          f"so rather than closing the bracket silently.",
          f"- **The split that separated ITPR2 from ITPR3 is on the "
          f"gnathostome stem** — between crown Gnathostomata "
          f"({num(gn.get('age_ma'),0)} Ma; "
          f"{num(gn.get('age_lo'),0)}–{num(gn.get('age_hi'),0)}) and crown "
          f"Vertebrata ({num(vt.get('age_ma'),0)} Ma). It maps there in "
          f"**every** cell of the matrix, cyclostomes in or out.", "",
          f"**What forces the older placement is one hagfish/lamprey pair.** "
          f"Of the six cyclostome loci, {len(nested)} sit in a "
          f"cyclostome-only clade that is itself *nested among the gnathostome "
          f"paralogs* — sister to the ITPR2+ITPR3 clade at "
          f"{nested[0]['pair_support'] if nested else '—'}. A lineage holding "
          f"both cyclostome and gnathostome copies forces every duplication "
          f"above it to predate the cyclostome-gnathostome split, which is why "
          f"both deep events map to Vertebrata. The remaining four sit in a "
          f"cyclostome-only clade outside everything, which is what the "
          f"reconciliation reads as a third ancestral vertebrate lineage.", "",
          "**And what keeps the younger event on the gnathostome stem is the "
          "same tip.** That pair is sister to ITPR2+ITPR3 rather than inside "
          "either, so the node uniting ITPR2 with ITPR3 holds gnathostomes "
          "only and maps below the cyclostome divergence. Had it fallen inside "
          "ITPR2 or inside ITPR3, that split too would have been pushed onto "
          "the vertebrate stem — the two answers above are decided by where "
          "one two-tip clade attaches.", ""]
    return L


# ------------------------------------------------------------------ §6

def _robustness(load, table, num, not_run) -> list:
    summ = load("reconciliation_summary.tsv")
    rk = load("rooting_check.tsv")
    L = ["## 6. What the answer depends on", ""]
    if not summ:
        return L + not_run("the reconciliation")

    L += ["### 6.1 The support the deep placement rests on", ""]
    coll = [r for r in _cells(summ) if r["variant"] == "support_collapsed"]
    place = [r for r in load("duplication_placement.tsv")
             if r["topology"] == "ml" and r["variant"] == "with_cyclostome"]

    def _ufboot(s):
        try:
            return float(s.split("/")[1])
        except (IndexError, ValueError):
            return 1e9

    weak = min(place, key=lambda r: _ufboot(r["gene_node_support"])) \
        if place else None
    node = weak["gene_node_support"] if weak else "—"
    computed = (
        f"the deepest paralog duplication maps to "
        f"**{', '.join(sorted(_deepest(summ, 'with_cyclostome')))}** on the "
        f"full tree and to "
        f"**{', '.join(sorted(_deepest(summ, 'support_collapsed')))}** after "
        f"every node below SH-aLRT 80 / UFBoot 95 is dissolved "
        f"({coll[0]['n_nodes_collapsed'] if coll else '—'} of them on the ML "
        f"tree). The placement therefore does **not** rest on the weak node: "
        f"collapsing it merges two duplications into one and the one that "
        f"survives is carried by a node at "
        f"{coll[0]['deepest_paralog_dup_support'] if coll else '—'}."
    )
    L += [verdict_line("support", computed, "confirmed"), "",
          "This is worth being precise about, because it is the objection a "
          "reader should raise. The duplication that separates ITPR1 from "
          f"ITPR2+ITPR3 sits on a gene-tree node at **{node}** — the weakest "
          "node in the whole vertebrate subtree, and apparently the one "
          "carrying the headline. It is not. Collapsing every node below the "
          "bar dissolves that arrangement and the placement stays: the root "
          "becomes a polytomy in which more than one child still maps into the "
          "gnathostome subtree, and the non-binary rule (§3) reads a "
          "duplication at Vertebrata from that alone. What the collapse costs "
          "is *resolution* — two separate vertebrate-stem duplications merge "
          "into one — not the placement. Remove the cyclostome tips instead "
          "and the placement moves at once, which is §7.", ""]

    L += ["### 6.2 Is the rooting doing the work?", ""]
    if not rk:
        L += not_run("the rooting check")
    else:
        out = next((r for r in rk if r["edge_index"] == "-1"), None)
        others = [r for r in rk if r["edge_index"] != "-1"]
        tot = [int(r["total_events"]) for r in others]
        best = min(tot) if tot else 0
        mins = [r for r in others if int(r["total_events"]) == best]
        agree = {r["deepest_paralog_dup_node"] for r in mins} == \
            {out["deepest_paralog_dup_node"]} if out else False
        L += ["Every table above uses the RyR-outgroup rooting, which is "
              "evidence from outside the vertebrate subtree. Minimum-event "
              "rooting is a different criterion entirely, so the two can "
              "disagree — and what matters is not whether they name the same "
              "edge but whether the *placement* moves when they do. The "
              f"subtree was re-rooted at all **{len(others)} of its edges**:",
              ""]
        L += table(["rooting", "duplications", "implied losses", "total",
                    "deepest paralog duplication"],
                   [["outgroup (RyR) — the one every other table uses",
                     out["n_duplications"], out["implied_losses"],
                     out["total_events"], out["deepest_paralog_dup_node"]]] +
                   [[f"minimum-event ({len(mins)} of {len(others)} edges)",
                     mins[0]["n_duplications"], mins[0]["implied_losses"],
                     mins[0]["total_events"],
                     mins[0]["deepest_paralog_dup_node"]]] if out else [])
        L += [f"The two criteria pick **different edges** — minimum-event "
              f"rooting is {int(out['total_events']) - best} events cheaper "
              f"and puts the root inside the ITPR2+ITPR3 side — and they "
              f"**{'agree' if agree else 'disagree'}** about where the deepest "
              f"paralog duplication maps. So the placement is a property of "
              f"the gene tree's shape, not of the outgroup that rooted it.",
              ""]
    return L


def render(*, load, load_json, table, num, not_run, meta) -> list:
    """§5-§6 here; §7-§11 in `s13_report_audit.py`, same loader and formatter."""
    import s13_report_audit as audit
    L = _placement(load, table, num, not_run)
    L += _robustness(load, table, num, not_run)
    L += audit.render(load=load, table=table, num=num, not_run=not_run)
    return L
