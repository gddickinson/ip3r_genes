"""S7 — the two checks on the tree itself, rendered.

Split out of `s7_report_results.py` to stay inside the 500-line budget,
and split *here* because these two sections ask the same kind of
question — not "what does the tree say" but "how much should this
particular tree be trusted":

  § the `--bnni` check — does any claim rest on support that UFBoot
    granted because the model is violated rather than because the data
    say so
  § is the reported tree the best tree found — a constrained search can
    reach a better optimum than the free one, and when it does, the tree
    every table is built from is not the global optimum

Both are gated on their own tables and both render `not run yet` rather
than nothing, so a check that was skipped is visible as skipped. Like
`s7_report_naming.py`, this takes the caller's loader, table renderer
and formatter so the halves cannot read a table differently.
"""

from __future__ import annotations

# ------------------------------------------------- the --bnni robustness

def _section_bnni(load, load_json, table, num) -> list[str]:
    """Does every claim survive UFBoot's own guard against model violation?"""
    rows = load("bnni_comparison.tsv")
    L = ["### § The `--bnni` check on the supports", "",
         "UFBoot is optimistic when the model is violated, and a 134-tip "
         "alignment spanning four kingdoms violates any single "
         "substitution model by construction. `--bnni` is UFBoot's own "
         "guard: an extra round of NNI optimisation on every bootstrap "
         "tree, which strips support that came from the model rather "
         "than from the data. It is a **check, not a second answer** — "
         "the reported tree is still the main search's — so what is "
         "asked here is whether the claims that tree carries survive "
         "the guard, on the same clade sets, read from "
         "`claim_members.tsv`.", ""]
    b = load_json("tree_stats.json").get("bnni_check") or {}
    if b:
        L += [f"The check ran under the same model (`{b.get('model','')}`), "
              f"the same seed ({b.get('seed','')}) and the same pinned "
              f"thread count ({b.get('threads','')}), for "
              f"{num(float(b.get('runtime_s') or 0) / 60, 1)} min, at "
              f"log-likelihood {num(b.get('log_likelihood'), 3)} against "
              f"the main search's "
              f"{num((load_json('tree_stats.json').get('iqtree') or {}).get('log_likelihood'), 3)}.",
              ""]
    if not rows:
        return L + ["*The `--bnni` re-run has not been compared yet — "
                    "`python3 scripts/s7_run.py bnni` then "
                    "`python3 scripts/s7_bnni.py`.*", ""]
    tested = [r for r in rows if r["verdict"] in
              ("held", "weakened", "lost", "strengthened",
               "unsupported in both")]
    n = {v: sum(1 for r in tested if r["verdict"] == v)
         for v in ("held", "weakened", "lost", "strengthened",
                   "unsupported in both")}
    L += table(["claim", "tips", "main SH-aLRT/UFBoot",
                "`--bnni` SH-aLRT/UFBoot", "verdict"],
               [[r["claim"], r["n_tips"],
                 f"{r['main_alrt'] or '—'}/{r['main_ufboot'] or '—'}",
                 ("not a clade" if r["alt_is_clade"] == "no" else
                  f"{r['alt_alrt'] or '—'}/{r['alt_ufboot'] or '—'}"),
                 ("**" + r["verdict"] + "**"
                  if r["verdict"] in ("weakened", "lost")
                  else r["verdict"])] for r in rows])
    bad = n["weakened"] + n["lost"]
    # "No claim weakened" and "every claim holds" are different
    # statements, and conflating them flatters the tree: a claim that was
    # never well supported cannot be weakened, so it lands in
    # `unsupported in both` and would be counted as holding. Report the
    # two separately and name the unsupported ones with their numbers.
    unsup = [r for r in rows if r["verdict"] == "unsupported in both"]
    if not bad and tested and not unsup:
        L += ["", f"**Every one of the {len(tested)} claims tested holds "
              f"under `--bnni`** — same clades, both thresholds still "
              f"cleared. The supports this task reports are not an "
              f"artefact of the model.", ""]
    elif not bad and tested:
        L += ["", f"**No claim is weakened or lost by the guard**: "
              f"{n['held']} of {len(tested)} clear both thresholds in the "
              f"`--bnni` tree as well, so the supports this task leans on "
              f"are not an artefact of the model.", "",
              f"The remaining {len(unsup)} was **not well supported in "
              f"the main tree either**, so `--bnni` took nothing away — "
              f"it was already outside what this report treats as "
              f"trustworthy, and it is listed here rather than folded "
              f"into the count: "
              + "; ".join(
                  f"**{r['claim']}** ({r['main_alrt']}/{r['main_ufboot']} "
                  f"in the reported tree, {r['alt_alrt']}/{r['alt_ufboot']} "
                  f"under `--bnni`)" for r in unsup)
              + ". The guard does move it further down, which is the "
                "honest reading: that clade is the census's labelled set "
                "taken bare, and the tree prefers a slightly different "
                "grouping (§5.1).", ""]
    else:
        wk = [r["claim"] for r in tested if r["verdict"] == "weakened"]
        ls = [r["claim"] for r in tested if r["verdict"] == "lost"]
        L += ["", f"**{bad} of {len(tested)} claims do not survive the "
              f"guard.** "
              + (f"Weakened (present, but the support does not hold): "
                 f"{', '.join(wk)}. " if wk else "")
              + (f"Lost (the clade is not in the `--bnni` tree at all, so "
                 f"this is a topology change and not a support change): "
                 f"{', '.join(ls)}. " if ls else "")
              + "A claim in either list must be stated with its `--bnni` "
                "number, not the main run's.", ""]
    return L


# ------------------------------------------- is the reported tree the best?

def _section_besttree(load, load_json, table, num) -> list[str]:
    """Did a constrained search beat the unconstrained one?

    Written because it happened, and it is the kind of thing a pipeline
    quietly does not mention. A constrained search explores a different
    path through tree space, and on a 134-tip alignment it can land on a
    better optimum than the free search did. When it does, the tree every
    downstream table is built from is **not** the best tree found, and
    the honest thing is to say so and to say which claims it changes —
    not to swap the tree in silently, and not to leave it out.
    """
    au = load("au_test.tsv")
    if not au:
        return []
    try:
        best = min(au, key=lambda r: float(r["deltaL"]))
        ml = next(r for r in au if r["tree"] == "ML")
        gap = float(ml["deltaL"])
    except (ValueError, KeyError, StopIteration):
        return []
    if best["tree"] == "ML" or gap <= 0:
        return ["### § Is the reported tree the best tree found?", "",
                "Yes. No constrained search reached a higher likelihood "
                "than the unconstrained one, so the tree every table and "
                "figure here is built from is the best of the four "
                "searched.", ""]

    L = ["### § Is the reported tree the best tree found?", "",
         f"**No — and this is not a rounding difference.** The "
         f"{best['hypothesis']} search reached "
         f"**{num(best['logL'], 3)}**, {num(gap, 1)} log-likelihood "
         f"units above the unconstrained search's {num(ml['logL'], 3)}. "
         f"A constraint that costs *negative* likelihood is a statement "
         f"about the search, not about the hypothesis: constraining part "
         f"of the topology narrowed the space enough for the search to "
         f"find an optimum the free run missed. It is ordinary on a "
         f"134-tip alignment and it is why the AU test is run over a set "
         f"of trees rather than read off one.", "",
         "Two things follow, and both are stated rather than resolved by "
         "swapping trees. **The tree reported here is still the "
         "unconstrained search's** — every table, every figure and the "
         "rooting are built from it, and substituting a tree found under "
         "a constraint would make the topology partly an assumption. "
         "**But it is not the global optimum**, so the question is which "
         "claims differ between the two.", ""]
    cmp_rows = load("best_tree_comparison.tsv")
    if not cmp_rows:
        return L + ["*The claim-by-claim comparison has not been run — "
                    "`python3 scripts/s7_bnni.py --tree <best> --label "
                    "best --out best_tree_comparison.tsv`.*", ""]
    tested = [r for r in cmp_rows if r["verdict"] in ("present", "absent",
                                                      "held", "weakened",
                                                      "lost")]
    absent = [r for r in tested if r["verdict"] in ("absent", "lost")]
    if not absent:
        L += [f"**None of them.** All {len(tested)} claim clades are "
              f"present in the better tree as well, so the two trees "
              f"agree everywhere this report makes a statement and "
              f"differ only where it does not. The likelihood gap is "
              f"real and the conclusions do not rest on it.", ""]
    else:
        L += [f"**{len(absent)} of {len(tested)} claim clades are not in "
              f"the better tree**: "
              + "; ".join(f"{r['claim']} ({r['n_tips']} tips)"
                          for r in absent)
              + ". Any statement resting on those is conditional on the "
                "reported tree rather than on the alignment, and is "
                "flagged as such.", ""]
    L += table(["claim", "tips", "in the reported tree",
                "in the better tree"],
               [[r["claim"], r["n_tips"],
                 f"{r['main_alrt'] or '—'}/{r['main_ufboot'] or '—'}",
                 ("**no**" if r["verdict"] in ("absent", "lost")
                  else "yes")] for r in cmp_rows])
    return L
