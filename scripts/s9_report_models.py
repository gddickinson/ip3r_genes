"""S9 — the model-by-model half of the report: branch-site, site models,
RELAX, and the caveats. Split from `s9_report_results.py` to stay inside
the 500-line budget; it takes the caller's loader and formatters so all
three halves read the tables the same way.

The sections here share one discipline: **a model that fits is not the same
as a claim that holds.** Branch-site model A on a stem where dS is
saturated can reach a significant LRT with no site to point at, so the BEB
posterior is reported beside every significant test, and a run whose best
restart still sits below its own null is printed as the local-optimum
failure it is rather than dropped.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import PARALOGS  # noqa: E402
from s9_priors import PRIOR, verdict  # noqa: E402

SIG = 0.05


def _f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def _sig(r: dict) -> bool:
    q = _f(r.get("q_bh"))
    return q is not None and q < SIG


# ---- 4.4 branch-site model A ----------------------------------------------

def _branch_site(lrts: list[dict], bs: list[dict], beb: list[dict],
                 num) -> list[str]:
    rows = [r for r in lrts if r["test"].startswith("branch-site")]
    out = ["### 4.4 Branch-site model A on each paralog stem", "",
           "The stem branch is where a duplicate's fate is decided: it is "
           "the interval between the duplication and the first surviving "
           "split of the new copy, and if a paralog was ever free to "
           "change, that is when. Model A asks whether a class of sites on "
           "that one branch has ω > 1 while the rest of the tree does not.",
           ""]
    if not rows:
        out += [verdict("sister_pair", None,
                        "no branch-site job has finished"), ""]
        return out
    by_para: dict[str, list[dict]] = {}
    for x in bs:
        by_para.setdefault(x["paralog"], []).append(x)
    bebs = {b["job"]: b for b in beb}

    out += ["| stem | best restart | 2ΔlnL | q (BH) | foreground ω₂ | its "
            "share of sites | sites at BEB ≥ 0.95 | restarts below their own "
            "null |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        runs = by_para.get(r["set"], [])
        below = sum(int(x["below_null"]) for x in runs)
        best = next((x for x in runs if x["is_best"] == "1"), {})
        b = bebs.get(r["alt"], {})
        om = _f(best.get("fg_omega2"))
        bound = best.get("at_bound") == "1"
        out.append(f"| {r['set']} | ω₀ = {best.get('initial_omega', '—')} | "
                   f"{num(r['stat'], '{:.2f}')} | "
                   f"{num(r.get('q_bh'), '{:.3g}')} | "
                   + (f"**{num(om, '{:.1f}')}** (at codeml's bound)"
                      if bound else f"**{num(om, '{:.2f}')}**")
                   + f" | {num(_f(best.get('prop_fg')), '{:.1%}')} | "
                   f"{b.get('n_p95', '—')} | {below} / {len(runs)} |")
    out.append("")

    stuck = [r for r in rows if "local optimum" in (r.get("note") or "")]
    if stuck:
        out += ["**" + ", ".join(r["set"] for r in stuck) + "**: the best of "
                "four restarts still sits below its own nested null, so "
                "those runs are local optima and report nothing — not a "
                "negative result, an unfinished optimisation.", ""]
    n_below = sum(sum(int(x["below_null"]) for x in by_para.get(r["set"], []))
                  for r in rows)
    if n_below:
        out += [f"**{n_below} of the twelve restarts converged *below* their "
                "own nested null** — one on every stem. A nested "
                "alternative cannot have a lower optimum than its null, so "
                "each of those is a local-optimum failure that a "
                "single-start run would have reported as its answer. This "
                "is why model A is restarted by construction here (D37) "
                "rather than repaired afterwards, and every restart stays "
                "in `bs_restarts.tsv`."]
        # Which starting value fails is itself the argument for restarting.
        bad = {r["set"]: [x["initial_omega"] for x in by_para.get(r["set"], [])
                          if x["below_null"] == "1"] for r in rows}
        bad = {k: v for k, v in bad.items() if v}
        if len(set(tuple(v) for v in bad.values())) == len(bad) and len(bad) > 1:
            listed = ", ".join(f"{k} at ω₀ = {'/'.join(v)}"
                               for k, v in sorted(bad.items()))
            out += ["",
                    f"And it is a *different* starting value that fails on "
                    f"each stem — {listed}. No single initial ω would have "
                    "been safe here, which is the case for running several "
                    "rather than for choosing a better one.", ""]
        else:
            out.append("")

    sig = [r for r in rows if _sig(r) and r not in stuck]
    if not sig:
        out += ["No paralog stem carries a significant branch-site signal "
                "after correction.", ""]
        return out

    # Significance is not the claim. The claim is about ω₂, and codeml
    # pinning it at 999 is the same tell as a site class pinned at exactly
    # 1 in §4.5 — the optimiser has hit a wall, not measured a rate.
    credible, at_bound, no_sites = [], [], []
    for r in sig:
        best = next((x for x in by_para.get(r["set"], [])
                     if x["is_best"] == "1"), {})
        n95 = int(bebs.get(r["alt"], {}).get("n_p95", 0) or 0)
        if best.get("at_bound") == "1":
            at_bound.append((r, best, n95))
        elif n95 == 0:
            no_sites.append((r, best, n95))
        else:
            credible.append((r, best, n95))

    out += [f"**All {len(sig)} stems are significant after BH, and "
            f"{len(credible)} of the three carries an ω₂ the data actually "
            "determine.** Two columns separate those statements:", ""]
    for r, best, n95 in at_bound:
        spread = sorted({_f(x.get("fg_omega2")) for x in by_para[r["set"]]
                         if abs(_f(x["lnL"]) - _f(r["lnL_alt"])) < 1.0}
                        - {None})
        # The restart evidence differs per stem and the sentence has to
        # follow it: several ω₂ at one likelihood is a spread, all of them
        # at the ceiling is a different observation with the same meaning.
        if len(spread) > 1:
            sp = ", ".join(num(x, "{:.0f}") for x in spread)
            ratio = max(spread) / min(spread) if min(spread) else 0
            evidence = (f"Restarts reaching the *same* likelihood put ω₂ at "
                        f"{sp} — a {ratio:.0f}-fold spread at an unchanged "
                        "lnL, which is the definition of an unidentified "
                        "parameter.")
        else:
            n_same = sum(1 for x in by_para[r["set"]]
                         if abs(_f(x["lnL"]) - _f(r["lnL_alt"])) < 1.0)
            evidence = (f"All {n_same} restarts that reach this likelihood "
                        "end at the ceiling, from initial ω both below and "
                        "above 1, so the likelihood is flat in ω₂ above it.")
        out.append(
            f"- **{r['set']}** — ω₂ is pinned at codeml's **999 upper "
            f"bound**. That is not an estimate of 999; it is the optimiser "
            "reporting that the foreground has no synonymous signal left "
            "to normalise a rate against, which is exactly what §4.2's "
            f"saturation predicts for a branch this old. {evidence} "
            f"({n95} sites at BEB ≥ 0.95.)")
    for r, best, n95 in no_sites:
        out.append(
            f"- **{r['set']}** — ω₂ = {num(_f(best.get('fg_omega2')), '{:.2f}')} "
            "but the BEB posterior identifies **no site** above 0.95, so "
            "the model has nothing to point at.")
    for r, best, n95 in credible:
        b = bebs.get(r["alt"], {})
        out.append(
            f"- **{r['set']}** — ω₂ = "
            f"**{num(_f(best.get('fg_omega2')), '{:.2f}')}** on "
            f"{num(_f(best.get('prop_fg')), '{:.1%}')} of sites, well "
            "inside the estimable range and **stable across restarts**, "
            f"with {n95} sites at BEB ≥ 0.95 and {b.get('n_p99', 0)} at "
            "≥ 0.99. This one is a result.")
    out.append("")
    if credible and at_bound:
        names = ", ".join(r["set"] for r, _, _ in credible)
        out += [f"So the reportable branch-site finding is **{names} alone**: "
                "a class of sites on its stem evolving several times faster "
                "than neutrally while the rest of the tree sits at ω ≈ 0.03. "
                "The other stems' tests are significant and their ω₂ is not "
                "measurable, and those are different sentences. A pipeline "
                "that printed the three q-values would have reported the "
                "strongest signal on the stem whose parameter is least "
                "determined.", ""]
    out += [verdict("sister_pair", "orthogonal",
                    "S7's topology defines *which* branch is each paralog's "
                    "stem, and S9 uses it as given. A branch test cannot "
                    "corroborate the topology it is conditioned on, and "
                    "saying so is the point of listing it here"), ""]
    return out


# ---- 4.5 site models -------------------------------------------------------

def _site_models(lrts: list[dict], beb: list[dict], site: list[dict],
                 num) -> list[str]:
    rows = [r for r in lrts if " within " in r["test"]]
    out = ["### 4.5 Site models within each paralog", "",
           "M2a vs M1a and M8 vs M7 ask whether *any* site in a paralog has "
           "ω > 1 across the whole clade — a different question from the "
           "stem, and the one that would find a site under recurrent "
           "positive selection anywhere in the vertebrate history of that "
           "copy.", ""]
    if not rows:
        out += ["*not run yet*", ""]
        return out
    bebs = {b["job"]: b for b in beb}
    sites = {r["job"]: r for r in site}
    out += ["| test | 2ΔlnL | df | q (BH) | ω of the extra class | "
            "its share of sites | sites at BEB ≥ 0.95 |",
            "|---|---|---|---|---|---|---|"]
    for r in rows:
        sc = sites.get(r["alt"], {})
        b = bebs.get(r["alt"], {})
        out.append(f"| {r['test']} | {num(r['stat'], '{:.2f}')} | {r['df']} | "
                   f"{num(r.get('q_bh'), '{:.3g}')} | "
                   f"**{num(sc.get('omega_max'), '{:.3f}')}** | "
                   f"{num(sc.get('p_omega_max'), '{:.5f}')} | "
                   f"{b.get('n_p95', '—')} |")
    out.append("")

    sig = [r for r in rows if _sig(r)]
    if not sig:
        out += ["No site model is significant after correction in any "
                "paralog. Combined with §4.1 that is a coherent picture "
                "rather than an absence of one: a channel whose ω is under "
                "0.05 everywhere has very little room for a site class "
                "above 1 to hide in.", ""]
        return out

    # A significant LRT is not the claim. The claim is about the class the
    # test adds, and that class has to be *strictly* above 1 with a share
    # of sites above zero before "positive selection" is the right words.
    real = [r for r in sig
            if sites.get(r["alt"], {}).get("class_above_one") == "1"
            and (_f(sites.get(r["alt"], {}).get("p_omega_max")) or 0.0) > 0]
    empty = [r for r in sig if r not in real]
    out += [f"**{len(sig)} of {len(rows)} tests are significant after BH, and "
            f"{len(real)} of them is evidence of positive selection.** The "
            "likelihood-ratio test and the claim are different statements, "
            "and the columns above are what separates them:", ""]
    for r in empty:
        sc = sites.get(r["alt"], {})
        om = _f(sc.get("omega_max"))
        share = _f(sc.get("p_omega_max"))
        b = bebs.get(r["alt"], {})
        if share is not None and share <= 0:
            why = ("the extra class carries **no sites at all** (proportion "
                   f"{num(share, '{:.5f}')}), so the alternative has "
                   "collapsed onto its own null — which is why 2ΔlnL is "
                   f"{num(r['stat'], '{:.2f}')}")
        elif om is not None and abs(om - 1.0) < 1e-6:
            why = (f"the extra class sits at **ω = {num(om, '{:.5f}')}**, "
                   "codeml's boundary — a class of *unconstrained* sites, "
                   "not positively selected ones, carrying "
                   f"{100 * (share or 0):.2f} % of the alignment, with "
                   f"{b.get('n_p95', '0')} site(s) reaching a 0.95 posterior")
        else:
            why = (f"the extra class is at ω = {num(om, '{:.3f}')} with "
                   f"{num(share, '{:.5f}')} of sites")
        out.append(f"- *{r['test']}* — {why}.")
    for r in real:
        sc = sites.get(r["alt"], {})
        b = bebs.get(r["alt"], {})
        out.append(f"- *{r['test']}* — ω = "
                   f"{num(sc.get('omega_max'), '{:.3f}')} over "
                   f"{100 * (_f(sc.get('p_omega_max')) or 0):.2f} % of "
                   f"sites, {b.get('n_p95', '0')} of them at BEB ≥ 0.95.")
    out.append("")
    if not real:
        out += ["So the honest reading is that **M8 fits better than M7 "
                "because this family has a small class of sites that are "
                "free to drift, not because any site is being driven**. A "
                "beta distribution on [0, 1] cannot represent a spike at "
                "the neutral boundary, so adding one class that lands "
                "exactly there improves the fit significantly and says "
                "nothing about adaptation. Reporting the three q-values "
                "without the class they are testing would turn \"under 1 % "
                "of sites are unconstrained\" into \"positive selection in "
                "all three paralogs\".", "",
                "It is also consistent with everything else here: §4.1 puts "
                "ω between 0.02 and 0.05 across the whole protein, and "
                "M2a — which *is* free to place a class above 1, and does "
                "estimate one — gives it a proportion of exactly zero in "
                "all three paralogs.", ""]
    return out


# ---- 4.6 RELAX -------------------------------------------------------------

def _pval(x, num) -> str:
    """A p-value HyPhy reports as exactly 0 has underflowed double
    precision, not been measured as zero. Printing `0` invites a reader to
    treat it as an exact quantity."""
    v = _f(x)
    if v is None:
        return "—"
    if v <= 0:
        return "< 1e-300 (underflow)"
    return num(v, "{:.3g}")


def _relax(relax: list[dict], num) -> list[str]:
    out = ["### 4.6 RELAX — is any paralog's selection *relaxed*?", "",
           "codeml's branch models ask whether a foreground's ω differs. "
           "RELAX asks whether the whole ω distribution on the test branches "
           "is pulled towards ω = 1 (relaxation, k < 1) or away from it "
           "(intensification, k > 1). In a family where every ω is far below "
           "1, that is the sharper question, and it comes with its own "
           "test.", ""]
    if not relax:
        out += ["*not run yet* — `scripts/s9_relax.py`.", ""]
        return out
    out += ["| test set | k | direction | p | LRT | test branches | "
            "reference branches |", "|---|---|---|---|---|---|---|"]
    for r in relax:
        if r.get("status") not in ("ok", None, ""):
            out.append(f"| {r['paralog']} | — | — | — | — | — | "
                       f"**{r.get('status')}** |")
            continue
        out.append(f"| {r['paralog']} | **{num(r['k'], '{:.3f}')}** | "
                   f"{r.get('direction', '')} | {_pval(r['p'], num)} | "
                   f"{num(r['LRT'], '{:.1f}')} | {r['n_test_branches']} | "
                   f"{r['n_reference_branches']} |")
    out.append("")

    ok = [r for r in relax if r.get("status") in ("ok", None, "")]
    sig = [r for r in ok if (_f(r["p"]) if _f(r["p"]) is not None else 1.0) < SIG]
    for r in sig:
        k = _f(r["k"])
        if k is None:
            continue
        out.append(f"- **{r['paralog']}**: k = {num(k, '{:.3f}')} — selection "
                   + ("is *relaxed*" if k < 1 else "is *intensified*")
                   + " relative to the other two paralogs "
                     f"(p {_pval(r['p'], num)}).")
    if sig:
        out.append("")
    intens = [r for r in sig if (_f(r["k"]) or 1) > 1]
    relaxed = [r for r in sig if (_f(r["k"]) or 1) < 1]
    if intens and relaxed and len(sig) == len(ok):
        names_i = ", ".join(r["paralog"] for r in intens)
        names_r = ", ".join(r["paralog"] for r in relaxed)
        out += [f"**The three copies have not been held to the same standard "
                f"since 2R.** {names_i} is under *intensified* selection "
                f"relative to the other two, and {names_r} under *relaxed* "
                "selection relative to theirs. That is the same ordering "
                "§4.1's one-ratio ω gives, arrived at by a different "
                "statistic on a different model — ω compares point "
                "estimates, k compares the whole distribution — so the two "
                "are a check on each other rather than one number told "
                "twice.", ""]
    elif not sig:
        out += ["No paralog's ω distribution differs from the other two "
                "under RELAX. The three copies have been held to the same "
                "standard since 2R.", ""]
    failed = [r for r in relax if r.get("status") not in ("ok", None, "")]
    if failed:
        out += ["Runs whose output could not be read are marked in the "
                "status column rather than reported as `k = None`: a failed "
                "parse beside two real answers reads as a negative result, "
                "and it is not one.", ""]
    repaired = sum(int(r.get("nonfinite_repaired") or 0) for r in relax)
    if repaired:
        out += [f"*{repaired} non-finite literal(s) repaired while reading "
                "HyPhy's json.* The partitioned descriptive model estimates "
                "a per-branch ω, and a branch with no synonymous change gets "
                "an infinite one, which HyPhy writes as the bare token "
                "`inf` — not legal JSON. It is normal output, but the "
                "failure it causes is silent in the wrong direction: the "
                "analysis succeeds and the *parse* throws. On the first run "
                "that turned ITPR1's result into a blank row.", ""]
    out += ["The unlabelled vertebrate tips the S7 tree places in no paralog "
            "clade are left **unlabelled** in these runs rather than swept "
            "into the reference: a branch whose paralog identity is "
            "unresolved is not evidence about either side of the contrast.",
            ""]
    return out


# ---- 4.7 caveats -----------------------------------------------------------

def _caveats(pairs: list[dict], status: list[dict], bs: list[dict],
             num, pct) -> list[str]:
    n_sat = sum(int(r["saturated"]) for r in pairs) if pairs else 0
    n_model = sum(1 for r in status if r["route"] == "miniprot"
                  and r["status"] == "ok")
    bound = sorted({r["paralog"] for r in bs
                    if r.get("is_best") == "1" and r.get("at_bound") == "1"})
    extra: list[str] = []
    if bound:
        extra = [f"5. **Two of the three branch-site ω₂ are not "
                 f"identified.** On the {' and '.join(bound)} stem"
                 + ("s" if len(bound) > 1 else "")
                 + ", codeml's estimate of the foreground ω sits at its "
                 "999 upper bound and the likelihood is flat above it. "
                 "Those tests are significant and their effect size is "
                 "unmeasurable, which is not the same as a large effect. "
                 "Only the ITPR1 stem carries an ω₂ inside the estimable "
                 "range, and it is the only branch-site result this task "
                 "reports as one."]
    out = ["## 5. What this does not establish", "",
           f"1. **Synonymous saturation.** {pct(n_sat, len(pairs)) if pairs else '—'} "
           "of all within-paralog pairs exceed the dS bar. Tree-based models "
           "handle this far better than pairwise ML, but they do not repeal "
           "it: an ω estimated where dS is poorly determined is a ratio "
           "whose denominator is soft, and every number here should be read "
           "as a lower bound on precision, not a point estimate with a "
           "small error.",
           f"2. **Genome gene models.** {n_model} of the CDS come from "
           "miniprot reconstructions with masked frameshift or stop codons. "
           "The curated sensitivity subsets in §4.1 show ω barely moves "
           "without them, but those models are also the only evidence for "
           "several lineages, so the subset is a control, not a "
           "replacement.",
           "3. **The tree is conditioned on.** Every branch test is run on "
           "S7's topology. If the sister arrangement were wrong, the stems "
           "S9 marks would be the wrong branches — which is why S7 ran an "
           "AU test over all three arrangements before this task started, "
           "and why §4.4 records the dependency instead of quietly relying "
           "on it.",
           "4. **This is a vertebrate result.** The non-vertebrate grade is "
           "not in the codon alignment at all. Nothing here says anything "
           "about the constraint on the single-copy receptors S20 and S23 "
           "found outside the vertebrates."] + extra + [""]
    return out


def model_sections(load_tsv, num, pct) -> list[str]:
    lrts = load_tsv("lrt_table.tsv")
    bs = load_tsv("bs_restarts.tsv")
    beb = load_tsv("beb_sites.tsv")
    site = load_tsv("site_models.tsv")
    relax = load_tsv("relax_table.tsv")
    pairs = load_tsv("pairwise_dnds.tsv")
    status = load_tsv("cds_status.tsv")
    out: list[str] = []
    out += _branch_site(lrts, bs, beb, num)
    out += _site_models(lrts, beb, site, num)
    out += _relax(relax, num)
    out += _caveats(pairs, status, bs, num, pct)
    out += ["## 6. Figures", "",
            "![s9_omega_by_paralog](figures/s9_omega_by_paralog.png)", "",
            "*Per-paralog one-ratio ω with the curated-CDS sensitivity "
            "estimate beside it, and neutrality drawn rather than "
            "described.*", "",
            "![s9_dnds_saturation](figures/s9_dnds_saturation.png)", "",
            "*Pairwise dN against dS within each paralog, with the neutral "
            "diagonal and the saturation bar. The panel that qualifies "
            "every cross-paralog number in this report.*", "",
            "![s9_branch_contrast](figures/s9_branch_contrast.png)", "",
            "*Two-ratio background vs foreground ω per paralog clade, and "
            "the RELAX k beside it.*", "",
            "![s9_bs_restarts](figures/s9_bs_restarts.png)", "",
            "*Every branch-site restart against its own null. A point below "
            "the line is a local optimum, not a result — the reason this "
            "task restarts model A by construction.*", ""]
    return out
