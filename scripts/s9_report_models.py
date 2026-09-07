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
           "split of the new copy, and if a paralog was ever free to change, "
           "that is when. Model A asks whether a class of sites on that one "
           "branch has ω > 1 while the rest of the tree does not.", ""]
    if not rows:
        out += [verdict("sister_pair", None,
                        "no branch-site job has finished"), ""]
        return out
    by_para = {r["paralog"]: [x for x in bs if x["paralog"] == r["paralog"]]
               for r in bs}
    out += ["| stem | best restart | lnL alt | lnL null | 2ΔlnL | p (½χ²₁) | "
            "q (BH) | restarts below their own null |",
            "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        runs = by_para.get(r["set"], [])
        below = sum(int(x["below_null"]) for x in runs)
        best = next((x for x in runs if x["is_best"] == "1"), {})
        out.append(f"| {r['set']} | ω₀ = {best.get('initial_omega', '—')} | "
                   f"{num(r['lnL_alt'], '{:.2f}')} | "
                   f"{num(r['lnL_null'], '{:.2f}')} | "
                   f"{num(r['stat'], '{:.2f}')} | "
                   f"{num(r['p_reported'], '{:.3g}')} | "
                   f"{num(r.get('q_bh'), '{:.3g}')} | "
                   f"{below} / {len(runs)} |")
    out.append("")

    stuck = [r for r in rows if "local optimum" in (r.get("note") or "")]
    if stuck:
        out += ["**" + ", ".join(r["set"] for r in stuck) + "**: the best of "
                "four restarts still sits below its own nested null. A "
                "nested alternative cannot do that, so those runs are local "
                "optima and report nothing — not a negative result, an "
                "unfinished optimisation. They are printed rather than "
                "dropped so the reader can see which stems the search "
                "failed on.", ""]

    sig = [r for r in rows if _sig(r) and r not in stuck]
    bebs = {b["job"]: b for b in beb}
    if sig:
        out += ["Where the test is significant after BH, the BEB posterior "
                "is what says whether the model has sites to point at:", "",
                "| stem | sites with P(ω>1) ≥ 0.95 | ≥ 0.99 | median P | "
                "max P |", "|---|---|---|---|---|"]
        for r in sig:
            b = bebs.get(r["alt"], {})
            out.append(f"| {r['set']} | {b.get('n_p95', '—')} | "
                       f"{b.get('n_p99', '—')} | "
                       f"{num(b.get('median_p'), '{:.3f}')} | "
                       f"{num(b.get('max_p'), '{:.3f}')} |")
        out.append("")
        empty = [r["set"] for r in sig
                 if int(bebs.get(r["alt"], {}).get("n_p95", 0) or 0) == 0]
        if empty:
            out += ["A significant LRT with **no site above the 0.95 "
                    "posterior** is the signature of a branch-site fit "
                    "driven by saturation or alignment error rather than by "
                    "identifiable adaptive substitutions, and §4.2 has just "
                    "shown how much saturation this alignment carries. "
                    + ", ".join(empty) + " is reported that way.", ""]
    else:
        out += ["No paralog stem carries a significant branch-site signal "
                "after correction. On a family this constrained that is the "
                "expected answer, and it is worth saying plainly: the "
                "duplicates that made ITPR1/2/3 were not followed by a "
                "detectable episode of positive selection on the stems — "
                "at the resolution 2,459 codons and a saturated dS allow.",
                ""]
    out += [verdict("sister_pair", "orthogonal",
                    "S7's topology defines *which* branch is each paralog's "
                    "stem, and S9 uses it as given. A branch test cannot "
                    "corroborate the topology it is conditioned on, and "
                    "saying so is the point of listing it here"), ""]
    return out


# ---- 4.5 site models -------------------------------------------------------

def _site_models(lrts: list[dict], beb: list[dict], num) -> list[str]:
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
    out += ["| test | 2ΔlnL | df | p | q (BH) | sites P ≥ 0.95 |",
            "|---|---|---|---|---|---|"]
    bebs = {b["job"]: b for b in beb}
    for r in rows:
        b = bebs.get(r["alt"], {})
        out.append(f"| {r['test']} | {num(r['stat'], '{:.2f}')} | {r['df']} | "
                   f"{num(r['p_reported'], '{:.3g}')} | "
                   f"{num(r.get('q_bh'), '{:.3g}')} | "
                   f"{b.get('n_p95', '—')} |")
    out.append("")
    sig = [r for r in rows if _sig(r)]
    if not sig:
        out += ["No site model is significant after correction in any "
                "paralog. Combined with §4.1 that is a coherent picture "
                "rather than an absence of one: a channel whose ω is under "
                "0.05 everywhere has very little room for a site class "
                "above 1 to hide in.", ""]
    else:
        out += ["Significant: " + ", ".join(r["test"] for r in sig)
                + ". Every one is qualified by §4.2 — the same saturated "
                "synonymous sites that make the pairwise matrix a "
                "diagnostic also inflate a site model's ability to find a "
                "high-ω class, and the BEB column is what distinguishes a "
                "handful of identifiable sites from a flat posterior.", ""]
    return out


# ---- 4.6 RELAX -------------------------------------------------------------

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
    out += ["| test set | k | p | LRT | test branches | reference branches |",
            "|---|---|---|---|---|---|"]
    for r in relax:
        out.append(f"| {r['paralog']} | **{num(r['k'], '{:.3f}')}** | "
                   f"{num(r['p'], '{:.3g}')} | {num(r['LRT'], '{:.2f}')} | "
                   f"{r['n_test_branches']} | {r['n_reference_branches']} |")
    out.append("")
    sig = [r for r in relax if (_f(r["p"]) or 1.0) < SIG]
    for r in sig:
        k = _f(r["k"])
        if k is None:
            continue
        out.append(f"- **{r['paralog']}**: k = {num(k, '{:.3f}')}, "
                   + ("selection is *relaxed* relative to the other two "
                      "paralogs" if k < 1 else
                      "selection is *intensified* relative to the other two "
                      "paralogs")
                   + f" (p = {num(r['p'], '{:.3g}')}).")
    if sig:
        out.append("")
    else:
        out += ["No paralog's ω distribution differs from the other two "
                "under RELAX. The three copies have been held to the same "
                "standard since 2R.", ""]
    out += ["The unlabelled vertebrate tips the S7 tree places in no paralog "
            "clade are left **unlabelled** in these runs rather than swept "
            "into the reference: a branch whose paralog identity is "
            "unresolved is not evidence about either side of the contrast.",
            ""]
    return out


# ---- 4.7 caveats -----------------------------------------------------------

def _caveats(pairs: list[dict], status: list[dict], num, pct) -> list[str]:
    n_sat = sum(int(r["saturated"]) for r in pairs) if pairs else 0
    n_model = sum(1 for r in status if r["route"] == "miniprot"
                  and r["status"] == "ok")
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
           "found outside the vertebrates.", ""]
    return out


def model_sections(load_tsv, num, pct) -> list[str]:
    lrts = load_tsv("lrt_table.tsv")
    bs = load_tsv("bs_restarts.tsv")
    beb = load_tsv("beb_sites.tsv")
    relax = load_tsv("relax_table.tsv")
    pairs = load_tsv("pairwise_dnds.tsv")
    status = load_tsv("cds_status.tsv")
    out: list[str] = []
    out += _branch_site(lrts, bs, beb, num)
    out += _site_models(lrts, beb, num)
    out += _relax(relax, num)
    out += _caveats(pairs, status, num, pct)
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
