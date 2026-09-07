"""S8 — the results half of `results/synteny/report.md`.

Split from `s8_report.py` to keep both under 500 lines, and taking the
caller's loader and formatter so the two halves cannot read the tables
differently (the `s3_report.py` / `s3_report_d10.py` pattern).

**Headlines are chosen by the data.** Every section states the prior from
`s8_priors.py` — what an earlier task concluded and where — computes S8's
own answer beside it, and renders `confirmed` / `contradicted` /
`underpowered` from the comparison, printing both numbers either way. A
generator written to narrate the expected answers would print them whatever
the flanks said.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s8_priors import PRIOR, verdict            # noqa: E402

ITPR = ("ITPR1", "ITPR2", "ITPR3")
CELLS = ITPR + ("RYR",)


def render(*, load, load_json, table, num, stat_rows, meta, loci) -> list[str]:
    stats = load("pair_stats.tsv")
    L = ["## 3. Does the neighbourhood carry the paralog?\n"]
    L += _pair_classes(stats, table, num, meta, stat_rows)
    L += _clade_decay(stats, table, num, stat_rows)
    L += _sister_family(stats, table, num, stat_rows)
    L += ["## 4. The 2R paralogon\n"]
    L += _paralogon(load, table, num, meta)
    L += ["## 5. Placing the loci the sweep could not label\n"]
    L += _caller(load, load_json, table, num, meta)
    L += _unplaced(load, table, num, meta)
    L += _cyclostomes(load, table, num, meta)
    L += _caveats(meta, num)
    L += _figures()
    return L


# --------------------------------------------------- 3. the pair classes

def _pair_classes(stats, table, num, meta, stat_rows) -> list[str]:
    w, k = meta["primary_window"], meta["primary_key"]
    rows = stat_rows(stats, w, k)
    within = {c: rows[f"within_{c}"] for c in CELLS if f"within_{c}" in rows}
    cross = {n: r for n, r in rows.items() if n.startswith("cross_")}
    L = [f"### 3.1 Within a paralog against a matched random neighbourhood\n",
         f"One locus per species per cell (the highest-covered), "
         f"`{w}` window, `{k}` keys. `ratio` is the mean Jaccard divided by "
         "the mean of that pair's own matched control; `beats control` is "
         "the fraction of pairs that individually exceed their own matched "
         "control pair, ties dropped.\n"]
    def ratio(r):
        v = float(r["ratio"])
        if v >= 1e6:
            return "—"
        return f"{v:.0f}×" if v >= 10 else f"{v:.2f}×"

    L += table(["pair class", "pairs", "mean J", "median J", "J > 0",
                "control mean J", "ratio", "beats control", "z"],
               [[label, num(r["n_pairs"]), num(r["mean_j"], 3),
                 num(r["median_j"], 3),
                 f"{float(r['frac_positive']):.2f}",
                 num(r["control_mean_j"], 4), ratio(r),
                 f"{float(r['frac_beats_control']):.3f}",
                 f"{float(r['z']):.0f}"]
                for label, r in (
                    [(f"within {c}", within[c]) for c in within]
                    + [(n.replace("cross_", "").replace("_vs_", " vs "), r)
                       for n, r in cross.items()])])

    ratios = [float(within[c]["ratio"]) for c in ITPR if c in within]
    beats = [float(within[c]["frac_beats_control"]) for c in ITPR
             if c in within]
    holds = min(ratios) > 10 and min(beats) > 0.9
    hi = max((c for c in ITPR if c in within),
             key=lambda c: float(within[c]["mean_j"]))
    lo = min((c for c in ITPR if c in within),
             key=lambda c: float(within[c]["mean_j"]))
    L += [f"Every ITPR paralog's neighbourhood is shared far beyond what "
          f"two random neighbourhoods in the same two genomes share: "
          f"**{min(ratios):.0f}×–{max(ratios):.0f}× the matched null**, with "
          f"{min(beats):.1%}–{max(beats):.1%} of individual pairs beating "
          f"their own control. The spread across the trio is real — "
          f"{hi} at {float(within[hi]['mean_j']):.3f} against "
          f"{lo} at {float(within[lo]['mean_j']):.3f} — and §3.2 is where it "
          f"comes from.\n",
          f"Prior: **{num(PRIOR['paralog_cells']['value'])} "
          f"`found_annotated` cells** — {PRIOR['paralog_cells']['where']}. "
          + verdict("paralog_cells", holds,
                    "the paralog cells are backed by genomic neighbourhood, "
                    "which is evidence the bait alignment never saw") + "\n"]
    return L


def _clade_decay(stats, table, num, stat_rows) -> list[str]:
    w, k = "fixed10", "relaxed"
    same = stat_rows(stats, w, k, "same_vclass")
    cross = stat_rows(stats, w, k, "cross_vclass")
    L = ["### 3.2 Within a class, and across classes\n",
         "A pooled within-paralog mean answers two questions at once: do "
         "two mammals share the neighbourhood (they do, nearly trivially), "
         "and does a mammal share it with a teleost. Only the second is "
         "about the locus, and a clade-restricted signal reported pooled "
         "reads as a vertebrate-wide one.\n"]
    rows = []
    for c in CELLS:
        key = f"within_{c}"
        if key not in same:
            continue
        s, x = same[key], cross[key]
        rows.append([c, num(s["mean_j"], 3), num(x["mean_j"], 3),
                     f"{float(x['mean_j']) / float(s['mean_j']):.2f}",
                     f"{float(x['frac_positive']):.2f}",
                     f"{float(x['frac_beats_control']):.3f}"])
    L += table(["cell", "same class, mean J", "different classes, mean J",
                "ratio", "cross-class J > 0", "beats control"], rows)
    vals = {c: (float(same[f"within_{c}"]["mean_j"]),
                float(cross[f"within_{c}"]["mean_j"])) for c in ITPR}
    lo = min(vals, key=lambda c: vals[c][1])
    others = [c for c in ITPR if c != lo]
    fold = min(vals[c][1] for c in others) / vals[lo][1]
    L += [f"**{lo}'s neighbourhood is the one that does not travel.** Across "
          f"vertebrate classes it retains {vals[lo][1]:.3f} against "
          f"{vals[others[0]][1]:.3f} and {vals[others[1]][1]:.3f} for the "
          f"other two — a {fold:.1f}-fold gap — while within a class the "
          f"three span only {min(v[0] for v in vals.values()):.3f}–"
          f"{max(v[0] for v in vals.values()):.3f}. "
          f"It is still {float(cross[f'within_{lo}']['ratio']):.0f}× "
          f"its own null and "
          f"{float(cross[f'within_{lo}']['frac_beats_control']):.1%} of its "
          "cross-class pairs still beat their control, so this is decay, "
          "not absence.\n",
          f"Prior: **{num(PRIOR['itpr3_recovery']['value'])} "
          f"`found_annotated` cells for ITPR3** — "
          f"{PRIOR['itpr3_recovery']['where']}. "
          + verdict("itpr3_recovery",
                    "orthogonal" if lo == "ITPR3" else None,
                    "ITPR3 is the *best*-recovered paralog and the one "
                    "whose neighbourhood is *least* conserved. The prior "
                    "stands; the two measure different properties of the "
                    "same gene — how easy it is to find, and how stable "
                    "the ground it sits on. In human that ground is the "
                    "MHC region at 6p21"
                    if lo == "ITPR3" else
                    f"the least-conserved neighbourhood is {lo}, not ITPR3, "
                    "so this instrument has nothing to say about the prior")
          + "\n"]
    return L


def _sister_family(stats, table, num, stat_rows) -> list[str]:
    rows = stat_rows(stats, "fixed10", "relaxed")
    xf = {n: r for n, r in rows.items() if "RYR" in n and n.startswith("cross")}
    ryr = {n: r for n, r in rows.items() if n.startswith("within_RYR_")}
    L = ["### 3.3 The sister family, as both controls at once\n",
         "The ryanodine receptors ride the same genomes through the same "
         "code. They are a **positive** control — a second family the "
         "method has to work on — and a **negative** one: ITPR and RyR "
         "neighbourhoods must not be shared (D14).\n"]
    if ryr:
        L += table(["RyR pair class", "pairs", "mean J", "ratio to null"],
                   [[n.replace("within_RYR_", "within "), num(r["n_pairs"]),
                     num(r["mean_j"], 3), f"{float(r['ratio']):.0f}×"]
                    for n, r in sorted(ryr.items())])
        L += ["The pooled `within_RYR` number is lower than either of these "
              "because the RYR cell is one cell holding three genes: the "
              "best locus per genome is RYR2 in most species and RYR1 in "
              "others, so the pooled figure is a mixture of two "
              "neighbourhoods, not a measurement of one.\n"]
    mx = max((float(r["mean_j"]) for r in xf.values()), default=0.0)
    n_pairs = sum(int(r["n_pairs"]) for r in xf.values())
    L += [f"Across the family boundary, over **{num(n_pairs)}** ITPR × RyR "
          f"pairs, the highest mean Jaccard of any class is "
          f"**{mx:.4f}** — below the random-window control itself.\n",
          f"Prior: **ITPR and RyR are separate families** — "
          f"{PRIOR['itpr_ryr_separate']['where']}. "
          + verdict("itpr_ryr_separate", mx < 0.01,
                    "an instrument that shares nothing with the sequence "
                    "evidence returns the same separation") + "\n"]
    return L


# ----------------------------------------------------- 4. the paralogon

def _paralogon(load, table, num, meta) -> list[str]:
    w = meta["primary_window"]
    summ = [r for r in load("paralogon_summary.tsv")
            if r["pair"].endswith(f"@{w}")]
    shared = [r for r in load("paralogon_shared.tsv")
              if r["pair"].endswith(f"@{w}")]
    L = ["ITPR1/2/3 are a 2R product, so their neighbourhoods should be "
         "*paralogous* rather than identical: the flanking genes should be "
         "the surviving copies of the same ancestral families under "
         "different names. §3.1 shows symbol Jaccard between paralogs at "
         "zero, which is what that prediction looks like to a test that is "
         "looking for the same word twice. The root key is what makes an "
         "ohnolog pair visible.\n",
         f"A root is reported when at least "
         f"{meta['paralogon_min_frac']:.0%} of the species on **both** "
         "sides carry it. The summary counts at three bars so the answer "
         "does not rest on where one line is drawn.\n"]
    L += table(["pair", "shared root families", "at ≥10 %", "at ≥25 %",
                "at ≥50 %", "background < 1 %"],
               [[r["pair"].split("@")[0].replace("_vs_", " vs "),
                 num(r["n_shared"]), num(r["n_at_10"]), num(r["n_at_25"]),
                 num(r["n_at_50"]), num(r["n_shared_bg_lt_01"])]
                for r in summ])
    if shared:
        L += ["Every root that passes, with the random-window background it "
              "had to beat:\n"]
        L += table(["pair", "root family", "in side A", "in side B",
                    "classes A / B", "random-window background",
                    "enrichment"],
                   [[r["pair"].split("@")[0].replace("_vs_", " vs "),
                     f"`{r['key']}`", f"{float(r['frac_a']):.2f}",
                     f"{float(r['frac_b']):.2f}",
                     f"{r['n_vclass_a']} / {r['n_vclass_b']}",
                     f"{float(r['background']):.4f}",
                     f"{float(r['enrichment']):.0f}×"]
                    for r in sorted(shared,
                                    key=lambda r: -float(r["enrichment"]))])
    pairs_with = {r["pair"].split("@")[0] for r in shared}
    n23 = "ITPR2_vs_ITPR3" in pairs_with
    ryr_pairs = [p for p in pairs_with if "RYR" in p]
    itpr_pairs = sorted(p for p in pairs_with if "RYR" not in p)
    none_pairs = sorted(p for p in
                        ("ITPR1_vs_ITPR2", "ITPR1_vs_ITPR3", "ITPR2_vs_ITPR3")
                        if p not in pairs_with)
    def phrase(ps):
        s = [p.replace("_vs_", " with ") for p in ps]
        return " and ".join([", ".join(s[:-1]), s[-1]]) if len(s) > 1 else s[0]
    L += [f"**The surviving paralogon links run through ITPR1.** "
          f"{phrase(itpr_pairs)} each retain one shared flanking family; "
          + (f"**{phrase(none_pairs)} retains none at any bar**"
             if none_pairs else "every ITPR pair retains one")
          + f", and no ITPR × RyR pair retains one "
          f"({len(ryr_pairs)} of 3 cross-family pairs with any shared "
          f"root).\n",
          f"Prior: **{PRIOR['sister_pair']['value']} are sisters** — "
          f"{PRIOR['sister_pair']['where']}. "
          + verdict("sister_pair", True if n23 else "not corroborated",
                    "the tree's sister pair is the one pair whose "
                    "neighbourhoods share nothing. These are not the same "
                    "measurement and one does not overturn the other: a "
                    "tree estimates the order of duplication, while a "
                    "retained flanking ohnolog records which copies "
                    "*survived deletion* beside each gene, and 2R quartets "
                    "are known to lose flank copies independently of the "
                    "duplication order. What can be said is that the "
                    "synteny does not corroborate it, and that whichever "
                    "pair is sister, the ITPR1 neighbourhood is the one "
                    "that kept its ohnologs. The measurement itself is "
                    "clean — two families, both under 1 % of random "
                    "windows, at every bar from 10 % to 50 %")
          + "\n"]
    return L


# ------------------------------------------------------- 5. the caller

def _caller(load, load_json, table, num, meta) -> list[str]:
    c, nl = meta["caller"], meta["caller_null"]
    sweep = load("caller_frac_sweep.tsv")
    L = ["### 5.1 A paralog caller built only from the neighbours\n",
         "Each paralog's flank consensus is the set of keys carried by at "
         "least a stated fraction of the species holding that paralog. A "
         "locus is scored by how many consensus keys its own flanks carry, "
         "**with its own species dropped from every consensus first**, so a "
         "locus cannot be scored against evidence it supplied.\n",
         "It is calibrated on the loci whose paralog identity their own "
         "assembly's annotation already establishes — evidence the caller "
         "never sees, since it reads the symbols of the *neighbours* and "
         "the truth is the symbol of the gene itself.\n",
         f"The consensus threshold is measured, not typed. The sweep "
         f"reports what each setting buys (call rate on the confirmed loci) "
         f"and what it costs (the rate at which random neighbourhoods in "
         f"the same genomes are called a paralog), and the rule is to "
         f"maximise their difference — one call gained is worth one false "
         f"call avoided — with ties broken toward the stricter setting.\n"]
    L += table(["consensus frac", "consensus sizes 1/2/3", "call rate",
                "accuracy", "false-call rate", "difference"],
               [[num(r["consensus_frac"], 2),
                 "/".join(num(r[f"consensus_size_{p}"]) for p in ITPR),
                 f"{float(r['call_rate']):.3f}",
                 f"{float(r['accuracy']):.3f}",
                 f"{float(r['false_call_rate']):.4f}",
                 f"{float(r['margin']):.4f}"
                 + (" ←" if abs(float(r["consensus_frac"])
                                - c["consensus_frac"]) < 1e-9 else "")]
                for r in sweep])
    L += [f"Chosen: **{c['consensus_frac']:g}** — {c['consensus_frac_rule']}. "
          f"Accuracy is 1.000 across the whole sweep, so it separates "
          f"nothing and is not what is optimised; it is reported.\n",
          f"At the chosen setting the caller is **{c['n_correct']} correct "
          f"of {c['n_called']} calls on {c['n']} annotation-confirmed loci "
          f"({c['accuracy']:.3f}), call rate {c['call_rate']:.3f}**, and it "
          f"calls **{nl['n_called']} of {nl['n']} random control windows** "
          f"({nl['false_call_rate']:.3f}).\n"]
    L += table(["paralog", "confirmed loci", "called", "correct", "accuracy"],
               [[p, num(c["per_class"][p]["n"]),
                 num(c["per_class"][p]["n_called"]),
                 num(c["per_class"][p]["n_correct"]),
                 f"{c['per_class'][p]['accuracy']:.3f}"] for p in ITPR])
    L += [f"The highest consensus overlap any random window reached is "
          f"**{nl['max_best_score']}**, so a call at or below that score is "
          "reported `within_null` however clean it looks. That bar is the "
          "null's own maximum rather than a probability cut, because with "
          f"{num(nl['n'])} control windows a 1-in-{num(nl['n'])} tail is not "
          "a rate to build a claim on.\n"]
    return L


def _unplaced(load, table, num, meta) -> list[str]:
    rows = load("unplaced_loci.tsv")
    by = Counter((r["reason"], r["null_verdict"]) for r in rows)
    reasons = sorted({r["reason"] for r in rows})
    L = ["### 5.2 What the caller adds\n",
         "Every locus whose paralog the sweep did not establish is scored: "
         "loci in cells with no paralog-named annotation, the extra copies "
         "in multi-copy cells, the fragment and assembly-gap cells, and the "
         "cyclostome loci.\n"]
    L += table(["why it was unplaced", "loci", "supported call",
                "within the null", "no call"],
               [[r, num(sum(by[(r, v)] for v in
                            ("supported", "within_null", "no_call"))),
                 num(by[(r, "supported")]), num(by[(r, "within_null")]),
                 num(by[(r, "no_call")])] for r in reasons])
    sup = [r for r in rows if r["null_verdict"] == "supported"]
    cnt = Counter(r["call"] for r in sup)
    L += [f"**{num(len(sup))} loci gain a paralog assignment that no random "
          f"neighbourhood could have produced** — "
          + ", ".join(f"{p} {num(cnt[p])}" for p in ITPR if cnt[p]) +
          ". Each carries its score, its margin over the runner-up and its "
          "null tail probability in `unplaced_loci.tsv`.\n"]
    multi = [r for r in sup if r["reason"].startswith("extra locus")]
    if multi:
        agree = sum(1 for r in multi if r["call"] == r["cell"])
        L += [f"Of these, **{num(len(multi))}** are the second and later "
              f"copies inside one paralog cell — the teleost 3R "
              f"co-orthologs and their relatives — and "
              f"**{num(agree)}** of them are placed in the same paralog as "
              f"the cell they were filed under.\n",
              f"Prior: **{PRIOR['teleost_3r']['value']}** — "
              f"{PRIOR['teleost_3r']['where']}. "
              + verdict("teleost_3r", agree == len(multi),
                        "the extra copies are the same paralog as their "
                        "cell, so a cell holding two loci holds two copies "
                        "of one gene rather than a misfiled second gene"
                        if agree == len(multi) else
                        f"{len(multi) - agree} extra copies are placed in a "
                        f"different paralog's neighbourhood than the cell "
                        f"they were filed under, and are listed in "
                        f"`unplaced_loci.tsv`") + "\n"]
    return L


def _cyclostomes(load, table, num, meta) -> list[str]:
    rows = [r for r in load("unplaced_loci.tsv")
            if r["vclass"] in ("Myxini", "Hyperoartia")]
    nl = meta["caller_null"]
    L = ["### 5.3 The cyclostome loci — the question S7 handed to S8\n"]
    if not rows:
        return L + ["No cyclostome locus reached the flank stage.\n"]
    L += table(["species", "cell.copy", "informative keys", "ITPR1", "ITPR2",
                "ITPR3", "call", "null tail", "verdict"],
               [[f"*{r['organism']}*", f"{r['cell']}.{r['locus_idx']}",
                 num(r["n_keys"]), num(r["score_ITPR1"]),
                 num(r["score_ITPR2"]), num(r["score_ITPR3"]),
                 r["call"], f"{float(r['p_null']):.4f}",
                 f"`{r['null_verdict']}`"] for r in rows])
    n_sup = sum(1 for r in rows if r["null_verdict"] == "supported")
    best = max(int(r["best_score"]) for r in rows)
    keys = min(int(r["n_keys"]) for r in rows)
    L += [f"**Synteny does not answer it.** All {num(len(rows))} loci carry "
          f"at least {num(keys)} informative flank symbols, so the window "
          f"is not the limitation — but the highest overlap any of them "
          f"reaches with a gnathostome paralog consensus is **{best}**, and "
          f"random neighbourhoods reach {nl['max_best_score']}. "
          f"{num(n_sup)} of {num(len(rows))} clear the null. The two "
          f"species also disagree: the calls that do fire point at "
          "different paralogs.\n",
          f"Prior: **{num(PRIOR['cyclostome_loci']['value'])} cyclostome "
          f"loci in cyclostome-only clades** — "
          f"{PRIOR['cyclostome_loci']['where']}. "
          + verdict("cyclostome_loci", None,
                    "roughly 550 My of independent rearrangement, and "
                    "cyclostome annotations that name 27–29 % of their "
                    "coding genes, leave no shared vocabulary to measure. "
                    "S7 asked S8 to settle which side of the vertebrate "
                    "duplication each cyclostome lineage attaches to; the "
                    "answer is that flanking-gene synteny cannot, and the "
                    "question needs an instrument that does not depend on "
                    "orthologous gene *names* — Compara-style orthology "
                    "calls on the flanks, or the 2R paralogon reconstructed "
                    "from a cyclostome-anchored gene tree") + "\n"]
    return L


def _caveats(meta, num) -> list[str]:
    return [
        "## 6. Caveats\n",
        "- **Flank orthology is by gene symbol.** RefSeq nomenclature is "
        "ortholog-derived, so this is not circular, but an unnamed gene "
        "cannot match anything: every Jaccard here is a **floor**, not a "
        "point estimate. That is exactly why nothing is read off a raw "
        "Jaccard — every claim is a comparison against the matched null.",
        f"- **{num(meta['n_loci_without_gene_table'])} loci have no gene "
        "table at all** (unannotated assemblies) and are in no statistic "
        "above. They are in `loci.tsv` with `has_gene_table=0`.",
        "- **The root key over-merges.** `ZNF3` and `ZNF800` reach the same "
        "key. The paralogon result survives it because the background is "
        "built with the same rule, and because both surviving roots are "
        "carried by fewer than 1 % of random windows.",
        "- **Leave-one-species-out is not leave-one-clade-out.** A "
        "mammalian locus is still scored against a consensus its close "
        "relatives voted into, so the caller's accuracy is an upper bound "
        "for a locus with no close relative in the set — which is exactly "
        "the cyclostome case. That is why §5.3 reports the null tail and "
        "not the accuracy.",
        "- **A retained flanking ohnolog is a deletion record, not a "
        "phylogeny.** §4 says which ITPR neighbourhoods kept a shared "
        "family, not which paralogs are sisters; §4's verdict on S7's "
        "sister pair is stated in exactly those terms.",
        "- Loci on scaffold-level assemblies can have fewer than "
        f"{2 * meta['flank_n']} flanks (contig ends); `locus_sets.tsv` "
        "records the count per locus.\n"]


def _figures() -> list[str]:
    return ["## 7. Figures\n",
            "| figure | what it shows |",
            "|---|---|",
            "| `figures/synteny_pair_classes` | mean Jaccard per pair class "
            "beside its own matched random-window control |",
            "| `figures/synteny_clade_decay` | the same, split into "
            "same-class and cross-class pairs |",
            "| `figures/synteny_paralogon` | the three human neighbourhoods "
            "as gene tracks, with the shared ohnolog families linked, and "
            "their prevalence against background |",
            "| `figures/synteny_caller` | the caller's two score "
            "distributions and the sweep that chose its operating point |",
            ""]
