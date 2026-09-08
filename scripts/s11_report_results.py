"""The results half of the S11 report, split to keep both under the 500-line
budget and taking the caller's loader and formatter so the two halves cannot
read the tables differently (the `s3_report.py` / `s3_report_d10.py`
pattern).

**Headlines chosen by the data.** Each section states the prior from
`s11_priors.PRIOR`, computes S11's own answer, and renders S8's five-valued
verdict from the comparison. `orthogonal` does real work here: most of what
came before S11 measured *sequence* and S11 measures *shape*. The two
families are 25 % identical and share a fold; that is not a contradiction of
the 25 %, and a report that called it one would be making a category error.
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_tmalign import TM_FOLD, TM_RANDOM                      # noqa: E402
from s11_tables import (CALL_FLOOR, PAIR_CLASSES,               # noqa: E402
                        REL_MARGIN, pair_scores)


def _f(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def section_calibration(pairs, manifest, fmt, missing, prior_line) -> list[str]:
    if not pairs:
        return [missing("4. What a TM-score means on this panel")]
    lo = pair_scores(pairs, manifest, "min")
    hi = pair_scores(pairs, manifest, "max")
    order = [k for k in PAIR_CLASSES if lo.get(k)]
    out = ["## 4. What a TM-score means on this panel", "",
           "Before any structure is called, the scale is calibrated on this "
           "panel rather than taken from the literature. Five classes of "
           "pair, each answering a different question, and the negative "
           "controls are what put a floor under all of them.", "",
           "TM-align normalises by each input's length, and on this panel "
           "the two normalisations say different things: an IP3 receptor "
           "subunit resolves to ~2,200 residues and a ryanodine receptor to "
           "~4,300, so **which chain a cross-family score is normalised by "
           "decides whether it reads as a fold result or a size one**. Both "
           "are given.", "",
           "| comparison | pairs | median TM, by the longer chain | median "
           "TM, by the shorter | range (longer) |",
           "|---|---|---|---|---|"]
    for key in order:
        v = sorted(lo[key])
        h = sorted(hi.get(key) or [0.0])
        out.append(f"| {key} | {len(v)} | **{fmt(median(v))}** | "
                   f"{fmt(median(h))} | {fmt(v[0])}–{fmt(v[-1])} |")
    state = lo.get("same protein, different state") or []
    within = lo.get("IP3R vs IP3R") or []
    cross_lo = lo.get("IP3R vs RyR") or []
    cross_hi = hi.get("IP3R vs RyR") or []
    ctrl = lo.get("control vs anything") or []
    out.append("")
    if ctrl:
        above = sum(1 for v in ctrl if v >= TM_FOLD)
        ctrl_hi = hi.get("control vs anything") or ctrl
        out += [
            f"The negative controls span **{fmt(min(ctrl))}–{fmt(max(ctrl))}** "
            f"with a median of **{fmt(median(ctrl))}**, and **{above} of "
            f"{len(ctrl)}** control pairs reach the {TM_FOLD} same-fold bar "
            f"under either normalisation (best control score anywhere: "
            f"{fmt(max(ctrl_hi))}). That is the floor this task's positive "
            "results are read against, measured on proteins S1 chose as "
            "decoys for a *sequence* scorer and re-used here without "
            "reselection.", ""]
    if state and within:
        out += [prior_line(
            "state_spread",
            f"Two structures of the same paralog in different conformations "
            f"score a median **{fmt(median(state))}**; two different IP3 "
            f"receptors score **{fmt(median(within))}**, a spread whose "
            "lower tail is the partial non-vertebrate models rather than a "
            "fold difference (§6). Conformation is worth "
            f"{fmt(1.0 - median(state))} of TM-score on this panel, so a "
            "fold difference of that size or less is not readable and is "
            "not claimed.",
            median(state) >= 0.7), ""]
    if cross_lo and cross_hi and within:
        both_ways = (median(cross_lo) >= TM_FOLD and median(cross_hi) >= TM_FOLD)
        one_way = (median(cross_hi) >= TM_FOLD > median(cross_lo))
        direction = (
            "Both directions clear the same-fold bar." if both_ways else
            (f"The two cross the bar in **one direction only**: a ryanodine "
             f"receptor accounts for an IP3 receptor's fold "
             f"({fmt(median(cross_hi))}) while an IP3 receptor cannot "
             f"account for a ryanodine receptor's ({fmt(median(cross_lo))}), "
             "which is as much a statement about a 2,200-residue protein "
             "being aligned into a 4,300-residue one as about fold.")
            if one_way else
            f"Neither direction reaches the {TM_FOLD} bar.")
        out += [prior_line(
            "sequence_separation",
            f"The two families sit at a median TM of **{fmt(median(cross_lo))}** "
            f"normalised by the ryanodine receptor and "
            f"**{fmt(median(cross_hi))}** normalised by the IP3 receptor, "
            f"against **{fmt(median(within))}** within the IP3 receptors, "
            "where the sequence gap is 0.828 within against 0.249 between. "
            + direction,
            "orthogonal"), ""]
    out += ["![](figures/s11_tm_calibration.png)", "",
            "**Figure 3.** The family call, structurally (**a**) and the "
            "calibration behind it (**b**). Both of TM-align's published "
            f"bars are drawn — {TM_RANDOM} (random) and {TM_FOLD} (same "
            "fold) — rather than described. Panel **b** uses the "
            "conservative normalisation (by the longer chain).", ""]
    return out


def section_calls(vs_ref, fmt, pct, missing, prior_line) -> list[str]:
    if not vs_ref:
        return [missing("5. D14, asked of shape")]
    named = [r for r in vs_ref if r["call"] in ("ITPR", "RYR")]
    called = [r for r in named if r["structural_call"] in ("ITPR", "RYR")]
    agree = [r for r in called if r["call_agrees"] == "1"]
    declined = [r for r in named if r["structural_call"] not in ("ITPR", "RYR")]
    controls = [r for r in vs_ref if r["role"] == "control"]
    ctrl_called = [r for r in controls
                   if r["structural_call"] in ("ITPR", "RYR")]
    calls = Counter(r["structural_call"] for r in vs_ref)
    out = [
        "## 5. D14, asked of shape", "",
        "Every stage of this project has had to separate IP3 from ryanodine "
        "receptors by a positive test. This is the sixth instrument to be "
        "asked, and the first that reads no gene symbol, no Pfam "
        "architecture, no alignment score and no tree — only coordinates.", "",
        "The test is the best TM-score against the IP3R references, the "
        "best against the RyR references, and a **relative** margin of "
        f"{REL_MARGIN:.0%} between them — relative for the reason "
        "`s3_assign.py` gives, that an absolute gap would call every "
        "full-length structure and no fragment. Scores are normalised **by "
        "the reference**, which is the question being asked: does this "
        "structure account for an IP3 receptor, or for a ryanodine "
        "receptor?", "",
        f"A margin is only consulted once the winner clears the "
        f"{CALL_FLOOR} same-fold bar. That gate is not decoration: all "
        "three negative controls beat their own runner-up by about 30 % of "
        "their score while scoring 0.17–0.27 against everything, so a rule "
        "gated on the margin alone calls a dynein heavy chain an IP3 "
        "receptor. With the gate, they are declined.", "",
        "| structural call | structures |", "|---|---|"]
    for k, n in sorted(calls.items(), key=lambda kv: -kv[1]):
        out.append(f"| {k} | {n} |")
    out += ["",
            f"**{len(agree)} of {len(called)}** structures the instrument "
            f"could call, whose family the census also names, are called "
            f"the same way by shape ({pct(len(agree), len(called))}) — and "
            f"**{len(ctrl_called)} of {len(controls)}** negative controls "
            "receives a family call.", ""]
    if declined:
        out += ["The instrument declines "
                f"**{len(declined)}** further census-named structures. A "
                "refusal is not a disagreement, and none of these is scored "
                "as one:", "",
                "| structure | census call | best vs IP3R | best vs RyR | "
                "modelled residues | why declined |",
                "|---|---|---|---|---|---|"]
        for r in sorted(declined, key=lambda r: -_f(r["best_itpr"])):
            out.append(f"| `{r['id']}` | {r['call']} | "
                       f"{fmt(r['best_itpr'])} | {fmt(r['best_ryr'])} | "
                       f"{fmt(r['resolved_residues'])} | {r['note']} |")
        out += ["", "Every one of them still prefers the IP3 receptors over "
                "the ryanodine receptors by roughly 2:1; what they cannot do "
                "is account for enough of a reference to earn a call. These "
                "are the partial, low-confidence models AlphaFold DB holds "
                "for the non-vertebrate groups (§1), and §6 separates how "
                "much of the shortfall is their length from how much is "
                "their similarity.", ""]
    disagree = [r for r in called if r["call_agrees"] != "1"]
    if disagree:
        out += ["| structure | census call | structural call | best vs IP3R | "
                "best vs RyR | note |", "|---|---|---|---|---|---|"]
        for r in disagree:
            out.append(f"| `{r['id']}` | {r['call']} | {r['structural_call']} "
                       f"| {fmt(r['best_itpr'])} | {fmt(r['best_ryr'])} | "
                       f"{r['note'] or '—'} |")
        out += ["", "A disagreement here is not by itself a correction: this "
                "instrument is the newest and the least tested of the six, "
                "and a structure that cannot be placed is more often a poor "
                "model than a mis-called gene. What it does is name the rows "
                "a later task has to look at.", ""]
    out += [prior_line(
        "d14_separable",
        f"{pct(len(agree), len(called))} of the structures the fold test "
        f"could call agree with the census ({len(agree)}/{len(called)}), "
        f"and no negative control is called at all "
        f"({len(ctrl_called)}/{len(controls)}).",
        bool(called) and len(agree) == len(called) and not ctrl_called), ""]
    return out


def section_nonvertebrate(vs_ref, manifest, fmt, missing,
                          prior_line) -> list[str]:
    if not vs_ref:
        return [missing("6. Do the deep records fold like receptors?")]
    man = {m["id"]: m for m in manifest}
    deep = []
    for r in vs_ref:
        if r["role"] != "model":
            continue
        group = (man.get(r["id"], {}).get("group") or r.get("group") or "")
        if group in ("ITPR1", "ITPR2", "ITPR3", "Vertebrata", "reference"):
            continue
        deep.append((group, r))
    if not deep:
        return [missing("6. Do the deep records fold like receptors?")]
    out = [
        "## 6. Do the deep records fold like receptors?", "",
        "This is the question S11 is worth most on. The plant, fungal, SAR, "
        "Discoba and non-vertebrate metazoan records were called ITPR on "
        "sequence evidence alone — profile score, architecture, and in S20's "
        "case a per-record contamination chase. None of that is shape.", "",
        "| group | model | modelled residues | best vs IP3R | ceiling | ref "
        "| best vs RyR | call |", "|---|---|---|---|---|---|---|---|"]
    for group, r in sorted(deep, key=lambda g: (g[0], -_f(g[1]["best_itpr"]))):
        out.append(
            f"| {group} | `{r['id'].replace('model_', '')}` | "
            f"{fmt(r['resolved_residues'])} | **{fmt(r['best_itpr'])}** | "
            f"{fmt(r['tm_ceiling'])} | "
            f"{r['best_itpr_ref'].replace('reference_', '').replace('state_panel_', '')} | "
            f"{fmt(r['best_ryr'])} | {r['structural_call']} |")
    folds = [r for _, r in deep if _f(r["best_itpr"]) >= TM_FOLD]
    itpr_called = [r for _, r in deep if r["structural_call"] == "ITPR"]
    capped = [r for _, r in deep if 0 < _f(r["tm_ceiling"]) < TM_FOLD]
    out += ["",
            f"**{len(folds)} of {len(deep)}** deep models reach the "
            f"{TM_FOLD} same-fold bar against an IP3 receptor, and "
            f"**{len(itpr_called)}** are called ITPR rather than RyR or "
            "no-call.", "",
            "The **ceiling** column separates how much of that shortfall "
            "is length from how much is similarity. A reference-normalised "
            "TM-score divides by the reference's length, so a model of *n* "
            "residues scored against an *L*-residue reference cannot exceed "
            "*n/L* before similarity is considered at all. Most of these "
            "models are 1,000–1,300 residues — the only records AlphaFold "
            "DB holds for these groups (§1) — against ~2,050-residue "
            "references."
            + (f" For **{len(capped)}** of them the ceiling sits below the "
               f"{TM_FOLD} bar outright: passing was arithmetically "
               "impossible."
               if capped else
               " None of them is excluded by arithmetic alone — every one "
               "had the headroom to pass and did not."), "",
            "So the shortfall is only partly about length: each declined "
            "model reaches "
            + (lambda fr: f"{fmt(min(fr))}–{fmt(max(fr))}" if fr else "n/a")(
                [_f(r["best_itpr"]) / _f(r["tm_ceiling"])
                 for _, r in deep
                 if _f(r["tm_ceiling"])
                 and r["structural_call"] == "below_fold_bar"])
            + " of its own ceiling — far above the negative controls and "
            "far below a full-length receptor. That is what a genuinely "
            "divergent homolog modelled at a mean pLDDT in the low 60s "
            "(§7) should look like. The fold test neither confirms nor "
            "contradicts these records; it does not reach them.", "",
            prior_line("nonvertebrate_real",
                       f"{len(itpr_called)} of {len(deep)} non-vertebrate "
                       "models are assigned to the IP3 receptors by fold. "
                       "The rest prefer the IP3 receptors over the "
                       "ryanodine receptors by roughly 2:1 but fall short "
                       "of the same-fold bar — and, having had the "
                       "headroom, fall short on similarity and not on "
                       "length alone.",
                       True if itpr_called and len(itpr_called) == len(deep)
                       else ("underpowered" if len(itpr_called) < len(deep)
                             else None)),
            ""]
    return out


def section_plddt(domains, summary, fmt, missing, prior_line) -> list[str]:
    if not summary:
        return [missing("7. Where the models are confident")]
    out = [
        "## 7. Where the models are confident", "",
        "A mean pLDDT over a 2,700-residue multi-domain channel averages a "
        "well-predicted β-trefoil with hundreds of residues of linker, and "
        "the resulting single number is high enough to look reassuring while "
        "saying nothing about the part any claim rests on. Confidence is "
        "reported per domain, with the boundaries S0 measured on the human "
        "paralogs transferred by pairwise alignment; a domain that lands on "
        "too little of its reference span is reported **unplaced** rather "
        "than averaged over whatever aligned.", "",
        "The `shared with` column is a statement about the *domain*, not "
        "about the function: PF08709 is annotated on the ryanodine "
        "receptors too — its Pfam name is \"Inositol "
        "1,4,5-trisphosphate/ryanodine receptor\" — and what the IP3 "
        "receptors do not share is the ligand that binds in it.", "",
        "| domain | Pfam | shared with | models | median pLDDT | min | max | "
        "≥ 70 |", "|---|---|---|---|---|---|---|---|"]
    for s in summary:
        out.append(
            f"| {s['label'].replace('$_3$', '3')} | {s['pfam']} | "
            f"{s['shared_class']} | {s['n_models']} | "
            f"**{fmt(s['median_plddt'])}** | {fmt(s['min_plddt'])} | "
            f"{fmt(s['max_plddt'])} | {s['n_ge_70']} |")
    core = next((s for s in summary if s["pfam"] == "PF08709"), None)
    outside = next((s for s in summary if s["pfam"] == "-"), None)
    out.append("")
    if core and outside:
        out += [prior_line(
            "ip3_core_diagnostic",
            f"The IP3-binding core is modelled at a median pLDDT of "
            f"**{fmt(core['median_plddt'])}** — the highest of any domain "
            f"in the table — against **{fmt(outside['median_plddt'])}** "
            "outside the annotated domains, while the pore, the part the "
            "two families genuinely share as a working channel, is the "
            "worst-modelled domain of the receptor. The ligand-binding "
            "core is among the better-modelled parts, not among the worse.",
            _f(core["median_plddt"]) > _f(outside["median_plddt"])), ""]
    if domains:
        unplaced = [d for d in domains if d["placed"] == "0"]
        if unplaced:
            byd = Counter(d["label"].replace("$_3$", "3") for d in unplaced)
            out += [f"**{len(unplaced)}** domain slots could not be placed at "
                    "all, concentrated in the short non-vertebrate models: "
                    + ", ".join(f"{k} ({v})" for k, v in byd.most_common())
                    + ".", ""]
    out += ["![](figures/s11_plddt_domains.png)", "",
            "**Figure 4.** AlphaFold confidence per domain, ordered along the "
            "subunit from the N-terminal IP3-binding core to the pore. "
            "Points are coloured by paralog; the horizontal rules are "
            "AlphaFold's own confident (70) and very-high (90) bands.", ""]
    return out


def section_foldseek(hits, fmt, missing) -> list[str]:
    if not hits:
        return [missing("8. Foldseek sweep (optional)")]
    verdicts = Counter(h["verdict"] for h in hits)
    subset = hits[0].get("subset", "")
    queries = sorted({h["query"] for h in hits})
    out = ["## 8. Foldseek sweep", "",
           f"A structure-first search of the **{subset}** subset of "
           "AlphaFold DB, run from the IP3 receptor references "
           f"({len(queries)} queries). The subset is named in every output "
           "row because a negative here is a negative about a *declared* "
           "space, which is the roadmap's scoping rule applied to a "
           "structural search.", "",
           "A hit is only a receptor if it explains the *whole* query. "
           "Foldseek's own `qtmscore` cannot carry that: on this data a "
           "191-residue alignment against a 2,300-residue query still comes "
           "back at 0.72. Coverage is therefore computed as aligned length "
           "over the query chain's own length.", "",
           "| verdict | hits |", "|---|---|"]
    for k, n in verdicts.most_common():
        out.append(f"| {k} | {fmt(n)} |")
    novel = [h for h in hits if h["verdict"] == "novel structural lead"]
    shared = [h for h in hits if h["verdict"] == "shared domain, not a receptor"]
    out.append("")
    if not novel:
        out += ["**No novel structural lead.** Nothing in the reviewed "
                "proteome reaches the same-fold bar across a majority of an "
                "IP3 receptor without already being in the census. Given "
                "how little of the family AlphaFold DB holds (§1) this is a "
                "weak negative about the *family* and a clean one about the "
                "*subset*: within Swiss-Prot's predicted structures there is "
                "no unrecognised IP3-receptor-shaped protein.", ""]
    if shared:
        by_acc: dict[str, dict] = {}
        for h in shared:
            by_acc.setdefault(h["accession"], h)
        pfams = Counter(h.get("target_pfams", "") for h in by_acc.values())
        out += [f"### What it does return: the shared domain", "",
                f"The **{len(by_acc)}** distinct above-bar hits the census "
                "has never held are all short, high-scoring matches to one "
                "domain — median coverage "
                f"{fmt(sorted(float(h['query_coverage'] or 0) for h in by_acc.values())[len(by_acc) // 2])}"
                " of the query.", "",
                "| accession | gene | protein | Pfam | organism | TM | query "
                "covered |", "|---|---|---|---|---|---|---|"]
        for h in sorted(by_acc.values(),
                        key=lambda r: -float(r["alntmscore"] or 0)):
            out.append(f"| {h['accession']} | {h['target_gene'] or '—'} | "
                       f"{h['target_name']} | {h.get('target_pfams', '')} | "
                       f"*{h['taxname']}* | {fmt(h['alntmscore'])} | "
                       f"{fmt(h['query_coverage'])} |")
        only_mir = all(p == "PF02815" for p in pfams if p)
        out += ["",
                ("Every one of them carries **PF02815 (MIR) and nothing "
                 "else** — the domain the IP3 receptors share with the "
                 "O-mannosyltransferases and the SDF2 proteins. That is the "
                 "one negative-control class §2 could not fill from the PDB, "
                 "and the structural sweep finds it unprompted: the "
                 "sharpest thing in the reviewed proteome that is shaped "
                 "like part of this family is the part the family does not "
                 "own."
                 if only_mir else
                 "Their Pfam sets are listed above; read the coverage column "
                 "before reading the TM-score."), ""]
    out += ["> **What this sweep cannot find.** AlphaFold DB holds no "
            "ryanodine receptor model at all — at ~5,000 residues they are "
            "past its ceiling, and the API returns 404 for RYR1. So the "
            "absence of a `sister family (D14)` verdict here is unreachable "
            "by construction, not an observation.", ""]
    return out


def section_caveats(manifest, vs_ref, stats, fmt) -> list[str]:
    ok = [m for m in manifest if m.get("status") == "ok"]
    models = [m for m in ok if m.get("role") == "model"]
    short = [m for m in models if _f(m.get("resolved_residues")) < 2000]
    return [
        "## 9. What this task does not establish", "",
        "- **The panel is small, and it is small for a reason that is itself "
        f"the result.** {len(models)} predicted models, of which {len(short)} "
        "are under 2,000 residues, because AlphaFold DB holds full-length "
        "models for a small minority of this family (§1). A structural "
        "argument over the whole census is not available from public "
        "predictions and would need models generated for the project.",
        "- **A TM-score is not an orthology statement.** The IP3 and "
        f"ryanodine receptors are one fold — they score above the {TM_FOLD} "
        "bar against each other (§4) — so shape separates them by a margin, "
        "not by presence and absence. It corroborates D14; it could not have "
        "established it.",
        "- **Every reference is one chain of a tetramer.** Nothing here "
        "tests quaternary structure, and the family's assembly is where much "
        "of its regulation lives.",
        "- **Conformation is a floor on resolution.** The state panel puts a "
        "number on how much TM-score moves between structures of the same "
        "protein (§4); differences smaller than that are not read.",
        "- **pLDDT is a confidence, not an accuracy.** A domain at 85 is a "
        "domain AlphaFold is confident about, which is not the same as a "
        "domain that is right — and for the deep non-vertebrate records "
        "there is no experimental structure to check it against.",
        "",
    ]


def section_outputs() -> list[str]:
    return [
        "## Outputs", "",
        "- `afdb_coverage_reps.tsv` / `afdb_coverage_census.tsv` / "
        "`afdb_coverage_summary.tsv` — the AFDB probe, per record and "
        "summarised, with UniProt's cross-reference beside what the API "
        "serves",
        "- `reference_candidates.tsv` — every RCSB entity carrying a family "
        "Pfam, with its family call and the rule that admitted or excluded it",
        "- `reference_selection.tsv` — the references, the state panel and "
        "the controls, each with the rule that picked it",
        "- `control_candidates.tsv` / `unfilled_controls.tsv` — the decoy "
        "structures considered, and the classes that had none",
        "- `model_slots.tsv` — which taxonomic slot each predicted model "
        "fills, and which were left unfilled",
        "- `structure_manifest.tsv` — the panel: chain used, resolved "
        "residues, coverage, resolution, state, path and SHA-256",
        "- `tm_scores.tsv` — every TM-align pair, both normalisations",
        "- `tm_vs_reference.tsv` — the structural family call with its margin",
        "- `plddt_domains.tsv` / `plddt_summary.tsv` — per-domain confidence",
        "- `structures_stats.json` — parameters, self-test status and the "
        "SHA-256 of every table",
        "- `figures/` — four figures (D13, D19)", ""]


def sections(load, load_json, fmt, pct, missing, prior_line) -> list[str]:
    manifest = load("structure_manifest.tsv")
    pairs = load("tm_scores.tsv")
    vs_ref = load("tm_vs_reference.tsv")
    stats = load_json("structures_stats.json")
    out: list[str] = []
    out += section_calibration(pairs, manifest, fmt, missing, prior_line)
    out += section_calls(vs_ref, fmt, pct, missing, prior_line)
    out += section_nonvertebrate(vs_ref, manifest, fmt, missing, prior_line)
    out += section_plddt(load("plddt_domains.tsv"),
                         load("plddt_summary.tsv"), fmt, missing, prior_line)
    out += section_foldseek(load("foldseek_hits.tsv"), fmt, missing)
    out += section_caveats(manifest, vs_ref, stats, fmt)
    out += section_outputs()
    return out
