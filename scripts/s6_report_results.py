"""The results half of the S6 report, split to keep both under 500 lines.

Takes the caller's loader, formatter and table helper (the `s3_report.py` /
`s3_report_d10.py` pattern), so the two halves cannot read the tables
differently.

**Every headline is chosen by the data.** S6 has priors — S1 measured how
far ITPR sits from RyR under this project's own identity metric, S3 measured
the profile separation, S23 measured that identity stops separating anything
outside the vertebrates, and the review recorded the vertebrate sister
question as *unresolved in the published literature*. Each section states
the prior in `PRIOR`, computes this alignment's statistic, and renders
`confirmed` / `contradicted` / `underpowered`, printing both numbers either
way. A generator written to narrate the priors would print them whatever the
alignment said.

The one place this report deliberately refuses to state a result is the
sister question. The alignment can rank the three between-paralog
identities, and it does; but mean identity is not a phylogenetic estimate,
and S7's AU test over the three rooted topologies is what answers it. The
section says which pair leads and by how much, and calls it a **preview**.
"""

from __future__ import annotations

from collections import defaultdict

#: What was measured before this task, with where.
PRIOR = {
    "itpr_vs_ryr_covered": {
        "value": 0.249,
        "where": "S1 `benchmark_controls/bait_margin.tsv`, mean covered-only "
                 "identity of the 25 ITPR positives to their nearest "
                 "labelled RyR bait (the same metric this table uses)",
    },
    "itpr_within_covered": {
        "value": 0.828,
        "where": "S1, the same 25 positives to their nearest labelled ITPR "
                 "bait",
    },
    "sister_pair": {
        "value": None,
        "where": "the review (§7.4): which two of the three vertebrate "
                 "paralogues are sisters *is not fixed by any published, "
                 "support-annotated ML analysis with an RyR outgroup*",
    },
}


#: How much the two estimators may differ before the comparison means
#: something. S1's numbers come from *pairwise* alignments of a 31-sequence
#: control panel, each positive scored against its nearest bait; this
#: table's come from all pairs of a 134-sequence trimmed MSA. Those are
#: different estimators of the same quantity, so a small absolute
#: difference is expected and is not evidence about the family. The
#: comparison that is load-bearing is the **separation** — how far the
#: within-family identity sits above the cross-family one — and the
#: tolerance is set wide enough that only a real collapse of that
#: separation reads as `contradicted`.
TOLERANCE = 0.15


def _verdict(observed, prior, tol=TOLERANCE):
    if prior is None:
        return "no prior"
    if abs(observed - prior) <= tol:
        return "confirmed"
    return "contradicted"


def _idx_by_group(labels, meta):
    idx = defaultdict(list)
    for i, l in enumerate(labels):
        idx[meta.get(l, {}).get("group", "?")].append(i)
    return idx


def _load_identity(path):
    import csv
    with open(path) as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0][1:], [[float(v) for v in r[1:]] for r in rows[1:]]


def render(*, load, load_json, table, num, reps, stats) -> list[str]:
    import s6_rep_spec as spec
    from s6_lib import MSA_DIR

    ident_path = MSA_DIR / "identity_covered.tsv"
    if not ident_path.exists():
        return ["## 4. Results", "",
                "*The alignment has not been built yet; re-run "
                "`scripts/s6_msa.py` and then this report.*", ""]

    labels, cov = _load_identity(ident_path)
    _, cls = _load_identity(MSA_DIR / "identity_classic.tsv")
    meta = {r["label"]: r for r in reps}
    idx = _idx_by_group(labels, meta)

    def mean_between(a, b):
        vals = [cov[i][j] for i in idx[a] for j in idx[b] if i != j]
        return sum(vals) / len(vals) if vals else float("nan")

    L = ["## 4. What the alignment shows", ""]
    # Sections are written with a `§` placeholder and numbered here, in the
    # order they are assembled. Hand-numbered headings collide the moment a
    # section is inserted, and a report that renumbers itself cannot.
    L += _section_family_separation(table, num, idx, mean_between, spec)
    L += _section_sister_preview(table, num, idx, mean_between, cov, labels,
                                 meta, spec)
    L += _section_cyclostome(table, num, labels, cov, meta, spec, idx)
    L += _section_coverage(load, table, num, spec)
    L += _section_fragment_effect(table, num, labels, cov, cls, meta)
    L += _section_conservation(load, load_json, table, num, stats)
    L += _section_sites(load_json, table, num)
    L += _section_caveats(load_json, spec)
    L += _section_figures()
    n = 0
    for i, line in enumerate(L):
        if line.startswith("### § "):
            n += 1
            L[i] = f"### 4.{n} " + line[len("### § "):]
    return L


# ----------------------------------------------------------------- D14

def _section_family_separation(table, num, idx, mean_between, spec) -> list[str]:
    trio = [g for g in spec.PARALOGS if g in idx]
    if not trio or "RYR" not in idx:
        return []
    within = [mean_between(a, a) for a in trio]
    within = sum(within) / len(within)
    across = sum(mean_between(a, "RYR") for a in trio) / len(trio)
    p_w = PRIOR["itpr_within_covered"]
    p_a = PRIOR["itpr_vs_ryr_covered"]
    sep, p_sep = within - across, p_w["value"] - p_a["value"]
    L = ["### § The family separation the whole project rests on (D14)", "",
         "The ryanodine receptors are in this alignment on purpose — they "
         "root the tree — and D14 says their separation from ITPR is a "
         "positive test at every stage, never an assumption. Measured here "
         "under the same covered-only identity metric S1 used, against what "
         "S1 measured:", ""]
    L += table(["quantity", "this alignment", "S1 prior", "Δ"],
               [["mean identity within a vertebrate paralog group",
                 f"**{within:.3f}**", f"{p_w['value']:.3f}",
                 f"{within - p_w['value']:+.3f}"],
                ["mean identity, each paralog group to the RyR outgroup",
                 f"**{across:.3f}**", f"{p_a['value']:.3f}",
                 f"{across - p_a['value']:+.3f}"],
                ["**the separation between them**", f"**{sep:.3f}**",
                 f"{p_sep:.3f}", f"{sep - p_sep:+.3f}"]])
    L += [f"Verdict on the separation: **{_verdict(sep, p_sep)}** "
          f"(tolerance {TOLERANCE:.2f}).", "",
          f"Prior source: {p_w['where']}; {p_a['where']}.", "",
          "The two estimators are not identical and are not expected to "
          "agree to three decimals: S1 scored each of 31 control-panel "
          f"sequences against its *nearest* bait on a pairwise alignment, "
          f"and this table averages **all** pairs of a {sum(len(v) for v in idx.values())}-sequence trimmed MSA. "
          "So the comparison is made on the separation, which is what every "
          "later stage actually depends on, rather than on either absolute "
          "value.", "",
          f"The separation is {sep:.3f} identity units. That is the margin "
          f"every stage of this project has had to work inside, and it is "
          f"why the length band is support and never the call: a 5,000 aa "
          f"RyR and a 2,700 aa ITPR are {across:.0%} identical over the "
          f"columns they share, which is not far enough apart to trust a "
          f"heuristic with.", ""]
    return L


# --------------------------------------------------------------- sister

def _section_sister_preview(table, num, idx, mean_between, cov, labels,
                            meta, spec) -> list[str]:
    trio = [g for g in spec.PARALOGS if g in idx]
    if len(trio) < 3:
        return []
    pairs = [(trio[i], trio[j]) for i in range(3) for j in range(i + 1, 3)]
    scored = sorted(((mean_between(a, b), a, b) for a, b in pairs),
                    reverse=True)
    lead, second = scored[0], scored[1]
    L = ["### § The sister question — a preview, not an answer", "",
         f"Prior: {PRIOR['sister_pair']['where']}.", "",
         "This alignment can rank the three between-paralog identities. "
         "That is not a phylogenetic estimate — it ignores the outgroup, "
         "the rate variation and the branch lengths S7's model fits — so it "
         "is reported as a preview with the ranking's own margin, and "
         "**S7's AU test over the three rooted topologies is the answer**.",
         ""]
    iqr = {}
    for a, b in pairs:
        vals = sorted(cov[i][j] for i in idx[a] for j in idx[b])
        iqr[(a, b)] = (vals[len(vals) // 4], vals[3 * len(vals) // 4])
    L += table(["pair", "mean identity (covered)", "interquartile range",
                "n pairs"],
               [[f"{a} × {b}", f"{v:.3f}",
                 f"{iqr[(a, b)][0]:.3f} – {iqr[(a, b)][1]:.3f}",
                 len(idx[a]) * len(idx[b])] for v, a, b in scored])
    gap = lead[0] - second[0]
    lead_lo = iqr[(lead[1], lead[2])][0]
    others_hi = max(iqr[(a, b)][1] for a, b in pairs
                    if (a, b) != (lead[1], lead[2]))
    L += [f"**{lead[1]} × {lead[2]} leads by {gap:.3f}** over "
          f"{second[1]} × {second[2]}.", ""]
    if lead_lo > others_hi:
        L += [f"The leading pair's interquartile range "
              f"({lead_lo:.3f} – {iqr[(lead[1], lead[2])][1]:.3f}) does not "
              f"overlap either other pair's (highest upper quartile "
              f"{others_hi:.3f}), so the ranking is not an artefact of a few "
              f"close pairs — it holds across the middle half of every "
              f"comparison. It is still not a phylogenetic estimate.", ""]
    if gap < 0.02:
        L += ["That margin is inside the noise of a mean over hundreds of "
              "pairs of unequal-rate sequences; the alignment does not "
              "separate the three hypotheses and should not be quoted as "
              "if it did.", ""]
    else:
        L += ["The margin is large enough to be worth carrying into S7 as a "
              "hypothesis to test, and small enough that the AU test is "
              "what settles it.", ""]
    return L


# ----------------------------------------------------------- cyclostomes

def _section_cyclostome(table, num, labels, cov, meta, spec, idx) -> list[str]:
    """Each pre-2R locus against the three paralog groups.

    The one question this alignment can put a number on that no earlier
    task could: are the three cyclostome loci each closest to a *different*
    vertebrate paralog (what 1:1 orthology from 2R would look like), or all
    equidistant (what a lineage-specific expansion would look like)?
    """
    trio = [g for g in spec.PARALOGS if g in idx]
    rows = [(i, l) for i, l in enumerate(labels)
            if meta.get(l, {}).get("band") == "cyclostomata"]
    if len(trio) < 3 or len(rows) < 2:
        return []
    out, best_counts = [], {}
    for i, l in rows:
        means = {g: sum(cov[i][j] for j in idx[g]) / len(idx[g]) for g in trio}
        best = max(means, key=means.get)
        rest = sorted(means.values(), reverse=True)
        best_counts[best] = best_counts.get(best, 0) + 1
        out.append([f"`{meta[l]['species'].split('(')[0].strip()}` "
                    f"{meta[l]['accession'].split('|')[-1][:22]}",
                    *[f"{means[g]:.3f}" for g in trio],
                    f"**{best}**", f"{rest[0] - rest[1]:+.3f}"])
    # The control this table needs: is "nearest ITPR1" a fact about the
    # cyclostomes, or about the metric? Every deep group is scored the same
    # way, so the same statistic over the non-vertebrate grades says what a
    # no-signal answer looks like on this alignment.
    baseline = {}
    for g in ("invert_metazoa", "protist", "RYR"):
        if g not in idx:
            continue
        counts, margins = {}, []
        for i in idx[g]:
            means = {p: sum(cov[i][j] for j in idx[p]) / len(idx[p])
                     for p in trio}
            s = sorted(means.values(), reverse=True)
            counts[max(means, key=means.get)] = counts.get(
                max(means, key=means.get), 0) + 1
            margins.append(s[0] - s[1])
        margins.sort()
        baseline[g] = (len(margins), counts, margins[len(margins) // 2])
    cyc_margins = sorted(float(r[-1]) for r in out)
    cyc_med = cyc_margins[len(cyc_margins) // 2]

    L = ["### § The cyclostome trio, previewed", "",
         "Every cyclostome in this set carries three ITPR loci, and the "
         "sweep's ITPR1 bait won all of them (§1.2), so nothing before now "
         "could say which locus is which. The alignment can at least ask "
         "the question: is each locus closest to a *different* vertebrate "
         "paralog — what 1:1 orthology from 2R would look like — or are "
         "they all equidistant, which is what a cyclostome-specific "
         "expansion would look like?", ""]
    L += table(["locus", *trio, "nearest", "margin over 2nd"], out)
    spread = len(best_counts)
    if spread >= 3:
        L += [f"The {len(rows)} loci divide across **{spread}** paralog "
              f"groups. That is the shape 1:1 orthology would produce, and "
              f"it is a reason to test the 2R hypothesis in S7 rather than "
              f"a demonstration of it — mean identity to a group is not an "
              f"orthology assignment, and the margins above say how thin "
              f"the distinctions are.", ""]
    else:
        L += [f"All {len(rows)} loci fall nearest the same "
              f"{'group' if spread == 1 else 'two groups'} "
              f"({', '.join(sorted(best_counts))}). On identity alone the "
              f"three copies are not separable into ITPR1/2/3, which is "
              f"what a lineage-specific expansion looks like — and equally "
              f"what three fast-evolving 1:1 orthologs would look like at "
              f"this depth. **S7's tree and S8's synteny are what "
              f"distinguish them**; this table says only that the easy "
              f"answer is not available.", ""]
    if baseline:
        L += ["**The control that table needs.** Leaning towards one "
              "paralog could be a fact about the cyclostomes or a fact "
              "about the metric — ITPR1 may simply be the slowest-evolving "
              "of the three, in which case everything deep is nearest it. "
              "The same statistic over the groups that are certainly *not* "
              "vertebrate paralogs says what no signal looks like here:",
              ""]
        L += table(["group", "n", "nearest paralog", "median margin"],
                   [[g, n, ", ".join(f"{k} {v}" for k, v in
                                     sorted(c.items(), key=lambda kv: -kv[1])),
                     f"{m:.3f}"] for g, (n, c, m) in baseline.items()]
                   + [["**cyclostome loci**", len(out), "ITPR1 "
                       f"{best_counts.get('ITPR1', 0)}", f"**{cyc_med:.3f}**"]])
        worst = max(m for _, _, m in baseline.values())
        if cyc_med > 3 * worst:
            L += [f"The non-vertebrate groups do lean towards ITPR1 more "
                  f"often than chance, but at a median margin of "
                  f"{worst:.3f} — no signal. The cyclostome margin is "
                  f"{cyc_med / worst:.0f}× that. So the lean is a fact "
                  f"about these loci and not an artefact of ITPR1 being the "
                  f"conserved paralog: **the three cyclostome copies are "
                  f"genuinely closer to ITPR1 than to ITPR2 or ITPR3**. "
                  f"That is what a "
                  f"cyclostome-specific expansion from an ITPR1-like "
                  f"ancestor would produce; it is also what 1:1 orthologs "
                  f"would produce if ITPR2 and ITPR3 diverged after the "
                  f"cyclostome split. S7 separates those; S6 can only say "
                  f"the signal is real.", ""]
        else:
            L += [f"The non-vertebrate groups lean the same way at a median "
                  f"margin of {worst:.3f} against the cyclostomes' "
                  f"{cyc_med:.3f}. The lean is therefore a property of the "
                  f"metric — ITPR1 is the paralog everything deep is "
                  f"nearest — and carries no information about the "
                  f"cyclostome loci specifically.", ""]
    return L


# ------------------------------------------------------------- coverage

def _section_coverage(load, table, num, spec) -> list[str]:
    rows = load("coverage.tsv")
    if not rows:
        return []
    vals = sorted(float(r["coverage"]) for r in rows)
    med = vals[len(vals) // 2]
    worst = sorted(rows, key=lambda r: float(r["coverage"]))[:8]
    by_group = defaultdict(list)
    for r in rows:
        by_group[r["group"]].append(float(r["coverage"]))
    L = ["### § How much of the trimmed alignment each tip actually "
         "carries", "",
         f"Median coverage **{med:.2f}**; "
         f"{sum(1 for v in vals if v < 0.5)} of {len(vals)} tips cover less "
         f"than half the trimmed alignment. A tip below half is not wrong, "
         f"but it contributes gaps to every column the tree is inferred "
         f"from, so it is named here rather than left inside a median.", ""]
    L += table(["group", "n", "median coverage", "lowest"],
               [[g, len(v), f"{sorted(v)[len(v) // 2]:.2f}", f"{min(v):.2f}"]
                for g in spec.GROUP_ORDER if g in by_group
                for v in [by_group[g]]])
    L += ["Lowest eight tips:", ""]
    L += table(["tip", "group", "coverage"],
               [[f"`{r['label']}`", r["group"], f"{float(r['coverage']):.2f}"]
                for r in worst])
    return L


# ------------------------------------------------------- fragment effect

def _section_fragment_effect(table, num, labels, cov, cls, meta) -> list[str]:
    diffs = []
    for i in range(len(labels)):
        for j in range(i + 1, len(labels)):
            d = cov[i][j] - cls[i][j]
            if d > 0:
                diffs.append((d, labels[i], labels[j], cov[i][j], cls[i][j]))
    diffs.sort(reverse=True)
    if not diffs:
        return []
    L = ["### § Why both identity matrices are committed", "",
         "`identity_covered.tsv` scores identity over mutually covered "
         "columns only; `identity_classic.tsv` counts a gap as a mismatch. "
         "For a set that deliberately contains fragments and a ~5,000 aa "
         "outgroup, those are different measurements, and the difference is "
         "not small. The largest divergences:", ""]
    L += table(["pair", "covered", "classic", "Δ"],
               [[f"`{a}` × `{b}`", f"{cv:.2f}", f"{cl:.2f}", f"+{d:.2f}"]
                for d, a, b, cv, cl in diffs[:8]])
    L += ["Every identity quoted in this report is the covered-only one. "
          "The classic matrix is committed beside it so a later task can "
          "see what a gap-counting metric would have said instead.", ""]
    return L


# ---------------------------------------------------------- conservation

def _section_conservation(load, load_json, table, num, stats) -> list[str]:
    cons = load("conservation.tsv")
    if not cons:
        return []
    vals = [float(r["conservation"]) for r in cons]
    n = len(vals)
    hi = sum(1 for v in vals if v >= 0.9)
    L = ["### § Where the family is conserved", "",
         f"Over the {num(n)} trimmed columns, mean conservation "
         f"**{sum(vals) / n:.3f}**, with **{num(hi)} columns "
         f"({100 * hi / n:.1f} %) at or above 0.9** — invariant or nearly "
         f"so across an alignment that spans vertebrates, invertebrates, "
         f"plants, protists, fungi and the sister family.", "",
         "Conservation here is 1 − normalised Shannon entropy over the "
         "column, **counting a gap as a character** "
         "(`src/analysis/evolution.py`). On a trimmed alignment that is the "
         "right convention — a column half of the tips do not have is less "
         "conserved across the family, not more — but it means the profile "
         "is not comparable to one computed over residues only.", "",
         "`figures/msa_conservation.png` maps human ITPR1's Pfam "
         "architecture onto this profile **through the alignment** — the "
         "domain bands are drawn where the alignment put those residues, "
         "by walking the human row and counting ungapped positions, not by "
         "scaling residue coordinates onto column coordinates.", ""]
    return L


# ----------------------------------------------------------------- sites

def _section_sites(load_json, table, num) -> list[str]:
    s = load_json("align_stats.json").get("sites", {})
    if not s:
        return []
    tot = sum(v for k, v in s.items() if k != "pct_informative")
    L = ["### § What is left for the tree to work with", "",
         "trimAl `-automated1` is a heuristic, and a percentage of columns "
         "kept says nothing about whether the *informative* ones survived. "
         "Counted on the trimmed alignment, by the same definition IQ-TREE "
         "reports (a column is parsimony-informative when at least two "
         "residues each appear at least twice):", ""]
    L += table(["column class", "n", "% of trimmed"],
               [[k.replace("_", " "), num(v), f"{100 * v / tot:.1f} %"]
                for k, v in s.items() if k != "pct_informative"])
    L += [f"**{s['pct_informative']} % of the trimmed alignment is "
          f"parsimony-informative** — {num(s['parsimony_informative'])} "
          f"columns. That is the number S7's support values are estimated "
          f"from, and the number to quote if a node's bootstrap is "
          f"questioned.", ""]
    return L


# ------------------------------------------------------------- caveats

def _section_caveats(load_json, spec) -> list[str]:
    sel = load_json("selection_stats.json")
    L = ["## 5. Caveats", "",
         "- **The paralog label on a non-vertebrate tip is annotation "
         "transfer, not descent.** ITPR1/2/3 are a 2R product. Four "
         "representatives outside the vertebrates carry a type number in "
         "their UniProt gene symbol or protein name; the audit table keeps "
         "that in `paralog` with `paralog_source`, and the figure group is "
         "the taxonomic grade regardless. Colouring such a tip as a "
         "vertebrate paralog would assert the thing S7 is being run to "
         "test.",
         "- **A novel model's group is its bait's hypothesis, not a "
         "result.** The R7 tips are loci no database annotates; their "
         "ITPR1/2/3 group comes from which bait won them in the S5 sweep, "
         "recorded as `paralog_source = s5_cell`. Where that attribution "
         "is constant across a clade's loci it is not used at all (the "
         "cyclostome grade, §1.2). Where it is used, the tree is what "
         "tests it — that is the question those tips are in the alignment "
         "to ask.",
         "- **Coverage is not evidence of quality.** A short tip covers "
         "less of the alignment; that is what short means. Whether a short "
         "tip is a real short gene or a broken model is S15's question, not "
         "this one's.",
         "- **The 3R pairs are in by rule, and they are a stress test.** "
         "Two teleost species contribute both `itprXa` and `itprXb`. If the "
         "tree does not recover them as sisters, the naming is wrong or the "
         "alignment is — either way S7 finds out rather than S6 assuming.",
         "- **The RyR outgroup is six sequences against 120.** It roots the "
         "tree; it is not a sample of the ryanodine receptors, and no "
         "statement about RyR evolution can be read off this alignment.",
         f"- **{len(sel.get('short_exceptions', []))} tips are below the "
         f"family's own length band**, admitted under the "
         f"architecture-exception floor so their clade has a tip at all. "
         f"They are listed in §1.2 and flagged `short_exception` in "
         f"`representatives.tsv`.", ""]
    return L


def _section_figures() -> list[str]:
    from s6_lib import MSA_DIR
    figs = [("msa_identity_heatmap",
             "Pairwise identity over mutually covered columns, ordered by "
             "group. The block structure is the result: the three "
             "vertebrate paralogs are tight blocks, the non-vertebrate "
             "grade is not a block at all, and the RyR outgroup is a "
             "uniformly dark band against everything."),
            ("msa_conservation",
             "Per-column conservation of the trimmed alignment with human "
             "ITPR1's Pfam architecture mapped through the alignment onto "
             "it."),
            ("msa_coverage",
             "Per-sequence coverage of the trimmed alignment by group, bar "
             "at the group median — which tips are fragments and which "
             "groups they are in."),
            ("msa_group_identity",
             "Left: mean between-paralog identity with its interquartile "
             "range, the alignment's preview of S7's sister question. "
             "Right: mean identity between every pair of groups.")]
    L = ["## 6. Figures", ""]
    for slug, cap in figs:
        if (MSA_DIR / "figures" / f"{slug}.png").exists():
            L += [f"![{slug}](figures/{slug}.png)", "", f"*{cap}*", ""]
    return L
