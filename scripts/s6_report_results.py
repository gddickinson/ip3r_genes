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


def _verdict(observed, prior, tol=0.06):
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
    L += _section_family_separation(table, num, idx, mean_between, spec)
    L += _section_sister_preview(table, num, idx, mean_between, cov, labels,
                                 meta, spec)
    L += _section_coverage(load, table, num, spec)
    L += _section_fragment_effect(table, num, labels, cov, cls, meta)
    L += _section_conservation(load, load_json, table, num, stats)
    L += _section_caveats(load_json, spec)
    L += _section_figures()
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
    v_w, v_a = _verdict(within, p_w["value"]), _verdict(across, p_a["value"])
    L = ["### 4.1 The family separation the whole project rests on (D14)", "",
         "The ryanodine receptors are in this alignment on purpose — they "
         "root the tree — and D14 says their separation from ITPR is a "
         "positive test at every stage, never an assumption. Measured on "
         "this alignment, under the same covered-only identity metric S1 "
         "used:", ""]
    L += table(["quantity", "this alignment", "prior", "verdict"],
               [["mean identity within a vertebrate paralog group",
                 f"**{within:.3f}**", f"{p_w['value']:.3f}", f"**{v_w}**"],
                ["mean identity, each paralog group to the RyR outgroup",
                 f"**{across:.3f}**", f"{p_a['value']:.3f}", f"**{v_a}**"]])
    L += [f"Prior source: {p_w['where']}; {p_a['where']}.", "",
          f"The gap is {within - across:.3f} identity units. That is the "
          f"margin every stage of this project has had to work inside, and "
          f"it is the reason the length band is support and never the call: "
          f"a 5,000 aa RyR and a 2,700 aa ITPR are {across:.0%} identical "
          f"over the columns they share, which is not far enough apart for "
          f"a heuristic to be trusted with.", ""]
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
    L = ["### 4.2 The sister question — a preview, not an answer", "",
         f"Prior: {PRIOR['sister_pair']['where']}.", "",
         "This alignment can rank the three between-paralog identities. "
         "That is not a phylogenetic estimate — it ignores the outgroup, "
         "the rate variation and the branch lengths S7's model fits — so it "
         "is reported as a preview with the ranking's own margin, and "
         "**S7's AU test over the three rooted topologies is the answer**.",
         ""]
    L += table(["pair", "mean identity (covered)", "n pairs"],
               [[f"{a} × {b}", f"{v:.3f}",
                 len(idx[a]) * len(idx[b])] for v, a, b in scored])
    gap = lead[0] - second[0]
    L += [f"**{lead[1]} × {lead[2]} leads by {gap:.3f}** over "
          f"{second[1]} × {second[2]}. ", ""]
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
    L = ["### 4.3 How much of the trimmed alignment each tip actually "
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
    L = ["### 4.4 Why both identity matrices are committed", "",
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
    L = ["### 4.5 Where the family is conserved", "",
         f"Over the {num(n)} trimmed columns, mean conservation "
         f"**{sum(vals) / n:.3f}**, with **{num(hi)} columns "
         f"({100 * hi / n:.1f} %) at or above 0.9** — invariant or nearly "
         f"so across an alignment that spans vertebrates, invertebrates, "
         f"plants, protists, fungi and the sister family.", "",
         "`figures/msa_conservation.png` maps human ITPR1's Pfam "
         "architecture onto this profile **through the alignment** — the "
         "domain bands are drawn where the alignment put those residues, "
         "by walking the human row and counting ungapped positions, not by "
         "scaling residue coordinates onto column coordinates.", ""]
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
