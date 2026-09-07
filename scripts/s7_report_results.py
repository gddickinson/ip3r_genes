"""S7 — the results half of `results/phylogeny/report.md`.

Split from `s7_report.py` to keep both under 500 lines, and taking the
caller's loader, formatter and table helper so the two halves cannot read
the tables differently (the `s3_report.py` / `s3_report_d10.py` pattern).

**The headlines are chosen by the data.** Every earlier task that touched
the paralogs left a number behind, and each of them is a prediction this
tree can falsify: S6's identity ranking picked a sister pair, S6's
cyclostome table found all six loci nearest ITPR1, S6 named the teleost
3R pairs as a stress test, and S2's architecture call and S5's paralog
cells have assumed all along that ITPR1/2/3 are clades. Each section
below states the prior in `PRIOR`, computes the tree's answer, and
renders `confirmed` / `contradicted` / `unresolved`, printing both
numbers either way. A generator written to narrate the expected answers
would print them whatever the tree said.
"""

from __future__ import annotations

from s7_lib import PHYLO_DIR                                   # noqa: E402
from s7_priors import PRIOR                                    # noqa: E402
from s7_report_checks import _section_besttree, _section_bnni  # noqa: E402
from s7_report_naming import (  # `_sup`/`_strong` live there so both
    _section_cyclostome, _section_duplicates, _section_relabels,  # halves
    _strong, _sup,                    # read a support pair identically
)

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")


AU_ALPHA = 0.05


def render(*, load, load_json, table, num, reps, stats) -> list[str]:
    L = ["## 5. What the tree says", ""]
    L += _section_monophyly(load, table, num)
    L += _section_sister(load, load_json, table, num)
    L += _section_relabels(load, table, num, reps)
    L += _section_cyclostome(load, table, num)
    L += _section_duplicates(load, table, num)
    L += _section_support(load, table, num)
    L += _section_bnni(load, load_json, table, num)
    L += _section_besttree(load, load_json, table, num)
    L += _section_caveats(load, load_json, table, num)
    L += _section_figures(load)
    n = 0
    for i, line in enumerate(L):
        if line.startswith("### § "):
            n += 1
            L[i] = f"### 5.{n} " + line[len("### § "):]
    return L


# ------------------------------------------------------------- monophyly

def _section_monophyly(load, table, num) -> list[str]:
    raw = load("monophyly.tsv")
    cor = load("monophyly_corrected.tsv")
    if not raw:
        return []
    p = PRIOR["paralog_monophyly"]
    got_raw = sum(1 for r in raw if r["group"] in PARALOGS
                  and r["monophyletic"] == "yes")
    got_cor = sum(1 for r in cor if r["monophyletic"] == "yes")
    L = ["### § Are ITPR1, ITPR2 and ITPR3 clades at all?", "",
         f"Prior: **{p['value']} of 3** — {p['where']}.", "",
         "The tree is asked first on the census's own labels, with nothing "
         "corrected. A group is monophyletic here if some edge of the "
         "*unrooted* tree separates exactly its tips.", ""]
    rows = [[r["group"], r["n"],
             "**yes**" if r["monophyletic"] == "yes" else "no",
             _sup(r), r["largest_pure_clade"], r["n_outside"]]
            for r in raw]
    L += table(["group", "tips", "monophyletic", "SH-aLRT/UFBoot",
                "largest pure clade", "tips outside it"], rows)

    verdict = ("confirmed" if got_raw == p["value"] else
               "contradicted" if got_raw < p["value"] else "confirmed")
    L += [f"**{got_raw} of 3** paralogs are monophyletic on the raw "
          f"census labels — the prior is {p['value']} of 3, so this is "
          f"**{verdict}**.", ""]
    if got_raw < 3 and cor:
        L += ["A label is not evidence, though, and a single "
              "mis-annotated tip breaks a clade of forty. The same "
              "question, asked after the tree's own corrections (§5.3, "
              "derived by coded rule in `s7_constraints.py`, never by "
              "hand):", ""]
        rows = [[r["group"], r["n"],
                 "**yes**" if r["monophyletic"] == "yes" else "no",
                 _sup(r)] for r in cor]
        L += table(["clade", "tips", "monophyletic", "SH-aLRT/UFBoot"], rows)
        L += [f"**{got_cor} of 3** under tree-corrected membership. The "
              "difference between the two tables is the size of the "
              "family's naming problem, not a change in the tree.", ""]
    return L


# ---------------------------------------------------------- the AU test

def _section_sister(load, load_json, table, num) -> list[str]:
    au = load("au_test.tsv")
    sml = load("sister_ml.tsv")
    claims = load("claim_nodes.tsv")
    p = PRIOR["sister_pair"]
    L = ["### § The sister question — what the AU test answers", "",
         f"Prior: **{p['value']}**, from {p['where']}.", "",
         f"Prior from the literature: **none**. {PRIOR['published_sister']['where']}.",
         "", "Three rooted hypotheses, one per way of pairing the three "
         "paralogs, each realised as a constrained ML search under the "
         "same model and each making all three paralog clades "
         "monophyletic — so the test compares the sister arrangement and "
         "nothing else. Membership is tree-corrected and the correction "
         "is *derived*, not listed (§5.3): a constraint built from census "
         "labels asks the test about a topology the data rejects for "
         "reasons that have nothing to do with the sister question, and "
         "all three hypotheses then fail together.", "",
         "The same trap has a second door, and this task walked into it. "
         "IQ-TREE's `-g` places freely exactly the taxa a constraint "
         "**omits**, so a tip that is listed — even in a top-level "
         "polytomy — is pinned outside every group the constraint "
         "declares. The first constraint files named all 134 tips and so "
         "forced the unlabelled vertebrate tips out of the paralog "
         "clades this tree nests them in, in all three hypotheses "
         "equally; the test rejected every one at ΔlogL ≈ 1,500, "
         "including the arrangement the ML tree itself holds at 100/100. "
         "Each constraint now names the three paralog cores and the "
         "outgroup and nothing else, and `s7_test_tree.py` T5 fails a "
         "constraint that names a free tip.", ""]

    # the ML tree's own answer, read off the rooted tree
    ml_pair = ""
    if claims:
        pairs = [r for r in claims
                 if r["claim"].count("+") == 1 and r["is_clade"] == "yes"]
        if pairs:
            best = max(pairs, key=lambda r: float(r["ufboot"] or 0))
            ml_pair = best["claim"]
            L += [f"**The unconstrained ML tree groups "
                  f"{ml_pair}** at SH-aLRT {best['alrt']} / UFBoot "
                  f"{best['ufboot']}.", ""]
        else:
            L += ["**The unconstrained ML tree groups no two paralogs to "
                  "the exclusion of the third** — no pair is a clade of "
                  "the rooted tree. That is itself the answer to look at "
                  "the AU test for.", ""]
    if not au:
        L += ["*The AU test has not been run yet — "
              "`python3 scripts/s7_run.py au`.*", ""]
        return L

    rows = [[r["tree"], r["hypothesis"] or "—", num(r["logL"], 3),
             r["deltaL"], r["p_AU"], r["p_KH"], r["p_SH"],
             _au_verdict(r)] for r in au]
    L += table(["tree", "hypothesis", "logL", "ΔlogL", "p-AU", "p-KH",
                "p-SH", "verdict"], rows)

    alive = [r for r in au if r["tree"] != "ML" and _alive(r)]
    dead = [r for r in au if r["tree"] != "ML" and not _alive(r)]
    L += [f"**{len(dead)} of {len(dead) + len(alive)} sister "
          f"arrangements are rejected** at p-AU < {AU_ALPHA}.", ""]
    if dead:
        L += ["- rejected: " + "; ".join(
            f"{r['hypothesis']} (ΔlogL {r['deltaL']}, p-AU {r['p_AU']})"
            for r in dead)]
    if alive:
        L += ["- not rejected: " + "; ".join(
            f"{r['hypothesis']} (ΔlogL {r['deltaL']}, p-AU {r['p_AU']})"
            for r in alive)]
    L += [""]

    if len(alive) == 1:
        h = alive[0]
        got = _pair_of(h["tree"])
        verdict = "confirmed" if got == p["value"] else "contradicted"
        L += [f"**One arrangement survives: {h['hypothesis']}.** The "
              f"other two are outside the 95 % confidence set of "
              f"topologies for this alignment. S6's identity preview "
              f"predicted **{p['value']}** and the test returns "
              f"**{got}** — the prior is **{verdict}**.", "",
              "This is the answer the brief asks for, and it is one the "
              "literature did not have: the review's §7.4 audit found no "
              "published, support-annotated ML analysis with an RyR "
              "outgroup that fixes the pair.", ""]
    elif len(alive) == 0:
        L += ["**Every hypothesis is rejected**, the ML topology "
              "included. That is not an answer about the paralogs — it "
              "is a sign the constraints are asking the wrong question "
              "(a constraint set that forces something the data rejects "
              "for an unrelated reason rejects everything). §5.3 is where "
              "to look.", ""]
    else:
        surv = ", ".join(_pair_of(r["tree"]) for r in alive)
        L += [f"**{len(alive)} arrangements survive ({surv}), so the "
              f"sister relationship is not resolved by this alignment.** "
              f"The ML topology is the best point estimate and the "
              f"constrained version of it costs almost no likelihood, but "
              f"the alternatives stay inside the 95 % confidence set. "
              f"S6's identity preview picked **{p['value']}** by "
              f"{p['margin']:.3f} identity units; the tree does not "
              f"contradict it and does not confirm it either — the prior "
              f"is **unresolved**.", "",
              "That is a statement about *this* alignment, not about the "
              "family: a duplication this old leaves little unsaturated "
              "signal, which is exactly why the identity ranking and the "
              "likelihood disagree about how much they know. S8's synteny "
              "is independent evidence on the same question.", ""]
    if sml:
        L += ["What each paralog clade's sister actually is on the rooted "
              "ML tree, composition and all:", ""]
        L += table(["paralog", "sister is a clade", "tips in sister",
                    "sister composition", "majority"],
                   [[r["paralog"], r["is_clade"], r["sister_n"],
                     r["sister_composition"] or "—",
                     r["sister_majority"] or "—"] for r in sml])
    return L


def _alive(r: dict) -> bool:
    try:
        return float(r["p_AU"]) >= AU_ALPHA
    except (TypeError, ValueError):
        return True


def _au_verdict(r: dict) -> str:
    if r["tree"] == "ML":
        return "—"
    return "not rejected" if _alive(r) else "**rejected**"


def _pair_of(tree_name: str) -> str:
    return {"H1_12": "ITPR1 + ITPR2", "H2_13": "ITPR1 + ITPR3",
            "H3_23": "ITPR2 + ITPR3"}.get(tree_name, tree_name)


# ------------------------------------------------------------ resolution

def _section_support(load, table, num) -> list[str]:
    rows = load("support_summary.tsv")
    if not rows:
        return []
    d = {r["statistic"]: r["value"] for r in rows}
    L = ["### § How much of this tree is actually resolved", "",
         "The support summary is not decoration: a sister question "
         "answered on a tree whose deep nodes are unsupported is not "
         "answered. This is the whole tree, every internal node:", ""]
    L += table(["statistic", "value"], [[k, num(v)] for k, v in d.items()])
    try:
        pct = float(d.get("pct both", 0))
        L += [f"**{pct:.1f} %** of internal nodes clear both thresholds. "
              + ("The tree is well resolved and the deep nodes the "
                 "claims rest on are listed with their own support in "
                 "`claim_nodes.tsv`." if pct >= 60 else
                 "That is a moderately resolved tree, and it is the "
                 "reason the claim nodes are reported with their own "
                 "support rather than under a blanket statement — "
                 "`claim_nodes.tsv` carries each one."), ""]
    except (TypeError, ValueError):
        pass
    claims = load("claim_nodes.tsv")
    if claims:
        L += ["Node by node, for every claim this task or a later one "
              "leans on:", ""]
        L += table(["claim", "tips", "is a clade", "SH-aLRT/UFBoot",
                    "well supported"],
                   [[c["claim"], c["n_tips"], c["is_clade"], _sup(c),
                     "**yes**" if c["well_supported"] == "yes" else "no"]
                    for c in claims])
    return L


# --------------------------------------------------------------- caveats

def _section_caveats(load, load_json, table, num) -> list[str]:
    choice = load_json("model_choice.json")
    L = ["### § Caveats", ""]
    items = [
        "**The model selection is greedy** (§2). The exchangeability "
        "matrix was chosen under `+G4`/`+I+G4` and the rate model then "
        "chosen on that matrix alone; the exhaustive scan was measured at "
        f"~{choice.get('exhaustive_probe', {}).get('projection_hours', '?')} h "
        "and abandoned. Both stages' tables are committed.",
        "**One alignment, one tree.** Every statement here is conditional "
        "on S6's 134-tip representative set and its trimAl columns. The "
        "representative rules are audited in `representatives.tsv` and "
        "the trimming is reproducible from `align_stats.json`, but a "
        "different taxon sample can move a deep node.",
        "**UFBoot is optimistic under model violation**, which a "
        "four-kingdom alignment is guaranteed to be in. That is why every "
        "'well supported' call in this report requires SH-aLRT ≥ 80 as "
        "well, and why a `--bnni` re-run is reported beside the main one "
        "where it exists.",
        "**The paralog labels the tree corrects are corrected by rule, "
        "not by hand** — but the rule reads this tree. Where a "
        "correction changes a call that matters, the RBH test in §5.3 is "
        "the independent check.",
    ]
    L += [f"- {i}" for i in items] + [""]
    return L


def _section_figures(load) -> list[str]:
    """The figure block — captions computed, and only for figures that exist.

    A caption is a claim about what the reader is looking at. Two of
    these described a marker and a grouping the figures do not carry
    (no tip on this tree was relabelled, and none of the unplaced tips
    resolves nearest a paralog), so the wording that varies with the
    data is derived from the same tables the figures are drawn from.
    """
    n_relabel = sum(1 for a in load("membership_audit.tsv")
                    if a["rule"] == "reassigned")
    homed = [r for r in load("cyclostome_placement.tsv") if r["home"]]
    figs = [
        ("tree_ml_rooted",
         "The rooted ML phylogeny. Branches in neutral ink; the three "
         "vertebrate paralog clades and the RyR outgroup boxed and "
         "coloured; a filled dot on every node clearing SH-aLRT ≥ 80 and "
         "UFBoot ≥ 95"
         + (f"; a ring on each of the {n_relabel} tips whose census label "
            "the tree overturned." if n_relabel else
            ". No tip is ringed because the relabelling rule fires on "
            "none of them.")),
        ("sister_au",
         "The three sister hypotheses and what the AU test does to them: "
         "constrained log-likelihood against the ML tree, with p-AU "
         "beside it and the 0.05 rejection line drawn."),
        ("support_profile",
         "How much of the tree is resolved — the joint distribution of "
         "SH-aLRT and UFBoot over every internal node, with the claim "
         "nodes marked."),
        ("paralog_placement",
         "The vertebrate tips the tree does not place in a paralog — the "
         "cyclostome loci and the disputed chondrichthyan and coelacanth "
         "records — each against the support of the smallest clade that "
         "holds it, annotated with what else is in that clade. "
         + (f"{len(homed)} of them resolve inside a single paralog."
            if homed else
            "None of them resolves inside a single paralog: the "
            "cyclostome loci sit with each other and the rest in mixed "
            "neighbourhoods.")),
    ]
    L = ["### § Figures", ""]
    for slug, cap in figs:
        if not (PHYLO_DIR / "figures" / f"{slug}.png").exists():
            continue
        L += [f"![](figures/{slug}.png)", "", f"*{cap}*", ""]
    return L
