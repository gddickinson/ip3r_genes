"""S7 — the AU test's three hypotheses, and the membership they rest on.

The sister question is about the *three vertebrate paralog clades*, so
each hypothesis is one split of the unrooted tree, with the ryanodine
receptors present to make the split a rooted statement:

  H1_12  ((ITPR1,ITPR2),ITPR3)   split {ITPR1 ∪ ITPR2} vs everything else
  H2_13  ((ITPR1,ITPR3),ITPR2)
  H3_23  ((ITPR2,ITPR3),ITPR1)

**Membership is tree-corrected, and the correction is derived rather
than listed.** A constraint set built from census labels would ask the
AU test about a topology the data rejects for a reason that has nothing
to do with the sister question — a single mis-annotated tip forced back
into the clade its name claims costs likelihood everywhere, and all
three hypotheses then fail together. (That is exactly what happened on
the PIEZO project's first pass, and the fix there was a hand-written
relabelling table. Here the tree writes the table itself.)

Three rules, each a positive test, each recorded in
`results/phylogeny/membership_audit.tsv` with the number it fired on:

  R1 `core`          the tip is inside its own group's largest pure
                     clade — its label and the tree agree.
  R2 `reassigned`    the tip is outside that core, and the **smallest**
                     clade containing it that holds at least
                     `MIN_NEIGHBOURS` other tips is made up *entirely* of
                     one other paralog's core — the tip is nested inside
                     that paralog — at a node clearing both of IQ-TREE's
                     own thresholds (SH-aLRT ≥ 80 **and** UFBoot ≥ 95).
                     The rule stops at that first qualifying clade and
                     does not keep climbing: a rule that climbs until
                     some group happens to dominate will hand a tip to
                     whichever paralog is largest three nodes up, which
                     is a fact about clade sizes and not about the tip.
                     (The self-test caught exactly that: `s7_test_tree.py`
                     T6/T7.)
  R3 `unconstrained` anything else — outside its core with no
                     well-supported home. Left out of all three
                     constraint sets rather than assumed either way,
                     which is the honest treatment of a tip the tree
                     cannot place and keeps it from deciding the test.

Nothing here reads a sequence, a length or a species.
"""

from __future__ import annotations

from scripts.s7_lib import Node, pure_clades, support_of

#: The three vertebrate paralogs, in the order their names imply.
PARALOGS = ("ITPR1", "ITPR2", "ITPR3")

#: Support a placement must clear before it may overrule a label
#: (IQ-TREE's own joint recommendation for a trustworthy clade).
MIN_ALRT = 80.0
MIN_UFBOOT = 95.0

#: How many members of another core group must share the smallest
#: containing clade before that clade counts as a home.
MIN_NEIGHBOURS = 3

HYPOTHESES = {
    "ML": "unconstrained ML tree",
    "H1_12": "ITPR1 + ITPR2 sisters, ITPR3 outside",
    "H2_13": "ITPR1 + ITPR3 sisters, ITPR2 outside",
    "H3_23": "ITPR2 + ITPR3 sisters, ITPR1 outside",
}

#: hypothesis -> the pair of paralogs it makes sisters
PAIRS = {"H1_12": ("ITPR1", "ITPR2"),
         "H2_13": ("ITPR1", "ITPR3"),
         "H3_23": ("ITPR2", "ITPR3")}


def cores(groups: dict[str, dict], tree: Node) -> dict[str, frozenset[str]]:
    """Each paralog's largest pure clade on the tree — R1's `core`."""
    all_leaves = tree.leaf_names()
    out: dict[str, frozenset[str]] = {}
    for g in PARALOGS:
        mem = {l for l, r in groups.items() if r["group"] == g} & all_leaves
        parts = pure_clades(tree, mem, all_leaves)
        out[g] = parts[0] if parts else frozenset()
    return out


def home_of(tree: Node, label: str, core: dict[str, frozenset[str]],
            min_k: int = None) -> dict:
    """Where the tree nests `label`, or nothing.

    Walks the clades containing `label` from smallest up, stops at the
    first that holds at least `min_k` *other* tips, and reports a home
    only if every one of those tips belongs to one paralog's core. An
    ambiguous neighbourhood returns `home=""` with the composition that
    made it ambiguous, so the audit can say why the tip was left free.
    """
    k = MIN_NEIGHBOURS if min_k is None else min_k
    best: dict = {"home": "", "clade_size": "", "n_home": "",
                  "alrt": "", "ufboot": "", "composition": ""}
    chosen = None
    for n in tree.walk():
        if n.is_leaf:
            continue
        ls = n.leaf_names()
        if label not in ls or len(ls) - 1 < k:
            continue
        if chosen is None or len(ls) < len(chosen.leaf_names()):
            chosen = n
    if chosen is None:
        return best
    others = chosen.leaf_names() - {label}
    counts = {g: len(others & c) for g, c in core.items()}
    inside = [g for g, m in counts.items() if m and others <= core[g]]
    alrt, ufb = support_of(chosen.name)
    # The residue is always reported. Without it a composition reading
    # "ITPR3:16" beside an `unconstrained` verdict looks like a rule that
    # failed to fire, when what actually happened is that the
    # neighbourhood also held tips belonging to no paralog core.
    residue = len(others) - sum(counts.values())
    parts = [f"{g}:{m}" for g, m in sorted(counts.items()) if m]
    if residue:
        parts.append(f"unlabelled:{residue}")
    best.update({
        "clade_size": len(chosen.leaf_names()),
        "alrt": alrt, "ufboot": ufb,
        "composition": ";".join(parts) or f"other:{len(others)}",
    })
    if len(inside) == 1:
        best["home"] = inside[0]
        best["n_home"] = counts[inside[0]]
    return best


def corrected_membership(groups: dict[str, dict], tree: Node,
                         ) -> tuple[dict[str, set[str]], list[dict]]:
    """The three constraint sets, plus the audit row behind every tip."""
    all_leaves = tree.leaf_names()
    core = cores(groups, tree)
    out: dict[str, set[str]] = {g: set(core[g]) for g in PARALOGS}
    audit: list[dict] = []

    for g in PARALOGS:
        labelled = {l for l, r in groups.items()
                    if r["group"] == g} & all_leaves
        for lab in sorted(labelled):
            if lab in core[g]:
                audit.append({"label": lab, "census_group": g,
                              "rule": "core", "assigned_to": g,
                              "clade_size": len(core[g]),
                              "neighbours": len(core[g]) - 1,
                              "alrt": "", "ufboot": "",
                              "why": f"inside the largest pure {g} clade"})
                continue
            h = home_of(tree, lab, core)
            home, alrt, ufb = h["home"], h["alrt"], h["ufboot"]
            strong = (alrt is not None and ufb is not None
                      and alrt >= MIN_ALRT and ufb >= MIN_UFBOOT)
            if home and home != g and strong:
                out[home].add(lab)
                audit.append({
                    "label": lab, "census_group": g, "rule": "reassigned",
                    "assigned_to": home, "clade_size": h["clade_size"],
                    "neighbours": h["n_home"], "alrt": alrt, "ufboot": ufb,
                    "why": (f"nested inside {home} — its nearest "
                            f"{h['n_home']}-tip neighbourhood is entirely "
                            f"{home}, at SH-aLRT {alrt}/UFBoot {ufb}")})
            else:
                if not home:
                    why = ("nearest neighbourhood is mixed ("
                           + (h["composition"] or "empty") + ")")
                elif home == g:
                    why = (f"nested inside its own group {g} but outside "
                           f"that group's largest pure clade")
                else:
                    why = (f"nested inside {home} but support "
                           f"{alrt}/{ufb} is under {MIN_ALRT}/{MIN_UFBOOT}")
                audit.append({
                    "label": lab, "census_group": g,
                    "rule": "unconstrained", "assigned_to": "",
                    "clade_size": h["clade_size"],
                    "neighbours": h["n_home"],
                    "alrt": alrt if alrt is not None else "",
                    "ufboot": ufb if ufb is not None else "",
                    "why": why})

    # A tip whose census group is not a paralog may still sit inside a
    # paralog core (it is, after all, the tree's own clade): it entered
    # `out` with the core and is recorded here rather than silently.
    for g in PARALOGS:
        for lab in sorted(core[g]):
            if groups.get(lab, {}).get("group") != g:
                audit.append({
                    "label": lab, "census_group": groups.get(lab, {}).get(
                        "group", "?"), "rule": "core_foreign",
                    "assigned_to": g, "clade_size": len(core[g]),
                    "neighbours": len(core[g]) - 1, "alrt": "", "ufboot": "",
                    "why": f"unlabelled tip inside the {g} core clade"})
    return out, audit


def au_constraints(core: dict[str, set[str]], outgroup: set[str],
                   all_leaves: frozenset[str]) -> dict[str, str]:
    """One `-g` constraint newick per hypothesis.

    **A constraint names only the taxa its hypothesis is about.**
    IQ-TREE's `-g` constraint may cover a subset of the alignment, and
    the taxa it leaves out are exactly the ones the search places
    freely — including *inside* a constrained clade. A taxon that is
    listed, even in a top-level polytomy, is thereby fixed **outside**
    every group the constraint declares.

    The first version of this function listed all 134 tips, on the
    reading that a top-level polytomy leaves them free. It does not: it
    pinned the ~14 unlabelled vertebrate tips the tree nests inside the
    paralog clades out of those clades, identically in all three
    hypotheses. The AU test then rejected every hypothesis at ΔlogL
    ≈ 1,500 — *including* H3_23, which the ML tree itself holds at
    SH-aLRT 100 / UFBoot 100. That is the PIEZO failure mode this
    module was written to avoid, arriving through the constraint's
    taxon set rather than through its labels: a test in which the true
    hypothesis is rejected along with the false ones is not measuring
    the sister question.

    So each constraint carries the three paralog cores and the
    outgroup, and nothing else. The outgroup is what makes the split a
    **rooted** statement — on three groups alone, `((A,B),C)` is an
    unrooted trifurcation and asserts only that each group is a clade,
    which every hypothesis shares.
    """
    assigned: set[str] = set()
    for v in core.values():
        assigned |= set(v)
    og = sorted(set(outgroup) & set(all_leaves))
    if not og:
        raise ValueError("no outgroup tips — the split would be unrooted")
    if set(og) & assigned:
        raise ValueError("outgroup is inside a paralog constraint set")
    free = all_leaves - assigned - set(og)
    if not free:
        raise ValueError("constraint would fix every tip — no free taxa")
    grp = {g: "(" + ",".join(sorted(core[g])) + ")" for g in PARALOGS}
    og_grp = "(" + ",".join(og) + ")"

    out: dict[str, str] = {}
    for name, (a, b) in PAIRS.items():
        c = next(g for g in PARALOGS if g not in (a, b))
        out[name] = f"((({grp[a]},{grp[b]}),{grp[c]}),{og_grp});"
    return out
