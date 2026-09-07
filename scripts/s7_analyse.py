"""S7 step 2 — everything the tree says, written out as tables.

The tree is read once here and turned into the committed TSVs that
`s7_report.py` and `s7_figure.py` render (D13: no report and no figure
recomputes anything from the Newick).

  rooted.nwk                the ML tree rooted on the RyR outgroup
  monophyly.tsv             every census group, monophyletic or not
  membership_audit.tsv      the tree-corrected paralog membership (R1-R3)
  claim_nodes.tsv           SH-aLRT / UFBoot at each node a claim rests on
  claim_members.tsv         the tip set behind each of those claims
  sister_ml.tsv             each paralog's sister group on the rooted tree
  support_summary.tsv       how much of the tree is actually resolved
  cyclostome_placement.tsv  S6's open question, asked of the tree
  duplicate_pairs.tsv       S6's 3R stress test: same species, same paralog
  naming_conflicts.tsv      tips whose census label the tree contradicts

Run:  python3 scripts/s7_analyse.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.s6_lib import write_tsv  # noqa: E402
from scripts.s7_lib import (  # noqa: E402
    PHYLO_DIR, Node, find_clade, find_clade_rooted, group_labels,
    load_groups, parse_newick,
    pure_clades, reroot, sister_of, support_of, to_newick,
)
from scripts.s7_constraints import (  # noqa: E402
    MIN_ALRT, MIN_UFBOOT, PARALOGS, cores, corrected_membership, home_of,
)

TREEFILE = PHYLO_DIR / "itpr_ml.treefile"
ROOTED = PHYLO_DIR / "rooted.nwk"
STATS = PHYLO_DIR / "tree_stats.json"

GROUP_ORDER = ["ITPR1", "ITPR2", "ITPR3", "vertebrate_basal",
               "invert_metazoa", "plant", "protist", "fungi", "RYR"]
OUTGROUP = "RYR"


def _sup(name: str) -> tuple[str, str]:
    a, u = support_of(name)
    return ("" if a is None else f"{a:g}", "" if u is None else f"{u:g}")


def main() -> int:
    groups = load_groups()
    tree = parse_newick(TREEFILE.read_text())
    all_leaves = tree.leaf_names()
    core = cores(groups, tree)

    # ---------------------------------------------------- 1. monophyly
    rows = []
    for g in GROUP_ORDER:
        mem = group_labels(groups, g) & all_leaves
        if not mem:
            continue
        mono, sup = find_clade(tree, mem)
        parts = pure_clades(tree, mem, all_leaves)
        outside = sorted(mem - parts[0]) if parts else sorted(mem)
        a, u = _sup(sup)
        rows.append({
            "group": g, "n": len(mem),
            "monophyletic": "yes" if mono else "no",
            "alrt": a, "ufboot": u,
            "largest_pure_clade": len(parts[0]) if parts else 0,
            "n_outside": 0 if mono else len(outside),
            "outside": "" if mono else ";".join(outside[:8]),
        })
    write_tsv(PHYLO_DIR / "monophyly.tsv",
              ["group", "n", "monophyletic", "alrt", "ufboot",
               "largest_pure_clade", "n_outside", "outside"], rows)

    # ------------------------------------- 2. tree-corrected membership
    corrected, audit = corrected_membership(groups, tree)
    write_tsv(PHYLO_DIR / "membership_audit.tsv",
              ["label", "census_group", "rule", "assigned_to", "clade_size",
               "neighbours", "alrt", "ufboot", "why"], audit)

    conflicts = [a for a in audit if a["rule"] in ("reassigned",
                                                   "unconstrained")]
    write_tsv(PHYLO_DIR / "naming_conflicts.tsv",
              ["label", "census_group", "rule", "assigned_to", "clade_size",
               "neighbours", "alrt", "ufboot", "why"], conflicts)

    # corrected monophyly of the three paralogs
    crows = []
    for g in PARALOGS:
        mono, sup = find_clade(tree, corrected[g])
        a, u = _sup(sup)
        crows.append({"group": g, "n": len(corrected[g]),
                      "monophyletic": "yes" if mono else "no",
                      "alrt": a, "ufboot": u})
    write_tsv(PHYLO_DIR / "monophyly_corrected.tsv",
              ["group", "n", "monophyletic", "alrt", "ufboot"], crows)

    # ------------------------------------------------------- 3. rooting
    out_lbls = group_labels(groups, OUTGROUP) & all_leaves
    rooted = reroot(parse_newick(TREEFILE.read_text()), out_lbls)
    ROOTED.write_text(to_newick(rooted) + "\n")
    parts_out = pure_clades(tree, out_lbls, all_leaves)
    root_row = [{
        "outgroup": OUTGROUP, "n": len(out_lbls),
        "monophyletic": "yes" if find_clade(tree, out_lbls)[0] else "no",
        "largest_pure_clade": len(parts_out[0]) if parts_out else 0,
        "rooted_file": ROOTED.name,
    }]
    write_tsv(PHYLO_DIR / "rooting.tsv",
              ["outgroup", "n", "monophyletic", "largest_pure_clade",
               "rooted_file"], root_row)

    # ------------- 3b. what each paralog clade actually contains
    # "Not monophyletic" is too blunt a verdict to act on, and it has to
    # be asked of the *rooted* tree: IQ-TREE writes an arbitrary root,
    # and a group straddling it has the whole tree as its MRCA. A clade
    # broken by another paralog's tips means the labels are wrong or the
    # tree is; a clade broken only by vertebrate tips carrying *no*
    # paralog label is a different statement entirely — it says the tree
    # places a locus S6 could not label, which is a result, not a
    # problem. So each paralog's **extended clade** is grown from its
    # core while no other paralog's labelled tip is swallowed, and what
    # it picked up is reported.
    rooted_leaves = rooted.leaf_names()
    labelled = {g: group_labels(groups, g) & rooted_leaves for g in PARALOGS}
    rparent: dict[int, Node] = {}
    for n in rooted.walk():
        for c in n.children:
            rparent[id(c)] = n
    extended: dict[str, Node] = {}
    mrows = []
    for g in PARALOGS:
        others = set().union(*[labelled[h] for h in PARALOGS if h != g])
        node = next((n for n in rooted.walk()
                     if not n.is_leaf and n.leaf_names() == corrected[g]),
                    None)
        if node is None:
            continue
        while True:
            up = rparent.get(id(node))
            if up is None or (up.leaf_names() & others):
                break
            node = up
        extended[g] = node
        ls = node.leaf_names()
        added = sorted(ls - labelled[g])
        by_group: dict[str, int] = {}
        for lab in added:
            k = groups.get(lab, {}).get("group", "?")
            by_group[k] = by_group.get(k, 0) + 1
        missing = sorted(labelled[g] - ls)
        a, u = _sup(node.name)
        if not added and not missing:
            verdict = "exactly the labelled set"
        elif missing:
            verdict = "some labelled tips fall outside"
        else:
            verdict = "labelled set plus unlabelled vertebrate tips"
        mrows.append({
            "paralog": g, "n_labelled": len(labelled[g]),
            "n_core": len(corrected[g]), "n_extended": len(ls),
            "alrt": a, "ufboot": u,
            "n_added": len(added), "added_groups":
                ";".join(f"{k}:{v}" for k, v in sorted(by_group.items())),
            "n_labelled_outside": len(missing),
            "verdict": verdict,
            "added": ";".join(added[:8]),
            "labelled_outside": ";".join(missing[:8]),
        })
    write_tsv(PHYLO_DIR / "paralog_clades.tsv",
              ["paralog", "n_labelled", "n_core", "n_extended", "alrt",
               "ufboot", "n_added", "added_groups", "n_labelled_outside",
               "verdict", "added", "labelled_outside"], mrows)

    # -------------------------------------------------- 4. claim nodes
    # The pair claims are asked of the **extended** clades: a pair test
    # on the bare cores would be answered "no" by the chondrichthyan tips
    # sitting between them, which is not what the sister question means.
    ext = {g: set(extended[g].leaf_names()) if g in extended
           else set(corrected[g]) for g in PARALOGS}
    claims: list[tuple[str, set[str]]] = []
    for g in PARALOGS:
        claims.append((f"{g} core clade", set(corrected[g])))
    for g in PARALOGS:
        claims.append((f"{g} clade incl. unlabelled tips", ext[g]))
    for a, b in (("ITPR1", "ITPR2"), ("ITPR1", "ITPR3"), ("ITPR2", "ITPR3")):
        claims.append((f"{a} + {b}", ext[a] | ext[b]))
    claims.append(("all three paralogs",
                   set().union(*[ext[g] for g in PARALOGS])))
    verts = group_labels(groups, *PARALOGS, "vertebrate_basal") & all_leaves
    claims.append(("vertebrate ITPRs (incl. unassigned)", verts))
    claims.append(("ITPR family (all non-RyR)", set(all_leaves) - out_lbls))
    claims.append(("RyR outgroup", set(out_lbls)))

    crow = []
    for desc, mem in claims:
        mem = set(mem) & rooted_leaves
        # the rooted test: "ITPR1 and ITPR2 are sisters" is a claim about
        # the root, and the unrooted test would answer a different question
        mono, sup = find_clade_rooted(rooted, mem)
        a, u = _sup(sup)
        crow.append({"claim": desc, "n_tips": len(mem),
                     "is_clade": "yes" if mono else "no",
                     "alrt": a, "ufboot": u,
                     "well_supported": ("yes" if (a and u
                                                  and float(a) >= MIN_ALRT
                                                  and float(u) >= MIN_UFBOOT)
                                        else "no")})
    # The membership behind every claim, one row per tip. `s7_bnni.py`
    # re-asks these same clade questions of the --bnni tree, and it must
    # ask about the *same sets* — recovering them from a clade size would
    # let the two modules disagree about what "ITPR2 + ITPR3" means.
    write_tsv(PHYLO_DIR / "claim_members.tsv", ["claim", "label"],
              [{"claim": desc, "label": lab}
               for desc, mem in claims
               for lab in sorted(set(mem) & rooted_leaves)])

    write_tsv(PHYLO_DIR / "claim_nodes.tsv",
              ["claim", "n_tips", "is_clade", "alrt", "ufboot",
               "well_supported"], crow)

    # ------------------------------------------- 5. the ML sister answer
    srows = []
    for g in PARALOGS:
        sis = sister_of(rooted, set(corrected[g]))
        comp = {}
        for h in GROUP_ORDER:
            k = len(sis & (group_labels(groups, h) & rooted_leaves))
            if k:
                comp[h] = k
        top = max(comp, key=comp.get) if comp else ""
        srows.append({
            "paralog": g, "is_clade": "yes" if sis else "no",
            "sister_n": len(sis),
            "sister_composition": ";".join(f"{k}:{v}" for k, v in
                                           sorted(comp.items(),
                                                  key=lambda t: -t[1])),
            "sister_majority": top,
        })
    write_tsv(PHYLO_DIR / "sister_ml.tsv",
              ["paralog", "is_clade", "sister_n", "sister_composition",
               "sister_majority"], srows)

    # ------------------------------------------------ 6. support summary
    alrts, ufs = [], []
    for n in tree.walk():
        if n.is_leaf or n is tree:
            continue
        a, u = support_of(n.name)
        if a is not None:
            alrts.append(a)
        if u is not None:
            ufs.append(u)
    both = sum(1 for n in tree.walk()
               if not n.is_leaf and n is not tree
               and (support_of(n.name)[0] or 0) >= MIN_ALRT
               and (support_of(n.name)[1] or 0) >= MIN_UFBOOT)
    tot = len(ufs)
    write_tsv(PHYLO_DIR / "support_summary.tsv",
              ["statistic", "value"], [
                  {"statistic": "internal nodes", "value": tot},
                  {"statistic": "UFBoot >= 95", "value":
                   sum(1 for u in ufs if u >= 95)},
                  {"statistic": "SH-aLRT >= 80", "value":
                   sum(1 for a in alrts if a >= 80)},
                  {"statistic": "both thresholds", "value": both},
                  {"statistic": "pct both", "value":
                   f"{100.0 * both / tot:.1f}" if tot else ""},
                  {"statistic": "median UFBoot", "value":
                   f"{sorted(ufs)[len(ufs) // 2]:g}" if ufs else ""},
                  {"statistic": "median SH-aLRT", "value":
                   f"{sorted(alrts)[len(alrts) // 2]:g}" if alrts else ""},
              ])

    # -------------------------------------------- 7. cyclostome question
    cyc = [l for l, r in groups.items()
           if r.get("band") == "cyclostomata" and l in all_leaves]
    crows2 = []
    for lab in sorted(cyc):
        h = home_of(tree, lab, core)
        crows2.append({"label": lab,
                       "species": groups[lab].get("species", ""),
                       "home": h["home"], "clade_size": h["clade_size"],
                       "n_home": h["n_home"],
                       "alrt": "" if h["alrt"] is None else f"{h['alrt']:g}",
                       "ufboot": "" if h["ufboot"] is None else f"{h['ufboot']:g}",
                       "composition": h["composition"]})
    write_tsv(PHYLO_DIR / "cyclostome_placement.tsv",
              ["label", "species", "home", "clade_size", "n_home",
               "alrt", "ufboot", "composition"], crows2)
    cyc_mono, cyc_sup = find_clade(tree, set(cyc)) if cyc else (False, "")
    by_sp: dict[str, list[str]] = {}
    for lab in cyc:
        by_sp.setdefault(groups[lab].get("species", "?"), []).append(lab)
    crows3 = [{"set": "all cyclostome loci", "n": len(cyc),
               "monophyletic": "yes" if cyc_mono else "no",
               "alrt": _sup(cyc_sup)[0], "ufboot": _sup(cyc_sup)[1],
               "members": ";".join(sorted(cyc))[:200]}]
    for sp, labs in sorted(by_sp.items()):
        m, s = find_clade(tree, set(labs))
        crows3.append({"set": sp, "n": len(labs),
                       "monophyletic": "yes" if m else "no",
                       "alrt": _sup(s)[0], "ufboot": _sup(s)[1],
                       "members": ";".join(sorted(labs))[:200]})
    # "Not one clade" is not the same as "scattered". The maximal pure
    # cyclostome pieces say which it is, and a piece spanning both
    # species is a duplication older than the hagfish/lamprey split.
    for piece in pure_clades(tree, set(cyc), all_leaves):
        if len(piece) < 2:
            continue
        m, s = find_clade(tree, set(piece))
        spp = {groups[l].get("species", "?").split(" (")[0] for l in piece}
        crows3.append({
            "set": f"cyclostome-only clade ({len(spp)} species)",
            "n": len(piece), "monophyletic": "yes" if m else "no",
            "alrt": _sup(s)[0], "ufboot": _sup(s)[1],
            "members": ";".join(sorted(piece))[:200]})
    write_tsv(PHYLO_DIR / "cyclostome_monophyly.tsv",
              ["set", "n", "monophyletic", "alrt", "ufboot", "members"],
              crows3)

    # --------------------------------------------- 8. the 3R stress test
    pairs: dict[tuple[str, str], list[str]] = {}
    for lab, r in groups.items():
        if lab not in all_leaves or r["group"] not in PARALOGS:
            continue
        pairs.setdefault((r["group"], r.get("species", "?")), []).append(lab)
    prows = []
    for (g, sp), labs in sorted(pairs.items()):
        if len(labs) < 2:
            continue
        m, s = find_clade(tree, set(labs))
        a, u = _sup(s)
        # A binary answer cannot tell one rogue copy from a real
        # scattering, and those call for opposite conclusions. The
        # largest pure piece and the MRCA's composition say which.
        parts = pure_clades(tree, set(labs), all_leaves)
        biggest = parts[0] if parts else frozenset()
        bm, bs = find_clade(tree, set(biggest)) if len(biggest) > 1 \
            else (False, "")
        mrca = None
        for n in rooted.walk():
            if n.is_leaf or not set(labs) <= n.leaf_names():
                continue
            if mrca is None or len(n.leaf_names()) < len(mrca.leaf_names()):
                mrca = n
        others: dict[str, int] = {}
        if mrca is not None:
            for lab in mrca.leaf_names() - set(labs):
                k = groups.get(lab, {}).get("group", "?")
                others[k] = others.get(k, 0) + 1
        ba, bu = _sup(bs)
        prows.append({"paralog": g, "species": sp, "n": len(labs),
                      "sisters": "yes" if m else "no", "alrt": a,
                      "ufboot": u,
                      "largest_pure_subset": len(biggest),
                      "subset_alrt": ba, "subset_ufboot": bu,
                      "mrca_size": len(mrca.leaf_names()) if mrca else "",
                      "mrca_intruders": ";".join(f"{k}:{v}" for k, v in
                                                 sorted(others.items())),
                      "labels": ";".join(sorted(labs))})
    write_tsv(PHYLO_DIR / "duplicate_pairs.tsv",
              ["paralog", "species", "n", "sisters", "alrt", "ufboot",
               "largest_pure_subset", "subset_alrt", "subset_ufboot",
               "mrca_size", "mrca_intruders", "labels"], prows)

    # ------------------------------------------------------- 9. summary
    summary = {
        "n_tips": len(all_leaves),
        "n_reassigned": sum(1 for a in audit if a["rule"] == "reassigned"),
        "n_unconstrained": sum(1 for a in audit
                               if a["rule"] == "unconstrained"),
        "paralogs_monophyletic_raw": sum(
            1 for r in rows if r["group"] in PARALOGS
            and r["monophyletic"] == "yes"),
        "paralogs_monophyletic_corrected": sum(
            1 for r in crows if r["monophyletic"] == "yes"),
        "paralogs_clean_extended_clade": sum(
            1 for r in mrows if r["n_labelled_outside"] == 0),
        "unlabelled_tips_placed_inside_a_paralog": sum(
            r["n_added"] for r in mrows),
        "outgroup_monophyletic": root_row[0]["monophyletic"],
        "pct_nodes_both_thresholds": (f"{100.0 * both / tot:.1f}"
                                      if tot else ""),
    }
    if STATS.exists():
        d = json.loads(STATS.read_text())
        d["analysis"] = summary
        STATS.write_text(json.dumps(d, indent=1) + "\n")
    print(json.dumps(summary, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
