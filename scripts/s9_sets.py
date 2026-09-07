"""S9 — which tips belong to which selection set, and why the tree decides.

codeml's branch models mark a *node*: a "foreground" that is not a clade
silently becomes something else entirely, and a two-ratio test on it is
meaningless rather than wrong. So the three paralog sets S9 estimates ω on
are not the census labels — they are S7's **extended paralog clades**, the
largest clade grown from each paralog's tree-corrected core while no other
paralog's labelled tip is swallowed.

Two consequences worth stating, because they are the reason this module
exists rather than a hand-written table:

  * A tip whose census label the tree declined to honour (S7 recorded five
    as `unconstrained`) is placed by the tree or by nothing. It is never
    put in a foreground because a database calls it ITPR2.
  * The unlabelled `vertebrate_basal` tips nested *inside* a paralog clade
    join that set. That is D30 read forwards: outside the vertebrates a
    paralog label means nothing, but inside a well-supported vertebrate
    paralog clade the tree has just supplied one.

The clades are re-derived here from `rooted.nwk` with S7's own rule and
then **cross-checked against `results/phylogeny/paralog_clades.tsv`**; a
disagreement in size or membership is a hard failure, because it would
mean S9 and S7 are describing different clades under the same name.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_constraints import PARALOGS, cores  # noqa: E402
from scripts.s7_lib import (Node, group_labels, load_groups,  # noqa: E402
                            parse_newick)
PHYLO_DIR = ROOT / "results" / "phylogeny"
ROOTED_NWK = PHYLO_DIR / "rooted.nwk"
PARALOG_CLADES_TSV = PHYLO_DIR / "paralog_clades.tsv"


def extended_clades(tree: Node) -> dict[str, tuple[frozenset[str], str]]:
    """S7's extended paralog clade per paralog: (leaf set, support label).

    Same rule as `s7_analyse.py` §3b — grow from the corrected core while
    the parent swallows no *other* paralog's labelled tip.
    """
    groups = load_groups()
    leaves = tree.leaf_names()
    labelled = {g: group_labels(groups, g) & leaves for g in PARALOGS}
    core = cores(groups, tree)
    parent: dict[int, Node] = {}
    for n in tree.walk():
        for c in n.children:
            parent[id(c)] = n
    out: dict[str, tuple[frozenset[str], str]] = {}
    for g in PARALOGS:
        others = set().union(*[labelled[h] for h in PARALOGS if h != g])
        node = next((n for n in tree.walk()
                     if not n.is_leaf and n.leaf_names() == core[g]), None)
        if node is None:
            continue
        while True:
            up = parent.get(id(node))
            if up is None or (up.leaf_names() & others):
                break
            node = up
        out[g] = (node.leaf_names(), node.name)
    return out


def check_against_s7(ext: dict[str, tuple[frozenset[str], str]]) -> list[str]:
    """Compare the re-derived clades with S7's committed table.

    Sizes are compared for all three; membership only for the tips S7's
    `added` column actually lists, since that column is truncated at eight
    entries and a truncated list is evidence of what *is* in the clade, not
    of what is not.
    """
    if not PARALOG_CLADES_TSV.exists():
        return [f"missing {PARALOG_CLADES_TSV}"]
    problems: list[str] = []
    groups = load_groups()
    with open(PARALOG_CLADES_TSV) as fh:
        rows = {r["paralog"]: r for r in csv.DictReader(fh, delimiter="\t")}
    for g in PARALOGS:
        if g not in rows:
            problems.append(f"{g}: no row in paralog_clades.tsv")
            continue
        if g not in ext:
            problems.append(f"{g}: no extended clade re-derived from the tree")
            continue
        members, _ = ext[g]
        want = int(rows[g]["n_extended"])
        if len(members) != want:
            problems.append(f"{g}: extended clade is {len(members)} tips here, "
                            f"{want} in paralog_clades.tsv")
        listed = [x for x in rows[g]["added"].split(";") if x]
        for lab in listed:
            if lab not in members:
                problems.append(f"{g}: S7 lists {lab} in the clade, "
                                "this derivation does not")
        labelled = group_labels(groups, g) & set(members)
        if len(members) - len(labelled) != len(listed) and len(listed) < 8:
            problems.append(f"{g}: {len(members) - len(labelled)} unlabelled "
                            f"tips in the clade, S7 lists {len(listed)}")
    return problems


def paralog_of(tips: list[str]) -> tuple[dict[str, str], dict[str, str]]:
    """(tip -> paralog set) for the tips in a paralog clade, and the reason.

    Tips in no extended clade are absent from the mapping: they stay in the
    whole-tree analyses as background and are in no per-paralog set. That
    is the honest handling of a locus the tree could not place — S7 found
    six of them, and dropping them from the tree entirely would change the
    branch lengths every other estimate is made on.
    """
    tree = parse_newick(ROOTED_NWK.read_text())
    ext = extended_clades(tree)
    problems = check_against_s7(ext)
    if problems:
        raise SystemExit("S9/S7 paralog clades disagree:\n  "
                         + "\n  ".join(problems))
    groups = load_groups()
    keep = set(tips)
    assign: dict[str, str] = {}
    why: dict[str, str] = {}
    for g in PARALOGS:
        if g not in ext:
            continue
        for lab in ext[g][0] & keep:
            census = groups.get(lab, {}).get("group", "?")
            assign[lab] = g
            why[lab] = ("census label, in the tree's clade" if census == g
                        else f"tree places a {census} tip in the {g} clade")
    return assign, why


def main() -> None:
    tree = parse_newick(ROOTED_NWK.read_text())
    ext = extended_clades(tree)
    problems = check_against_s7(ext)
    for g in PARALOGS:
        members, sup = ext.get(g, (frozenset(), ""))
        print(f"{g}: {len(members)} tips  support={sup or '-'}")
    print("cross-check vs S7:", "OK" if not problems else "FAILED")
    for p in problems:
        print("  ", p)
    raise SystemExit(1 if problems else 0)


if __name__ == "__main__":
    main()
