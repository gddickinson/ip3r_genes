"""s15b_dollo.py — Dollo parsimony, and the bound the tree's resolution
puts on any count it returns.

Dollo's assumption is the one this character can carry: a complex gene is
gained once and may be lost many times.  So the single gain is placed at
the MRCA of every tip that has the gene, and a loss is an edge below it
whose whole subtree has lost it.

Two properties matter more than the algorithm.

**The count must be reachable.**  On S15a's matrix it is zero, and a
routine that returns zero because it cannot return anything else is not a
measurement.  `s15b_test_counts.py` puts a constructed loss on a known
edge and requires this module to find it, and puts two on sister edges and
requires them merged into their parent.

**A count of *independent* losses is bounded by the tree's resolution.**
Under a polytomy of degree m, k sibling losses might be k separate events
or one event on a branch the polytomy does not resolve, so the honest
answer is an interval: `n_max` counts every loss edge, `n_min` counts one
per parent that carries any.  They differ only where a polytomy carries
more than one loss, which is exactly the case where a single number would
be an artefact of the tree rather than a result.

**Where the gain goes is a decision, and it is made per character.**
Unpinned, Dollo places the single gain at the MRCA of the present tips,
and a clade entirely without the gene is then *ancestrally* absent rather
than a loss — which is the right reading for ITPR2 and ITPR3, whose
duplications S13 places on the gnathostome stem, so no cyclostome ever had
them.  It is the wrong reading for the family as a whole: S20 and S23
found IP3 receptors across the eukaryotes, so the family was present at
the root of any vertebrate tree and a vertebrate clade without one has
lost it.  `pin_gain_at_root=True` says so explicitly, and every row of
every S15b table records which rule it was counted under.

`?` (undecided) tips are neither present nor absent.  A subtree carrying
only `?` supports no loss, and a `?` tip inside an otherwise-lost subtree
does not rescue it — the loss is placed and the tip is recorded as
uninformative for it.
"""

from __future__ import annotations

import s15b_lib as lib


def _subtree_counts(root: lib.Node, char: dict[str, int | None]):
    """Per node: (#present tips, #absent tips, #unknown tips) below it."""
    counts: dict[int, tuple[int, int, int]] = {}
    for n in lib.postorder(root):
        if n.is_tip:
            v = char.get(n.name)
            counts[id(n)] = ((1, 0, 0) if v == 1 else
                             (0, 1, 0) if v == 0 else (0, 0, 1))
        else:
            p = a = u = 0
            for c in n.children:
                cp, ca, cu = counts[id(c)]
                p, a, u = p + cp, a + ca, u + cu
            counts[id(n)] = (p, a, u)
    return counts


def mrca_of_present(root: lib.Node, char: dict[str, int | None]):
    """The deepest node whose subtree holds every present tip."""
    counts = _subtree_counts(root, char)
    total = counts[id(root)][0]
    if total == 0:
        return None, counts
    node = root
    while True:
        nxt = next((c for c in node.children if counts[id(c)][0] == total),
                   None)
        if nxt is None:
            return node, counts
        node = nxt


def dollo(root: lib.Node, char: dict[str, int | None],
          pin_gain_at_root: bool = False) -> dict:
    """Losses under Dollo parsimony, with the resolution bound.

    Returns the loss edges (each named by the node below it), the naive
    count `n_max`, the resolution-collapsed minimum `n_min`, and the gain
    node the whole reconstruction hangs from.

    `pin_gain_at_root` places the gain above the root instead of at the
    MRCA of the present tips.  Use it when an outside result establishes
    that the character was present at the root — for this project, S20's
    eukaryote-wide range of the family — because otherwise a clade with no
    copy is scored as never having had one.
    """
    mrca, counts = mrca_of_present(root, char)
    gain = root if pin_gain_at_root else mrca
    n_tips = len(lib.tips(root))
    scored = sum(1 for t in lib.tips(root) if char.get(t.name) is not None)
    if gain is None or (pin_gain_at_root and mrca is None):
        return dict(losses=[], n_max=0, n_min=0, gain_node="",
                    gain_rule="pinned at root" if pin_gain_at_root
                    else "MRCA of present tips",
                    n_tips=n_tips, n_scored=scored, n_present=0,
                    n_absent=sum(1 for t in lib.tips(root)
                                 if char.get(t.name) == 0),
                    n_unknown=n_tips - scored,
                    note="no tip carries the character; no gain to place")

    losses: list[dict] = []
    stack = list(gain.children)
    while stack:
        n = stack.pop()
        p, a, u = counts[id(n)]
        if p == 0 and a > 0:
            parent = n.parent
            losses.append(dict(
                node=n.name or f"<unnamed:{a + u} tips>",
                parent=parent.name if parent is not None else "",
                parent_degree=len(parent.children) if parent else 0,
                n_tips_lost=a, n_tips_unknown=u,
                tips=";".join(sorted(t.name for t in lib.tips(n)))))
            continue
        stack.extend(n.children)

    by_parent: dict[str, int] = {}
    for ls in losses:
        by_parent[ls["parent"]] = by_parent.get(ls["parent"], 0) + 1
    n_min = len(by_parent)
    for ls in losses:
        ls["siblings_lost_under_parent"] = by_parent[ls["parent"]]
        ls["parent_is_polytomy"] = int(ls["parent_degree"] >= 3)
    pres, abst, unk = counts[id(root)]
    return dict(losses=losses, n_max=len(losses), n_min=n_min,
                gain_node=gain.name or "<unnamed>",
                gain_rule=("pinned at root" if pin_gain_at_root
                           else "MRCA of present tips"),
                mrca_of_present=(mrca.name or "<unnamed>") if mrca else "",
                n_tips=n_tips, n_scored=scored, n_present=pres,
                n_absent=abst, n_unknown=unk,
                note=("" if len(losses) == n_min else
                      f"{len(losses)} loss edges collapse to {n_min} "
                      f"independent events if every polytomy carrying "
                      f"more than one is resolved against them"))


def polytomy_profile(root: lib.Node) -> dict:
    """The resolution a count would be read against, whatever the count."""
    degrees = [len(n.children) for n in lib.preorder(root) if n.children]
    poly = [d for d in degrees if d >= 3]
    unary = [d for d in degrees if d == 1]
    return dict(n_internal=len(degrees), n_polytomies=len(poly),
                n_unary=len(unary),
                max_degree=max(degrees, default=0),
                n_edges_below_polytomies=sum(poly),
                degree_counts={str(d): degrees.count(d)
                               for d in sorted(set(degrees))})
