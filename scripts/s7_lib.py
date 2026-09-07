"""S7 shared helpers — Newick parsing, rooting, splits, group metadata.

Stdlib-only tree handling for the IQ-TREE outputs: a small recursive
Newick parser (handles IQ-TREE's `)aLRT/UFBoot:length` internal labels),
bipartition extraction for monophyly tests on unrooted trees, and
outgroup re-rooting for the publication figure.

Ported from the PIEZO project's `s7_lib.py`, which carries no family
constants — the group vocabulary is read out of S6's
`representatives.tsv`, so this module never names a gene.

Used by `s7_run.py` (ML search driver), `s7_report.py` and `s7_figure.py`.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MSA_DIR = ROOT / "results" / "msa_v2"
PHYLO_DIR = ROOT / "results" / "phylogeny"
REPS_TSV = MSA_DIR / "representatives.tsv"


# ------------------------------------------------------------------ tree

@dataclass
class Node:
    name: str = ""              # tip label, or internal support label
    length: float = 0.0
    children: list["Node"] = field(default_factory=list)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list["Node"]:
        if self.is_leaf:
            return [self]
        out: list[Node] = []
        stack = [self]
        while stack:
            n = stack.pop()
            if n.is_leaf:
                out.append(n)
            else:
                stack.extend(n.children)
        return out

    def leaf_names(self) -> frozenset[str]:
        return frozenset(l.name for l in self.leaves())

    def walk(self):
        """Yield every node, preorder."""
        stack = [self]
        while stack:
            n = stack.pop()
            yield n
            stack.extend(n.children)


def parse_newick(text: str) -> Node:
    """Parse a Newick string (single tree, ';'-terminated)."""
    s = text.strip().rstrip(";")
    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse_node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
        # label (tip name or internal support), then optional :length
        start = pos
        while pos < len(s) and s[pos] not in ",():;":
            pos += 1
        node.name = s[start:pos]
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",():;":
                pos += 1
            node.length = float(s[start:pos])
        return node

    return parse_node()


def to_newick(node: Node, with_support: bool = True) -> str:
    def rec(n: Node) -> str:
        if n.is_leaf:
            return f"{n.name}:{n.length:.10g}"
        inner = ",".join(rec(c) for c in n.children)
        lab = n.name if with_support else ""
        return f"({inner}){lab}:{n.length:.10g}"
    inner = ",".join(rec(c) for c in node.children)
    return f"({inner}){node.name if with_support else ''};"


# ------------------------------------------------------------- bipartitions

def clade_map(tree: Node) -> dict[frozenset[str], Node]:
    """Leafset -> node for every internal node below the (arbitrary) root."""
    out: dict[frozenset[str], Node] = {}
    for n in tree.walk():
        if not n.is_leaf and n is not tree:
            out[n.leaf_names()] = n
    return out


def find_clade(tree: Node, labels: set[str]) -> tuple[bool, str]:
    """Is `labels` a clade of the *unrooted* tree? Returns (mono, support).

    On an unrooted tree a set is monophyletic iff some edge splits it
    exactly; with the parsed (arbitrarily rooted) tree that means some
    node's leafset equals the set, or equals its complement.
    """
    want = frozenset(labels)
    all_leaves = tree.leaf_names()
    comp = all_leaves - want
    for leafset, node in clade_map(tree).items():
        if leafset == want or leafset == comp:
            return True, node.name
    return False, ""


def find_clade_rooted(tree: Node, labels: set[str]) -> tuple[bool, str]:
    """Is `labels` a clade of *this rooted* tree? Returns (mono, support).

    The rooted question, deliberately separate from `find_clade`: on an
    unrooted tree a set and its complement are the same split, which is
    what monophyly means there — but "ITPR1 and ITPR2 are sisters" is a
    statement about the root, and answering it with the unrooted test
    would call a pair sisters whenever the *rest* of the tree happened
    to form a clade somewhere below the root.
    """
    want = frozenset(labels)
    for leafset, node in clade_map(tree).items():
        if leafset == want:
            return True, node.name
    return False, ""


def pure_clades(tree: Node, members: set[str],
                all_leaves: frozenset[str]) -> list[frozenset[str]]:
    """Maximal splits containing only `members` (unrooted sense)."""
    found: list[frozenset[str]] = []
    for n in tree.walk():
        if n is tree or n.is_leaf:
            continue
        for side in (n.leaf_names(), all_leaves - n.leaf_names()):
            if side <= members and side not in found:
                if not any(side < f for f in found):
                    found = [f for f in found if not f < side] + [side]
    # single tips can be maximal too
    for m in members:
        if not any(m in f for f in found):
            found.append(frozenset([m]))
    return sorted(found, key=len, reverse=True)


def reroot(tree: Node, outgroup: set[str]) -> Node:
    """Root on the edge that best separates `outgroup` from the rest.

    Picks the largest pure-outgroup clade edge (works even when the
    full outgroup set is not monophyletic in the gene tree); the root
    is placed at the midpoint of that edge.
    """
    best: Node | None = None
    best_n = -1
    parent: dict[int, Node] = {}
    for n in tree.walk():
        for c in n.children:
            parent[id(c)] = n
    for n in tree.walk():
        if n is tree:
            continue
        ls = n.leaf_names()
        if ls <= outgroup and len(ls) > best_n:
            best, best_n = n, len(ls)
    if best is None:  # no pure-outgroup edge: return unchanged
        return tree

    # Re-hang the tree at the midpoint of the edge above `best`.
    # path = [best, p1, p2, ..., old_root]; the edge (p_i, p_{i+1}) is
    # stored on p_i (its .length/.name) while p_i is the child — after
    # reversal it hangs above p_{i+1}, so p_{i+1} inherits p_i's attrs.
    path = [best]
    while id(path[-1]) in parent:
        path.append(parent[id(path[-1])])
    orig_len = {id(n): n.length for n in path}
    orig_name = {id(n): n.name for n in path}
    for i in range(1, len(path)):
        path[i].children.remove(path[i - 1])
    for i in range(1, len(path) - 1):
        path[i].children.append(path[i + 1])
        path[i + 1].length = orig_len[id(path[i])]
        path[i + 1].name = orig_name[id(path[i])]
    half = orig_len[id(best)] / 2.0
    best.length = half
    path[1].length = half
    # the root edge keeps its support label (unless best is a leaf)
    path[1].name = orig_name[id(best)] if not best.is_leaf else ""
    new_root = Node(children=[best, path[1]])
    # collapse any degree-2 node left behind (e.g. a rooted input tree)
    _suppress_unary(new_root)
    return new_root


def _suppress_unary(node: Node) -> None:
    for c in list(node.children):
        _suppress_unary(c)
        if len(c.children) == 1:
            gc = c.children[0]
            gc.length += c.length
            node.children[node.children.index(c)] = gc


def sister_of(tree: Node, clade: set[str]) -> frozenset[str]:
    """On a *rooted* tree, the leafset of `clade`'s sister group.

    Empty if `clade` is not a clade of this rooted tree.
    """
    want = frozenset(clade)
    for n in tree.walk():
        for i, c in enumerate(n.children):
            if c.leaf_names() == want:
                other = [d for j, d in enumerate(n.children) if j != i]
                out: set[str] = set()
                for d in other:
                    out |= d.leaf_names()
                return frozenset(out)
    return frozenset()


# ------------------------------------------------------------------ groups

def load_groups() -> dict[str, dict]:
    """label -> representatives.tsv row (group, species, rule, ...)."""
    out: dict[str, dict] = {}
    with open(REPS_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["label"] and r["label"] != "(dropped)":
                out[r["label"]] = r
    return out


def group_labels(groups: dict[str, dict], *names: str) -> set[str]:
    return {l for l, r in groups.items() if r["group"] in names}


def support_of(name: str) -> tuple[float | None, float | None]:
    """IQ-TREE internal label 'aLRT/UFBoot' -> (aLRT, UFBoot)."""
    if not name:
        return None, None
    parts = name.split("/")
    try:
        if len(parts) == 2:
            return float(parts[0]), float(parts[1])
        return None, float(parts[0])
    except ValueError:
        return None, None
