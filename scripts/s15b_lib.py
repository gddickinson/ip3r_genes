"""s15b_lib.py — S15b's data layer: the tree, its branch lengths, and the
paths every S15b module shares.

S15b adds no new measurement.  Every number it reports is read off a table
S15a committed, which is why this module is plumbing and a Newick parser
and nothing else.

Three things are worth naming.

`parse_newick()` keeps **unary nodes**.  S15a's tree wraps each accession
tip in its species node — `(GCA_033238685.1)Lates_japonicus` — because the
sweep's unit of observation is an assembly and the species name is still
the most exclusive taxon containing it.  A parser that suppressed those
would silently change the edge a loss is placed on.

`branch_lengths()` is the brief's branch-length axis, and it exists to be
shown **not** to matter to the primary count: Dollo parsimony is a count of
edges and does not read a length at all.  Lengths enter only the Mk
section, so all three schemes are declared here and the report says which
result each one could possibly move.  `calibrated` joins S13's committed
node ages by name (23 of the 29 calibrated nodes are present in this tree)
and spaces the uncalibrated nodes evenly between their nearest calibrated
ancestor and their descendants — the ages are an **input** (D15), not
something S15b estimates.

`OUT` is a directory of its own rather than more files in
`results/loss_dynamics/`, because S15a's directory is the instrument and
S15b's is the count, and a reader who wants to know what was measured
should not have to sort one from the other.
"""

from __future__ import annotations

import json
from pathlib import Path

import s15_lib as base

PROJECT = base.PROJECT
RESULTS = base.RESULTS
S15A = RESULTS / "loss_dynamics"
OUT = RESULTS / "loss_counts"
FIGS = OUT / "figures"

MATRIX = S15A / "character_matrix.tsv"
CANDIDATES = S15A / "loss_candidates.tsv"
POLYTOMIES = S15A / "tree_polytomies.tsv"
NEWICK = S15A / "species_tree_309.nwk"
INTEGRITY_LOCI = S15A / "integrity_loci.tsv"
INTEGRITY_PAIRS = S15A / "integrity_pairs.tsv"
INTEGRITY_TESTS = S15A / "integrity_tests.tsv"
S15A_STATS = S15A / "loss_dynamics_stats.json"
CALIBRATIONS = RESULTS / "reconciliation" / "species_tree_calibrations.tsv"

read_tsv = base.read_tsv
write_tsv = base.write_tsv
sha256 = base.sha256
median = base.median
sign_test = base.sign_test

ITPR_CELLS = base.ITPR_CELLS

#: branch-length schemes, in the order the sensitivity table walks them
BL_SCHEMES = ("unit", "ultrametric", "calibrated")


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    """Dashboard live panel, tagged S15b rather than S15a."""
    payload = dict(task="S15b", stage=stage, workers=1,
                   steps=[dict(label=s, done=bool(d)) for s, d in steps])
    p = RESULTS / "session_live.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(payload, indent=1), encoding="utf-8")


# ----------------------------------------------------------------- the tree

class Node:
    __slots__ = ("name", "children", "parent", "length")

    def __init__(self, name: str = ""):
        self.name = name
        self.children: list[Node] = []
        self.parent: Node | None = None
        self.length: float = 1.0

    @property
    def is_tip(self) -> bool:
        return not self.children

    def __repr__(self) -> str:                          # pragma: no cover
        return f"Node({self.name!r}, {len(self.children)} children)"


def parse_newick(text: str) -> Node:
    """Parse a Newick string, keeping unary nodes and internal labels."""
    s = text.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = 0

    def label() -> tuple[str, float]:
        nonlocal pos
        start = pos
        while pos < len(s) and s[pos] not in "(),":
            pos += 1
        tok = s[start:pos].strip()
        if ":" in tok:
            nm, _, ln = tok.partition(":")
            try:
                return nm.strip(), float(ln)
            except ValueError:
                return nm.strip(), 1.0
        return tok, 1.0

    def node() -> Node:
        nonlocal pos
        n = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                child = node()
                child.parent = n
                n.children.append(child)
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                break
            if pos >= len(s) or s[pos] != ")":
                raise ValueError(f"unbalanced newick at {pos}")
            pos += 1
        n.name, n.length = label()
        return n

    root = node()
    if pos != len(s):
        raise ValueError(f"trailing newick text at {pos}: {s[pos:pos + 40]!r}")
    return root


def load_tree(path: Path | None = None) -> Node:
    return parse_newick(Path(path or NEWICK).read_text(encoding="utf-8"))


def postorder(root: Node):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(n.children)
    return reversed(out)


def preorder(root: Node):
    stack, out = [root], []
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(reversed(n.children))
    return out


def tips(root: Node) -> list[Node]:
    return [n for n in preorder(root) if n.is_tip]


def to_newick(root: Node, lengths: bool = False) -> str:
    def enc(n: Node) -> str:
        inner = ("(" + ",".join(enc(c) for c in n.children) + ")"
                 if n.children else "")
        suffix = f":{n.length:g}" if lengths and n.parent is not None else ""
        return f"{inner}{n.name}{suffix}"
    return enc(root) + ";"


# ------------------------------------------------------------ branch lengths

def _node_ages(root: Node, cal: dict[str, float]) -> dict[int, float]:
    """Age in Ma for every node: calibrated where named, interpolated else.

    A node with no calibration is placed evenly between its nearest
    calibrated ancestor and 0 (the present), by how many edges separate
    them, so the interpolation cannot invert an age it was given.
    """
    ages: dict[int, float] = {}
    for n in preorder(root):
        if n.name in cal:
            ages[id(n)] = cal[n.name]
    root_age = ages.get(id(root)) or max(cal.values(), default=1.0)
    ages[id(root)] = root_age
    for n in preorder(root):
        if id(n) in ages:
            continue
        anc, up = n.parent, 1
        while anc is not None and id(anc) not in ages:
            anc, up = anc.parent, up + 1
        top = ages[id(anc)] if anc is not None else root_age
        # depth in edges from this node down to its farthest tip
        down = _height_edges(n)
        ages[id(n)] = top * down / (down + up) if (down + up) else 0.0
    for n in tips(root):
        ages[id(n)] = 0.0
    return ages


def _height_edges(n: Node) -> int:
    if n.is_tip:
        return 0
    return 1 + max(_height_edges(c) for c in n.children)


def load_calibrations(path: Path | None = None) -> dict[str, float]:
    rows = read_tsv(Path(path or CALIBRATIONS))
    out = {}
    for r in rows:
        try:
            out[r["node"]] = float(r["age_ma"])
        except (KeyError, TypeError, ValueError):
            continue
    return out


def branch_lengths(root: Node, scheme: str,
                   cal: dict[str, float] | None = None) -> dict[int, float]:
    """Edge length per node id under one of `BL_SCHEMES`.

    `unit`         every edge 1 — the scheme parsimony implicitly uses.
    `ultrametric`  every root-to-tip path the same total length, edges
                   spaced by the number of nodes below them.
    `calibrated`   S13's committed node ages (D15), interpolated where this
                   tree carries a node S13 did not name.
    """
    if scheme not in BL_SCHEMES:
        raise ValueError(f"unknown branch-length scheme {scheme!r}")
    if scheme == "unit":
        return {id(n): (0.0 if n.parent is None else 1.0)
                for n in preorder(root)}
    if scheme == "ultrametric":
        h = {id(n): _height_edges(n) for n in preorder(root)}
        total = h[id(root)] or 1
        return {id(n): (0.0 if n.parent is None else
                        (h[id(n.parent)] - h[id(n)]) / total)
                for n in preorder(root)}
    ages = _node_ages(root, cal if cal is not None else load_calibrations())
    return {id(n): (0.0 if n.parent is None else
                    max(0.0, ages[id(n.parent)] - ages[id(n)]))
            for n in preorder(root)}


def apply_lengths(root: Node, lengths: dict[int, float]) -> None:
    for n in preorder(root):
        n.length = lengths.get(id(n), 1.0)


# -------------------------------------------------------------- small stats

def bh(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg, stdlib (S9's rule: one family, one correction)."""
    m = len(pvals)
    if not m:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [0.0] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, pvals[i] * m / k)
        q[i] = min(1.0, prev)
    return q


def pfmt(p: float | str) -> str:
    """A p-value at four decimal places renders 2.1e-07 as `0.0000`."""
    try:
        v = float(p)
    except (TypeError, ValueError):
        return str(p)
    if v == 0:
        return "0"
    return f"{v:.3g}"
