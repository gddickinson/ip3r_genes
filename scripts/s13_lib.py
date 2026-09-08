"""S13 — reconciliation helpers: tip mapping, LCA reconciliation, losses.

The reconciliation itself is the classic LCA algorithm (Goodman et al. 1979;
Page 1994; loss counting after Zmasek & Eddy 2001), implemented here in the
stdlib rather than pulled in as a dependency, for the same reason S7 carries
its own Newick code: the trees are small, the algorithm is exact and
deterministic, and a hand-checkable implementation is worth more than a black
box for a claim this load-bearing.

Four project-specific rules are baked in.

* **Tips are labelled by the S7 tree, not the census.** S9 established that
  for the selection sets (D30 read forwards); the same applies here. A tip's
  paralog is the extended paralog clade of `rooted.nwk` that holds it, and a
  tip inside none of the three is `unplaced` — never assigned by its database
  name.
* **A rooted gene tree is required.** Reconciliation is meaningless on an
  unrooted tree, so the vertebrate subtree is extracted from the
  RyR-outgroup-rooted S7 tree and its own root is the vertebrate ITPR MRCA.
* **Species come through `s6_lib.binomial`.** The census carries
  `Latimeria chalumnae` from the genome sweep and `Latimeria chalumnae
  (Coelacanth)` from UniProt; a raw-string key makes that two species and the
  reconciliation then reports a duplication that is a string difference.
* **Polytomies are first-class.** `collapse_unsupported` turns nodes below
  S7's own thresholds into multifurcations, and the LCA mapping handles a
  node with any number of children — which matters because the node the
  deepest duplication placement rests on sits at SH-aLRT 17.4 / UFBoot 54.
"""

from __future__ import annotations

import csv
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s6_lib import binomial                                     # noqa: E402
from s7_lib import (Node, parse_newick, support_of,             # noqa: E402,F401
                    to_newick, load_groups, group_labels)

RECON_DIR = ROOT / "results" / "reconciliation"
PHYLO_DIR = ROOT / "results" / "phylogeny"
REPS_TSV = ROOT / "results" / "msa_v2" / "representatives.tsv"
PARALOG_CLADES_TSV = PHYLO_DIR / "paralog_clades.tsv"

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
VERT_GROUPS = ("ITPR1", "ITPR2", "ITPR3", "vertebrate_basal")

# S7's own "well supported" bar, reused unchanged so a collapse here and a
# claim there mean the same thing.
MIN_ALRT, MIN_UFBOOT = 80.0, 95.0


# ------------------------------------------------------------------ mapping

def tip_metadata() -> dict[str, dict]:
    """gene-tree tip label -> representative row (species, group, class...)."""
    return load_groups()


def species_of(label: str, meta: dict[str, dict]) -> str:
    """Species tag as it appears in the species tree (underscored binomial)."""
    row = meta.get(label)
    if not row:
        raise KeyError(f"tip {label!r} has no representatives.tsv row")
    return binomial(row["species"]).replace(" ", "_")


def vertebrate_labels(meta: dict[str, dict]) -> set[str]:
    return {lab for lab, r in meta.items() if r["group"] in VERT_GROUPS}


# ------------------------------------------------------- species-tree indices

@dataclass
class SpeciesIndex:
    """Depth, parent and leafset for every node of the species tree."""
    root: Node
    depth: dict[int, int] = field(default_factory=dict)
    parent: dict[int, "Node | None"] = field(default_factory=dict)
    leafset: dict[int, frozenset] = field(default_factory=dict)
    by_name: dict[str, Node] = field(default_factory=dict)

    @classmethod
    def build(cls, root: Node) -> "SpeciesIndex":
        idx = cls(root=root)

        def rec(n: Node, parent: "Node | None", d: int) -> frozenset:
            idx.depth[id(n)] = d
            idx.parent[id(n)] = parent
            if n.name:
                idx.by_name[n.name] = n
            ls = (frozenset({n.name}) if n.is_leaf
                  else frozenset().union(*(rec(c, n, d + 1) for c in n.children)))
            idx.leafset[id(n)] = ls
            return ls

        rec(root, None, 0)
        return idx

    def lca(self, a: Node, b: Node) -> Node:
        pa, pb = self.path_to_root(a), self.path_to_root(b)
        seen = {id(n) for n in pa}
        for n in pb:
            if id(n) in seen:
                return n
        return self.root

    def path_to_root(self, n: Node) -> list[Node]:
        out, cur = [], n
        while cur is not None:
            out.append(cur)
            cur = self.parent[id(cur)]
        return out

    def leaf(self, species: str) -> Node:
        node = self.by_name.get(species)
        if node is None or not node.is_leaf:
            raise KeyError(f"{species!r} is not a tip of the species tree")
        return node

    def label(self, n: Node) -> str:
        return n.name or f"<unnamed:{sorted(self.leafset[id(n)])[:2]}>"


# ------------------------------------------------------------ reconciliation

@dataclass
class Event:
    gene_node: Node
    kind: str                 # "duplication" | "speciation"
    species_node: Node
    species_label: str
    n_leaves: int
    leaf_groups: tuple
    losses: int
    support: str


def _is_duplication(m: Node, child_maps: list, sp: SpeciesIndex) -> bool:
    """Non-binary LCA duplication test (Vernot et al. 2008).

    True when a child maps to M(v) itself, or when two children map into the
    same child-subtree of M(v) — i.e. two descendant lineages that coexisted
    in the same ancestral species.
    """
    if any(cm is m for cm in child_maps):
        return True
    slots = []
    for cm in child_maps:
        path = sp.path_to_root(cm)
        try:
            slot = path[path.index(m) - 1]      # the child of m on cm's path
        except (ValueError, IndexError):
            return True                          # cm is not under m: coexistence
        slots.append(id(slot))
    return len(set(slots)) < len(slots)


def reconcile(gene_root: Node, sp: SpeciesIndex, tip_species: dict[str, str],
              group_of=None) -> tuple[dict, list[Event]]:
    """LCA-reconcile a rooted gene tree against a species tree.

    A gene node is a **duplication** when two of its descendant lineages
    coexist in the same ancestral species, and a speciation otherwise. On a
    binary tree that is the familiar test "the node's species mapping equals a
    child's". On a **multifurcation it is not**, and the difference is not
    cosmetic: the support-collapsed variant turns the root into a four-way
    polytomy whose children map to Cyclostomata, Gnathostomata, Cyclostomata
    and Gnathostomata, and under the binary rule no child maps to Vertebrata so
    the node reads as a *speciation* — the collapse would have
    been reported as the duplications disappearing when in fact the rule could
    not see them. The test used is therefore the standard non-binary
    generalisation (Vernot et al. 2008): map each child into the child-subtree
    of M(v) that contains it, and call a duplication when a child maps to M(v)
    itself or when two children land in the same subtree.

    Losses are counted per Zmasek & Eddy 2001 from the depth gap each child
    crosses.
    """
    mapping: dict[int, Node] = {}

    def rec_map(g: Node) -> Node:
        if g.is_leaf:
            node = sp.leaf(tip_species[g.name])
        else:
            kids = [rec_map(c) for c in g.children]
            node = kids[0]
            for k in kids[1:]:
                node = sp.lca(node, k)
        mapping[id(g)] = node
        return node

    rec_map(gene_root)

    events: list[Event] = []
    for g in gene_root.walk():
        if g.is_leaf:
            continue
        m = mapping[id(g)]
        child_maps = [mapping[id(c)] for c in g.children]
        is_dup = _is_duplication(m, child_maps, sp)
        gap = 0
        for cm in child_maps:
            d = sp.depth[id(cm)] - sp.depth[id(m)]
            gap += d if is_dup else max(d - 1, 0)
        groups = tuple(sorted({group_of(l.name) for l in g.leaves()})) \
            if group_of else ()
        events.append(Event(gene_node=g,
                            kind="duplication" if is_dup else "speciation",
                            species_node=m, species_label=sp.label(m),
                            n_leaves=len(g.leaves()), leaf_groups=groups,
                            losses=gap, support=g.name or ""))
    return mapping, events


# ------------------------------------------------------------- subtree tools

def extract_clade(root: Node, wanted: set) -> Node:
    """Smallest clade of `root` containing every label in `wanted`.

    Raises if that clade also carries labels outside `wanted` — a silent
    superset would put non-vertebrate tips into a vertebrate reconciliation
    and every duplication count downstream would be wrong.
    """
    best: Node | None = None
    for n in root.walk():
        if n.is_leaf:
            continue
        labels = set(n.leaf_names())
        if wanted <= labels and (best is None or len(labels) < len(best.leaf_names())):
            best = n
    if best is None:
        raise ValueError("no clade contains all requested labels")
    extra = set(best.leaf_names()) - wanted
    if extra:
        raise ValueError(f"smallest containing clade also holds {sorted(extra)}")
    return best


def prune_tips(root: Node, drop: set) -> Node:
    """Remove tips and collapse the degree-2 nodes left behind.

    Dropping a tip by filtering the *wanted* set is not enough: the tip is
    still physically in the tree, so the smallest containing clade still holds
    it and it still maps to a species. A sensitivity analysis that "excludes"
    a taxon has to actually remove it.
    """
    def rec(n: Node) -> "Node | None":
        if n.is_leaf:
            return None if n.name in drop else n
        kept = [c for c in (rec(c) for c in n.children) if c is not None]
        if not kept:
            return None
        if len(kept) == 1:                     # suppress the unary node
            kept[0].length += n.length
            return kept[0]
        n.children = kept
        return n

    out = rec(root)
    if out is None:
        raise ValueError("pruning removed every tip")
    return out


def collapse_unsupported(root: Node, min_alrt: float = MIN_ALRT,
                         min_ufboot: float = MIN_UFBOOT) -> tuple[Node, int]:
    """Collapse internal edges below S7's own support bar into polytomies.

    A resolution the tree does not support is not evidence, and reconciliation
    reads resolution as event structure: a poorly supported node is silently
    promoted into a duplication with an age bracket. Collapsing it moves the
    question to "does the placement survive when the unsupported resolution is
    removed", which is answerable. The root is never collapsed (it has no
    edge), and support is read with S7's own parser.

    **Refuses a tree with no support labels.** A constrained IQ-TREE search
    writes none, and every node of such a tree then reads as unsupported: the
    first run of this variant dissolved the constrained topologies down to a
    single polytomy of all 134 tips and reported it as a collapse of 40 nodes.
    Grading an unmeasured support as a failing one is the same error
    `s7_bnni.py` refuses when it compares a constrained tree on topology
    alone.
    """
    labelled = sum(1 for n in root.walk()
                   if not n.is_leaf and n is not root
                   and support_of(n.name)[1] is not None)
    if not labelled:
        raise ValueError("tree carries no support labels: nothing to collapse "
                         "against (a constrained search computes none)")
    n_collapsed = 0

    def rec(n: Node) -> Node:
        if n.is_leaf:
            return n
        n.children = [rec(c) for c in n.children]
        kids: list[Node] = []
        for c in n.children:
            if c.is_leaf:
                kids.append(c)
                continue
            alrt, ufboot = support_of(c.name)
            ok = (alrt is not None and alrt >= min_alrt
                  and ufboot is not None and ufboot >= min_ufboot)
            if ok:
                kids.append(c)
            else:
                nonlocal n_collapsed
                n_collapsed += 1
                kids.extend(c.children)         # dissolve the edge
        n.children = kids
        return n

    return rec(root), n_collapsed


def edges(root: Node) -> list[Node]:
    """Every node that has a parent, i.e. every rerootable edge."""
    return [n for n in root.walk() if n is not root]


def reroot_at(root: Node, target: Node) -> Node:
    """Re-root the tree on the edge above `target`, returning a new root.

    Mutates `root`, so the caller re-parses a fresh tree per candidate edge.
    Used only by the minimum-event rooting check: it asks where a
    reconciliation *would* put the root if the RyR outgroup had not already
    decided, which is an independent test of the rooting the answer uses.
    """
    parent: dict[int, Node] = {}
    for n in root.walk():
        for c in n.children:
            parent[id(c)] = n
    if id(target) not in parent:
        return root
    chain: list[Node] = []
    cur = parent[id(target)]
    while True:
        chain.append(cur)
        nxt = parent.get(id(cur))
        if nxt is None:
            break
        cur = nxt
    chain[0].children = [c for c in chain[0].children if c is not target]
    for i in range(len(chain) - 1):                 # detach every chain edge
        par = chain[i + 1]
        par.children = [c for c in par.children if c is not chain[i]]
    for i in range(len(chain) - 1):                 # re-hang them inverted
        chain[i].children.append(chain[i + 1])
    return Node(name="", length=0.0,
                children=[target, _suppress_unary(chain[0])])


def _suppress_unary(n: Node) -> Node:
    """Drop degree-2 internal nodes, which rerooting always leaves behind."""
    while (not n.is_leaf) and len(n.children) == 1:
        child = n.children[0]
        child.length += n.length
        n = child
    if not n.is_leaf:
        n.children = [_suppress_unary(c) for c in n.children]
    return n


# --------------------------------------------------------------- paralog sets

def paralog_sets(tree: Node) -> dict[str, frozenset]:
    """S7's extended paralog clades, read from its committed table.

    Read, not re-derived: `s9_sets.py` re-derives them and cross-checks, and a
    third derivation that disagreed would be a silent fork in what the project
    means by "ITPR2". The table lists the *added* unlabelled tips explicitly,
    so the set is the census-labelled tips of that paralog (restricted to this
    tree) plus those.
    """
    groups = load_groups()
    leaves = set(tree.leaf_names())
    out: dict[str, frozenset] = {}
    with PARALOG_CLADES_TSV.open() as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            p = r["paralog"]
            base = group_labels(groups, p) & leaves
            added = {x for x in r["added"].split(";") if x} & leaves
            out[p] = frozenset(base | added)
    return out


def group_resolver(tree: Node):
    """label -> tree-resolved paralog, or 'unplaced'."""
    sets = paralog_sets(tree)
    lookup: dict[str, str] = {}
    for p, labels in sets.items():
        for lab in labels:
            lookup[lab] = p

    def group_of(label: str) -> str:
        return lookup.get(label, "unplaced")

    return group_of, sets


# ---------------------------------------------------------------------- io

def load_tree(path: Path) -> Node:
    return parse_newick(path.read_text())


def load_species_tree() -> SpeciesIndex:
    return SpeciesIndex.build(load_tree(RECON_DIR / "species_tree.nwk"))


def calibrations() -> dict[str, dict]:
    path = RECON_DIR / "species_tree_calibrations.tsv"
    with path.open() as fh:
        return {r["node"]: r for r in csv.DictReader(fh, delimiter="\t")}


def write_tsv(path: Path, header: list, rows: list) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def read_tsv(path: Path) -> list[dict]:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ------------------------------------------------- report-side matrix helpers
# Both halves of the report read the matrix through these, for the same reason
# they share a loader: a second definition that filtered the skipped cells
# differently would make the two halves disagree about how many
# reconciliations ran.

def matrix_cells(summary_rows: list) -> list:
    """The cells that actually ran — a skipped cell carries a `note`."""
    return [r for r in summary_rows if not r.get("note")]


def deepest_nodes(summary_rows: list, variant: str) -> set:
    """Species nodes the deepest paralog duplication maps to, for one variant."""
    return {r["deepest_paralog_dup_node"]
            for r in matrix_cells(summary_rows) if r["variant"] == variant}
