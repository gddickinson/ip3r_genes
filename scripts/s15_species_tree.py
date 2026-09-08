"""s15_species_tree.py — the tree the loss count is placed on, declared as
an input (D15) and built where the losses are.

S13 curated a 31-species, literature-calibrated species tree because a
reconciliation needs ages.  S15 cannot use it: the losses this task has to
count are cells of the **309-genome** sweep, and 278 of those species are
not in S13's tree at all.  A character matrix over 309 genomes needs a
topology over 309 genomes.

So the topology here is NCBI taxonomy, taken from the same archived
`datasets` dumps S4 built the manifest from, plus one archived call for
the 991 ancestor names.  Three consequences, all of them stated rather
than smoothed over.

**It is an input, not a result.**  Nothing in S15 estimates this tree, and
no loss placement is evidence about it.  It is the same discipline D15
applies to S13's tree, with a different source: there, curated ages; here,
a curated taxonomy.

**Its polytomies are real and are not resolved.**  A taxonomy tree
collapses every arrangement the taxonomy does not name, so a clade of
forty passerine genomes is a forty-way polytomy.  For Dollo parsimony
that is not a defect — a loss is placed at the shallowest node whose
whole subtree lacks the gene, and a polytomy makes that placement *less*
confident and never wrong — but a count of independent losses under a
polytomy is an upper bound on resolution, so `polytomy_report()` commits
the degree distribution and the report quotes it.

**It is checked against S13's tree rather than assumed compatible.**
`compare_with_s13()` asks whether every clade S13's curated tree names is
recovered as a clade here on the species the two share.  A disagreement
would mean the two halves of the project place a loss on differently
shaped branches, and it has to be visible in a table.
"""

from __future__ import annotations

import json
from pathlib import Path

import s15_lib as lib

VERT_TAX = "vertebrata_taxonomy.jsonl"
ANC_TAX = "s15_ancestor_taxonomy.jsonl"

#: taxids above the family's own root that carry no information here
_TRIM = {1, 131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 89593}


def _archive() -> Path:
    return lib.data_root() / "raw_api" / "ncbi_datasets"


def load_taxonomy() -> tuple[dict[int, dict], dict[int, str], dict[int, str]]:
    """(record by taxid, name by taxid, rank by taxid) from the archives."""
    recs: dict[int, dict] = {}
    names: dict[int, str] = {}
    ranks: dict[int, str] = {}
    for fn in (VERT_TAX, ANC_TAX):
        path = _archive() / fn
        if not path.exists():
            raise FileNotFoundError(
                f"{path} is missing; S15's tree is built offline from S4's "
                f"archived datasets dumps")
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            t = json.loads(line).get("taxonomy") or {}
            tid = t.get("tax_id")
            if tid is None:
                continue
            tid = int(tid)
            recs[tid] = t
            nm = (t.get("current_scientific_name") or {}).get("name")
            if nm:
                names[tid] = nm
            if t.get("rank"):
                ranks[tid] = str(t["rank"])
            for rank, v in (t.get("classification") or {}).items():
                if isinstance(v, dict) and v.get("id") and v.get("name"):
                    names.setdefault(int(v["id"]), v["name"])
                    ranks.setdefault(int(v["id"]), rank.upper())
    return recs, names, ranks


class Node:
    __slots__ = ("taxid", "name", "rank", "children", "parent", "tips")

    def __init__(self, taxid: int, name: str, rank: str):
        self.taxid, self.name, self.rank = taxid, name, rank
        self.children: list["Node"] = []
        self.parent: "Node | None" = None
        self.tips: list[str] = []

    @property
    def is_tip(self) -> bool:
        return not self.children

    @property
    def label(self) -> str:
        return self.name or f"taxid{self.taxid}"

    def walk(self):
        yield self
        for c in self.children:
            yield from c.walk()

    def leaves(self) -> list["Node"]:
        return [n for n in self.walk() if n.is_tip]


def build(manifest_rows: list[dict]) -> tuple[Node, dict[str, Node], list[dict]]:
    """The taxonomy tree over the manifest's genomes.

    One tip per **genome accession**, not per species: the sweep's unit of
    observation is an assembly, two accessions of one species are two
    independent observations of its gene complement, and collapsing them
    would silently average a `found` cell with a `trace` one.
    """
    recs, names, ranks = load_taxonomy()
    nodes: dict[int, Node] = {}
    audit: list[dict] = []

    def node(tid: int) -> Node:
        n = nodes.get(tid)
        if n is None:
            n = nodes[tid] = Node(tid, names.get(tid, ""), ranks.get(tid, ""))
        return n

    root: Node | None = None
    tips: dict[str, Node] = {}
    for r in manifest_rows:
        tid = int(r["taxid"])
        rec = recs.get(tid)
        if rec is None:
            audit.append(dict(accession=r["accession"], taxid=tid,
                              organism=r["organism"], status="no_taxonomy",
                              path=""))
            continue
        chain = [x for x in rec.get("parents", []) if x not in _TRIM]
        chain.append(tid)
        prev = None
        for x in chain:
            cur = node(int(x))
            if prev is not None and cur.parent is None and cur is not prev:
                cur.parent = prev
                prev.children.append(cur)
            if prev is None:
                root = root or cur
            prev = cur
        # the genome is a tip hanging off its species node
        leaf = Node(tid, r["accession"], "ASSEMBLY")
        leaf.parent = prev
        prev.children.append(leaf)
        tips[r["accession"]] = leaf
        audit.append(dict(accession=r["accession"], taxid=tid,
                          organism=r["organism"], status="placed",
                          path=" > ".join(names.get(int(x), f"taxid{x}")
                                          for x in chain)))
    if root is None:
        raise RuntimeError("no genome could be placed on the taxonomy")
    _suppress_unary(root)
    for n in root.walk():
        n.tips = [l.name for l in ([n] if n.is_tip else n.leaves())]
    return root, tips, audit


def _suppress_unary(root: Node) -> None:
    """Collapse chains of single-child internal nodes, keeping the deepest
    named one — a taxonomy path is mostly unary and every unary node would
    otherwise become a branch a loss could be placed on."""
    for n in list(root.walk()):
        while len(n.children) == 1 and not n.children[0].is_tip:
            only = n.children[0]
            n.children = only.children
            for c in n.children:
                c.parent = n
            if only.name and only.rank not in ("NO_RANK", ""):
                n.name, n.rank, n.taxid = only.name, only.rank, only.taxid


def newick(n: Node) -> str:
    if n.is_tip:
        return n.name.replace("(", "_").replace(")", "_").replace(",", "_")
    inner = ",".join(newick(c) for c in n.children)
    lab = n.label.replace(" ", "_").replace("(", "_").replace(")", "_")
    return f"({inner}){lab}"


def polytomy_report(root: Node) -> list[dict]:
    rows = []
    for n in root.walk():
        if n.is_tip:
            continue
        rows.append(dict(taxid=n.taxid, node=n.label, rank=n.rank,
                         n_children=len(n.children), n_tips=len(n.tips),
                         is_polytomy=int(len(n.children) > 2)))
    rows.sort(key=lambda r: (-r["n_children"], r["node"]))
    return rows


def _clades(root: Node) -> dict[str, frozenset]:
    return {n.label: frozenset(n.tips) for n in root.walk() if not n.is_tip}


def compare_with_s13(root: Node, tips: dict[str, Node],
                     manifest_rows: list[dict]) -> list[dict]:
    """Is every clade S13's curated tree names a clade here too?

    Restricted to the species the two trees share, and asked so that the
    two trees' different *sampling* cannot register as a disagreement:
    S13's tree carries 31 species and this one 309 assemblies, so the
    question is whether the smallest clade here containing a curated
    clade's members contains any assembly of a species S13's tree places
    **outside** that clade.  Assemblies of species S13 never sampled are
    not intruders — they are unsampled taxa, and counting them was the
    first version's error: it reported 21 of 29 curated clades as
    unrecovered while every matched node carried the right name.
    """
    import s13_lib
    sp_by_acc = {r["accession"]: r["organism"] for r in manifest_rows}
    try:
        s13_tree = s13_lib.load_tree(
            lib.RESULTS / "reconciliation" / "species_tree.nwk")
    except Exception as exc:                        # pragma: no cover
        return [dict(clade="", n_s13=0, n_shared=0, recovered=0,
                     note=f"S13 tree unreadable: {exc}")]
    import s6_lib
    acc_of_binomial: dict[str, list[str]] = {}
    for acc, org in sp_by_acc.items():
        acc_of_binomial.setdefault(s6_lib.binomial(org), []).append(acc)
    s13_binomials = {s6_lib.binomial(m) for m in _tip_names(s13_tree)}
    shared_universe = frozenset(
        a for b in s13_binomials for a in acc_of_binomial.get(b, []))
    out = []
    stack = [s13_tree]
    while stack:
        nd = stack.pop()
        kids = getattr(nd, "children", []) or []
        stack.extend(kids)
        if not kids:
            continue
        lbl = (getattr(nd, "name", "") or "").split(":")[0]
        if not lbl:
            continue
        members = [t for t in _tip_names(nd)]
        shared = [a for m in members
                  for a in acc_of_binomial.get(s6_lib.binomial(m), [])]
        if len(shared) < 2:
            out.append(dict(clade=lbl, n_s13=len(members),
                            n_shared=len(shared), recovered=-1,
                            note="fewer than two shared assemblies"))
            continue
        want = frozenset(shared)
        hit = None
        for label, tipset in _clades(root).items():
            if want <= tipset:
                if hit is None or len(tipset) < len(hit[1]):
                    hit = (label, tipset)
        universe = frozenset(a for m in members for a in
                             acc_of_binomial.get(s6_lib.binomial(m), []))
        extra = (hit[1] & shared_universe) - universe if hit else frozenset()
        out.append(dict(clade=lbl, n_s13=len(members), n_shared=len(shared),
                        recovered=int(bool(hit) and not extra),
                        matched_node=hit[0] if hit else "",
                        n_intruders=len(extra),
                        intruders=";".join(sorted(sp_by_acc[a]
                                                  for a in extra)[:4]),
                        note="" if hit and not extra else
                        ("no containing clade" if not hit
                         else f"{len(extra)} assemblies S13 places outside "
                              f"this clade are inside it here")))
    out.sort(key=lambda r: r["clade"])
    return out


def _tip_names(node) -> list[str]:
    kids = getattr(node, "children", []) or []
    if not kids:
        return [(getattr(node, "name", "") or "").split(":")[0]
                .replace("_", " ")]
    out = []
    for k in kids:
        out.extend(_tip_names(k))
    return out


AUDIT_COLS = ["accession", "organism", "taxid", "status", "path"]
POLY_COLS = ["taxid", "node", "rank", "n_children", "n_tips", "is_polytomy"]
CMP_COLS = ["clade", "n_s13", "n_shared", "recovered", "matched_node",
            "n_intruders", "intruders", "note"]
