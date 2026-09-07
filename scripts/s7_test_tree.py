"""Negative controls for the tree utilities, run on every build.

The pattern S5 and S6 use (`s5_bait_screen.self_test()`,
`s6_test_selection.py`): constructed cases, each rejected or accepted by
the rule responsible, plus the properties the step needs to have for the
result downstream of it to mean anything.

What is checked, and why each one is here rather than assumed:

  T1  Newick round-trip — the parser handles IQ-TREE's
      `)aLRT/UFBoot:length` internal labels. A parser that silently ate
      the support label would make every "well supported" call in the
      report read as missing rather than as wrong.
  T2  `find_clade` is an **unrooted** test: a set and its complement are
      the same split. Getting this wrong makes paralog monophyly depend
      on where the arbitrary root happens to sit.
  T3  `reroot` preserves the leafset, leaves no unary node, and puts the
      outgroup on one side of the root. A rerooting that drops a tip
      would change every count in the report and change nothing visible.
  T4  `sister_of` on a rooted tree returns the sister group — the
      function the ML sister answer is read off.
  T5  **The constraint newicks say what they claim, and no more.** Each
      of the three is parsed back and required to hold its pair as a
      clade with the third paralog outside it; to carry the outgroup,
      without which the split is unrooted and says nothing about
      sisters; to name each constrained tip exactly once; and — the
      check that would have caught the first version — to name **no
      free tip at all**. A taxon listed in a `-g` constraint is fixed
      outside every group the constraint declares, even from a
      top-level polytomy, so a constraint that names all 134 tips
      silently forbids the unlabelled vertebrate tips from sitting
      inside the paralog clades where this tree puts them. That cost is
      identical across the three hypotheses and rejected all of them,
      the true one included. A constraint file that is merely
      well-formed can still encode a different hypothesis from the one
      its name promises, and the AU test would answer that one.
  T6  **A low-support placement does not overrule a census label.** The
      whole point of R3 (`unconstrained`): on a constructed tree where a
      mislabelled tip sits inside another paralog at weak support, the
      rule must decline to move it. Without this the AU test would be
      asked about a membership the tree does not actually assert.
  T7  A well-supported placement *does* overrule the label (R2 fires) —
      the positive half of the same rule, so T6 cannot pass by the rule
      simply never firing.

Run:  python3 scripts/s7_test_tree.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts.s7_lib import (  # noqa: E402
    Node, find_clade, find_clade_rooted, parse_newick, pure_clades, reroot,
    sister_of, support_of, to_newick,
)
from scripts.s7_constraints import (  # noqa: E402
    PARALOGS, au_constraints, corrected_membership,
)

FAILS: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok  ' if ok else 'FAIL'} {name}" + (f" — {detail}" if detail
                                                    else ""))
    if not ok:
        FAILS.append(name)


# A binary tree (what IQ-TREE returns) with three four-tip paralog
# clades, a two-tip outgroup and a free tip; support labels in
# IQ-TREE's `aLRT/UFBoot` form.
A = "((a1:0.1,a2:0.1)99/100:0.1,(a3:0.1,a4:0.1)99/100:0.1)98/100"
B = "((b1:0.1,b2:0.1)99/100:0.1,(b3:0.1,b4:0.1)99/100:0.1)97/99"
C = "((c1:0.1,c2:0.1)99/100:0.1,(c3:0.1,c4:0.1)99/100:0.1)99/100"
TREE = (f"(((({A}:0.2,{B}:0.2)95/98:0.3,{C}:0.4)90/96:0.5,"
        f"(x1:0.2,x2:0.2)99/100:0.6)80/95:0.1,free:0.9);")

GROUPS = {f"a{i}": {"group": "ITPR1"} for i in range(1, 5)}
GROUPS |= {f"b{i}": {"group": "ITPR2"} for i in range(1, 5)}
GROUPS |= {f"c{i}": {"group": "ITPR3"} for i in range(1, 5)}
GROUPS |= {"x1": {"group": "RYR"}, "x2": {"group": "RYR"},
           "free": {"group": "invert_metazoa"}}


def t1_roundtrip() -> None:
    t = parse_newick(TREE)
    back = parse_newick(to_newick(t))
    check("T1 round-trip preserves the leafset",
          t.leaf_names() == back.leaf_names(),
          f"{len(t.leaf_names())} tips")
    labs = [n.name for n in t.walk() if not n.is_leaf and n.name]
    check("T1 support labels survive the parse",
          any(support_of(l) == (98.0, 100.0) for l in labs),
          f"{len(labs)} internal labels")


def t2_unrooted() -> None:
    t = parse_newick(TREE)
    mono, sup = find_clade(t, {"a1", "a2", "a3", "a4"})
    check("T2 a real clade is found", mono and support_of(sup) == (98.0, 100.0))
    comp = set(t.leaf_names()) - {"a1", "a2", "a3", "a4"}
    check("T2 its complement is the same split", find_clade(t, comp)[0])
    check("T2 a non-clade is rejected", not find_clade(t, {"a1", "c1"})[0])
    parts = pure_clades(t, {"a1", "a2", "c1"}, t.leaf_names())
    check("T2 pure_clades returns maximal pieces",
          sorted(len(p) for p in parts) == [1, 2],
          f"{[sorted(p) for p in parts]}")


def t3_reroot() -> None:
    t = parse_newick(TREE)
    r = reroot(t, {"x1", "x2"})
    check("T3 rerooting preserves every tip",
          r.leaf_names() == parse_newick(TREE).leaf_names())
    check("T3 no unary node survives",
          all(len(n.children) != 1 for n in r.walk()))
    sides = [c.leaf_names() for c in r.children]
    check("T3 the outgroup is one side of the root",
          frozenset({"x1", "x2"}) in sides, f"{[sorted(s) for s in sides]}")


def t4_sister() -> None:
    r = reroot(parse_newick(TREE), {"x1", "x2"})
    sis = sister_of(r, {"a1", "a2", "a3", "a4"})
    check("T4 sister_of returns the sister group",
          sis == frozenset({"b1", "b2", "b3", "b4"}), f"got {sorted(sis)}")
    check("T4 a non-clade has no sister", sister_of(r, {"a1", "c1"}) ==
          frozenset())
    # The rooted test must not accept a set merely because its
    # complement is a clade. On the rerooted tree the three paralogs
    # together form a node, so their complement {x1,x2,free} is an
    # unrooted split — and is not a rooted clade.
    split = {"x1", "x2", "free"}
    check("T4 rooted test accepts a real rooted clade",
          find_clade_rooted(r, {"a1", "a2", "a3", "a4"})[0])
    check("T4 the unrooted test accepts the complement-of-a-clade",
          find_clade(r, split)[0])
    check("T4 the rooted test rejects it",
          not find_clade_rooted(r, split)[0],
          "complement-of-a-clade must not pass the rooted test")


def t5_constraints() -> None:
    t = parse_newick(TREE)
    core, _ = corrected_membership(GROUPS, t)
    outgroup = {"x1", "x2"}
    cons = au_constraints(core, outgroup, t.leaf_names())
    check("T5 three constraints are built", len(cons) == 3, str(list(cons)))
    want = {"H1_12": ("ITPR1", "ITPR2", "ITPR3"),
            "H2_13": ("ITPR1", "ITPR3", "ITPR2"),
            "H3_23": ("ITPR2", "ITPR3", "ITPR1")}
    constrained = set().union(*(set(core[g]) for g in PARALOGS)) | outgroup
    free = set(t.leaf_names()) - constrained
    for name, nwk in cons.items():
        ct = parse_newick(nwk)
        a, b, c = want[name]
        pair = set(core[a]) | set(core[b])
        check(f"T5 {name} holds {a}+{b} as a clade",
              find_clade(ct, pair)[0])
        check(f"T5 {name} leaves {c} outside that clade",
              not (set(core[c]) & pair))
        names = [n.name for n in ct.leaves()]
        check(f"T5 {name} carries the outgroup (the split is rooted)",
              outgroup <= set(names), f"{len(outgroup)} tips")
        check(f"T5 {name} names each constrained tip exactly once",
              sorted(names) == sorted(constrained),
              f"{len(names)} of {len(constrained)}")
        # The check that would have caught the first version: a tip
        # named in a `-g` constraint is pinned outside every group the
        # constraint declares, so naming a free tip turns the
        # hypothesis into a different, stricter one.
        check(f"T5 {name} names no free tip",
              not (set(names) & free),
              f"{len(free)} free tips withheld")


def _mislabel(support: str) -> tuple[dict, Node]:
    """`b1` is nested inside the ITPR1 clade but is labelled ITPR2, at
    the given support on the node that puts it there."""
    nwk = (f"(((((({A}):0.1,b1:0.1){support}:0.2,"
           f"((b2:0.1,b3:0.1)99/100:0.1,b4:0.1)97/99:0.2)95/98:0.3,"
           f"{C}:0.4)90/96:0.5,(x1:0.2,x2:0.2)99/100:0.6)80/95:0.1,"
           f"free:0.9);")
    return dict(GROUPS), parse_newick(nwk)


def t6_weak_placement_does_not_relabel() -> None:
    g, t = _mislabel("40/70")
    _, audit = corrected_membership(g, t)
    row = next(a for a in audit if a["label"] == "b1")
    check("T6 a weakly-placed tip is left unconstrained",
          row["rule"] == "unconstrained",
          f"rule={row['rule']} support={row['alrt']}/{row['ufboot']}")


def t7_strong_placement_relabels() -> None:
    g, t = _mislabel("99/100")
    core, audit = corrected_membership(g, t)
    row = next(a for a in audit if a["label"] == "b1")
    check("T7 a strongly-placed tip is reassigned",
          row["rule"] == "reassigned" and row["assigned_to"] == "ITPR1",
          f"rule={row['rule']} -> {row['assigned_to']}")
    check("T7 and it lands in that constraint set", "b1" in core["ITPR1"]
          and "b1" not in core["ITPR2"])


def main() -> int:
    print("S7 tree-utility self-test")
    for fn in (t1_roundtrip, t2_unrooted, t3_reroot, t4_sister,
               t5_constraints, t6_weak_placement_does_not_relabel,
               t7_strong_placement_relabels):
        fn()
    if FAILS:
        print(f"\n{len(FAILS)} FAILED: " + "; ".join(FAILS))
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
