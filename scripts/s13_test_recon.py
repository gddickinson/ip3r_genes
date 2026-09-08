"""Negative controls for the S13 reconciliation rules, run on every build.

`s5_bait_screen.self_test()`'s pattern: every rule here returns a perfectly
plausible number when it is wrong, so the tests are checks on **refusal** and
on the cases that separate a correct rule from a rule that merely runs.

Two of these were written because the rule failed them on first use, and both
would have been invisible in the output:

* **T2/T3** — the classic LCA duplication test ("the node maps where a child
  maps") is only valid on a *binary* tree. On the support-collapsed variant the
  root is a four-way polytomy whose children map to Cyclostomata,
  Gnathostomata, Cyclostomata and Gnathostomata; under the binary rule no child
  maps to Vertebrata, so the node reads as a speciation
  and the collapse would have been reported as "the duplications disappear".
* **T6** — a constrained IQ-TREE search writes no support values, so every node
  of such a tree reads as unsupported and the collapse dissolved all three
  constrained topologies into a single 134-tip polytomy, reporting it as a
  collapse of 40 nodes.

    python scripts/s13_test_recon.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s7_lib import parse_newick                                    # noqa: E402
import s13_lib as L                                                # noqa: E402
import s13_losses as LO                                            # noqa: E402
import s13_species_tree as ST                                      # noqa: E402

# A miniature species tree: ((a,b)AB,(c,d)CD)ROOT
SP_NWK = "((a:1,b:1)AB:2,(c:1,d:1)CD:2)ROOT;"


def sp_index():
    return L.SpeciesIndex.build(parse_newick(SP_NWK))


def recon(gene_nwk: str, tip_species: dict):
    sp = sp_index()
    tree = parse_newick(gene_nwk)
    return sp, tree, L.reconcile(tree, sp, tip_species)[1]


def _events_by_size(events):
    return {e.n_leaves: e for e in events}


def t1_speciation_and_duplication() -> list:
    """A clean speciation and a clean duplication, on a binary gene tree."""
    out = []
    _, _, ev = recon("((g_a:1,g_b:1):1,(g_c:1,g_d:1):1);",
                     {"g_a": "a", "g_b": "b", "g_c": "c", "g_d": "d"})
    by = _events_by_size(ev)
    if by[4].kind != "speciation" or by[4].species_label != "ROOT":
        out.append(f"T1: 1:1:1:1 tree root should be a speciation at ROOT, "
                   f"got {by[4].kind}@{by[4].species_label}")
    _, _, ev = recon("((g1_a:1,g1_b:1):1,(g2_a:1,g2_b:1):1);",
                     {"g1_a": "a", "g1_b": "b", "g2_a": "a", "g2_b": "b"})
    by = _events_by_size(ev)
    if by[4].kind != "duplication" or by[4].species_label != "AB":
        out.append(f"T1: two copies per species should duplicate at AB, "
                   f"got {by[4].kind}@{by[4].species_label}")
    return out


def t2_polytomy_duplication() -> list:
    """A polytomy whose children coexist in one ancestor is a duplication.

    Under the binary rule no child maps to ROOT, so this reads as a speciation.
    """
    out = []
    _, _, ev = recon("(g1_a:1,g2_a:1,g_c:1);",
                     {"g1_a": "a", "g2_a": "a", "g_c": "c"})
    root = _events_by_size(ev)[3]
    if root.kind != "duplication":
        out.append("T2: a polytomy with two lineages in the same species must "
                   f"be a duplication, got {root.kind}")
    if root.species_label != "ROOT":
        out.append(f"T2: mapping should be ROOT, got {root.species_label}")
    return out


def t3_polytomy_speciation() -> list:
    """A polytomy whose children occupy distinct subtrees is a speciation."""
    out = []
    _, _, ev = recon("(g_a:1,g_c:1);", {"g_a": "a", "g_c": "c"})
    root = _events_by_size(ev)[2]
    if root.kind != "speciation":
        out.append(f"T3: disjoint children must be a speciation, got {root.kind}")
    return out


def t4_extract_clade_refuses_superset() -> list:
    out = []
    tree = parse_newick("((x:1,y:1):1,z:1);")
    try:
        L.extract_clade(tree, {"x", "z"})
        out.append("T4: extract_clade accepted a clade holding an extra tip")
    except ValueError:
        pass
    got = L.extract_clade(tree, {"x", "y"})
    if set(got.leaf_names()) != {"x", "y"}:
        out.append("T4: extract_clade returned the wrong clade")
    return out


def t5_prune_tips() -> list:
    out = []
    tree = L.prune_tips(parse_newick("((x:1,y:1):1,(z:1,w:1):1);"), {"y"})
    if set(tree.leaf_names()) != {"x", "z", "w"}:
        out.append(f"T5: pruning left {sorted(tree.leaf_names())}")
    if any((not n.is_leaf) and len(n.children) == 1 for n in tree.walk()):
        out.append("T5: pruning left a unary node")
    return out


def t6_collapse_refuses_unlabelled() -> list:
    """A tree with no support values must refuse the collapse, not dissolve."""
    out = []
    try:
        L.collapse_unsupported(parse_newick("((x:1,y:1):1,(z:1,w:1):1);"))
        out.append("T6: collapse accepted a tree carrying no support labels")
    except ValueError:
        pass
    return out


def t7_collapse_is_selective() -> list:
    """Exactly the nodes below the bar collapse; those above survive."""
    out = []
    nwk = "(((x:1,y:1)100/100:1,(z:1,w:1)20/50:1)99/99:1,v:1);"
    tree, n = L.collapse_unsupported(parse_newick(nwk))
    if n != 1:
        out.append(f"T7: expected 1 collapsed node, got {n}")
    sets = {frozenset(nd.leaf_names()) for nd in tree.walk() if not nd.is_leaf}
    if frozenset({"z", "w"}) in sets:
        out.append("T7: the unsupported clade survived the collapse")
    if frozenset({"x", "y"}) not in sets:
        out.append("T7: a supported clade was collapsed")
    return out


def t8_reroot_at() -> list:
    out = []
    src = "((a:1,b:1)n1:1,(c:1,(d:1,e:1)n3:1)n2:1);"
    n_edges = len(L.edges(parse_newick(src)))
    for i in range(n_edges):
        tree = parse_newick(src)
        r = L.reroot_at(tree, L.edges(tree)[i])
        if set(r.leaf_names()) != {"a", "b", "c", "d", "e"}:
            out.append(f"T8: rerooting at edge {i} lost or gained a tip")
            break
        if any((not n.is_leaf) and len(n.children) == 1 for n in r.walk()):
            out.append(f"T8: rerooting at edge {i} left a unary node")
            break
    return out


def t9_species_key_is_binomial() -> list:
    """UniProt's common-name suffix must not make a second species."""
    out = []
    meta = {"t1": {"species": "Latimeria chalumnae"},
            "t2": {"species": "Latimeria chalumnae (Coelacanth)"}}
    if L.species_of("t1", meta) != L.species_of("t2", meta):
        out.append("T9: the same species under two spellings mapped to two "
                   "species-tree tips")
    return out


def t10_loss_verdicts() -> list:
    """Both halves of the rule that stops a bait-panel limit reading as loss."""
    out = []
    cases = [
        (False, "found_annotated", 0, 0, "sampling_artefact"),
        (False, "found_no_annotation", 0, 0, "sampling_artefact"),
        (False, "absent", 0, 0, "corroborated_loss"),
        (False, "absent", 2, 0, "paralog_unassignable"),
        (False, "absent", 0, 3, "paralog_unassignable"),
        (False, "tblastn_trace", 0, 0, "loss_with_remnant"),
        (False, "fragment", 0, 0, "undecidable"),
        (False, "assembly_gap", 0, 0, "undecidable"),
        (False, None, 0, 0, "no_genome_in_manifest"),
        (True, "absent", 0, 0, "sampled_in_gene_tree"),
    ]
    for in_tree, status, spare, unplaced, want in cases:
        got = LO.verdict(in_tree, status, spare, unplaced)
        if got != want:
            out.append(f"T10: verdict({in_tree},{status},{spare},{unplaced}) "
                       f"= {got}, expected {want}")
    return out


def t11_loss_counting() -> list:
    """Zmasek & Eddy: a speciation skipping a species node implies a loss.

    The four cases below are hand-derived on the miniature species tree
    ((a,b)AB,(c,d)CD)ROOT, because a loss count is the one number in a
    reconciliation that is easy to get plausibly wrong: the formula runs, and
    a factor-of-two error in it looks exactly like a real result.
    """
    out = []
    cases = [
        # gene tree, tip->species, expected (kind, species, losses) at the root
        ("(g_a:1,g_b:1);", {"g_a": "a", "g_b": "b"},
         ("speciation", "AB", 0)),
        # one copy in a, one in c: lost in b and in d
        ("(g_a:1,g_c:1);", {"g_a": "a", "g_c": "c"},
         ("speciation", "ROOT", 2)),
        # two copies of the same species duplicate *in that species*
        ("(g1_a:1,g2_a:1);", {"g1_a": "a", "g2_a": "a"},
         ("duplication", "a", 0)),
        # a duplication at ROOT whose second copy survives only in a
        ("((g1_a:1,g_c:1):1,g2_a:1);",
         {"g1_a": "a", "g_c": "c", "g2_a": "a"},
         ("duplication", "ROOT", 2)),
    ]
    for nwk, tips, (kind, label, losses) in cases:
        ev = recon(nwk, tips)[2]
        root = max(ev, key=lambda e: e.n_leaves)
        got = (root.kind, root.species_label, root.losses)
        if got != (kind, label, losses):
            out.append(f"T11: {nwk} -> {got}, expected "
                       f"{(kind, label, losses)}")
    return out


def t12_species_tree_check_can_fail() -> list:
    """The validator must be able to fail — on each thing it claims to catch."""
    out = []
    if ST.check():
        out.append(f"T12: the committed species tree does not validate: "
                   f"{ST.check()}")
    real = ST.CALIBRATIONS["Gnathostomata"]
    ST.CALIBRATIONS["Gnathostomata"] = (700.0, 421.0, 468.0, 563.0, "broken")
    if not ST.check():
        out.append("T12: an age older than its parent was not caught")
    ST.CALIBRATIONS["Gnathostomata"] = (462.0, 421.0, 468.0, 400.0, "broken")
    if not ST.check():
        out.append("T12: a stem age younger than the crown age was not caught")
    ST.CALIBRATIONS["Gnathostomata"] = real
    dropped = ST.CALIBRATIONS.pop("Amniota")
    if not ST.check():
        out.append("T12: an uncalibrated internal node was not caught")
    ST.CALIBRATIONS["Amniota"] = dropped
    if ST.check():
        out.append("T12: the validator did not recover after the mutations")
    return out


def t13_determinism() -> list:
    out = []
    a = [(e.kind, e.species_label, e.n_leaves, e.losses)
         for e in recon("((g1_a:1,g1_b:1):1,(g2_a:1,g_c:1):1);",
                        {"g1_a": "a", "g1_b": "b", "g2_a": "a", "g_c": "c"})[2]]
    b = [(e.kind, e.species_label, e.n_leaves, e.losses)
         for e in recon("((g1_a:1,g1_b:1):1,(g2_a:1,g_c:1):1);",
                        {"g1_a": "a", "g1_b": "b", "g2_a": "a", "g_c": "c"})[2]]
    if a != b:
        out.append("T13: two reconciliations of the same input disagreed")
    return out


def t14_paralog_sets_match_s7() -> list:
    """The clade sets used here must be the ones S7 committed."""
    out = []
    if not L.PARALOG_CLADES_TSV.exists():
        return ["T14: S7's paralog_clades.tsv is missing"]
    tree = L.load_tree(L.PHYLO_DIR / "rooted.nwk")
    sets = L.paralog_sets(tree)
    for r in L.read_tsv(L.PARALOG_CLADES_TSV):
        want = int(r["n_extended"])
        got = len(sets.get(r["paralog"], ()))
        if got != want:
            out.append(f"T14: {r['paralog']} extended clade is {got} tips here "
                       f"and {want} in S7's table")
    all_labels = set().union(*sets.values()) if sets else set()
    if len(all_labels) != sum(len(v) for v in sets.values()):
        out.append("T14: a tip belongs to two paralog clades")
    return out


TESTS = [t1_speciation_and_duplication, t2_polytomy_duplication,
         t3_polytomy_speciation, t4_extract_clade_refuses_superset,
         t5_prune_tips, t6_collapse_refuses_unlabelled,
         t7_collapse_is_selective, t8_reroot_at, t9_species_key_is_binomial,
         t10_loss_verdicts, t11_loss_counting,
         t12_species_tree_check_can_fail, t13_determinism,
         t14_paralog_sets_match_s7]


def self_test(verbose: bool = True) -> list:
    problems: list = []
    for fn in TESTS:
        got = fn()
        problems.extend(got)
        if verbose:
            print(f"  {'FAIL' if got else 'ok  '}  {fn.__name__}")
    return problems


def main() -> None:
    print(f"S13 negative controls ({len(TESTS)} groups)")
    problems = self_test()
    for p in problems:
        print(f"  ! {p}")
    if problems:
        raise SystemExit(f"{len(problems)} failure(s)")
    print("all S13 self-tests passed")


if __name__ == "__main__":
    main()
