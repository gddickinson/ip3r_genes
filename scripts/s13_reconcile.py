"""S13 step 2 — reconcile the S7 gene trees against the species tree.

Run over a **matrix**, not a single tree, because three things S7 left open
each move the answer and none can be argued away:

  * *topology* — the AU test left two of the three sister arrangements
    standing (`ml` and `H3_23`, which agree on ITPR2+ITPR3) and rejected two.
    All four constrained/unconstrained topologies are reconciled, and the
    `--bnni` re-search is added as a fifth: it is the same data under the
    model-violation guard, so a placement that moves between them is a
    placement that depends on the model.
  * *the cyclostome tips* — six loci in two cyclostome-only clades, one of
    them sister to the ITPR2+ITPR3 clade. They are the **only** tips that can
    put a duplication below the cyclostome-gnathostome split, i.e. the entire
    weight of a "these are 2R ohnologs" reading rests on them.
  * *resolution the tree does not support* — the node separating ITPR1 from
    (ITPR2+ITPR3+cyclostome) sits at SH-aLRT 17.4 / UFBoot 54. Reconciliation
    reads resolution as event structure and cannot tell a supported node from
    an unsupported one, so a third variant collapses every node below S7's own
    bar into a polytomy and asks the same question of what is left.

Every cell of that matrix is reported. Where the cells agree, the result is
robust; where they disagree, the disagreement *is* the result.

A fourth check runs alongside: the reported answer uses the RyR-outgroup
rooting, so `--rooting-check` re-roots the vertebrate subtree at every edge
in turn and reports where the minimum-event root would sit. If the outgroup
rooting and the minimum-event rooting name the same edge, the placement does
not depend on which of the two decided it.

    python scripts/s13_reconcile.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from s13_lib import (PARALOGS, PHYLO_DIR, RECON_DIR, SpeciesIndex,  # noqa: E402
                     calibrations, collapse_unsupported, edges, extract_clade,
                     group_resolver, load_species_tree, load_tree, prune_tips,
                     reconcile, reroot_at, species_of, tip_metadata, to_newick,
                     vertebrate_labels, write_tsv)
from s7_lib import group_labels, load_groups, reroot                # noqa: E402

TOPOLOGIES = {
    "ml":     (PHYLO_DIR / "rooted.nwk", "unconstrained ML tree (the reported tree)"),
    "H1_12":  (PHYLO_DIR / "au" / "H1_12.treefile", "ITPR1+ITPR2 sisters (AU-rejected, p 1.8e-05)"),
    "H2_13":  (PHYLO_DIR / "au" / "H2_13.treefile", "ITPR1+ITPR3 sisters (AU-rejected, p 1.65e-05)"),
    "H3_23":  (PHYLO_DIR / "au" / "H3_23.treefile", "ITPR2+ITPR3 sisters (not rejected, p 0.476)"),
    "bnni":   (PHYLO_DIR / "itpr_ml_bnni.treefile", "the --bnni model-violation guard re-search"),
}

VARIANTS = ("with_cyclostome", "without_cyclostome", "support_collapsed")


def vertebrate_subtree(path: Path, wanted: set, drop: set, collapse: bool):
    """Rooted vertebrate ITPR subtree of a tree file.

    The AU constraint searches wrote *unrooted* trees. Reconciliation on an
    unrooted tree is meaningless — the root decides which node is the
    ancestral duplication — so each is re-rooted on the ryanodine receptors,
    the same outgroup S7 used, before the subtree is taken.
    """
    tree = load_tree(path)
    groups = load_groups()
    n_collapsed = 0
    if path.name != "rooted.nwk":
        tree = reroot(tree, group_labels(groups, "RYR"))
    if collapse:
        tree, n_collapsed = collapse_unsupported(tree)
    if drop:
        tree = prune_tips(tree, drop)
    return extract_clade(tree, wanted), n_collapsed


def bracket(species_label: str, sp: SpeciesIndex, cal: dict) -> tuple:
    """Age bracket for a duplication mapped to `species_label`.

    A duplication mapped to species node X happened on the branch *above* X:
    after the split that gave X its stem, before the split at X. The bracket
    is therefore [age(X), stem_age(X)] — and `stem_age` is a literature age
    for the real sister split, not the sampled parent's age, so a duplication
    landing on a node whose neighbours went unsampled is not handed a bracket
    hundreds of millions of years wider than the evidence (see
    `s13_species_tree.py`). At the root it is open-ended and reported as such.
    """
    node = sp.by_name.get(species_label)
    if node is None:
        return ("", "", "")
    row = cal.get(species_label, {})
    young = row.get("age_ma", "")
    if sp.parent[id(node)] is None:
        return (str(young), "unbounded",
                "root of the sampled species tree: the duplication is on the "
                "vertebrate stem and nothing in this tree bounds it from above")
    return (str(young), row.get("stem_age_ma", ""), row.get("source", ""))


def paralog_ancestors(events, group_of):
    """The duplication nodes that separate the paralog subfamilies.

    Ranked oldest-first by how many paralogs their descendants span, then by
    size. A duplication whose descendants carry only `unplaced` tips is not a
    paralog-creating event as far as this project can tell, and is kept in the
    table with its groups shown rather than filtered out.
    """
    dups = [e for e in events if e.kind == "duplication"]
    spanning = [e for e in dups
                if len(set(e.leaf_groups) & set(PARALOGS)) >= 2]
    return sorted(spanning, key=lambda e: (-len(set(e.leaf_groups) & set(PARALOGS)),
                                           -e.n_leaves))


def rooting_check(path: Path, wanted: set, sp: SpeciesIndex, meta: dict,
                  group_of) -> list:
    """Duplications+losses over every rooting of the vertebrate subtree.

    The reported answer uses the RyR-outgroup rooting, which is evidence from
    outside this subtree. Minimum-event rooting is a different criterion
    entirely, and the check is worth running because the two can disagree —
    what matters is not whether they name the same edge but whether the
    *duplication placement* moves when they do. The outgroup rooting is in the
    table as `edge_index -1` so the two are read off the same scale.
    """
    def score(rooted):
        tip_species = {lab: species_of(lab, meta) for lab in rooted.leaf_names()}
        _, ev = reconcile(rooted, sp, tip_species, group_of)
        n_dup = sum(1 for e in ev if e.kind == "duplication")
        n_loss = sum(e.losses for e in ev)
        root_ev = next((e for e in ev if e.gene_node is rooted), None)
        anc = paralog_ancestors(ev, group_of)
        return (n_dup, n_loss, root_ev.kind if root_ev else "",
                root_ev.species_label if root_ev else "",
                anc[0].species_label if anc else "")

    rows = []
    base, _ = vertebrate_subtree(path, wanted, set(), False)
    d, l, k, rs, deep = score(base)
    rows.append([-1, len(base.leaf_names()), d, l, d + l, k, rs, deep,
                 "outgroup rooting (RyR) - the one every other table uses"])
    n_edges = len(edges(base))
    for i in range(n_edges):
        sub, _ = vertebrate_subtree(path, wanted, set(), False)
        target = edges(sub)[i]
        size = len(target.leaf_names())
        rooted = reroot_at(sub, target)
        d, l, k, rs, deep = score(rooted)
        rows.append([i, size, d, l, d + l, k, rs, deep,
                     ";".join(sorted(target.leaf_names())[:2])])
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=RECON_DIR)
    ap.add_argument("--rooting-check", action="store_true", default=True)
    ap.add_argument("--no-rooting-check", dest="rooting_check",
                    action="store_false")
    args = ap.parse_args()

    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    meta = tip_metadata()
    sp = load_species_tree()
    cal = calibrations()
    vert = vertebrate_labels(meta)

    ml_tree = load_tree(PHYLO_DIR / "rooted.nwk")
    group_of, psets = group_resolver(ml_tree)
    cyclo = sorted(l for l in vert
                   if meta[l]["species"].split()[0] in ("Myxine", "Petromyzon"))
    print(f"vertebrate tips: {len(vert)}; cyclostome tips: {len(cyclo)}")
    for p in PARALOGS:
        print(f"  {p}: {len(psets[p])} tips (S7 extended clade)")

    event_rows, summary_rows, matrix_rows = [], [], []

    for topo, (path, note) in TOPOLOGIES.items():
        if not path.exists():
            print(f"  ! {topo}: {path} missing - skipped")
            continue
        for variant in VARIANTS:
            drop = set(cyclo) if variant == "without_cyclostome" else set()
            collapse = variant == "support_collapsed"
            wanted = vert - drop
            try:
                clade, n_collapsed = vertebrate_subtree(path, wanted, drop,
                                                        collapse)
            except ValueError as exc:
                # A constrained search writes no support values, so the
                # collapse variant has nothing to read. Recorded as skipped
                # with its reason rather than run on absent evidence.
                matrix_rows.append([topo, note, variant] + [""] * 10
                                   + [f"skipped: {exc}"])
                print(f"  {topo:8s} {variant:19s} skipped - {exc}")
                continue
            tip_species = {lab: species_of(lab, meta) for lab in clade.leaf_names()}
            _, events = reconcile(clade, sp, tip_species, group_of)
            n_dup = sum(1 for e in events if e.kind == "duplication")
            n_loss = sum(e.losses for e in events)
            anc = paralog_ancestors(events, group_of)
            root_ev = next((e for e in events if e.gene_node is clade), None)

            for e in events:
                event_rows.append([topo, variant, e.kind, e.species_label,
                                   e.n_leaves, "+".join(e.leaf_groups),
                                   e.losses, e.support,
                                   len(e.gene_node.children)])
                # Annotate the tree in place so the committed Newick carries
                # the reconciliation, not just a table beside it.
                e.gene_node.name = (("D@" if e.kind == "duplication" else "S@")
                                    + e.species_label)
            if topo == "ml":
                (out_dir / f"reconciled_ml_{variant}.nwk").write_text(
                    to_newick(clade) + "\n")
            for rank, e in enumerate(anc[:6], start=1):
                lo, hi, src = bracket(e.species_label, sp, cal)
                summary_rows.append([topo, variant, rank, e.species_label,
                                     "+".join(e.leaf_groups), e.n_leaves,
                                     lo, hi, e.support, src])
            matrix_rows.append([
                topo, note, variant, len(clade.leaf_names()), n_collapsed,
                n_dup, n_loss,
                root_ev.species_label if root_ev else "",
                root_ev.kind if root_ev else "",
                anc[0].species_label if anc else "",
                "+".join(anc[0].leaf_groups) if anc else "",
                anc[0].support if anc else "",
                len(anc), ""])
            print(f"  {topo:8s} {variant:19s} tips {len(clade.leaf_names()):3d}  "
                  f"collapsed {n_collapsed:3d}  dups {n_dup:3d}  "
                  f"losses {n_loss:4d}  root->{root_ev.species_label if root_ev else '?'}"
                  f"  deepest paralog dup->{anc[0].species_label if anc else '-'}")

    write_tsv(out_dir / "reconciliation_events.tsv",
              ["topology", "variant", "event", "species_node", "n_gene_leaves",
               "paralog_groups", "implied_losses", "gene_node_support",
               "n_children"], event_rows)
    write_tsv(out_dir / "duplication_placement.tsv",
              ["topology", "variant", "rank", "species_node", "paralog_groups",
               "n_gene_leaves", "age_young_ma", "age_old_ma",
               "gene_node_support", "age_source"], summary_rows)
    write_tsv(out_dir / "reconciliation_summary.tsv",
              ["topology", "topology_note", "variant", "n_tips",
               "n_nodes_collapsed", "n_duplications", "implied_losses",
               "root_species_node", "root_event", "deepest_paralog_dup_node",
               "deepest_paralog_dup_groups", "deepest_paralog_dup_support",
               "n_paralog_duplications", "note"], matrix_rows)

    if args.rooting_check:
        rows = rooting_check(PHYLO_DIR / "rooted.nwk", vert, sp, meta, group_of)
        best = min(r[4] for r in rows)
        for r in rows:
            r.append("minimum" if r[4] == best else "")
        write_tsv(out_dir / "rooting_check.tsv",
                  ["edge_index", "clade_size", "n_duplications",
                   "implied_losses", "total_events", "root_event",
                   "root_species_node", "deepest_paralog_dup_node",
                   "example_tips", "is_minimum"], rows)
        mins = [r for r in rows if r[-1] == "minimum"]
        out = rows[0]
        agree = {r[7] for r in mins} == {out[7]}
        print(f"\nrooting check: {len(rows)-1} edges + the outgroup rooting; "
              f"minimum total events {best} at {len(mins)} rooting(s); "
              f"outgroup rooting {out[4]}; deepest paralog duplication "
              f"{'agrees' if agree else 'DIFFERS'} "
              f"({out[7]} vs {sorted({r[7] for r in mins})})")

    print(f"\nwrote reconciliation tables to {out_dir}")


if __name__ == "__main__":
    main()
