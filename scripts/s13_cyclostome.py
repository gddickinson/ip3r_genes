"""S13 step 4 — the six cyclostome loci, and what they are load-bearing for.

The whole difference between "the ITPR1 / ITPR2+ITPR3 duplication predates the
cyclostome-gnathostome split" and "it happened on the gnathostome stem" is six
tips: three *Myxine glutinosa* loci and three *Petromyzon marinus* loci, in two
cyclostome-only clades, one of them sister to the ITPR2+ITPR3 clade. S7 said so
and handed the question to S8; S8 measured it and came back **underpowered** —
550 My of rearrangement and annotations that name 27-29 % of coding genes leave
no shared flank vocabulary — and named the instrument that could answer it: "the
2R paralogon reconstructed from a cyclostome-anchored gene tree", which is this
task.

So this module does not re-ask S8's question. It asks the two questions that
decide how much weight the answer can carry:

**1. Does the tree's pairing of hagfish and lamprey loci survive an instrument
that never saw the alignment?** Each of the three cyclostome clades pairs one
*Myxine* locus with one *Petromyzon* locus — an orthology claim. S8's consensus
paralog caller read only the flanking gene symbols, so its per-locus call is
independent of every sequence in this project. The join is exact and offline:
the *Myxine* tips are S5 genome models and carry their own coordinates, and the
*Petromyzon* tips are UniProt records whose `gene` field in census v6 is the
same `LOC` symbol S8 recorded as `annot_gene`.

**2. Are these tips long branches?** Cyclostome sequences are the classic
long-branch attraction risk in vertebrate phylogeny, and long branches are
attracted to the root — which is exactly where the answer sits. A tip whose
root-to-tip distance is an outlier among the 57 is a tip whose deep placement
should be discounted, and the number is measured here rather than assumed
either way.

    python scripts/s13_cyclostome.py
"""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s6_lib import binomial                                        # noqa: E402
from s13_lib import (PHYLO_DIR, RECON_DIR, group_resolver,         # noqa: E402
                     load_tree, read_tsv, tip_metadata,
                     vertebrate_labels, write_tsv)

CENSUS = ROOT / "results" / "census_v6" / "census_v6.tsv"
S8_LOCI = ROOT / "results" / "synteny" / "loci.tsv"
S8_UNPLACED = ROOT / "results" / "synteny" / "unplaced_loci.tsv"
CYCLO_GENERA = ("Myxine", "Petromyzon")


def parse_model_accession(acc: str) -> tuple:
    """`GCF_x|ITPR1|NC_1:100-200-` -> ('NC_1', 100, 200); else ('', 0, 0)."""
    parts = acc.split("|")
    if len(parts) < 3 or ":" not in parts[2]:
        return ("", 0, 0)
    contig, span = parts[2].split(":", 1)
    span = span.rstrip("+-")
    try:
        lo, hi = (int(x) for x in span.split("-"))
    except ValueError:
        return ("", 0, 0)
    return (contig, lo, hi)


def census_gene(acc: str) -> str:
    for r in read_tsv(CENSUS):
        if r["accession"] == acc:
            return r.get("gene", "")
    return ""


def s8_rows() -> list:
    """Every S8 locus row for a cyclostome, with its caller verdict merged."""
    calls = {}
    for r in read_tsv(S8_UNPLACED):
        calls[(r["organism"], r["contig"], r["start"], r["end"])] = r
    out = []
    for r in read_tsv(S8_LOCI):
        if r["organism"].split()[0] not in CYCLO_GENERA:
            continue
        key = (r["organism"], r["contig"], r["start"], r["end"])
        r["_call"] = calls.get(key, {})
        out.append(r)
    return out


def join_tip(label: str, meta: dict, loci: list) -> dict:
    """Match a gene-tree tip to its S8 locus row. Exact or not at all."""
    row = meta[label]
    acc = row["accession"]
    contig, lo, hi = parse_model_accession(acc)
    if contig:
        for l in loci:
            if (l["contig"] == contig and int(l["start"]) == lo
                    and int(l["end"]) == hi):
                return {"how": "coordinates", "locus": l}
        return {"how": "no_match_by_coordinates", "locus": {}}
    gene = census_gene(acc)
    if gene:
        for l in loci:
            if l["annot_gene"] and l["annot_gene"] == gene:
                return {"how": f"census gene symbol {gene}", "locus": l}
        return {"how": f"no_match_for_gene_symbol {gene}", "locus": {}}
    return {"how": "no_join_key", "locus": {}}


def root_to_tip(tree) -> dict:
    """Path length from the tree's root to every tip."""
    out: dict = {}

    def rec(n, acc: float):
        acc += n.length
        if n.is_leaf:
            out[n.name] = acc
            return
        for c in n.children:
            rec(c, acc)

    for c in tree.children:                 # skip the root's own stub length
        rec(c, 0.0)
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=RECON_DIR)
    args = ap.parse_args()

    meta = tip_metadata()
    ml = load_tree(PHYLO_DIR / "rooted.nwk")
    group_of, _ = group_resolver(ml)
    vert = vertebrate_labels(meta)
    cyclo = sorted(l for l in vert
                   if meta[l]["species"].split()[0] in CYCLO_GENERA)
    loci = s8_rows()

    # --- which cyclostome clade each tip is in, from the tree, not by hand
    clades: dict = {}
    for n in ml.walk():
        if n.is_leaf:
            continue
        labels = set(n.leaf_names())
        if labels and labels <= set(cyclo):
            for lab in labels:
                # keep the *largest* pure-cyclostome clade holding the tip
                if lab not in clades or len(labels) > len(clades[lab][1]):
                    clades[lab] = (n, labels)

    # The orthology claim the tree makes is the *smallest* pure-cyclostome
    # clade holding both species — the hagfish/lamprey pair — not the largest
    # pure-cyclostome piece, which merges two such pairs and would report an
    # agreement between loci the tree never paired.
    pairs: dict = {}
    for n in ml.walk():
        if n.is_leaf:
            continue
        labels = set(n.leaf_names())
        if not labels or not labels <= set(cyclo):
            continue
        spp = {binomial(meta[x]["species"]) for x in labels}
        if len(spp) < 2:
            continue
        for lab in labels:
            if lab not in pairs or len(labels) < len(pairs[lab][1]):
                pairs[lab] = (n, labels)

    def s8_call_of(lab: str) -> str:
        j = join_tip(lab, meta, loci)
        return (j["locus"].get("_call", {}) or {}).get("call", "") if j["locus"] else ""

    rows = []
    for lab in cyclo:
        node, labels = clades.get(lab, (None, frozenset()))
        pnode, plabels = pairs.get(lab, (None, frozenset()))
        partners = sorted(plabels - {lab})
        pcalls = {s8_call_of(x) for x in plabels}
        if not pcalls or "" in pcalls or {"no_call"} & pcalls:
            agree = "uninformative"
        elif len(pcalls) == 1:
            agree = "agree"
        else:
            agree = "disagree"
        j = join_tip(lab, meta, loci)
        l = j["locus"]
        call = l.get("_call", {}) if l else {}
        rows.append([
            lab, binomial(meta[lab]["species"]), meta[lab]["accession"],
            group_of(lab),
            len(labels), node.name if node is not None else "",
            ";".join(sorted(binomial(meta[x]["species"]) for x in labels)),
            j["how"], l.get("cell", ""), l.get("annot_gene", ""),
            call.get("call", ""), call.get("best_score", ""),
            call.get("null_verdict", ""), call.get("reason", ""),
            len(plabels), pnode.name if pnode is not None else "",
            ";".join(p[:46] for p in partners), agree,
        ])
    write_tsv(args.out / "cyclostome_loci.tsv",
              ["label", "species", "accession", "tree_paralog",
               "cyclostome_clade_size", "clade_support", "clade_species",
               "s8_join", "s5_cell", "s8_annot_gene", "s8_paralog_call",
               "s8_call_score", "s8_null_verdict", "s8_reason",
               "pair_size", "pair_support", "pair_partner",
               "s8_pair_agreement"], rows)

    # --- long-branch check
    d = root_to_tip(ml)
    vals = sorted(d[l] for l in vert)
    med = statistics.median(vals)
    lb_rows = []
    for lab in sorted(vert):
        v = d[lab]
        rank = sum(1 for x in vals if x > v) + 1
        lb_rows.append([lab, group_of(lab), round(v, 4), round(v / med, 3),
                        rank, len(vals),
                        "cyclostome" if lab in cyclo else ""])
    lb_rows.sort(key=lambda r: -r[2])
    write_tsv(args.out / "branch_lengths.tsv",
              ["label", "tree_paralog", "root_to_tip", "vs_median", "rank",
               "n_tips", "is_cyclostome"], lb_rows)

    cy_ranks = [r[4] for r in lb_rows if r[6]]
    cy_ratio = [r[3] for r in lb_rows if r[6]]
    print(f"{len(cyclo)} cyclostome tips in {len({id(clades[l][0]) for l in cyclo if l in clades})} pure clades")
    for r in rows:
        print(f"  {r[0][:46]:46s} clade={r[4]} pair={r[14]}@{r[15] or '-':9s} "
              f"s8={r[10] or '-':7s} ({r[12] or 'n/a'}) -> {r[17]}")
    print(f"\nroot-to-tip: median {med:.3f}; cyclostome tips rank "
          f"{sorted(cy_ranks)} of {len(vals)}, "
          f"{min(cy_ratio):.2f}-{max(cy_ratio):.2f}x the median")
    print(f"wrote cyclostome_loci.tsv + branch_lengths.tsv to {args.out}")


if __name__ == "__main__":
    main()
