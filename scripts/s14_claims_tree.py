"""S14 claims: what the phylogeny figure draws (added 2026-09-23).

The tree figure boxes each paralogue's whole clade, and its legend states the
clade sizes, the unlabelled members the tree places inside them, the support
the labelled core has on its own, and where the six jawless-fish genes sit.
Each of those numbers is declared here once, against the committed table it
comes from. The thesis and the origin paper carry them by identifier.
Imported and concatenated by `s14_claims.py`.
"""

from __future__ import annotations

_NODES = "results/phylogeny/claim_nodes.tsv"
_CLADES = "results/phylogeny/paralog_clades.tsv"
_CYC = "results/phylogeny/cyclostome_placement.tsv"
_EXT = ("ITPR1 clade incl. unlabelled tips,ITPR2 clade incl. unlabelled "
        "tips,ITPR3 clade incl. unlabelled tips")

CLAIMS: list[dict] = [
    dict(id="C277", section="Fig3",
         claim="the whole ITPR1 clade holds 19 proteins",
         source=_NODES, op="cell", column="n_tips",
         where={"claim": "ITPR1 clade incl. unlabelled tips"}, expect="19"),
    dict(id="C278", section="Fig3",
         claim="the whole ITPR2 clade holds 13 proteins",
         source=_NODES, op="cell", column="n_tips",
         where={"claim": "ITPR2 clade incl. unlabelled tips"}, expect="13"),
    dict(id="C279", section="Fig3",
         claim="the whole ITPR3 clade holds 19 proteins",
         source=_NODES, op="cell", column="n_tips",
         where={"claim": "ITPR3 clade incl. unlabelled tips"}, expect="19"),
    dict(id="C280", section="Fig3",
         claim="all three whole paralogue clades are supported at 100/100",
         source=_NODES, op="count",
         where={"claim__in": _EXT, "alrt": "100", "ufboot": "100"},
         expect="3"),
    dict(id="C281", section="Fig3",
         claim="seven unlabelled vertebrate proteins are placed inside a "
               "paralogue clade by the tree",
         source=_CLADES, op="sum", column="n_added", expect="7"),
    dict(id="C282", section="Fig3",
         claim="the name-labelled ITPR1 core clade alone holds 13 proteins",
         source=_NODES, op="cell", column="n_tips",
         where={"claim": "ITPR1 core clade"}, expect="13"),
    dict(id="C283", section="Fig3",
         claim="the name-labelled ITPR1 core clade has SH-aLRT 47.8",
         source=_NODES, op="cell", column="alrt",
         where={"claim": "ITPR1 core clade"}, expect="47.8"),
    dict(id="C284", section="Fig3",
         claim="all 57 vertebrate proteins form one clade",
         source=_NODES, op="cell", column="n_tips",
         where={"claim": "vertebrate ITPRs (incl. unassigned)"},
         expect="57"),
    dict(id="C285", section="Fig3",
         claim="four hagfish and lamprey genes form their own clade",
         source=_CYC, op="count", where={"clade_size": "4"}, expect="4"),
    dict(id="C286", section="Fig3",
         claim="that clade has SH-aLRT 99.5",
         source=_CYC, op="nunique", column="alrt",
         where={"clade_size": "4", "alrt": "99.5"}, expect="1"),
    dict(id="C287", section="Fig3",
         claim="two hagfish and lamprey genes join the ITPR2 + ITPR3 clade "
               "at SH-aLRT 83.1",
         source=_CYC, op="count", where={"clade_size": "34", "alrt": "83.1",
                                         "ufboot": "77"}, expect="2"),
]
