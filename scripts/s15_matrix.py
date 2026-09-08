"""s15_matrix.py — the character matrix S15b counts losses on.

One row per genome x paralog cell, carrying both evidence axes the brief
asks for side by side: the alignment-derived state from `s15_states.py`
and the neighbourhood-derived state from `s15_synteny_reach.py`.  Two
axes, never merged into one, because they answer different questions and
— measured here — one of them almost never arrives.

`implied_copies()` is the per-genome statistic that makes the matrix
countable.  Reference *coverage* cannot count copies: the paralogs are
61-68 % identical, so a genome holding only ITPR1 recovers most of the
ITPR2 reference too, which is exactly why `s15_calibrate_recon.py` needed
a decoy at all.  Genomic sequence can: a locus and a reconstruction
occupy different places in the assembly, so the per-cell contributions
add.  A cell contributes 1.0 if the aligner placed a whole gene, its own
coverage if it placed a partial one, and its reconstruction's
gene-equivalents if it placed none — and a genome's total is what a copy
number would be if the assembly were contiguous.

`wide()` is the same matrix one row per genome, which is the shape a
Dollo pass reads, and `state_counts()` is the summary the report renders
from (D13: nothing downstream recomputes a state).
"""

from __future__ import annotations

import collections

import s15_lib as lib
import s15_states as st


def build(ledger: list[dict], summaries: dict[str, dict],
          recon_rows: list[dict], synteny_calls: list[dict],
          reach_rows: list[dict], bar: float,
          window: str = "informative10") -> list[dict]:
    recon = {(r["accession"], r["cell"]): r for r in recon_rows
             if r["scope"] == "own_clade"}
    calls: dict[tuple, list[dict]] = collections.defaultdict(list)
    for c in synteny_calls:
        calls[(c["accession"], c["cell"])].append(c)
    reach: dict[tuple, list[dict]] = collections.defaultdict(list)
    for r in reach_rows:
        if r.get("window") == window:
            reach[(r["accession"], r["cell"])].append(r)
    rows: list[dict] = []
    for r in ledger:
        if r["class"] not in lib.ITPR_CELLS:
            continue
        s = summaries.get(r["accession"])
        if s is None:
            continue
        row = st.classify(r, s, recon.get((r["accession"], r["class"])), bar)
        key = (r["accession"], r["class"])
        row.update(st.synteny_state(r["class"], calls.get(key, []),
                                    reach.get(key, [])))
        rows.append(row)
    rows.sort(key=lambda r: (r["vclass"], r["organism"], r["cell"]))
    return rows


def implied_copies(matrix: list[dict]) -> list[dict]:
    """Per genome: how many ITPR genes' worth of sequence the assembly holds.

    Additive because the contributions are disjoint pieces of genome, not
    overlapping pieces of a reference.
    """
    by_acc: dict[str, list[dict]] = collections.defaultdict(list)
    for r in matrix:
        by_acc[r["accession"]].append(r)
    out = []
    for acc, cells in by_acc.items():
        contrib = {}
        for r in cells:
            if r["state"] == "present_single_locus":
                contrib[r["cell"]] = 1.0
            elif r["state"] in ("present_truncated", "present_partial"):
                contrib[r["cell"]] = float(r["best_coverage"])
            elif r["state"] == "present_fragmented":
                contrib[r["cell"]] = float(r["recon_gene_equiv"])
            else:
                contrib[r["cell"]] = 0.0
        any_cell = cells[0]
        # the spare loci are ITPR genes the bait panel could not file into a
        # cell (D45) — real copies, and leaving them out would report the
        # cyclostomes as carrying one gene when the sweep placed three
        spare = int(any_cell["n_spare_itpr_loci"] or 0)
        out.append(dict(
            accession=acc, organism=any_cell["organism"],
            vclass=any_cell["vclass"],
            contig_n50=any_cell["contig_n50"],
            contig_spans_gene=any_cell["contig_spans_gene"],
            control_ok=any_cell["control_ok"],
            n_spare_itpr_loci=spare,
            implied_copies=round(sum(contrib.values()) + spare, 3),
            n_present=sum(1 for r in cells if r["state"] in st.PRESENT),
            n_undecided=sum(1 for r in cells if r["state"] in st.UNDECIDED),
            n_absent=sum(1 for r in cells if r["state"] == "absent"),
            **{f"copies_{c}": round(contrib.get(c, 0.0), 3)
               for c in lib.ITPR_CELLS}))
    out.sort(key=lambda r: (r["vclass"], r["organism"]))
    return out


def wide(matrix: list[dict]) -> list[dict]:
    by_acc: dict[str, dict] = collections.defaultdict(dict)
    meta: dict[str, dict] = {}
    for r in matrix:
        by_acc[r["accession"]][r["cell"]] = r
        meta[r["accession"]] = r
    out = []
    for acc, cells in by_acc.items():
        m = meta[acc]
        row = dict(accession=acc, organism=m["organism"], vclass=m["vclass"],
                   contig_n50=m["contig_n50"],
                   contig_spans_gene=m["contig_spans_gene"],
                   control_ok=m["control_ok"])
        for c in lib.ITPR_CELLS:
            cr = cells.get(c)
            row[f"state_{c}"] = cr["state"] if cr else ""
            row[f"synteny_{c}"] = cr["synteny_state"] if cr else ""
        out.append(row)
    out.sort(key=lambda r: (r["vclass"], r["organism"]))
    return out


def state_counts(matrix: list[dict]) -> list[dict]:
    rows = []
    for scope, key in (("all", lambda r: "all"),
                       ("by_cell", lambda r: r["cell"]),
                       ("by_class", lambda r: r["vclass"])):
        c = collections.Counter((key(r), r["state"]) for r in matrix)
        for (grp, state), n in sorted(c.items()):
            rows.append(dict(scope=scope, group=grp, state=state, n=n))
    for (grp, state), n in sorted(collections.Counter(
            (r["cell"], r["synteny_state"]) for r in matrix).items()):
        rows.append(dict(scope="synteny_by_cell", group=grp, state=state, n=n))
    return rows


def loss_candidates(matrix: list[dict]) -> list[dict]:
    """Every cell that reaches `absent`, and every cell that nearly did.

    The near-misses are committed beside the losses because the near-miss
    list *is* the sensitivity of the count: each row names the rule that
    stopped it, so S15b's sensitivity matrix can be read as "which of
    these would a different filter have promoted".
    """
    out = []
    for r in matrix:
        if r["state"] == "absent":
            out.append(dict(r, candidate_kind="loss"))
        elif r["state"] in st.UNDECIDED:
            out.append(dict(r, candidate_kind="blocked_by_" + r["rule"]))
    return out


COPY_COLS = ["accession", "organism", "vclass", "contig_n50",
             "contig_spans_gene", "control_ok", "n_spare_itpr_loci",
             "implied_copies",
             "copies_ITPR1", "copies_ITPR2", "copies_ITPR3",
             "n_present", "n_undecided", "n_absent"]
WIDE_COLS = ["accession", "organism", "vclass", "contig_n50",
             "contig_spans_gene", "control_ok",
             "state_ITPR1", "state_ITPR2", "state_ITPR3",
             "synteny_ITPR1", "synteny_ITPR2", "synteny_ITPR3"]
COUNT_COLS = ["scope", "group", "state", "n"]
CAND_COLS = st.COLS + ["candidate_kind"]
