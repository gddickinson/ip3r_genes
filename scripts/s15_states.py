"""s15_states.py — the character matrix's vocabulary, and the ordered rules
that assign one state to one genome x paralog cell.

A loss count is only as good as the state that licenses it, and the S5
ledger's statuses were not built to license one.  `absent` there means
"the rescue attributed no region to this cell", which in a shattered
assembly is a statement about contig lengths; `tblastn_trace` means "some
regions were attributed", which is not a presence call either.  So S15
re-states every cell, and the vocabulary is designed around one
constraint: **exactly one state may be counted as a loss**, and it is
reached only by passing every positive test that could explain the
absence otherwise.

The rules run in order and the first that fires wins.  Each writes the
number it fired on into the row, so a state is always readable back to
its evidence.

  R0 `no_control`            the genome's RyR positive control did not
                             fire, so the search is not demonstrated to
                             have reached this assembly and no cell in it
                             supports any claim (S5's own rule, D14).
  R1 `present_single_locus`  a locus at or above the sweep's own coverage
                             bar.
  R2 `present_truncated`     a locus below it whose truncation the
                             assembly explains — a contig edge or an
                             N-run under the model.
  R3 `present_partial`       a locus below it with no assembly excuse.
                             Present, and the shortfall is unexplained.
  R4 `present_fragmented`    no locus, but the reference reassembles from
                             sequence outside every locus the aligner
                             found, above the bar `s15_calibrate_recon.py`
                             measured against genes that are accounted
                             for elsewhere in the same genome.
  R5 `paralog_unassignable`  no locus and no reconstruction, but the
                             genome carries family loci no cell claimed.
                             D45: the cyclostomes read `absent` for ITPR2
                             and ITPR3 while carrying three ITPR loci
                             apiece, because the bait panel has no
                             cyclostome-labelled bait.  A count that read
                             those as losses would score two per genome
                             that never happened.
  R6 `undecidable_contiguity` nothing found, and the assembly's contig N50
                             is below D4's bar for this gene — an
                             assembly that cannot represent the gene
                             cannot be evidence that it is missing.
  R7 `absent`                nothing found, in a controlled assembly that
                             could have held the gene.  The only state
                             S15b may count.

`SYNTENY_STATES` is the second axis the brief asks for, and it is
deliberately not a duplicate of the first: it records what the
neighbourhood says, including — for 424 of 432 trace regions — that there
is no neighbourhood to ask.
"""

from __future__ import annotations

import s15_lib as lib
import s5_calibration as s5cal

#: the sweep's own coverage bar for a full locus (s5_classify.COV_FOUND)
COV_FOUND = 0.70

STATES = ("present_single_locus", "present_truncated", "present_partial",
          "present_fragmented", "paralog_unassignable",
          "undecidable_contiguity", "absent", "no_control")

#: the states that assert the gene is in the assembly
PRESENT = ("present_single_locus", "present_truncated", "present_partial",
           "present_fragmented")
#: the states that assert nothing either way
UNDECIDED = ("paralog_unassignable", "undecidable_contiguity", "no_control")

SYNTENY_STATES = ("corroborated", "contradicted", "no_call",
                  "no_neighbourhood", "not_asked")


def spare_loci(summary: dict) -> int:
    """ITPR loci in this genome beyond the number of cells they filled.

    S13's definition (`s13_losses.ledger_by_species`), recomputed here
    from the per-genome summary rather than the ledger because the ledger
    holds only each cell's best locus.  The count is what separates a
    bait-panel limit from a biological absence, and the cyclostomes are
    why it exists: each carries **three** ITPR loci, all of them filed in
    the ITPR1 cell as secondaries because the panel has no
    cyclostome-labelled bait, and both `absent` cells sit beside two spare
    loci.  `other_loci` — the loci no cell claimed at all — is added,
    since it is the same evidence recorded in a second place.
    """
    cells = summary.get("cells") or {}
    total = sum(len(cells.get(c, {}).get("loci") or []) for c in lib.ITPR_CELLS)
    filled = sum(1 for c in lib.ITPR_CELLS
                 if (cells.get(c, {}).get("loci") or []))
    unclaimed = sum(1 for loc in (summary.get("other_loci") or [])
                    if str(loc.get("family") or loc.get("locus_family") or "")
                    .upper().startswith("ITPR"))
    return max(0, total - filled) + unclaimed


def classify(cell_row: dict, summary: dict, recon_row: dict | None,
             bar: float) -> dict:
    """One cell, one state, with the number the rule fired on."""
    acc = cell_row["accession"]
    n50 = int(cell_row.get("contig_n50") or 0)
    status = cell_row["status"]
    cov = _f(cell_row.get("best_coverage"))
    control_ok = bool(summary.get("control_ok"))
    spare = spare_loci(summary)
    contiguous = s5cal.spans_a_gene(n50)
    rc = _f(recon_row.get("coverage")) if recon_row else 0.0
    rg = _f(recon_row.get("gene_equiv")) if recon_row else 0.0
    n_contigs = int(recon_row.get("n_contigs") or 0) if recon_row else 0

    def out(state, rule, why):
        return dict(accession=acc, organism=cell_row["organism"],
                    vclass=cell_row["vclass"], cell=cell_row["class"],
                    ledger_status=status, state=state, rule=rule, reason=why,
                    best_coverage=cov, contig_n50=n50,
                    contig_spans_gene=int(contiguous),
                    control_ok=int(control_ok), n_spare_itpr_loci=spare,
                    recon_coverage=round(rc, 4),
                    recon_gene_equiv=round(rg, 4),
                    recon_n_contigs=n_contigs)

    if not control_ok:
        return out("no_control", "R0",
                   f"the RyR positive control did not fire in {acc}")
    if status.startswith("found"):
        return out("present_single_locus", "R1",
                   f"one locus at coverage {cov:.2f} >= {COV_FOUND}")
    if status == "assembly_gap":
        return out("present_truncated", "R2",
                   f"locus at coverage {cov:.2f}, truncated at a contig "
                   f"edge or an N-run")
    if status == "fragment":
        return out("present_partial", "R3",
                   f"locus at coverage {cov:.2f} with no assembly excuse")
    if rc >= bar:
        return out("present_fragmented", "R4",
                   f"{rc:.2f} of the reference reassembled from "
                   f"{n_contigs} contigs outside every locus the aligner "
                   f"found, bar {bar:.3f}")
    if spare > 0:
        return out("paralog_unassignable", "R5",
                   f"no locus in this cell, but {spare} family loci in "
                   f"{acc} that no cell claimed")
    if not contiguous:
        return out("undecidable_contiguity", "R6",
                   f"contig N50 {n50:,} bp is below D4's bar "
                   f"({int(s5cal.itpr_span_stats()['median']):,} bp) for "
                   f"this gene")
    return out("absent", "R7",
               f"nothing found in a controlled assembly with contig N50 "
               f"{n50:,} bp, reconstruction {rc:.3f} < bar {bar:.3f}")


def synteny_state(cell: str, calls: list[dict],
                  reach_rows: list[dict] | None = None) -> dict:
    """The second axis: what the neighbourhood says about this cell.

    `reach_rows` is what stops a cell the caller could not reach from
    reading as a cell nobody asked about.  A cell with rescue regions and
    no reachable one is `no_neighbourhood` — the measurement that answers
    the brief's first step — and only a cell with no regions at all
    (because the aligner placed its gene) is `not_asked`.
    """
    reach_rows = reach_rows or []
    if not calls:
        if not reach_rows:
            return dict(synteny_state="not_asked", synteny_call="",
                        synteny_n_regions=0, synteny_n_keys=0,
                        synteny_reason="the aligner placed this gene; "
                                       "synteny was not asked to")
        keys = [int(r.get("n_keys") or 0) for r in reach_rows]
        state = "no_neighbourhood" if max(keys, default=0) < 4 else "no_call"
        return dict(
            synteny_state=state, synteny_call="",
            synteny_n_regions=len(reach_rows), synteny_n_keys=max(keys or [0]),
            synteny_reason=(
                f"{sum(1 for k in keys if k == 0)} of {len(keys)} regions "
                f"have no flanking gene at all; the best has {max(keys or [0])} "
                f"informative neighbours against the {4} the caller needs"))
    called = [c for c in calls if c["call"] != "no_call"]
    if not called:
        return dict(synteny_state="no_call", synteny_call="",
                    synteny_n_regions=len(calls),
                    synteny_n_keys=max((int(c.get("n_keys") or 0)
                                        for c in calls), default=0),
                    synteny_reason="every reachable region returned no_call")
    agree = [c for c in called if c["call"] == cell]
    state = "corroborated" if len(agree) >= len(called) - len(agree) \
        else "contradicted"
    return dict(synteny_state=state,
                synteny_call=";".join(sorted({c["call"] for c in called})),
                synteny_n_regions=len(calls),
                synteny_n_keys=max(int(c.get("n_keys") or 0) for c in called),
                synteny_reason=f"{len(agree)} of {len(called)} reachable "
                               f"regions call {cell}")


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


COLS = ["accession", "organism", "vclass", "cell", "ledger_status", "state",
        "rule", "reason", "best_coverage", "contig_n50",
        "contig_spans_gene", "control_ok", "n_spare_itpr_loci",
        "recon_coverage", "recon_gene_equiv", "recon_n_contigs",
        "synteny_state", "synteny_call", "synteny_n_regions",
        "synteny_n_keys", "synteny_reason"]
