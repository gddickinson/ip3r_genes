"""s15b_coding.py — re-running S15a's rule chain under a stated setting.

S15a assigned one state per genome x paralog cell with an ordered chain of
positive tests (`s15_states.py`, R0-R7).  S15b's deliverable is the
question of which *settings of that chain* manufacture a loss, so this
module re-runs the chain rather than relabelling its output — a setting
that disables R4 has to let the cell fall through to R5, R6 and R7 in
order, which a relabelling cannot do.

Nothing here re-measures anything.  Every input is a column S15a
committed: the ledger status, the locus coverage, the reconstruction
coverage, the genome's spare family loci, its contiguity against D4's bar
and its control.

**The evidence ladder is ordered and its rungs are named after what they
refuse**, from the loosest bar the calibration allows to accepting nothing
but a complete locus.  Two of the rungs are the calibration's own gap
edges rather than round numbers, because those are the two values the
S15a measurement actually licenses (D46) and a reader should be able to
see that moving the bar across its own uncertainty changes nothing.

`use_r5` and `use_r6` are the two *rules* the brief's axes name: D45's
`paralog_unassignable` (a cell the bait panel cannot fill is not an
absence) and D4's contiguity bar (an assembly that cannot represent the
gene is not evidence the gene is missing).  Turning either off is not a
sensitivity setting anybody should adopt; it is the measurement of what
that decision is worth, in losses.
"""

from __future__ import annotations

#: the reconstruction calibration's committed operating point and the two
#: edges of the gap it was placed in (S15a `recon_calibration`)
BAR_CALIBRATED = 0.13685
BAR_GAP_LO = 0.100
BAR_GAP_HI = 0.1737

PRESENT_STATES = ("present_single_locus", "present_truncated",
                  "present_partial", "present_fragmented")
UNDECIDED_STATES = ("paralog_unassignable", "undecidable_contiguity",
                    "no_control")

#: (name, recon bar or None to disable R4, the locus rules accepted)
EVIDENCE_LEVELS: tuple[tuple[str, float | None, tuple[str, ...]], ...] = (
    ("gap_lo", BAR_GAP_LO, ("R1", "R2", "R3")),
    ("calibrated", BAR_CALIBRATED, ("R1", "R2", "R3")),
    ("gap_hi", BAR_GAP_HI, ("R1", "R2", "R3")),
    ("recon_half", 0.50, ("R1", "R2", "R3")),
    ("recon_strict", 0.90, ("R1", "R2", "R3")),
    ("no_recon", None, ("R1", "R2", "R3")),
    ("locus_only", None, ("R1", "R2")),
    ("full_locus_only", None, ("R1",)),
)
EVIDENCE_NAMES = tuple(e[0] for e in EVIDENCE_LEVELS)

#: the setting S15a itself ran at, named once so the report cannot drift
BASE_SETTING = dict(evidence="calibrated", use_r5=True, use_r6=True)

_DESCRIPTIONS = {
    "gap_lo": "reconstruction bar at the calibration gap's lower edge "
              "(the decoy's maximum)",
    "calibrated": "S15a's operating point, the gap's midpoint",
    "gap_hi": "reconstruction bar at the gap's upper edge "
              "(the lowest candidate)",
    "recon_half": "reconstruction bar at half the reference",
    "recon_strict": "reconstruction bar at 0.90 of the reference",
    "no_recon": "a reference reassembled across contigs is not accepted "
                "as presence at all",
    "locus_only": "and a locus below the coverage bar with no assembly "
                  "excuse is not accepted either",
    "full_locus_only": "nothing but a locus at or above the sweep's own "
                       "coverage bar counts as presence",
}


def describe(evidence: str) -> str:
    return _DESCRIPTIONS[evidence]


def _f(v) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def recode(row: dict, evidence: str = "calibrated", use_r5: bool = True,
           use_r6: bool = True) -> dict:
    """One committed matrix row, re-run through the chain under a setting.

    Returns `state` and the `rule` that fired, exactly as `s15_states.py`
    would have written them had it run with these knobs.
    """
    bar, locus_rules = None, ()
    for name, b, rules in EVIDENCE_LEVELS:
        if name == evidence:
            bar, locus_rules = b, rules
            break
    else:
        raise ValueError(f"unknown evidence level {evidence!r}")

    status = row.get("ledger_status", "")
    if not int(_f(row.get("control_ok", 1))):
        return dict(state="no_control", rule="R0")
    if "R1" in locus_rules and status.startswith("found"):
        return dict(state="present_single_locus", rule="R1")
    if "R2" in locus_rules and status == "assembly_gap":
        return dict(state="present_truncated", rule="R2")
    if "R3" in locus_rules and status == "fragment":
        return dict(state="present_partial", rule="R3")
    if bar is not None and _f(row.get("recon_coverage")) >= bar:
        return dict(state="present_fragmented", rule="R4")
    if use_r5 and _f(row.get("n_spare_itpr_loci")) > 0:
        return dict(state="paralog_unassignable", rule="R5")
    if use_r6 and not int(_f(row.get("contig_spans_gene"))):
        return dict(state="undecidable_contiguity", rule="R6")
    return dict(state="absent", rule="R7")


def recode_all(rows: list[dict], evidence: str = "calibrated",
               use_r5: bool = True, use_r6: bool = True) -> list[dict]:
    out = []
    for r in rows:
        d = dict(r)
        d.update(recode(r, evidence, use_r5, use_r6))
        out.append(d)
    return out


# ------------------------------------------------------------- the codings

def paralog_states(rows: list[dict]) -> dict[str, dict[str, str]]:
    """`{cell: {accession: state}}` — the paralog-resolved coding."""
    out: dict[str, dict[str, str]] = {}
    for r in rows:
        out.setdefault(r["cell"], {})[r["accession"]] = r["state"]
    return out


def family_states(rows: list[dict]) -> dict[str, str]:
    """One state per genome: does the assembly carry the family at all?

    D46's coding, and the primary one.  A genome is `present` if any cell
    is present **or** it carries spare family loci no cell claimed — the
    spare loci are family evidence whatever the paralog chain does with
    them.  It is `absent` only if every cell reached `absent`; anything
    else is `undecided`, naming the state that stopped it.
    """
    by_genome: dict[str, list[dict]] = {}
    for r in rows:
        by_genome.setdefault(r["accession"], []).append(r)
    out = {}
    for acc, cells in by_genome.items():
        states = [c["state"] for c in cells]
        spare = max(_f(c.get("n_spare_itpr_loci")) for c in cells)
        if any(s in PRESENT_STATES for s in states) or spare > 0:
            out[acc] = "present"
        elif all(s == "absent" for s in states):
            out[acc] = "absent"
        else:
            out[acc] = "undecided"
    return out


def to_character(states: dict[str, str]) -> dict[str, int | None]:
    """State -> Dollo character: 1 present, 0 absent, None undecided."""
    out: dict[str, int | None] = {}
    for acc, st in states.items():
        if st in PRESENT_STATES or st == "present":
            out[acc] = 1
        elif st == "absent":
            out[acc] = 0
        else:
            out[acc] = None
    return out


def settings_grid(evidence_levels: tuple[str, ...] | None = None):
    """Every (evidence, use_r5, use_r6) the sensitivity matrix walks."""
    for ev in (evidence_levels or EVIDENCE_NAMES):
        for use_r5 in (True, False):
            for use_r6 in (True, False):
                yield dict(evidence=ev, use_r5=use_r5, use_r6=use_r6)
