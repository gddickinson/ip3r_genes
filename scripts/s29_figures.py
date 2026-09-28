"""S29's thesis figure, from committed tables only (D13): the pore across the
type-3 state panel (`s29_state_pores.py`).

Panel a overlays every deposition's pore profile with z measured from its
own filter, so the two constrictions line up and the one state that differs
is seen to differ at the gate. Panel b gives the two numbers per state,
ordered by the gate, on one radius axis, so the filter's constancy and the
gate's single jump are read off the same scale.

    python scripts/s29_figures.py
"""

from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
from matplotlib.lines import Line2D                            # noqa: E402

import figstyle as fs                                          # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "structures"
FIGS = OUT / "figures"
OPEN_STATE = "activated"


def _rows(name: str) -> list[dict]:
    with (OUT / name).open(encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def fig_state_pores() -> None:
    cons = sorted(_rows("state_constrictions.tsv"),
                  key=lambda r: float(r["gate_radius_A"]))
    prof = defaultdict(list)
    for r in _rows("state_pore_profiles.tsv"):
        prof[r["pdb_id"]].append((float(r["z_from_filter_A"]),
                                  float(r["min_heavy_atom_radius_A"])))
    closed = [r for r in cons if r["state"] != OPEN_STATE]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(fs.W_FULL, 3.1),
                                 width_ratios=[1.45, 1.0])
    for r in cons:
        is_open = r["state"] == OPEN_STATE
        z, rad = zip(*sorted(prof[r["pdb_id"]]))
        ax.plot(z, rad, color=fs.ACCENT if is_open else fs.BLUES[3],
                lw=1.4 if is_open else 0.8, zorder=3 if is_open else 2)
    ax.axvline(0, color=fs.FAINT, lw=0.6, ls=":", zorder=1)
    ax.text(0.8, 12.6, "filter", fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax.set_xlim(-25, 40)
    ax.set_ylim(0, 14)
    ax.set_xlabel("distance from the filter along the axis (Å)",
                  fontsize=fs.FS_LABEL)
    ax.set_ylabel("minimum heavy-atom radius (Å)", fontsize=fs.FS_LABEL)
    ax.legend(handles=[
        Line2D([], [], color=fs.ACCENT, lw=1.4, label="activated (8TKF)"),
        Line2D([], [], color=fs.BLUES[3], lw=0.8,
               label=f"the other {len(closed)} states")],
        loc="lower left", fontsize=fs.FS_TICK, frameon=False)
    fs.panel(ax, "a", "Only the activated pore widens past the filter")
    fs.despine(ax)

    y = range(len(cons))
    for yi, r in zip(y, cons):
        g, f = float(r["gate_radius_A"]), float(r["filter_radius_A"])
        bx.plot([g, f], [yi, yi], color=fs.GRID, lw=1.0, zorder=1)
        bx.plot(f, yi, "s", ms=4, color=fs.MUTED, zorder=2)
        bx.plot(g, yi, "o", ms=4.5, zorder=3,
                color=fs.ACCENT if r["state"] == OPEN_STATE else fs.BLUES[4])
    bx.set_yticks(list(y), [f"{r['pdb_id']} {r['state']}" for r in cons],
                  fontsize=fs.FS_TICK)
    bx.set_xlim(1.5, 6.5)
    bx.set_xlabel("radius at the constriction (Å)", fontsize=fs.FS_LABEL)
    bx.legend(handles=[
        Line2D([], [], marker="o", ls="none", color=fs.BLUES[4], ms=4.5,
               label="gate"),
        Line2D([], [], marker="s", ls="none", color=fs.MUTED, ms=4,
               label="filter")],
        loc="lower right", fontsize=fs.FS_TICK, frameon=False)
    fs.panel(bx, "b", "The gate moves and the filter does not")
    fs.despine(bx)
    fs.hgrid(bx, "x")
    fig.tight_layout()
    fs.save(fig, FIGS / "s29_state_pores")


def main() -> int:
    fs.use()
    fig_state_pores()
    print("[s29] figure state_pores")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
