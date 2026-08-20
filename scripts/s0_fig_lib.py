"""Shared helpers for the review figures.

Every figure in `docs/figures/` is drawn by an `s0_figs_*.py` module, goes
through `figstyle`, and renders **only** from the committed tables under
`results/s0_baseline/review_figures/` (plus the S0/S1 tables that were already
committed). That is Decision D13 applied to figures: a figure cannot disagree
with the data behind it, because it has no other source.

This module holds what those modules share — table loading, the domain palette,
the domain-track primitive, and the provenance tag every caption carries.

**Provenance.** The review tags each statement `[db]` / `[lit]` / `[open]`; the
figures carry the same discipline in a visible corner tag:

    measured   drawn from a committed table of live database / structure output
    computed   derived here from committed sequences (alignment, identity)
    schematic  a drawing of a cited mechanism — no data, and not to scale
    curated    positions and categories transcribed from the cited literature

A schematic is not a lesser figure, but it must never be mistaken for a
measurement, so the tag is drawn on the canvas rather than left to the caption.
"""

from __future__ import annotations

import json
from pathlib import Path

import figstyle

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "results" / "s0_baseline" / "review_figures"
S0 = ROOT / "results" / "s0_baseline"
S1 = ROOT / "results" / "benchmark_controls"
FIG_DIR = ROOT / "docs" / "figures"

#: The **review's** text block: A4 less the 2.4 cm margins set in the LaTeX
#: preamble of `s0_review_build.py`. `figstyle.W_FULL` is the *manuscript's*
#: text width, and the review is a different document with different geometry
#: — drawing at 6.7 in and placing in a 6.38 in block would scale every label
#: down by 5%, which is exactly the failure figstyle exists to prevent.
W_REVIEW = (21.0 - 2 * 2.4) / 2.54

# ------------------------------------------------------------------ tables

def load_tsv(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    head = lines[0].split("\t")
    return [dict(zip(head, ln.split("\t")))
            for ln in lines[1:] if ln.strip()]


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def num(row: dict, key: str, cast=float):
    return cast(row[key])


# ------------------------------------------------------------------ colour

#: Domain colour encodes the §7.1 argument rather than domain identity for its
#: own sake: blues are the signatures IP3R shares with RyR, amber is the pore
#: domain the whole cation-channel world shares, violet is RyR-only — the same
#: violet `figstyle.GROUP` already reserves for the sister family.
DOMAIN_COLOUR = {
    "PF08709": "#184f95",
    "PF02815": "#2a78d6",
    "PF01365": "#6da7ec",
    "PF08454": "#b3d0f7",
    "PF00520": "#c08a2b",
    "PF02026": "#4a3aa7",
    "PF06459": "#6f62c4",
    "PF21119": "#948bd6",
    "PF00622": "#bab4e8",
    "PF13499": "#d5d1f0",
}
DOMAIN_SHORT = {
    "PF08709": "IP$_3$-binding core", "PF02815": "MIR", "PF01365": "RIH",
    "PF08454": "RIH-assoc", "PF00520": "pore", "PF02026": "RyR repeat",
    "PF06459": "RyR TM4-6", "PF21119": "RyR jsol", "PF00622": "SPRY",
    "PF13499": "EF-hand pair",
}
CLASS_LABEL = {
    "shared": "shared with the ryanodine receptors",
    "generic": "generic cation-channel pore",
    "ryr_only": "ryanodine-receptor–specific",
}
CLASS_SWATCH = {"shared": "#2a78d6", "generic": "#c08a2b",
                "ryr_only": "#4a3aa7"}

#: Disease mechanism — an ordered severity-free categorical set, kept distinct
#: from `figstyle.CLINICAL`, which encodes pathogenicity rather than mechanism.
MECHANISM_COLOUR = {
    "haploinsufficiency": "#2a78d6",
    "gain of function": "#c08a2b",
    "dominant negative": "#b3261e",
    "recessive loss": "#4a3aa7",
    "unresolved": "#8a897f",
}

PROV_COLOUR = {"measured": "#184f95", "computed": "#2a78d6",
               "schematic": "#8a897f", "curated": "#c08a2b"}


# ------------------------------------------------------------------ pieces

def provenance(fig, kind: str, note: str = "") -> None:
    """The corner tag. See the module docstring for what each word promises."""
    if kind not in PROV_COLOUR:
        raise ValueError(f"unknown provenance {kind!r}")
    txt = kind if not note else f"{kind} · {note}"
    # Anchored a fixed number of points *below* the canvas, so it clears the
    # lowest axis label whatever the figure height; the tight bounding box
    # then grows to include it. A figure-fraction offset would not: the same
    # fraction is a different gap on a 2-inch and a 5-inch figure.
    from matplotlib.transforms import offset_copy
    tr = offset_copy(fig.transFigure, fig=fig, x=0, y=-7, units="points")
    fig.text(0.995, 0.0, txt, ha="right", va="top", transform=tr,
             fontsize=figstyle.FS_NOTE - 0.4, color=PROV_COLOUR[kind],
             style="italic")


def domain_track(ax, domains: list[dict], y: float, length: int,
                 height: float = 0.30, label_min_aa: int = 130,
                 edge: str = "#ffffff") -> None:
    """One protein as a scaled bar with its Pfam domains drawn on it.

    `domains` are rows of `domain_coords.tsv`. The backbone is drawn full
    length so proteins of different sizes are visually comparable — the whole
    point of putting a 2,700 aa IP3R next to a 5,000 aa RyR.
    """
    ax.add_patch(_rect(0, y - height / 5, length, height * 2 / 5,
                       "#e5e4df", lw=0))
    for d in domains:
        s, e = int(d["start"]), int(d["end"])
        col = DOMAIN_COLOUR.get(d["pfam"], "#8a897f")
        ax.add_patch(_rect(s, y - height / 2, e - s, height, col,
                           lw=0.5, ec=edge))
        if e - s >= label_min_aa:
            ax.text((s + e) / 2, y, DOMAIN_SHORT.get(d["pfam"], d["pfam"]),
                    ha="center", va="center", fontsize=figstyle.FS_NOTE - 0.6,
                    color="#ffffff", zorder=5)


def _rect(x, y, w, h, colour, lw=0.0, ec="none"):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, h, facecolor=colour, linewidth=lw,
                     edgecolor=ec)


def domain_legend(ax, classes=("shared", "generic", "ryr_only"), **kw):
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=CLASS_SWATCH[c], label=CLASS_LABEL[c])
               for c in classes]
    return ax.legend(handles=handles, **kw)


def save(fig, slug: str):
    figstyle.use()
    return figstyle.save(fig, FIG_DIR / slug)
