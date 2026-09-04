"""One figure style for every publication figure in this project.

Import this and call `figstyle.use()` before building a figure. It fixes the
three things that made the figure set inconsistent:

**Size.** Figures used to be drawn 11-17 inches wide and then scaled to the
6.7-inch text block of the page, so 7 pt labels arrived at 3-4 pt. Everything is
now drawn at its final printed width (`W_FULL`), so a point is a point and the
PDF does no scaling.

**Colour.** One validated palette. Paralog identity is categorical and always
the same three hues; evidence and annotation quality are *ordered*, so they use
a diverging blue -> grey -> red scale instead of more categorical hues. The
values are the reference palette's slots 1/2/3 and its blue/red ramps; the trio
clears all-pairs CVD separation and each arm of the diverging scale is a
monotone single-hue ramp (checked with the palette validator, not by eye).

**Chrome.** Hairline solid axes, no top/right spines, recessive grid, sans face
matching the manuscript's headings, and text saved as text (`pdf.fonttype 42`).
"""

from __future__ import annotations

import subprocess
import warnings
from pathlib import Path

import matplotlib
from matplotlib import font_manager

# ----------------------------------------------------------------- geometry

# Page geometry has one definition, in the manuscript module that also places
# the figures; importing it here keeps drawing and placement in step.
try:
    from s14_lib import H_MAX, W_FULL, W_HALF
except ImportError:                                # standalone use
    W_FULL, W_HALF, H_MAX = 6.7, 3.25, 8.6

# ------------------------------------------------------------------- colour

#: Categorical: paralog identity. Reference palette slots 1, 2, 3 — the three
#: that clear the all-pairs CVD floor. Never reassign these.
PARALOG = {
    "ITPR1": "#2a78d6",   # blue
    "ITPR2": "#eb6834",   # orange
    "ITPR3": "#1baf7a",   # aqua
}
PARALOG_ORDER = ["ITPR1", "ITPR2", "ITPR3"]

#: The representative-alignment sequence groups. Categorical colour is spent on
#: paralog identity and nothing else, so the trio is `PARALOG`, the
#: non-vertebrate grade the tree nests it inside is neutral, and the sister
#: family — the ryanodine receptors, which every family-wide search returns and
#: which root the tree — gets the accent violet so it is never mistaken for an
#: ITPR. Used by the alignment figures and by every figure that colours a
#: representative-alignment row.
GROUP = {
    "ITPR1": PARALOG["ITPR1"], "ITPR2": PARALOG["ITPR2"],
    "ITPR3": PARALOG["ITPR3"], "invert_metazoa": "#52514e",
    "plant": "#a9a79e", "protist": "#c9c3b0", "fungi": "#8a897f",
    "RYR": "#4a3aa7",
}
GROUP_ORDER = ["ITPR1", "ITPR2", "ITPR3", "invert_metazoa", "plant",
               "protist", "fungi", "RYR"]
GROUP_LABEL = {"ITPR1": "ITPR1", "ITPR2": "ITPR2", "ITPR3": "ITPR3",
               "invert_metazoa": "invertebrates", "plant": "plants",
               "protist": "protists", "fungi": "fungi",
               "RYR": "ryanodine receptors"}

#: Diverging: how good the evidence for the gene is, present -> absent.
#: Blue arm (found), neutral greys (cannot tell), red arm (dead or gone).
STATUS = {
    "found_annotated": "#184f95",
    "found_unannotated": "#2a78d6",
    "found_no_annotation": "#86b6ef",
    "fragment": "#c9d8e8",
    "assembly_gap": "#d6d5cf",
    "tblastn_trace_ambiguous": "#a9a79e",
    "tblastn_trace": "#ef9a90",
    "absent": "#b3261e",
}
#: Reading order for stacked bars and legends: best evidence first.
STATUS_ORDER = [
    "found_annotated", "found_unannotated", "found_no_annotation",
    "fragment", "assembly_gap", "tblastn_trace_ambiguous", "tblastn_trace",
    "absent",
]
#: Short labels — the raw status strings are too long for a legend at 7 pt.
STATUS_LABEL = {
    "found_annotated": "found, annotated",
    "found_unannotated": "found, not annotated",
    "found_no_annotation": "found, no gene set",
    "fragment": "partial locus",
    "assembly_gap": "assembly gap",
    "tblastn_trace_ambiguous": "trace, ambiguous",
    "tblastn_trace": "remnant",
    "absent": "absent",
}

#: Diverging: annotation quality, same semantic direction as STATUS.
QUALITY = {
    "complete": "#184f95",
    "split": "#86b6ef",
    "fragmentary": "#d6d5cf",
    "noncoding": "#ef9a90",
    "unannotated": "#b3261e",
}

#: Clinical significance — a status encoding, never reused for a series.
#: The ITPR family has a real clinical variant set (SCA15/SCA29 and
#: Gillespie syndrome in ITPR1, anhidrosis in ITPR2, neuropathy in ITPR3),
#: so this scale carries weight rather than being decorative.
CLINICAL = {"pathogenic": "#b3261e", "uncertain": "#a9a79e",
            "benign": "#0f7d3d"}

#: Sequential blue ramp (light -> dark) for magnitude.
BLUES = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95"]

INK = "#14140f"        # primary text
MUTED = "#52514e"      # secondary text
FAINT = "#8a897f"      # annotation, de-emphasised marks
GRID = "#e5e4df"       # gridlines and hairline rules
SURFACE = "#ffffff"
ACCENT = "#4a3aa7"     # callouts, highlight boxes (violet, slot 7)
HILITE = "#f2f1ec"     # region shading behind a highlighted span

# ------------------------------------------------------------------ typography

FS_TICK = 6.6
FS_LABEL = 7.4
FS_TITLE = 8.2         # panel title
FS_LETTER = 9.0        # panel letter
FS_SUPTITLE = 9.4
FS_NOTE = 6.4          # in-plot annotation

_FONT_STACK = ["TeX Gyre Heros", "Helvetica Neue", "Helvetica", "Arial",
               "DejaVu Sans"]


def _register_document_sans() -> None:
    """Make the manuscript's sans face available to matplotlib.

    The PDF's headings are TeX Gyre Heros, which lives in the TeX tree and is
    not on matplotlib's font path. Registering it means figure text and heading
    text are the same typeface; if TeX is absent the stack falls back.
    """
    try:
        out = subprocess.run(
            ["kpsewhich", "texgyreheros-regular.otf", "texgyreheros-bold.otf",
             "texgyreheros-italic.otf", "texgyreheros-bolditalic.otf"],
            capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return
    for line in out.stdout.split("\n"):
        path = line.strip()
        if path and Path(path).exists():
            try:
                font_manager.fontManager.addfont(path)
            except Exception:                       # noqa: BLE001 - optional
                pass


_REGISTERED = False


def use() -> None:
    """Apply the project figure style. Safe to call more than once."""
    global _REGISTERED
    if not _REGISTERED:
        _register_document_sans()
        _REGISTERED = True
    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": _FONT_STACK,
        "font.size": FS_LABEL,
        "axes.titlesize": FS_TITLE,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.titlepad": 4.0,
        "axes.labelsize": FS_LABEL,
        "axes.labelcolor": INK,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.6,
        "axes.labelpad": 2.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "grid.linestyle": "-",
        "xtick.labelsize": FS_TICK,
        "ytick.labelsize": FS_TICK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.labelcolor": INK,
        "ytick.labelcolor": INK,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.4,
        "ytick.major.size": 2.4,
        "xtick.major.pad": 1.8,
        "ytick.major.pad": 1.8,
        "legend.fontsize": FS_TICK,
        "legend.frameon": False,
        "legend.handlelength": 1.1,
        "legend.handletextpad": 0.5,
        "legend.columnspacing": 1.0,
        "legend.labelspacing": 0.35,
        "legend.borderaxespad": 0.2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.0,
        "patch.linewidth": 0.0,
        "figure.titlesize": FS_SUPTITLE,
        "figure.titleweight": "bold",
        "figure.dpi": 120,
        "figure.facecolor": SURFACE,
        "savefig.dpi": 400,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "text.color": INK,
    })


# ------------------------------------------------------------------- helpers

def panel(ax, letter: str, title: str = "", pad: float = 4.0) -> None:
    """Bold panel letter at the top-left, description in roman beside it."""
    ax.set_title("  " + title if title else "", loc="left", pad=pad,
                 fontsize=FS_TITLE, color=INK)
    ax.annotate(letter, xy=(0.0, 1.0), xycoords="axes fraction",
                xytext=(-1.0, pad + 1.0), textcoords="offset points",
                fontsize=FS_LETTER, fontweight="bold", color=INK,
                ha="right", va="baseline", annotation_clip=False)


def despine(ax, keep=("left", "bottom")) -> None:
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def hgrid(ax, axis: str = "y") -> None:
    """A recessive solid grid, drawn behind the marks."""
    ax.set_axisbelow(True)
    ax.grid(True, axis=axis, color=GRID, linewidth=0.5, linestyle="-")


def paralog_handles(labels=None):
    """Legend handles for the paralog trio, in fixed order."""
    from matplotlib.patches import Patch
    labels = labels or PARALOG_ORDER
    return [Patch(facecolor=PARALOG[p], label=p) for p in labels]


def save(fig, stem, formats=("png", "pdf"), dpi: int = 400) -> list[Path]:
    """Save one figure to every format, and report if it exceeds the page."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    w, h = fig.get_size_inches()
    if w > W_FULL + 0.01:
        print(f"  [figstyle] {stem.name}: {w:.2f} in wide — wider than the "
              f"{W_FULL} in text block, it will be scaled down in the PDF")
    if h > H_MAX:
        print(f"  [figstyle] {stem.name}: {h:.2f} in tall — taller than "
              f"{H_MAX} in, it will be scaled down in the PDF")
    # A single artist placed in the wrong coordinate system can push the
    # tight bounding box hundreds of axes-heights off-canvas, and savefig
    # will happily write a several-hundred-megapixel file. Catch it here.
    try:
        fig.canvas.draw()
        tb = fig.get_tightbbox(fig.canvas.get_renderer())
        if tb.width > 3 * w or tb.height > 3 * h:
            raise RuntimeError(
                f"{stem.name}: tight bbox is {tb.width:.1f}x{tb.height:.1f} in "
                f"for a {w:.1f}x{h:.1f} in figure — an artist is drawn far "
                f"off-canvas (check transform= on any annotation)")
    except RuntimeError:
        raise
    except Exception:                              # noqa: BLE001 - best effort
        pass
    out = []
    for ext in formats:
        path = stem.with_suffix(f".{ext}")
        # A character the figure font lacks is dropped from the output with
        # only a warning, so a label can silently lose glyphs. Promote it.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            fig.savefig(path, dpi=dpi)
        missing = {str(w.message) for w in caught
                   if "missing from font" in str(w.message)}
        if missing:
            raise RuntimeError(f"{stem.name}: {len(missing)} glyph(s) have no "
                               f"outline in the figure font and would be "
                               f"dropped — {sorted(missing)[0]}")
        out.append(path)
    return out
