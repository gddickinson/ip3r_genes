"""The figure-text checks every publication figure passes at save time (S28).

`figstyle.save()` calls `audit()` after the figure is drawn and refuses to
write a figure the audit fails. The checks exist because the fault they catch
was reported three times before it was checked mechanically (D78), and all
run against the *drawn* figure rather than the source, because a title built
from an f-string is not reachable by a regex over the module.

**R5 — no text overlaps other text, and no text is clipped.** Every visible
`Text` artist's extent is collected after a draw and every pair is tested for
intersection beyond a small tolerance. Rotated labels are tested as rotated
rectangles (separating-axis test), because the axis-aligned box of a 45°
tick label overlaps its neighbour's long before the glyphs do. A text whose
clip box cuts it, or an annotation whose anchor lies outside its axes and is
therefore not drawn at all, is reported as clipped.

**R5, for the panel furniture and the legends** — added after the first
page-by-page read of the typeset thesis found defects the text-against-text
test cannot see. A panel title may not run into another panel's area (one
was hidden behind the structure image beside it), must leave `TITLE_GAP_PT`
clear of the next panel's title and letter (three abutted with no visible
gap and no overlap), and when it wraps its letter must sit beside the first
line (`figstyle.panel` raises the letter to do so). A legend may not be drawn
over the data it keys: every line vertex (densified along each segment),
scatter offset and every bar (as an area, not its corners) of its axes is
tested against the legend's box, unlabelled reference lines included.

**R2 — no bare noun phrase where a reader reads a label.** A panel title, a
suptitle or a free annotation of `MIN_LABEL_WORDS` words or more must carry a
finite verb (`prose_lex.has_finite_verb`). Tick labels, axis labels and legend
entries are quantities and categories and are exempt; so is anything matching
`LABEL_EXEMPT`, each entry of which carries its reason.

Deliberate exceptions to R5 go in `ALLOWED_OVERLAPS` with a reason, keyed by
figure stem, in the same spirit as `s25_figmap.UNPLACED`.

Two environment variables exist for measurement, not for production:
`FIGSTYLE_LENIENT=1` reports every finding instead of raising, and
`FIGSTYLE_TEXT_LOG=<path>` appends one row per text artist to a TSV so the
starting state can be recorded before anything is edited.
"""

from __future__ import annotations

import math
import os
import re
from pathlib import Path

from matplotlib.text import Annotation, Text

import prose_lex as P

#: Two texts whose boxes overlap by less than this (in points) are touching,
#: not overlapping: a font box carries ascender and descender space that a
#: string of digits never inks, and the panel letter's descender space meets
#: the top tick label's ascender space on almost every axes.
TOL_PT = 1.2
#: A panel title and the next panel's letter or title closer than this (in
#: points) read as one run of text: the S28 page read found three titles
#: abutting the next letter with no visible gap, none of them overlapping.
TITLE_GAP_PT = 3.0
#: A free annotation with at least this many words is prose, not a label.
MIN_LABEL_WORDS = 5
#: A panel title with at least this many words is a description, not a name.
MIN_TITLE_WORDS = 3

#: (pattern, reason). A text matching the pattern is exempt from R2.
LABEL_EXEMPT: list[tuple[re.Pattern, str]] = [
    (re.compile(r"^(measured|computed|schematic|curated)( [·—] .*)?$"),
     "the provenance corner tag is a classification, not prose"),
    (re.compile(r"^[A-Za-z]$"), "a panel letter"),
    (re.compile(r"^(n|N|p|q|r|ρ|ω|J|AUC|TM)\s*[=<>≤≥]"),
     "a statistic printed beside its mark"),
]

#: figure stem -> [(text a, text b, reason)]; both texts matched by prefix.
ALLOWED_OVERLAPS: dict[str, list[tuple[str, str, str]]] = {}


# ------------------------------------------------------------- inventory

def _role_map(fig) -> dict[int, str]:
    roles: dict[int, str] = {}

    def tag(artists, role):
        for a in artists:
            roles.setdefault(id(a), role)

    for ax in fig.axes:
        tag([ax.title, ax._left_title, ax._right_title], "title")
        for axis in (ax.xaxis, ax.yaxis):
            # `ax.axis("off")` and `axis.set_visible(False)` leave every tick
            # label reporting itself visible while nothing draws it, and the
            # locator keeps Tick objects beyond the view limits whose labels
            # are never drawn either. `_update_ticks` is the list the axis
            # actually draws.
            off = not ax.axison or not axis.get_visible()
            tag([axis.label], "off" if off else "axis")
            tag([axis.get_offset_text()], "off" if off else "tick")
            drawn = axis._update_ticks() if not off else []
            for tick in drawn:
                tag([tick.label1, tick.label2], "tick")
            for tick in axis.majorTicks + axis.minorTicks:
                tag([tick.label1, tick.label2], "off")
        leg = ax.get_legend()
        if leg is not None:
            tag(leg.get_texts() + [leg.get_title()], "legend")
        tag(ax.texts, "text")
    for leg in fig.legends:
        tag(leg.get_texts() + [leg.get_title()], "legend")
    if fig._suptitle is not None:
        tag([fig._suptitle], "suptitle")
    tag(fig.texts, "figtext")
    return roles


def _visible(t: Text) -> bool:
    if not t.get_visible() or not t.get_text().strip():
        return False
    alpha = t.get_alpha()
    return alpha is None or alpha > 0


def _extent(t: Text, renderer):
    """The text's own box. `Annotation.get_window_extent` unions the arrow
    in, which would make every arrowed label overlap whatever it points at."""
    return Text.get_window_extent(t, renderer)


def _oriented(t: Text, renderer):
    """(cx, cy, half_w, half_h, angle_rad) of the drawn text, in pixels."""
    bb = _extent(t, renderer)
    angle = t.get_rotation() % 360.0
    if angle in (0.0, 90.0, 180.0, 270.0):
        return (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2, \
            bb.width / 2, bb.height / 2, 0.0
    saved = t.get_rotation()
    t.set_rotation(0)
    flat = _extent(t, renderer)
    t.set_rotation(saved)
    return (bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2, \
        flat.width / 2, flat.height / 2, math.radians(angle)


def _overlap(a, b, tol_px: float) -> bool:
    """Separating-axis test on two oriented rectangles shrunk by `tol_px`."""
    rects = []
    for cx, cy, hw, hh, ang in (a, b):
        hw, hh = max(hw - tol_px, 0.0), max(hh - tol_px, 0.0)
        c, s = math.cos(ang), math.sin(ang)
        corners = [(cx + x * c - y * s, cy + x * s + y * c)
                   for x, y in ((-hw, -hh), (hw, -hh), (hw, hh), (-hw, hh))]
        rects.append((corners, [(c, s), (-s, c)]))
    for _corners, axes in rects:
        for ax_x, ax_y in axes:
            proj = [[x * ax_x + y * ax_y for x, y in corners]
                    for corners, _ in rects]
            if max(proj[0]) <= min(proj[1]) or max(proj[1]) <= min(proj[0]):
                return False
    return True


def inventory(fig) -> list[dict]:
    """One row per visible text artist: role, text, extent, rotation."""
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    roles = _role_map(fig)
    rows = []
    for t in fig.findobj(Text):
        if not _visible(t):
            continue
        role = roles.get(id(t), "text")
        if role == "off":
            continue
        drawn = True
        if isinstance(t, Annotation) and not t._check_xy(renderer):
            drawn = False
        rows.append({"artist": t, "role": role, "drawn": drawn,
                     "text": " ".join(t.get_text().split()),
                     "rect": _oriented(t, renderer) if drawn else None})
    return rows


# ---------------------------------------------------------------- checks

def _exempt(text: str) -> bool:
    return any(pat.match(text) for pat, _r in LABEL_EXEMPT)


def _short(text: str, n: int = 50) -> str:
    return text if len(text) <= n else text[:n - 1] + "…"


def words(text: str) -> int:
    """Tokens carrying a letter. Numbers and symbols are not words, so a
    clade box reading `ITPR1 / n = 13 / 48/95` is a label, not a phrase."""
    return sum(1 for t in P._tokens(P.clean(text))
               if any(c.isalpha() for c in t))


def label_problems(rows: list[dict]) -> list[str]:
    out = []
    for r in rows:
        text, role = r["text"], r["role"]
        if role in ("tick", "axis", "legend") or _exempt(text):
            continue
        n_words = words(text)
        floor = MIN_TITLE_WORDS if role in ("title", "suptitle") \
            else MIN_LABEL_WORDS
        if n_words >= floor and not P.has_finite_verb(text):
            out.append(f"{role} has no finite verb: '{_short(text)}'")
    return out


def overlap_problems(fig, rows: list[dict], stem: str) -> list[str]:
    out = []
    tol_px = TOL_PT * fig.dpi / 72.0
    drawn = [r for r in rows if r["drawn"]]
    allowed = ALLOWED_OVERLAPS.get(stem, [])
    for i in range(len(drawn)):
        for j in range(i + 1, len(drawn)):
            a, b = drawn[i], drawn[j]
            if a["artist"] is b["artist"]:
                continue
            if not _overlap(a["rect"], b["rect"], tol_px):
                continue
            ta, tb = a["text"], b["text"]
            if any((ta.startswith(x) and tb.startswith(y)) or
                   (ta.startswith(y) and tb.startswith(x))
                   for x, y, _r in allowed):
                continue
            out.append(f"{a['role']} '{_short(ta)}' overlaps {b['role']} "
                       f"'{_short(tb)}' {_where(fig, a['rect'], b['rect'])}")
    for r in rows:
        if not r["drawn"]:
            out.append(f"{r['role']} '{_short(r['text'])}' is not drawn: its "
                       f"anchor lies outside the axes")
            continue
        t = r["artist"]
        if t.get_clip_on() and t.get_clip_box() is not None:
            bb = _extent(t, fig.canvas.get_renderer())
            cb = t.get_clip_box()
            if (bb.x0 < cb.x0 - tol_px or bb.x1 > cb.x1 + tol_px
                    or bb.y0 < cb.y0 - tol_px or bb.y1 > cb.y1 + tol_px):
                out.append(f"{r['role']} '{_short(r['text'])}' is clipped "
                           f"by its axes")
    return out


def _panel_letter(r: dict) -> bool:
    return r["role"] == "text" and bool(re.fullmatch(r"[a-z]", r["text"]))


def panel_problems(fig, rows: list[dict]) -> list[str]:
    """R5 for the panel furniture: a title stays on one line, inside its own
    column, clear of the next panel's letter and title."""
    out = []
    renderer = fig.canvas.get_renderer()
    gap_px = TITLE_GAP_PT * fig.dpi / 72.0
    titles = [r for r in rows if r["role"] == "title" and r["drawn"]]
    heads = titles + [r for r in rows if _panel_letter(r) and r["drawn"]]
    boxes = [(ax, ax.get_window_extent(renderer)) for ax in fig.axes
             if ax.get_visible()]
    for r in titles:
        t = r["artist"]
        bb = _extent(t, renderer)
        own = t.axes
        n_lines = t.get_text().strip().count("\n") + 1
        if n_lines > 1:
            line_h = bb.height / n_lines
            for o in heads:
                if o["artist"].axes is own and _panel_letter(o):
                    lb = _extent(o["artist"], renderer)
                    if lb.y1 < bb.y1 - line_h / 2:
                        out.append(f"title '{_short(r['text'])}' wraps and "
                                   f"its panel letter sits beside a lower "
                                   f"line")
        own_bb = own.get_window_extent(renderer) if own is not None else None
        for ax, ab in boxes:
            if ax is own or (own_bb is not None and
                             abs(ab.x0 - own_bb.x0) < 1 and
                             abs(ab.x1 - own_bb.x1) < 1):
                continue            # the title's own axes, or a twin of it
            if (bb.x1 > ab.x0 + gap_px and bb.x0 < ab.x1 and
                    bb.y1 > ab.y0 and bb.y0 < ab.y1):
                out.append(f"title '{_short(r['text'])}' runs into another "
                           f"panel's area")
                break
        for o in heads:
            if o is r or o["artist"].axes is own and _panel_letter(o):
                continue            # its own letter sits beside it by design
            if _overlap(r["rect"], o["rect"], -gap_px):
                out.append(f"title '{_short(r['text'])}' has less than "
                           f"{TITLE_GAP_PT:g} pt clear of {o['role']} "
                           f"'{_short(o['text'])}'")
    return out


def _mark_points(ax, renderer):
    """Display-space points the data marks of `ax` pass through: line
    vertices densified along each segment, scatter offsets, and the corners
    and centre of every bar."""
    import numpy as np
    pts = []
    for ln in ax.get_lines():
        if not ln.get_visible() or ln.get_zorder() <= 0:
            continue            # zorder <= 0 is a background guide, as a grid
                                # is; unlabelled reference lines still count
        xy = ln.get_transform().transform(ln.get_xydata())
        xy = xy[np.isfinite(xy).all(axis=1)]
        for a, b in zip(xy[:-1], xy[1:]):
            n = max(2, int(np.hypot(*(b - a)) // 3))
            pts.extend(a + (b - a) * f for f in np.linspace(0, 1, n))
        pts.extend(xy)
    for coll in ax.collections:
        if not coll.get_visible() or not hasattr(coll, "get_offsets"):
            continue
        offs = coll.get_offsets()
        if len(offs):
            pts.extend(coll.get_offset_transform().transform(offs))
    return pts


def _bar_boxes(ax, renderer):
    """Display-space extents of the bars drawn in `ax`."""
    return [p.get_window_extent(renderer) for p in ax.patches
            if p.get_visible() and type(p).__name__ == "Rectangle"
            and p.get_width() and p.get_height()]


def legend_problems(fig) -> list[str]:
    """R5 for legends: a legend is not drawn over the data it keys."""
    out = []
    renderer = fig.canvas.get_renderer()
    for ax in fig.axes:
        leg = ax.get_legend()
        if leg is None or not leg.get_visible():
            continue
        lb = leg.get_window_extent(renderer)
        hits = sum(1 for x, y in _mark_points(ax, renderer)
                   if lb.x0 < x < lb.x1 and lb.y0 < y < lb.y1)
        # a bar is an area: a legend inside its body touches no corner
        hits += sum(1 for e in _bar_boxes(ax, renderer)
                    if e.x0 < lb.x1 and e.x1 > lb.x0 and
                    e.y0 < lb.y1 and e.y1 > lb.y0)
        if hits:
            first = leg.get_texts()[0].get_text() if leg.get_texts() else ""
            out.append(f"legend '{_short(first)}' is drawn over {hits} "
                       f"point(s) of the data it keys")
    return out


def _where(fig, a, b) -> str:
    """The collision's position as a fraction of the canvas, x right and y
    up, so a reader can find it on the rendered figure."""
    r = fig.canvas.get_renderer()
    x = (a[0] + b[0]) / 2 / max(r.width, 1)
    y = (a[1] + b[1]) / 2 / max(r.height, 1)
    return f"(at x={x:.2f}, y={y:.2f} of the canvas)"


def audit(fig, stem) -> list[str]:
    """Every R2 and R5 finding for `fig`; logs the inventory if asked."""
    stem = Path(stem)
    rows = inventory(fig)
    problems = (overlap_problems(fig, rows, stem.name) + label_problems(rows)
                + panel_problems(fig, rows) + legend_problems(fig))
    log = os.environ.get("FIGSTYLE_TEXT_LOG")
    if log:
        with open(log, "a", encoding="utf-8") as fh:
            for r in rows:
                fh.write("\t".join([stem.name, r["role"], str(words(r["text"])),
                                    "1" if P.has_finite_verb(r["text"])
                                    else "0", "1" if r["drawn"] else "0",
                                    r["text"]]) + "\n")
            for p in problems:
                fh.write("\t".join([stem.name, "PROBLEM", "", "", "", p])
                         + "\n")
    return problems


def disabled() -> bool:
    """`FIGSTYLE_NOCHECK=1` skips the audit; for the byte-identity check."""
    return os.environ.get("FIGSTYLE_NOCHECK", "") not in ("", "0")


def lenient() -> bool:
    return os.environ.get("FIGSTYLE_LENIENT", "") not in ("", "0")
