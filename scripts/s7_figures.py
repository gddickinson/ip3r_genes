"""S7 figures 2–4, from the committed S7 tables only (D13, D19).

Figure 1 (the rooted phylogram) is `s7_figure.py` — it is large enough to
own a module. This one draws the three that read off tables:

  sister_au.png         the AU test: how much likelihood each constrained
                        topology costs, and what the test does to it
  support_profile.png   the joint SH-aLRT × UFBoot distribution over every
                        internal node, with the claim nodes marked
  paralog_placement.png where the vertebrate tips that carry no paralog
                        label resolve, and at what support

Form follows the job, not habit. The AU panel does **not** put ΔlogL and
p-AU on two y-scales — they are different measures on different scales, so
they are two panels sharing one row of categories. The support figure is a
scatter rather than two histograms because the claim is about the *joint*
condition (both thresholds), which a pair of marginals cannot show.

Run:  python3 scripts/s7_figures.py [--only <slug>]
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from s6_lib import read_tsv                                   # noqa: E402
from s7_lib import (                                          # noqa: E402
    MSA_DIR, PHYLO_DIR, parse_newick, support_of,
)

FIGS = PHYLO_DIR / "figures"
AU_ALPHA = 0.05
MIN_ALRT, MIN_UFBOOT = 80.0, 95.0

#: What a tip resolves with when no single paralog owns its
#: neighbourhood — a category, not a missing value.
NO_HOME = "no single paralog"

PAIR_LABEL = {"H1_12": "ITPR1 + ITPR2", "H2_13": "ITPR1 + ITPR3",
              "H3_23": "ITPR2 + ITPR3", "ML": "unconstrained ML"}


def load(name: str) -> list[dict]:
    p = PHYLO_DIR / name
    return read_tsv(p) if p.exists() else []


def _f(v, default=float("nan")) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# --------------------------------------------------------------- fig 2

def fig_sister_au(fs, plt) -> str | None:
    rows = load("au_test.tsv")
    if not rows:
        return None
    rows = sorted(rows, key=lambda r: _f(r["deltaL"], 0))
    names = [PAIR_LABEL.get(r["tree"], r["tree"]) for r in rows]
    y = list(range(len(rows)))[::-1]
    dl = [_f(r["deltaL"], 0.0) for r in rows]
    au = [_f(r["p_AU"]) for r in rows]
    rejected = [a < AU_ALPHA for a in au]

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(fs.W_FULL, 0.42 * len(rows) + 1.35),
        gridspec_kw={"width_ratios": [1.15, 1.0], "wspace": 0.08})

    # A — likelihood cost of each constrained topology
    fs.hgrid(ax1, "x")
    ax1.barh(y, dl, height=0.5, color=fs.INK, zorder=3)
    for yi, v in zip(y, dl):
        ax1.annotate(f"{v:,.1f}", (v, yi), xytext=(4, 0),
                     textcoords="offset points", va="center", ha="left",
                     fontsize=fs.FS_NOTE, color=fs.MUTED)
    ax1.set_yticks(y, names, fontsize=fs.FS_TICK)
    ax1.set_xlabel("log-likelihood cost against the best tree (ΔlogL)",
                   fontsize=fs.FS_LABEL, color=fs.MUTED)
    ax1.set_xlim(0, max(max(dl) * 1.28, 1e-9))
    fs.despine(ax1)
    fs.panel(ax1, "a", "what each topology costs")

    # B — the test's own verdict, on its own scale
    fs.hgrid(ax2, "x")
    ax2.axvline(AU_ALPHA, color=fs.CLINICAL["pathogenic"], lw=0.9,
                zorder=2)
    ax2.annotate(f"p-AU = {AU_ALPHA}\nrejection", (AU_ALPHA, len(rows) - 0.35),
                 xytext=(4, 0), textcoords="offset points", va="top",
                 ha="left", fontsize=fs.FS_NOTE,
                 color=fs.CLINICAL["pathogenic"])
    for yi, a, rej in zip(y, au, rejected):
        col = fs.CLINICAL["pathogenic"] if rej else fs.INK
        ax2.plot([0, a], [yi, yi], color=fs.GRID, lw=1.4, zorder=2)
        ax2.plot([a], [yi], "o", ms=5.0, color=col, mew=0, zorder=4)
        # A rejected hypothesis sits at p-AU ~ 0, which is *left* of the
        # 0.05 line, so a label offset from the point starts underneath
        # the line and loses its first character — "1.24e-51" reads as
        # ".24e-51", which is a different number. Anchor those labels
        # past the line instead.
        ax2.annotate(f"{a:.3g}  {'rejected' if rej else 'not rejected'}",
                     (max(a, AU_ALPHA), yi), xytext=(7, 0),
                     textcoords="offset points",
                     va="center", ha="left", fontsize=fs.FS_NOTE,
                     color=col)
    ax2.set_yticks(y, ["" for _ in y])
    ax2.set_xlim(0, 1.0)
    ax2.set_xlabel("p-AU (10,000 RELL replicates)",
                   fontsize=fs.FS_LABEL, color=fs.MUTED)
    fs.despine(ax2, keep=("bottom",))
    fs.panel(ax2, "b", "what the AU test does to it")

    fig.suptitle("The three sister hypotheses for ITPR1/2/3, tested",
                 fontsize=fs.FS_SUPTITLE, fontweight="bold", color=fs.INK,
                 x=0.005, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fs.save(fig, FIGS / "sister_au")
    plt.close(fig)
    return "sister_au"


# --------------------------------------------------------------- fig 3

def fig_support_profile(fs, plt) -> str | None:
    tree_path = PHYLO_DIR / "itpr_ml.treefile"
    if not tree_path.exists():
        return None
    tree = parse_newick(tree_path.read_text())
    pts = []
    for n in tree.walk():
        if n.is_leaf or n is tree:
            continue
        a, u = support_of(n.name)
        if a is not None and u is not None:
            pts.append((a, u))
    if not pts:
        return None
    claims = [c for c in load("claim_nodes.tsv")
              if c["is_clade"] == "yes" and c["alrt"] and c["ufboot"]]

    fig, ax = plt.subplots(figsize=(fs.W_HALF * 1.55, 2.95))
    fs.hgrid(ax, "both")
    ax.axvspan(MIN_ALRT, 100, ymin=0, ymax=1, color=fs.GRID, alpha=0.45,
               zorder=0, lw=0)
    ax.axhline(MIN_UFBOOT, color=fs.FAINT, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.axvline(MIN_ALRT, color=fs.FAINT, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=3.0, mew=0,
            color=fs.MUTED, alpha=0.55, zorder=3, label="internal node")
    ax.plot([_f(c["alrt"]) for c in claims], [_f(c["ufboot"]) for c in claims],
            "o", ms=6.0, mfc="none", mec=fs.PARALOG["ITPR1"], mew=1.2,
            zorder=5, label="node a claim rests on")
    ok = sum(1 for a, u in pts if a >= MIN_ALRT and u >= MIN_UFBOOT)
    ax.annotate(f"both thresholds\n{ok} of {len(pts)} nodes "
                f"({100.0 * ok / len(pts):.0f} %)",
                xy=(0.985, 0.045), xycoords="axes fraction", ha="right",
                va="bottom", fontsize=fs.FS_NOTE, color=fs.MUTED,
                multialignment="right")
    ax.set_xlabel("SH-aLRT (%)", fontsize=fs.FS_LABEL, color=fs.MUTED)
    ax.set_ylabel("ultrafast bootstrap (%)", fontsize=fs.FS_LABEL,
                  color=fs.MUTED)
    ax.set_xlim(-3, 103)
    ax.set_ylim(-3, 103)
    fs.despine(ax)
    ax.legend(loc="lower left", frameon=False, fontsize=fs.FS_NOTE,
              handletextpad=0.4, borderaxespad=0.3)
    ax.set_title("How much of the tree is resolved", loc="left",
                 fontsize=fs.FS_SUPTITLE, fontweight="bold", color=fs.INK,
                 pad=5)
    fig.tight_layout()
    fs.save(fig, FIGS / "support_profile")
    plt.close(fig)
    return "support_profile"


# --------------------------------------------------------------- fig 4

def _read(path: Path) -> list[dict]:
    return read_tsv(path) if path.exists() else []


def _paren(text: str) -> str:
    """The composition inside `nearest neighbourhood is mixed (...)`."""
    m = re.search(r"\(([^)]*)\)\s*$", text or "")
    return m.group(1) if m else ""


def _neighbourhood(comp: str) -> str:
    """A clade composition, read as a phrase.

    `ITPR2:10;ITPR3:16;unlabelled:7` is exact but unreadable at figure
    size, and it is the *variable* here — every tip on this panel has the
    same empty `home`, so annotating the home would repeat one word down
    the whole axis and show nothing.
    """
    if not comp:
        return NO_HOME
    parts = []
    for piece in comp.split(";"):
        name, _, n = piece.partition(":")
        if not name:
            continue
        parts.append(f"{n} {'unlabelled' if name == 'unlabelled' else name}")
    if not parts:
        return NO_HOME
    if len(parts) == 1 and parts[0].endswith("unlabelled"):
        return f"{parts[0]} tips, no paralog"
    return " + ".join(parts)


def fig_paralog_placement(fs, plt) -> str | None:
    rows = load("cyclostome_placement.tsv")
    audit = {a["label"]: a for a in load("membership_audit.tsv")}
    extra = [a for a in load("naming_conflicts.tsv")
             if a["label"] not in {r["label"] for r in rows}]
    items = [{"label": r["label"], "species": r["species"],
              "home": r["home"], "alrt": r["alrt"], "ufboot": r["ufboot"],
              "n": r.get("clade_size", ""),
              "what": _neighbourhood(r.get("composition", "")),
              "kind": "cyclostome locus"} for r in rows]
    # The species comes from S6's own table, not from splitting the tip
    # label on underscores: the label's shape varies with where the
    # record came from, and a parse that works on UniProt tips silently
    # mangles the genome-model ones.
    sp = {r["label"]: r.get("species", "")
          for r in _read(MSA_DIR / "representatives.tsv")}
    items += [{"label": a["label"], "species": sp.get(a["label"], ""),
               "home": a["assigned_to"],
               "alrt": a["alrt"], "ufboot": a["ufboot"],
               "n": a.get("clade_size", ""),
               "what": _neighbourhood(_paren(a.get("why", ""))),
               "kind": "census label overturned"} for a in extra]
    # A tip with no `home` is not a tip with nothing to show. On this
    # tree *every* one of them lacks a home — the loci sit in
    # cyclostome-only clades and the disputed names in mixed
    # neighbourhoods — and dropping them left the figure empty while the
    # skip line blamed a missing table. That "none of them lands in a
    # paralog" is the result, so the unplaced tips get their own
    # category and the panel draws the support each one does have.
    for i in items:
        i["home"] = i["home"] or NO_HOME
    if not items:
        return None
    items.sort(key=lambda i: (i["home"], -_f(i["ufboot"], 0)))
    y = list(range(len(items)))[::-1]

    fig, ax = plt.subplots(figsize=(fs.W_FULL,
                                    0.28 * len(items) + 1.5))
    fs.hgrid(ax, "x")
    ax.axvline(MIN_UFBOOT, color=fs.FAINT, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.annotate("UFBoot 95", (MIN_UFBOOT, len(items) - 0.4),
                xytext=(3, 0), textcoords="offset points", va="top",
                ha="left", fontsize=fs.FS_NOTE, color=fs.MUTED)
    seen: set[str] = set()
    for yi, it in zip(y, items):
        u = _f(it["ufboot"], 0.0)
        col = fs.PARALOG.get(it["home"], fs.FAINT)
        ax.plot([0, u], [yi, yi], color=fs.GRID, lw=1.3, zorder=2)
        ax.plot([u], [yi], "o", ms=5.0, color=col, mew=0, zorder=4,
                label=it["home"] if it["home"] not in seen else None)
        seen.add(it["home"])
        # direct label: identity never rests on colour alone (and the
        # aqua slot's contrast against the surface obliges a visible one)
        # not ", 4 tips" after a phrase that already counted them
        n = ("" if (not it["n"] or it["what"].endswith("tips")
                    or it["what"].endswith("paralog"))
             else f", {it['n']} tips")
        ax.annotate(f"{it['what']}{n}", (u, yi), xytext=(7, 0),
                    textcoords="offset points", va="center", ha="left",
                    fontsize=fs.FS_NOTE, color=fs.MUTED)
    names = [(i["species"] or "").split(" (")[0] or i["label"][:34]
             for i in items]
    ax.set_yticks(y, [f"{n}" for n in names], fontsize=fs.FS_TICK)
    ax.set_xlim(0, 168)
    ax.set_xticks([0, 20, 40, 60, 80, 95, 100])
    ax.set_xlabel("ultrafast bootstrap of the clade placing the tip (%)",
                  fontsize=fs.FS_LABEL, color=fs.MUTED)
    fs.despine(ax)
    if len(seen) > 1:
        ax.legend(loc="lower right", frameon=False, fontsize=fs.FS_NOTE,
                  handletextpad=0.4, borderaxespad=0.4,
                  title="resolves with", title_fontsize=fs.FS_NOTE)
    ax.set_title("Vertebrate tips the tree does not place in a paralog, "
                 "and what sits with them instead",
                 loc="left", fontsize=fs.FS_SUPTITLE,
                 fontweight="bold", color=fs.INK, pad=5)
    fig.tight_layout()
    fs.save(fig, FIGS / "paralog_placement")
    plt.close(fig)
    return "paralog_placement"


FIGURES = {"sister_au": fig_sister_au,
           "support_profile": fig_support_profile,
           "paralog_placement": fig_paralog_placement}


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import figstyle as fs
    fs.use()

    only = []
    args = sys.argv[1:]
    while args:
        a = args.pop(0)
        if a == "--only":
            only.append(args.pop(0))
        elif a == "--list":
            print("\n".join(FIGURES))
            return 0
    FIGS.mkdir(parents=True, exist_ok=True)
    made, skipped = [], []
    for slug, fn in FIGURES.items():
        if only and slug not in only:
            continue
        got = fn(fs, plt)
        (made if got else skipped).append(slug)
    for s in made:
        print(f"figure: {FIGS / s}.png (+.pdf)")
    for s in skipped:
        print(f"skipped {s}: it has nothing to draw — the table it reads "
              f"is missing or holds no row it can plot")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
