"""The four S22 figures, from committed S22 tables only (D13, D19).

Three drawing decisions are the figures' own arguments and are made here
rather than described in a caption.

**Figure 1 draws the boundary the answer turns on.**  The core-versus-pore
result reverses depending on whether a fifty-residue luminal loop is left
inside the pore module, so the panel plots the paired per-tip difference
under *both* definitions on one axis with zero marked.  A figure showing
only the primary definition would be asserting the choice instead of
showing what it costs.

**Figure 2 is a scatter against distance and not a bar of shells.**  The
claim is that constraint is elevated across the whole pocket and does not
fall off sharply at the contact residues; four bars cannot show the absence
of a step, and a scatter with the shell edges drawn can.

**Figure 4 draws the power curve beside the result.**  A lineage test that
finds nothing is only readable next to what it could have found, so the
minimum detectable shift and the shift the RyR control actually produces
are on the same axis as the observed strata.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import figstyle as FS  # noqa: E402
import s22_lib as L  # noqa: E402
import s22_shells as SH  # noqa: E402

FS.use()
# The ligand core is a *region*, not a group, so it takes a neutral dark
# rather than a hue: `figstyle.GROUP` reserves violet for the ryanodine
# receptors and figure 4 draws them, so a violet region here would collide
# with the one control the task depends on.
CORE_INK = "#3d3b34"
INK = FS.INK if hasattr(FS, "INK") else "#22211d"


def _f(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------

def fig1_modules() -> None:
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")
    contrast = L.read_tsv(L.OUT_DIR / "module_contrast.tsv")
    tips = L.read_tsv(L.OUT_DIR / "tip_divergence.tsv")
    dm = L.domain_map()

    fig, axes = plt.subplots(2, 1, figsize=(FS.W_FULL, 4.4),
                             gridspec_kw={"height_ratios": [1.0, 1.35]})
    ax = axes[0]
    order = ["ITPR1", "ITPR2", "ITPR3"]
    for i, p in enumerate(order):
        y = len(order) - 1 - i
        length = max(r["resi"] for r in L.constraint(p))
        ax.plot([1, length], [y, y], color="#c9c7bf", lw=3,
                solid_capstyle="butt", zorder=1)
        for r in dm:
            if r["paralog"] != p or r["kind"] != "pfam":
                continue
            ax.plot([r["start"], r["end"]], [y, y], color="#9c9a92", lw=7,
                    solid_capstyle="butt", zorder=2)
        for m in mods:
            if m["paralog"] != p or m["is_primary"] != "True":
                continue
            col = CORE_INK if m["module"] == "ligand_core" else FS.PARALOG[p]
            ax.plot([int(m["start"]), int(m["end"])], [y, y], color=col, lw=9,
                    solid_capstyle="butt", zorder=3)
        lum = next(r for r in dm if r["paralog"] == p
                   and r["element"] == "luminal_loop")
        ax.plot([lum["start"], lum["end"]], [y, y], color="#d9534f", lw=9,
                solid_capstyle="butt", zorder=4)
        ax.text(-40, y, p, ha="right", va="center", fontsize=7)
    ax.set_ylim(-0.9, len(order) - 0.05)
    ax.set_xlim(-420, 2900)
    ax.set_yticks([])
    ax.set_xlabel("residue")
    FS.despine(ax, keep=("bottom",))
    FS.panel(ax, "a", "the two modules on the receptor")
    # The pore module is drawn in each paralogue's own colour, so it is
    # annotated on the track rather than given a legend swatch that could
    # only show one of the three.
    top = len(order) - 1
    core = next(m for m in mods if m["paralog"] == order[0]
                and m["module"] == "ligand_core" and m["is_primary"] == "True")
    pore = next(m for m in mods if m["paralog"] == order[0]
                and m["module"] == "pore_module" and m["is_primary"] == "True")
    for m, text in ((core, "ligand core"), (pore, "pore module")):
        mid = (int(m["start"]) + int(m["end"])) / 2
        ax.annotate(text, xy=(mid, top), xytext=(0, 11),
                    textcoords="offset points", ha="center", fontsize=6.2,
                    color="#52514e")
    from matplotlib.lines import Line2D
    ax.legend(handles=[
        Line2D([], [], color=CORE_INK, lw=5, label="ligand core"),
        Line2D([], [], color="#d9534f", lw=5, label="luminal loop"),
        Line2D([], [], color="#9c9a92", lw=4, label="Pfam domain")],
        loc="lower right", frameon=False, fontsize=6, ncol=3,
        handlelength=1.4, columnspacing=1.0, borderaxespad=0.2)

    ax = axes[1]
    xs, labels, colours = [], [], []
    for pore_def, tag in (("channel_minus_luminal", "loop out\n(primary)"),
                          ("channel_all", "loop in")):
        for p in order:
            vals = []
            for t in tips:
                if t["paralog"] != p:
                    continue
                a, b = _f(t["contact_span_identity"]), _f(f"{t[pore_def+'_identity']}")
                ca, cb = _f(t["contact_span_coverage"], 0), _f(t[pore_def + "_coverage"], 0)
                if a is None or b is None or ca < 0.5 or cb < 0.5:
                    continue
                vals.append(a - b)
            xs.append(vals)
            labels.append(f"{p}\n{tag}")
            colours.append(FS.PARALOG[p])
    bp = ax.boxplot(xs, patch_artist=True, widths=0.6, showfliers=False,
                    medianprops={"color": "#22211d", "lw": 1.1})
    for patch, col, i in zip(bp["boxes"], colours, range(len(xs))):
        patch.set_facecolor(col)
        patch.set_alpha(0.85 if i < 3 else 0.35)
        patch.set_edgecolor("#6f6d66")
    ax.axhline(0, color="#22211d", lw=0.8, ls=(0, (3, 2)))
    ax.axvline(3.5, color="#c9c7bf", lw=0.8)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("core − pore identity, per orthologue")
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "b", "the paired difference, under both pore definitions")
    fig.tight_layout()
    FS.save(fig, L.FIG_DIR / "s22_fig1_modules")
    plt.close(fig)


def fig2_shells() -> None:
    shells = L.read_tsv(L.OUT_DIR / "shell_index.tsv")
    sc = L.read_tsv(L.OUT_DIR / "shell_constraint.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.7))
    ax = axes[0]
    for p in ("ITPR1", "ITPR2", "ITPR3"):
        by_resi = {r["resi"]: r for r in L.constraint(p)}
        xs, ys, mark = [], [], []
        for s in shells:
            if s["paralog"] != p or s["resi"] in ("", "None"):
                continue
            r = by_resi.get(int(s["resi"]))
            if r is None or r["deep_jsd"] is None:
                continue
            xs.append(_f(s["median_distance_A"]))
            ys.append(r["deep_jsd"])
            mark.append(s["is_s0_contact"] == "True")
        ax.scatter([x for x, m in zip(xs, mark) if not m],
                   [y for y, m in zip(ys, mark) if not m],
                   s=6, color=FS.PARALOG[p], alpha=0.55, lw=0)
        ax.scatter([x for x, m in zip(xs, mark) if m],
                   [y for y, m in zip(ys, mark) if m],
                   s=26, facecolor="none", edgecolor=FS.PARALOG[p], lw=1.0)
    prot = L.mean([r["deep_jsd"] for r in L.constraint("ITPR1")
                   if r["deep_jsd"] is not None])
    ax.axhline(prot, color="#6f6d66", lw=0.9, ls=(0, (4, 2)))
    ax.annotate("whole-protein mean", xy=(14.6, prot), xytext=(0, -8),
                textcoords="offset points", ha="right", fontsize=6,
                color="#6f6d66")
    for _lo, hi, _name in SH.SHELL_EDGES[:-1]:
        ax.axvline(hi, color="#e0ded6", lw=0.7, zorder=0)
    ax.set_xlabel("distance to IP3, Å (median over 6 structures)")
    ax.set_ylabel("conservation (JSD, deep layer)")
    ax.set_xlim(0, 15.2)
    FS.despine(ax)
    FS.panel(ax, "a", "the pocket, residue by residue")
    ax.legend(handles=FS.paralog_handles(), loc="lower left", frameon=False,
              fontsize=6, ncol=3, handlelength=1.0, columnspacing=0.8)

    ax = axes[1]
    order = ["contact", "second", "third", "fourth"]
    width = 0.26
    for i, p in enumerate(("ITPR1", "ITPR2", "ITPR3")):
        vals = [_f(next((r["mean_jsd"] for r in sc
                         if r["paralog"] == p and r["shell"] == s), None))
                for s in order]
        ax.bar([x + (i - 1) * width for x in range(len(order))], vals,
               width=width, color=FS.PARALOG[p], edgecolor="none")
    ax.axhline(prot, color="#6f6d66", lw=0.9, ls=(0, (4, 2)))
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["contact\n≤4.5 Å", "second\n4.5–8", "third\n8–11.5",
                        "fourth\n11.5–15"], fontsize=6)
    ax.set_ylabel("mean conservation")
    ax.set_ylim(0.70, 0.86)
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "b", "no step at the contact shell")
    fig.tight_layout()
    FS.save(fig, L.FIG_DIR / "s22_fig2_shells")
    plt.close(fig)


def fig3_omega() -> None:
    bm = L.read_tsv(L.OUT_DIR / "omega_by_module.tsv")
    bs = L.read_tsv(L.OUT_DIR / "omega_by_shell.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(FS.W_FULL, 2.6))
    ax = axes[0]
    regions = ["ligand_core", "pore_module", "WHOLE_PROTEIN"]
    width = 0.26
    for i, p in enumerate(("ITPR1", "ITPR2", "ITPR3")):
        vals = [_f(next((r["mean_beta"] for r in bm
                         if r["paralog"] == p and r["region"] == g), None), 0)
                for g in regions]
        ax.bar([x + (i - 1) * width for x in range(len(regions))], vals,
               width=width, color=FS.PARALOG[p], edgecolor="none")
    ax.set_xticks(range(len(regions)))
    ax.set_xticklabels(["ligand core", "pore module", "whole protein"],
                       fontsize=6.5)
    ax.set_ylabel("mean β (non-synonymous rate)")
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "a", "substitution rate by module")

    ax = axes[1]
    order = ["contact", "second", "third", "fourth"]
    for i, p in enumerate(("ITPR1", "ITPR2", "ITPR3")):
        vals = [_f(next((r["frac_purifying_q05"] for r in bs
                         if r["paralog"] == p and r["shell"] == s), None), 0)
                for s in order]
        # A deterministic offset by paralogue index, not a jitter: all
        # three sit exactly on 1.0 at the contact shell and one marker
        # would otherwise hide the other two (D24 applied to a drawing).
        ax.plot([x + (i - 1) * 0.045 for x in range(len(order))], vals,
                marker="o", ms=3.5, lw=1.2, color=FS.PARALOG[p])
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(["contact", "second", "third", "fourth"], fontsize=6.5)
    ax.set_ylabel("share of sites purifying (FEL q ≤ 0.05)")
    ax.set_ylim(0, 1.02)
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "b", "and by shell")
    ax.legend(handles=FS.paralog_handles(), loc="lower left", frameon=False,
              fontsize=6, ncol=3, handlelength=1.0, columnspacing=0.8)
    fig.tight_layout()
    FS.save(fig, L.FIG_DIR / "s22_fig3_omega")
    plt.close(fig)


def fig4_lineage() -> None:
    co = L.read_tsv(L.OUT_DIR / "pathway_cooccurrence.tsv")
    panel = L.read_tsv(L.OUT_DIR / "deep_lineage_panel.tsv")
    ryr = L.read_tsv(L.OUT_DIR / "deep_lineage_ryr_control.tsv")
    matched = L.read_tsv(L.OUT_DIR / "deep_lineage_matched.tsv")
    mp = L.read_tsv(L.OUT_DIR / "deep_lineage_matched_power.tsv")
    fig, axes = plt.subplots(1, 3, figsize=(FS.W_FULL, 2.7),
                             gridspec_kw={"width_ratios": [0.95, 1.25, 1.0]})

    ax = axes[0]
    # The bar is the *fraction* of ITPR-carrying proteomes with no PI-PLC,
    # with the counts annotated.  Overlaying two counts on a log axis — the
    # obvious way to fit 745 vertebrate and 15 plant proteomes on one scale —
    # misreads proportions by construction (s19_figures' rule), and the
    # proportion is the result.
    rows = sorted((g for g in co if int(g["n_itpr"]) > 0),
                  key=lambda g: int(g["n_itpr_no_plc"]) / int(g["n_itpr"]))
    groups = [g["group"] for g in rows]
    frac = [int(g["n_itpr_no_plc"]) / int(g["n_itpr"]) for g in rows]
    ax.barh(range(len(groups)), frac, color="#b3261e", height=0.6)
    for i, g in enumerate(rows):
        ax.annotate(f"{g['n_itpr_no_plc']}/{g['n_itpr']}", xy=(frac[i], i),
                    xytext=(4, 0), textcoords="offset points", va="center",
                    fontsize=6, color="#52514e")
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels([g.replace("_", " ") for g in groups], fontsize=6)
    ax.set_xlabel("share of ITPR-carrying proteomes\nwith no PI-PLC")
    ax.set_xlim(0, 1.0)
    FS.despine(ax)
    FS.panel(ax, "a", "receptor without enzyme")

    # The confound and the control, on one axis: the paired difference
    # against how far the whole protein has moved from the reference.  The
    # two PLC cells occupy different halves of the x-axis, which is the
    # result §9.3 reports, and the RyR control sits at the test group's own
    # divergence rather than at the reference group's.
    ax = axes[1]
    live = [r for r in panel if r["in_test"] == "True"]
    for status, colour, size in (("present", "#9c9a92", 7),
                                 ("absent", "#b3261e", 16)):
        xs = [_f(r["pore_identity"]) for r in live if r["plc_status"] == status]
        ys = [_f(r["delta_core_minus_pore"]) for r in live
              if r["plc_status"] == status]
        ax.scatter(xs, ys, s=size, color=colour, alpha=0.7, lw=0,
                   label=f"PI-PLC {status} (n={len(xs)})")
    rx = [_f(r["pore_identity"]) for r in ryr if r["tip"] != "SUMMARY"]
    ry = [_f(r["delta_core_minus_pore"]) for r in ryr if r["tip"] != "SUMMARY"]
    ax.scatter(rx, ry, s=34, marker="D", facecolor="none",
               edgecolor=FS.GROUP["RYR"], lw=1.1,
               label=f"RyR control (n={len(rx)})")
    ax.axhline(0, color="#22211d", lw=0.7, ls=(0, (3, 2)))
    ax.set_xlabel("pore identity to human ITPR3 (divergence)")
    ax.set_ylabel("core − pore identity")
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "b", "the confound, and the control")
    ax.legend(loc="lower right", frameon=False, fontsize=5.4)

    ax = axes[2]
    xs = [_f(r["shift"]) for r in mp if r["shift"] != ""]
    ys = [_f(r["power"]) for r in mp if r["shift"] != ""]
    ax.plot(xs, ys, marker="o", ms=3.5, lw=1.2, color="#2a78d6")
    ax.axhline(0.8, color="#6f6d66", lw=0.9, ls=(0, (4, 2)))
    summ = next((r for r in ryr if r["tip"] == "SUMMARY"), {})
    cs = _f(summ.get("shift_vs_itpr", ""))
    if cs is not None:
        ax.axvline(abs(cs), color=FS.GROUP["RYR"], lw=1.1)
        ax.annotate("RyR control", xy=(abs(cs), 0.07), xytext=(-3, 0),
                    textcoords="offset points", fontsize=5.8, ha="right",
                    color=FS.GROUP["RYR"])
    obs = [abs(_f(r["difference"])) for r in matched if r["difference"] != ""]
    med = L.median([_f(r["difference"]) for r in matched
                    if r["difference"] != ""])
    if med is not None:
        ax.axvline(abs(med), color="#b3261e", lw=1.1)
        ax.annotate("observed", xy=(abs(med), 0.9), xytext=(3, 0),
                    textcoords="offset points", fontsize=5.8,
                    color="#b3261e")
    ax.set_xlabel("true shift in core − pore identity")
    ax.set_ylabel("power of the matched test, α = 0.05")
    ax.set_ylim(0, 1.03)
    FS.despine(ax)
    FS.hgrid(ax)
    FS.panel(ax, "c", "what the test could have seen")
    fig.tight_layout()
    FS.save(fig, L.FIG_DIR / "s22_fig4_lineage")
    plt.close(fig)


FIGURES = {"fig1": fig1_modules, "fig2": fig2_shells, "fig3": fig3_omega,
           "fig4": fig4_lineage}


def main(argv: list[str] | None = None) -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", action="append", choices=sorted(FIGURES))
    a = ap.parse_args(argv)
    L.FIG_DIR.mkdir(parents=True, exist_ok=True)
    for name in (a.only or sorted(FIGURES)):
        L.log(f"drawing {name}")
        FIGURES[name]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
