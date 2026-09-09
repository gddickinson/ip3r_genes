"""The four S18 figures, from committed S18 tables only (D13, D19).

Three choices are worth naming.

**The state palette is `figstyle.QUALITY`, not a fresh one.** Its five colours
are this project's ordered evidence scale for exactly these five states, so a
locus drawn red here is red for the same reason it is red in S5's ledger
figures. Nothing in S18 invents a colour.

**The by-source panel draws the contiguity control beside the raw contrast,
not instead of it.** D9's difference is the figure's subject and D4's bar is
its most obvious confounder — GenBank assemblies are less contiguous, and a
locus on a contig too short to hold the gene cannot be annotated completely by
anyone. Drawing only the controlled version would hide how much of the gap the
control removes; drawing only the raw one would overstate it.

**The calibration panel draws the bar rather than captioning it**, and draws
the population's own distribution behind it, because the claim being made is
that the inherited threshold sits in that distribution's far lower tail and a
number in a legend cannot show a tail.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib                                                 # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                   # noqa: E402

import figstyle as F                                              # noqa: E402
import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402

STATES = ["complete", "split", "fragmentary", "noncoding", "unannotated"]
SHORT = {"complete": "complete", "split": "split", "fragmentary": "fragmentary",
         "noncoding": "non-coding only", "unannotated": "unannotated"}


def _stacked(ax, labels, fracs, totals):
    """One horizontal stacked bar per label, states in evidence order."""
    y = range(len(labels))
    left = [0.0] * len(labels)
    for st in STATES:
        vals = [f.get(st, 0.0) for f in fracs]
        ax.barh(list(y), vals, left=left, height=0.62,
                color=F.QUALITY[st], edgecolor="white", linewidth=0.5,
                label=SHORT[st])
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(list(y))
    ax.set_yticklabels([f"{lab}\nn = {n}" for lab, n in zip(labels, totals)])
    ax.set_xlim(0, 1)
    ax.set_xlabel("share of loci")
    ax.invert_yaxis()
    F.despine(ax, keep=("bottom",))


# --------------------------------------------------------------------------
def fig1_by_source():
    src = L.read_tsv(L.OUT / "by_source.tsv")
    scopes = ["all scorable loci",
              "loci on a contig that spans the gene (D4)"]
    fig, axes = plt.subplots(2, 1, figsize=(F.W_FULL, 3.3))
    for ax, scope, letter in zip(axes, scopes, "ab"):
        rows = [r for r in src if r["scope"] == scope]
        by_state = {r["state"]: r for r in rows}
        fr = [{st: float(by_state[st]["frac_refseq"]) for st in STATES},
              {st: float(by_state[st]["frac_genbank"]) for st in STATES}]
        tot = [int(rows[0]["n_refseq_total"]), int(rows[0]["n_genbank_total"])]
        _stacked(ax, ["RefSeq (GCF_)", "GenBank (GCA_)"], fr, tot)
        title = ("every scorable locus" if letter == "a"
                 else "the same contrast above D4's contiguity bar")
        F.panel(ax, letter, title)
    axes[0].legend(loc="upper center", bbox_to_anchor=(0.5, 1.62), ncol=5,
                   frameon=False, fontsize=F.FS_TICK, handlelength=1.1,
                   columnspacing=1.0, handletextpad=0.4)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    F.save(fig, L.FIGS / "s18_fig1_by_source")
    plt.close(fig)


def fig2_calibration():
    rows = L.read_tsv(L.OUT / "locus_audit.tsv")
    for r in rows:
        r["coverage"] = float(r.get("coverage") or 0)
    pop = R.calibration_population(rows)
    vals = sorted(float(r["best_model_frac"]) for r in pop)
    cal = L.read_tsv(L.OUT / "complete_calibration.tsv")[0]
    bar, over = float(cal["bar"]), float(cal["overcredit_frac"])

    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 2.35),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    ax = axes[0]
    bins = [i / 40 for i in range(41)]
    ax.hist(vals, bins=bins, color=F.QUALITY["complete"], edgecolor="white",
            linewidth=0.3)
    ax.set_yscale("log")
    ax.axvline(bar, color=F.CLINICAL["pathogenic"], lw=1.1)
    ax.axvline(over, color=F.FAINT, lw=1.0, ls=(0, (3, 2)))
    ax.annotate(f"inherited bar {bar:g}\n({float(cal['bar_percentile']):.1%} "
                f"of this population)", xy=(bar, 0.9), xycoords=("data",
                                                                 "axes fraction"),
                xytext=(4, 0), textcoords="offset points", ha="left",
                va="top", fontsize=F.FS_NOTE, color=F.CLINICAL["pathogenic"])
    ax.annotate(f"{over:g}", xy=(over, 0.35), xycoords=("data",
                                                        "axes fraction"),
                xytext=(3, 0), textcoords="offset points", ha="left",
                fontsize=F.FS_NOTE, color=F.MUTED)
    ax.set_xlabel("best single annotated model's share of the gene")
    ax.set_ylabel("loci (log)")
    ax.set_xlim(0, 1)
    F.panel(ax, "a", f"loci the annotation names correctly (n = {len(vals)})")
    F.despine(ax)

    ax = axes[1]
    sens = L.read_tsv(L.OUT / "state_sensitivity.tsv")
    sens = [r for r in sens if float(r["min_piece"]) == R.MIN_PIECE_FRAC]
    xs = [float(r["bar"]) for r in sens]
    n = float(sens[0]["n_loci"])
    bottom = [0.0] * len(xs)
    for st in STATES:
        ys = [int(r[f"n_{st}"]) / n for r in sens]
        ax.fill_between(xs, bottom, [a + b for a, b in zip(bottom, ys)],
                        color=F.QUALITY[st], linewidth=0)
        bottom = [a + b for a, b in zip(bottom, ys)]
    ax.axvline(bar, color=F.CLINICAL["pathogenic"], lw=1.1)
    ax.set_xlabel("completeness bar")
    ax.set_ylabel("share of ITPR loci")
    ax.set_xlim(min(xs), max(xs))
    ax.set_ylim(0, 1)
    F.panel(ax, "b", "what moving the bar changes")
    F.despine(ax)
    fig.tight_layout()
    F.save(fig, L.FIGS / "s18_fig2_calibration")
    plt.close(fig)


def fig3_family_vs_control():
    cells = L.read_tsv(L.OUT / "state_by_cell.tsv")
    order = ["ITPR1", "ITPR2", "ITPR3", "RYR"]
    cells = sorted([r for r in cells if r["cell"] in order],
                   key=lambda r: order.index(r["cell"]))
    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 2.45),
                             gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axes[0]
    _stacked(ax, [r["cell"] for r in cells],
             [{st: float(r[f"frac_{st}"]) for st in STATES} for r in cells],
             [int(r["n_loci"]) for r in cells])
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.45), ncol=3,
              frameon=False, fontsize=F.FS_TICK, handlelength=1.1,
              columnspacing=0.9, handletextpad=0.4)
    F.panel(ax, "a", "the three paralogs and the sister-family control")

    ax = axes[1]
    fvc = L.read_tsv(L.OUT / "family_vs_control.tsv")
    scopes = ["all scorable loci",
              "loci on a contig that spans the gene (D4)"]
    labels = ["all scorable", "above D4's bar"]
    x = [0, 1]
    # Not `QUALITY["unannotated"]` for the ITPR bar: that red means
    # *unannotated* in panel A, and a colour cannot mean a state on one side
    # of a figure and a family on the other. The RyR bar keeps the accent
    # violet `figstyle.GROUP` already reserves for the outgroup.
    tops = []
    for i, key in enumerate(("frac_itpr_failing", "frac_ryr_failing")):
        ys = [float(next(r for r in fvc if r["scope"] == s
                         and r["measure"] == "any annotation failure")[key])
              for s in scopes]
        tops.append(ys)
        ax.bar([v + (i - 0.5) * 0.34 for v in x], ys, width=0.32,
               color=(F.BLUES[4] if i == 0 else F.ACCENT),
               label=("ITPR1-3" if i == 0 else "RyR control"))
    for j, (s, xi) in enumerate(zip(scopes, x)):
        q = float(next(r for r in fvc if r["scope"] == s
                       and r["measure"] == "any annotation failure")["q_bh"])
        ax.annotate(f"q = {q:.2f}", xy=(xi, max(tops[0][j], tops[1][j])),
                    xytext=(0, 3), textcoords="offset points", ha="center",
                    fontsize=F.FS_NOTE, color=F.MUTED)
    ax.set_ylim(0, max(max(tops[0]), max(tops[1])) * 1.28)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("share of loci with any failure")
    ax.legend(frameon=False, fontsize=F.FS_TICK, loc="upper right")
    F.hgrid(ax)
    F.despine(ax)
    F.panel(ax, "b", "is the family worse than its sister?")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    F.save(fig, L.FIGS / "s18_fig3_family_vs_control")
    plt.close(fig)


def fig4_protein_side():
    rows = L.read_tsv(L.OUT / "protein_audit.tsv")
    # Nine verdicts need a three-column key, and a key that large does not fit
    # above a 2.4 in figure: the first version's `tight_layout` gave the legend
    # the whole panel and drew both axes as hairlines. It goes at the foot of
    # the figure instead, with the rect reserving the space.
    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 2.30),
                             gridspec_kw={"width_ratios": [1.6, 1]})
    ax = axes[0]
    order = ["correct_paralog", "paralog_unspecified",
             "correct_family_wrong_paralog", "paralog_not_callable",
             "family_ambiguous", "wrong_family", "family_unnamed",
             "placeholder", "absent"]
    lab = {"correct_paralog": "names the paralog, correctly",
           "paralog_unspecified": "names the family, not the paralog",
           "correct_family_wrong_paralog": "names a different paralog",
           "paralog_not_callable": "paralog outside the panel's reach",
           "family_ambiguous": "names both families",
           "wrong_family": "names the sister family",
           "family_unnamed": "names something else",
           "placeholder": "a locus tag", "absent": "nothing"}
    colour = {"correct_paralog": F.QUALITY["complete"],
              "paralog_unspecified": F.BLUES[2],
              "correct_family_wrong_paralog": F.QUALITY["noncoding"],
              "paralog_not_callable": F.FAINT,
              "family_ambiguous": F.BLUES[0],
              "wrong_family": F.QUALITY["unannotated"],
              "family_unnamed": F.QUALITY["fragmentary"],
              "placeholder": "#e8e7e1", "absent": "#f4f3ee"}
    fields = [("gene symbol", "symbol_verdict"),
              ("protein name", "name_verdict")]
    left = [0.0, 0.0]
    n = len(rows)
    for st in order:
        vals = [sum(1 for r in rows if r[f] == st) / n for _, f in fields]
        ax.barh([0, 1], vals, left=left, height=0.5, color=colour[st],
                edgecolor="white", linewidth=0.5, label=lab[st])
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks([0, 1])
    ax.set_yticklabels([f for f, _ in fields])
    ax.set_xlim(0, 1)
    ax.set_ylim(1.6, -0.6)
    ax.set_xlabel(f"share of the {n:,} full-length family protein records")
    F.panel(ax, "a", "what the protein databases call these records")
    F.despine(ax, keep=("bottom",))

    ax = axes[1]
    z = L.read_tsv(L.OUT / "zero_hit_summary.tsv")
    names = {"gene_caller_missed_it": "the genome\nhas it",
             "assembly_cannot_carry_it": "assembly\ntoo broken",
             "genome_also_empty": "genome\nempty too",
             "undecidable_no_genome": "no genome\nin scope"}
    ys = [int(r["n_proteomes"]) for r in z]
    # One colour: the verdicts are named on the axis, so colour would carry
    # nothing here — and every hue in this figure's key already means a naming
    # verdict in panel A.
    ax.bar(range(len(z)), ys, color=F.BLUES[4], width=0.62)
    for i, v in enumerate(ys):
        ax.annotate(str(v), xy=(i, v), xytext=(0, 2),
                    textcoords="offset points", ha="center",
                    fontsize=F.FS_NOTE, color=F.INK)
    ax.set_xticks(range(len(z)))
    ax.set_xticklabels([names.get(r["verdict"], r["verdict"]) for r in z],
                       fontsize=F.FS_TICK)
    ax.set_ylabel("reference proteomes")
    ax.set_ylim(0, max(ys) * 1.20)
    ax.set_yticks(list(range(0, max(ys) + 1, 5)))
    F.hgrid(ax)
    F.despine(ax)
    F.panel(ax, "b", "the 15 proteomes S3 found nothing in")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.tight_layout(rect=(0, 0.235, 1, 1))
    fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False,
               fontsize=F.FS_TICK, handlelength=1.1, columnspacing=1.0,
               handletextpad=0.4, bbox_to_anchor=(0.5, 0.015))
    F.save(fig, L.FIGS / "s18_fig4_protein_side")
    plt.close(fig)


def main() -> None:
    F.use()
    L.FIGS.mkdir(parents=True, exist_ok=True)
    for fn in (fig1_by_source, fig2_calibration, fig3_family_vs_control,
               fig4_protein_side):
        print(f"  [s18 fig] {fn.__name__}")
        fn()


if __name__ == "__main__":
    main()
