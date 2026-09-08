"""The four S17 figures, from committed S17 tables only (D13, D19).

Four decisions worth naming, each because the obvious version of the panel says
something the data does not.

* **The per-element panel plots the composition-free metric beside the JSD.**
  JSD is a divergence from a background frequency table, so an element built of
  common amino acids scores lower at equal conservation, and the channel is
  built of exactly those. Drawing one bar per element would report a metric
  property as a biological one; the panel draws both and the reader can see
  they agree.
* **The channel profile is drawn at residue resolution, not as a mean.** The
  channel's mean conservation is the average of the most conserved stretch in
  the protein and the least, and a bar of that mean is a number no residue has.
  The panel draws the profile and marks the filter, the gate and the
  geometrically-defined luminal loop on it.
* **The classifier panel is a full ROC, not a bar of AUCs.** Four layers are
  compared and the claim is about which separates the labelled variants; an AUC
  bar hides that the layers cross.
* **The variant panel counts positions, not alleles, and prints the labelled
  n on the axis.** 88 % of this family's ClinVar missense record is uncertain,
  and a panel that showed only P/LP and B/LB would imply a labelled set an
  order of magnitude larger than the one the test ran on.
"""

from __future__ import annotations

import argparse
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import figstyle as F                                           # noqa: E402
import s17_lib as L                                            # noqa: E402

FIG_DIR = L.OUT_DIR / "figures"

#: Reading order for the elements: N to C along the protein, controls last.
ELEMENT_ORDER = ["nterm_trefoil", "MIR", "RIH_N", "RIH_C", "RIH_assoc",
                 "channel", "selectivity_filter", "gate", "luminal_loop"]
ELEMENT_LABEL = {
    "nterm_trefoil": "β-trefoil (PF08709)", "MIR": "MIR", "RIH_N": "RIH (N)",
    "RIH_C": "RIH (C)", "RIH_assoc": "RIH-assoc", "channel": "channel (TM)",
    "selectivity_filter": "filter", "gate": "gate",
    "luminal_loop": "luminal loop",
}


def _rows(name: str) -> list[dict]:
    return L.read_tsv(L.OUT_DIR / name)


def _f(x, default=None):
    return float(x) if x not in ("", None) else default


# --------------------------------------------------------------- figure 1

def fig_elements(out: Path) -> None:
    """Per-element constraint in all three paralogs, on two metrics."""
    import matplotlib.pyplot as plt

    rows = _rows("metric_controls.tsv")
    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 3.5), sharey=True)
    order = [e for e in ELEMENT_ORDER
             if any(r["element"] == e for r in rows)]
    y = list(range(len(order)))[::-1]

    for ax, key, lab in ((axes[0], "mean_jsd", "Jensen–Shannon divergence"),
                         (axes[1], "mean_frac_modal",
                          "modal-residue fraction (composition-free)")):
        for i, paralog in enumerate(L.PARALOGS):
            vals, ys = [], []
            for e, yy in zip(order, y):
                hit = [r for r in rows if r["paralog"] == paralog
                       and r["element"] == e]
                if hit:
                    vals.append(_f(hit[0][key]))
                    ys.append(yy + (i - 1) * 0.24)
            ax.scatter(vals, ys, s=16, color=F.PARALOG[paralog],
                       label=paralog if ax is axes[0] else None, zorder=3)
        ctrl = [_f(r[key]) for r in rows if r["is_control"] == "True"]
        if ctrl:
            ax.axvline(st.mean(ctrl), color=F.MUTED, lw=0.8, ls="--", zorder=1)
            # offset in points from the axes floor, not in data units: at
            # min(y) - 0.55 the label lands on the spine and is clipped by it.
            ax.annotate("linker mean", xy=(st.mean(ctrl), 0.0),
                        xycoords=("data", "axes fraction"),
                        xytext=(0, -22), textcoords="offset points",
                        fontsize=6, color=F.MUTED, ha="center", va="top",
                        annotation_clip=False)
        ax.set_xlabel(lab)
        F.despine(ax)
        F.hgrid(ax, "x")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels([ELEMENT_LABEL.get(e, e) for e in order])
    axes[0].legend(frameon=False, fontsize=6, loc="center left")
    F.panel(axes[0], "a", "constraint by element")
    F.panel(axes[1], "b", "the same, on a metric with no background model")
    fig.tight_layout()
    F.save(fig, out / "s17_elements")
    plt.close(fig)


# --------------------------------------------------------------- figure 2

def fig_channel(out: Path) -> None:
    """The pore-forming half at bin resolution, with the measured landmarks.

    **Binned, and the bins never straddle an element boundary** — that is
    `s5_figures.py`'s rule, and here it is load-bearing twice over. A sliding
    window across the luminal loop's edge would report a conservation no
    residue has, drawing the curve straight through the very boundary the panel
    exists to show; and a per-residue bar at 2,700 residues in 6.7 inches is a
    tenth of a point wide, so the "profile" would be aliasing. Bins are 10
    residues, cut at every element boundary, so a short element is one bin and
    no bin mixes two elements.
    """
    import matplotlib.pyplot as plt

    BIN = 10
    REGION = ("linker_RIH_assoc_channel", "channel", "selectivity_filter",
              "gate", "luminal_loop", "cterm")
    fig, axes = plt.subplots(3, 1, figsize=(F.W_FULL, 4.6))
    for ax, paralog in zip(axes, L.PARALOGS):
        acc = L.REFERENCES[paralog][1]
        rows = [r for r in _rows(f"constraint_{paralog}_{acc}.tsv")
                if r["element"] in REGION]
        rows.sort(key=lambda r: int(r["resi"]))
        # group into runs of one element, then bin inside each run
        runs: list[tuple[str, list[dict]]] = []
        for r in rows:
            if runs and runs[-1][0] == r["element"] \
                    and int(r["resi"]) == int(runs[-1][1][-1]["resi"]) + 1:
                runs[-1][1].append(r)
            else:
                runs.append((r["element"], [r]))
        xs, ys, cols = [], [], []
        for elem, res in runs:
            for i in range(0, len(res), BIN):
                chunk = res[i:i + BIN]
                v = [_f(c["deep_frac_modal"]) for c in chunk
                     if c["deep_frac_modal"] != ""]
                if not v:
                    continue
                xs.append((int(chunk[0]["resi"]), int(chunk[-1]["resi"])))
                ys.append(st.mean(v))
                cols.append(elem)
        base = F.PARALOG[paralog]
        accent = {"luminal_loop": "#b3261e", "selectivity_filter": "#184f95",
                  "gate": "#4a3aa7"}
        for (lo, hi), v, elem in zip(xs, ys, cols):
            ax.bar((lo + hi) / 2, v - 0.2, bottom=0.2, width=hi - lo + 1,
                   color=accent.get(elem, base), lw=0,
                   align="center", zorder=3)
        ctrl = [_f(r["mean_frac_modal"]) for r in _rows("metric_controls.tsv")
                if r["paralog"] == paralog and r["is_control"] == "True"]
        if ctrl:
            ax.axhline(st.mean(ctrl), color=F.MUTED, lw=0.7, ls="--", zorder=4)
        ax.set_ylim(0.2, 1.0)
        ax.set_xlim(min(lo for lo, _h in xs), max(h for _l, h in xs))
        ax.set_ylabel(paralog, fontsize=7, color=base)
        ax.tick_params(labelsize=6)
        F.despine(ax)
        if paralog == "ITPR3":
            from matplotlib.patches import Patch
            ax.legend(handles=[Patch(color=accent["luminal_loop"],
                                     label="luminal loop (geometric)"),
                               Patch(color=accent["selectivity_filter"],
                                     label="selectivity filter"),
                               Patch(color=accent["gate"], label="gate")],
                      frameon=False, fontsize=6, ncol=3, loc="lower center",
                      bbox_to_anchor=(0.5, -0.02))
    axes[-1].set_xlabel("residue (each paralog's own UniProt numbering); "
                        f"{BIN}-residue bins, never crossing an element boundary")
    F.panel(axes[0], "a", "modal-residue fraction across the pore-forming half; "
                          "dashed line is that protein's own linker mean")
    fig.tight_layout()
    F.save(fig, out / "s17_channel_profile")
    plt.close(fig)


# --------------------------------------------------------------- figure 3

def fig_classifier(out: Path) -> None:
    """ROC of every layer on the pooled ClinVar labelling, and the AUC table."""
    import matplotlib.pyplot as plt

    tests = _rows("variant_constraint_test.tsv")
    variants = _rows("variants.tsv")
    sites = {p: {int(r["resi"]): r for r in
                 _rows(f"constraint_{p}_{L.REFERENCES[p][1]}.tsv")}
             for p in L.PARALOGS}

    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 3.1))
    ax = axes[0]
    colours = {"deep": "#184f95", "shallow": "#86b6ef", "vert": "#eb6834",
               "family": "#1baf7a"}
    # Every curve is drawn on the **same** positions — the ones where all four
    # layers have a reliable score — because the layers do not cover the same
    # residues and four curves built on four different variant sets would be
    # comparing the sets. It is also the contrast the report's §7.2 table
    # quotes, so the figure and the table cannot show different AUCs for the
    # same comparison (D13).
    LAYERS = ("deep", "shallow", "vert", "family")

    def _scorable(gene: str, resi: int) -> bool:
        s = sites[gene].get(resi)
        if not s:
            return False
        for ly in LAYERS:
            if _f(s.get(f"{ly}_jsd", "")) is None:
                return False
            occ = s.get(f"{ly}_occupancy", "")
            if occ not in ("", None) and float(occ) < 0.50:
                return False
        return True

    for layer in LAYERS:
        pos, neg = [], []
        seen = set()
        for v in variants:
            if v["source"] != "clinvar":
                continue
            resi = int(v["resi"])
            key = (v["gene"], resi, v["class_bucket"])
            if key in seen or not _scorable(v["gene"], resi):
                continue
            seen.add(key)
            x = _f(sites[v["gene"]][resi].get(f"{layer}_jsd", ""))
            if x is None:
                continue
            if v["class_bucket"] == "P/LP":
                pos.append(x)
            elif v["class_bucket"] == "B/LB":
                neg.append(x)
        if len(pos) < 3 or len(neg) < 3:
            continue
        cuts = sorted(set(pos + neg), reverse=True)
        tpr = [0.0] + [sum(1 for p in pos if p >= c) / len(pos) for c in cuts] + [1.0]
        fpr = [0.0] + [sum(1 for n in neg if n >= c) / len(neg) for c in cuts] + [1.0]
        auc = [r for r in tests if r["gene"] == "POOLED" and r["layer"] == layer
               and r["contrast"] == "P/LP vs B/LB, all layers scorable"]
        lab = f"{layer} (AUC {auc[0]['auc']})" if auc else layer
        ax.plot(fpr, tpr, color=colours[layer], lw=1.4, label=lab)
    ax.plot([0, 1], [0, 1], color=F.MUTED, lw=0.7, ls=":")
    ax.set_xlabel("false-positive rate (benign called constrained)")
    ax.set_ylabel("true-positive rate")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    F.despine(ax)
    F.panel(ax, "a", "pathogenic vs benign, pooled; one fixed set of positions")

    ax2 = axes[1]
    strat = _rows("vus_stratification.tsv")
    genes = [g for g in L.PARALOGS if any(r["gene"] == g for r in strat)]
    width = 0.26
    for i, layer in enumerate(("deep", "vert", "family")):
        rows_l = [next((r for r in strat if r["gene"] == g
                        and r["layer"] == layer), None) for g in genes]
        vals = [_f(r["frac_vus_above_pathogenic_median"]) if r else 0
                for r in rows_l]
        xs = [j + (i - 1) * width for j in range(len(genes))]
        ax2.bar(xs, vals, width=width, color=colours[layer], label=layer)
        # each layer scores a different number of VUS — the msa_v2 layers lose
        # gappy columns the deep alignment fills — so the n goes on the bar it
        # belongs to, not on the gene's tick label.
        for x, v, r in zip(xs, vals, rows_l):
            if r:
                ax2.annotate(r["n_vus"], xy=(x, v), fontsize=5, color=F.MUTED,
                             ha="center", va="bottom", rotation=90)
    ax2.set_xticks(range(len(genes)))
    ax2.set_xticklabels(genes, fontsize=7)
    ax2.set_ylabel("VUS at or above the P/LP median")
    ax2.legend(frameon=False, fontsize=6)
    F.despine(ax2)
    F.hgrid(ax2)
    F.panel(ax2, "b", "where the uncertain variants sit")
    fig.tight_layout()
    F.save(fig, out / "s17_variant_classifier")
    plt.close(fig)


# --------------------------------------------------------------- figure 4

def fig_sites(out: Path) -> None:
    """The measured functional residues against their own element."""
    import matplotlib.pyplot as plt

    rows = _rows("functional_site_constraint.tsv")
    # 3.8 in rather than 3.0: panel b carries ten rows of three offset points
    # each, and at 3.0 in the triplets are 34 px apart in a 70 px row, which is
    # close enough to read a point onto its neighbour's label.
    fig, axes = plt.subplots(1, 2, figsize=(F.W_FULL, 3.8))

    ax = axes[0]
    classes = ["ip3_contact", "filter_lining", "gate_lining"]
    labels = {"ip3_contact": "IP3 contacts (10)", "filter_lining": "filter (2)",
              "gate_lining": "gate (2)"}
    width = 0.26
    for i, paralog in enumerate(L.PARALOGS):
        site = [next((_f(r["mean_jsd"]) for r in rows
                      if r["paralog"] == paralog and r["site_class"] == c), 0)
                for c in classes]
        ax.bar([j + (i - 1) * width for j in range(len(classes))], site,
               width=width, color=F.PARALOG[paralog], label=paralog)
    prot = [_f(r["whole_protein_mean"]) for r in rows]
    own = [_f(r["own_element_mean"]) for r in rows if r["own_element_mean"] != ""]
    ax.axhline(st.mean(prot), color=F.MUTED, lw=0.8, ls="--")
    ax.annotate("whole-protein mean", xy=(len(classes) - 0.5, st.mean(prot)),
                fontsize=6, color=F.MUTED, ha="right", va="bottom")
    ax.axhline(st.mean(own), color="#b3261e", lw=0.8, ls=":")
    ax.annotate("their own elements' mean", xy=(-0.45, st.mean(own)),
                fontsize=6, color="#b3261e", ha="left", va="bottom")
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels([labels[c] for c in classes], fontsize=6)
    ax.set_ylabel("mean JSD")
    ax.set_ylim(0.6, 0.92)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    F.despine(ax)
    F.hgrid(ax)
    F.panel(ax, "a", "the measured residues against two controls")

    ax2 = axes[1]
    ident = _rows("paralog_identity_by_element.tsv")
    order = [e for e in ELEMENT_ORDER + ["ip3_contact"]
             if any(r["element"] == e for r in ident)]
    pair_colour = {"ITPR1_vs_ITPR2": "#184f95", "ITPR1_vs_ITPR3": "#eb6834",
                   "ITPR2_vs_ITPR3": "#1baf7a"}
    y = list(range(len(order)))[::-1]
    for i, (pair, colour) in enumerate(pair_colour.items()):
        xs, ys = [], []
        for e, yy in zip(order, y):
            hit = [r for r in ident if r["pair"] == pair and r["element"] == e]
            if hit:
                xs.append(_f(hit[0]["identity"]))
                ys.append(yy + (i - 1) * 0.24)
        ax2.scatter(xs, ys, s=15, color=colour, label=pair.replace("_vs_", " × "),
                    zorder=3)
    whole = [_f(r["whole_protein_identity"]) for r in ident
             if r["whole_protein_identity"] != ""]
    if whole:
        ax2.axvline(st.mean(whole), color=F.MUTED, lw=0.8, ls="--")
    ax2.set_yticks(y)
    ax2.set_yticklabels([ELEMENT_LABEL.get(e, "IP3 contacts") for e in order],
                        fontsize=6)
    ax2.set_xlabel("identity between paralogs (covered columns)")
    ax2.legend(frameon=False, fontsize=6, loc="center left")
    F.despine(ax2)
    F.hgrid(ax2, "x")
    F.panel(ax2, "b", "what the three copies have kept alike")
    fig.tight_layout()
    F.save(fig, out / "s17_functional_sites")
    plt.close(fig)


FIGURES = [("elements", fig_elements), ("channel", fig_channel),
           ("classifier", fig_classifier), ("sites", fig_sites)]


def run(out_dir: Path = L.OUT_DIR) -> dict:
    F.use()
    out = out_dir / "figures"
    out.mkdir(parents=True, exist_ok=True)
    made = []
    for name, fn in FIGURES:
        try:
            fn(out)
            made.append(name)
            print(f"[s17] figure {name}")
        except FileNotFoundError as e:
            print(f"[s17] figure {name}: skipped — {e}")
    return {"figures": made}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    run(ap.parse_args().out)


if __name__ == "__main__":
    main()
