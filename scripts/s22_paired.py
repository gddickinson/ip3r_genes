"""S22 stage 3 — the ligand core against the pore, paired inside one alignment.

The comparison the brief asks for is easy to get wrong in one specific way:
measure the core on one set of sequences and the pore on another and the
difference reported is a difference in alignment depth, ortholog quality or
taxon sampling.  Both modules sit in the same protein, so the fix is to
never let the sets differ.

Three pairings, in increasing order of what they control:

* `column_contrast` — every core column against every pore column, in one
  paralogue's own deep alignment.  The sequences are identical by
  construction; what is not controlled is that the two modules have
  different amino-acid composition and different gap structure, so the
  occupancy of each column is reported beside its score.
* `tip_divergence` / `module_contrast` — **one core number and one pore
  number per ortholog**, so the pairing is at the level of the sequence:
  a tip that is generally divergent contributes a divergent core *and* a
  divergent pore, and only the difference between them enters the test.
  This is the design with the power — 249-265 pairs per paralogue.
* `paralog_contrast` — the same module in two paralogues, paired column by
  column through the msa_v2 alignment the three share, so a difference
  cannot be a difference in where the module's boundaries fall.

Every tip must cover both modules to enter the paired test, and the tips
that do not are written out with which module lost them.  A tip covering
the pore and not the core would otherwise enter the test as an extreme
core divergence, which is the exact artefact a gene model truncated at the
N-terminus produces.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402

MIN_MODULE_COVERAGE = 0.50   # a tip must resolve half of each module
DEFINITIONS = (("ligand_core", "contact_span"),
               ("ligand_core", "ibc_literature"),
               ("pore_module", "channel_minus_luminal"),
               ("pore_module", "channel_all"))


def _module_columns(mods: list[dict], paralog: str,
                    definition: str) -> list[int]:
    """The deep-alignment columns a module occupies in one paralogue."""
    want = M.residues(mods, paralog, definition)
    return [r["deep_col"] for r in L.constraint(paralog)
            if r["resi"] in want and r["deep_col"] is not None]


def tip_divergence(mods: list[dict]) -> list[dict]:
    """Per ortholog tip: covered identity to the reference in each module."""
    rows: list[dict] = []
    for p in L.PARALOGS:
        aln = L.deep_alignment(p)
        ref_label = L.deep_reference_label(p)
        ref = aln[ref_label]
        cols = {d: _module_columns(mods, p, d) for _, d in DEFINITIONS}
        for label, seq in aln.items():
            if label == ref_label:
                continue
            rec = {"paralog": p, "tip": label}
            for _mod, d in DEFINITIONS:
                cc = cols[d]
                n_ref = sum(1 for c in cc if ref[c] != "-")
                covered = [c for c in cc if ref[c] != "-" and seq[c] != "-"]
                ident = sum(1 for c in covered if seq[c] == ref[c])
                rec[f"{d}_n_ref"] = n_ref
                rec[f"{d}_n_covered"] = len(covered)
                rec[f"{d}_coverage"] = len(covered) / n_ref if n_ref else 0.0
                rec[f"{d}_identity"] = ident / len(covered) if covered else None
            rows.append(rec)
    return rows


def _pairs(tips: list[dict], paralog: str, def_a: str, def_b: str):
    """(tip, identity_a, identity_b) for tips resolving both modules."""
    out, dropped = [], []
    for t in tips:
        if t["paralog"] != paralog:
            continue
        ca, cb = t[f"{def_a}_coverage"], t[f"{def_b}_coverage"]
        ia, ib = t[f"{def_a}_identity"], t[f"{def_b}_identity"]
        if ia is None or ib is None or ca < MIN_MODULE_COVERAGE or cb < MIN_MODULE_COVERAGE:
            lost = def_a if (ia is None or ca < MIN_MODULE_COVERAGE) else def_b
            dropped.append((t["tip"], lost, ca, cb))
            continue
        out.append((t["tip"], ia, ib))
    return out, dropped


def module_contrast(tips: list[dict]) -> tuple[list[dict], list[dict]]:
    """The paired test, run for the primary pair and both sensitivity pairs."""
    combos = [("contact_span", "channel_minus_luminal", True),
              ("ibc_literature", "channel_minus_luminal", False),
              ("contact_span", "channel_all", False),
              ("ibc_literature", "channel_all", False)]
    rows: list[dict] = []
    drops: list[dict] = []
    for core_def, pore_def, primary in combos:
        for p in L.PARALOGS:
            pairs, dropped = _pairs(tips, p, core_def, pore_def)
            for tip, lost, ca, cb in dropped:
                drops.append({"paralog": p, "core_definition": core_def,
                              "pore_definition": pore_def, "tip": tip,
                              "lost_on": lost, "core_coverage": L.fmt(ca),
                              "pore_coverage": L.fmt(cb)})
            diffs = [a - b for _, a, b in pairs]
            w = L.wilcoxon(diffs)
            s = L.sign_test(diffs)
            lo, hi = L.bootstrap_ci(diffs)
            rows.append({
                "paralog": p, "core_definition": core_def,
                "pore_definition": pore_def, "is_primary": primary,
                "n_tips": len(pairs), "n_dropped": len(dropped),
                "mean_core_identity": L.fmt(L.mean([a for _, a, _ in pairs])),
                "mean_pore_identity": L.fmt(L.mean([b for _, _, b in pairs])),
                "mean_difference": L.fmt(L.mean(diffs)),
                "median_difference": L.fmt(L.median(diffs)),
                "ci95_lo": L.fmt(lo), "ci95_hi": L.fmt(hi),
                "n_core_more_conserved": s["n_pos"],
                "n_pore_more_conserved": s["n_neg"], "n_ties": s["n_ties"],
                "p_sign": L.fmt(s["p"], 6), "p_wilcoxon": L.fmt(w["p"], 6),
                "direction": ("core > pore" if (L.mean(diffs) or 0) > 0
                              else "pore > core"),
            })
    qs = L.benjamini_hochberg([float(r["p_wilcoxon"]) if r["p_wilcoxon"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_wilcoxon"] = L.fmt(q, 6)
    return rows, drops


def column_contrast(mods: list[dict]) -> list[dict]:
    """Per-column constraint, core against pore, in each conservation layer."""
    rows: list[dict] = []
    for p in L.PARALOGS:
        by_resi = {r["resi"]: r for r in L.constraint(p)}
        for core_def, pore_def, primary in (
                ("contact_span", "channel_minus_luminal", True),
                ("ibc_literature", "channel_all", False)):
            core = M.residues(mods, p, core_def)
            pore = M.residues(mods, p, pore_def)
            for layer in L.LAYERS:
                a = [by_resi[i][f"{layer}_jsd"] for i in sorted(core)
                     if i in by_resi]
                b = [by_resi[i][f"{layer}_jsd"] for i in sorted(pore)
                     if i in by_resi]
                mw = L.mann_whitney(a, b)
                oa = [by_resi[i][f"{layer}_occupancy"] for i in sorted(core)
                      if i in by_resi]
                ob = [by_resi[i][f"{layer}_occupancy"] for i in sorted(pore)
                      if i in by_resi]
                # The composition control.  JSD is a divergence from a
                # background amino-acid table, so a transmembrane module
                # scores lower than a soluble one at equal conservation
                # (S17 measured this).  `frac_modal` — the weighted share
                # of sequences carrying the column's commonest residue —
                # has no background term, so the two metrics disagreeing
                # is itself the reading.
                fa = fb = ea = eb = None
                if layer == "deep":
                    fa = L.mean([float(by_resi[i]["deep_frac_modal"])
                                 for i in sorted(core) if i in by_resi])
                    fb = L.mean([float(by_resi[i]["deep_frac_modal"])
                                 for i in sorted(pore) if i in by_resi])
                    ea = L.mean([float(by_resi[i]["deep_entropy"])
                                 for i in sorted(core) if i in by_resi])
                    eb = L.mean([float(by_resi[i]["deep_entropy"])
                                 for i in sorted(pore) if i in by_resi])
                rows.append({
                    "paralog": p, "layer": layer, "is_primary": primary,
                    "core_definition": core_def, "pore_definition": pore_def,
                    "n_core": mw["n_a"], "n_pore": mw["n_b"],
                    "mean_core_jsd": L.fmt(L.mean(a)),
                    "mean_pore_jsd": L.fmt(L.mean(b)),
                    "median_core_jsd": L.fmt(L.median(a)),
                    "median_pore_jsd": L.fmt(L.median(b)),
                    "difference": L.fmt((L.mean(a) or 0) - (L.mean(b) or 0)),
                    "cles_core_gt_pore": L.fmt(mw["cles"]),
                    "p_mannwhitney": L.fmt(mw["p"], 6),
                    "mean_core_occupancy": L.fmt(L.mean(oa)),
                    "mean_pore_occupancy": L.fmt(L.mean(ob)),
                    "mean_core_frac_modal": L.fmt(fa),
                    "mean_pore_frac_modal": L.fmt(fb),
                    "mean_core_entropy": L.fmt(ea),
                    "mean_pore_entropy": L.fmt(eb),
                })
    qs = L.benjamini_hochberg([float(r["p_mannwhitney"]) if r["p_mannwhitney"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_mannwhitney"] = L.fmt(q, 6)
    return rows


def paralog_contrast(mods: list[dict]) -> list[dict]:
    """The same module in two paralogues, paired column by column.

    The frame is msa_v2's alignment, which is what makes this a pairing:
    without it the two paralogues' modules have different lengths and the
    comparison would be of two unrelated distributions.
    """
    frames = {}
    for p in L.PARALOGS:
        frames[p] = {int(r["msa_col"]): r for r in L.constraint(p)
                     if r.get("msa_col") not in ("", None)}
    rows: list[dict] = []
    pairs = [("ITPR1", "ITPR2"), ("ITPR1", "ITPR3"), ("ITPR2", "ITPR3")]
    for module, definition in (("ligand_core", "contact_span"),
                               ("pore_module", "channel_minus_luminal")):
        for a, b in pairs:
            ra = M.residues(mods, a, definition)
            rb = M.residues(mods, b, definition)
            shared = [(c, frames[a][c], frames[b][c])
                      for c in sorted(set(frames[a]) & set(frames[b]))
                      if frames[a][c]["resi"] in ra and frames[b][c]["resi"] in rb]
            diffs = [x["deep_jsd"] - y["deep_jsd"] for _, x, y in shared
                     if x["deep_jsd"] is not None and y["deep_jsd"] is not None]
            s = L.sign_test(diffs)
            w = L.wilcoxon(diffs)
            ident = sum(1 for _, x, y in shared if x["aa"] == y["aa"])
            rows.append({
                "module": module, "definition": definition,
                "paralog_a": a, "paralog_b": b,
                "n_shared_columns": len(shared),
                "n_identical_residues": ident,
                "frac_identical": L.fmt(ident / len(shared) if shared else None),
                "mean_jsd_a": L.fmt(L.mean([x["deep_jsd"] for _, x, _ in shared])),
                "mean_jsd_b": L.fmt(L.mean([y["deep_jsd"] for _, _, y in shared])),
                "mean_difference": L.fmt(L.mean(diffs)),
                "n_a_higher": s["n_pos"], "n_b_higher": s["n_neg"],
                "n_ties": s["n_ties"],
                "p_sign": L.fmt(s["p"], 6), "p_wilcoxon": L.fmt(w["p"], 6),
            })
    qs = L.benjamini_hochberg([float(r["p_wilcoxon"]) if r["p_wilcoxon"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_wilcoxon"] = L.fmt(q, 6)
    return rows


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = [("tip divergence", False), ("module contrast", False),
             ("column contrast", False), ("paralogue contrast", False)]
    L.live("paired", steps)
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")
    for r in mods:
        r["start"] = int(r["start"]); r["end"] = int(r["end"])

    tips = tip_divergence(mods)
    cols = ["paralog", "tip"]
    for _m, d in DEFINITIONS:
        cols += [f"{d}_n_ref", f"{d}_n_covered", f"{d}_coverage", f"{d}_identity"]
    L.write_tsv(L.OUT_DIR / "tip_divergence.tsv", tips, cols)
    steps[0] = ("tip divergence", True); L.live("paired", steps)

    contrast, drops = module_contrast(tips)
    L.write_tsv(L.OUT_DIR / "module_contrast.tsv", contrast, list(contrast[0]))
    L.write_tsv(L.OUT_DIR / "module_contrast_dropped.tsv", drops,
                ["paralog", "core_definition", "pore_definition", "tip",
                 "lost_on", "core_coverage", "pore_coverage"])
    steps[1] = ("module contrast", True); L.live("paired", steps)

    colrows = column_contrast(mods)
    L.write_tsv(L.OUT_DIR / "column_contrast.tsv", colrows, list(colrows[0]))
    steps[2] = ("column contrast", True); L.live("paired", steps)

    par = paralog_contrast(mods)
    L.write_tsv(L.OUT_DIR / "paralog_contrast.tsv", par, list(par[0]))
    steps[3] = ("paralogue contrast", True); L.live("paired", steps)

    for r in contrast:
        if r["is_primary"]:
            L.log(f"{r['paralog']}: core {r['mean_core_identity']} vs pore "
                  f"{r['mean_pore_identity']}  diff {r['mean_difference']}  "
                  f"{r['n_core_more_conserved']}+/{r['n_pore_more_conserved']}- "
                  f"of {r['n_tips']}  q={r['q_wilcoxon']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
