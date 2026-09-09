"""S22 stage 5 — the same two modules asked of the substitution rate.

Conservation and omega are both called constraint in prose and they are not
the same quantity: a column's dispersion is a statement about which amino
acids are seen, a site's beta is a statement about the rate at which they
change on a tree.  The two can disagree, and a task comparing modules has
to run both or it cannot tell a real difference from a property of one
metric (S17's `orthogonal` verdict, applied inside one task).

The rates come from S17's committed FEL run on S9's three codon alignments;
nothing here re-fits a model.  Three cautions carried forward from S9/S17
and applied here rather than restated:

* the aggregate `sum(beta)/sum(alpha)` is dominated by sites whose
  synonymous rate is unidentifiable, so a per-site omega is reported only
  over sites where alpha is not at the bound, and the count of excluded
  sites is a column;
* `frac_purifying` is read at FEL's own q <= 0.05, not at a p-value;
* the comparison between modules is a rank test on beta, because beta is
  zero-inflated (most sites of this protein have no non-synonymous change
  at all) and a difference of means over a zero-inflated variable is a
  statement about the zeros.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402
import s22_contacts as C  # noqa: E402


def _summarise(sites: list[dict], label: dict) -> dict:
    betas = [s["beta"] for s in sites if s["beta"] is not None]
    usable = [s for s in sites if not s["alpha_at_bound"]
              and s["omega"] is not None and s["alpha"] is not None]
    om = [s["omega"] for s in usable]
    pur = [s for s in sites if s["q_value"] is not None
           and s["q_value"] <= 0.05 and s["verdict"] == "purifying"]
    row = dict(label)
    row.update({
        "n_sites": len(sites),
        "median_beta": L.fmt(L.median(betas)),
        "mean_beta": L.fmt(L.mean(betas)),
        "n_beta_zero": sum(1 for b in betas if b == 0.0),
        "n_usable_alpha": len(usable),
        "n_alpha_at_bound": sum(1 for s in sites if s["alpha_at_bound"]),
        "median_omega_usable": L.fmt(L.median(om)),
        "frac_purifying_q05": L.fmt(len(pur) / len(sites) if sites else None),
    })
    return row


def by_module(mods: list[dict]) -> list[dict]:
    rows: list[dict] = []
    fel = L.fel_sites()
    for p in L.PARALOGS:
        sites = fel.get(p, [])
        by_resi = {s["resi"]: s for s in sites}
        regions = {
            "ligand_core": M.residues(mods, p, "contact_span"),
            "ligand_core_ibc": M.residues(mods, p, "ibc_literature"),
            "pore_module": M.residues(mods, p, "channel_minus_luminal"),
            "pore_module_all": M.residues(mods, p, "channel_all"),
        }
        for name, resis in regions.items():
            sub = [by_resi[i] for i in sorted(resis) if i in by_resi]
            rows.append(_summarise(sub, {"paralog": p, "region": name,
                                         "is_primary": name in ("ligand_core",
                                                                "pore_module")}))
        rows.append(_summarise(sites, {"paralog": p, "region": "WHOLE_PROTEIN",
                                       "is_primary": False}))
    return rows


def module_test(mods: list[dict]) -> list[dict]:
    """Ligand core against pore module, on beta, per paralogue."""
    rows: list[dict] = []
    fel = L.fel_sites()
    for core_def, pore_def, primary in (
            ("contact_span", "channel_minus_luminal", True),
            ("contact_span", "channel_all", False),
            ("ibc_literature", "channel_minus_luminal", False)):
        for p in L.PARALOGS:
            by_resi = {s["resi"]: s for s in fel.get(p, [])}
            a = [by_resi[i]["beta"] for i in sorted(M.residues(mods, p, core_def))
                 if i in by_resi and by_resi[i]["beta"] is not None]
            b = [by_resi[i]["beta"] for i in sorted(M.residues(mods, p, pore_def))
                 if i in by_resi and by_resi[i]["beta"] is not None]
            mw = L.mann_whitney(a, b)
            rows.append({
                "paralog": p, "core_definition": core_def,
                "pore_definition": pore_def, "is_primary": primary,
                "n_core_sites": mw["n_a"], "n_pore_sites": mw["n_b"],
                "median_core_beta": L.fmt(L.median(a)),
                "median_pore_beta": L.fmt(L.median(b)),
                "mean_core_beta": L.fmt(L.mean(a)),
                "mean_pore_beta": L.fmt(L.mean(b)),
                "cles_core_gt_pore": L.fmt(mw["cles"]),
                "p_mannwhitney": L.fmt(mw["p"], 6),
                "direction": ("core evolves faster" if (L.mean(a) or 0) > (L.mean(b) or 0)
                              else "pore evolves faster"),
            })
    qs = L.benjamini_hochberg([float(r["p_mannwhitney"]) if r["p_mannwhitney"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_mannwhitney"] = L.fmt(q, 6)
    return rows


def by_shell(shells: list[dict]) -> list[dict]:
    rows: list[dict] = []
    fel = L.fel_sites()
    for p in L.PARALOGS:
        by_resi = {s["resi"]: s for s in fel.get(p, [])}
        prot = [s["beta"] for s in fel.get(p, []) if s["beta"] is not None]
        for shell in C.SHELL_ORDER:
            resis = [int(s["resi"]) for s in shells
                     if s["paralog"] == p and s["shell"] == shell
                     and s["resi"] is not None]
            sub = [by_resi[i] for i in sorted(resis) if i in by_resi]
            row = _summarise(sub, {"paralog": p, "shell": shell})
            mw = L.mann_whitney([s["beta"] for s in sub if s["beta"] is not None],
                                prot, alternative="less")
            row["whole_protein_median_beta"] = L.fmt(L.median(prot))
            row["p_less_than_protein"] = L.fmt(mw["p"], 6)
            rows.append(row)
    qs = L.benjamini_hochberg([float(r["p_less_than_protein"])
                               if r["p_less_than_protein"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_less_than_protein"] = L.fmt(q, 6)
    return rows


def contact_omega(shells: list[dict], mods: list[dict]) -> list[dict]:
    """The contact residues' own rates, against both backgrounds."""
    rows: list[dict] = []
    fel = L.fel_sites()
    for p in L.PARALOGS:
        by_resi = {s["resi"]: s for s in fel.get(p, [])}
        pocket = {int(s["resi"]) for s in shells
                  if s["paralog"] == p and s["resi"] is not None}
        contacts = {int(s["resi"]) for s in shells
                    if s["paralog"] == p and s["resi"] is not None
                    and str(s["is_s0_contact"]) == "True"}
        for background, pool in (("rest_of_core",
                                  M.residues(mods, p, "contact_span")),
                                 ("rest_of_pocket", pocket)):
            vals = [by_resi[i]["beta"] for i in sorted(pool)
                    if i in by_resi and by_resi[i]["beta"] is not None]
            labs = [i in contacts for i in sorted(pool)
                    if i in by_resi and by_resi[i]["beta"] is not None]
            perm = L.permutation_label(vals, labs, iters=C.PERM_ITERS,
                                       statistic=lambda xs: -sum(xs) / len(xs))
            inside = [v for v, l in zip(vals, labs) if l]
            outside = [v for v, l in zip(vals, labs) if not l]
            rows.append({
                "paralog": p, "background": background,
                "n_background": len(vals), "n_contacts": sum(labs),
                "mean_contact_beta": L.fmt(L.mean(inside)),
                "mean_background_beta": L.fmt(L.mean(outside)),
                "n_contact_beta_zero": sum(1 for v in inside if v == 0.0),
                "p_permutation_lower": L.fmt(perm["p"], 6),
            })
    qs = L.benjamini_hochberg([float(r["p_permutation_lower"])
                               if r["p_permutation_lower"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_permutation_lower"] = L.fmt(q, 6)
    return rows


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    L.live("omega", [("by module", False), ("module test", False),
                     ("by shell", False)])
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")
    shells = L.read_tsv(L.OUT_DIR / "shell_index.tsv")
    shells = [s for s in shells if s["resi"] not in ("", "None")]

    bm = by_module(mods)
    L.write_tsv(L.OUT_DIR / "omega_by_module.tsv", bm, list(bm[0]))
    mt = module_test(mods)
    L.write_tsv(L.OUT_DIR / "omega_module_test.tsv", mt, list(mt[0]))
    bs = by_shell(shells)
    L.write_tsv(L.OUT_DIR / "omega_by_shell.tsv", bs, list(bs[0]))
    co = contact_omega(shells, mods)
    L.write_tsv(L.OUT_DIR / "omega_contacts.tsv", co, list(co[0]))
    L.live("omega", [("by module", True), ("module test", True),
                     ("by shell", True)])

    for r in mt:
        if r["is_primary"]:
            L.log(f"{r['paralog']} beta core {r['mean_core_beta']} vs pore "
                  f"{r['mean_pore_beta']}  q={r['q_mannwhitney']}  "
                  f"{r['direction']}")
    for r in co:
        L.log(f"{r['paralog']} contact beta {r['mean_contact_beta']} vs "
              f"{r['background']} {r['mean_background_beta']} "
              f"q={r['q_permutation_lower']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
