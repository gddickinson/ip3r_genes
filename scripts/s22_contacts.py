"""S22 stage 4 — the contacts against the rest of the core, and the gradient.

Two questions, and the second is the one with the power.

**Are the ten residues that touch IP3 more constrained than the rest of the
binding core?**  S17 asked a version of this against MIR and RIH_N, the
Pfam domains the contacts happen to fall in.  Asked against the *core* the
test is sharper and the null is stated: permute which of the core's ~300
positions carry the contact label, keeping the constraint values where they
are.  The unit of resampling is the label, because the hypothesis is about
these particular residues and not about this particular set of scores.

**Does constraint decay with distance from the ligand?**  A binary
contact/not label throws away the only thing the structure can say that an
alignment cannot.  The shells measured in stage 2 turn it into a trend, and
a trend over four bins and 125 residues is a much harder thing to produce
by accident than a difference between 10 residues and 295.

The shells are measured on human ITPR3 and carried to the other two
paralogues through S17's `transfer_positions`, the same routine every other
cross-paralogue coordinate in this project went through.  A shell residue
that does not transfer is dropped and counted; it is not assigned to the
nearest position that did.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402
import s17_lib as S17  # noqa: E402

SHELL_ORDER = ("contact", "second", "third", "fourth")
PERM_ITERS = 200000
# A residue is called a consensus contact when a majority of the six IP3-bound
# depositions place it within 4.5 A.  S0's ten all clear 6/6; the bar exists
# so a residue resolved as a contact in one map alone cannot join them.
CONSENSUS_MAJORITY = 4


_TRANSFER: dict[str, dict[int, int]] = {}


def transfer(paralog: str) -> dict[int, int]:
    """ITPR3 (structure) residue -> this paralogue's residue."""
    if paralog == L.STRUCTURE_REF:
        return {}
    if paralog not in _TRANSFER:
        src = S17.uniprot_fasta(L.STRUCTURE_ACC)
        dst = S17.uniprot_fasta(L.REFERENCES[paralog][1])
        _TRANSFER[paralog] = S17.transfer_positions(src, dst)
    return _TRANSFER[paralog]


def shell_index() -> list[dict]:
    """Every measured shell residue, in all three paralogues' numbering."""
    shells = L.read_tsv(L.OUT_DIR / "ligand_shells.tsv")
    out: list[dict] = []
    for s in shells:
        resi3 = int(s["resi"])
        n_contact = int(s["n_structures_contact"])
        for p in L.PARALOGS:
            if p == L.STRUCTURE_REF:
                resi = resi3
                moved = "measured_here"
            else:
                m = transfer(p)
                if resi3 not in m:
                    out.append({"paralog": p, "structure_resi": resi3,
                                "resi": None, "shell": s["shell"],
                                "median_distance_A": s["median_distance_A"],
                                "n_structures_contact": n_contact,
                                "is_s0_contact": s["is_s0_contact"],
                                "is_consensus_contact": n_contact >= CONSENSUS_MAJORITY,
                                "transfer": "unaligned"})
                    continue
                resi = m[resi3]
                moved = "transferred"
            out.append({"paralog": p, "structure_resi": resi3, "resi": resi,
                        "shell": s["shell"],
                        "median_distance_A": s["median_distance_A"],
                        "n_structures_contact": n_contact,
                        "is_s0_contact": s["is_s0_contact"],
                        "is_consensus_contact": n_contact >= CONSENSUS_MAJORITY,
                        "transfer": moved})
    return out


# ---------------------------------------------------------------------------

def contact_test(mods: list[dict], shells: list[dict]) -> list[dict]:
    """Contacts against two backgrounds, by permutation.

    `rest_of_core` is the brief's question: are the contacts special inside
    the binding core.  `rest_of_pocket` is the harder one and the reason
    stage 2 measured shells at all — the contacts against every *other*
    residue the structure puts within 15 A of the ligand.  A test that is
    significant against the core and null against the pocket says the
    constrained thing is the pocket, not the ten residues.
    """
    rows: list[dict] = []
    for label, key in (("s0_contact", "is_s0_contact"),
                       ("consensus_contact", "is_consensus_contact")):
      for background in ("rest_of_core", "rest_of_pocket"):
        for p in L.PARALOGS:
            pocket = {int(s["resi"]) for s in shells
                      if s["paralog"] == p and s["resi"] is not None}
            core = (M.residues(mods, p, M.PRIMARY["ligand_core"])
                    if background == "rest_of_core" else pocket)
            marked = {int(s["resi"]) for s in shells
                      if s["paralog"] == p and s["resi"] is not None
                      and str(s[key]) == "True"}
            by_resi = {r["resi"]: r for r in L.constraint(p)}
            for layer in L.LAYERS:
                pos = sorted(core)
                vals = [by_resi[i][f"{layer}_jsd"] for i in pos if i in by_resi]
                labs = [i in marked for i in pos if i in by_resi]
                perm = L.permutation_label(vals, labs, iters=PERM_ITERS)
                inside = [v for v, l in zip(vals, labs) if l and v is not None]
                outside = [v for v, l in zip(vals, labs) if not l and v is not None]
                rows.append({
                    "contact_set": label, "background": background,
                    "paralog": p, "layer": layer,
                    "n_background": len(pos), "n_contacts": sum(labs),
                    "mean_contact_jsd": L.fmt(L.mean(inside)),
                    "mean_background_jsd": L.fmt(L.mean(outside)),
                    "difference": L.fmt((L.mean(inside) or 0) - (L.mean(outside) or 0)),
                    "p_permutation": L.fmt(perm["p"], 6),
                    "iters": perm.get("iters", 0),
                    "n_invariant_contacts": sum(
                        1 for i in sorted(marked)
                        if i in by_resi and by_resi[i].get(f"{layer}_jsd") is not None
                        and by_resi[i].get("deep_frac_modal") not in ("", None)
                        and float(by_resi[i]["deep_frac_modal"]) >= 0.999),
                })
    qs = L.benjamini_hochberg([float(r["p_permutation"]) if r["p_permutation"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_permutation"] = L.fmt(q, 6)
    return rows


def shell_constraint(shells: list[dict]) -> list[dict]:
    """Constraint by shell, and the trend across them."""
    rows: list[dict] = []
    for p in L.PARALOGS:
        by_resi = {r["resi"]: r for r in L.constraint(p)}
        mine = [s for s in shells if s["paralog"] == p and s["resi"] is not None]
        prot = [r["deep_jsd"] for r in L.constraint(p) if r["deep_jsd"] is not None]
        for shell in SHELL_ORDER:
            vals, fmods = [], []
            for s in mine:
                if s["shell"] != shell:
                    continue
                r = by_resi.get(int(s["resi"]))
                if r is None:
                    continue
                if r["deep_jsd"] is not None:
                    vals.append(r["deep_jsd"])
                if r.get("deep_frac_modal") not in ("", None):
                    fmods.append(float(r["deep_frac_modal"]))
            mw = L.mann_whitney(vals, prot, alternative="greater")
            rows.append({
                "paralog": p, "shell": shell, "n_residues": len(vals),
                "mean_jsd": L.fmt(L.mean(vals)),
                "median_jsd": L.fmt(L.median(vals)),
                "mean_frac_modal": L.fmt(L.mean(fmods)),
                "whole_protein_mean_jsd": L.fmt(L.mean(prot)),
                "p_greater_than_protein": L.fmt(mw["p"], 6),
            })
    qs = L.benjamini_hochberg([float(r["p_greater_than_protein"])
                               if r["p_greater_than_protein"] else None
                               for r in rows])
    for r, q in zip(rows, qs):
        r["q_greater_than_protein"] = L.fmt(q, 6)
    return rows


def shell_trend(shells: list[dict]) -> list[dict]:
    """Spearman of constraint against distance, over the shell residues.

    Reported per paralogue on the residues the structure actually resolves.
    A negative rho is the decay the ligand hypothesis predicts: further from
    IP3, less constrained.
    """
    from s15_lib import spearman
    rows: list[dict] = []
    for p in L.PARALOGS:
        by_resi = {r["resi"]: r for r in L.constraint(p)}
        xs, ys, fm = [], [], []
        for s in shells:
            if s["paralog"] != p or s["resi"] is None:
                continue
            r = by_resi.get(int(s["resi"]))
            if r is None or r["deep_jsd"] is None:
                continue
            xs.append(float(s["median_distance_A"]))
            ys.append(r["deep_jsd"])
            fm.append(float(r["deep_frac_modal"]) if r.get("deep_frac_modal")
                      not in ("", None) else None)
        rho, pv, n = spearman(xs, ys)
        pairs = [(a, b) for a, b in zip(xs, fm) if b is not None]
        rho2, pv2, n2 = spearman([a for a, _ in pairs], [b for _, b in pairs])
        rows.append({
            "paralog": p, "n_residues": n,
            "rho_distance_vs_jsd": L.fmt(rho),
            "p_distance_vs_jsd": L.fmt(pv, 6),
            "rho_distance_vs_frac_modal": L.fmt(rho2),
            "p_distance_vs_frac_modal": L.fmt(pv2, 6),
            "direction": ("constraint falls with distance" if (rho or 0) < 0
                          else "constraint rises with distance"),
        })
    return rows


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = [("shell index", False), ("contact test", False),
             ("shell constraint", False), ("trend", False)]
    L.live("contacts", steps)
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")

    shells = shell_index()
    L.write_tsv(L.OUT_DIR / "shell_index.tsv", shells,
                ["paralog", "structure_resi", "resi", "shell",
                 "median_distance_A", "n_structures_contact", "is_s0_contact",
                 "is_consensus_contact", "transfer"])
    n_un = sum(1 for s in shells if s["transfer"] == "unaligned")
    L.log(f"shell index: {len(shells)} rows, {n_un} unaligned on transfer")
    steps[0] = ("shell index", True); L.live("contacts", steps)

    ct = contact_test(mods, shells)
    L.write_tsv(L.OUT_DIR / "contact_test.tsv", ct, list(ct[0]))
    steps[1] = ("contact test", True); L.live("contacts", steps)

    sc = shell_constraint(shells)
    L.write_tsv(L.OUT_DIR / "shell_constraint.tsv", sc, list(sc[0]))
    steps[2] = ("shell constraint", True); L.live("contacts", steps)

    tr = shell_trend(shells)
    L.write_tsv(L.OUT_DIR / "shell_trend.tsv", tr, list(tr[0]))
    steps[3] = ("trend", True); L.live("contacts", steps)

    for r in ct:
        if r["contact_set"] == "s0_contact" and r["layer"] == "deep":
            L.log(f"{r['paralog']} contacts {r['mean_contact_jsd']} vs "
                  f"{r['background']} {r['mean_background_jsd']}  "
                  f"q={r['q_permutation']}")
    for r in tr:
        L.log(f"{r['paralog']} rho(distance, jsd) = {r['rho_distance_vs_jsd']} "
              f"p={r['p_distance_vs_jsd']} n={r['n_residues']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
