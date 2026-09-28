"""S29 — the gating transition at the pore, across the type-3 state panel.

S0 measured the pore on one structure, 6DQN, and found a wide filter
(5.08 Å) and a narrow gate (2.55 Å) (the review's Figure 3c). S11 then
chose one ITPR3 deposition per conformational state by stated rules. This
stage measures every one of them with S0's own code (`s0_figdata_structure.
measure_pore`: the same axis, orientation, profile and constriction rule),
so the question "does the pore change between states, and where" is
answered by one instrument rather than read from six papers.

Structures: 6DQN (S0's, IP3-bound) and the six ITPR3 rows of S11's
`structure_manifest.tsv` (resting 8TKG, activated 8TKF, apo 6DQJ,
higher-order inhibited 8TLA, labile resting 8TKH, preactivated 7T3P), as
full tetramers from the data root (S11's mmCIF cache).

Outputs (results/structures/):
    state_pore_profiles.tsv   each deposit's profile, z also given from its filter
    state_constrictions.tsv   filter and gate: z, radius, lining residues

Waters are removed first (see `without_water`). Positive control (raises):
6DQN must reproduce S0's committed profile and
constrictions exactly, or this is not S0's instrument. Each deposit's filter
must also land on the GGGVGD motif, as S0 requires of 6DQN.

The companion simulator measures the same panel with its own code
(`ip3r_simulation`, `structure.states`); its gate and filter radii agree
with these to the rounding (docs/s29_simulator_review.md, Figure S29.5).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s0_figdata_structure as S0  # noqa: E402
import s11_lib as S11  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "structures"
MANIFEST = OUT / "structure_manifest.tsv"
S0_DIR = ROOT / "results" / "s0_baseline" / "review_figures"
S0_PDB, S0_STATE = "6DQN", "IP3-bound"


class StatePoreRefusal(RuntimeError):
    pass


def panel() -> list[tuple[str, str]]:
    """(pdb_id, state): S0's structure, then S11's ITPR3 depositions."""
    with MANIFEST.open(encoding="utf-8") as fh:
        rows = [r for r in csv.DictReader(fh, delimiter="\t")
                if r["paralog"] == "ITPR3" and r["source"] == "PDB"]
    return [(S0_PDB, S0_STATE)] + [(r["source_id"], r["state"]) for r in rows]


def without_water(st: dict) -> dict:
    """Drop modelled waters. S0's rule counts every heavy atom, which on 6DQN
    (no waters in the pore) is the protein wall; 8TKG's 2.5 Å map models
    waters inside the pore, and its luminal minimum then lands on them, off
    the filter motif (the motif check below refused it). The wall is the
    protein's, so waters are removed; 6DQN is unchanged (the S0 control)."""
    keep = [i for i, (_c, _s, comp) in enumerate(st["heavy_id"]) if comp != "HOH"]
    return {**st, "heavy": st["heavy"][keep],
            "heavy_id": [st["heavy_id"][i] for i in keep]}


def measure(pdb_id: str) -> dict:
    cif = S11.structures_dir() / "reference" / f"{pdb_id.lower()}.cif"
    if not cif.exists():
        cif = S11.pdb_fetch(pdb_id)
    st = without_water(S0.parse(cif))
    axis, centre, _ = S0.four_fold_frame(st["ca"])
    pore = S0.measure_pore(st, axis, centre)
    heavy_f = pore["heavy_f"]
    (fz, fr), (gz, gr) = pore["filter"], pore["gate"]
    filt = S0.lining_residues(heavy_f, st["heavy_id"], fz, fr)
    gate = S0.lining_residues(heavy_f, st["heavy_id"], gz, gr)
    motif = S0.motif_position()
    nums = [int("".join(c for c in r if c.isdigit())) for r in filt]
    if motif and not any(motif - 2 <= n <= motif + len(S0.FILTER_MOTIF)
                         for n in nums):
        raise StatePoreRefusal(f"{pdb_id}: the luminal constriction ({filt}) "
                               f"is not on the {S0.FILTER_MOTIF} filter")
    return {"profile": pore["profile"], "filter": (fz, fr), "gate": (gz, gr),
            "filter_lining": filt, "gate_lining": gate}


def s0_control(m: dict) -> None:
    """6DQN through the shared code must be S0's committed measurement."""
    with (S0_DIR / "structure_pore.tsv").open(encoding="utf-8") as fh:
        ref = [(float(r["z_along_axis"]), float(r["min_heavy_atom_radius"]))
               for r in csv.DictReader(fh, delimiter="\t")]
    got = [(round(z, 2), round(r, 3)) for z, r in m["profile"]]
    meta = json.loads((S0_DIR / "structure_meta.json").read_text())
    same = (got == ref
            and round(m["filter"][1], 2) == meta["filter_min_radius_A"]
            and round(m["gate"][1], 2) == meta["gate_min_radius_A"])
    if not same:
        raise StatePoreRefusal("6DQN does not reproduce S0's committed pore "
                               "profile and constrictions")


def main() -> int:
    prof_rows, cons_rows = [], []
    for pdb_id, state in panel():
        m = measure(pdb_id)
        if pdb_id == S0_PDB:
            s0_control(m)
        fz = m["filter"][0]
        for z, r in m["profile"]:
            prof_rows.append({"pdb_id": pdb_id, "state": state,
                              "z_A": f"{z:.2f}", "z_from_filter_A": f"{z - fz:.2f}",
                              "min_heavy_atom_radius_A": f"{r:.3f}"})
        cons_rows.append({
            "pdb_id": pdb_id, "state": state,
            "filter_z_A": f"{m['filter'][0]:.1f}",
            "filter_radius_A": f"{m['filter'][1]:.2f}",
            "gate_z_A": f"{m['gate'][0]:.1f}",
            "gate_radius_A": f"{m['gate'][1]:.2f}",
            "gate_above_filter_A": f"{m['gate'][0] - fz:.1f}",
            "filter_lining": ";".join(m["filter_lining"]),
            "gate_lining": ";".join(m["gate_lining"])})
        print(f"[s29] {pdb_id} {state}: filter {m['filter'][1]:.2f} Å, "
              f"gate {m['gate'][1]:.2f} Å ({', '.join(m['gate_lining'])})")
    for name, rows in (("state_pore_profiles.tsv", prof_rows),
                       ("state_constrictions.tsv", cons_rows)):
        with (OUT / name).open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, list(rows[0]), delimiter="\t",
                               lineterminator="\n")
            w.writeheader()
            w.writerows(rows)
    print(f"[s29] {len(cons_rows)} depositions -> results/structures/"
          "state_pore_profiles.tsv, state_constrictions.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
