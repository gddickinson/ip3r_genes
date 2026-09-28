"""S29 — the IP3 contact control under both contact rules.

S22's positive control is that all six IP3-bound depositions recover S0's
ten IP3 contacts at <= 4.5 A.  S22 measures over *all* atoms, and a
deposition that models hydrogens then measures some distances to a
hydrogen.  S0 defined the ten on 6DQN, which models none, so S0's rule is
the heavy-atom one.  The companion simulator (`../ip3r_simulation`, check
`P6.contacts_heavy_atom`) found that the two rules disagree in two
depositions.  This stage measures that with this project's own reader,
from this project's own structure files, so the papers can state it from a
committed table.

Output: `results/ligand_site/contact_rule.tsv`, one row per deposition and
S0 contact: the minimum distance to the chain's own IP3 over all atoms and
over heavy atoms, and whether each is inside the cutoff.

Positive control (raises): the all-atom column must reproduce S22's
committed `shell_agreement.tsv` count of recovered contacts in every
deposition, or this reader is not S22's and nothing here means anything.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_shells as SH  # noqa: E402

CUTOFF_A = 4.5                       # S0's contact cutoff, as S22 uses it
OUT = L.OUT_DIR / "contact_rule.tsv"
FIELDS = ["pdb_id", "resi", "aa3", "d_all_atoms_A", "d_heavy_atoms_A",
          "contact_all_atoms", "contact_heavy_atoms"]


class ContactRuleRefusal(RuntimeError):
    pass


def committed_recovery() -> dict[str, int]:
    return {r["pdb_id"]: int(r["n_s0_contacts_recovered"])
            for r in L.read_tsv(L.OUT_DIR / "shell_agreement.tsv")}


def rows_for(pdb_id: str, want: list[int]) -> list[dict]:
    every = SH.measure(pdb_id)
    heavy = SH.measure(pdb_id, heavy_only=True)
    out = []
    for resi in want:
        a = every[resi]["d_same"] if resi in every else float("inf")
        h = heavy[resi]["d_same"] if resi in heavy else float("inf")
        out.append({"pdb_id": pdb_id, "resi": resi,
                    "aa3": (every.get(resi) or heavy.get(resi) or {}).get("resn", ""),
                    "d_all_atoms_A": round(a, 3), "d_heavy_atoms_A": round(h, 3),
                    "contact_all_atoms": a <= CUTOFF_A,
                    "contact_heavy_atoms": h <= CUTOFF_A})
    return out


def main() -> int:
    want = SH.s0_contacts()
    committed = committed_recovery()
    rows: list[dict] = []
    for pdb_id in (SH.S0_STRUCTURE,) + SH.REPLICATES:
        got = rows_for(pdb_id, want)
        n_all = sum(r["contact_all_atoms"] for r in got)
        if n_all != committed[pdb_id]:
            raise ContactRuleRefusal(
                f"{pdb_id}: all-atom recovery {n_all} differs from S22's "
                f"committed {committed[pdb_id]}")
        n_heavy = sum(r["contact_heavy_atoms"] for r in got)
        L.log(f"{pdb_id}: {n_all}/{len(want)} all atoms, {n_heavy}/{len(want)} "
              "heavy atoms")
        rows.extend(got)
    L.write_tsv(OUT, rows, FIELDS)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
