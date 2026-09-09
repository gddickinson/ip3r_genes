"""S22 stage 2 — how far every residue is from the ligand, measured.

S17 knows ten residues contact IP3, because S0 measured them on 6DQN at
<= 4.5 A.  Ten positions is a thin basis for a constraint claim, and a
binary contact/not-contact label throws away the one thing a structure can
say that an alignment cannot: constraint should *decay* with distance from
the ligand if the ligand is what holds those residues in place.  So this
stage measures the whole gradient.

Three things make the measurement checkable.

**All-atom, not CA.**  A CA-only trace puts an arginine's CA 8 A from a
phosphate its guanidinium is hydrogen-bonded to.  The reader here is
deliberately small — S11 wrote one for CA coordinates and this is the same
idea for the atoms a distance needs.

**Five independent depositions, not one.**  6DQN is S0's structure; 8TKG,
8TKF, 8TKH, 7T3P and 8TLA are four further human ITPR3 depositions with IP3
bound, from a different group and different conformational states.  A shell
assignment that only holds in one of them is a fact about that map.  The
consensus is committed with the per-structure disagreement beside it.

**S0's contact set is the positive control.**  The <= 4.5 A set recovered
here must contain S0's ten residues in the structure S0 used, or the reader
is wrong and everything downstream of it is worthless.  That check raises.
"""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s11_lib as S11  # noqa: E402

LIGAND = "I3P"                       # inositol 1,4,5-trisphosphate
S0_STRUCTURE = "6DQN"
REPLICATES = ("8TKG", "8TKF", "8TKH", "7T3P", "8TLA")

# Shell edges in angstrom.  4.5 is S0's contact cutoff, kept so the first
# shell is that set and not a new one; the rest are equal 3.5 A steps out to
# 15 A, past which no side chain of a 300-residue module can reach a ligand.
SEARCH_RADIUS_A = 15.0
SHELL_EDGES = ((0.0, 4.5, "contact"),
               (4.5, 8.0, "second"),
               (8.0, 11.5, "third"),
               (11.5, 15.0, "fourth"))

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


class ShellRefusal(RuntimeError):
    pass


# ---------------------------------------------------------------------------
# a deliberately small all-atom mmCIF reader
# ---------------------------------------------------------------------------

def read_atoms(path: Path) -> tuple[list[tuple], list[tuple]]:
    """(protein atoms, ligand atoms) as (chain, resi, resn, x, y, z).

    Only the first model is read and only altloc `.`/`A` is kept, for the
    reason S11's reader gives: a second model or a B conformer would enter
    the distance minimum as if it were a second copy of the residue.
    """
    prot: list[tuple] = []
    lig: list[tuple] = []
    cols: dict[str, int] = {}
    in_loop = False
    with open(path) as fh:
        for line in fh:
            if line.startswith("_atom_site."):
                cols[line.strip().split(".", 1)[1]] = len(cols)
                in_loop = True
                continue
            if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
                f = line.split()
                if len(f) < len(cols):
                    continue
                if f[cols["pdbx_PDB_model_num"]] not in ("1", "."):
                    continue
                alt = f[cols["label_alt_id"]]
                if alt not in (".", "?", "A"):
                    continue
                resn = f[cols["auth_comp_id"]]
                try:
                    xyz = (float(f[cols["Cartn_x"]]), float(f[cols["Cartn_y"]]),
                           float(f[cols["Cartn_z"]]))
                except ValueError:
                    continue
                chain = f[cols["auth_asym_id"]]
                try:
                    resi = int(f[cols["auth_seq_id"]])
                except ValueError:
                    continue
                rec = (chain, resi, resn) + xyz
                if resn == LIGAND:
                    lig.append(rec)
                elif resn in AA3:
                    prot.append(rec)
            elif in_loop and line.startswith("#"):
                if prot or lig:
                    break
    return prot, lig


def _grid(atoms: list[tuple], cell: float):
    g = defaultdict(list)
    for a in atoms:
        key = (int(a[3] // cell), int(a[4] // cell), int(a[5] // cell))
        g[key].append(a)
    return g


def min_distances(prot: list[tuple], lig: list[tuple],
                  cutoff: float = 15.0) -> dict[tuple[str, int], dict]:
    """Per protein residue: min distance to the ligand, same chain and any.

    A uniform grid over the ligand atoms keeps this linear in protein atoms;
    a 11,000-residue tetramer against 132 ligand atoms is otherwise 5.8 M
    pair distances per structure and five structures of it.
    """
    cell = cutoff
    grid = _grid(lig, cell)
    out: dict[tuple[str, int], dict] = {}
    for ch, resi, resn, x, y, z in prot:
        key = (ch, resi)
        rec = out.get(key)
        if rec is None:
            rec = out[key] = {"resn": resn, "d_same": math.inf,
                              "d_any": math.inf, "lig_chain_any": ""}
        gx, gy, gz = int(x // cell), int(y // cell), int(z // cell)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for lc, _li, _ln, lx, ly, lz in grid.get((gx + dx, gy + dy, gz + dz), ()):
                        d = math.dist((x, y, z), (lx, ly, lz))
                        if d < rec["d_any"]:
                            rec["d_any"] = d
                            rec["lig_chain_any"] = lc
                        if lc == ch and d < rec["d_same"]:
                            rec["d_same"] = d
    return out


def shell_of(d: float) -> str:
    """Shell name, or `beyond_radius` past the search radius.

    Residues past `SEARCH_RADIUS_A` are absent from `ligand_shells.tsv`
    entirely rather than pooled into a last open-ended bin: the grid search
    scores a residue at 20 A only when it happens to land in a neighbouring
    cell, so an open bin would be a population defined by the geometry of
    the search and not by the protein.
    """
    for lo, hi, name in SHELL_EDGES:
        if lo <= d < hi:
            return name
    return "beyond_radius"


# ---------------------------------------------------------------------------
# the stage
# ---------------------------------------------------------------------------

def measure(pdb_id: str) -> dict[int, dict]:
    """Per residue number of human ITPR3, the best distance over subunits.

    The tetramer carries four equivalent sites, so a residue's distance is
    the *minimum over its own chain's* ligand where one is present — the
    same residue in a chain whose site is empty is not evidence of anything.
    """
    path = S11.pdb_fetch(pdb_id)
    prot, lig = read_atoms(path)
    if not lig:
        raise ShellRefusal(f"{pdb_id}: no {LIGAND} atoms found")
    d = min_distances(prot, lig, cutoff=SEARCH_RADIUS_A)
    per_resi: dict[int, dict] = {}
    lig_chains = {a[0] for a in lig}
    for (ch, resi), rec in d.items():
        if rec["d_same"] > SEARCH_RADIUS_A:
            continue
        cur = per_resi.get(resi)
        if cur is None or rec["d_same"] < cur["d_same"]:
            per_resi[resi] = {"resn": rec["resn"], "d_same": rec["d_same"],
                              "d_any": rec["d_any"], "chain": ch}
    L.log(f"{pdb_id}: {len(prot)} protein atoms, {len(lig)} {LIGAND} atoms in "
          f"{len(lig_chains)} chain(s), {len(per_resi)} residues scored")
    return per_resi


def s0_contacts() -> list[int]:
    return sorted(r["source_resi"] for r in L.functional_sites()
                  if r["paralog"] == L.STRUCTURE_REF
                  and r["site_class"] == "ip3_contact")


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    L.live("shells", [("measure structures", False), ("consensus", False)])
    want = s0_contacts()
    per_struct: dict[str, dict[int, dict]] = {}
    rows: list[dict] = []
    for pdb_id in (S0_STRUCTURE,) + REPLICATES:
        m = measure(pdb_id)
        per_struct[pdb_id] = m
        got = sorted(r for r, v in m.items() if v["d_same"] <= 4.5)
        missing = [r for r in want if r not in got]
        note = "" if not missing else "missing:" + ";".join(map(str, missing))
        if pdb_id == S0_STRUCTURE and missing:
            raise ShellRefusal(
                f"{S0_STRUCTURE}: the reader does not recover S0's contact "
                f"set — missing {missing}")
        rows.append({
            "pdb_id": pdb_id, "search_radius_A": SEARCH_RADIUS_A,
            "n_residues_scored": len(m),
            "n_contact_le_4.5A": len(got),
            "n_s0_contacts_recovered": len(want) - len(missing),
            "n_s0_contacts": len(want),
            "extra_contacts": ";".join(str(r) for r in got if r not in want),
            "note": note,
        })
    L.write_tsv(L.OUT_DIR / "shell_agreement.tsv", rows,
                ["pdb_id", "search_radius_A", "n_residues_scored",
                 "n_contact_le_4.5A",
                 "n_s0_contacts_recovered", "n_s0_contacts",
                 "extra_contacts", "note"])
    L.live("shells", [("measure structures", True), ("consensus", False)])

    # consensus: the median distance over the depositions that resolve the
    # residue, and how many of them put it in the contact shell.
    all_resi = sorted({r for m in per_struct.values() for r in m})
    out: list[dict] = []
    for resi in all_resi:
        ds = [m[resi]["d_same"] for m in per_struct.values() if resi in m]
        med = L.median(ds)
        n_contact = sum(1 for d in ds if d <= 4.5)
        aa = next((m[resi]["resn"] for m in per_struct.values() if resi in m), "")
        out.append({
            "resi": resi, "aa3": aa, "aa": AA3.get(aa, "X"),
            "n_structures": len(ds),
            "min_distance_A": min(ds), "median_distance_A": med,
            "max_distance_A": max(ds),
            "shell": shell_of(med),
            "n_structures_contact": n_contact,
            "is_s0_contact": resi in want,
            "distances": ";".join(f"{p}:{m[resi]['d_same']:.2f}"
                                  for p, m in per_struct.items() if resi in m),
        })
    L.write_tsv(L.OUT_DIR / "ligand_shells.tsv", out,
                ["resi", "aa3", "aa", "n_structures", "min_distance_A",
                 "median_distance_A", "max_distance_A", "shell",
                 "n_structures_contact", "is_s0_contact", "distances"])
    by_shell: dict[str, int] = defaultdict(int)
    for r in out:
        by_shell[r["shell"]] += 1
    L.log("consensus shells: " + "  ".join(
        f"{n}={by_shell[n]}" for _, _, n in SHELL_EDGES))
    L.live("shells", [("measure structures", True), ("consensus", True)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
