"""S17 stage 6 — constraint painted onto the S11 structures.

A per-site table is not a structure-mapped result until it is on a structure.
This writes the per-site score into the **B-factor column** of a CA trace for
each ITPR structure in the S11 panel, which makes the result usable by anyone
with PyMOL or ChimeraX and no code:

    load constraint_reference_ITPR3_8TKG.pdb
    spectrum b, blue_white_red, minimum=40, maximum=100

Two layers are written per structure, on **one scale read the same way round**
so they can be coloured with the same command: `constraint` is `deep_jsd x 100`,
and `selection` is the FEL non-synonymous rate inverted, `(1 - min(beta, 1)) x
100`, so 100 means "no non-synonymous change" in both. Per-site omega is
deliberately not painted — the FEL stage documents why alpha is unidentifiable
at a share of sites, and a ratio with an unidentified denominator would be
painted as a colour.

Residues with no scored counterpart — unmodelled loops, unmapped termini — are
written as **-1, never 0**, so an unscored residue cannot be read as an
unconstrained one. The luminal loop is the case that makes this matter: it is
the least conserved element in the protein *and* the one the cryo-EM maps
resolve worst, so the two must be distinguishable on the coloured structure.

Mapping is by pairwise alignment of each structure's own sequence (read out of
its CA records) to the reference protein, so a rat ITPR1 cryo-EM structure with
its own numbering and its own disordered gaps is handled without any assumption
that it shares residue numbers with human ITPR1.

    python scripts/s17_paint.py
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_lib as S11                                          # noqa: E402
import s11_struct_io as SIO                                    # noqa: E402
import s17_lib as L                                            # noqa: E402

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M", "SEC": "U",
}

#: Minimum share of a structure's residues that must map to the reference for
#: the painting to be written. A model that maps a tenth of itself is being
#: painted with somebody else's profile.
MIN_MAPPED_FRAC = 0.50


def targets() -> list[dict]:
    """Every ITPR structure in S11's panel, from its committed manifest.

    Read, never listed: S11 chose the panel by seven stated rules and a second
    hand-written list here could silently disagree with it. Rows with no file
    on disk are reported as unpainted rather than skipped — S11's human ITPR2
    AFDB entry is the 181-residue isoform and has no path, and that absence is
    a result S11 already reported.
    """
    out = []
    for r in L.read_tsv(L.STRUCT_DIR / "structure_manifest.tsv"):
        if r["call"] != "ITPR" or r["paralog"] not in L.PARALOGS:
            continue
        out.append(r)
    return out


def chain_sequence(chain) -> tuple[str, list]:
    residues = [x for x in chain.residues if x.name.upper() in THREE_TO_ONE]
    return "".join(THREE_TO_ONE[x.name.upper()] for x in residues), residues


def _scores(out_dir: Path) -> tuple[dict, dict]:
    constraint: dict[str, dict[int, float]] = {}
    for paralog in L.PARALOGS:
        acc = L.REFERENCES[paralog][1]
        constraint[paralog] = {
            int(r["resi"]): float(r["deep_jsd"])
            for r in L.read_tsv(out_dir / f"constraint_{paralog}_{acc}.tsv")
            if r["deep_jsd"] not in ("", None) and r["deep_reliable"] == "True"}
    selection: dict[str, dict[int, float]] = {p: {} for p in L.PARALOGS}
    fel = out_dir / "fel_sites.tsv"
    if fel.exists():
        for r in L.read_tsv(fel):
            if not r["resi"]:
                continue
            selection[r["paralog"]][int(r["resi"])] = max(
                0.0, 1.0 - min(float(r["beta"]), 1.0))
    return constraint, selection


def paint(out_dir: Path = L.OUT_DIR) -> list[dict]:
    dest_dir = out_dir / "painted"
    dest_dir.mkdir(parents=True, exist_ok=True)
    constraint, selection = _scores(out_dir)

    rows = []
    for t in targets():
        label, paralog, path = t["id"], t["paralog"], t["path"]
        if not path or not Path(path).exists():
            rows.append({"structure": label, "layer": "", "paralog": paralog,
                         "role": t["role"], "file": "", "n_residues": 0,
                         "n_painted": 0, "frac_painted": 0.0, "n_mapped": 0,
                         "note": t.get("note") or "no structure file on disk"})
            continue
        chain = SIO.largest_chain(SIO.read_structure(Path(path)))
        seq, residues = chain_sequence(chain)
        ref_seq = L.uniprot_fasta(L.REFERENCES[paralog][1])
        xfer = L.transfer_positions(seq, ref_seq)
        frac_mapped = len(xfer) / len(residues) if residues else 0.0
        if frac_mapped < MIN_MAPPED_FRAC:
            rows.append({"structure": label, "layer": "", "paralog": paralog,
                         "role": t["role"], "file": "", "n_residues": len(residues),
                         "n_painted": 0, "frac_painted": 0.0, "n_mapped": len(xfer),
                         "note": f"only {frac_mapped:.0%} of the chain maps to "
                                 f"{paralog} — refused"})
            continue
        for layer, table, note in (
                ("constraint", constraint, "b-factor = deep_jsd x 100"),
                ("selection", selection, "b-factor = (1 - min(dN rate, 1)) x 100")):
            for res in chain.residues:
                res.bfactor = -1.0
            painted = 0
            for i, res in enumerate(residues, start=1):
                resi = xfer.get(i)
                val = table[paralog].get(resi) if resi else None
                if val is not None:
                    res.bfactor = round(val * 100, 2)
                    painted += 1
            dest = dest_dir / f"{layer}_{label}.pdb"
            SIO.write_ca_pdb(chain, dest)
            rows.append({
                "structure": label, "layer": layer, "paralog": paralog,
                "role": t["role"], "file": dest.name, "n_residues": len(residues),
                "n_painted": painted, "n_mapped": len(xfer),
                "frac_painted": round(painted / len(residues), 4) if residues else 0.0,
                "note": note + "; -1 = no score",
            })
        done = {r["layer"]: r["n_painted"] for r in rows if r["structure"] == label}
        print(f"[s17] painted {label}: constraint {done.get('constraint', 0)}, "
              f"selection {done.get('selection', 0)} of {len(residues)} residues")
    L.write_tsv(out_dir / "painted_structures.tsv", rows)
    return rows


def run(out_dir: Path = L.OUT_DIR) -> dict:
    rows = paint(out_dir)
    return {"n_structures": len({r["structure"] for r in rows if r["n_painted"]})}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    run(ap.parse_args().out)


if __name__ == "__main__":
    main()
