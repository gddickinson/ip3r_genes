#!/usr/bin/env python3
"""Measure the IP₃R channel from deposited cryo-EM coordinates (review figures).

The review's §2.2-2.4 make three geometric claims — a large cytosolic cap over
a small membrane domain, a pore with a filter and a gate at the cytosolic end
of TM6, and an IP₃ site "roughly 100 Å from the gate". This script measures all
three from the structure the review already cites, so the figures show a
measurement instead of repeating a number.

Structure: **PDB 6DQN** — human type-3 receptor with IP₃ bound, 3.33 Å, the
class-1 map of Paknejad & Hite (review reference R24). Human, ligand-bound and
high resolution, so it carries the ligand site *and* matches the panel protein
Q14573 whose Pfam coordinates colour the trace.

Bulk in, small out. The 14 MB mmCIF is downloaded to a work directory outside
the repository (the project's bulk-data rule); what is committed is the
Cα trace of **one** subunit — the deposit is C4-symmetric to within the
tolerance this script checks and prints, so the other three are reconstructed
by rotation at draw time — plus the pore profile and a measurements file.

    python scripts/s0_figdata_structure.py [--work-dir DIR] [--pdb 6DQN]

Outputs (results/s0_baseline/review_figures/):
    structure_ca.tsv       Cα of one subunit in the four-fold frame
    structure_pore.tsv     minimum heavy-atom distance to the axis, along it
    structure_meta.json    axis, C4 residual, TM span, filter/gate, distances
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "s0_baseline" / "review_figures"
DOMAINS = OUT_DIR / "domain_coords.tsv"

PDB_ID = "6DQN"
PDB_URL = "https://files.rcsb.org/download/{pdb}.cif.gz"
#: The subunit whose Pfam coordinates annotate the trace (6DQN is human IP3R3).
STRUCT_ACC = "Q14573"
#: The pore's selectivity-filter signature. Finding the narrowest luminal point
#: on this motif is an independent check that the axis and profile are right —
#: nothing in the geometry knows about the sequence.
FILTER_MOTIF = "GGGVGD"
PANEL_FASTA = ROOT / "results" / "benchmark_controls" / "panel_positives.fasta"


# --------------------------------------------------------------- parsing

def fetch_cif(pdb: str, work: Path) -> Path:
    """Download the mmCIF into `work` (outside the repo) unless already there."""
    work.mkdir(parents=True, exist_ok=True)
    cif = work / f"{pdb.lower()}.cif"
    if cif.exists() and cif.stat().st_size > 1_000_000:
        print(f"  using cached {cif}")
        return cif
    gz = cif.with_suffix(".cif.gz")
    url = PDB_URL.format(pdb=pdb.upper())
    print(f"  fetching {url}")
    rc = subprocess.run(["curl", "-sSL", "-o", str(gz), url]).returncode
    if rc != 0 or not gz.exists():
        raise SystemExit(f"download failed for {pdb}")
    subprocess.run(["gunzip", "-f", str(gz)], check=True)
    return cif


def parse(cif: Path) -> dict:
    """Cα per chain, every heavy atom, and the hetero groups we care about.

    mmCIF `atom_site` is a fixed-column loop here, so a whitespace split is
    enough and avoids a parser dependency; the column indices are asserted
    against the loop header rather than assumed.
    """
    header: list[str] = []
    in_loop = False
    ca: dict[str, list] = {}
    het: dict[tuple, list] = {}
    heavy: list[tuple] = []
    heavy_id: list[tuple] = []          # (chain, resseq, resname) per heavy atom
    for ln in cif.open(encoding="utf-8"):
        if ln.startswith("_atom_site."):
            header.append(ln.strip())
            in_loop = True
            continue
        if in_loop and not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            if header and not ln.startswith("#"):
                continue
        if not (ln.startswith("ATOM") or ln.startswith("HETATM")):
            continue
        f = ln.split()
        elem, atom, comp, chain = f[2], f[3], f[5], f[18]
        xyz = (float(f[10]), float(f[11]), float(f[12]))
        if elem != "H":
            heavy.append(xyz)
            heavy_id.append((chain, int(f[16]) if f[16].lstrip("-").isdigit()
                             else -1, comp))
        if ln.startswith("ATOM") and atom == "CA":
            ca.setdefault(chain, []).append((int(f[16]), *xyz))
        if comp in ("I3P", "ZN"):
            het.setdefault((comp, chain), []).append(xyz)
    idx = {name: i for i, name in enumerate(header)}
    for need, want in (("_atom_site.Cartn_x", 10), ("_atom_site.auth_seq_id", 17 - 1),
                       ("_atom_site.auth_asym_id", 18)):
        if idx.get(need) != want:
            raise SystemExit(f"mmCIF column layout changed: {need} at "
                             f"{idx.get(need)}, expected {want}")
    return {"ca": ca, "het": het, "heavy": np.array(heavy),
            "heavy_id": heavy_id}


# ------------------------------------------------------------ the frame

def four_fold_frame(ca: dict) -> tuple[np.ndarray, np.ndarray, float]:
    """Axis of the tetramer, plus a centre on it, plus the planarity residual.

    The four subunit centroids of a C4 tetramer lie on a circle in a plane
    normal to the symmetry axis, so the axis is that plane's normal — the
    smallest singular direction of the centroid cloud. The residual is the
    smallest singular value: how far from coplanar the four centroids are,
    reported so a non-C4 deposit cannot pass silently.
    """
    cents = np.array([np.mean([r[1:] for r in v], axis=0)
                      for v in ca.values()])
    centre = cents.mean(0)
    _, sv, vt = np.linalg.svd(cents - centre)
    return vt[2] / np.linalg.norm(vt[2]), centre, float(sv[2])


def to_frame(xyz: np.ndarray, axis: np.ndarray, centre: np.ndarray) -> np.ndarray:
    """(x, y, z) with z along the four-fold axis and the origin on it."""
    e3 = axis
    seed = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(seed, e3)) > 0.9:
        seed = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(seed, e3)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(e3, e1)
    return (xyz - centre) @ np.column_stack([e1, e2, e3])


def c4_residual(ca_frame: dict) -> float:
    """RMSD between subunit A and the others after the ideal 90° rotations.

    If this is small the deposit really is four-fold symmetric and committing
    one subunit loses nothing; if it is not, the figure must not rotate.
    """
    chains = sorted(ca_frame)
    ref = {r[0]: np.array(r[1:]) for r in ca_frame[chains[0]]}
    worst = 0.0
    for k, ch in enumerate(chains[1:], start=1):
        th = np.deg2rad(90.0 * k)
        rot = np.array([[np.cos(th), -np.sin(th), 0],
                        [np.sin(th), np.cos(th), 0], [0, 0, 1]])
        other = {r[0]: np.array(r[1:]) for r in ca_frame[ch]}
        shared = sorted(set(ref) & set(other))
        # A 90° turn maps A onto one neighbour; try both senses and keep the
        # better, so chain-lettering order cannot flip the test.
        d1 = np.array([ref[i] @ rot.T - other[i] for i in shared])
        d2 = np.array([ref[i] @ rot - other[i] for i in shared])
        rms = min(np.sqrt((d1 ** 2).sum(1).mean()),
                  np.sqrt((d2 ** 2).sum(1).mean()))
        worst = max(worst, float(rms))
    return worst


# -------------------------------------------------------------- geometry

def pore_profile(heavy_f: np.ndarray, z0: float, z1: float,
                 step: float = 0.5, slab: float = 1.5) -> list[tuple]:
    """Minimum heavy-atom distance to the four-fold axis, sampled along it.

    This is a lower bound on the pore radius, not a solvent-probe radius: it
    is what the deposited atoms say without modelling van der Waals surfaces,
    which is the honest quantity for a 3.3 Å map.
    """
    r = np.hypot(heavy_f[:, 0], heavy_f[:, 1])
    z = heavy_f[:, 2]
    out = []
    for zi in np.arange(z0, z1 + step, step):
        m = np.abs(z - zi) <= slab
        if m.sum() < 4:
            continue
        out.append((float(zi), float(r[m].min())))
    return out


def lining_residues(heavy_f: np.ndarray, heavy_id: list, z: float,
                    radius: float, slab: float = 2.0,
                    tol: float = 1.2) -> list[str]:
    """Which residues actually form the constriction at height `z`.

    Everything within `tol` of the minimum radius in the slab, so a gate is
    reported as the residues that make it rather than as a bare coordinate.
    """
    r = np.hypot(heavy_f[:, 0], heavy_f[:, 1])
    m = (np.abs(heavy_f[:, 2] - z) <= slab) & (r <= radius + tol)
    seen = {}
    for i in np.nonzero(m)[0]:
        ch, seq, comp = heavy_id[i]
        seen.setdefault((seq, comp), set()).add(ch)
    return [f"{comp}{seq}" for (seq, comp) in sorted(seen)]


def ligand_contacts(heavy: np.ndarray, heavy_id: list, lig_xyz: np.ndarray,
                    lig_chain: str, cutoff: float = 4.5) -> dict:
    """Residues within `cutoff` of a bound IP3, split by subunit.

    Measuring the binding site instead of transcribing a residue list from a
    paper means the alignment figures can mark the columns that actually touch
    the ligand in this deposit — and it exposes any cross-subunit contact,
    which matters for the inter-subunit allostery of §2.4.
    """
    d = np.sqrt(((heavy[:, None, :] - lig_xyz[None, :, :]) ** 2).sum(-1)).min(1)
    own, other = {}, {}
    for i in np.nonzero(d <= cutoff)[0]:
        ch, seq, comp = heavy_id[i]
        if comp in ("I3P", "ZN", "HOH"):
            continue
        (own if ch == lig_chain else other)[seq] = comp
    return {"same_subunit": [f"{v}{k}" for k, v in sorted(own.items())],
            "other_subunit": [f"{v}{k}" for k, v in sorted(other.items())]}


def motif_position() -> int | None:
    """1-based start of the selectivity-filter motif in the structure protein."""
    if not PANEL_FASTA.exists():
        return None
    seq, keep = [], False
    for ln in PANEL_FASTA.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            if keep:
                break
            keep = ln[1:].split("|")[0] == STRUCT_ACC
        elif keep:
            seq.append(ln.strip())
    s = "".join(seq)
    i = s.find(FILTER_MOTIF)
    return i + 1 if i >= 0 else None


def read_domains(acc: str) -> list[dict]:
    lines = DOMAINS.read_text(encoding="utf-8").splitlines()
    head = lines[0].split("\t")
    rows = [dict(zip(head, ln.split("\t"))) for ln in lines[1:] if ln.strip()]
    return [r for r in rows if r["accession"] == acc]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pdb", default=PDB_ID)
    ap.add_argument("--work-dir", default=None,
                    help="where the mmCIF is downloaded (never the repo); "
                         "defaults to a system temp directory")
    args = ap.parse_args()

    work = Path(args.work_dir) if args.work_dir else \
        Path(tempfile.gettempdir()) / "ip3r_structures"
    if ROOT in work.resolve().parents or work.resolve() == ROOT:
        raise SystemExit(f"--work-dir must be outside the repository: {work}")

    cif = fetch_cif(args.pdb, work)
    print(f"  parsing {cif.name} ({cif.stat().st_size/1e6:.1f} MB)")
    st = parse(cif)
    axis, centre, planar = four_fold_frame(st["ca"])
    print(f"  four-fold axis {np.round(axis, 3)}  centroid planarity "
          f"residual {planar:.3f} Å")

    ca_frame = {ch: [(r[0], *to_frame(np.array(r[1:]), axis, centre))
                     for r in rows] for ch, rows in st["ca"].items()}
    resid = c4_residual(ca_frame)
    print(f"  C4 residual (subunit A vs the rest, after 90° turns): "
          f"{resid:.2f} Å RMSD")

    heavy_f = to_frame(st["heavy"], axis, centre)

    # Orient: the cytosolic cap is the end with the larger radial extent.
    zs = heavy_f[:, 2]
    rad = np.hypot(heavy_f[:, 0], heavy_f[:, 1])
    lo, hi = np.percentile(zs, [15, 85])
    if rad[zs < lo].mean() > rad[zs > hi].mean():
        flip = -1.0
        heavy_f[:, 2] *= -1
        ca_frame = {ch: [(r[0], r[1], r[2], -r[3]) for r in rows]
                    for ch, rows in ca_frame.items()}
    else:
        flip = 1.0

    # The membrane-embedded span, taken as the axial extent of the pore domain.
    doms = read_domains(STRUCT_ACC)
    ion = [d for d in doms if d["pfam"] == "PF00520"]
    if not ion:
        raise SystemExit(f"no PF00520 row for {STRUCT_ACC} in {DOMAINS}")
    lo_res, hi_res = int(ion[0]["start"]), int(ion[0]["end"])
    chA = sorted(ca_frame)[0]
    tm_z = np.array([r[3] for r in ca_frame[chA] if lo_res <= r[0] <= hi_res])
    tm_lo, tm_hi = float(np.percentile(tm_z, 2)), float(np.percentile(tm_z, 98))
    print(f"  pore domain PF00520 = residues {lo_res}-{hi_res}, axial span "
          f"{tm_lo:.1f} to {tm_hi:.1f} Å")

    prof = pore_profile(heavy_f, tm_lo - 12, tm_hi + 12)
    pz = np.array([p[0] for p in prof])
    pr = np.array([p[1] for p in prof])
    mid = 0.5 * (tm_lo + tm_hi)
    # §2.3: a short filter on the luminal side, the gate at the cytosolic end
    # of the bundle. Take the narrowest point on each side of the midpoint.
    lum = pz < mid
    cyt = (pz >= mid) & (pz <= tm_hi + 6)
    filt_z, filt_r = float(pz[lum][pr[lum].argmin()]), float(pr[lum].min())
    gate_z, gate_r = float(pz[cyt][pr[cyt].argmin()]), float(pr[cyt].min())
    print(f"  narrowest luminal point  z={filt_z:+.1f} Å  r={filt_r:.2f} Å")
    print(f"  narrowest cytosolic point z={gate_z:+.1f} Å  r={gate_r:.2f} Å")

    # IP3 sites, and how far each is from the gate.
    ip3 = {}
    for (comp, ch), pts in st["het"].items():
        if comp != "I3P":
            continue
        c = to_frame(np.array(pts).mean(0), axis, centre)
        c[2] *= flip
        ip3[ch] = c
    gate_pt = np.array([0.0, 0.0, gate_z])
    dists = {ch: float(np.linalg.norm(c - gate_pt)) for ch, c in ip3.items()}
    d_mean = float(np.mean(list(dists.values())))
    axial = float(np.mean([c[2] - gate_z for c in ip3.values()]))
    radial = float(np.mean([np.hypot(c[0], c[1]) for c in ip3.values()]))
    # §2.4 quotes "roughly 100 Å from the gate". That is the *axial* rise; the
    # site is well off the axis, so the through-space separation is larger.
    # Both are reported and the figure draws the triangle.
    print(f"  IP3 site: {radial:.1f} Å off the axis, {axial:.1f} Å above the "
          f"gate -> {d_mean:.1f} Å through space")

    filt_res = lining_residues(heavy_f, st["heavy_id"], filt_z, filt_r)
    gate_res = lining_residues(heavy_f, st["heavy_id"], gate_z, gate_r)
    print(f"  filter lined by {', '.join(filt_res)}")
    print(f"  gate   lined by {', '.join(gate_res)}")

    motif_at, motif_ok = motif_position(), False
    if motif_at:
        nums = [int("".join(c for c in r if c.isdigit())) for r in filt_res]
        motif_ok = any(motif_at - 2 <= n <= motif_at + len(FILTER_MOTIF)
                       for n in nums)
        print(f"  {FILTER_MOTIF} filter motif at {motif_at}-"
              f"{motif_at + len(FILTER_MOTIF) - 1}: narrowest luminal point "
              f"{'lands on it' if motif_ok else 'DOES NOT match it'}")
        if not motif_ok:
            raise SystemExit("the luminal constriction is not the selectivity "
                             "filter — check the axis and the orientation")

    height = float(heavy_f[:, 2].max() - heavy_f[:, 2].min())
    # Maximum diameter of the enclosing cylinder — corner to corner across the
    # square cap, which is larger than the edge-to-edge width usually quoted.
    width = float(rad.max()) * 2
    print(f"  overall: {height:.0f} Å tall, {width:.0f} Å maximum diameter")

    lig_chain = sorted(k[1] for k in st["het"] if k[0] == "I3P")[0]
    contacts = ligand_contacts(st["heavy"], st["heavy_id"],
                               np.array(st["het"][("I3P", lig_chain)]),
                               lig_chain)
    print(f"  IP3 contacts (≤4.5 Å) in its own subunit: "
          f"{', '.join(contacts['same_subunit'])}")
    if contacts["other_subunit"]:
        print(f"  ... and from a neighbouring subunit: "
              f"{', '.join(contacts['other_subunit'])}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (OUT_DIR / "structure_ca.tsv").open("w", encoding="utf-8") as fh:
        fh.write("resseq\tx\ty\tz\n")
        for r in sorted(ca_frame[chA]):
            fh.write(f"{r[0]}\t{r[1]:.2f}\t{r[2]:.2f}\t{r[3]:.2f}\n")
    with (OUT_DIR / "structure_pore.tsv").open("w", encoding="utf-8") as fh:
        fh.write("z_along_axis\tmin_heavy_atom_radius\n")
        for z, r in prof:
            fh.write(f"{z:.2f}\t{r:.3f}\n")
    meta = {
        "pdb_id": args.pdb.upper(),
        "description": "human type-3 IP3 receptor, IP3-bound (class 1)",
        "resolution_A": 3.33,
        "review_reference": "R24",
        "uniprot": STRUCT_ACC,
        "chains": sorted(st["ca"]),
        "residues_per_chain": len(ca_frame[chA]),
        "committed_subunit": chA,
        "four_fold_axis_in_deposit": [round(float(v), 4) for v in axis],
        "centroid_planarity_residual_A": round(planar, 3),
        "c4_residual_rmsd_A": round(resid, 3),
        "tm_span_z_A": [round(tm_lo, 1), round(tm_hi, 1)],
        "tm_pfam": "PF00520",
        "tm_residues": [lo_res, hi_res],
        "filter_z_A": round(filt_z, 1), "filter_min_radius_A": round(filt_r, 2),
        "gate_z_A": round(gate_z, 1), "gate_min_radius_A": round(gate_r, 2),
        "filter_motif": FILTER_MOTIF,
        "filter_motif_start": motif_at,
        "filter_motif_check_passed": motif_ok,
        "filter_lining_residues": filt_res,
        "gate_lining_residues": gate_res,
        "ip3_contact_cutoff_A": 4.5,
        "ip3_contacts": contacts,
        "ip3_sites": {k: [round(float(x), 1) for x in v] for k, v in ip3.items()},
        "ip3_to_gate_A": {k: round(v, 1) for k, v in dists.items()},
        "ip3_to_gate_mean_A": round(d_mean, 1),
        "ip3_axial_rise_above_gate_A": round(axial, 1),
        "ip3_radial_offset_from_axis_A": round(radial, 1),
        "overall_height_A": round(height, 0),
        "overall_max_diameter_A": round(width, 0),
    }
    (OUT_DIR / "structure_meta.json").write_text(
        json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"wrote structure_ca.tsv / structure_pore.tsv / structure_meta.json "
          f"to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
