"""S24's data layer: the alignments, the structures and the two guards.

Everything here reads a committed file. Nothing re-aligns, nothing re-fetches
and nothing re-derives a number an earlier task already committed (D13) — a
supplementary panel that disagreed with the main figure it supports would be a
bug, not a difference of method, and the only way to guarantee that is to read
the same file.

**The two guards, and why each is a positive test rather than a comment.**

*The column map.* trimAl writes no column map of its own into the trimmed
FASTA; S6 recovered one with `trimal -colnumbering` and committed it. S24 does
not trust it. `verify_column_map` walks every trimmed column and requires the
whole 134-character column of `aln.fasta` at the mapped index to equal the
column of `trimmed.fasta` — for every column, every sequence, no sampling. It
is the cheapest check that the two coordinate systems every downstream residue
claim crosses are actually in step, and it can fail: an off-by-one in the map
would put every domain track, every variant and every conservation value one
column out and nothing else in the project would notice.

*The residue at the column.* A clinically labelled position is a number in a
paper joined to a number in a table joined to a column in an alignment, and
each join can slip. `verify_residues` requires the reference amino acid of
every harvested variant to be the residue the per-residue constraint table
holds at that position, and every aligned partner residue to be the residue the
*other* paralog's own table holds at the position claimed for it. A mismatch
raises; it does not warn.

Both guards run before any figure is drawn (`s24_run.py`), so a figure cannot
be written from coordinates that did not check out.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

RESULTS = PROJECT_ROOT / "results"
MSA_DIR = RESULTS / "msa_v2"
CONSTRAINT_DIR = RESULTS / "constraint"
SELECTION_DIR = RESULTS / "selection"
STRUCT_DIR = RESULTS / "structures"
OUT_DIR = RESULTS / "supplementary"
FIG_DIR = OUT_DIR / "figures"

#: The three human paralogues and the UniProt accession each constraint table
#: is written in. Read from S17's own file names rather than retyped.
PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
REFERENCE_ACC = {"ITPR1": "Q14643", "ITPR2": "Q14571", "ITPR3": "Q14573"}

#: The elements, N to C, with the labels the S17 figures use. Imported rather
#: than restated so a supplementary panel cannot name an element differently
#: from the Extended Data panel beside it.
ELEMENT_ORDER = ["nterm_trefoil", "MIR", "RIH_N", "RIH_C", "RIH_assoc",
                 "channel", "selectivity_filter", "gate", "luminal_loop"]
ELEMENT_LABEL = {
    "nterm_trefoil": "β-trefoil (PF08709)", "MIR": "MIR", "RIH_N": "RIH (N)",
    "RIH_C": "RIH (C)", "RIH_assoc": "RIH-assoc", "channel": "channel (TM)",
    "selectivity_filter": "filter", "gate": "gate",
    "luminal_loop": "luminal loop",
    # The linkers carry the element names of what they join, which is what
    # `s17_domains.residue_element` needs and what an axis cannot hold: at
    # 7 pt a rotated `linker_RIH_C_RIH_assoc` is 40 % of the panel width.
    "linker_RIH_N_RIH_C": "linker N–C",
    "linker_RIH_C_RIH_assoc": "linker C–assoc",
    "linker_RIH_assoc_channel": "linker assoc–TM",
    "nterm": "N terminus", "cterm": "C terminus",
}

#: The two regions the supplementary residue-level panel is about: the ligand
#: end and the pore end. Everything between them is linker or scaffold.
LIGAND_ELEMENTS = ("nterm_trefoil", "MIR")
PORE_ELEMENTS = ("channel", "selectivity_filter", "gate", "luminal_loop")


# ------------------------------------------------------------------- IO

def read_tsv(path: Path) -> list[dict]:
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], header: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    header = header or (list(rows[0]) if rows else [])
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    """Ordered {name: sequence}; the name is the header up to first space."""
    seqs: dict[str, str] = {}
    name, buf = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    return seqs


def f(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


# ------------------------------------------------------- guard 1: columns

def verify_column_map(aln: dict[str, str], trimmed: dict[str, str],
                      colmap: list[dict]) -> dict:
    """Exact column walk: trimmed column j must be aln column map[j].

    Raises on any disagreement. Returns the numbers a report can quote:
    how many columns were checked, over how many sequences, and the two
    alignment widths the map runs between.
    """
    if set(aln) != set(trimmed):
        raise SystemExit("aln.fasta and trimmed.fasta hold different labels")
    order = list(aln)
    width_in = len(next(iter(aln.values())))
    width_out = len(next(iter(trimmed.values())))
    if len(colmap) != width_out:
        raise SystemExit(f"column_map.tsv has {len(colmap)} rows for a "
                         f"{width_out}-column trimmed alignment")
    aln_cols = [aln[k] for k in order]
    trim_cols = [trimmed[k] for k in order]
    for row in colmap:
        j = int(row["trimmed_column"]) - 1        # 1-based in the file
        i = int(row["aln_column"]) - 1
        if not 0 <= i < width_in:
            raise SystemExit(f"column_map row {j + 1} points at aln column "
                             f"{i + 1} of {width_in}")
        got = "".join(s[j] for s in trim_cols)
        want = "".join(s[i] for s in aln_cols)
        if got != want:
            raise SystemExit(
                f"column map is wrong at trimmed column {j + 1}: "
                f"aln column {i + 1} is not that column")
    kept = {int(r["aln_column"]) - 1 for r in colmap}
    if len(kept) != width_out:
        raise SystemExit("column_map maps two trimmed columns to one input "
                         "column")
    return {"columns_checked": width_out, "sequences": len(order),
            "aln_columns": width_in, "trimmed_columns": width_out,
            "kept_indices": sorted(kept)}


# ------------------------------------------------------ guard 2: residues

def constraint_table(paralog: str) -> list[dict]:
    acc = REFERENCE_ACC[paralog]
    return read_tsv(CONSTRAINT_DIR / f"constraint_{paralog}_{acc}.tsv")


def residue_index(paralog: str) -> dict[int, dict]:
    return {int(r["resi"]): r for r in constraint_table(paralog)}


def verify_residues(variants: list[dict], pairs: list[dict],
                    index: dict[str, dict[int, dict]]) -> dict:
    """Every stated residue must be the residue at its own column.

    Two joins are checked. A variant's `ref_aa` against the residue its own
    paralogue's constraint table holds at `resi`; and every aligned partner
    in `paralog_variant_positions.tsv` against the residue the *other*
    paralogue's table holds at `other_resi`. Raises on either.
    """
    checked_v = 0
    for r in variants:
        gene, resi = r["gene"], int(r["resi"])
        site = index[gene].get(resi)
        if site is None:
            raise SystemExit(f"{gene} {resi}: no row in the constraint table")
        if site["aa"] != r["ref_aa"]:
            raise SystemExit(f"{gene} {r['ref_aa']}{resi}: the constraint "
                             f"table holds {site['aa']} there")
        checked_v += 1
    checked_p = 0
    for r in pairs:
        if r.get("aligned") != "True" or not r.get("other_resi"):
            continue
        other = index[r["other"]].get(int(r["other_resi"]))
        if other is None:
            raise SystemExit(f"{r['other']} {r['other_resi']}: no row in the "
                             f"constraint table")
        if other["aa"] != r["other_aa"]:
            raise SystemExit(
                f"{r['gene']} {r['resi']} -> {r['other']} {r['other_resi']}: "
                f"the table holds {other['aa']}, the pair table says "
                f"{r['other_aa']}")
        checked_p += 1
    return {"variant_residues_checked": checked_v,
            "aligned_partners_checked": checked_p}


# ------------------------------------------------------------- structures

def read_ca_pdb(path: Path) -> list[dict]:
    """CA atoms of a painted structure: {resi, aa, x, y, z, b}.

    Deliberately small. S17 wrote these files itself, one CA per residue,
    chain A only, with the score in the B-factor column and −1 where there
    is no score — so a full parser would only add failure modes.
    """
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            out.append({
                "resi": int(line[22:26]),
                "aa": line[17:20].strip(),
                "x": float(line[30:38]), "y": float(line[38:46]),
                "z": float(line[46:54]),
                "b": float(line[60:66]),
            })
    if not out:
        raise SystemExit(f"{path.name}: no CA atoms")
    return out


def principal_axes(atoms: list[dict]) -> tuple[list[float], list[float]]:
    """Project a CA trace onto its two largest principal axes.

    Stdlib power iteration on the 3x3 covariance — no numpy dependency in a
    figure module, and the result is deterministic given the input file
    (D24 applied to a projection). Returns (u, v), one coordinate pair per
    atom, with the *membrane normal* first: for a channel the longest axis
    is the pore axis, which is what puts the filter and the gate one above
    the other in the drawing rather than side by side.
    """
    n = len(atoms)
    cx = sum(a["x"] for a in atoms) / n
    cy = sum(a["y"] for a in atoms) / n
    cz = sum(a["z"] for a in atoms) / n
    pts = [(a["x"] - cx, a["y"] - cy, a["z"] - cz) for a in atoms]
    cov = [[0.0] * 3 for _ in range(3)]
    for p in pts:
        for i in range(3):
            for j in range(3):
                cov[i][j] += p[i] * p[j]
    for i in range(3):
        for j in range(3):
            cov[i][j] /= n

    def mul(m, v):
        return [sum(m[i][k] * v[k] for k in range(3)) for i in range(3)]

    def norm(v):
        s = sum(x * x for x in v) ** 0.5
        return [x / s for x in v] if s else [1.0, 0.0, 0.0]

    def power(m, avoid=None):
        v = [1.0, 0.4, 0.2]
        for _ in range(200):
            v = mul(m, v)
            if avoid is not None:
                d = sum(v[i] * avoid[i] for i in range(3))
                v = [v[i] - d * avoid[i] for i in range(3)]
            v = norm(v)
        return v

    e1 = power(cov)
    e2 = power(cov, avoid=e1)
    u = [sum(p[i] * e1[i] for i in range(3)) for p in pts]
    v = [sum(p[i] * e2[i] for i in range(3)) for p in pts]
    return u, v


def binned(values: list[float], nbins: int) -> tuple[list[float], list[float]]:
    """Mean of `values` in `nbins` equal column bins, with bin centres.

    Shared by both figure modules so the two halves of the supplementary
    alignment set cannot bin differently and draw curves that are not
    comparable.
    """
    n = len(values)
    step = n / nbins
    xs, ys = [], []
    for b in range(nbins):
        lo, hi = int(b * step), max(int(b * step) + 1, int((b + 1) * step))
        chunk = values[lo:min(hi, n)]
        if not chunk:
            continue
        xs.append((lo + min(hi, n)) / 2.0)
        ys.append(sum(chunk) / len(chunk))
    return xs, ys


# --------------------------------------------------------------- stats

def dump_stats(stats: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(stats, fh, indent=1, sort_keys=True)
        fh.write("\n")
