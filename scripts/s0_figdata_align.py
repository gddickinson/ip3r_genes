#!/usr/bin/env python3
"""Alignments behind the review's sequence figures.

Two alignments, because the review makes two different sequence arguments and
they need different sequence sets.

**Family alignment** — the 25 committed control positives (ITPR1/2/3 across the
vertebrate panel plus the invertebrate / non-metazoan single-*itpr* grade).
Gives a per-column conservation profile mapped onto human IP₃R1 numbering, so
§2-§3 domain structure and the §9 disease positions can be read against
constraint.

**Superfamily alignment** — human IP₃R1/2/3 with human RyR1/2/3 and two
outgroup-grade invertebrates. This is what §7.1 asserts and this is what shows
it: at the pore the two families align residue for residue, and at the
IP₃-binding site — the one module RyR does not use — they do not.

The windows drawn in the figures are anchored on the *measured* sites from
`s0_figdata_structure.py` (IP₃ contacts, selectivity filter, gate), not on a
residue list copied from a paper.

    python scripts/s0_figdata_align.py [--threads 4]

Outputs (results/s0_baseline/review_figures/):
    align_conservation.tsv   per-position conservation on human IP3R1 numbering
    align_identity.tsv       all-pairs identity, family and superfamily sets
    align_blocks.tsv         the alignment windows the zoom figures print
    align_meta.json          MAFFT version + return codes, set membership
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = ROOT / "results" / "s0_baseline" / "review_figures"
POSITIVES = ROOT / "results" / "benchmark_controls" / "panel_positives.fasta"
DECOYS = ROOT / "results" / "benchmark_controls" / "panel_decoys.fasta"
META = OUT_DIR / "structure_meta.json"

REF = "Q14643"          # human ITPR1 — the numbering the profile is mapped to
STRUCT = "Q14573"       # human ITPR3 — the numbering the measured sites use

#: The superfamily set: both families, plus the invertebrate grade that shows
#: the ITPR side is not just three near-identical human paralogues.
SUPERFAMILY = ["Q14643", "Q14571", "Q14573",        # human ITPR1/2/3
               "P29993", "Q9Y0A1",                   # fly Itpr, worm itr-1
               "P21817", "Q92736", "Q15413"]         # human RYR1/2/3

AA = "ACDEFGHIKLMNPQRSTVWY"


# ------------------------------------------------------------------ io

def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """accession -> (header, sequence)."""
    out, key, buf, head = {}, None, [], ""
    for ln in path.read_text(encoding="utf-8").splitlines():
        if ln.startswith(">"):
            if key:
                out[key] = (head, "".join(buf))
            head = ln[1:]
            key = head.split("|")[0]
            buf = []
        else:
            buf.append(ln.strip())
    if key:
        out[key] = (head, "".join(buf))
    return out


def label_of(header: str) -> str:
    """`Q14643|ITPR1|Homo_sapiens|ITPR1` -> `ITPR1 H. sapiens`."""
    parts = header.split("|")
    sym = parts[1] if len(parts) > 1 else parts[0]
    sp = parts[2].replace("_", " ") if len(parts) > 2 else ""
    if sp:
        g, _, rest = sp.partition(" ")
        sp = f"{g[0]}. {rest}"
    return f"{sym} {sp}".strip()


def group_of(header: str) -> str:
    parts = header.split("|")
    return parts[3] if len(parts) > 3 else "?"


def run_mafft(seqs: dict[str, str], tag: str, threads: int) -> tuple[dict, dict]:
    """Align with MAFFT, and record that it really ran.

    S1 established that a silent MAFFT failure falls back to a star alignment
    without saying so, so the return code and version are captured and written
    to the meta file rather than assumed.
    """
    exe = shutil.which("mafft")
    if not exe:
        raise SystemExit("mafft not on PATH — required for the sequence figures")
    with tempfile.TemporaryDirectory() as td:
        fin = Path(td) / "in.fa"
        fin.write_text("".join(f">{k}\n{v}\n" for k, v in seqs.items()),
                       encoding="utf-8")
        cmd = [exe, "--auto", "--anysymbol", "--quiet",
               "--thread", str(threads), str(fin)]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(proc.stderr[-2000:], file=sys.stderr)
            raise SystemExit(f"mafft failed on {tag} (rc={proc.returncode})")
        aln, key, buf = {}, None, []
        for ln in proc.stdout.splitlines():
            if ln.startswith(">"):
                if key:
                    aln[key] = "".join(buf)
                key, buf = ln[1:].strip(), []
            else:
                buf.append(ln.strip())
        if key:
            aln[key] = "".join(buf)
    ver = subprocess.run([exe, "--version"], capture_output=True, text=True)
    trace = {"tag": tag, "returncode": proc.returncode, "n_seqs": len(aln),
             "n_columns": len(next(iter(aln.values()))),
             "version": (ver.stderr or ver.stdout).strip().split()[0],
             "command": " ".join(cmd[:-1]) + " <fasta>"}
    print(f"  mafft [{tag}] rc={proc.returncode} {trace['n_seqs']} seqs "
          f"x {trace['n_columns']} columns  ({trace['version']})")
    return aln, trace


# ----------------------------------------------------------- statistics

def column_conservation(col: str) -> tuple[float, float, str]:
    """(modal fraction, normalised 1 - Shannon entropy, modal residue).

    Gaps are excluded from the composition and reported separately, so a column
    that is 90% gap does not read as conserved on the strength of two residues.
    """
    res = [c for c in col if c not in "-."]
    if not res:
        return 0.0, 0.0, "-"
    counts = Counter(res)
    modal, n_modal = counts.most_common(1)[0]
    total = len(res)
    h = -sum((n / total) * math.log(n / total) for n in counts.values())
    hmax = math.log(min(len(AA), total)) or 1.0
    return n_modal / total, max(0.0, 1.0 - h / hmax), modal


def pairwise_identity(aln: dict[str, str], covered_only: bool) -> dict:
    """All-pairs percent identity.

    `covered_only` scores over mutually covered columns only — the
    fragment-aware metric S1 measures every margin under, kept here so the
    ITPR-vs-RyR numbers in the figures cannot be an artefact of length.
    """
    keys = list(aln)
    out = {}
    for i, a in enumerate(keys):
        for b in keys[i:]:
            sa, sb = aln[a], aln[b]
            match = n = 0
            for ca, cb in zip(sa, sb):
                ga, gb = ca in "-.", cb in "-."
                if ga and gb:
                    continue
                if covered_only and (ga or gb):
                    continue
                n += 1
                if ca == cb and not ga:
                    match += 1
            out[(a, b)] = match / n if n else 0.0
    return out


def map_columns(aln: dict[str, str], ref: str) -> list[int | None]:
    """Column index -> 1-based residue number in `ref` (None where ref gaps)."""
    pos, out = 0, []
    for c in aln[ref]:
        if c in "-.":
            out.append(None)
        else:
            pos += 1
            out.append(pos)
    return out


def ref_position_to_column(aln: dict[str, str], ref: str) -> dict[int, int]:
    return {p: i for i, p in enumerate(map_columns(aln, ref)) if p is not None}


# --------------------------------------------------------------- blocks

def measured_windows() -> list[dict]:
    """The windows the zoom figure prints, anchored on measured sites.

    Read from `structure_meta.json`; nothing here is a residue range typed in
    by hand, so re-running the structure script moves the figure with it.
    """
    if not META.exists():
        raise SystemExit(f"{META} missing — run s0_figdata_structure.py first")
    m = json.loads(META.read_text(encoding="utf-8"))
    nums = sorted(int("".join(c for c in r if c.isdigit()))
                  for r in m["ip3_contacts"]["same_subunit"])
    # Cluster the contacts into the runs they fall in, then pad each run.
    runs, cur = [], [nums[0]]
    for n in nums[1:]:
        if n - cur[-1] <= 12:
            cur.append(n)
        else:
            runs.append(cur)
            cur = [n]
    runs.append(cur)
    wins = [{"name": f"IP$_3$ contact {i+1}", "kind": "ligand",
             "start": r[0] - 4, "end": r[-1] + 4,
             "marks": r} for i, r in enumerate(runs)]
    f0 = m["filter_motif_start"]
    wins.append({"name": "selectivity filter", "kind": "pore",
                 "start": f0 - 5, "end": f0 + len(m["filter_motif"]) + 3,
                 "marks": list(range(f0, f0 + len(m["filter_motif"])))})
    g = sorted(int("".join(c for c in r if c.isdigit()))
               for r in m["gate_lining_residues"])
    wins.append({"name": "gate", "kind": "pore",
                 "start": g[0] - 5, "end": g[-1] + 5, "marks": g})
    return wins


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--threads", type=int, default=4)
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pos = read_fasta(POSITIVES)
    dec = read_fasta(DECOYS)
    traces = []

    # ---- family alignment -> conservation on human ITPR1 numbering
    fam_seqs = {k: v[1] for k, v in pos.items()}
    fam, tr = run_mafft(fam_seqs, "family", args.threads)
    traces.append(tr)
    colmap = map_columns(fam, REF)
    rows = []
    for i, refpos in enumerate(colmap):
        if refpos is None:
            continue
        col = "".join(fam[k][i] for k in fam)
        modal, cons, aa = column_conservation(col)
        gap = sum(1 for c in col if c in "-.") / len(col)
        rows.append((refpos, fam[REF][i], modal, cons, gap, aa))
    with (OUT_DIR / "align_conservation.tsv").open("w", encoding="utf-8") as fh:
        fh.write("itpr1_pos\titpr1_aa\tmodal_fraction\tconservation\t"
                 "gap_fraction\tmodal_aa\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]:.4f}\t{r[3]:.4f}\t{r[4]:.4f}\t"
                     f"{r[5]}\n")
    print(f"  conservation profile: {len(rows)} positions on {REF}")

    # ---- superfamily alignment -> identity + the zoom blocks
    pool = {**{k: v for k, v in pos.items()}, **{k: v for k, v in dec.items()}}
    missing = [a for a in SUPERFAMILY if a not in pool]
    if missing:
        raise SystemExit(f"not in the committed panels: {missing}")
    sup_seqs = {a: pool[a][1] for a in SUPERFAMILY}
    sup, tr = run_mafft(sup_seqs, "superfamily", args.threads)
    traces.append(tr)

    with (OUT_DIR / "align_identity.tsv").open("w", encoding="utf-8") as fh:
        fh.write("set\ta\tb\tlabel_a\tlabel_b\tgroup_a\tgroup_b\t"
                 "identity\tidentity_covered\n")
        for name, aln, src in (("family", fam, pos), ("superfamily", sup, pool)):
            full = pairwise_identity(aln, covered_only=False)
            cov = pairwise_identity(aln, covered_only=True)
            for (a, b), v in full.items():
                fh.write(f"{name}\t{a}\t{b}\t{label_of(src[a][0])}\t"
                         f"{label_of(src[b][0])}\t{group_of(src[a][0])}\t"
                         f"{group_of(src[b][0])}\t{v:.4f}\t{cov[(a, b)]:.4f}\n")
    print("  identity matrices written for both sets")

    p2c = ref_position_to_column(sup, STRUCT)
    # The conservation profile is drawn on ITPR1 numbering but the measured
    # sites come from an ITPR3 structure, so each window column also carries
    # the ITPR1 residue it aligns to — the two figures can then share an axis.
    c2ref = {i: p for i, p in enumerate(map_columns(sup, REF))}
    with (OUT_DIR / "align_blocks.tsv").open("w", encoding="utf-8") as fh:
        fh.write("window\tkind\taccession\tlabel\tgroup\titpr3_pos\t"
                 "itpr1_pos\tcolumn\tresidue\tis_marked\n")
        for w in measured_windows():
            cols = [(p, p2c[p]) for p in range(w["start"], w["end"] + 1)
                    if p in p2c]
            for acc in SUPERFAMILY:
                for p, c in cols:
                    fh.write(f"{w['name']}\t{w['kind']}\t{acc}\t"
                             f"{label_of(pool[acc][0])}\t{group_of(pool[acc][0])}\t"
                             f"{p}\t{c2ref[c] if c2ref[c] else ''}\t{c}\t"
                             f"{sup[acc][c]}\t{int(p in w['marks'])}\n")
    print(f"  {len(measured_windows())} measured windows written")

    (OUT_DIR / "align_meta.json").write_text(json.dumps({
        "reference_numbering": REF,
        "window_numbering": STRUCT,
        "family_set": sorted(fam_seqs),
        "superfamily_set": SUPERFAMILY,
        "mafft": traces,
        "windows": measured_windows(),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"wrote alignment tables to {OUT_DIR.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
