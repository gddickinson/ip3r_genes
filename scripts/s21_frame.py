"""The common coordinate frame, and the anchor test that can refuse it.

An exon boundary in the sweep's output is a position in the **bait's** residue
numbering, and the sweep used 38 baits. Comparing an ITPR1 gene's boundaries
with an ITPR2 gene's therefore needs one frame both reach, and that frame is
S6's committed alignment: all three human paralogues and human RYR2 are tips of
it, so a boundary travels bait → its cell's human reference → alignment column.

Two steps, each with a way to fail.

* **bait → human reference**, by pairwise MAFFT at `--thread 1` (D24: MAFFT is
  not reproducible at `--thread -1`, and every boundary in S21 is transferred
  through this map). A bait that *is* the reference maps to itself and no
  aligner runs. The identity and the fraction of the bait that lands on the
  reference are recorded per bait, because a boundary transferred through a
  poor pair is a weaker claim than one transferred through a good one.
* **human reference → alignment column**, by walking the reference's own
  aligned row. This is where the frame becomes shared, and it is checked rather
  than assumed: `anchor_test()` requires all 14 residues S0 *measured* on the
  6DQN structure — the ten IP3 contacts, the two selectivity-filter residues
  and the two gate residues, in each paralogue's own numbering from S17's
  committed `functional_sites.tsv` — to land in the **same column** in all
  three paralogues. They do, all 14. A frame that had slipped anywhere in the
  pore or the ligand core would fail there, and a slipped frame is exactly the
  error that makes a shared-intron count look like a result.

The unlabelled `vertebrate_basal` baits have no paralogue, so a locus won by
one has no reference of its own. Its boundaries are transferred through the
reference of the **cell** it filled and the row is flagged `frame_via_cell`;
S7 declined to place those tips, so the cross-paralogue tests are run with and
without them and both counts are reported.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as S17                                             # noqa: E402
import s21_lib as L                                               # noqa: E402

FUNCTIONAL_SITES = L.RESULTS / "constraint" / "functional_sites.tsv"

#: A bait whose pairwise alignment to its reference puts less than this
#: fraction of its own residues on the reference is refused a frame: its
#: boundaries would be transferred through gaps. Both families are within-family
#: comparisons here (61-68 % identical for ITPR, more for RyR), so this is a
#: floor against a broken alignment, not a similarity filter.
MIN_BAIT_ON_REF = 0.70


def res_to_col(aligned_row: str) -> dict[int, int]:
    """Residue number (1-based) -> alignment column (1-based)."""
    out, r = {}, 0
    for i, ch in enumerate(aligned_row):
        if ch != "-":
            r += 1
            out[r] = i + 1
    return out


def col_to_res(aligned_row: str) -> dict[int, int]:
    return {c: r for r, c in res_to_col(aligned_row).items()}


def cache_path() -> Path:
    d = L.data_root() / "s21"
    d.mkdir(parents=True, exist_ok=True)
    return d / "bait_frames.json"


def build(refresh: bool = False) -> dict:
    """`{bait_acc: {cell: {"map": {bait_res: column}, ...}}}`, cached.

    A bait is only ever asked about the cells it actually won, but the map is
    built per (bait, reference) pair and the reference is a property of the
    cell, so the cache is keyed that way. 38 baits and four references is at
    most 152 pairwise alignments and in practice 41.
    """
    path = cache_path()
    if path.exists() and not refresh:
        blob = json.loads(path.read_text())
        return {b: {c: {"map": {int(k): int(v) for k, v in d["map"].items()},
                        **{k: v for k, v in d.items() if k != "map"}}
                    for c, d in cells.items()}
                for b, cells in blob.items()}
    baits = L.bait_seqs()
    msa = L.msa_rows()
    labels = L.bait_labels()
    out: dict[str, dict] = {}
    for bait_acc, seq in sorted(baits.items()):
        meta = labels.get(bait_acc) or {}
        family = meta.get("family", "")
        cells = ["RYR"] if family == "RYR" else list(L.PARALOGS)
        for cell in cells:
            ref_acc = L.FRAME_ACC[cell]
            row = msa.get(ref_acc)
            if row is None:
                continue
            ref_seq = row.replace("-", "")
            r2c = res_to_col(row)
            if bait_acc == ref_acc:
                res_map = {i: r2c[i] for i in range(1, len(ref_seq) + 1)
                           if i in r2c}
                ident = 1.0
            else:
                b2r = S17.transfer_positions(seq, ref_seq)
                res_map = {b: r2c[r] for b, r in b2r.items() if r in r2c}
                ident = _pair_identity(seq, ref_seq, b2r)
            out.setdefault(bait_acc, {})[cell] = {
                "map": res_map, "ref": ref_acc,
                "bait_len": len(seq),
                "on_ref": round(len(res_map) / len(seq), 4) if seq else 0.0,
                "identity": round(ident, 4)}
    path.write_text(json.dumps(out, indent=0))
    return out


def _pair_identity(seq_a: str, seq_b: str, amap: dict[int, int]) -> float:
    """Identity over the residues the pair aligns (covered-only, S6's metric)."""
    if not amap:
        return 0.0
    same = sum(1 for a, b in amap.items()
               if seq_a[a - 1] == seq_b[b - 1])
    return same / len(amap)


def frame_rows(frames: dict) -> list[dict]:
    """One committed row per (bait, cell): what the frame is and what it cost."""
    labels = L.bait_labels()
    rows = []
    for bait, cells in sorted(frames.items()):
        meta = labels.get(bait) or {}
        for cell, d in sorted(cells.items()):
            rows.append({
                "bait": bait, "bait_paralog": meta.get("paralog", ""),
                "bait_clade": meta.get("clade", ""),
                "bait_family": meta.get("family", ""),
                "cell": cell, "reference": d["ref"], "bait_len": d["bait_len"],
                "mapped_residues": len(d["map"]),
                "frac_bait_on_ref": d["on_ref"], "pair_identity": d["identity"],
                "usable": int(d["on_ref"] >= MIN_BAIT_ON_REF),
                "frame_via_cell": int(meta.get("paralog", "") == ""
                                      and meta.get("family") == "ITPR")})
    return rows


# --------------------------------------------------------------------------
# the anchor test
# --------------------------------------------------------------------------
def anchor_test() -> list[dict]:
    """Every measured functional residue, in all three paralogues, one column.

    Read off S17's committed `functional_sites.tsv`, which carries S0's 6DQN
    measurements in each paralogue's own numbering with its own transfer
    verdict. Nothing here re-derives a coordinate; the test is only whether the
    alignment frame agrees with them.
    """
    msa = L.msa_rows()
    maps = {a: res_to_col(msa[a]) for a in
            (L.FRAME_ACC["ITPR1"], L.FRAME_ACC["ITPR2"], L.FRAME_ACC["ITPR3"])}
    by: dict[tuple[str, str], dict] = {}
    with open(FUNCTIONAL_SITES) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            by.setdefault((r["site_class"], r["source_resi"]), {})[
                r["paralog"]] = (r["ref_acc"], int(r["resi"]), r["aa"])
    rows = []
    for (site_class, src), d in sorted(by.items(), key=lambda kv: (kv[0][0],
                                                                  int(kv[0][1]))):
        cols = {}
        for paralog, (acc, resi, _aa) in d.items():
            cols[paralog] = maps.get(acc, {}).get(resi, 0)
        distinct = set(cols.values())
        rows.append({
            "site_class": site_class, "source_resi": int(src),
            "n_paralogs": len(d),
            **{f"resi_{p}": d[p][1] for p in sorted(d)},
            **{f"aa_{p}": d[p][2] for p in sorted(d)},
            **{f"col_{p}": cols[p] for p in sorted(cols)},
            "one_column": int(len(distinct) == 1 and 0 not in distinct)})
    return rows


def require_anchors() -> list[dict]:
    """The anchor test as a gate: a frame that fails it is not written."""
    rows = anchor_test()
    bad = [r for r in rows if not r["one_column"]]
    if bad or not rows:
        raise SystemExit(
            f"[s21] alignment frame refused: {len(bad)} of {len(rows)} measured "
            "functional residues do not share a column across the paralogues")
    return rows


def column_of(frames: dict, bait: str, cell: str, resi: int) -> int:
    """The alignment column a bait residue reaches, 0 when it reaches none."""
    d = (frames.get(bait) or {}).get(cell)
    if not d:
        return 0
    return d["map"].get(resi, 0)
