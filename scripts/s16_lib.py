"""S16 shared helpers — duplication history: loci, copies, Ensembl.

S16 asks where ITPR1/2/3 came from and what 3R did to them. Everything it
needs is already on disk except one thing:

* `<data_root>/genome_sweep/<acc>/summary.json` — **every modelled locus**
  per genome, not just the one the ledger kept. The ledger holds one row per
  genome x cell and therefore only each cell's best locus, which is exactly
  the information a duplication question needs (S8 and S15a both hit this).
* `<data_root>/genome_sweep/<acc>/genes_slim.tsv` — the annotation's gene
  table, read through `s8_flank_lib` so S16's windows are S8's windows.
* `results/loss_dynamics/integrity_loci.tsv` — S15a's calibrated lesion bar,
  read back rather than re-derived, so a copy that is a corpse can be told
  from a copy that is a gene.
* Ensembl BioMart — the human paralogy map the 2R test runs on. Fetched by
  `s16_ensembl_map.py` and cached permanently under
  `<data_root>/raw_api/s16/`; every later run is offline.

IO and the copy rule only. The analyses live in `s16_copy_number.py`,
`s16_paralogon.py`, `s16_quartet.py`, `s16_teleost.py`.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.utils.data_root import get_data_root                  # noqa: E402
import s5_calibration as s5cal                                 # noqa: E402

RESULTS = PROJECT_ROOT / "results"
LEDGER = RESULTS / "genome_ledger" / "genome_ledger.tsv"
MANIFEST = RESULTS / "genome_manifest.tsv"
INTEGRITY = RESULTS / "loss_dynamics" / "integrity_loci.tsv"
OUT_DIR = RESULTS / "duplication"

#: The three vertebrate paralog cells. `RYR` is the sister family's cell and
#: is *not* a paralog here — it is this task's positive control, because the
#: ryanodine receptors are a three-member vertebrate family of the same age
#: and the same 2R candidacy. D14's hazard, used as an instrument.
PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
CONTROL_CELL = "RYR"

#: A locus counts as a copy when its own model covers this much of its bait.
#: Not a tuning knob: two loci each covering a complete bait are two genes,
#: whereas two covering complementary halves are one gene the aligner split,
#: which `merge_split_models` folds before the bar is ever applied. The
#: sensitivity of every count to this bar is measured in `copy_sensitivity`.
COV_FULL = 0.50
#: Below this the alignment is a scrap, whatever its coverage of a short bait.
MIN_COPY_ALIGNED_AA = 500
#: The sweep already refused to record a cluster below 0.40 identity
#: (`min_locus_identity` in every summary.json), so this floor is the sweep's
#: own, read back rather than retyped, and is here to be *stated*.
MIN_COPY_IDENTITY = 0.40

#: Two alignments of one gene that the aligner emitted separately. Merged
#: when they sit on one strand within this gap and their bait spans are
#: complementary rather than repeated (`s23_copy_number.merge_split_loci`'s
#: rule, applied to the S5 sweep's own records).
MERGE_GAP_BP = 100_000
MERGE_MAX_QUERY_OVERLAP = 0.20


def data_root() -> Path:
    return get_data_root()


def sweep_dir(acc: str = "") -> Path:
    d = data_root() / "genome_sweep"
    return d / acc if acc else d


def cache_dir() -> Path:
    d = data_root() / "raw_api" / "s16"
    d.mkdir(parents=True, exist_ok=True)
    return d


def out_dir() -> Path:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    return OUT_DIR


# --------------------------------------------------------------- table IO

def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, header: list[str], rows: list[list]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fnum(v, cast=float, default=0):
    try:
        return cast(v)
    except (TypeError, ValueError):
        return default


def fmt(v) -> str:
    """Floats at %.6g — a p-value written to four decimal places is 0.0000,
    which hides how strong a claim is rather than how weak (S15a's rule)."""
    if isinstance(v, float):
        return f"{v:.6g}"
    return "" if v is None else str(v)


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """BH-adjusted q-values, stdlib. S9's rule: this task asks the same
    paralogy question of many pairs, windows and level sets, and the
    smallest of many uncorrected p-values is the error the self-tests
    exist to avoid."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank, i in enumerate(reversed(order), 1):
        k = n - rank + 1
        prev = min(prev, pvals[i] * n / k)
        q[i] = min(1.0, prev)
    return q


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    p = RESULTS / "session_live.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        "task": "S16", "stage": stage,
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "steps": [{"label": l, "done": bool(d)} for l, d in steps],
    }, indent=1))


# --------------------------------------------------------------- genomes

_GIDX: dict[str, dict] | None = None


def genome_index() -> dict[str, dict]:
    """accession -> genome-level columns, plus D4's contiguity flag.

    `contig_spans_gene` is `s5_calibration.spans_a_gene` and not a threshold
    of this task's own: a copy count in an assembly whose contigs cannot
    carry the gene is a different measurement from one above the bar, and
    the two tasks must not disagree about where the bar is.
    """
    global _GIDX
    if _GIDX is not None:
        return _GIDX
    idx: dict[str, dict] = {}
    for r in read_tsv(LEDGER):
        acc = r["accession"]
        if acc in idx:
            continue
        n50 = fnum(r.get("contig_n50"), int, 0)
        idx[acc] = {
            "accession": acc,
            "organism": r["organism"],
            "vclass": r["vclass"],
            "vorder": r["vorder"],
            "annotated": r.get("annotated_assembly", ""),
            "assembly_level": r.get("assembly_level", ""),
            "contig_n50": n50,
            "genome_bp": fnum(r.get("genome_bp"), int, 0),
            "contig_spans_gene": bool(s5cal.spans_a_gene(n50)),
        }
    _GIDX = idx
    return idx


def lesion_bar() -> tuple[float, str]:
    """S15a's committed lesion-density bar, recovered from its own table.

    S15a scored each locus `intact` / `lesion_rich` against a bar measured
    on full-coverage, contiguous, family-named loci. Recomputing the bar
    here would let two tasks drift; instead the bar is read back off the
    verdicts S15a committed — the largest density it still called intact.
    """
    if not INTEGRITY.exists():
        return 0.0, "S15a integrity table absent"
    intact = [fnum(r.get("lesion_density"), float, 0.0)
              for r in read_tsv(INTEGRITY) if r.get("verdict") == "intact"]
    if not intact:
        return 0.0, "S15a table has no intact verdict"
    return max(intact), (f"largest density S15a still called intact, over "
                         f"{len(intact)} loci")


# ----------------------------------------------------------- per-locus IO

def _locus_record(acc: str, g: dict, cell: str, cell_status: str,
                  slot: str, idx: int, loc: dict) -> dict:
    aligned = fnum(loc.get("aligned_aa"), int, 0)
    if not aligned:
        span = loc.get("q_span") or []
        if len(span) == 2:
            aligned = max(0, fnum(span[1], int, 0) - fnum(span[0], int, 0) + 1)
    fs = fnum(loc.get("frameshifts"), int, 0)
    st = fnum(loc.get("stop_codons"), int, 0)
    lesions = fs + st
    annot = loc.get("annot_gene") or {}
    qs = loc.get("q_span") or [0, 0]
    return {
        "accession": acc,
        "organism": g.get("organism", ""),
        "vclass": g.get("vclass", ""),
        "vorder": g.get("vorder", ""),
        "assembly_level": g.get("assembly_level", ""),
        "contig_n50": g.get("contig_n50", 0),
        "contig_spans_gene": g.get("contig_spans_gene", False),
        "cell": cell,
        "cell_status": cell_status,
        "slot": slot,
        "locus_idx": idx,
        "locus_clade": loc.get("locus_top_clade", ""),
        "locus_family": loc.get("locus_family", ""),
        "assigned_via": loc.get("assigned_via", ""),
        "contig": loc.get("contig", ""),
        "start": fnum(loc.get("start"), int, 0),
        "end": fnum(loc.get("end"), int, 0),
        "strand": loc.get("strand", ""),
        "bait": loc.get("bait", ""),
        "identity": fnum(loc.get("identity"), float, 0.0),
        "coverage": fnum(loc.get("coverage"), float, 0.0),
        "aligned_aa": aligned,
        "q_start": fnum(qs[0], int, 0),
        "q_end": fnum(qs[1], int, 0),
        "bait_len": fnum(loc.get("bait_len"), int, 0),
        "frameshifts": fs,
        "stop_codons": st,
        "lesions": lesions,
        "lesion_density": round(1000.0 * lesions / aligned, 4) if aligned else 0.0,
        "family_margin": fnum(loc.get("family_margin"), float, 0.0),
        "paralog_margin": fnum(loc.get("paralog_margin"), float, 0.0),
        "contig_edge": bool(loc.get("contig_edge")),
        "longest_n_run": fnum(loc.get("longest_n_run"), int, 0),
        "annot_gene": (annot.get("name") if isinstance(annot, dict) else "") or "",
        "annot_frac_cds": (fnum(annot.get("frac_cds"), float, 0.0)
                           if isinstance(annot, dict) else 0.0),
        "annot_paralog_matches": bool(loc.get("annot_paralog_matches")),
        "merged_models": 1,
        "merged_from": "",
    }


def iter_loci(accessions: list[str] | None = None,
              cells: tuple[str, ...] = PARALOGS + (CONTROL_CELL,)
              ) -> list[dict]:
    """Every modelled locus in every swept genome, one row each.

    `slot` records where it came from: `cell` = inside a bait-class cell (the
    sweep's own loci, including the secondary ones the ledger drops), `other`
    = the sweep's `other_loci`, alignments that joined no cell.
    """
    gidx = genome_index()
    accs = accessions if accessions is not None else sorted(gidx)
    out: list[dict] = []
    for acc in accs:
        path = sweep_dir(acc) / "summary.json"
        if not path.exists():
            continue
        with open(path) as fh:
            summ = json.load(fh)
        g = gidx.get(acc, {})
        for cell, c in (summ.get("cells") or {}).items():
            if cell not in cells:
                continue
            for i, loc in enumerate(c.get("loci") or []):
                out.append(_locus_record(acc, g, cell, c.get("status", ""),
                                         "cell", i, loc))
        for i, loc in enumerate(summ.get("other_loci") or []):
            cell = loc.get("locus_top_clade", "") or "unassigned"
            out.append(_locus_record(acc, g, cell, "other", "other", i, loc))
    return out


# ------------------------------------------------------------- the merge

def _query_overlap_frac(a: dict, b: dict) -> float:
    """How much of the shorter bait span the two alignments share."""
    lo = max(a["q_start"], b["q_start"])
    hi = min(a["q_end"], b["q_end"])
    ov = max(0, hi - lo + 1)
    shortest = min(a["q_end"] - a["q_start"] + 1, b["q_end"] - b["q_start"] + 1)
    return ov / shortest if shortest > 0 else 0.0


def merge_split_models(recs: list[dict], gap: int = MERGE_GAP_BP,
                       max_overlap: float = MERGE_MAX_QUERY_OVERLAP
                       ) -> tuple[list[dict], list[dict]]:
    """Fold neighbouring same-cell models that look like one split gene.

    This is the brief's requirement and it is load-bearing in exactly one
    direction: a split model manufactures a duplication, which is the claim
    S16 is testing. Two models are one gene when they sit on one strand of
    one contig within `gap`, **and** their bait spans are complementary
    rather than repeated — two alignments each covering residues 1-2748 of
    the same bait are two genes however close they are, and two covering
    1-1300 and 1310-2748 are one gene however far apart.

    Returns (merged records, merge records). Every merge writes down the two
    numbers it was made on, so a copy number can be checked rather than
    trusted.
    """
    by_key: dict[tuple, list[dict]] = {}
    for r in recs:
        by_key.setdefault((r["accession"], r["cell"], r["contig"],
                           r["strand"]), []).append(r)
    kept: list[dict] = []
    merges: list[dict] = []
    for key, group in sorted(by_key.items()):
        group.sort(key=lambda r: (r["start"], r["end"]))
        cur: dict | None = None
        for r in group:
            if cur is None:
                cur = dict(r)
                continue
            dist = r["start"] - cur["end"]
            ov = _query_overlap_frac(cur, r)
            if dist <= gap and ov <= max_overlap:
                merges.append({
                    "accession": key[0], "cell": key[1], "contig": key[2],
                    "strand": key[3],
                    "kept_start": cur["start"], "kept_end": cur["end"],
                    "merged_start": r["start"], "merged_end": r["end"],
                    "gap_bp": dist, "query_overlap": round(ov, 4),
                    "kept_q": f"{cur['q_start']}-{cur['q_end']}",
                    "merged_q": f"{r['q_start']}-{r['q_end']}",
                    "why": (f"{dist:,} bp apart on one strand and their bait "
                            f"spans overlap {ov:.0%} <= "
                            f"{max_overlap:.0%} — one gene, two alignments"),
                })
                best = cur if (cur["coverage"], cur["aligned_aa"]) >= \
                    (r["coverage"], r["aligned_aa"]) else r
                new = dict(best)
                new["start"] = min(cur["start"], r["start"])
                new["end"] = max(cur["end"], r["end"])
                new["q_start"] = min(cur["q_start"], r["q_start"])
                new["q_end"] = max(cur["q_end"], r["q_end"])
                new["aligned_aa"] = cur["aligned_aa"] + r["aligned_aa"]
                new["coverage"] = min(1.0, cur["coverage"] + r["coverage"])
                new["lesions"] = cur["lesions"] + r["lesions"]
                new["frameshifts"] = cur["frameshifts"] + r["frameshifts"]
                new["stop_codons"] = cur["stop_codons"] + r["stop_codons"]
                new["lesion_density"] = (
                    round(1000.0 * new["lesions"] / new["aligned_aa"], 4)
                    if new["aligned_aa"] else 0.0)
                new["merged_models"] = cur["merged_models"] + r["merged_models"]
                new["merged_from"] = ";".join(
                    x for x in (cur["merged_from"],
                                f"{r['contig']}:{r['start']}-{r['end']}") if x)
                cur = new
            else:
                kept.append(cur)
                cur = dict(r)
        if cur is not None:
            kept.append(cur)
    kept.sort(key=lambda r: (r["accession"], r["cell"], r["contig"],
                             r["start"]))
    return kept, merges


# ---------------------------------------------------------- the copy rule

def copy_exclusion(rec: dict, cov_full: float = COV_FULL) -> str:
    """Why this locus is not counted as a gene copy — "" when it is one."""
    if rec["cell"] not in PARALOGS + (CONTROL_CELL,):
        return f"cell_{rec['cell'] or 'unassigned'}"
    if rec["identity"] < MIN_COPY_IDENTITY:
        return "identity_below_sweep_floor"
    if rec["coverage"] < cov_full:
        return "coverage_below_full"
    if rec["aligned_aa"] < MIN_COPY_ALIGNED_AA:
        return "too_short"
    return ""


def is_copy(rec: dict, cov_full: float = COV_FULL) -> bool:
    return copy_exclusion(rec, cov_full) == ""


def integrity_call(rec: dict, bar: float) -> str:
    """intact / lesion_rich, on S15a's own bar. One-sided by decision:
    S10 established that zero lesions falsifies a pseudogene call and a
    handful does not establish one, so `lesion_rich` is a flag, not a
    verdict of death."""
    if rec["lesions"] == 0:
        return "intact"
    return "lesion_rich" if rec["lesion_density"] > bar else "intact"
