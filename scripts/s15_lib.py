"""s15_lib.py — S15's data layer: paths, IO, the sweep's own per-genome
evidence, and the two statistics the loss argument is read through.

Nothing here decides a state; `s15_states.py` does that.  What lives here
is the plumbing every S15 module shares, plus three things worth naming.

`load_summaries()` reads the sweep's archived `summary.json` per genome
rather than the ledger.  The ledger holds one row per genome x cell and
therefore only the *best* locus, and S15 has to reason about every locus a
genome carries — the spare loci are what tell `paralog_unassignable` from
an absence (D45).

`known_loci()` reproduces `s5_run_sweep.filter_hsps_outside`'s exclusion
set: every locus the aligner clustered anywhere in the genome, any bait,
padded.  Every reconstruction in S15 is computed on HSPs *outside* that
set, so a reference cannot be reassembled out of the genome's other
paralogs' genes.  This is the cross-paralog control, and it is inherited
from the sweep rather than re-invented.

`spearman()` and `sign_test()` are stdlib because S15's confounder
controls are a rank correlation and a paired test, and adding a scipy
dependency to a task whose whole point is a negative result would put a
third-party version number inside a claim.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from pathlib import Path

PROJECT = Path(__file__).resolve().parent.parent
RESULTS = PROJECT / "results"
LEDGER = RESULTS / "genome_ledger" / "genome_ledger.tsv"
RESCUE = RESULTS / "genome_ledger" / "rescue_regions.tsv"
MANIFEST = RESULTS / "genome_manifest.tsv"
OUT = RESULTS / "loss_dynamics"
FIGS = OUT / "figures"

ITPR_CELLS = ("ITPR1", "ITPR2", "ITPR3")
CONTROL_CELL = "RYR"
ALL_CELLS = ITPR_CELLS + (CONTROL_CELL,)

#: the sweep's own known-locus pad, from s5_run_sweep.filter_hsps_outside
KNOWN_PAD = 5000
#: the sweep's own rescue significance cut, from s5_rescue.RESCUE_E
RESCUE_E = 1e-5


def data_root() -> Path:
    import sys
    sys.path.insert(0, str(PROJECT))
    from src.utils.data_root import require_data_root
    return require_data_root()


def sweep_dir() -> Path:
    return data_root() / "genome_sweep"


# ------------------------------------------------------------------- TSV IO

def read_tsv(path: Path) -> list[dict]:
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if not lines:
        return []
    cols = lines[0].split("\t")
    return [dict(zip(cols, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def write_tsv(path: Path, rows: list[dict], cols: list[str]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c, "")) for c in cols) + "\n")
    return path


def _fmt(v) -> str:
    """Six significant figures, not four decimal places.

    A fixed 4-dp format writes a p-value of 2.1e-07 as `0.0000`, which
    the report then reads back as zero and renders as "p < 1e-300".  It
    did exactly that on this task's first build, so floats are written at
    `%.6g` — faithful across the eleven orders of magnitude the p-values
    and the densities span, and still short.
    """
    if v is None:
        return ""
    if isinstance(v, bool):
        return "1" if v else "0"
    if isinstance(v, float):
        if v != v or v in (float("inf"), float("-inf")):
            return "nan"
        return f"{v:.6g}"
    return str(v).replace("\t", " ").replace("\n", " ")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    p = RESULTS / "session_live.json"
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({
        "task": "S15a", "stage": stage,
        "updated": time.strftime("%Y-%m-%d %H:%M:%S"),
        "steps": [{"label": l, "done": bool(d)} for l, d in steps],
    }, indent=1))


# ------------------------------------------------------- the sweep's output

def load_summaries(accessions: list[str] | None = None) -> dict[str, dict]:
    """The archived per-genome summary.json, keyed by accession."""
    root = sweep_dir()
    out: dict[str, dict] = {}
    accs = accessions if accessions is not None else sorted(
        p.name for p in root.iterdir() if p.is_dir())
    for acc in accs:
        p = root / acc / "summary.json"
        if p.exists():
            out[acc] = json.loads(p.read_text())
    return out


def known_loci(summary: dict) -> list[tuple[str, int, int]]:
    """Every locus the aligner clustered in this genome, any bait.

    The same set `s5_run_sweep.filter_hsps_outside` excluded before the
    rescue HSPs were clustered — so a reconstruction built on the filtered
    HSPs cannot reuse the genome's other family genes.
    """
    out = []
    for cell in (summary.get("cells") or {}).values():
        for loc in cell.get("loci", []) or []:
            out.append((loc["contig"], int(loc["start"]), int(loc["end"])))
    for loc in summary.get("other_loci") or []:
        if loc.get("contig"):
            out.append((loc["contig"], int(loc["start"]), int(loc["end"])))
    return out


def outside_known(contig: str, lo: int, hi: int,
                  known: list[tuple[str, int, int]],
                  pad: int = KNOWN_PAD) -> bool:
    return not any(c == contig and lo <= e + pad and hi >= s - pad
                   for c, s, e in known)


def bait_lengths() -> dict[str, int]:
    """Reference lengths from the committed bait panel, keyed by header id."""
    path = RESULTS / "s5_baits" / "baits.faa"
    if not path.exists():                       # the sweep's archived copy
        path = sweep_dir() / "baits.faa"
    out, name, seq = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = len("".join(seq))
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip())
    if name:
        out[name] = len("".join(seq))
    return out


def bait_clade(header: str) -> str:
    """The bait's clade band / paralog label — the last '|' field."""
    return header.split("|")[-1]


# ------------------------------------------------------------- geometry

def union_length(intervals: list[tuple[int, int]]) -> int:
    """Total length covered by a set of closed 1-based intervals."""
    if not intervals:
        return 0
    total, cs, ce = 0, None, None
    for a, b in sorted((min(x, y), max(x, y)) for x, y in intervals):
        if cs is None:
            cs, ce = a, b
        elif a <= ce + 1:
            ce = max(ce, b)
        else:
            total += ce - cs + 1
            cs, ce = a, b
    if cs is not None:
        total += ce - cs + 1
    return total


def overlaps(a_lo: int, a_hi: int, b_lo: int, b_hi: int) -> bool:
    return a_lo <= b_hi and a_hi >= b_lo


# ------------------------------------------------------------- statistics

def _ranks(xs: list[float]) -> list[float]:
    order = sorted(range(len(xs)), key=lambda i: xs[i])
    ranks = [0.0] * len(xs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and xs[order[j + 1]] == xs[order[i]]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> tuple[float, float, int]:
    """Spearman rho with a t-approximation p-value. Returns (rho, p, n)."""
    pairs = [(x, y) for x, y in zip(xs, ys)
             if x is not None and y is not None]
    n = len(pairs)
    if n < 4:
        return (float("nan"), float("nan"), n)
    rx = _ranks([p[0] for p in pairs])
    ry = _ranks([p[1] for p in pairs])
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx)
                    * sum((b - my) ** 2 for b in ry))
    if den == 0:
        return (float("nan"), float("nan"), n)
    rho = num / den
    if abs(rho) >= 1.0:
        return (rho, 0.0, n)
    t = rho * math.sqrt((n - 2) / (1 - rho * rho))
    return (rho, _t_sf(abs(t), n - 2) * 2, n)


def _t_sf(t: float, df: int) -> float:
    """Upper tail of Student's t, via the incomplete beta function."""
    x = df / (df + t * t)
    return 0.5 * _betainc(df / 2.0, 0.5, x)


def _betainc(a: float, b: float, x: float) -> float:
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lbeta = (math.lgamma(a) + math.lgamma(b) - math.lgamma(a + b))
    front = math.exp(math.log(x) * a + math.log(1 - x) * b - lbeta) / a
    f, c, d = 1.0, 1.0, 0.0
    for i in range(0, 300):
        m = i // 2
        if i == 0:
            num = 1.0
        elif i % 2 == 0:
            num = (m * (b - m) * x) / ((a + 2 * m - 1) * (a + 2 * m))
        else:
            num = -((a + m) * (a + b + m) * x) / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        d = 1e-30 if abs(d) < 1e-30 else d
        d = 1.0 / d
        c = 1.0 + num / c
        c = 1e-30 if abs(c) < 1e-30 else c
        f *= c * d
        if abs(1.0 - c * d) < 1e-12:
            break
    return front * (f - 1.0)


def sign_test(diffs: list[float]) -> dict:
    """Two-sided exact sign test; ties dropped and counted (S8's rule)."""
    pos = sum(1 for d in diffs if d > 0)
    neg = sum(1 for d in diffs if d < 0)
    ties = sum(1 for d in diffs if d == 0)
    n = pos + neg
    if n == 0:
        return dict(n=0, n_pos=0, n_neg=0, n_ties=ties, p=float("nan"))
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(0, k + 1)) / (2.0 ** n)
    return dict(n=n, n_pos=pos, n_neg=neg, n_ties=ties,
                p=min(1.0, 2.0 * tail))


def median(xs: list[float]) -> float:
    s = sorted(x for x in xs if x is not None)
    if not s:
        return float("nan")
    m = len(s) // 2
    return s[m] if len(s) % 2 else (s[m - 1] + s[m]) / 2.0


def youden(pos: list[float], neg: list[float]) -> dict:
    """The threshold on a score that best separates two populations.

    Reported with its J statistic so a pair of populations that does *not*
    separate is visible as such — `s23_calibrate_loci.separation()`'s rule:
    a threshold is only worth quoting beside the separation it achieves.
    """
    if not pos or not neg:
        return dict(threshold=float("nan"), j=float("nan"),
                    sens=float("nan"), spec=float("nan"),
                    n_pos=len(pos), n_neg=len(neg))
    best = None
    for t in sorted(set(pos) | set(neg)):
        sens = sum(1 for x in pos if x >= t) / len(pos)
        spec = sum(1 for x in neg if x < t) / len(neg)
        j = sens + spec - 1.0
        if best is None or j > best[1]:
            best = (t, j, sens, spec)
    return dict(threshold=best[0], j=best[1], sens=best[2], spec=best[3],
                n_pos=len(pos), n_neg=len(neg))
