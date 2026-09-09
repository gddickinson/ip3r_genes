"""S22 shared layer — paths, loaders, and the statistics the tests are read from.

Three conventions every S22 stage inherits.

**Nothing here re-derives what an earlier task committed.**  The per-residue
constraint layers come from S17's `constraint_<gene>_<acc>.tsv`, the per-site
selection from its `fel_sites.tsv`, the deep alignments from its
`aln_<gene>.fasta`, and the element vocabulary from its `domain_map.tsv`.
S22's contribution is a *module* vocabulary on top of that, and the paired
tests that vocabulary makes possible.

**Every test is paired inside one alignment where it can be.**  The ligand
core and the pore module sit in the same protein, so a difference between
them measured on the same sequences cannot be a difference in taxon
sampling, alignment depth or ortholog quality.  `paired_by_tip` is the shape
that buys: one core number and one pore number per ortholog, ~250 of them
per paralogue.

**scipy is used for the two rank tests and for nothing else.**  The
permutation nulls are written out here because their unit of resampling is
the thing being tested (which residues carry the contact label), and a
library test cannot express that.
"""

from __future__ import annotations

import hashlib
import json
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as S17  # noqa: E402

PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results" / "ligand_site"
FIG_DIR = OUT_DIR / "figures"
CONSTRAINT_DIR = PROJECT_ROOT / "results" / "constraint"
SELECTION_DIR = PROJECT_ROOT / "results" / "selection"
MSA_DIR = PROJECT_ROOT / "results" / "msa_v2"
S20_DIR = PROJECT_ROOT / "results" / "s20_sweep"

PARALOGS = ("ITPR1", "ITPR2", "ITPR3")
REFERENCES = S17.REFERENCES
STRUCTURE_REF = "ITPR3"          # 6DQN / 8TK* are human ITPR3
STRUCTURE_ACC = "Q14573"

# The conservation layer every module comparison is primarily read on.
# `deep` is the per-paralogue sweep-ortholog layer (249-265 sequences);
# it is the only layer with enough depth for a per-tip paired test.
PRIMARY_LAYER = "deep"
LAYERS = ("deep", "vert", "family", "shallow")

read_fasta = S17.read_fasta
write_fasta = S17.write_fasta
read_tsv = S17.read_tsv
write_tsv = S17.write_tsv
fmt = S17.fmt


def log(*a) -> None:
    print("[s22]", *a, flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def live(stage: str, steps: list[tuple[str, bool]]) -> None:
    """Dashboard live panel (session protocol)."""
    payload = {
        "task": "S22",
        "stage": stage,
        "steps": [{"label": lab, "done": bool(done)} for lab, done in steps],
    }
    try:
        (PROJECT_ROOT / "results" / "session_live.json").write_text(
            json.dumps(payload, indent=2))
    except OSError:
        pass


# ---------------------------------------------------------------------------
# loaders — every one of these reads a committed table
# ---------------------------------------------------------------------------

_CONSTRAINT: dict[str, list[dict]] = {}


def constraint(paralog: str) -> list[dict]:
    """S17's per-residue table for one paralogue, in its own numbering."""
    if paralog not in _CONSTRAINT:
        acc = REFERENCES[paralog][1]
        rows = read_tsv(CONSTRAINT_DIR / f"constraint_{paralog}_{acc}.tsv")
        for r in rows:
            r["resi"] = int(r["resi"])
            r["deep_col"] = int(r["deep_col"]) if r["deep_col"] != "" else None
            for lay in LAYERS:
                v = r.get(f"{lay}_jsd", "")
                r[f"{lay}_jsd"] = float(v) if v not in ("", "nan") else None
                o = r.get(f"{lay}_occupancy", "")
                r[f"{lay}_occupancy"] = float(o) if o not in ("", "nan") else None
            for flag in ("ip3_contact", "filter_lining", "gate_lining"):
                r[flag] = r.get(flag, "") == "True"
        _CONSTRAINT[paralog] = rows
    return _CONSTRAINT[paralog]


_FEL: dict[str, list[dict]] | None = None


def fel_sites() -> dict[str, list[dict]]:
    """S17's FEL per-site rates, keyed by paralogue, on reference numbering."""
    global _FEL
    if _FEL is None:
        out: dict[str, list[dict]] = {p: [] for p in PARALOGS}
        for r in read_tsv(CONSTRAINT_DIR / "fel_sites.tsv"):
            if r["resi"] in ("", "None"):
                continue
            r["resi"] = int(r["resi"])
            for k in ("alpha", "beta", "omega", "p_value", "q_value"):
                v = r.get(k, "")
                r[k] = float(v) if v not in ("", "nan", "None") else None
            r["alpha_at_bound"] = r.get("alpha_at_bound", "") == "True"
            out.setdefault(r["paralog"], []).append(r)
        _FEL = out
    return _FEL


_DEEP: dict[str, dict[str, str]] = {}


def deep_alignment(paralog: str) -> dict[str, str]:
    """S17's per-paralogue deep ortholog alignment (the `deep` layer)."""
    if paralog not in _DEEP:
        _DEEP[paralog] = read_fasta(CONSTRAINT_DIR / f"aln_{paralog}.fasta")
    return _DEEP[paralog]


def deep_reference_label(paralog: str) -> str:
    return f"REF|{paralog}|{REFERENCES[paralog][1]}"


def domain_map() -> list[dict]:
    rows = read_tsv(CONSTRAINT_DIR / "domain_map.tsv")
    for r in rows:
        r["start"] = int(r["start"])
        r["end"] = int(r["end"])
    return rows


def functional_sites() -> list[dict]:
    rows = read_tsv(CONSTRAINT_DIR / "functional_sites.tsv")
    for r in rows:
        r["resi"] = int(r["resi"])
        r["source_resi"] = int(r["source_resi"])
    return rows


def msa_v2() -> tuple[dict[str, str], list[dict]]:
    return (read_fasta(MSA_DIR / "aln.fasta"),
            read_tsv(MSA_DIR / "representatives.tsv"))


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def median(xs) -> float | None:
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    n = len(xs)
    return xs[n // 2] if n % 2 else 0.5 * (xs[n // 2 - 1] + xs[n // 2])


def mean(xs) -> float | None:
    xs = [x for x in xs if x is not None]
    return sum(xs) / len(xs) if xs else None


def sign_test(diffs: list[float]) -> dict:
    """Two-sided exact sign test.  Ties are dropped and counted (S8's rule)."""
    pos = sum(1 for d in diffs if d > 0)
    neg = sum(1 for d in diffs if d < 0)
    ties = sum(1 for d in diffs if d == 0)
    n = pos + neg
    if n == 0:
        return {"n": 0, "n_pos": 0, "n_neg": 0, "n_ties": ties, "p": None}
    k = min(pos, neg)
    tail = sum(math.comb(n, i) for i in range(k + 1)) / (2 ** n)
    return {"n": n, "n_pos": pos, "n_neg": neg, "n_ties": ties,
            "p": min(1.0, 2 * tail)}


def wilcoxon(diffs: list[float]) -> dict:
    """Paired Wilcoxon signed-rank, two-sided, zeros dropped."""
    nz = [d for d in diffs if d != 0]
    if len(nz) < 6:
        return {"n": len(nz), "stat": None, "p": None,
                "note": "too few non-zero pairs"}
    from scipy.stats import wilcoxon as _w
    r = _w(nz, alternative="two-sided", zero_method="wilcox")
    return {"n": len(nz), "stat": float(r.statistic), "p": float(r.pvalue),
            "note": ""}


def mann_whitney(a: list[float], b: list[float],
                 alternative: str = "two-sided") -> dict:
    a = [x for x in a if x is not None]
    b = [x for x in b if x is not None]
    if not a or not b:
        return {"n_a": len(a), "n_b": len(b), "u": None, "p": None,
                "cles": None}
    from scipy.stats import mannwhitneyu
    r = mannwhitneyu(a, b, alternative=alternative)
    # Common-language effect size: P(a > b) + 0.5 P(a == b)
    cles = float(r.statistic) / (len(a) * len(b))
    return {"n_a": len(a), "n_b": len(b), "u": float(r.statistic),
            "p": float(r.pvalue), "cles": cles}


def permutation_label(values: list[float], labelled: list[bool],
                      iters: int = 100000, seed: int = 20220422,
                      statistic=None) -> dict:
    """Permute *which* positions carry the label, keeping the values fixed.

    The unit of resampling is the label, because the null being tested is
    "these particular residues are no more constrained than any other
    residue of the same module".  A test that resampled values instead
    would be asking a different question.
    """
    vals = [v for v, _ in zip(values, labelled) if v is not None]
    labs = [l for v, l in zip(values, labelled) if v is not None]
    k = sum(labs)
    if k == 0 or k == len(vals):
        return {"n": len(vals), "k": k, "observed": None, "p": None,
                "note": "label is empty or universal"}
    statistic = statistic or (lambda xs: sum(xs) / len(xs))
    obs = statistic([v for v, l in zip(vals, labs) if l])
    rng = random.Random(seed)
    idx = list(range(len(vals)))
    ge = 0
    for _ in range(iters):
        pick = rng.sample(idx, k)
        if statistic([vals[i] for i in pick]) >= obs:
            ge += 1
    return {"n": len(vals), "k": k, "observed": obs,
            "p": (ge + 1) / (iters + 1), "iters": iters, "note": ""}


def bootstrap_ci(xs: list[float], iters: int = 10000, seed: int = 20220422,
                 stat=None, alpha: float = 0.05) -> tuple[float, float] | tuple[None, None]:
    xs = [x for x in xs if x is not None]
    if len(xs) < 3:
        return (None, None)
    stat = stat or (lambda v: sum(v) / len(v))
    rng = random.Random(seed)
    n = len(xs)
    draws = sorted(stat([xs[rng.randrange(n)] for _ in range(n)])
                   for _ in range(iters))
    lo = draws[int(alpha / 2 * iters)]
    hi = draws[min(iters - 1, int((1 - alpha / 2) * iters))]
    return (lo, hi)


def benjamini_hochberg(pvals: list[float]) -> list[float]:
    """BH q-values; None passes through as None and takes no rank."""
    idx = [i for i, p in enumerate(pvals) if p is not None]
    m = len(idx)
    out: list[float | None] = [None] * len(pvals)
    if m == 0:
        return out
    order = sorted(idx, key=lambda i: pvals[i])
    prev = 1.0
    for rank in range(m, 0, -1):
        i = order[rank - 1]
        q = min(prev, pvals[i] * m / rank)
        out[i] = q
        prev = q
    return out


def power_binomial(n: int, p0: float, p1: float, alpha: float = 0.05) -> float:
    """Power of a two-sided sign test at n pairs to detect a shift to p1.

    Exact: the rejection region is computed from the null binomial and its
    probability taken under the alternative.  Returned so a lineage test
    with four taxa can say what it could and could not have seen.
    """
    if n <= 0:
        return 0.0
    from math import comb
    null = [comb(n, k) * p0 ** k * (1 - p0) ** (n - k) for k in range(n + 1)]
    # two-sided rejection region: smallest tails summing under alpha
    order = sorted(range(n + 1), key=lambda k: abs(k - n * p0), reverse=True)
    reject, acc = set(), 0.0
    for k in order:
        if acc + null[k] > alpha:
            break
        reject.add(k)
        acc += null[k]
    return sum(comb(n, k) * p1 ** k * (1 - p1) ** (n - k) for k in reject)
