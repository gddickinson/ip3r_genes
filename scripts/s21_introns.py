"""Intron positions: what is conserved, and what is shared between paralogues.

An intron's *position* here is a pair — the alignment column its upstream exon
ends in, and its classical phase (0, 1 or 2: how many bases of the interrupted
codon lie upstream). Both halves are needed. Two paralogues can carry an intron
between the same two residues in different frames, which is not one ancestral
intron; and a column on its own would call any two introns in the same region
shared.

**Within a paralogue** the question is whether the ~58-exon architecture is one
architecture: how many positions occur in what fraction of the genomes carrying
that paralogue, measured over the sweep rather than asserted from human.

**Between paralogues** the question is how many of those positions are
ancestral, and it needs a null, because two genes with ~58 introns each spread
over ~2,700 aligned residues will share some positions by chance. Every
comparison is **paired within genome** (D16) and scored two ways, which is the
point:

* an **exact Poisson-binomial** null — each of B's introns placed
  independently and uniformly on the columns *both* loci have residues in,
  keeping its own phase, so the match probability differs per intron and the
  distribution of the match count is the Poisson-binomial of those
  probabilities. No RNG, no seed, no replicate count: the p-value is exact
  under the stated null.
* a **seeded permutation** null that relocates B's introns without
  replacement, which a real gene obeys and the analytic null does not.

Both are committed. They agree to within the permutation's own resolution, and
where they would not, the permutation is the one to believe.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as S16                                             # noqa: E402
import s21_lib as L                                               # noqa: E402

#: Column tolerances a shared position is tested at. 0 is the operating point;
#: the others say whether the answer rests on exact column identity, which in a
#: 11,777-column alignment of four kingdoms it should not have to.
TOLERANCES = (0, 1, 2)
OPERATING_TOL = 0

#: Replicates for the permutation cross-check, and its seed (D24).
N_PERM = 1000
SEED = 20260908

#: A position is called conserved within a paralogue at these prevalences.
PREVALENCE_BARS = (0.50, 0.90, 0.99)


def available_columns(frames: dict, bait: str, cell: str,
                      q_start: int, q_end: int) -> set[int]:
    """Every alignment column this locus has a residue in.

    The columns the null draws from: a locus that aligns residues 200-2,600 of
    its bait cannot carry an intron at a column outside that range, and a null
    drawn from the whole alignment would make every real match look surprising.
    """
    m = (frames.get(bait) or {}).get(cell)
    if not m:
        return set()
    return {c for r, c in m["map"].items() if q_start <= r <= q_end}


def positions(irows: list[dict]) -> dict[tuple[str, str, int], set]:
    """`(accession, cell, locus_idx) -> {(column, phase)}`, columns known only."""
    out: dict[tuple[str, str, int], set] = {}
    for r in irows:
        if r["col_left"] <= 0 or r["intron_phase"] < 0:
            continue
        out.setdefault((r["accession"], r["cell"], r["locus_idx"]),
                       set()).add((r["col_left"], r["intron_phase"]))
    return out


# --------------------------------------------------------------------------
# within a paralogue
# --------------------------------------------------------------------------
def conservation(pos: dict, arows: list[dict]) -> list[dict]:
    """One row per (cell, position): in how many of that cell's genomes it is.

    One locus per (genome, cell) — the best-covered — so a genome with two
    copies does not vote twice.
    """
    best: dict[tuple[str, str], tuple] = {}
    cov: dict[tuple[str, str], float] = {}
    for r in arows:
        if not r["in_scope"]:
            continue
        k = (r["accession"], r["cell"])
        if k not in cov or r["coverage"] > cov[k]:
            cov[k] = r["coverage"]
            best[k] = (r["accession"], r["cell"], r["locus_idx"])
    per_cell: dict[str, list[set]] = {}
    for (acc, cell), key in best.items():
        if key in pos:
            per_cell.setdefault(cell, []).append(pos[key])
    rows = []
    for cell, sets in sorted(per_cell.items()):
        n = len(sets)
        counts: dict[tuple[int, int], int] = {}
        for s in sets:
            for p in s:
                counts[p] = counts.get(p, 0) + 1
        for (col, phase), c in sorted(counts.items()):
            rows.append({"cell": cell, "column": col, "phase": phase,
                         "n_loci_with": c, "n_loci": n,
                         "prevalence": round(c / n, 4) if n else 0.0})
    return rows


def conservation_summary(rows: list[dict]) -> list[dict]:
    """How many positions clear each prevalence bar, per paralogue."""
    out = []
    for cell in sorted({r["cell"] for r in rows}):
        sub = [r for r in rows if r["cell"] == cell]
        n_loci = sub[0]["n_loci"] if sub else 0
        d = {"cell": cell, "n_loci": n_loci, "n_positions_seen": len(sub)}
        for bar in PREVALENCE_BARS:
            d[f"n_at_{int(bar * 100)}pc"] = sum(1 for r in sub
                                                if r["prevalence"] >= bar)
        d["median_positions_per_locus"] = round(
            sum(r["n_loci_with"] for r in sub) / n_loci, 1) if n_loci else 0.0
        out.append(d)
    return out


# --------------------------------------------------------------------------
# between paralogues
# --------------------------------------------------------------------------
def _targets(pos_a: set, tol: int) -> dict[int, set[int]]:
    """Phase -> the columns within `tol` of an A intron of that phase."""
    out: dict[int, set[int]] = {}
    for col, phase in pos_a:
        s = out.setdefault(phase, set())
        for d in range(-tol, tol + 1):
            s.add(col + d)
    return out


def _poisson_binomial_tail(probs: list[float], k: int) -> float:
    """P(X >= k) for independent Bernoullis with these probabilities, exactly.

    `k <= 0` returns exactly 1.0 rather than the sum of the whole distribution:
    summing 60 terms leaves 1 - 9e-16, and a p-value printed as 0.9999999999999991
    where the answer is "no signal at all" invites the reading that something
    was measured.
    """
    if k <= 0:
        return 1.0
    dist = [1.0]
    for p in probs:
        nxt = [0.0] * (len(dist) + 1)
        for i, v in enumerate(dist):
            nxt[i] += v * (1.0 - p)
            nxt[i + 1] += v * p
        dist = nxt
    return min(1.0, max(0.0, sum(dist[min(k, len(dist) - 1):])))


def shared_pair(pos_a: set, pos_b: set, shared_cols: set[int], tol: int,
                rng: random.Random | None = None) -> dict:
    """One paralogue pair in one genome: observed matches and both nulls."""
    tgt = _targets(pos_a, tol)
    b_in = [(c, p) for c, p in pos_b if c in shared_cols]
    observed = sum(1 for c, p in b_in if c in tgt.get(p, ()))
    n_shared = len(shared_cols)
    probs = []
    for _c, phase in b_in:
        hits = len(tgt.get(phase, set()) & shared_cols)
        probs.append(hits / n_shared if n_shared else 0.0)
    exp_analytic = sum(probs)
    p_analytic = (_poisson_binomial_tail(probs, observed)
                  if probs else float("nan"))
    perm_ge = perm_sum = 0
    ran_perm = rng is not None and bool(b_in) and n_shared >= len(b_in)
    if ran_perm:
        cols = sorted(shared_cols)
        phases = [p for _c, p in b_in]
        for _ in range(N_PERM):
            draw = rng.sample(cols, len(b_in))
            m = sum(1 for c, ph in zip(draw, phases) if c in tgt.get(ph, ()))
            perm_sum += m
            perm_ge += m >= observed
    return {"n_introns_a": len(pos_a), "n_introns_b": len(pos_b),
            "n_introns_b_in_shared": len(b_in), "n_shared_columns": n_shared,
            "observed": observed,
            "expected_analytic": round(exp_analytic, 3),
            "p_analytic": p_analytic,
            "n_perm": N_PERM if ran_perm else 0,
            "expected_perm": round(perm_sum / N_PERM, 3) if ran_perm
            else float("nan"),
            "p_perm": round((perm_ge + 1) / (N_PERM + 1), 5) if ran_perm
            else float("nan")}


def shared_table(pos: dict, arows: list[dict], frames: dict,
                 tolerances=TOLERANCES, permute: bool = True) -> list[dict]:
    """Every paralogue pair in every genome that carries both, at each tolerance."""
    best: dict[tuple[str, str], dict] = {}
    for r in arows:
        if not r["in_scope"]:
            continue
        k = (r["accession"], r["cell"])
        if k not in best or r["coverage"] > best[k]["coverage"]:
            best[k] = r
    cells = list(L.PARALOGS) + [L.CONTROL_CELL]
    cols_cache: dict[tuple[str, str], set[int]] = {}

    def avail(row: dict) -> set[int]:
        k = (row["accession"], row["cell"])
        if k not in cols_cache:
            cols_cache[k] = available_columns(
                frames, row["bait"], row["cell"],
                *_q_range(row))
        return cols_cache[k]

    rows = []
    for tol in tolerances:
        rng = random.Random(SEED + tol) if permute else None
        for i, a in enumerate(cells):
            for b in cells[i + 1:]:
                for acc in sorted({k[0] for k in best}):
                    ra, rb = best.get((acc, a)), best.get((acc, b))
                    if not ra or not rb:
                        continue
                    ka = (acc, a, ra["locus_idx"])
                    kb = (acc, b, rb["locus_idx"])
                    if ka not in pos or kb not in pos:
                        continue
                    shared = avail(ra) & avail(rb)
                    if not shared:
                        continue
                    d = shared_pair(pos[ka], pos[kb], shared, tol, rng)
                    d.update(accession=acc, organism=ra["organism"],
                             vclass=ra["vclass"], cell_a=a, cell_b=b,
                             tolerance=tol,
                             is_operating_point=int(tol == OPERATING_TOL),
                             frame_via_cell=int(bool(ra["frame_via_cell"]
                                                     or rb["frame_via_cell"])))
                    rows.append(d)
    return rows


def _q_range(row: dict) -> tuple[int, int]:
    """The bait residues this locus actually aligned.

    Not the whole bait: a locus aligning residues 200-2,600 cannot carry an
    intron outside that range, and drawing the null from the bait's full column
    image would make the shared set too large, the expected match count too
    small and every real match look more surprising than it is.
    """
    return (int(row["q_start"]), int(row["q_end"]))


def shared_summary(rows: list[dict], drop_via_cell: bool = False) -> list[dict]:
    """Per (pair, tolerance): the pooled counts and the BH-corrected verdict."""
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        if drop_via_cell and r["frame_via_cell"]:
            continue
        groups.setdefault((r["cell_a"], r["cell_b"], r["tolerance"]),
                          []).append(r)
    out = []
    for (a, b, tol), sub in sorted(groups.items()):
        obs = sum(r["observed"] for r in sub)
        exp = sum(r["expected_analytic"] for r in sub)
        n_sig = sum(1 for r in sub if r["p_analytic"] == r["p_analytic"]
                    and r["p_analytic"] < 0.05)
        out.append({
            "cell_a": a, "cell_b": b, "tolerance": tol, "n_genomes": len(sub),
            "observed_total": obs, "expected_total": round(exp, 1),
            "enrichment": round(obs / exp, 3) if exp else float("nan"),
            "median_observed": round(L.quantile(
                [float(r["observed"]) for r in sub], 0.5), 1),
            "median_expected": round(L.quantile(
                [float(r["expected_analytic"]) for r in sub], 0.5), 2),
            "n_genomes_p_lt_0.05": n_sig,
            "frac_genomes_significant": round(n_sig / len(sub), 4) if sub else 0.0,
            "excludes_frame_via_cell": int(drop_via_cell),
            "is_operating_point": int(tol == OPERATING_TOL)})
    ps = []
    for r in out:
        sub = groups[(r["cell_a"], r["cell_b"], r["tolerance"])]
        worst = max((x["p_analytic"] for x in sub
                     if x["p_analytic"] == x["p_analytic"]), default=1.0)
        ps.append(worst)
    for r, q in zip(out, S16.benjamini_hochberg(ps)):
        r["q_worst_genome"] = round(q, 6)
    return out
