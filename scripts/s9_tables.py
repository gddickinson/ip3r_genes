"""S9 — every committed table, written from the codeml/HyPhy output once.

D13 applies to reports; this module is what makes that possible for S9.
The report and the figures read these TSVs and never re-parse an `mlc`, so
a number in the text cannot differ from the number in the table beside it.

Four things here are decisions rather than bookkeeping:

* **Branch-site model A is tested against a 50:50 mixture**, not a plain
  χ²₁. The null fixes ω₂ = 1 at the boundary of the parameter space, so the
  asymptotic distribution of 2ΔlnL is ½χ²₀ + ½χ²₁ and the naive p-value is
  twice too small. Both are written out; the halved one is the reported
  one, which is codeml's own recommendation.
* **An alternative that lands below its own null is not a result.** Nested
  models cannot do that, so such a run is a local optimum. Every restart is
  written to `bs_restarts.tsv` with a flag, and the LRT uses the best.
* **Pairwise dS carries a saturation flag.** Between 2R paralogs dS is
  expected to be past the point where it can be estimated, and an ω built
  on a saturated dS is a ratio of a small number to an unknown one.
* **BH is applied across the whole LRT family**, because S9 runs the same
  test three times — once per paralog — and reporting the smallest of three
  p-values without correction is the multiple-testing error this project
  spends its self-tests avoiding.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import OUT_DIR, PARALOGS  # noqa: E402
from s9_codeml_lib import benjamini_hochberg, lrt, parse_mlc  # noqa: E402
from s9_jobs import BS_INIT_OMEGAS, CODEML_DIR  # noqa: E402

#: dS above which a pairwise ω is a ratio to an unknown quantity. The
#: conventional bar; the table carries dS itself so a reader can move it.
SATURATED_DS = 1.5


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c, "")) for c in cols) + "\n")


def _fmt(v) -> str:
    if v is None or v == "":
        return ""
    if isinstance(v, float):
        return f"{v:.6g}"
    return str(v)


def load_results() -> dict:
    """The driver's json, back-filled by a scan of the job directories.

    A job run outside the driver lives only on disk, and a stale json must
    not be able to hide it — nor an `mlc` be preferred to a recorded run.
    """
    out: dict = {}
    path = OUT_DIR / "codeml_results.json"
    if path.exists():
        out.update(json.loads(path.read_text()))
    # The per-job files win over the shared json: two drivers may have run
    # at once, and only the per-job file is written by the driver that owns
    # that job.
    for d in sorted(p for p in CODEML_DIR.glob("*") if p.is_dir()):
        rj = d / "result.json"
        if rj.exists():
            try:
                out[d.name] = json.loads(rj.read_text())
                continue
            except Exception:  # noqa: BLE001
                pass
    for d in sorted(p for p in CODEML_DIR.glob("*") if p.is_dir()):
        if out.get(d.name, {}).get("lnL"):
            continue
        res = parse_mlc(d / "mlc", d.name)
        if res:
            out[d.name] = {"lnL": res.lnL, "np": res.np, "omegas": res.omegas,
                           "kappa": res.kappa, "tree_length": res.tree_length,
                           "site_classes": res.site_classes,
                           "description": "(recovered from mlc)"}
    return out


# ---- ω table ---------------------------------------------------------------

OMEGA_COLS = ["job", "set", "model", "lnL", "np", "kappa", "tree_length",
              "omega", "omegas", "minutes", "description"]


def _set_of(job: str) -> str:
    for p in PARALOGS:
        if job.endswith(p) or f"_{p}_" in job or job.endswith(f"{p}_curated"):
            return p
    return "all" if job.endswith("_all") else ""


def omega_rows(res: dict) -> list[dict]:
    rows = []
    for job, r in sorted(res.items()):
        oms = r.get("omegas") or []
        rows.append({
            "job": job, "set": _set_of(job),
            "model": job.split("_")[0],
            "lnL": r.get("lnL"), "np": r.get("np"), "kappa": r.get("kappa"),
            "tree_length": r.get("tree_length"),
            "omega": oms[0] if oms else "",
            "omegas": ";".join(f"{x:.4f}" for x in oms),
            "minutes": r.get("minutes", ""),
            "description": r.get("description", "")})
    return rows


# ---- pairwise dN/dS --------------------------------------------------------

def read_paml_matrix(path: Path) -> tuple[list[str], list[list[float]]]:
    """PAML lower-triangular distance matrix (2ML.dN / 2ML.dS)."""
    if not path.exists():
        return [], []
    lines = [l for l in path.read_text().split("\n") if l.strip()]
    names: list[str] = []
    rows: list[list[float]] = []
    for line in lines[1:]:
        parts = line.split()
        if not parts:
            continue
        names.append(parts[0])
        rows.append([float(x) for x in parts[1:]])
    return names, rows


PAIR_COLS = ["set", "a", "b", "dN", "dS", "omega", "saturated"]


def pairwise_rows() -> list[dict]:
    out: list[dict] = []
    for para in PARALOGS:
        wd = CODEML_DIR / f"pair_{para}"
        n_dn, m_dn = read_paml_matrix(wd / "2ML.dN")
        n_ds, m_ds = read_paml_matrix(wd / "2ML.dS")
        if not n_dn or n_dn != n_ds:
            continue
        for i, name_i in enumerate(n_dn):
            for j in range(len(m_dn[i])):
                dn, ds = m_dn[i][j], m_ds[i][j]
                out.append({"set": para, "a": name_i, "b": n_dn[j],
                            "dN": dn, "dS": ds,
                            "omega": (dn / ds) if ds > 1e-6 else "",
                            "saturated": int(ds > SATURATED_DS)})
    return out


# ---- branch-site restarts --------------------------------------------------

BS_COLS = ["paralog", "run", "initial_omega", "lnL", "np",
           "lnL_null", "below_null", "is_best", "foreground_omega"]


def foreground_omega(job: str) -> str:
    mlc = CODEML_DIR / job / "mlc"
    if not mlc.exists():
        return ""
    m = re.search(r"^foreground w\s+(.*)$", mlc.read_text(), re.M)
    return ";".join(m.group(1).split()) if m else ""


def bs_rows(res: dict) -> list[dict]:
    rows: list[dict] = []
    for para in PARALOGS:
        null = res.get(f"bs_null_{para}")
        if not null:
            continue
        runs = [(f"bs_alt_{para}_w{w:g}", w) for w in BS_INIT_OMEGAS]
        runs = [(n, w) for n, w in runs if n in res]
        if not runs:
            continue
        best = max(runs, key=lambda t: res[t[0]]["lnL"])[0]
        for name, w0 in runs:
            r = res[name]
            rows.append({"paralog": para, "run": name, "initial_omega": w0,
                         "lnL": r["lnL"], "np": r.get("np"),
                         "lnL_null": null["lnL"],
                         "below_null": int(r["lnL"] < null["lnL"]),
                         "is_best": int(name == best),
                         "foreground_omega": foreground_omega(name)})
    return rows


def best_bs_alt(res: dict, para: str) -> str | None:
    runs = [f"bs_alt_{para}_w{w:g}" for w in BS_INIT_OMEGAS]
    runs = [n for n in runs if n in res]
    return max(runs, key=lambda n: res[n]["lnL"]) if runs else None


# ---- BEB sites -------------------------------------------------------------

BEB_COLS = ["job", "set", "n_sites", "n_p95", "n_p99", "median_p", "max_p"]


def beb_row(job: str) -> dict | None:
    """Spread of the BEB posterior P(ω > 1) over the alignment.

    A localised adaptive signal puts a handful of sites near 1 and leaves
    the rest low. A posterior that is flat and mid-valued everywhere means
    the model found no site to point at — the signature of a fit driven by
    saturation or alignment error rather than by identifiable substitutions.
    """
    mlc = CODEML_DIR / job / "mlc"
    if not mlc.exists():
        return None
    text = mlc.read_text()
    i = text.find("Bayes Empirical Bayes (BEB)")
    if i < 0:
        return None
    probs: list[float] = []
    for line in text[i:].split("\n")[2:]:
        m = re.match(r"^\s*\d+\s+[A-Z-]\s+([0-9.]+)", line)
        if not m:
            if probs:
                break
            continue
        probs.append(float(m.group(1)))
    if not probs:
        return None
    probs.sort()
    return {"job": job, "set": _set_of(job), "n_sites": len(probs),
            "n_p95": sum(p >= 0.95 for p in probs),
            "n_p99": sum(p >= 0.99 for p in probs),
            "median_p": probs[len(probs) // 2], "max_p": probs[-1]}


# ---- LRTs ------------------------------------------------------------------

LRT_COLS = ["test", "set", "null", "alt", "lnL_null", "lnL_alt", "stat",
             "df", "p", "p_reported", "q_bh", "mixture", "note"]


def lrt_rows(res: dict) -> list[dict]:
    rows: list[dict] = []

    def add(test: str, para: str, null: str, alt: str, df: int,
            mixture: bool = False, note: str = "") -> None:
        if null not in res or alt not in res:
            return
        n, a = res[null], res[alt]
        from types import SimpleNamespace
        st = lrt(SimpleNamespace(lnL=n["lnL"], np=n.get("np", 0)),
                 SimpleNamespace(lnL=a["lnL"], np=a.get("np", 0)), df)
        p_rep = st["p"] / 2 if mixture else st["p"]
        rows.append({"test": test, "set": para, "null": null, "alt": alt,
                     "lnL_null": n["lnL"], "lnL_alt": a["lnL"],
                     "stat": st["stat"], "df": df, "p": st["p"],
                     "p_reported": p_rep,
                     "mixture": int(mixture), "note": note})

    for para in PARALOGS:
        add(f"two-ratio ({para} clade) vs one-ratio", para,
            "m0_all", f"two_ratio_{para}", 1)
        best = best_bs_alt(res, para)
        if best:
            below = res[best]["lnL"] < res.get(f"bs_null_{para}", {}).get(
                "lnL", float("-inf"))
            add(f"branch-site model A on the {para} stem", para,
                f"bs_null_{para}", best, 1, mixture=True,
                note=("best restart still below its own null — local "
                      "optimum, not evidence" if below else
                      f"best of {len(BS_INIT_OMEGAS)} restarts"))
        add(f"M2a vs M1a within {para}", para, f"m1a_{para}", f"m2a_{para}", 2)
        add(f"M8 vs M7 within {para}", para, f"m7_{para}", f"m8_{para}", 2)

    if rows:
        qs = benjamini_hochberg([r["p_reported"] for r in rows])
        for r, q in zip(rows, qs):
            r["q_bh"] = q
    return rows


# ---- driver ----------------------------------------------------------------

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    res = load_results()
    if not res:
        raise SystemExit("no codeml results yet — run scripts/s9_codeml.py")

    write_tsv(OUT_DIR / "omega_table.tsv", OMEGA_COLS, omega_rows(res))
    write_tsv(OUT_DIR / "pairwise_dnds.tsv", PAIR_COLS, pairwise_rows())
    write_tsv(OUT_DIR / "bs_restarts.tsv", BS_COLS, bs_rows(res))
    beb = [r for r in (beb_row(f"{m}_{p}") for p in PARALOGS
                       for m in ("m2a", "m8")) if r]
    beb += [r for r in (beb_row(best_bs_alt(res, p) or "") for p in PARALOGS)
            if r]
    write_tsv(OUT_DIR / "beb_sites.tsv", BEB_COLS, beb)
    rows = lrt_rows(res)
    write_tsv(OUT_DIR / "lrt_table.tsv", LRT_COLS, rows)

    stats = {"n_jobs": len(res),
             "n_lrt": len(rows),
             "saturated_ds_bar": SATURATED_DS,
             "bs_init_omegas": list(BS_INIT_OMEGAS),
             "sha256": {p.name: sha256(p) for p in sorted(OUT_DIR.glob("*.tsv"))}}
    (OUT_DIR / "selection_stats.json").write_text(
        json.dumps(stats, indent=2, sort_keys=True))
    print(f"{len(res)} codeml jobs, {len(rows)} LRTs -> {OUT_DIR}")
    for r in rows:
        print(f"  {r['test']:46} 2dlnL={r['stat']:9.3f} df={r['df']} "
              f"p={r['p_reported']:.3g} q={r.get('q_bh', float('nan')):.3g}")


if __name__ == "__main__":
    main()
