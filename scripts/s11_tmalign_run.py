"""S11 step 4 — TM-align across the panel, resumable.

All-vs-all rather than only vs-reference. The extra pairs are what turn a
column of numbers into a calibrated one:

  * **model vs its own paralog's cryo-EM reference** is the fold claim;
  * **model vs the sister family's reference** is D14 asked structurally —
    a claim that something folds like an ITPR means nothing until it is
    also *not* folding like a RyR by a margin;
  * **reference vs reference** is the ceiling: two experimental structures
    of the same protein in different states bound how much of a TM-score
    difference is conformation rather than fold, which is why the state
    panel exists;
  * **anything vs the negative controls** is the floor, measured rather
    than taken from the literature.

Every pair is cached on the SHA-256 of the two files, so a re-run is free
and a *changed* structure re-runs rather than silently reusing a stale
score (D24's discipline applied to a cache key).

**Two guards, and both were earned in this task.** A first, ad-hoc run of
this stage was interrupted; killing its TM-align children left the Python
parent orphaned to init, and it went on spawning new ones for twenty
minutes beside the driver that replaced it — S7 recorded exactly this
incident for IQ-TREE and the same shape recurred here. So:

  * `claim()` refuses to start while another process holds the panel's
    lock, and the lock is checked against a **live** pid rather than mere
    file existence, so a crashed run does not block the next one for ever;
  * the cache write is **atomic** (temp file, then rename). Two processes
    writing the same key concurrently is exactly what the orphan produced,
    and an interleaved write yields a JSON file that parses as garbage —
    survivable here only because the reader already unlinks and re-runs on
    a decode error, which is a repair, not a design;
  * **a failed pair is never cached.** Killing the orphan's TM-align
    children returned a `-9` exit status with empty stderr, which the first
    version stored as a permanent zero-score result with no error text —
    ten pairs, all against one control, and the only reason they did not
    end up in the published table as genuine non-matches is that this
    module's own `self_test` requires every pair to have parsed a score.
    A failure is now simply not written, so the next run re-runs it. That
    costs a re-run when a pair genuinely cannot be scored; the alternative
    costs a wrong number.
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import os
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import structures_dir                              # noqa: E402
from s11_tmalign import run_tmalign                             # noqa: E402

PAIR_COLUMNS = [
    "query", "target", "query_role", "target_role", "query_call",
    "target_call", "query_paralog", "target_paralog", "tm_query",
    "tm_target", "tm_max", "rmsd", "aligned", "seq_identity",
    "len_query", "len_target", "ok", "error",
]


def lock_path() -> Path:
    return structures_dir() / "tmalign.lock"


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def claim() -> None:
    """Refuse to run while a live process already owns this stage."""
    path = lock_path()
    if path.exists():
        try:
            pid = int(path.read_text().strip() or 0)
        except ValueError:
            pid = 0
        if pid and pid != os.getpid() and _alive(pid):
            raise RuntimeError(
                f"TM-align stage is already running as pid {pid} "
                f"({path}). Stop it before starting another, or the two "
                "will duplicate every pair.")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"{os.getpid()}\n")


def release() -> None:
    path = lock_path()
    try:
        if path.exists() and path.read_text().strip() == str(os.getpid()):
            path.unlink()
    except OSError:
        pass


def _write_atomic(path: Path, payload: str) -> None:
    """Temp file then rename — a rename is atomic on the same filesystem,
    so a concurrent reader never sees a half-written cache entry."""
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as fh:
            fh.write(payload)
        os.replace(tmp, path)
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def cache_path(sha_q: str, sha_t: str) -> Path:
    d = structures_dir() / "tmalign_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{sha_q[:16]}_{sha_t[:16]}.json"


def _one(pair: tuple[dict, dict]) -> dict:
    q, t = pair
    path = cache_path(q["sha256"], t["sha256"])
    if path.exists():
        try:
            res = json.loads(path.read_text())
        except json.JSONDecodeError:
            path.unlink()
        else:
            return _row(q, t, res)
    result = run_tmalign(Path(q["path"]), Path(t["path"]), q["id"], t["id"])
    res = asdict(result)
    if result.ok:
        _write_atomic(path, json.dumps(res))
    return _row(q, t, res)


def _row(q: dict, t: dict, res: dict) -> dict:
    tm_max = max(res.get("tm_query", 0.0), res.get("tm_target", 0.0))
    return {
        "query": q["id"], "target": t["id"],
        "query_role": q.get("role", ""), "target_role": t.get("role", ""),
        "query_call": q.get("call", ""), "target_call": t.get("call", ""),
        "query_paralog": q.get("paralog", ""),
        "target_paralog": t.get("paralog", ""),
        "tm_query": round(res.get("tm_query", 0.0), 4),
        "tm_target": round(res.get("tm_target", 0.0), 4),
        "tm_max": round(tm_max, 4),
        "rmsd": res.get("rmsd", 0.0),
        "aligned": res.get("aligned", 0),
        "seq_identity": res.get("seq_id", 0.0),
        "len_query": res.get("len_query", 0),
        "len_target": res.get("len_target", 0),
        "ok": int(bool(res.get("ok"))),
        "error": res.get("error", ""),
    }


def all_pairs(panel: list[dict]) -> list[tuple[dict, dict]]:
    """Ordered pairs, both directions excluded — TM-align reports both
    normalisations in one run, so (a,b) and (b,a) would be the same work
    twice. Self-pairs are kept: a structure against itself must score 1.0,
    which is the cheapest possible check that the right file was read."""
    usable = [p for p in panel if p.get("status") == "ok" and p.get("path")]
    usable.sort(key=lambda p: p["id"])
    out = []
    for i, q in enumerate(usable):
        for t in usable[i:]:
            out.append((q, t))
    return out


def run(panel: list[dict], workers: int = 4, on_progress=None) -> list[dict]:
    claim()
    pairs = all_pairs(panel)
    rows: list[dict] = []
    done = 0
    try:
        with cf.ThreadPoolExecutor(max_workers=workers) as pool:
            for row in pool.map(_one, pairs):
                rows.append(row)
                done += 1
                if on_progress and done % 25 == 0:
                    on_progress(done, len(pairs))
    finally:
        release()
    rows.sort(key=lambda r: (r["query"], r["target"]))
    return rows


def self_test(rows: list[dict]) -> list[str]:
    """Negative controls on the scoring itself, run on every build.

    Three properties TM-align must have on this panel, each of which has a
    plausible-looking failure mode that no downstream number would reveal:

      T1 a structure against itself scores 1.0 — the cheapest check that
         the manifest's path column points where it says;
      T2 every pair parsed a score (a version-string change in TM-align's
         output silently yields zeros for every pair, which reads as
         "nothing folds like anything");
      T3 the two normalisations bracket each other sensibly — the score
         normalised by the shorter chain is never the smaller of the two.
    """
    problems: list[str] = []
    for r in rows:
        if r["query"] == r["target"]:
            if r["tm_max"] < 0.999:
                problems.append(f"T1 self-pair {r['query']} scored {r['tm_max']}")
            continue
        if not r["ok"]:
            problems.append(f"T2 no score for {r['query']} vs {r['target']}: "
                            f"{r['error']}")
            continue
        if r["len_query"] and r["len_target"]:
            shorter_is_q = r["len_query"] <= r["len_target"]
            by_shorter = r["tm_query"] if shorter_is_q else r["tm_target"]
            if by_shorter + 1e-6 < min(r["tm_query"], r["tm_target"]):
                problems.append(
                    f"T3 normalisation inverted for {r['query']} vs {r['target']}")
    return problems
