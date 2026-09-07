"""S9 step 1 — fetch + validate a CDS for every vertebrate representative.

Writes `results/selection/cds.fasta` (nucleotide, headers = the S6/S7 tip
labels, internal-stop codons masked to NNN) and `cds_status.tsv` (one row
per tip: route, source id, validation stats, kept/failed).

Tips are fetched concurrently — the REST routes are network-bound and
rest.ensembl.org is slow in bursts — with every raw response cached, so a
re-run only repeats work that previously failed. The miniprot route is
cached at the *result* level as well: a locus realignment costs a genome
index plus an alignment, which the HTTP cache cannot help with.

Run:  python scripts/s9_fetch_cds.py [--workers N] [--miniprot-threads N]
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import (ALN_FASTA, CACHE_DIR, OUT_DIR, ROOT,  # noqa: E402
                        cds_from_miniprot_model, ensembl_candidates,
                        load_vertebrate_tips, read_fasta, route_for,
                        uniprot_candidates, validate_and_mask, write_fasta)

#: Ceiling on how many CDS candidates one tip may be offered. A UniProt
#: entry can list a dozen transcripts; past this the route is not going to
#: produce the aligned isoform and the tip is better recorded as failed.
MAX_CANDIDATES = 12

LIVE_JSON = ROOT / "results" / "session_live.json"
_lock = threading.Lock()

STATUS_COLS = ["label", "group", "route", "source", "status", "cds_len",
               "prot_len", "n_mismatch", "n_stop_masked", "n_masked",
               "n_candidates", "rejected", "note"]


def live(stage: str, done: int, total: int, note: str = "") -> None:
    LIVE_JSON.write_text(json.dumps({
        "task": "S9", "stage": stage,
        "steps": [{"label": stage, "done": done, "total": total}],
        "note": note, "ts": time.strftime("%H:%M:%S")}))


def _miniprot_candidates(label: str, kw: dict, expected: str, threads: int):
    """The miniprot route, cached per tip — a realignment is minutes, not
    milliseconds, and the HTTP cache does not cover it. It offers at most
    one candidate: the locus rerun either reproduces the S5 protein exactly
    or it has found a different gene model."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in label)
    cache = CACHE_DIR / f"miniprot_{safe}.json"
    if cache.exists():
        data = json.loads(cache.read_text())
        if data.get("cds"):
            yield data["cds"], data.get("note", "cached")
            return
    cds, note = cds_from_miniprot_model(kw["assembly"], expected,
                                        threads=threads,
                                        locus=kw.get("locus"))
    if cds:
        cache.write_text(json.dumps({"cds": cds, "note": note}))
        yield cds, note


def fetch_one(row: dict, expected: str | None, mp_threads: int) -> dict:
    """Resolve one representative to a validated CDS.

    Every candidate a route offers is translated and checked against the S6
    protein, and the first that *matches* is kept — not the first that
    downloads. The ones rejected on the way are recorded, because "this tip
    has three transcripts and only one encodes the aligned isoform" is a
    fact about the record, not noise.
    """
    label = row["label"]
    route, kw = route_for(row)
    base = {"label": label, "group": row["group"], "route": route}
    if expected is None:
        return {**base, "source": "-", "status": "no_msa_seq", "cds": None,
                "stats": {}, "rejected": "", "note": "label not in aln.fasta"}
    if route == "ensembl":
        cands = ensembl_candidates(kw["ensembl_id"])
    elif route == "uniprot":
        cands = uniprot_candidates(kw["acc"])
    else:
        cands = _miniprot_candidates(label, kw, expected, mp_threads)

    tried = 0
    rejected: list[str] = []
    last_stats: dict = {"prot_len": len(expected)}
    for cds, note in cands:
        tried += 1
        masked, st = validate_and_mask(cds, expected)
        if masked:
            return {**base, "source": note, "status": "ok", "cds": masked,
                    "stats": {**st, "n_candidates": tried},
                    "rejected": ";".join(rejected), "note": st.get("note", "")}
        last_stats = {**st, "n_candidates": tried}
        rejected.append(f"{note}({st.get('note') or 'mismatch'},"
                        f"{st.get('trans_len', 0)}aa)")
        if tried >= MAX_CANDIDATES:
            break
    status = "fetch_failed" if not tried else "validation_failed"
    return {**base, "source": rejected[-1] if rejected else "-",
            "status": status, "cds": None, "stats": last_stats,
            "rejected": ";".join(rejected),
            "note": last_stats.get("note", "no candidate matched")}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--miniprot-threads", type=int, default=3)
    args = ap.parse_args()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    tips = load_vertebrate_tips()
    aln = read_fasta(ALN_FASTA)
    # The ungapped alignment row *is* the protein each CDS must encode.
    # Taking it from representatives.fasta instead would be a second copy
    # that could differ from the one the tree was built on.
    prots = {lab: seq.replace("-", "").upper() for lab, seq in aln.items()}

    done = 0
    total = len(tips)

    def work(row: dict) -> dict:
        nonlocal done
        out = fetch_one(row, prots.get(row["label"]), args.miniprot_threads)
        with _lock:
            done += 1
            live("fetch_cds", done, total, out["label"])
            st = out["stats"]
            print(f"[{done:2}/{total}] {out['status']:18} {out['label'][:56]:56} "
                  f"mism={st.get('n_mismatch', '-')} "
                  f"stops={st.get('n_stop_masked', '-')} "
                  f"n={st.get('n_candidates', 0)} ({out['source'][:56]})",
                  flush=True)
        return out

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        results = list(pool.map(work, tips))

    kept = [(r["label"], r["cds"]) for r in results if r["cds"]]
    rows = ["\t".join(STATUS_COLS)]
    for r in results:
        st = r["stats"]
        rows.append("\t".join(str(x) for x in [
            r["label"], r["group"], r["route"], r["source"], r["status"],
            st.get("cds_len", 0), st.get("prot_len", 0),
            st.get("n_mismatch", "-"), st.get("n_stop_masked", "-"),
            st.get("n_masked", "-"), st.get("n_candidates", 0),
            r.get("rejected", ""), r["note"]]))
    write_fasta(OUT_DIR / "cds.fasta", kept)
    (OUT_DIR / "cds_status.tsv").write_text("\n".join(rows) + "\n")
    live("fetch_cds", total, total, f"{len(kept)}/{total} CDS validated")

    by_status: dict[str, int] = {}
    for r in results:
        by_status[r["status"]] = by_status.get(r["status"], 0) + 1
    print(f"\nkept {len(kept)}/{total} -> {OUT_DIR / 'cds.fasta'}")
    for k, v in sorted(by_status.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
