"""S11 step 2 — how much of this family does AlphaFold DB actually hold?

Before comparing anything we ask what already exists, because the answer
sets the shape of the session: if AFDB covers the representatives, S11 is a
comparison exercise; if it does not, the absence *is* the result and the
comparison runs on whatever remains.

This family makes the question sharp. A vertebrate ITPR subunit is
~2,700 residues, sitting right on AFDB's monomer ceiling, and the sister
family at ~5,000 is past it outright. So coverage is measured, never
assumed, and it is measured in three ways that do not agree:

  * **`in_alphafold`** — UniProt's own cross-reference, already a column on
    every census v6 row. Cheap, family-wide, and second-hand.
  * **the API probe** — asking AFDB itself, per accession. Authoritative,
    and the only one that can see *what* is served.
  * **model coverage** — the modelled span against the length the census
    holds. This is the one that matters, and the one the first two cannot
    see: AFDB answers an accession query with whatever record it has,
    including an **isoform**. Two of the three human paralogs come back as
    isoform models, and for ITPR2 the served model is 181 residues of a
    2,701-residue protein. Counting that as "covered" would overstate the
    family's structural coverage by an order of magnitude on that row.

Two probe sets:
  * `reps`   — the 135 S6 representatives, S11's actual targets.
  * `census` — every census v6 record the family call assigns, restricted
               to UniProt-shaped accessions (a genome gene model has no
               AFDB entry by construction and is counted as *not
               applicable*, never as a miss).

    python scripts/s11_afdb_probe.py --set reps
    python scripts/s11_afdb_probe.py --set census
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import (PROJECT_ROOT, afdb_probe, load_census,     # noqa: E402
                     load_representatives, results_dir, write_tsv)
from s6_lib import normalise_group                              # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT))
from src.utils import family                                    # noqa: E402

#: The census carries its group column in two spellings (S23 left
#: `Viridiplantae` and `viridiplantae` both in place); counted naively that
#: is two groups. `s6_lib.normalise_group` is the project's one answer.

#: UniProt accession shape. Anything else — an Ensembl id, a RefSeq id, a
#: genome gene-model id from S5/S23 — cannot have an AFDB entry, so it is
#: `not_applicable` and is kept out of the coverage denominator.
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$")


def accession_kind(acc: str) -> str:
    return "uniprot" if UNIPROT_RE.match(acc or "") else "other"


def rep_targets() -> list[dict]:
    out = []
    for r in load_representatives():
        acc = (r.get("accession") or "").strip()
        out.append({
            "accession": acc,
            "label": r.get("label", ""),
            "kind": accession_kind(acc),
            "group": normalise_group(r.get("group", "")),
            "paralog": r.get("paralog", ""),
            "species": r.get("species", ""),
            "kingdom": r.get("kingdom", ""),
            "phylum": r.get("phylum", ""),
            "length": int(r.get("length") or 0),
            "source": r.get("source", ""),
            "call": r.get("family_call", "") or "ITPR",
            "in_alphafold": "",
        })
    return out


def census_targets(call: str = "ITPR") -> list[dict]:
    out, seen = [], set()
    for r in load_census():
        if r.get("call") != call:
            continue
        acc = (r.get("accession") or "").strip()
        if not acc or acc in seen:
            continue
        seen.add(acc)
        out.append({
            "accession": acc,
            "label": r.get("gene", "") or r.get("protein_name", "")[:40],
            "kind": accession_kind(acc),
            "group": normalise_group(r.get("group", "")),
            "paralog": r.get("gene", ""),
            "species": r.get("species", ""),
            "kingdom": r.get("kingdom", ""),
            "phylum": r.get("phylum", ""),
            "length": int(r.get("length") or 0),
            "source": r.get("source", ""),
            "call": r.get("call", ""),
            "in_alphafold": r.get("in_alphafold", ""),
        })
    return out


def probe_all(targets: list[dict], workers: int = 8,
              on_progress=None) -> list[dict]:
    """Probe every UniProt-shaped target, concurrently and cached."""
    rows: list[dict] = []
    todo = [t for t in targets if t["kind"] == "uniprot"]
    done = 0
    with cf.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(afdb_probe, t["accession"]): t for t in todo}
        for fut in cf.as_completed(futures):
            t = futures[fut]
            entry = fut.result()
            rows.append(_row(t, entry))
            done += 1
            if on_progress and done % 250 == 0:
                on_progress(done, len(todo))
    for t in targets:
        if t["kind"] != "uniprot":
            rows.append(_row(t, None))
    rows.sort(key=lambda r: r["accession"])
    return rows


def _row(t: dict, entry) -> dict:
    if entry is None:
        return {**t, "afdb_status": "not_applicable", "model_accession": "",
                "afdb_entry": "", "model_length": 0, "model_coverage": 0.0,
                "is_canonical": "", "afdb_version": "", "seq_version_date": "",
                "error": "accession is not a UniProt accession"}
    if not entry.found:
        return {**t, "afdb_status": "absent", "model_accession": "",
                "afdb_entry": "", "model_length": 0, "model_coverage": 0.0,
                "is_canonical": "", "afdb_version": "", "seq_version_date": "",
                "error": entry.error}
    cov = entry.coverage(t["length"])
    status = ("full" if entry.is_canonical and cov >= 0.95
              else "isoform" if not entry.is_canonical
              else "partial")
    return {**t, "afdb_status": status,
            "model_accession": entry.model_accession,
            "afdb_entry": entry.entry_id,
            "model_length": entry.model_length,
            "model_coverage": cov,
            "is_canonical": int(entry.is_canonical),
            "afdb_version": entry.version,
            "seq_version_date": entry.seq_version_date,
            "error": ""}


COLUMNS = ["accession", "label", "kind", "call", "group", "paralog", "species",
           "kingdom", "phylum", "length", "source", "in_alphafold",
           "afdb_status", "model_accession", "afdb_entry", "model_length",
           "model_coverage", "is_canonical", "afdb_version",
           "seq_version_date", "error"]


def summarise(rows: list[dict], scope: str) -> list[list]:
    """Coverage by group, with the three measurements side by side.

    `xref_agrees` is the point: UniProt's cross-reference against what the
    API serves. A disagreement is not a bookkeeping detail — it is the
    difference between the coverage a table claims and the coverage that
    exists.
    """
    out = []
    by_group: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_group[r.get("group") or "unassigned"].append(r)
    for group in sorted(by_group) + ["ALL"]:
        rs = rows if group == "ALL" else by_group[group]
        applicable = [r for r in rs if r["afdb_status"] != "not_applicable"]
        status = Counter(r["afdb_status"] for r in rs)
        usable = [r for r in applicable if r["model_coverage"] >= 0.95]
        xref = [r for r in applicable if r["in_alphafold"] in ("0", "1")]
        agree = sum(1 for r in xref
                    if (r["in_alphafold"] == "1") ==
                    (r["afdb_status"] != "absent"))
        out.append([
            scope, group, len(rs), len(applicable),
            status["full"], status["partial"], status["isoform"],
            status["absent"], status["not_applicable"],
            len(usable),
            round(len(usable) / len(applicable), 4) if applicable else 0.0,
            len(xref), agree,
            round(agree / len(xref), 4) if xref else "",
        ])
    return out


SUMMARY_COLUMNS = ["scope", "group", "n_records", "n_applicable", "full",
                   "partial", "isoform", "absent", "not_applicable",
                   "usable_ge_95pct", "usable_fraction", "n_with_xref",
                   "xref_agrees", "xref_agreement"]


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--set", dest="scope", choices=["reps", "census", "both"],
                    default="both")
    ap.add_argument("--call", default="ITPR",
                    help="census call to probe (default ITPR)")
    ap.add_argument("--workers", type=int, default=8)
    args = ap.parse_args(argv)

    out = results_dir()
    scopes = ["reps", "census"] if args.scope == "both" else [args.scope]
    summary: list[list] = []
    for scope in scopes:
        targets = rep_targets() if scope == "reps" else census_targets(args.call)
        n_up = sum(1 for t in targets if t["kind"] == "uniprot")
        print(f"[s11] {scope}: {len(targets)} targets, {n_up} UniProt-shaped")
        rows = probe_all(
            targets, workers=args.workers,
            on_progress=lambda d, n: print(f"[s11]   {scope} {d}/{n}", flush=True))
        write_tsv(out / f"afdb_coverage_{scope}.tsv", COLUMNS,
                  [[r.get(c, "") for c in COLUMNS] for r in rows])
        summary += summarise(rows, scope)
        tail = [r for r in rows if r["afdb_status"] != "not_applicable"]
        print(f"[s11] {scope}: "
              + ", ".join(f"{k}={v}" for k, v in
                          sorted(Counter(r["afdb_status"] for r in tail).items())))
    write_tsv(out / "afdb_coverage_summary.tsv", SUMMARY_COLUMNS, summary)
    print(f"[s11] wrote {out/'afdb_coverage_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
