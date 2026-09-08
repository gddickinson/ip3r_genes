"""S11 — every committed table, written once so the report and the figures
never re-derive anything (D13).

The one table that is not a dump of a previous step is
`tm_vs_reference.tsv`: **D14 asked structurally.** Every stage of this
project has had to separate ITPR from RyR by a positive test, and at this
stage the test is a fold comparison with a margin — the best TM-score
against the ITPR references, the best against the RyR references, and how
far apart they are. A structure that scores 0.71 against an IP3R means
nothing until it is also shown not to score 0.70 against a ryanodine
receptor.

The margin is **relative**, for the reason `s3_assign.py` gives: an
absolute gap would call every full-length structure and no fragment.
"""

from __future__ import annotations

import json
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import results_dir, sha256, write_tsv              # noqa: E402
from s11_tmalign import TM_FOLD, TM_RANDOM, tmalign_version     # noqa: E402

#: D7's number, reused. A structural family call needs the winner to beat
#: the runner-up by this fraction of its own score.
REL_MARGIN = 0.10

#: **The bar a family call needs before the margin is even consulted.**
#: A margin between two low scores is a margin between two non-matches: on
#: this panel all three negative controls beat their runner-up by 30 % of
#: their own score while scoring 0.17-0.27 against everything, so a rule
#: gated on the margin alone calls a dynein an IP3 receptor. The gate is
#: TM-align's own same-fold bar, applied to the winner, and it is what
#: makes the controls come back **declined** — which is the only thing
#: that makes the positives worth reading.
CALL_FLOOR = TM_FOLD

CANDIDATE_COLUMNS = [
    "entry_id", "entity_id", "call", "census_accession", "gene", "paralog",
    "method", "resolution", "released", "sample_length", "organism",
    "state", "ligands", "ligand_note", "pfams", "found_by_pfam", "chains",
    "n_chains", "deposited_residues", "eligible", "eligibility",
    "description", "title",
]
CONTROL_CANDIDATE_COLUMNS = [
    "entry_id", "entity_id", "control_class", "decoy_gene",
    "decoy_accession", "uniprot_length", "sample_length", "sample_fraction",
    "method", "resolution", "released", "organism", "state", "chains",
    "eligible", "eligibility", "title",
]
SELECTION_COLUMNS = [
    "role", "entry_id", "call", "paralog", "control_class", "decoy_gene",
    "decoy_accession", "organism", "resolution", "sample_length",
    "uniprot_length", "sample_fraction", "state", "ligand_note", "released",
    "selection_rule", "title",
]
SLOT_COLUMNS = ["slot", "outcome", "filled_with", "n_census_records",
                "stage_lost"]
UNFILLED_COLUMNS = ["control_class", "n_decoys", "stage_lost"]
VS_REF_COLUMNS = [
    "id", "role", "call", "paralog", "group", "source", "resolved_residues",
    "model_coverage", "best_itpr", "best_itpr_ref", "best_ryr",
    "best_ryr_ref", "tm_ceiling", "margin", "rel_margin", "structural_call",
    "call_agrees", "note",
]


def _rows(dicts: list[dict], columns: list[str]) -> list[list]:
    return [[d.get(c, "") for c in columns] for d in dicts]


def write_candidates(rows: list[dict]) -> Path:
    path = results_dir() / "reference_candidates.tsv"
    write_tsv(path, CANDIDATE_COLUMNS, _rows(rows, CANDIDATE_COLUMNS))
    return path


def write_control_candidates(rows: list[dict]) -> Path:
    path = results_dir() / "control_candidates.tsv"
    write_tsv(path, CONTROL_CANDIDATE_COLUMNS,
              _rows(rows, CONTROL_CANDIDATE_COLUMNS))
    return path


def write_selection(rows: list[dict]) -> Path:
    path = results_dir() / "reference_selection.tsv"
    write_tsv(path, SELECTION_COLUMNS, _rows(rows, SELECTION_COLUMNS))
    return path


def write_unfilled(rows: list[dict]) -> Path:
    path = results_dir() / "unfilled_controls.tsv"
    write_tsv(path, UNFILLED_COLUMNS, _rows(rows, UNFILLED_COLUMNS))
    return path


def write_slots(rows: list[dict]) -> Path:
    path = results_dir() / "model_slots.tsv"
    write_tsv(path, SLOT_COLUMNS, _rows(rows, SLOT_COLUMNS))
    return path


def write_manifest(rows: list[dict], columns: list[str]) -> Path:
    path = results_dir() / "structure_manifest.tsv"
    write_tsv(path, columns, _rows(rows, columns))
    return path


def write_tm_scores(rows: list[dict], columns: list[str]) -> Path:
    path = results_dir() / "tm_scores.tsv"
    write_tsv(path, columns, _rows(rows, columns))
    return path


def write_plddt(rows: list[dict], columns: list[str],
                summary: list[list], summary_columns: list[str]) -> list[Path]:
    a = results_dir() / "plddt_domains.tsv"
    b = results_dir() / "plddt_summary.tsv"
    write_tsv(a, columns, _rows(rows, columns))
    write_tsv(b, summary_columns, summary)
    return [a, b]


# --------------------------------------------------------------------------
# What a pair *is* — one classifier, used by the figure and by the report
# --------------------------------------------------------------------------

#: Reading order for the calibration classes, best-case first. The figure
#: and the report share this list so a class cannot appear in one and not
#: the other, and cannot be ordered differently in the two.
PAIR_CLASSES = ["same protein, different state", "IP3R vs IP3R",
                "IP3R vs RyR", "RyR vs RyR", "control vs anything"]


def pair_class(q: dict, t: dict) -> str:
    """Which calibration class a scored pair belongs to.

    Lives here rather than in the figure *and* the report, which is where
    it started: two copies classified identically but labelled differently,
    so a class silently vanished from one panel when a label was edited in
    the other.
    """
    roles = {q.get("role"), t.get("role")}
    if "control" in roles:
        return "control vs anything"
    if q.get("call") == "RYR" and t.get("call") == "RYR":
        return "RyR vs RyR"
    if {q.get("call"), t.get("call")} == {"ITPR", "RYR"}:
        return "IP3R vs RyR"
    if q.get("paralog") and q.get("paralog") == t.get("paralog"):
        return "same protein, different state"
    return "IP3R vs IP3R"


def pair_scores(pairs: list[dict], manifest: list[dict], which: str = "min"
                ) -> dict[str, list[float]]:
    """Every scored pair bucketed by `pair_class`.

    `which` selects the normalisation, and **both are reported** because on
    this panel they say different things. TM-align normalises by each
    input's length; `min` is the score against the *longer* chain and `max`
    against the shorter. An IP3 receptor subunit is ~2,200 resolved
    residues and a ryanodine receptor ~4,300, so a cross-family pair scores
    0.40 by `min` and 0.75 by `max` — the first says an IP3 receptor cannot
    account for a ryanodine receptor, the second that a ryanodine receptor
    accounts for an IP3 receptor. Reporting only one of those numbers would
    make a size relationship read as a fold result.
    """
    man = {m["id"]: m for m in manifest}
    pick = min if which == "min" else max
    out: dict[str, list[float]] = {}
    for p in pairs:
        if p["query"] == p["target"] or not int(p.get("ok") or 0):
            continue
        q, t = man.get(p["query"]), man.get(p["target"])
        if not q or not t:
            continue
        try:
            score = pick(float(p["tm_query"]), float(p["tm_target"]))
        except (TypeError, ValueError):
            continue
        out.setdefault(pair_class(q, t), []).append(score)
    return out


# --------------------------------------------------------------------------
# D14, structurally
# --------------------------------------------------------------------------

def _best_against(pairs: dict[tuple[str, str], dict], ident: str,
                  refs: list[dict]) -> tuple[float, str]:
    """Best TM-score of `ident` against a set of references.

    Normalised **by the reference**, not by the query. That choice is the
    whole difference between a fold statement and a size statement: a
    1,100-residue fragment aligned into a 2,700-residue receptor scores
    high normalised by its own length whatever it is, because it only has
    to explain itself. Normalised by the reference it has to explain the
    receptor, which is the claim being made.
    """
    best, who = 0.0, ""
    for ref in refs:
        rid = ref["id"]
        if rid == ident:
            continue
        row = pairs.get((ident, rid)) or pairs.get((rid, ident))
        if not row or not row.get("ok"):
            continue
        score = (row["tm_query"] if row["query"] == rid else row["tm_target"])
        if score > best:
            best, who = score, rid
    return round(best, 4), who


def tm_ceiling(query_residues: int, reference_residues: int) -> float:
    """The highest reference-normalised TM-score a structure of this size
    could possibly reach.

    TM-score normalised by the reference sums over *aligned* pairs and
    divides by the reference's length, so a structure with `n` residues
    scored against an `L`-residue reference is bounded above by `n / L`
    — before any question of similarity. On this panel the partial
    non-vertebrate models are 1,000-1,300 residues against ~2,050-residue
    references, which puts several of them within a few hundredths of the
    same-fold bar by arithmetic alone. Printing the ceiling is what stops
    the report saying "could not reach the bar however good it is", which
    is true of some of those rows and not of others.
    """
    if not query_residues or not reference_residues:
        return 0.0
    return round(min(1.0, query_residues / reference_residues), 4)


def vs_reference(manifest: list[dict], tm_rows: list[dict]) -> list[dict]:
    pairs = {(r["query"], r["target"]): r for r in tm_rows}
    residues = {m["id"]: int(m.get("resolved_residues") or 0)
                for m in manifest}
    refs = [m for m in manifest
            if m.get("role") in ("reference", "state_panel")
            and m.get("status") == "ok"]
    itpr_refs = [r for r in refs if r.get("call") == "ITPR"]
    ryr_refs = [r for r in refs if r.get("call") == "RYR"]
    out = []
    for m in manifest:
        if m.get("status") != "ok":
            continue
        b_i, who_i = _best_against(pairs, m["id"], itpr_refs)
        b_r, who_r = _best_against(pairs, m["id"], ryr_refs)
        top = max(b_i, b_r)
        margin = round(abs(b_i - b_r), 4)
        rel = round(margin / top, 4) if top else 0.0
        if top < TM_RANDOM:
            call = "no_fold_match"
        elif top < CALL_FLOOR:
            call = "below_fold_bar"
        elif rel < REL_MARGIN:
            call = "no_call"
        else:
            call = "ITPR" if b_i > b_r else "RYR"
        expected = m.get("call", "")
        # `call_agrees` is blank, not 0, where the instrument declined to
        # call at all. A refusal is not a disagreement, and scoring it as
        # one would turn "this model is too partial to place" into "shape
        # contradicts the census" (S7's `unconstrained`, S23's
        # `uncontrolled`).
        agrees = ("" if expected not in ("ITPR", "RYR")
                  or call not in ("ITPR", "RYR")
                  else int(call == expected))
        winner_ref = who_i if b_i >= b_r else who_r
        out.append({
            "id": m["id"], "role": m.get("role", ""), "call": expected,
            "paralog": m.get("paralog", ""), "group": m.get("group", ""),
            "source": m.get("source", ""),
            "resolved_residues": m.get("resolved_residues", ""),
            "model_coverage": m.get("model_coverage", ""),
            "best_itpr": b_i, "best_itpr_ref": who_i,
            "best_ryr": b_r, "best_ryr_ref": who_r,
            "tm_ceiling": tm_ceiling(residues.get(m["id"], 0),
                                     residues.get(winner_ref, 0)),
            "margin": margin, "rel_margin": rel, "structural_call": call,
            "call_agrees": agrees,
            "note": (
                "below the random-structure bar" if call == "no_fold_match"
                else f"top score {top} is below the {CALL_FLOOR} same-fold "
                     "bar — no family call"
                if call == "below_fold_bar"
                else f"within {REL_MARGIN:.0%} of both references"
                if call == "no_call" else ""),
        })
    return out


def write_vs_reference(rows: list[dict]) -> Path:
    path = results_dir() / "tm_vs_reference.tsv"
    write_tsv(path, VS_REF_COLUMNS, _rows(rows, VS_REF_COLUMNS))
    return path


# --------------------------------------------------------------------------
# The run's own record
# --------------------------------------------------------------------------

def write_stats(paths: list[Path], parameters: dict,
                self_tests: dict) -> Path:
    stats = {
        "task": "S11",
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tmalign_version": tmalign_version(),
        "parameters": {
            "rel_margin": REL_MARGIN,
            "call_floor": CALL_FLOOR,
            "tm_random_bar": TM_RANDOM,
            "tm_fold_bar": TM_FOLD,
            **parameters,
        },
        "self_tests": self_tests,
        "tables": {p.name: {"sha256": sha256(p), "bytes": p.stat().st_size}
                   for p in sorted(paths) if p.exists()},
    }
    path = results_dir() / "structures_stats.json"
    path.write_text(json.dumps(stats, indent=2, sort_keys=True) + "\n")
    return path
