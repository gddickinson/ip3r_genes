"""S11 — which structures are the references, and which are the controls.

Seven rules, written out and enforced rather than described. Every one is a
**positive test**, and none of them reads an entry title to decide what
family a structure belongs to.

  **R1 — family by this project's own call, never by the title.** An
  entity's UniProt accession must be a census v6 row whose `call` is the
  family in question. The titles here would be a trap: "ryanodine receptor"
  and "inositol trisphosphate receptor" entries share every diagnostic
  Pfam, and D14 has required a positive test at every previous stage.
  Entities with no census row are reported, not guessed at; on the current
  candidate set there are none carrying a family signature at length, so
  the profile fallback is never reached — which is the check, not an
  assumption.

  **R2 — cryo-EM with a recorded resolution.** Every X-ray ITPR entity in
  the candidate set is a truncated construct: they run 226–2,217 aa at
  1.9–7.4 Å against full-length cryo-EM entities of 2,633–2,771 aa, and the
  two longest (the "large cytosolic domain", 2,217 aa) are the two lowest
  resolution structures in the whole family — 7.3 and 7.4 Å. So "cryo-EM"
  is a measurement about what exists, not a preference. The report computes
  those ranges rather than quoting them.

  **R3 — full length.** The entity's sample sequence must sit in the
  family's declared size band (`family.MIN_LENGTH_AA` … `MAX_LENGTH_AA`)
  for an ITPR, and above `RYR_MIN_LENGTH_AA` for a RyR — the band is
  defined to *exclude* RyRs, so a RyR reference cannot be selected by the
  ITPR band. A binding-core construct would make every TM-score a
  statement about one domain.

  **R4 — the state must be recorded.** A state is parsed from the title
  against a stated vocabulary and corroborated against the deposited
  ligands. An entry whose state cannot be read stays a candidate and is
  never a reference.

  **R5 — one primary reference per paralog**, ranked by resolution, ties
  broken on resolved residues then entry id, so the pick is deterministic
  (D24).

  **R6 — the state panel.** For the paralog offering the most distinct
  states, one reference per state. Without it a predicted model is scored
  against a single conformation and the fold question is confounded with
  the gating question — the panel measures how far TM-score moves between
  two structures of the *same* protein.

  **R7 — negative controls come from S1's committed decoy panel**, one per
  decoy class, and are held to R2 and R3 as well: a control must be a
  full-length experimental structure or it is not a control, it is a
  fragment scoring low for reasons of size. A class with no qualifying
  structure is reported as unfilled rather than quietly dropped
  (`s5_build_baits.py`'s `unfilled_slots.tsv` rule).

The tetramer is reduced to one chain before anything is scored; see
`s11_struct_io.largest_chain`.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_rcsb as rcsb                                       # noqa: E402
from s11_lib import PROJECT_ROOT, load_census                 # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT))
from src.utils import family                                  # noqa: E402

#: A RyR is defined by being *above* the family band, which is where the
#: band's upper bound was put on purpose (`src/utils/family.py`).
RYR_MIN_LENGTH_AA = 4_000

#: R7's source — S1's decoy panel, with the class each decoy was chosen for.
DECOY_PANEL = (PROJECT_ROOT / "results" / "benchmark_controls" /
               "panel_decoys.json")

#: The sister family is a *reference* here, not a control: the brief asks
#: for the cryo-EM RyR entry, and D14 makes it the sharpest thing an ITPR
#: model can be scored against. S1's RyR decoys are therefore excluded
#: from the control classes.
EXCLUDED_DECOY_CLASSES = ("RyR (sister family)",)

#: R4's vocabulary. Ordered: the first match wins, so the compound states
#: are tested before the words they contain ("inactive-like" before
#: "active", "preactivated" before "activated").
STATE_PATTERNS: list[tuple[str, str]] = [
    ("inactive-like", r"inactive[- ]like"),
    ("labile resting", r"labile resting"),
    ("higher-order inhibited", r"higher[- ]order inhibited"),
    ("preactivated", r"pre[- ]?activ"),
    ("activated", r"\bactivat|\bactive\b"),
    ("primed", r"\bprimed\b"),
    ("resting", r"\bresting\b"),
    ("apo", r"\bapo\b|ligand[- ]free|\bligand free\b"),
    ("open", r"\bopen\b"),
    ("closed", r"\bclosed\b"),
    ("inactivated", r"\binactivat"),
]

#: Ligand chemical-component ids that corroborate a state claim. These are
#: what the depositor modelled, so they are evidence; the title is not.
LIGAND_MEANING = {"I3P": "IP3", "CA": "Ca2+", "ATP": "ATP",
                  "JYP": "adenophostin A", "CFF": "caffeine",
                  "RYA": "ryanodine"}


def parse_state(title: str) -> str:
    low = (title or "").lower()
    for name, pattern in STATE_PATTERNS:
        if re.search(pattern, low):
            return name
    return ""


def ligand_note(ligands: str) -> str:
    got = [LIGAND_MEANING[l] for l in ligands.split(";")
           if l in LIGAND_MEANING]
    return "+".join(got)


# --------------------------------------------------------------------------
# R1 — the family call
# --------------------------------------------------------------------------

def census_index() -> dict[str, dict]:
    return {r["accession"]: r for r in load_census()}


def call_entity(row: dict, census: dict[str, dict]) -> tuple[str, str, str]:
    """(call, accession, gene) for one polymer entity, by R1.

    Returns `no_census_row` when the entity maps to UniProt but the census
    has never seen that accession, and `no_uniprot` when RCSB records no
    UniProt mapping at all. Both are reported; neither is guessed at.
    """
    accs = [a for a in (row.get("uniprot") or "").split(";") if a]
    if not accs:
        return "no_uniprot", "", ""
    for acc in accs:
        rec = census.get(acc)
        if rec:
            return rec.get("call", ""), acc, rec.get("gene", "")
    return "no_census_row", accs[0], ""


def paralog_of(gene: str, call: str) -> str:
    """The paralog label, upper-cased so rat `Itpr1` and human `ITPR1` are
    one group.

    A gene the census leaves blank gets **no** label rather than falling
    back to the family name. Such a row is a real receptor structure and
    stays in the candidate table, but it is not a paralog reference: a
    group called `RYR` sitting beside RYR1 and RYR2 in a figure would
    assert a paralog that nothing named.
    """
    return (gene or "").upper()


# --------------------------------------------------------------------------
# R2/R3 — eligibility
# --------------------------------------------------------------------------

def length_ok(call: str, sample_length: int) -> bool:
    if call == "ITPR":
        return family.MIN_LENGTH_AA <= sample_length <= family.MAX_LENGTH_AA
    if call == "RYR":
        return sample_length >= RYR_MIN_LENGTH_AA
    return False


def eligibility(row: dict) -> tuple[bool, str]:
    """R2 + R3 + R4 as a single positive test, naming the rule that failed."""
    if row.get("call") not in ("ITPR", "RYR"):
        return False, "R1 not called family"
    if row.get("method") != "ELECTRON MICROSCOPY":
        return False, "R2 not cryo-EM"
    if not row.get("resolution"):
        return False, "R2 no resolution"
    if not length_ok(row["call"], int(row.get("sample_length") or 0)):
        return False, "R3 not full length"
    if not row.get("state"):
        return False, "R4 state not readable"
    if not row.get("paralog"):
        return False, "R5 paralog not named by the census"
    return True, "eligible"


# --------------------------------------------------------------------------
# R5/R6 — selection
# --------------------------------------------------------------------------

def rank_key(row: dict) -> tuple:
    """Lexicographic, so the audit can name the component that decided it."""
    return (float(row.get("resolution") or 99.0),
            -int(row.get("deposited_residues") or 0),
            row.get("entry_id", ""))


def build_candidates() -> list[dict]:
    """Every family-Pfam RCSB entity, flattened and called (R1), with the
    eligibility verdict on each row. This is the committed denominator."""
    ids, found_by = rcsb.family_candidate_entries()
    entries = rcsb.fetch_entries(ids)
    census = census_index()
    rows: list[dict] = []
    for entry_id in sorted(entries):
        for row in rcsb.entity_rows(entries[entry_id]):
            call, acc, gene = call_entity(row, census)
            row["call"] = call
            row["census_accession"] = acc
            row["gene"] = gene
            row["paralog"] = paralog_of(gene, call) if call in ("ITPR", "RYR") else ""
            row["state"] = parse_state(row.get("title", ""))
            row["ligand_note"] = ligand_note(row.get("ligands", ""))
            row["found_by_pfam"] = ";".join(found_by.get(entry_id, []))
            ok, reason = eligibility(row)
            row["eligible"] = int(ok)
            row["eligibility"] = reason
            rows.append(row)
    return rows


def select_primary(rows: list[dict]) -> list[dict]:
    """R5 — the best entry per paralog, both families."""
    best: dict[str, dict] = {}
    for r in rows:
        if not r["eligible"]:
            continue
        key = r["paralog"]
        if key not in best or rank_key(r) < rank_key(best[key]):
            best[key] = r
    out = []
    for key in sorted(best, key=lambda k: (best[k]["call"], k)):
        sel = dict(best[key])
        sel["role"] = "reference"
        sel["selection_rule"] = "R5 best resolution for paralog"
        out.append(sel)
    return out


def select_state_panel(rows: list[dict], family_call: str = "ITPR") -> list[dict]:
    """R6 — one entry per state, for the paralog offering the most states.

    The paralog is chosen by measurement (most distinct eligible states,
    ties on best resolution then name), so the conformational control is
    taken where the data actually supports one.
    """
    by_paralog: dict[str, list[dict]] = {}
    for r in rows:
        if r["eligible"] and r["call"] == family_call:
            by_paralog.setdefault(r["paralog"], []).append(r)
    if not by_paralog:
        return []
    def score(p: str) -> tuple:
        rs = by_paralog[p]
        return (-len({x["state"] for x in rs}),
                min(float(x["resolution"]) for x in rs), p)
    chosen = min(by_paralog, key=score)
    best: dict[str, dict] = {}
    for r in by_paralog[chosen]:
        st = r["state"]
        if st not in best or rank_key(r) < rank_key(best[st]):
            best[st] = r
    out = []
    for st in sorted(best):
        sel = dict(best[st])
        sel["role"] = "state_panel"
        sel["selection_rule"] = f"R6 best resolution for {chosen} state '{st}'"
        out.append(sel)
    return out


# --------------------------------------------------------------------------
# R7 — negative controls from S1's decoy panel
# --------------------------------------------------------------------------

def decoy_classes() -> dict[str, list[dict]]:
    """S1's decoy panel grouped by the class it was chosen for.

    The class is read back out of the committed description S1 wrote
    (`DECOY [<class>] <name>`) rather than re-declared here, so S11's
    controls cannot drift from the panel S1 validated the scorer on.
    """
    if not DECOY_PANEL.exists():
        return {}
    out: dict[str, list[dict]] = {}
    for rec in json.loads(DECOY_PANEL.read_text()):
        m = re.search(r"DECOY \[(.*?)\]", rec.get("description") or "")
        cls = m.group(1) if m else ""
        if not cls or cls in EXCLUDED_DECOY_CLASSES:
            continue
        out.setdefault(cls, []).append(rec)
    return out


def control_candidates(min_fraction: float = 0.80) -> list[dict]:
    """Every decoy's experimental entries, flattened and screened by R2/R3.

    R3 for a control is "full length" relative to *its own* UniProt length,
    not to the family band: a control is a different protein, and the point
    of the size match is that the chain being scored is the whole thing.
    """
    rows: list[dict] = []
    for cls, recs in sorted(decoy_classes().items()):
        for rec in recs:
            acc = rec["accession"]
            up_len = int(rec.get("length_aa") or 0)
            ids = rcsb.entries_for_uniprot(acc)
            if not ids:
                continue
            for entry_id, entry in sorted(rcsb.fetch_entries(ids).items()):
                for row in rcsb.entity_rows(entry):
                    if acc not in (row.get("uniprot") or ""):
                        continue
                    sample = int(row.get("sample_length") or 0)
                    frac = (sample / up_len) if up_len else 0.0
                    row.update({
                        "control_class": cls,
                        "decoy_gene": rec.get("gene_symbol", ""),
                        "decoy_accession": acc,
                        "uniprot_length": up_len,
                        "sample_fraction": round(frac, 3),
                        "call": "control",
                        "state": parse_state(row.get("title", "")),
                    })
                    ok = (row.get("method") == "ELECTRON MICROSCOPY"
                          and row.get("resolution")
                          and frac >= min_fraction)
                    row["eligible"] = int(bool(ok))
                    row["eligibility"] = (
                        "eligible" if ok else
                        "R2 not cryo-EM" if row.get("method") != "ELECTRON MICROSCOPY"
                        else "R2 no resolution" if not row.get("resolution")
                        else f"R3 fragment ({frac:.2f} of UniProt length)")
                    rows.append(row)
    return rows


def select_controls(rows: list[dict], target_length: int) -> tuple[list[dict], list[dict]]:
    """R7 — one control per class, size-matched to the family reference.

    "Closest in length to the reference chain" and not "largest": TM-score
    is length-normalised but the score two *unrelated* structures can reach
    still depends on their size ratio, so a size-matched control is what
    makes the measured floor mean anything.

    Returns (selected, unfilled) — a class with no qualifying structure is
    named with the stage that lost it, never silently dropped.
    """
    by_class: dict[str, list[dict]] = {}
    for r in rows:
        if r["eligible"]:
            by_class.setdefault(r["control_class"], []).append(r)
    selected = []
    for cls in sorted(by_class):
        best = min(by_class[cls], key=lambda r: (
            abs(int(r["sample_length"]) - target_length),
            float(r.get("resolution") or 99.0), r["entry_id"]))
        sel = dict(best)
        sel["role"] = "control"
        sel["selection_rule"] = (
            f"R7 closest to {target_length} aa in class '{cls}'")
        selected.append(sel)
    unfilled = []
    for cls, recs in sorted(decoy_classes().items()):
        if cls in by_class:
            continue
        had_any = any(r["control_class"] == cls for r in rows)
        unfilled.append({
            "control_class": cls,
            "n_decoys": len(recs),
            "stage_lost": "no experimental entry" if not had_any
                          else "no full-length cryo-EM entry",
        })
    return selected, unfilled
