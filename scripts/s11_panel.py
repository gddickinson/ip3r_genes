"""S11 step 3 — assemble the structure panel and write its manifest.

Four rules, and the fourth is the one that keeps the table honest.

  **P1 — every reference and control selected in `s11_refs`.** Reduced to
  one chain first (`s11_struct_io.largest_chain`): both families are
  homotetramers and an IP3R assembly is ~11,000 residues, so scoring a
  monomer model against the deposited assembly would answer a question
  about quaternary structure.

  **P2 — every S6 representative AFDB holds a model for**, at whatever
  coverage. Coverage is a *column*, never a filter, because the reason a
  representative has no usable model is itself the result S11 reports.

  **P3 — group fill.** Where an S6 group (or a vertebrate paralog) has no
  representative modelled at `FULL_COVERAGE`, the best-covered census ITPR
  record in that group stands in, and `fill_rule` records that the slot was
  filled rather than met. Without this the panel is six proteins.

  **P4 — nothing enters the panel unread.** Every file is parsed and its
  largest chain measured; `model_coverage` in the manifest is *resolved
  residues over the length the census holds*, not what the API advertised.
  A chain below `MIN_CHAIN_RESIDUES` is rejected with the reason recorded.
  AFDB's advertised span and the residues actually in the file are two
  different numbers, and only the second one is scored.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_refs as refs                                        # noqa: E402
from s11_lib import (afdb_fetch, afdb_probe, load_tsv, pdb_fetch,  # noqa: E402
                     results_dir, sha256, structures_dir)
from s11_struct_io import largest_chain, read_structure, write_ca_pdb  # noqa: E402

#: A model covering at least this much of the protein is "full"; below it
#: the row is kept and reported, never silently used as if it were whole.
FULL_COVERAGE = 0.95

#: Under this many resolved residues a chain is not a channel subunit and
#: cannot carry a fold claim about one.
MIN_CHAIN_RESIDUES = 300

MANIFEST_COLUMNS = [
    "id", "role", "call", "paralog", "group", "species", "organism",
    "source", "source_id", "chain", "resolved_residues", "expected_length",
    "model_coverage", "resolution", "state", "ligands", "method",
    "mean_plddt", "mean_bfactor", "fill_rule", "path", "sha256", "status",
    "note",
]


def panel_dir() -> Path:
    d = structures_dir() / "panel"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _mean_bfactor(chain) -> float:
    if not chain.residues:
        return 0.0
    return round(sum(r.bfactor for r in chain.residues) / len(chain.residues), 2)


def _install(src: Path, ident: str, expected_length: int,
             predicted: bool) -> dict:
    """P4 — parse, reduce to one chain, write the CA trace, measure it.

    `predicted` decides what the B-factor column *means*. AFDB writes
    pLDDT there; a cryo-EM deposition writes an atomic displacement
    parameter, and on this panel those run to 122. Reporting either as the
    other would put a confidence score of 121.88 in a table — so the mean
    is written to `mean_plddt` only for a predicted model, and to
    `mean_bfactor` otherwise.
    """
    try:
        chains = read_structure(src)
    except Exception as exc:                                    # noqa: BLE001
        return {"status": "unreadable", "note": str(exc)[:120]}
    if not chains:
        return {"status": "unreadable", "note": "no CA atoms parsed"}
    chain = largest_chain(chains)
    if len(chain) < MIN_CHAIN_RESIDUES:
        return {"status": "too_short", "chain": chain.chain_id,
                "resolved_residues": len(chain),
                "note": f"{len(chain)} residues < {MIN_CHAIN_RESIDUES}"}
    dest = panel_dir() / f"{ident}.pdb"
    write_ca_pdb(chain, dest)
    cov = round(len(chain) / expected_length, 4) if expected_length else 0.0
    mean = _mean_bfactor(chain)
    return {"status": "ok", "chain": chain.chain_id,
            "resolved_residues": len(chain), "model_coverage": cov,
            "mean_plddt": mean if predicted else "",
            "mean_bfactor": "" if predicted else mean,
            "path": str(dest), "sha256": sha256(dest), "note": ""}


# --------------------------------------------------------------------------
# P1 — references and controls
# --------------------------------------------------------------------------

def install_experimental(selected: list[dict]) -> list[dict]:
    rows = []
    for sel in selected:
        entry_id = sel["entry_id"]
        role = sel["role"]
        ident = f"{role}_{sel.get('paralog') or sel.get('decoy_gene', '')}_{entry_id}"
        ident = ident.replace(" ", "_")
        try:
            cif = pdb_fetch(entry_id)
        except Exception as exc:                                # noqa: BLE001
            rows.append({"id": ident, "role": role, "source": "PDB",
                         "source_id": entry_id, "status": "fetch_failed",
                         "note": str(exc)[:120]})
            continue
        expected = int(sel.get("sample_length") or 0)
        measured = _install(cif, ident, expected, predicted=False)
        rows.append({
            "id": ident, "role": role,
            "call": sel.get("call", ""),
            "paralog": sel.get("paralog", "") or sel.get("control_class", ""),
            "group": sel.get("control_class", "") or "reference",
            "species": sel.get("organism", ""),
            "organism": sel.get("organism", ""),
            "source": "PDB", "source_id": entry_id,
            "expected_length": expected,
            "resolution": sel.get("resolution", ""),
            "state": sel.get("state", ""),
            "ligands": sel.get("ligand_note", "") or sel.get("ligands", ""),
            "method": sel.get("method", ""),
            "fill_rule": sel.get("selection_rule", ""),
            **measured,
        })
    return rows


# --------------------------------------------------------------------------
# P2/P3 — predicted models
# --------------------------------------------------------------------------

def _model_rows(coverage_tsv: Path) -> list[dict]:
    path = results_dir() / coverage_tsv
    return load_tsv(path) if path.exists() else []


def choose_models() -> tuple[list[dict], list[dict]]:
    """P2 + P3 — which predicted models the panel takes, and why each.

    Returns (chosen, slot_audit). The audit names every group slot and
    whether it was met by a representative, filled from the census, or
    left empty — `s5_build_baits.py`'s unfilled-slot rule.
    """
    all_reps = _model_rows(Path("afdb_coverage_reps.tsv"))
    all_census = _model_rows(Path("afdb_coverage_census.tsv"))
    has_model = ("absent", "not_applicable")
    reps = [r for r in all_reps if r.get("afdb_status") not in has_model]
    census = [r for r in all_census if r.get("afdb_status") not in has_model]

    def cov(r: dict) -> float:
        try:
            return float(r.get("model_coverage") or 0)
        except ValueError:
            return 0.0

    def census_slot(r: dict) -> str:
        """The slot a census row belongs to.

        Vertebrates are split by paralog, everything else by taxonomic
        group. The paralog question is a vertebrate question (2R), and
        labelling a protist row with a paralog would assert what S7 was run
        to test — `s6_rep_spec.group_of`'s rule, applied here.
        """
        group = (r.get("group") or "").strip()
        para = (r.get("paralog") or "").upper()
        if group == "Vertebrata" and para.startswith("ITPR"):
            return para[:5]
        return group or "unassigned"

    # One vocabulary, and it is the census's. S6's representative table
    # names its groups differently (`plant`, `protist`, `invert_metazoa`),
    # and S6's `protist` is four census groups at once — so a slot audit
    # built from both vocabularies at face value reports `plant` unfilled
    # beside `Viridiplantae` filled, which is one group counted twice under
    # two names (the `Viridiplantae`/`viridiplantae` problem S23 left, in a
    # different key). Every representative is a census accession by
    # construction, so its slot is looked up rather than translated.
    slot_by_acc = {r["accession"]: census_slot(r) for r in all_census}

    def slot_of(r: dict) -> str:
        got = slot_by_acc.get(r.get("accession", ""))
        if got:
            return got
        return census_slot(r)

    chosen: dict[str, dict] = {}
    for r in reps:
        r["fill_rule"] = "P2 S6 representative"
        chosen[r["accession"]] = r

    met = {slot_of(r) for r in reps if cov(r) >= FULL_COVERAGE}
    # The slot universe comes from **every** census and representative row,
    # not only the ones that have a model. Built from the modelled rows
    # alone, a group AFDB holds nothing for would vanish from the audit
    # instead of appearing in it as unfilled — which is the one outcome the
    # audit exists to record.
    slots = sorted({slot_of(r) for r in all_reps if r.get("call") != "RYR"}
                   | {slot_of(r) for r in all_census})
    audit = []
    for slot in slots:
        if slot in met:
            audit.append({"slot": slot, "outcome": "met by representative",
                          "filled_with": "", "n_census_records":
                          sum(1 for r in all_census if slot_of(r) == slot),
                          "stage_lost": ""})
            continue
        pool = [r for r in census if slot_of(r) == slot
                and cov(r) >= FULL_COVERAGE
                and r["accession"] not in chosen]
        if not pool:
            in_slot = [r for r in census if slot_of(r) == slot]
            best = max(in_slot, key=cov, default=None)
            n_records = sum(1 for r in all_census if slot_of(r) == slot)
            audit.append({
                "slot": slot, "outcome": "unfilled", "filled_with": "",
                "n_census_records": n_records,
                "stage_lost": (
                    "no AFDB model at any coverage" if not in_slot else
                    f"no model at or above {FULL_COVERAGE:.0%} coverage "
                    f"(best {cov(best):.2f})")})
            continue
        pick = max(pool, key=lambda r: (cov(r), int(r.get("length") or 0),
                                        r["accession"]))
        pick["fill_rule"] = f"P3 census fill for slot '{slot}'"
        chosen[pick["accession"]] = pick
        audit.append({"slot": slot, "outcome": "filled from census",
                      "filled_with": pick["accession"], "n_census_records":
                      sum(1 for r in all_census if slot_of(r) == slot),
                      "stage_lost": ""})
    return sorted(chosen.values(), key=lambda r: r["accession"]), audit


def install_models(chosen: list[dict]) -> list[dict]:
    rows = []
    for r in chosen:
        acc = r["accession"]
        ident = f"model_{acc}"
        entry = afdb_probe(acc)
        expected = int(r.get("length") or 0)
        try:
            src = afdb_fetch(entry, structures_dir() / "afdb" / f"{acc}.pdb")
        except Exception as exc:                                # noqa: BLE001
            src = None
            note = str(exc)[:120]
        if src is None:
            rows.append({"id": ident, "role": "model", "source": "AFDB",
                         "source_id": entry.entry_id, "status": "fetch_failed",
                         "call": r.get("call", "") or "ITPR",
                         "paralog": r.get("paralog", ""),
                         "group": r.get("group", ""),
                         "species": r.get("species", ""),
                         "expected_length": expected,
                         "note": locals().get("note", "no pdb url")})
            continue
        measured = _install(src, ident, expected, predicted=True)
        rows.append({
            "id": ident, "role": "model", "call": r.get("call", "") or "ITPR",
            "paralog": r.get("paralog", ""), "group": r.get("group", ""),
            "species": r.get("species", ""), "organism": r.get("species", ""),
            "source": "AFDB", "source_id": entry.entry_id,
            "expected_length": expected, "resolution": "", "state": "predicted",
            "ligands": "", "method": "PREDICTION",
            "fill_rule": r.get("fill_rule", ""),
            **measured,
        })
    return rows


def build_panel() -> tuple[list[dict], list[dict], list[dict], list[dict],
                           list[dict], list[dict]]:
    """Everything step 3 produces, in one pass.

    Returns (manifest, rcsb_candidates, control_candidates, unfilled_control
    classes, model slot audit, selection) — the selection separately from
    the manifest because a selected reference that then failed to download
    still has to appear as *selected*, or the table would show the rule
    picking whatever survived rather than what it picked.
    """
    candidates = refs.build_candidates()
    primary = refs.select_primary(candidates)
    states = refs.select_state_panel(candidates)
    itpr_lengths = sorted(int(p["sample_length"]) for p in primary
                          if p["call"] == "ITPR")
    target = (itpr_lengths[len(itpr_lengths) // 2] if itpr_lengths else 2700)
    ctrl_cands = refs.control_candidates()
    controls, unfilled = refs.select_controls(ctrl_cands, target)
    # A state-panel entry that is also the primary reference is one file.
    seen = {p["entry_id"] for p in primary}
    states = [s for s in states if s["entry_id"] not in seen]
    selection = primary + states + controls
    manifest = install_experimental(selection)
    chosen, slot_audit = choose_models()
    manifest += install_models(chosen)
    return (manifest, candidates, ctrl_cands, unfilled, slot_audit, selection)
