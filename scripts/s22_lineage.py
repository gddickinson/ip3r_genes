"""S22 stage 7 — the comparative question, and what it could have seen.

The brief: in lineages where the upstream PLC/IP3 pathway is reduced or
absent, is the binding core relaxed relative to the pore?  Two halves.

**The taxon list is derived, not read.**  Stage 6 measured the PI-PLC
repertoire of every eukaryotic reference proteome this project swept.  A
lineage qualifies when its proteome carries an ITPR and no PI-PLC.  The
list is committed whatever its length, because a list of length zero is an
answer to the question and a list of length four is a statement about
power.

**The statistic is the paired difference, not the core.**  A distant tip is
divergent everywhere, so its core identity alone says nothing; what the
hypothesis predicts is that the *gap* between core and pore widens.  Each
tip therefore contributes one number, core identity minus pore identity,
measured on the same alignment columns for every tip.

**The positive control is the ryanodine receptors.**  They are in the same
alignment, they have the pore, and they do not bind IP3 — which is the
condition the lineage test is looking for, arrived at by evolution rather
than by a missing enzyme.  Without it a null result cannot be told from a
test that could never fire, and the power calculation is anchored on the
effect the control actually produces rather than on a number chosen to make
the study look adequate.

Alignment caveat, stated rather than buried: msa_v2 spans four kingdoms and
the distant tips sit at 20-25 % identity to the human reference, where the
placement of an individual core column is not certain.  Coverage per tip
per module is a column of the panel, and a tip resolving less than half of
either module is dropped and counted.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402

MIN_MODULE_COVERAGE = 0.50
REF_TIP = "ITPR3_Homo_sapiens_Human_Q14573"
S20_PRESENCE = L.S20_DIR / "proteome_presence.tsv"
#: The groups the PI-PLC sweep covered.  Imported rather than retyped so the
#: co-occurrence table can never name a group that was not searched.
from s22_plc import EUK_GROUPS as PLC_GROUPS                    # noqa: E402
S3_ASSIGN = L.PROJECT_ROOT / "results" / "hmm_sweep" / "hmmsearch_assignments.tsv"
VERT_MANIFEST = L.PROJECT_ROOT / "results" / "hmm_sweep" / "proteome_manifest.tsv"


# ---------------------------------------------------------------------------
# the pathway table
# ---------------------------------------------------------------------------

def cooccurrence() -> tuple[list[dict], list[dict]]:
    """ITPR presence against PI-PLC presence, per proteome and per group.

    Vertebrate ITPR presence is keyed on **taxid**, because S3's sweep
    committed assignments and not an accession-to-proteome attribution; 5 of
    763 vertebrate proteomes share a taxid with another, and the key is
    recorded in the table so the looseness is visible rather than implied.
    The four non-vertebrate groups use S20's measured per-proteome call.
    """
    plc = {(r["group"], r["upid"]): r
           for r in L.read_tsv(L.OUT_DIR / "plc_repertoire.tsv")}
    itpr_taxids = {r["taxon_id"] for r in L.read_tsv(S3_ASSIGN)
                   if r["assignment"] == "ITPR"}
    per: list[dict] = []
    for r in L.read_tsv(S20_PRESENCE):
        # Only the groups the PI-PLC sweep covered may appear.  Archaea and
        # the bacterial sample were not searched, and a row reading "0 carry
        # a PI-PLC" for a domain that demonstrably has phospholipases would
        # be a claim this task never made.
        if r["group"] not in PLC_GROUPS:
            continue
        key = (r["group"], r["upid"])
        p = plc.get(key)
        per.append({"group": r["group"], "upid": r["upid"],
                    "organism": r["organism"], "taxid": r["taxid"],
                    "itpr_status": r["itpr_status"],
                    "itpr_key": "proteome",
                    "n_itpr": r["n_itpr"],
                    "n_pi_plc": p["n_pi_plc"] if p else "",
                    "plc_status": p["plc_status"] if p else "not_searched"})
    for r in L.read_tsv(VERT_MANIFEST):
        upid, taxid = r["Proteome Id"], r["Organism Id"]
        p = plc.get(("vertebrata", upid))
        per.append({"group": "vertebrata", "upid": upid,
                    "organism": r["Organism"], "taxid": taxid,
                    "itpr_status": ("present" if taxid in itpr_taxids else "absent"),
                    "itpr_key": "taxid",
                    "n_itpr": "",
                    "n_pi_plc": p["n_pi_plc"] if p else "",
                    "plc_status": p["plc_status"] if p else "not_searched"})
    groups: dict[str, dict] = {}
    for r in per:
        g = groups.setdefault(r["group"], {
            "group": r["group"], "n_proteomes": 0, "n_itpr": 0, "n_plc": 0,
            "n_both": 0, "n_itpr_no_plc": 0, "n_plc_no_itpr": 0,
            "n_neither": 0, "itpr_key": r["itpr_key"]})
        i = r["itpr_status"] == "present"
        p = r["plc_status"] == "present"
        g["n_proteomes"] += 1
        g["n_itpr"] += i
        g["n_plc"] += p
        g["n_both"] += (i and p)
        g["n_itpr_no_plc"] += (i and not p)
        g["n_plc_no_itpr"] += (p and not i)
        g["n_neither"] += (not i and not p)
    return per, [groups[g] for g in sorted(groups)]


# ---------------------------------------------------------------------------
# the panel
# ---------------------------------------------------------------------------

def _msa_columns(mods: list[dict]) -> tuple[list[int], list[int]]:
    """msa_v2 columns of the two primary modules, from the ITPR3 reference."""
    rows = {r["resi"]: r for r in L.constraint(L.STRUCTURE_REF)}
    core = M.residues(mods, L.STRUCTURE_REF, M.PRIMARY["ligand_core"])
    pore = M.residues(mods, L.STRUCTURE_REF, M.PRIMARY["pore_module"])
    cc = sorted(int(rows[i]["msa_col"]) for i in core
                if i in rows and rows[i].get("msa_col") not in ("", None))
    pc = sorted(int(rows[i]["msa_col"]) for i in pore
                if i in rows and rows[i].get("msa_col") not in ("", None))
    return cc, pc


def _plc_index() -> tuple[dict[str, dict], dict[str, dict]]:
    """PLC status by taxid and by lowercase binomial."""
    by_tax: dict[str, dict] = {}
    by_name: dict[str, dict] = {}
    for r in L.read_tsv(L.OUT_DIR / "plc_repertoire.tsv"):
        by_tax.setdefault(r["taxid"], r)
        binom = " ".join(r["organism"].split()[:2]).lower().strip("(),")
        by_name.setdefault(binom, r)
    return by_tax, by_name


def panel(mods: list[dict]) -> list[dict]:
    aln, reps = L.msa_v2()
    by_label = {r["label"]: r for r in reps}
    ref = aln[REF_TIP]
    cc, pc = _msa_columns(mods)
    by_tax, by_name = _plc_index()
    rows: list[dict] = []
    for label, seq in aln.items():
        if label == REF_TIP:
            continue
        r = by_label[label]
        rec = {"tip": label, "group": r["group"], "paralog": r["paralog"],
               "species": r["species"], "taxid": r["taxon_id"],
               "kingdom": r["kingdom"], "phylum": r["phylum"]}
        for name, cols in (("core", cc), ("pore", pc)):
            n_ref = sum(1 for c in cols if ref[c] != "-")
            cov = [c for c in cols if ref[c] != "-" and seq[c] != "-"]
            ident = sum(1 for c in cov if seq[c] == ref[c])
            rec[f"{name}_n_ref"] = n_ref
            rec[f"{name}_n_covered"] = len(cov)
            rec[f"{name}_coverage"] = L.fmt(len(cov) / n_ref if n_ref else 0)
            rec[f"{name}_identity"] = L.fmt(ident / len(cov) if cov else None)
        ok = (rec["core_identity"] != "" and rec["pore_identity"] != ""
              and float(rec["core_coverage"]) >= MIN_MODULE_COVERAGE
              and float(rec["pore_coverage"]) >= MIN_MODULE_COVERAGE)
        rec["delta_core_minus_pore"] = (
            L.fmt(float(rec["core_identity"]) - float(rec["pore_identity"]))
            if ok else "")
        rec["in_test"] = ok
        rec["drop_reason"] = "" if ok else "module coverage below 0.50"

        hit = by_tax.get(r["taxon_id"]) if r["taxon_id"] else None
        if hit is None:
            hit = by_name.get(" ".join(r["species"].split()[:2]).lower())
        rec["plc_upid"] = hit["upid"] if hit else ""
        rec["n_pi_plc"] = hit["n_pi_plc"] if hit else ""
        rec["plc_status"] = hit["plc_status"] if hit else "no_reference_proteome"
        if r["group"] == "RYR":
            rec["stratum"] = "ryr_control"
        elif r["kingdom"] == "Metazoa" and r["group"] in ("ITPR1", "ITPR2",
                                                          "ITPR3",
                                                          "vertebrate_basal"):
            rec["stratum"] = "vertebrate_itpr"
        else:
            rec["stratum"] = f"nonvert_itpr_plc_{rec['plc_status']}"
        rows.append(rec)
    return rows


# ---------------------------------------------------------------------------
# tests and power
# ---------------------------------------------------------------------------

def _deltas(rows: list[dict], stratum: str) -> list[float]:
    return [float(r["delta_core_minus_pore"]) for r in rows
            if r["in_test"] and r["stratum"] == stratum]


def tests(rows: list[dict]) -> list[dict]:
    strata = sorted({r["stratum"] for r in rows})
    out: list[dict] = []
    itpr_all = [float(r["delta_core_minus_pore"]) for r in rows
                if r["in_test"] and r["stratum"] != "ryr_control"]
    comparisons = [("ryr_control", "every ITPR tip", itpr_all,
                    "positive control — the pore without the ligand")]
    absent = _deltas(rows, "nonvert_itpr_plc_absent")
    present = _deltas(rows, "nonvert_itpr_plc_present")
    comparisons.append(("nonvert_itpr_plc_absent",
                        "nonvert_itpr_plc_present", present,
                        "the brief's lineage test"))
    for stratum, against, base, note in comparisons:
        a = _deltas(rows, stratum)
        mw = L.mann_whitney(a, base)
        out.append({
            "test": f"{stratum} vs {against}", "note": note,
            "n_test": len(a), "n_reference": len(base),
            "median_test_delta": L.fmt(L.median(a)),
            "median_reference_delta": L.fmt(L.median(base)),
            "shift": L.fmt((L.median(a) or 0) - (L.median(base) or 0))
                     if a and base else "",
            "cles": L.fmt(mw["cles"]), "p_mannwhitney": L.fmt(mw["p"], 6),
        })
    for s in strata:
        d = _deltas(rows, s)
        out.append({
            "test": f"{s} (description)", "note": "no comparison — the stratum itself",
            "n_test": len(d), "n_reference": "",
            "median_test_delta": L.fmt(L.median(d)),
            "median_reference_delta": "", "shift": "", "cles": "",
            "p_mannwhitney": "",
        })
    return out


def power(rows: list[dict], iters: int = 4000, seed: int = 20220422) -> list[dict]:
    """What shift the lineage test could have detected, by simulation.

    The reference distribution is the PLC-present non-vertebrate tips; a
    synthetic test group of the observed size is drawn from it with a shift
    added and the same Mann-Whitney run.  The smallest shift reaching 80 %
    rejection is the minimum detectable effect, and it is reported beside
    the shift the RyR control actually produces — which is the largest
    effect this alignment is known to be able to show.
    """
    ref = _deltas(rows, "nonvert_itpr_plc_present")
    n_abs = len(_deltas(rows, "nonvert_itpr_plc_absent"))
    ryr = _deltas(rows, "ryr_control")
    itpr = [float(r["delta_core_minus_pore"]) for r in rows
            if r["in_test"] and r["stratum"] != "ryr_control"]
    control_shift = ((L.median(ryr) or 0) - (L.median(itpr) or 0)) if ryr else None
    out: list[dict] = []
    if n_abs == 0 or len(ref) < 5:
        out.append({"n_test_group": n_abs, "n_reference": len(ref),
                    "shift": "", "power": "",
                    "control_shift_ryr": L.fmt(control_shift),
                    "note": "no test group, or reference too small — "
                            "the test could not be run at any effect size"})
        return out
    rng = random.Random(seed)
    for shift in (0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.50):
        rej = 0
        for _ in range(iters):
            a = [rng.choice(ref) + shift for _ in range(n_abs)]
            b = [rng.choice(ref) for _ in range(len(ref))]
            p = L.mann_whitney(a, b)["p"]
            if p is not None and p < 0.05:
                rej += 1
        out.append({"n_test_group": n_abs, "n_reference": len(ref),
                    "shift": L.fmt(shift), "power": L.fmt(rej / iters),
                    "control_shift_ryr": L.fmt(control_shift),
                    "note": ""})
    return out


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    steps = [("co-occurrence", False), ("panel", False), ("tests", False),
             ("power", False)]
    L.live("lineage", steps)
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")

    per, by_group = cooccurrence()
    L.write_tsv(L.OUT_DIR / "pathway_by_proteome.tsv", per, list(per[0]))
    L.write_tsv(L.OUT_DIR / "pathway_cooccurrence.tsv", by_group,
                list(by_group[0]))
    taxa = [r for r in per if r["itpr_status"] == "present"
            and r["plc_status"] == "absent"]
    L.write_tsv(L.OUT_DIR / "plc_absent_taxa.tsv", taxa,
                list(per[0]) if taxa else
                ["group", "upid", "organism", "taxid", "itpr_status",
                 "itpr_key", "n_itpr", "n_pi_plc", "plc_status"])
    L.log(f"ITPR-carrying proteomes with no PI-PLC: {len(taxa)}")
    steps[0] = ("co-occurrence", True); L.live("lineage", steps)

    rows = panel(mods)
    L.write_tsv(L.OUT_DIR / "lineage_panel.tsv", rows, list(rows[0]))
    steps[1] = ("panel", True); L.live("lineage", steps)

    tr = tests(rows)
    L.write_tsv(L.OUT_DIR / "lineage_test.tsv", tr, list(tr[0]))
    steps[2] = ("tests", True); L.live("lineage", steps)

    pw = power(rows)
    L.write_tsv(L.OUT_DIR / "lineage_power.tsv", pw, list(pw[0]))
    steps[3] = ("power", True); L.live("lineage", steps)

    for r in tr:
        if r["p_mannwhitney"]:
            L.log(f"{r['test']}: n={r['n_test']} shift={r['shift']} "
                  f"p={r['p_mannwhitney']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
