"""S3 step 4 — merge the HMM sweep into census v3.

Census v3 is census v2 with a **second, independent verdict on every row**
and the sweep's novel hits appended. Two instruments now score the family:

  * the architecture rule of D14b (S2) — which Pfam signatures a record
    carries, and which RyR-specific ones it does not;
  * best-profile assignment with a bit-score margin (s3_assign.py) — which
    of `itpr.hmm` / `ryr.hmm` scores the sequence higher, and by how much.

They see different evidence. The architecture rule reads InterPro's
annotation of a record; the profiles read the residues. So the merged call
is stated as a rule over both, and every row says which instruments spoke:

  agree            → the call, confidence high
  one speaks only  → that call, at that instrument's confidence
  disagree         → `conflict`, kept and reported, never silently resolved
  neither          → `unassigned`

That third case is the one worth watching, and the second is the one that
does work: S2 left 2,181 records `unassigned` because a partial architecture
cannot be called by a rule that reads absence as evidence, and a profile has
no such problem.

Outputs (results/census_v3/):
  census_v3.tsv             every record, both verdicts, the merged call
  novel_hits.tsv            sweep hits absent from census v2
  call_changes.tsv          rows whose call differs from census v2
  conflicts.tsv             the two instruments disagreeing
  missed_by_hmm.tsv         v2 ITPR records in swept proteomes with no hit
  proteomes_without_hits.tsv  swept proteomes with no ITPR call (S5 leads)
  convergence.tsv           per-round jackhmmer counts + the D10 verdict
  census_v3_stats.json

Run:  python3 scripts/s3_census_v3.py
"""

from __future__ import annotations

import json
import re
import sys
import time
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_assign import ASSIGN_FIELDS  # noqa: E402
from scripts.s3_hmm_lib import (  # noqa: E402
    CENSUS_V3_DIR, HMM_SWEEP_DIR, acc_key, load_census_v2, read_tsv,
    write_tsv,
)
from scripts.s3_kill import accepted_targets, evaluate  # noqa: E402
from scripts.s3_seed_spec import JACKHMMER_SEEDS  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

# A gene symbol that names this family or its sister, in any of UniProt's
# casings, with the numeric suffixes UniProt appends when a proteome holds
# several entries for one gene name (`Itpr1_0`, `Ryr3_2`).
FAMILY_NAME_RE = re.compile(r"^(itpr|itr-?\d*|ryr|unc-68)", re.IGNORECASE)

V2_FIELDS = [
    "accession", "gene", "protein_name", "species", "taxon_id", "group",
    "kingdom", "phylum", "class", "order", "family", "length", "in_band",
    "fragment", "reviewed", "protein_existence", "in_alphafold",
    "seed_pfams", "pfams", "n_itpr_arch", "ryr_pfams",
]
V3_FIELDS = (["accession", "call", "confidence", "reason", "instruments",
              "arch_call", "profile_call", "profile_confidence",
              "itpr_score", "ryr_score", "rel_margin", "searches", "source"]
             + V2_FIELDS[1:])


def log(msg: str) -> None:
    print(f"[census_v3] {msg}", flush=True)


# ------------------------------------------------------------ the merge rule
def merge_call(arch: str, profile: str, profile_conf: str,
               arch_conf: str) -> tuple[str, str, str, str]:
    """(call, confidence, reason, instruments) from the two verdicts."""
    arch_spoke = arch in ("ITPR", "RYR")
    prof_spoke = profile in ("ITPR", "RYR")
    if arch_spoke and prof_spoke:
        if arch == profile:
            return (arch, "high",
                    f"architecture and profile both call {arch}", "both")
        return ("conflict", "none",
                f"architecture says {arch}, profile says {profile}", "both")
    if arch_spoke:
        return (arch, arch_conf or "medium",
                "architecture only — no profile verdict", "architecture")
    if prof_spoke:
        return (profile, profile_conf,
                "profile only — architecture incomplete", "profile")
    return ("unassigned", "none",
            "neither instrument calls it", "neither")


# ------------------------------------------------------------------- inputs
def load_assignments(path: Path) -> dict[str, dict]:
    if not path.exists():
        raise SystemExit(f"missing {path} — run the sweep first")
    return {acc_key(r["accession"]): r for r in read_tsv(path)}


def load_convergence() -> tuple[list[dict], dict]:
    """Per-round jackhmmer rows + the D10 verdict per seed tag.

    The kill criterion needs to know which family each included target
    belongs to; that comes from the sweep's own profile assignments, so the
    criterion is evaluated on the same evidence the census is built from.
    """
    raw_dir = require_data_root() / "hmmer"
    sweep = load_assignments(HMM_SWEEP_DIR / "hmmsearch_assignments.tsv")

    def call_of(target_name: str) -> str:
        """jackhmmer target names are full UniProt names (sp|ACC|ID)."""
        parts = target_name.split("|")
        acc = acc_key(parts[1] if len(parts) >= 3 else target_name)
        return sweep.get(acc, {}).get("assignment", "unassigned")


    def evidence_of(target_name: str) -> str:
        parts = target_name.split("|")
        acc = acc_key(parts[1] if len(parts) >= 3 else target_name)
        return sweep.get(acc, {}).get("evidence", "no sweep hit")

    rows: list[dict] = []
    verdicts: dict[str, dict] = {}
    composition: list[dict] = []
    for tag in sorted(JACKHMMER_SEEDS):
        path = raw_dir / f"jackhmmer_{tag}_rounds.json"
        if not path.exists():
            log(f"  no jackhmmer run for {tag} — skipped")
            continue
        data = json.loads(path.read_text())
        rounds = data["rounds"]
        by_name = {t: call_of(t) for rd in rounds for t in rd.get("included", [])}
        verdict = evaluate(rounds, by_name, own_family="ITPR",
                           converged=bool(data["converged"]))
        verdict.update({"tag": tag, "accession": data["accession"],
                        "converged": data["converged"],
                        "n_rounds": len(rounds)})
        verdict["accepted_targets"] = len(
            accepted_targets(rounds, verdict["accepted_rounds"]))
        # What an unconverged run is actually accreting. K1 watches the
        # sister family only, on purpose — the "neither family" bucket also
        # holds real fragmentary members the D22 span gate declined — so
        # the composition of the final model is reported instead of being
        # folded into a pass/fail.
        final = rounds[-1].get("included", []) if rounds else []
        comp = Counter((call_of(t), evidence_of(t)) for t in final)
        for (call, ev), n in sorted(comp.items(), key=lambda kv: -kv[1]):
            composition.append({
                "seed_tag": tag, "round": rounds[-1]["round"] if rounds else 0,
                "profile_call": call, "sweep_evidence": ev, "targets": n,
                "share": round(n / len(final), 4) if final else 0.0})
        verdicts[tag] = verdict
        for pr in verdict["per_round"]:
            rows.append({"seed_tag": tag, "accession": data["accession"],
                         **pr, "converged": data["converged"],
                         "verdict": verdict["verdict"],
                         "kill_rule": verdict["rule"]})
        log(f"  {tag}: {len(rounds)} rounds, converged={data['converged']}, "
            f"D10 {verdict['verdict']} ({verdict['reason']})")
    if composition:
        write_tsv(CENSUS_V3_DIR / "jackhmmer_model_composition.tsv",
                  ["seed_tag", "round", "profile_call", "sweep_evidence",
                   "targets", "share"], composition)
    return rows, verdicts


# -------------------------------------------------------------------- build
def build(v2_rows: list[dict], calib: dict[str, dict],
          sweep: dict[str, dict]) -> tuple[list[dict], list[dict]]:
    """census v3 rows + the novel rows, from v2 plus the two assignment sets."""
    out: list[dict] = []
    for r in v2_rows:
        acc = acc_key(r["accession"])
        # A record's profile verdict comes from the calibration run (which
        # scored the whole seeded space); the sweep confirms it for the
        # subset that is also in a reference proteome.
        p = calib.get(acc) or sweep.get(acc) or {}
        searches = [s for s, d in (("interpro", True), ("calibration",
                    acc in calib), ("sweep", acc in sweep)) if d]
        call, conf, why, instruments = merge_call(
            r["call"], p.get("assignment", ""), p.get("confidence", ""),
            r["confidence"])
        out.append({
            **{k: r.get(k, "") for k in V2_FIELDS},
            "call": call, "confidence": conf, "reason": why,
            "instruments": instruments,
            "arch_call": r["call"], "profile_call": p.get("assignment", ""),
            "profile_confidence": p.get("confidence", ""),
            "itpr_score": p.get("itpr_score", ""),
            "ryr_score": p.get("ryr_score", ""),
            "rel_margin": p.get("rel_margin", ""),
            "searches": "+".join(searches), "source": "InterPro",
        })

    v2_accs = {acc_key(r["accession"]) for r in v2_rows}
    novel: list[dict] = []
    module_only: list[dict] = []
    for acc, p in sorted(sweep.items()):
        if acc in v2_accs:
            continue
        # D22: a hit that spans only a shared module is not a candidate
        # family member, it is a protein carrying one of the family's
        # domains. Those are recorded and counted, not enrolled — otherwise
        # the census fills with troponins and SPRY proteins.
        if p["evidence"] == "module":
            module_only.append(p)
            continue
        call = p["assignment"]
        novel.append({
            "accession": acc, "gene": p["gene"],
            "protein_name": p["protein_name"], "species": p["species"],
            "taxon_id": p["taxon_id"],
            # The sweep DB is Vertebrata by construction; the lineage
            # columns InterPro filled for v2 rows are not available here,
            # so only what the sweep actually knows is written.
            "group": "Vertebrata", "kingdom": "Metazoa",
            "phylum": "Chordata", "class": "", "order": "", "family": "",
            "length": p["length"], "in_band": "", "fragment": "",
            "reviewed": p["reviewed"], "protein_existence": "",
            "in_alphafold": "", "seed_pfams": "", "pfams": "",
            "n_itpr_arch": "", "ryr_pfams": "",
            "call": call, "confidence": p["confidence"],
            "reason": f"profile only — not in the InterPro census; {p['reason']}",
            "instruments": "profile", "arch_call": "",
            "profile_call": call, "profile_confidence": p["confidence"],
            "itpr_score": p["itpr_score"], "ryr_score": p["ryr_score"],
            "rel_margin": p["rel_margin"],
            "searches": "sweep", "source": "HMM",
        })
    return out + novel, novel, module_only


def missed_by_hmm(v3: list[dict], sweep: dict[str, dict],
                  swept_taxa: set[int], db_path: Path) -> list[dict]:
    """v2 ITPR records from a swept proteome that the sweep did not hit.

    A real miss would break the completeness argument, so it is measured
    rather than assumed — and the measurement has to separate two very
    different things:

      * the record is **not in the sweep database**. A UniProt reference
        proteome is one canonical protein per gene; the InterPro census
        enumerated all of UniProtKB, which holds far more entries for the
        same taxa. An accession absent from the DB was never searched, so
        it is a statement about the database's scope, not about the
        profiles' sensitivity.
      * the record **is in the DB and the sweep did not find it**. That is
        a sensitivity failure and the only kind that bears on completeness.

    Telling them apart needs one streaming pass over the 8.8 GB DB, which
    is why it is done here once for the whole candidate set rather than
    guessed from gene names.
    """
    hit_accs = set(sweep)
    hit_gene_species = {(r["gene"].upper(), r["species"])
                        for r in sweep.values() if r["gene"]}
    candidates = []
    for r in v3:
        if r["source"] != "InterPro" or r["arch_call"] != "ITPR":
            continue
        if not r["taxon_id"] or int(r["taxon_id"]) not in swept_taxa:
            continue
        if r["accession"] in hit_accs:
            continue
        candidates.append(r)
    if not candidates:
        return []

    wanted = {r["accession"] for r in candidates}
    in_db: set[str] = set()
    log(f"  checking {len(wanted)} missed accessions against the sweep DB…")
    with db_path.open() as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].split(None, 1)[0].split("|")
                if len(parts) >= 3 and acc_key(parts[1]) in wanted:
                    in_db.add(acc_key(parts[1]))
    log(f"  {len(in_db)} of them are in the sweep DB")

    out = []
    for r in candidates:
        if r["accession"] in in_db:
            verdict = "IN THE DB AND NOT FOUND — a sensitivity failure"
        elif (r["gene"].upper(), r["species"].split(" (")[0]) in \
                hit_gene_species:
            verdict = "not in the DB; another entry of the same gene was hit"
        else:
            verdict = "not in the sweep DB — never searched"
        out.append({**r, "verdict": verdict})
    return out


def main() -> int:
    t0 = time.time()
    CENSUS_V3_DIR.mkdir(parents=True, exist_ok=True)
    v2_rows = load_census_v2()
    calib = load_assignments(HMM_SWEEP_DIR / "calibration_assignments.tsv")
    sweep = load_assignments(HMM_SWEEP_DIR / "hmmsearch_assignments.tsv")
    log(f"v2 {len(v2_rows)} rows | calibration {len(calib)} | sweep {len(sweep)}")

    v3, novel, module_only = build(v2_rows, calib, sweep)
    write_tsv(CENSUS_V3_DIR / "census_v3.tsv", V3_FIELDS, v3)
    write_tsv(CENSUS_V3_DIR / "novel_hits.tsv", V3_FIELDS, novel)
    # The full declined list is not written out: every one of these rows is
    # already in `hmm_sweep/hmmsearch_assignments.tsv`, recoverable with
    # `evidence == "module"`, and duplicating 12,624 rows of troponins to
    # say "not this family" is not worth 3 MB of repository.
    #
    # A module-only hit whose own gene symbol is a family name is not a
    # SPRY protein that happens to score — it is a piece of a split or
    # truncated gene model of a real family gene, sitting below the D22
    # gate because there is not enough of it left to span a family domain.
    # Those are an annotation-quality lead (S18/S10), so they are separated
    # out rather than left in the declined pile.
    fragments = [r for r in module_only
                 if FAMILY_NAME_RE.match(r["gene"] or "")]
    write_tsv(CENSUS_V3_DIR / "subthreshold_family_fragments.tsv",
              ASSIGN_FIELDS, fragments)
    log(f"declined {len(module_only)} sweep hits as module-only matches "
        f"(D22); {len(fragments)} of them carry a family gene symbol and "
        "are kept as annotation leads")

    v2_call = {acc_key(r["accession"]): r["call"] for r in v2_rows}
    changes = [r for r in v3 if r["source"] == "InterPro"
               and r["call"] != v2_call.get(r["accession"])]
    conflicts = [r for r in v3 if r["call"] == "conflict"]
    write_tsv(CENSUS_V3_DIR / "call_changes.tsv",
              V3_FIELDS + ["v2_call"],
              [{**r, "v2_call": v2_call.get(r["accession"], "")}
               for r in changes])
    write_tsv(CENSUS_V3_DIR / "conflicts.tsv", V3_FIELDS, conflicts)

    manifest = read_tsv(HMM_SWEEP_DIR / "proteome_manifest.tsv")
    swept_taxa = {int(m["Organism Id"]) for m in manifest}
    db_path = Path(json.loads(
        (HMM_SWEEP_DIR / "proteome_db_stats.json").read_text())["db_path"])
    missed = missed_by_hmm(v3, sweep, swept_taxa, db_path)
    write_tsv(CENSUS_V3_DIR / "missed_by_hmm.tsv", V3_FIELDS + ["verdict"],
              missed)

    itpr_taxa = {int(r["taxon_id"]) for r in v3
                 if r["call"] == "ITPR" and r["taxon_id"]}
    no_hits = [m for m in manifest
               if int(m["Organism Id"]) not in itpr_taxa]
    write_tsv(CENSUS_V3_DIR / "proteomes_without_hits.tsv",
              list(manifest[0].keys()), no_hits)

    conv_rows, verdicts = load_convergence()
    for tag, v in verdicts.items():
        pr = v["per_round"]
        if len(pr) >= 2:
            v["offfamily_growth"] = round(
                pr[-1]["n_offfamily"] / pr[0]["n_offfamily"], 2) \
                if pr[0]["n_offfamily"] else None
    if conv_rows:
        write_tsv(CENSUS_V3_DIR / "convergence.tsv",
                  ["seed_tag", "accession", "round", "new_targets",
                   "n_included", "n_own", "n_sister", "n_offfamily",
                   "sister_frac", "sister_rise", "growth", "converged",
                   "verdict", "kill_rule"], conv_rows)

    counts = Counter(r["call"] for r in v3)
    inst = Counter(r["instruments"] for r in v3)
    stats = {
        "v2_rows": len(v2_rows), "novel_rows": len(novel),
        "v3_rows": len(v3),
        "calls": dict(counts), "instruments": dict(inst),
        "call_changes": len(changes),
        "changes_by_transition": dict(Counter(
            f"{v2_call.get(r['accession'], '')} → {r['call']}"
            for r in changes)),
        "conflicts": len(conflicts),
        "missed_by_hmm": {
            "total": len(missed),
            "by_verdict": dict(Counter(m["verdict"] for m in missed)),
            "sensitivity_failures": sum(
                1 for m in missed if m["verdict"].startswith("IN THE DB")),
        },
        "swept_proteomes": len(manifest),
        "proteomes_without_itpr": len(no_hits),
        "novel_calls": dict(Counter(r["call"] for r in novel)),
        "novel_evidence": dict(Counter(
            sweep[r["accession"]]["evidence"] for r in novel)),
        "module_only_declined": len(module_only),
        "subthreshold_family_fragments": len(fragments),
        "subthreshold_fragment_genes": dict(Counter(
            r["gene"] for r in fragments).most_common(12)),
        "module_only_top_genes": dict(Counter(
            r["gene"] for r in module_only if r["gene"]).most_common(12)),
        "jackhmmer": {t: {k: v for k, v in d.items() if k != "per_round"}
                      for t, d in verdicts.items()},
        "elapsed_s": round(time.time() - t0, 1),
    }
    (CENSUS_V3_DIR / "census_v3_stats.json").write_text(
        json.dumps(stats, indent=1))
    # Dashboard live panel (session protocol): the merge is the last step,
    # so this is where S3 reports itself finished.
    (PROJECT_ROOT / "results" / "session_live.json").write_text(json.dumps({
        "task": "S3",
        "steps": [
            {"label": "profiles built (itpr.hmm + ryr.hmm)", "done": True},
            {"label": "calibration vs the S2 architecture call", "done": True},
            {"label": "hmmsearch sweep, 763 proteomes", "done": True},
            {"label": f"jackhmmer ({len(verdicts)}/3 seeds parsed)",
             "done": len(verdicts) == 3},
            {"label": f"census v3 ({stats['v3_rows']:,} records)",
             "done": True},
        ],
    }, indent=1))
    log(f"census v3: {len(v3)} rows ({len(novel)} novel), calls {dict(counts)}")
    log(f"call changes vs v2: {len(changes)}; conflicts: {len(conflicts)}")
    log(f"done in {stats['elapsed_s']}s → {CENSUS_V3_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
