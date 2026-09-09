"""S19 — what each search channel actually contributed.

Four questions, four denominators, kept apart on purpose because conflating
them is how a methods section overstates a method.

1. **Records.** Per channel: how many records it returned, how many *only* it
   returned, and how the census accumulated in the order the channels ran.
2. **Head to head, inside one database.** The InterPro/Pfam enumeration and
   the profile-HMM sweep are comparable only in the sequence space both could
   see, so the comparison is restricted to the accessions actually present in
   the swept reference-proteome FASTAs (`s19_lib.build_universe`). Without
   that restriction "the enumeration holds a record the profile HMM did not
   return" cannot be told apart from "that record was never in the database
   the profile HMM searched". Split by length band, because the question is
   *where* in the length distribution a family profile earns its place over
   domain annotation.
3. **Genes.** The genome sweep is the ground truth for where the genes are,
   so "what would a proteome-only search have missed?" is asked per genome x
   cell, not per record. Three ways to miss a gene: the species has no
   reference proteome at all, it has one and the gene is not in it, or it is
   in it and no protein record resolves to that paralog.
4. **Cost.** Recorded seconds per channel against the records its census step
   added — read from each session's own committed run statistics, and left
   blank rather than read as zero where a session did not record one.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402
import s19_recovery as R                                        # noqa: E402

#: The census versions in run order, with the channel each one added.
CENSUS_STEPS = [
    ("v1", "targeted database search (the app: NCBI/Ensembl/UniProt/Compara)"),
    ("v2", "InterPro/Pfam exhaustive enumeration"),
    ("v3", "profile HMM + jackhmmer, vertebrate reference proteomes"),
    ("v4", "genome sweep (miniprot, 309 vertebrate assemblies)"),
    ("v5", "profile HMM + jackhmmer, 6,928 other reference proteomes"),
    ("v6", "genome sweep (miniprot, 194 non-vertebrate assemblies)"),
]

LENGTH_BANDS = [
    ("<1000 aa (profile tail)", 0, 1000),
    ("1000-1999 aa (partial)", 1000, 2000),
    (">=2000 aa (gene-scale)", 2000, 10 ** 9),
]


# ------------------------------------------------------------------ 1. records

def record_sets(log=S.log) -> tuple[dict, list[list]]:
    sets = S.method_accessions()
    n_app = sets.pop("app_search_n", 0)
    rows = []
    for m in S.METHOD_ORDER:
        mine = sets.get(m, set())
        others = set().union(*[s for k, s in sets.items() if k != m])
        n = len(mine) if mine else (n_app if m == "app_search" else 0)
        rows.append([m, S.METHOD_LABEL[m], S.METHOD_UNIVERSE[m], n,
                     len(mine - others) if mine else "",
                     "" if mine else "membership not retained; count only"])
    S.write_tsv(S.out_dir() / "method_sets.tsv",
                ["method", "label", "universe", "n_records", "n_unique",
                 "note"], rows)
    log(f"method_sets: {len(rows)} channels")
    return sets, rows


def census_growth(log=S.log) -> list[list]:
    """Census size at each version, and the records the step added."""
    rows = []
    prev: set[str] = set()
    v1 = S.read_json(S.V1_DELTA)
    rows.append(["v1", CENSUS_STEPS[0][1], v1.get("v1_accessions", 0),
                 v1.get("v1_accessions", 0), 0,
                 "v1 is the app's four search bundles; S2 verdicted every "
                 "accession of it that v2 does not hold"])
    prev = set()
    for version, channel in CENSUS_STEPS[1:]:
        accs = {S.acc_key(r["accession"]) for r in S.census_rows(version)}
        note = ""
        if version == "v2":
            note = (f"{v1.get('shared', 0)} of v1's "
                    f"{v1.get('v1_accessions', 0)} carried forward; "
                    f"{v1.get('only_v1', 0)} verdicted as outside the "
                    "seeded space")
        rows.append([version, channel, len(accs), len(accs - prev),
                     len(prev - accs), note])
        prev = accs
    S.write_tsv(S.out_dir() / "census_growth.tsv",
                ["census", "channel_added", "n_records", "n_added",
                 "n_dropped", "note"], rows)
    log(f"census_growth: v1 {rows[0][2]} -> v6 {rows[-1][2]} records")
    return rows


def pairwise(sets: dict, log=S.log) -> list[list]:
    keys = [m for m in S.METHOD_ORDER if sets.get(m)]
    rows = []
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            inter = sets[a] & sets[b]
            union = sets[a] | sets[b]
            rows.append([a, b, len(sets[a]), len(sets[b]), len(inter),
                         round(len(inter) / len(union), 4) if union else 0.0,
                         int(S.METHOD_UNIVERSE[a] == S.METHOD_UNIVERSE[b])])
    S.write_tsv(S.out_dir() / "method_pairwise.tsv",
                ["method_a", "method_b", "n_a", "n_b", "n_shared", "jaccard",
                 "same_universe"], rows)
    log(f"method_pairwise: {len(rows)} pairs")
    return rows


# -------------------------------------------------- 2. head to head, one DB

def _record_lengths() -> dict[str, int]:
    out: dict[str, int] = {}
    for version in ("v2", "v6"):
        for r in S.census_rows(version):
            acc = S.acc_key(r["accession"])
            if acc not in out:
                out[acc] = int(S.fnum(r.get("length"), float, 0))
    for path in [S.HMM_DIR / "hmmsearch_assignments.tsv"] + \
            [S.S20_DIR / f"assignments_{g}.tsv" for g in S.NONVERT_GROUPS]:
        for r in S.read_tsv(path):
            acc = S.acc_key(r["accession"])
            out.setdefault(acc, int(S.fnum(r.get("length"), float, 0)))
    return out


def _hmm_called(group: str) -> tuple[set[str], set[str]]:
    """(targets scored, targets the profile pair actually called family).

    The distinction is not cosmetic. `hmmsearch` returns every protein either
    profile scores above the E-value, which in the vertebrate database is
    18,501 targets of which 13,371 are `unassigned` — SPRY-domain proteins
    the RyR profile reaches and D22's 200-match-state gate declines. Scoring
    the enumeration against the raw target list would credit the profile
    sweep with 13,000 records it never claimed. Both are carried, and the
    family-called set is what the comparison is made on.
    """
    path = (S.HMM_DIR / "hmmsearch_assignments.tsv" if group == "vertebrata"
            else S.S20_DIR / f"assignments_{group}.tsv")
    scored, called = set(), set()
    for r in S.read_tsv(path):
        acc = S.acc_key(r["accession"])
        scored.add(acc)
        if r.get("assignment") in ("ITPR", "RYR"):
            called.add(acc)
    return scored, called


def head_to_head(sets: dict, log=S.log) -> tuple[list[list], dict]:
    """Pfam enumeration vs profile HMM inside the searched sequence space.

    Everything outside the swept FASTAs is excluded from **both** sides, so a
    miss here is a genuine failure to return a record that was in the
    database rather than an accounting artefact of two different search
    spaces. Both channels are counted by what they **call family** — the
    enumeration by carrying a family signature, the sweep by clearing D22's
    gate — because a comparison of a curated set against a raw hit list is
    not a comparison of two methods. The jackhmmer column is what iteration
    added over one pass in the same database, counting only the rounds D10
    accepted.
    """
    pfam = sets["pfam_enumeration"]
    length = _record_lengths()
    rows: list[list] = []
    summary: dict = {}
    scopes = [("vertebrata", "profile_hmm_vert", "jackhmmer_vert")]
    scopes += [(g, "profile_hmm_nonvert", "jackhmmer_nonvert")
               for g in S.NONVERT_GROUPS]
    for group, hmm_key, _jack_key in scopes:
        scored, hmm = _hmm_called(group)
        jack = S.jackhmmer_targets(group)
        query = pfam | hmm | jack
        in_db = S.universe_intersect(group, query)
        pfam_db, hmm_db, jack_db = pfam & in_db, hmm & in_db, jack & in_db
        for label, lo, hi in LENGTH_BANDS:
            def band(s):
                return {a for a in s if lo <= length.get(a, 0) < hi}
            p, h = band(pfam_db), band(hmm_db)
            rows.append([group, label, len(p | h), len(p), len(h), len(p & h),
                         len(h - p), len(p - h), ""])
        # jackhmmer's extra targets are **not** banded. Their length is not
        # recorded anywhere this task can read — they are outside the census
        # by construction — and defaulting an unknown length to zero would
        # file all 19,969 vertebrate accretions in the "<1000 aa" band and
        # print a measurement that was never made.
        rows.append([group, "(all lengths; accreted targets have no "
                     "recorded length)", "", "", "", "", "", "",
                     len(jack_db - hmm_db)])
        summary[group] = {
            "db_accessions": S.universe_size(group),
            "hmm_targets_scored": len(scored),
            "hmm_targets_declined_by_gate": len(scored - hmm),
            "pfam_in_db": len(pfam_db), "hmm_in_db": len(hmm_db),
            "shared": len(pfam_db & hmm_db),
            "hmm_only": len(hmm_db - pfam_db),
            "pfam_only": len(pfam_db - hmm_db),
            "jackhmmer_beyond_hmmsearch": len(jack_db - hmm_db),
            "pfam_total_outside_db": len(pfam - in_db),
        }
        log(f"head_to_head {group}: pfam {len(pfam_db)}, hmm {len(hmm_db)}, "
            f"hmm-only {len(hmm_db - pfam_db)}, pfam-only "
            f"{len(pfam_db - hmm_db)}")
    S.write_tsv(S.out_dir() / "head_to_head.tsv",
                ["database", "length_band", "n_union", "n_pfam", "n_hmm",
                 "n_shared", "n_hmm_only", "n_pfam_only",
                 "n_jackhmmer_beyond_hmmsearch"], rows)
    return rows, summary


# --------------------------------------------------------------- 4. cost/yield

#: Which census step each channel belongs to. Two channels ran inside one
#: step (hmmsearch and jackhmmer over the same DB), so the records a step
#: added are attributed to the **step**, never split between its channels —
#: a split would need a counterfactual neither session ran.
CHANNEL_STEP = {
    "app_search": CENSUS_STEPS[0][1],
    "pfam_enumeration": CENSUS_STEPS[1][1],
    "profile_hmm_vert": CENSUS_STEPS[2][1],
    "jackhmmer_vert": CENSUS_STEPS[2][1],
    "genome_miniprot_vert": CENSUS_STEPS[3][1],
    "profile_hmm_nonvert": CENSUS_STEPS[4][1],
    "jackhmmer_nonvert": CENSUS_STEPS[4][1],
    "genome_miniprot_nonvert": CENSUS_STEPS[5][1],
}


def _sweep_seconds(accessions) -> float:
    total = 0.0
    for acc in accessions:
        total += S.fnum(S.read_json(S.sweep_dir(acc) / "summary.json")
                        .get("total_s"), float, 0.0)
    return total


def _s23_seconds() -> float:
    root = S.data_root() / "s23_sweep"
    if not root.exists():
        return 0.0
    total = 0.0
    for path in root.glob("*/summary.json"):
        total += S.fnum(S.read_json(path).get("total_s"), float, 0.0)
    return total


def cost_yield(sets: dict, growth: list[list], log=S.log) -> list[list]:
    """Elapsed seconds per channel against the records its step added.

    Where a session did not record a wall clock the cell says so rather than
    reading zero: the InterPro enumeration's own statistics were written on a
    re-parse from archived pages, and the app's four search bundles carry no
    timing at all.
    """
    added = {r[1]: r[3] for r in growth}
    secs: dict[str, float] = {}
    note: dict[str, str] = {}

    blob = S.read_json(S.RESULTS / "census_v2" / "interpro_enumeration_stats.json")
    walk = blob.get("signatures", blob)
    secs["pfam_enumeration"] = sum(
        S.fnum((v or {}).get("elapsed_s"), float, 0.0)
        for v in walk.values() if isinstance(v, dict))
    note["pfam_enumeration"] = ("measured" if secs["pfam_enumeration"] > 1
                                else "not recorded (statistics written on a "
                                     "re-parse from archived pages)")

    hs = S.read_json(S.HMM_DIR / "sweep_stats_hmmsearch.json")
    secs["profile_hmm_vert"] = sum(
        S.fnum(v.get("elapsed_s"), float, 0.0)
        for v in ((hs.get("hmmsearch") or {}).get("profiles") or {}).values())
    note["profile_hmm_vert"] = "measured"

    jv = 0.0
    for tag in ("itpr1_human", "itpr_fly", "itpr_acanthamoeba"):
        for run in (S.read_json(
                S.HMM_DIR / f"sweep_stats_jackhmmer_{tag}.json")
                .get("jackhmmer") or {}).values():
            jv += S.fnum(run.get("elapsed_s"), float, 0.0)
    secs["jackhmmer_vert"] = jv
    note["jackhmmer_vert"] = "measured"

    s20 = S.read_json(S.S20_DIR / "sweep_stats.json")
    nv_h = nv_j = 0.0
    for g in S.NONVERT_GROUPS:
        blob = S.read_json(S.S20_DIR / f"proteome_db_stats_{g}.json")
        nv_h += S.fnum(blob.get("hmmsearch_elapsed_s"), float, 0.0)
    for key, val in _walk_elapsed(s20):
        if "jackhmmer" in key:
            nv_j += val
        elif "hmmsearch" in key:
            nv_h += val
    secs["profile_hmm_nonvert"] = nv_h
    secs["jackhmmer_nonvert"] = nv_j
    note["profile_hmm_nonvert"] = "measured" if nv_h else "not recorded"
    note["jackhmmer_nonvert"] = "measured" if nv_j else "not recorded"

    accs = {c["accession"] for c in S.ledger_cells()}
    secs["genome_miniprot_vert"] = _sweep_seconds(accs)
    note["genome_miniprot_vert"] = "measured (sum of per-genome sweep totals)"
    secs["genome_miniprot_nonvert"] = _s23_seconds()
    note["genome_miniprot_nonvert"] = ("measured (sum of per-genome sweep "
                                       "totals)" if secs[
                                           "genome_miniprot_nonvert"]
                                       else "not recorded")
    secs["app_search"] = 0.0
    note["app_search"] = "not recorded"

    step_secs: dict[str, float] = {}
    for m in S.METHOD_ORDER:
        step_secs[CHANNEL_STEP[m]] = (step_secs.get(CHANNEL_STEP[m], 0.0)
                                      + secs.get(m, 0.0))
    rows = []
    for m in S.METHOD_ORDER:
        step = CHANNEL_STEP[m]
        n_new = added.get(step, 0)
        sec = secs.get(m, 0.0)
        measured = note.get(m, "").startswith("measured")
        rows.append([m, S.METHOD_LABEL[m], step, len(sets.get(m, ())), n_new,
                     round(sec, 1) if measured else "",
                     round(sec / 3600, 2) if measured else "",
                     round(step_secs[step] / n_new, 1)
                     if n_new and measured else "", note.get(m, "")])
    S.write_tsv(S.out_dir() / "method_cost.tsv",
                ["method", "label", "census_step", "n_records_returned",
                 "n_records_step_added_to_census", "elapsed_s", "elapsed_h",
                 "step_seconds_per_new_record", "timing_source"], rows)
    total = sum(v for k, v in secs.items()
                if note.get(k, "").startswith("measured"))
    log(f"method_cost: {total / 3600:.1f} recorded compute-hours over "
        f"{len(rows)} channels")
    return rows


def _walk_elapsed(blob, prefix=""):
    """Every `elapsed_s` in a nested statistics blob, with its key path."""
    if isinstance(blob, dict):
        for k, v in blob.items():
            if k == "elapsed_s":
                yield (prefix, S.fnum(v, float, 0.0))
            else:
                yield from _walk_elapsed(v, f"{prefix}/{k}")
    elif isinstance(blob, list):
        for v in blob:
            yield from _walk_elapsed(v, prefix)


def run(log=S.log) -> dict:
    sets, _ = record_sets(log)
    growth = census_growth(log)
    pairwise(sets, log)
    _, h2h = head_to_head(sets, log)
    _, genes = R.gene_recovery(log)
    R.db_status(log)
    cost_yield(sets, growth, log)
    summary = {"records": {m: len(s) for m, s in sets.items()},
               "head_to_head": h2h, "genes": genes,
               "census_v6_records": len(S.census_rows("v6"))}
    S.write_json(S.out_dir() / "contribution_summary.json", summary)
    return summary


if __name__ == "__main__":
    run()
