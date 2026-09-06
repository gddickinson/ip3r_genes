"""s23_bait_select.py — measuring, ranking and filling the S23 bait slots.

Split out of `s23_build_baits.py` to keep both under the project's 500-line
limit (the `s3_report.py` / `s3_report_d10.py` pattern). The driver owns the
outputs and the rule table; this owns the decisions:

  `measure()`          B3/B5 — the length band and the architecture-exception
                       floor, **stratified per band**. Measured over the whole
                       population they are the arthropod family's numbers.
  `build_candidates()` every in-scope record judged by the shape rules, with
                       the reason it failed kept.
  `load_s3_seeds()`    B6, enforced: a seed whose census row has drifted
                       aborts the build.
  `fill()`             the slots, with B3's spread rule shaping the shortlist
                       rather than filtering its output — the distinction that
                       cost this panel its *Chlamydomonas* bait once.
"""

from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import s23_bait_spec as spec                                   # noqa: E402
import s5_bait_screen as screen_lib                            # noqa: E402
from s3_hmm_lib import acc_key                                 # noqa: E402
from s5_build_baits import (extract_from_archive, read_fasta,   # noqa: E402
                            read_tsv)

HMM_DIR = PROJECT_ROOT / "results" / "hmm_sweep"

#: How many ranked candidates per slot go to the (expensive) chimera screen.
DEFAULT_SHORTLIST = 6


def _f(x) -> float:
    try:
        return float(x or 0.0)
    except (TypeError, ValueError):
        return 0.0


def _i(x) -> int:
    try:
        return int(x or 0)
    except (TypeError, ValueError):
        return 0


# ------------------------------------------------------------- measurement

def measure(rows: list[dict]) -> dict:
    """B3/B5 — the length band and the architecture-exception floor.

    Both from the same population: in-scope, complete-architecture,
    non-vertebrate ITPR records, **capped per band** (B3's
    `MEASURE_CAP_PER_BAND`). Measuring them together matters — the floor is a
    percentile of the scores *of the records that define the band*, so a band
    measured on one population and a floor on another would not compose.
    """
    complete = [r for r in rows if str(r.get("n_itpr_arch")) == "5"
                and _i(r.get("length"))]
    by_band_len: dict[str, list[int]] = defaultdict(list)
    for r in complete:
        by_band_len[spec.band_of(r)].append(_i(r["length"]))
    band = spec.measure_band(spec.stratify(by_band_len))
    lo, hi = band
    by_band_score: dict[str, list[float]] = defaultdict(list)
    for r in complete:
        if lo <= _i(r["length"]) <= hi and _f(r.get("itpr_score")) > 0:
            by_band_score[spec.band_of(r)].append(_f(r["itpr_score"]))
    floor = spec.measure_arch_exception_floor(spec.stratify(by_band_score))
    spec.set_measured_band(band, floor)
    unstratified = spec.measure_band([_i(r["length"]) for r in complete])
    return {"band_aa": list(band), "arch_floor_bits": round(floor, 1),
            "measured_on_records": len(complete),
            "band_aa_unstratified": list(unstratified),
            "measure_cap_per_band": spec.MEASURE_CAP_PER_BAND,
            "bands_in_measurement": len(by_band_len)}


# -------------------------------------------------------------- candidates

def bait_id(row: dict, band: str, role: str) -> str:
    sp = (row.get("species") or row.get("organism") or "unknown")
    sp = sp.split("(")[0].strip().replace(" ", "_")[:38]
    gene = (row.get("gene") or "unnamed").replace(" ", "_")[:20]
    return f"{row['accession']}|{gene}|{sp}|{band}|{role}"


def build_candidates(census: list[dict]) -> tuple[dict, dict, list, dict]:
    """(slot -> ranked candidates, band -> record count, rejected rows).

    A slot is `(band, role)`. Every in-scope record is judged by the rules and
    the reason it failed is kept, so `bait_candidates.tsv` shows what the panel
    passed over rather than only what it took.
    """
    by_band_len: dict[str, list[int]] = defaultdict(list)
    slots: dict[tuple, list[dict]] = defaultdict(list)
    rejected: list[dict] = []
    in_band: dict[tuple, int] = defaultdict(int)
    for r in census:
        band = spec.band_of(r)
        if not band:
            continue
        call = r.get("call")
        if call not in (spec.FAMILY_ITPR, spec.FAMILY_RYR):
            continue
        length = _i(r.get("length"))
        if call == spec.FAMILY_ITPR:
            by_band_len[band].append(length)
        if call == spec.FAMILY_RYR and not spec.expects_ryr(band):
            continue                      # B2: RyR is a control only in metazoa
        in_band[(band, call)] += 1
        if call == spec.FAMILY_ITPR:
            ok, why, is_exc = spec.passes_shape_with_exception(
                call, length, r.get("n_itpr_arch"), r.get("itpr_score"))
        else:
            ok, why = spec.passes_shape(call, length, r.get("n_itpr_arch"))
            is_exc = False
        row = dict(r, band=band, role=call, arch_exception=int(is_exc),
                   shape_reason=why or "in band")
        if not ok:
            rejected.append(dict(row, verdict="shape", detail=why))
            continue
        # B1: the record's own source must be the census's, not a bare sweep
        # hit with no annotation behind it (B5 already excludes those for
        # ITPR; this catches the RyR side).
        slots[(band, call)].append(row)
    ref_len = {}
    for band, lens in by_band_len.items():
        ins = sorted(lens)
        ref_len[band] = ins[len(ins) // 2] if ins else 0
    for slot, rows in slots.items():
        band = slot[0]
        rows.sort(key=lambda r: spec.rank_candidate(r, ref_len.get(band, 2800)))
    return slots, ref_len, rejected, dict(in_band)


def load_s3_seeds(census_by_acc: dict[str, dict]) -> tuple[dict, list]:
    """B6 — the S3 seeds that fall in an S23 band, with their census rows.

    Enforced, not trusted: a seed whose census row has since changed call,
    band or length aborts the build. S3's own builder found a mis-transcribed
    architecture count this way on its first run.
    """
    seqs = read_fasta(HMM_DIR / "itpr_seed.faa")
    seqs.update(read_fasta(HMM_DIR / "ryr_seed.faa"))
    manifest = read_tsv(HMM_DIR / "seed_manifest.tsv")
    kept, problems = [], []
    for row in manifest:
        acc = acc_key(row["accession"])
        crow = census_by_acc.get(acc)
        if crow is None:
            continue
        band = spec.band_of(crow)
        if not band:
            continue                       # vertebrate seeds: B7 leaves them out
        role = spec.FAMILY_RYR if row["profile"] == "ryr" else spec.FAMILY_ITPR
        if role == spec.FAMILY_RYR and not spec.expects_ryr(band):
            continue
        if row["id"] not in seqs:
            problems.append(f"{row['id']}: no sequence in the S3 seed FASTAs")
            continue
        seq = seqs[row["id"]]
        if _i(row["length"]) and abs(len(seq) - _i(row["length"])) > 0:
            problems.append(
                f"{row['id']}: manifest length {row['length']} != "
                f"{len(seq)} in the FASTA")
            continue
        kept.append(dict(crow, band=band, role=role, seed_id=row["id"],
                         sequence=seq, s3_clade=row["clade"],
                         s3_note=row.get("note", ""), arch_exception=0,
                         shape_reason="S3 seed, reused under B6"))
    if problems:
        raise SystemExit("S3 seed reuse (B6) failed its own checks:\n  "
                         + "\n  ".join(problems))
    return {r["seed_id"]: r for r in kept}, kept


# ------------------------------------------------------------------ filling

def fill(slots: dict, seeds: list, shortlist: int, work: Path,
         in_band_counts: dict | None = None) -> tuple[list, list, list]:
    """Fill every slot; returns (panel rows, unfilled slots, screen rows).

    The screen is run once over the union of every shortlist rather than per
    slot: one hmmsearch per profile over a few hundred sequences, and one set
    of scores that every slot's decision is read off.
    """
    seeds_by_slot: dict[tuple, list[dict]] = defaultdict(list)
    for s in seeds:
        seeds_by_slot[(s["band"], s["role"])].append(s)

    wanted: dict[tuple, list[dict]] = {}
    for band in spec.BAND_ORDER:
        for role in (spec.FAMILY_ITPR, spec.FAMILY_RYR):
            if role == spec.FAMILY_RYR and not spec.expects_ryr(band):
                continue
            # The shortlist must itself be **spread**, not merely the top N by
            # score. A band here is a whole phylum: Chlorophyta's nine
            # shape-passing records are led by six *Cymbomonas* ones, so a
            # score-ranked shortlist never reached *Chlamydomonas reinhardtii*
            # — and applying the spread rule afterwards only removed a bait
            # instead of replacing it (85 baits -> 78, Chlorophyta down to
            # one). Take the best candidate per spread key first, in score
            # order, then fill any remaining slots from the rest.
            ranked = slots.get((band, role), [])
            first_of_key, rest, seen_keys = [], [], set()
            for r in ranked:
                key = spec.spread_key(r)
                if key and key not in seen_keys:
                    seen_keys.add(key)
                    first_of_key.append(r)
                else:
                    rest.append(r)
            wanted[(band, role)] = (first_of_key + rest)[:shortlist]

    # sequences for every shortlisted candidate, in one pass over the archive
    need = {r["accession"] for rows in wanted.values() for r in rows}
    need -= {s["accession"] for s in seeds}
    archive = extract_from_archive({acc_key(a) for a in need})

    cand_seqs: dict[str, str] = {}
    cand_meta: dict[str, dict] = {}
    for slot, rows in wanted.items():
        for r in rows:
            bid = bait_id(r, r["band"], r["role"])
            seq = archive.get(acc_key(r["accession"]))
            if not seq:
                continue
            cand_seqs[bid] = seq
            cand_meta[bid] = dict(r, id=bid, origin="census v5",
                                  sequence=seq)
    for s in seeds:
        bid = bait_id(s, s["band"], s["role"])
        cand_seqs[bid] = s["sequence"]
        cand_meta[bid] = dict(s, id=bid, origin=f"S3 seed ({s['seed_id']})")

    verdicts = screen_lib.screen(cand_seqs, work)
    screen_rows = []
    for bid, v in verdicts.items():
        m = cand_meta[bid]
        ok, why = screen_lib.verdict(v, m["role"], spec)
        screen_rows.append({"id": bid, "accession": m["accession"],
                            "band": m["band"], "role": m["role"],
                            "length": len(cand_seqs[bid]), "passed": int(ok),
                            "reason": why, **v})
        cand_meta[bid]["screen_ok"] = ok
        cand_meta[bid]["screen_reason"] = why

    panel, unfilled = [], []
    for band in spec.BAND_ORDER:
        for role in (spec.FAMILY_ITPR, spec.FAMILY_RYR):
            if role == spec.FAMILY_RYR and not spec.expects_ryr(band):
                continue
            seed_pool = [cand_meta[bait_id(s, band, role)]
                         for s in seeds_by_slot.get((band, role), [])
                         if bait_id(s, band, role) in cand_meta]
            census_pool = [cand_meta[bait_id(r, band, role)]
                           for r in wanted.get((band, role), [])
                           if bait_id(r, band, role) in cand_meta]
            pool = seed_pool + census_pool
            # B6: a reused S3 seed is *added to* the band's quota, not
            # subtracted from it. Letting a seed consume the slot cost the
            # first build its Chlamydomonas bait, and the pilot then found
            # nothing in the Chlamydomonas genome.
            need_n = spec.quota(band, role) + len(
                [m for m in seed_pool if m.get("screen_ok")])
            seen, taken, spreads = set(), [], set()
            for m in pool:
                if m["id"] in seen or not m.get("screen_ok"):
                    continue
                # B3's spread rule: a band's second and later baits must come
                # from a different order (or genus, where no order is ranked)
                # than the ones already taken. Without it Chlorophyta filled
                # both its slots with the same *Cymbomonas* genus.
                key = spec.spread_key(m)
                if taken and key and key in spreads:
                    continue
                seen.add(m["id"])
                spreads.add(key)
                taken.append(m)
                if len(taken) >= need_n:
                    break
            if len(taken) < spec.quota(band, role):
                # relax the spread rather than leave a band unbaited: a bait
                # from a sister genus beats no bait at all.
                for m in pool:
                    if m["id"] in seen or not m.get("screen_ok"):
                        continue
                    seen.add(m["id"])
                    taken.append(m)
                    if len(taken) >= need_n:
                        break
            panel.extend(taken)
            # A slot is unfilled when it has fewer baits than the *census*
            # quota asks for, not fewer than quota+seeds: a band carrying only
            # its reused S3 seed has a bait, and reporting it unfilled would
            # make the seed reuse look like a gap.
            if len(taken) < spec.quota(band, role):
                # Say which stage lost the slot. The first version of this
                # message reported "no census record in this band at all" for
                # seven bands, four of which have dozens of records that fail
                # the *shape* rules — a reason that points at the census when
                # the truth is a threshold. An unfilled slot is a limitation
                # of the panel and has to name its own cause.
                n_records = (in_band_counts or {}).get((band, role), 0)
                need_n = spec.quota(band, role)   # report the census quota
                n_shape = len(slots.get((band, role), []))
                n_seq = len(pool)
                if not n_records:
                    why = "no census record in this band at all"
                elif not n_shape:
                    why = (f"{n_records} census record(s), none passing the "
                           "shape rules (length band / architecture); the "
                           "band's records are fragmentary or divergent")
                elif not n_seq:
                    why = (f"{n_shape} record(s) pass the shape rules but "
                           "none has a sequence in the archived FASTAs")
                else:
                    why = (f"{n_seq} shortlisted candidate(s), none passing "
                           "the chimera screen")
                unfilled.append({
                    "band": band, "role": role, "wanted": need_n,
                    "filled": len(taken), "records": n_records,
                    "shape_ok": n_shape, "with_sequence": n_seq, "why": why})
    return panel, unfilled, screen_rows


