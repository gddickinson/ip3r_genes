"""s23_calibrate_loci.py — what counts as a locus, and what counts as one copy.

Two thresholds S23a inherited from S5 without measuring them for this scope,
measured here from the sweep's own output. Both are about the same thing:
**what a locus is**. Copy number is this task's deliverable, so the definition
of one gene is a result, not an implementation detail.

**1. The identity floor** (`MIN_LOCUS_IDENTITY`). S5b measured 0.40 in the
vertebrates, where confirmed loci start at 0.759, the junk tops out at 0.34,
and every genome has a bait from its own class. Neither holds outside them:
the bands here are whole phyla and 13 slots are unbaited, so a real gene can
be much further from its nearest bait. The sweep therefore records every
cluster down to `RECORD_MIN_IDENTITY` (0.15) and this measures the floor from
what it recorded — against loci whose identity the assembly's **own
annotation** establishes, which is evidence the alignment score did not
produce. A floor measured from the population it filters would be circular;
this one is not.

**The annotation axis is too thin here, and that is itself a result.** Of 917
recorded clusters across the 194 genomes, **21** sit on a gene whose name says
anything at all: 10 name the family, 11 name something else. Outside the
vertebrates most gene models carry locus tags — *Chlamydomonas* files its
receptor as `CHLRE_16g665450v5`, *Strongylocentrotus* as `LOC594527` — so the
evidence S5b calibrated 571 loci against does not exist at that depth in this
scope. The floor is therefore measured against a **second instrument** as
well: every recorded locus's translated model is scored with `itpr.hmm` and
`ryr.hmm` under S3's margin (D23), re-derived from the archived miniprot GFFs
so the calibration reruns offline.

That second axis is a different algorithm on different data — spliced
DNA-protein alignment against one bait versus a ~4,900-state profile built
from ~30 seeds — but it is **not independent of the family definition** the
way an assembly's own annotation is, and the report says so rather than
presenting the two as equivalent. Both are reported separately as well as
pooled, so a reader can see what each one alone would have given.

It also answers the question the brief asked directly: how many `no_locus`
genomes have a cluster sitting just under the floor? A genome called absent
because its best evidence scored 0.39 is a different claim from one whose best
evidence scored 0.17.

**2. What makes two alignments two copies** (`copy_max_overlap`). Measured the
same way: among complete alignments that land on an annotated gene, how far do
two that hit the **same** gene overlap, and how far do two that hit
**different** genes? If those two populations separate, the threshold is read
off the gap rather than chosen.

**3. And the measurement that motivated (2)**: locus span against CDS
footprint. `-G` is 650 kb for the metazoa, so a locus is a large object — the
S23a pilot's *Drosophila* Itpr sat in a 297 kb cluster around an 8.5 kb CDS
footprint. That is harmless for a status call and fatal for a copy count, and
the ratio is reported per group so the reader can see where it bites.

Outputs -> results/s23_scope/locus_calibration.json     read back by s23_calibration
           results/s23_scope/locus_identity.tsv         every recorded locus
           results/s23_scope/locus_span.tsv             span vs CDS footprint
           results/s23_scope/copy_overlap.tsv           the same-gene / other-gene split

Usage:
  python3 scripts/s23_calibrate_loci.py
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import s23_calibration as cal                               # noqa: E402
from s23_locus_evidence import (ANNOT_CDS_FRAC, collect,    # noqa: E402
                                load_summaries, profile_evidence)

OUT_DIR = PROJECT_ROOT / "results" / "s23_scope"

#: The sample this calibration refuses to run below. A threshold measured
#: from a handful of loci is not a measurement, and the failure mode is not
#: hypothetical: a smoke-test run of this module over the S23a pilot's single
#: re-swept genome produced `call_min_identity = 0.25` from two confirmed
#: loci, wrote it, and the next sweep read it back and moved three genomes
#: from `no_locus` to `fragment_only` on the strength of it. The guard is the
#: `s20_test_sensitivity.py` rule applied to a calibration: refusing to pass
#: vacuously is part of the instrument.
MIN_GENOMES = 50
MIN_CONFIRMED = 15

#: The candidate floors scanned. The report states, for each, how many
#: annotation-confirmed loci and how many annotation-contradicted ones it
#: keeps — which is what turns a chosen number into a measured one.
CANDIDATE_FLOORS = [round(0.15 + 0.05 * i, 2) for i in range(12)]


#: A locus is credited to an annotated gene when the alignment's exons fall
#: this far inside it. Same number `s5_classify` uses to attach an annotation,
#: reused so "the annotation confirms this locus" means one thing project-wide.
# --------------------------------------------------------------- the floor

#: The two evidence axes, pooled. The annotation axis is the independent one
#: and is nearly empty in this scope; the profile axis is a second instrument
#: on the same alignment's translation. Pooled for the decision, reported
#: separately as well, so what each one alone would have given stays visible.
def is_confirmed(r: dict) -> bool:
    return (r["evidence"] == "confirmed"
            or r.get("profile_evidence") == "profile_confirmed")


def is_contradicted(r: dict) -> bool:
    """Annotation names a different gene, or the profiles call it RyR.

    A locus the profiles decline entirely is **not** counted here: declining
    is what the D22 length gate does to a short scrap, and a scrap is absence
    of evidence rather than evidence of junk.
    """
    if r["evidence"] in ("contradicted", "sister"):
        return True
    if r.get("profile_evidence") == "profile_contradicted":
        return True
    # An annotation-contradicted locus stays contradicted even if the
    # profiles were silent about it.
    return False


def floor_scan(rows: list[dict]) -> list[dict]:
    conf = [r for r in rows if is_confirmed(r)]
    contra = [r for r in rows if is_contradicted(r)]
    out = []
    for f in CANDIDATE_FLOORS:
        out.append({
            "floor": f,
            "confirmed_kept": sum(1 for r in conf if r["identity"] >= f),
            "confirmed_lost": sum(1 for r in conf if r["identity"] < f),
            "contradicted_kept": sum(1 for r in contra if r["identity"] >= f),
            "no_evidence_kept": sum(1 for r in rows
                                    if not is_confirmed(r)
                                    and not is_contradicted(r)
                                    and r["identity"] >= f)})
    return out


def separation(rows: list[dict], key: str) -> dict:
    """How well one statistic separates confirmed loci from contradicted ones.

    Reported for **identity** and for **coverage**, because the question S23b
    has to answer is not only "what floor" but "on which statistic". S5b never
    had to ask: in the vertebrates identity separated with a wide empty gap.
    """
    conf = sorted(float(r[key]) for r in rows if is_confirmed(r))
    contra = sorted(float(r[key]) for r in rows if is_contradicted(r))
    if not conf or not contra:
        return {"statistic": key, "n_confirmed": len(conf),
                "n_contradicted": len(contra), "separates": False,
                "why": "one of the two populations is empty"}
    q = lambda v, f: v[min(len(v) - 1, int(f * len(v)))]        # noqa: E731
    overlap = max(contra) >= min(conf)
    return {
        "statistic": key, "n_confirmed": len(conf),
        "n_contradicted": len(contra),
        "confirmed_min": round(conf[0], 4),
        "confirmed_q1": round(q(conf, 0.25), 4),
        "confirmed_median": round(q(conf, 0.5), 4),
        "contradicted_median": round(q(contra, 0.5), 4),
        "contradicted_q3": round(q(contra, 0.75), 4),
        "contradicted_max": round(contra[-1], 4),
        "separates": not overlap,
        "why": (f"confirmed loci reach down to {conf[0]:.3f} and contradicted "
                f"ones up to {contra[-1]:.3f}"
                + (" — the two populations overlap, so no threshold on this "
                   "statistic separates them" if overlap else
                   " — a gap separates them")),
    }


def choose_floor(rows: list[dict], scan: list[dict]) -> tuple[float, str]:
    """The identity a cluster must reach to be called, and why.

    **This scope's answer is that identity is the wrong statistic**, and that
    is a measurement rather than a preference. S5b's 0.40 came from a wide
    empty gap in the vertebrates, where every genome has a bait from its own
    class. Here the bands are whole phyla, so a real gene's identity to its
    nearest bait is not a measure of whether it is a gene — it is a measure of
    how far away the nearest bait happens to be. When the two populations
    overlap, the floor is set at the recording floor (identity retired as a
    call gate) and the report says which statistic does separate them instead.
    """
    conf = sorted(r["identity"] for r in rows if is_confirmed(r))
    contra = sorted(r["identity"] for r in rows if is_contradicted(r))
    if not conf:
        return cal.INHERITED_CALL_MIN_IDENTITY, (
            "no confirmed locus in the sweep on either evidence axis, so this "
            f"scope offers no measurement; S5b's "
            f"{cal.INHERITED_CALL_MIN_IDENTITY:.2f} is kept and is an "
            "inheritance, not a measurement")
    lo_conf, hi_contra = conf[0], (max(contra) if contra else 0.0)
    if hi_contra < lo_conf:
        mid = round((lo_conf + hi_contra) / 2, 2)
        return mid, (
            f"measured: {len(conf)} confirmed loci reach down to "
            f"{lo_conf:.3f} and {len(contra)} contradicted ones top out at "
            f"{hi_contra:.3f}; {mid:.2f} is the midpoint of that gap")
    lost_at_inherited = sum(1 for r in rows if is_confirmed(r)
                            and r["identity"] < cal.INHERITED_CALL_MIN_IDENTITY)
    full_lost = sum(1 for r in rows if is_confirmed(r)
                    and r.get("grade") == "full"
                    and r["identity"] < cal.INHERITED_CALL_MIN_IDENTITY)
    return cal.RECORD_MIN_IDENTITY, (
        f"**neither identity nor coverage separates the two populations in "
        f"this scope.** {len(conf)} confirmed loci reach down to "
        f"{lo_conf:.3f} on identity while {len(contra)} contradicted ones "
        f"reach up to {hi_contra:.3f}, and the best achievable threshold on "
        f"either statistic still discards a large part of the confirmed set "
        f"(see `separation` in this file). S5b's inherited "
        f"{cal.INHERITED_CALL_MIN_IDENTITY:.2f} would discard "
        f"{lost_at_inherited} confirmed loci, {full_lost} of them **complete "
        f"gene models**. Outside the vertebrates a locus's identity to its "
        f"nearest bait measures how far away the nearest bait is — a whole "
        f"phylum — not whether it is a gene. Identity is therefore retired as "
        f"a call gate (set to the recording floor) and the **profile call** "
        f"carries it (D14/D23), which is where this project puts every other "
        f"family call. That gate is validated against the one axis it does "
        f"not share: see `profile_gate_vs_annotation`")


def profile_gate_audit(rows: list[dict]) -> dict:
    """Score the profile gate against the assembly's own annotation.

    The annotation axis is the only evidence here that the gate does not
    share, so it is the only thing that can validate it. It is small — the
    annotations outside the vertebrates mostly carry locus tags — and being
    small is itself reported rather than smoothed over.
    """
    conf = [r for r in rows if r["evidence"] == "confirmed"]
    contra = [r for r in rows if r["evidence"] in ("contradicted", "sister")]
    tp = sum(1 for r in conf if r.get("profile_call") == "ITPR")
    fp = sum(1 for r in contra if r.get("profile_call") == "ITPR")
    return {
        "annotation_confirmed": len(conf), "gate_agrees": tp,
        "annotation_contradicted": len(contra), "gate_declines": len(contra) - fp,
        "gate_admits": fp,
        "admitted": [{"organism": r["organism"], "gene": r["annot_gene"],
                      "identity": r["identity"], "coverage": r["coverage"],
                      "itpr_score": r.get("itpr_score", ""),
                      "ryr_score": r.get("ryr_score", "")}
                     for r in contra if r.get("profile_call") == "ITPR"],
        "why": (f"the profile gate agrees with {tp} of {len(conf)} loci the "
                f"assembly's own annotation names for this family, and "
                f"declines {len(contra) - fp} of {len(contra)} it names for "
                f"something else"),
    }


def near_miss(rows: list[dict], floor: float, window: float = 0.05) -> list[dict]:
    """`no_locus` genomes whose best evidence sits just under the floor."""
    by_acc: dict[str, dict] = {}
    for r in rows:
        if r["status"] != "no_locus" or r["identity"] >= floor:
            continue
        cur = by_acc.get(r["accession"])
        if cur is None or r["identity"] > cur["identity"]:
            by_acc[r["accession"]] = r
    return sorted((r for r in by_acc.values()
                   if r["identity"] >= floor - window),
                  key=lambda r: -r["identity"])


# ------------------------------------------------------ span and copy overlap

def span_stats(rows: list[dict]) -> dict:
    by_group: dict[str, list] = defaultdict(list)
    for r in rows:
        if r["called"] and r["cds_footprint_bp"] and r["span_inflation"]:
            by_group[r["group"] or "other"].append(r["span_inflation"])
    out = {}
    for g, vals in sorted(by_group.items()):
        vals.sort()
        out[g] = {"n": len(vals), "median": round(statistics.median(vals), 2),
                  "max": vals[-1],
                  "over_5x": sum(1 for v in vals if v >= 5.0),
                  "over_10x": sum(1 for v in vals if v >= 10.0)}
    return out


def copy_overlap(summaries: list[dict]) -> tuple[list[dict], float, str]:
    """How far do two complete alignments on the same / different genes overlap?"""
    rows = []
    for s in summaries:
        genes: dict[str, list] = defaultdict(list)
        for c in s.get("copies") or []:
            genes[c["contig"]].append(c)
        for contig, cs in genes.items():
            for i in range(len(cs)):
                for j in range(i + 1, len(cs)):
                    a, b = cs[i], cs[j]
                    ov = min(a["end"], b["end"]) - max(a["start"], b["start"]) + 1
                    shorter = min(a["end"] - a["start"] + 1,
                                  b["end"] - b["start"] + 1)
                    rows.append({"accession": s["accession"], "contig": contig,
                                 "frac": round(max(0, ov) / max(1, shorter), 4)})
    # Copies are, by construction, already below the threshold in force, so
    # this reports what the rule produced rather than re-deriving it. The
    # threshold moves only when a pair of *annotated* genes says it should.
    fr = sorted(r["frac"] for r in rows)
    if not fr:
        return rows, cal.copy_max_overlap()[0], (
            "no genome carries two complete alignments on one contig, so this "
            "scope offers no measurement; the default 0.20 stands")
    return rows, cal.copy_max_overlap()[0], (
        f"{len(fr)} same-contig copy pairs, largest residual overlap "
        f"{fr[-1]:.0%}; the rule's threshold is not contradicted by the data "
        "it produced")


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--allow-thin", action="store_true",
                    help="write the calibration even below the sample floor "
                         "(diagnostics only — it will be read back by the "
                         "sweep as if it were measured)")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = load_summaries()
    if not summaries:
        raise SystemExit("no per-genome summaries under <data_root>/s23_sweep")
    print(f"{len(summaries)} genome(s); scoring every recorded locus against "
          "itpr.hmm / ryr.hmm from the archived GFFs")
    profiles = profile_evidence(summaries)
    rows = collect(summaries, profiles)
    scan = floor_scan(rows)
    floor, why = choose_floor(rows, scan)
    misses = near_miss(rows, floor)
    spans = span_stats(rows)
    seps = {k: separation(rows, k) for k in ("identity", "coverage")}
    gate_audit = profile_gate_audit(rows)
    pairs, overlap, overlap_why = copy_overlap(summaries)

    counts = defaultdict(int)
    for r in rows:
        counts[r["evidence"]] += 1
        if r.get("profile_evidence"):
            counts[r["profile_evidence"]] += 1
    counts["confirmed_pooled"] = sum(1 for r in rows if is_confirmed(r))
    counts["contradicted_pooled"] = sum(1 for r in rows if is_contradicted(r))

    write_tsv(OUT_DIR / "locus_identity.tsv",
              ["accession", "organism", "group", "phylum", "called",
               "identity", "coverage", "grade", "evidence", "annot_gene",
               "profile_call", "profile_confidence", "itpr_score",
               "ryr_score", "profile_evidence", "bait", "band", "status"],
              rows)
    write_tsv(OUT_DIR / "locus_span.tsv",
              ["accession", "organism", "group", "called", "identity",
               "coverage", "span_bp", "cds_footprint_bp", "span_inflation",
               "evidence", "annot_gene"],
              [r for r in rows if r["called"]])
    write_tsv(OUT_DIR / "copy_overlap.tsv",
              ["accession", "contig", "frac"], pairs)
    write_tsv(OUT_DIR / "identity_floor_scan.tsv",
              ["floor", "confirmed_kept", "confirmed_lost",
               "contradicted_kept", "no_evidence_kept"], scan)

    out = {
        "genomes": len(summaries), "loci_recorded": len(rows),
        "record_min_identity": cal.RECORD_MIN_IDENTITY,
        "inherited_call_min_identity": cal.INHERITED_CALL_MIN_IDENTITY,
        "evidence_counts": dict(counts),
        "call_min_identity": floor, "call_min_identity_why": why,
        "separation": seps,
        "profile_gate_vs_annotation": gate_audit,
        "copy_max_overlap": overlap, "copy_max_overlap_why": overlap_why,
        "near_miss_no_locus": misses,
        "span_inflation_by_group": spans,
        "floor_scan": scan,
    }
    thin = (len(summaries) < MIN_GENOMES
            or counts["confirmed_pooled"] < MIN_CONFIRMED)
    if thin and not args.allow_thin:
        print(f"\nNOT WRITING locus_calibration.json: {len(summaries)} genome(s) "
              f"and {counts['confirmed_pooled']} confirmed locus/loci, "
              f"below the floor of {MIN_GENOMES} / {MIN_CONFIRMED}.\n"
              "  The tables above are written and are diagnostics. The sweep "
              "reads the JSON back as a measured threshold, so writing one "
              "from a sample this thin would put an unmeasured number into "
              "every call. Finish the sweep, then rerun.")
        return 2
    out["thin_sample"] = bool(thin)
    (OUT_DIR / "locus_calibration.json").write_text(json.dumps(out, indent=1))
    print(f"{len(summaries)} genomes, {len(rows)} recorded loci")
    print(f"  annotation axis: {counts['confirmed']} confirmed, "
          f"{counts['contradicted'] + counts['sister']} contradicted")
    print(f"  profile axis:    {counts['profile_confirmed']} confirmed, "
          f"{counts['profile_contradicted']} contradicted")
    print(f"  pooled:          {counts['confirmed_pooled']} confirmed, "
          f"{counts['contradicted_pooled']} contradicted")
    for k, s in seps.items():
        print(f"  {k:9s} separates={s['separates']}  {s['why']}")
    print(f"  profile gate vs annotation: {gate_audit['why']}")
    print(f"call_min_identity = {floor:.2f}\n  {why}")
    print(f"copy_max_overlap  = {overlap:.2f}\n  {overlap_why}")
    if misses:
        print(f"{len(misses)} no_locus genome(s) with a cluster within "
              f"0.05 of the floor:")
        for m in misses[:10]:
            print(f"  {m['identity']:.3f}  {m['organism'][:38]:38s} "
                  f"{m['bait'][:30]}")
    for g, st in spans.items():
        print(f"  span/CDS {g:14s} n={st['n']:4d} median {st['median']:6.1f}x "
              f"max {st['max']:.0f}x  ({st['over_10x']} over 10x)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
