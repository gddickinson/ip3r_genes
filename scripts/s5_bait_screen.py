"""s5_bait_screen.py — D5's chimera screen, run with this project's own profiles.

R4 of `s5_bait_spec.py`. A bait that is a fusion, a read-through or a
mis-joined gene model passes every rule a table can check — its census call is
right, its architecture count is 5/5, its length is in band — and then drags a
whole clade's loci onto the wrong coordinates. What catches it is the shape of
its alignment to a profile: a real family member aligns to `itpr.hmm` in one
long envelope, while a chimera aligns over part of itself and leaves a large
interior or terminus unexplained.

So the screen asks two things of every candidate:

  1. **Family, positively (D14).** Score against `itpr.hmm` *and* `ryr.hmm`
     and assign under S3's own relative margin. A bait must be assigned to
     the family its label claims. This is the same instrument, with the same
     thresholds, that called 16,039 census records.

  2. **Envelope.** The winning profile must cover at least
     `MIN_ENVELOPE_FRAC` of the bait, with no unaligned run inside that
     envelope longer than `MAX_INTERNAL_GAP_AA`.

`self_test()` is a negative control for the screen itself, run on every
build: a chimera is synthesised from two panel baits and the screen must
reject it. A screen that has never rejected anything is not evidence.
"""

from __future__ import annotations

import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

import s5_bait_spec as spec                                # noqa: E402
from s3_assign import REL_MARGIN, MIN_SCORE, assign        # noqa: E402
from s3_hmm_lib import best_hits_by_target, parse_domtblout  # noqa: E402

HMM_DIR = PROJECT_ROOT / "results" / "hmm_sweep"

NO_HIT = {"assignment": "no_hit", "rel_margin": 0.0, "itpr_score": 0.0,
          "ryr_score": 0.0, "profile_reason": "no profile hit",
          "envelope_frac": 0.0, "max_internal_gap": 0, "n_segments": 0}


def write_fasta(path: Path, seqs: dict[str, str]) -> None:
    with path.open("w") as out:
        for k, v in seqs.items():
            out.write(f">{k}\n")
            for i in range(0, len(v), 60):
                out.write(v[i:i + 60] + "\n")


def hmmsearch(profile: Path, faa: Path, domtbl: Path) -> None:
    proc = subprocess.run(
        ["hmmsearch", "--domtblout", str(domtbl), "--cpu", "4",
         "-o", os.devnull, str(profile), str(faa)],
        capture_output=True, text=True)
    if proc.returncode != 0:
        raise SystemExit(f"hmmsearch failed on {profile.name}: "
                         f"{proc.stderr.strip()[:300]}")


def profile_hits(domtbl: Path) -> dict[str, dict]:
    """bait id -> best-hit summary, keyed on the FASTA name itself.

    Deliberately *not* `s3_assign.load_profile_hits`, which re-parses each
    header as a UniProt one and keys on the accession it finds there. These
    headers are bait ids (`ACC|PARALOG|Species`), so that parse takes the
    wrong field; keying on the target name is exact. `assign()` reads a few
    metadata fields off the hit, so they are supplied empty — the screen
    scores sequences, and every label it could read is already known.
    """
    return {target: {**hit, "accession": target, "gene": "", "species": "",
                     "taxon_id": "", "reviewed": False, "protein_name": ""}
            for target, hit in
            best_hits_by_target(parse_domtblout(domtbl)).items()}


def _merge(spans: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not spans:
        return []
    spans = sorted(spans)
    out = [list(spans[0])]
    for s, e in spans[1:]:
        if s <= out[-1][1] + 1:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(a, b) for a, b in out]


def envelope_of(domtbl: Path) -> dict[str, dict]:
    """target -> merged alignment envelope on the *bait*, and its worst gap."""
    spans: dict[str, list[tuple[int, int]]] = defaultdict(list)
    tlen: dict[str, int] = {}
    for r in parse_domtblout(domtbl):
        spans[r["target_name"]].append((r["ali_from"], r["ali_to"]))
        tlen[r["target_name"]] = r["tlen"]
    out = {}
    for target, sp in spans.items():
        merged = _merge(sp)
        covered = sum(b - a + 1 for a, b in merged)
        gaps = [b1 - a2 - 1 for (_, a2), (b1, _) in zip(merged, merged[1:])]
        out[target] = {"envelope_frac": round(covered / max(1, tlen[target]), 4),
                       "max_internal_gap": max(gaps) if gaps else 0,
                       "n_segments": len(merged), "length": tlen[target]}
    return out


def screen(seqs: dict[str, str], work: Path) -> dict[str, dict]:
    """Score every candidate against both S3 profiles; return the verdict."""
    faa = work / "candidates.faa"
    write_fasta(faa, seqs)
    hits, envs = {}, {}
    for fam, profile in (("itpr", HMM_DIR / "itpr.hmm"),
                         ("ryr", HMM_DIR / "ryr.hmm")):
        domtbl = work / f"{fam}.domtblout"
        hmmsearch(profile, faa, domtbl)
        hits[fam] = profile_hits(domtbl)
        envs[fam] = envelope_of(domtbl)
    by_id = {r["accession"]: r
             for r in assign(hits["itpr"], hits["ryr"],
                             rel_margin=REL_MARGIN, min_score=MIN_SCORE)}
    verdicts: dict[str, dict] = {}
    for name in seqs:
        a = by_id.get(name)
        v = {"assignment": (a or {}).get("assignment", "no_hit"),
             "rel_margin": (a or {}).get("rel_margin", 0.0),
             "itpr_score": (a or {}).get("itpr_score", 0.0),
             "ryr_score": (a or {}).get("ryr_score", 0.0),
             "profile_reason": (a or {}).get("reason", "no profile hit")}
        win = "itpr" if v["itpr_score"] >= v["ryr_score"] else "ryr"
        env = envs[win].get(name) or {}
        v.update({"envelope_frac": env.get("envelope_frac", 0.0),
                  "max_internal_gap": env.get("max_internal_gap", 0),
                  "n_segments": env.get("n_segments", 0)})
        verdicts[name] = v
    return verdicts


def verdict(v: dict, want_family: str, spec_mod=None) -> tuple[bool, str]:
    """D5's verdict on one candidate.

    `spec_mod` lets a second panel be screened against its own thresholds and
    its own length band without copying this module. S23 passes
    `s23_bait_spec`; everything else gets S5's. The screen itself is
    deliberately *not* parameterised — one screen, one instrument.
    """
    spec = spec_mod or globals()["spec"]
    if v["assignment"] != want_family:
        return False, (f"profile assignment {v['assignment']} != {want_family} "
                       f"({v['profile_reason']})")
    if v["envelope_frac"] < spec.MIN_ENVELOPE_FRAC:
        return False, (f"{want_family.lower()}.hmm covers only "
                       f"{v['envelope_frac']:.0%} of the bait "
                       f"(< {spec.MIN_ENVELOPE_FRAC:.0%}) — partial or chimeric")
    if v["max_internal_gap"] > spec.MAX_INTERNAL_GAP_AA:
        return False, (f"{v['max_internal_gap']} aa unaligned inside the "
                       f"envelope (> {spec.MAX_INTERNAL_GAP_AA}) — the shape "
                       "of a fusion or a mis-joined model")
    return True, (f"{want_family.lower()}.hmm wins at {v['rel_margin']:.0%} "
                  f"margin over an envelope covering {v['envelope_frac']:.0%} "
                  f"in {v['n_segments']} segment(s)")


def _build_cases(panel: dict[str, str], meta: dict[str, dict]
                 ) -> tuple[dict[str, str], dict[str, tuple[str, str, str]]]:
    """Synthesise the negative controls from real panel baits.

    One per rule, aimed at what that rule is actually responsible for. The
    first version of this test conflated them: it asked the *envelope* test
    to reject a truncation, which it cannot and should not — a protein
    truncated to 40 % of its length aligns 100 % of *itself* to the profile
    in one clean segment. Truncation is the length band's job (R2), fusion
    is the envelope's, and cross-family joining is the margin's.
    """
    def longest(family: str) -> tuple[str, str]:
        ids = [i for i in panel if meta.get(i, {}).get("family") == family]
        best = max(ids, key=lambda i: len(panel[i]))
        return best, panel[best]

    itpr_id, itpr = longest("ITPR")
    ryr_id, ryr = longest("RYR")
    half_i, half_r = len(itpr) // 2, len(ryr) // 2

    cases = {
        "fusion_itpr_ryr": itpr[:int(len(itpr) * 0.45)] + ryr[-int(len(ryr) * 0.55):],
        "insertion_in_itpr": itpr[:half_i] + ryr[half_r:half_r + 600] + itpr[half_i:],
        "truncated_itpr": itpr[:int(len(itpr) * 0.40)],
    }
    described = {
        "fusion_itpr_ryr": (
            "screen", f"{itpr_id}[:45%] + {ryr_id}[-55%:]",
            "an ITPR N-terminus joined to a RyR C-terminus — the two profiles "
            "score it too closely to call, so D14's margin rejects it"),
        "insertion_in_itpr": (
            "screen", f"{itpr_id} with 600 aa of {ryr_id} inserted mid-sequence",
            "a mis-joined gene model — the profile envelope splits around the "
            "insertion, which is what the internal-gap test is for"),
        "truncated_itpr": (
            "length_band", f"{itpr_id}[:40%]",
            "a truncated model — aligns cleanly to the profile over 100 % of "
            "itself, so only R2's length band can catch it (S2 found seven "
            "real records of exactly this shape)"),
    }
    return cases, described


def self_test(panel: dict[str, str], meta: dict[str, dict],
              work: Path, spec_mod=None) -> list[dict]:
    """Negative controls for the selection rules. Every case must be rejected.

    Run on every build. Without it the screen is machinery no result depends
    on: across this panel it rejects nothing, which is either a clean
    shortlist or a screen that cannot fire, and those look identical from the
    outside.
    """
    spec_here = spec_mod or spec
    cases, described = _build_cases(panel, meta)
    verdicts = screen(cases, work)
    rows = []
    for name, seq in cases.items():
        by_rule, built, note = described[name]
        if by_rule == "length_band":
            ok, why = spec_here.passes_shape("ITPR", len(seq), "5")
            passed = ok
            v = dict(NO_HIT)
            why = why or "inside the bait length band"
        else:
            v = verdicts.get(name, dict(NO_HIT))
            passed, why = verdict(v, "ITPR", spec_here)
        rows.append({"case": name, "rejected_by": by_rule, "built_from": built,
                     "description": note, "rejected": int(not passed),
                     "reason": why, "length": len(seq), **v})
        if passed:
            raise SystemExit(
                f"self-test FAILED: the synthetic {name} ({note}; built from "
                f"{built}) was accepted by rule '{by_rule}' — {why}.\n"
                "  That rule is not discriminating; do not trust a panel it "
                "approved.")
    return rows
