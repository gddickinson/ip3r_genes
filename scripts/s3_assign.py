"""S3 — best-profile assignment with a bit-score margin (D14 at HMM scale).

D14 says ITPR vs RYR is a positive test at every stage; D7 says a best-hit
assignment without a margin is not a call. Here that is: score every target
against **both** profiles and assign it to the better-scoring one, but only
when the margin between the two clears a threshold expressed as a fraction
of the winning score. Inside the band the target is `unassigned` and counted
— exactly as S2 leaves a partial architecture uncalled rather than guessing.

The margin is **relative**, not absolute. A 2,700-aa ITPR and a 600-aa
fragment of one cannot be held to the same bit-score gap: bit scores scale
with alignable length, so a fixed gap would call every long protein
confidently and no short one. `rel_margin = (win - lose) / win` is
length-normalised the way D7's identity margin is.

Nothing here consults a gene symbol, a protein name or a length. Those
columns travel with each row so the assignment can be *audited* against
them — the same separation S2's `rule_audit.tsv` keeps.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_hmm_lib import (  # noqa: E402
    acc_key, best_hits_by_target, parse_domtblout, parse_uniprot_header,
)

# D7's 10 % no-call band, carried over from the identity margin to the
# bit-score margin. A target whose two profile scores are within 10 % of
# the winner is not called by this instrument.
REL_MARGIN = 0.10

# A target scoring against only one profile is still only called when that
# score is meaningful on its own: a handful of bits over a 2,645-state
# profile is not a family assignment.
MIN_SCORE = 30.0

# The gate that matters, and the one the first run of this sweep did not
# have. Both profiles are full-length channel models, so a protein sharing
# a single small module with one of them scores against it — and `ryr.hmm`
# carries SPRY (PF00622), which sits in ~114,000 UniProt proteins and is
# "in no sense RyR-specific" in D14b's own words. Without a gate the sweep
# called 14,981 vertebrate proteins RYR, of which 12,874 matched under a
# tenth of the profile and were named Tnnc2, Rspry1, Cabp1, Rnf123, Ash2l —
# EF-hand and SPRY proteins, not receptors.
#
# So a profile match is family evidence only if it spans at least
# MIN_PROFILE_POSITIONS match states. The number is measured, not tuned:
# in this project's own S0 domain coordinates
# (`results/s0_baseline/review_figures/domain_coords.tsv`) the shortest
# observed instance of PF08709 — the IP3-binding core that *names* the
# family — is 200 aa, while the longest observed SPRY is 137 aa. The floor
# sits in that gap, so it admits any match spanning a family domain and
# excludes any match spanning only the shared module that causes the
# problem. Split gene models survive it: the 200-300 position band is
# almost entirely `Itpr1_0` / `Itpr2_1` / `Ryr3_0` — pieces of real genes.
#
# Cost, measured on the calibration set: 86 of 12,020 records the
# architecture rule called lose their profile verdict (0.7 %).
MIN_PROFILE_POSITIONS = 200

# Above this share of the profile the match spans the whole architecture
# rather than part of it — reported as an evidence class, not used as a gate.
COMPLETE_COVERAGE = 0.50

ASSIGN_FIELDS = [
    "accession", "assignment", "confidence", "evidence", "reason",
    "itpr_score", "ryr_score", "score_margin", "rel_margin",
    "itpr_evalue", "ryr_evalue", "best_evalue",
    "itpr_coverage", "ryr_coverage", "itpr_positions", "ryr_positions",
    "target_coverage",
    "length", "gene", "species", "taxon_id", "reviewed", "protein_name",
]


def _positions(hit: dict | None) -> int:
    """Match states of the profile the alignment actually spans."""
    return round(hit["hmm_coverage"] * hit["qlen"]) if hit else 0


def load_profile_hits(domtbl: Path) -> dict[str, dict]:
    """acc → best-hit summary for one profile's domtblout."""
    out: dict[str, dict] = {}
    for target, hit in best_hits_by_target(parse_domtblout(domtbl)).items():
        meta = parse_uniprot_header(f"{target} {hit['description']}")
        acc = acc_key(meta["accession"])
        prev = out.get(acc)
        if prev is None or hit["full_score"] > prev["full_score"]:
            out[acc] = {**hit, **meta}
    return out


def assign(itpr_hits: dict[str, dict], ryr_hits: dict[str, dict],
           rel_margin: float = REL_MARGIN,
           min_score: float = MIN_SCORE,
           min_positions: int = MIN_PROFILE_POSITIONS) -> list[dict]:
    """Merge the two profiles' hits into one assignment per accession."""
    rows: list[dict] = []
    for acc in sorted(set(itpr_hits) | set(ryr_hits)):
        i = itpr_hits.get(acc)
        r = ryr_hits.get(acc)
        meta = i or r
        i_score = i["full_score"] if i else 0.0
        r_score = r["full_score"] if r else 0.0
        i_pos, r_pos = _positions(i), _positions(r)
        win = max(i_score, r_score)
        lose = min(i_score, r_score)
        margin = win - lose
        rel = round(margin / win, 4) if win > 0 else 0.0
        # The winner is picked as a *record*, not by score alone. A hit whose
        # full-sequence bit score is exactly 0.0 with no hit from the other
        # profile made `i_score > r_score` false and sent `winner` to the
        # profile that has no record at all — which S23b hit the first time
        # these profiles were run over genomic gene models, where marginal
        # alignments at or below 0 bits are normal (the same population
        # `s20_test_sensitivity.py` found). Ties still resolve to `unassigned`
        # through the margin test below, so no existing call changes.
        if i is not None and (r is None or i_score >= r_score):
            winner, win_rec, win_pos = "ITPR", i, i_pos
        else:
            winner, win_rec, win_pos = "RYR", r, r_pos
        win_cov = win_rec["hmm_coverage"] if win_rec else 0.0
        evidence = ("architecture" if win_cov >= COMPLETE_COVERAGE
                    else "partial" if win_pos >= min_positions
                    else "module")

        if win < min_score:
            call, conf, why = "unassigned", "none", \
                f"best profile score {win:.1f} < {min_score:.0f} bits"
        elif win_pos < min_positions:
            call, conf, why = "unassigned", "none", (
                f"{winner.lower()}.hmm matched only {win_pos} profile "
                f"positions (< {min_positions}) — a shared module, not a "
                "family assignment")
        elif rel < rel_margin:
            call, conf, why = "unassigned", "none", (
                f"profiles within {rel:.1%} of each other "
                f"(< {rel_margin:.0%} no-call band, D7)")
        else:
            call = winner
            conf = ("high" if rel >= 0.30 and evidence == "architecture"
                    else "medium")
            why = (f"{call.lower()}.hmm wins by {margin:.1f} bits "
                   f"({rel:.1%} of {win:.1f}) over {win_pos} profile "
                   "positions")
            if lose == 0.0:
                why += f"; no {'ryr' if call == 'ITPR' else 'itpr'}.hmm hit"

        rows.append({
            "accession": acc, "assignment": call, "confidence": conf,
            "evidence": evidence, "reason": why,
            "itpr_score": round(i_score, 1), "ryr_score": round(r_score, 1),
            "score_margin": round(margin, 1), "rel_margin": rel,
            "itpr_evalue": f"{i['full_evalue']:.3g}" if i else "",
            "ryr_evalue": f"{r['full_evalue']:.3g}" if r else "",
            "best_evalue": f"{min(h['full_evalue'] for h in (i, r) if h):.3g}",
            "itpr_coverage": i["hmm_coverage"] if i else 0.0,
            "ryr_coverage": r["hmm_coverage"] if r else 0.0,
            "itpr_positions": i_pos, "ryr_positions": r_pos,
            "target_coverage": win_rec["target_coverage"] if win_rec else 0.0,
            "length": meta["tlen"], "gene": meta["gene"],
            "species": meta["species"],
            "taxon_id": meta["taxon_id"] or "",
            "reviewed": int(bool(meta["reviewed"])),
            "protein_name": meta["protein_name"].replace("\t", " "),
        })
    return rows


def audit_against(rows: list[dict], truth: dict[str, str],
                  label: str) -> dict:
    """Score the profile assignment against a call it never sees.

    `truth` maps accession → ITPR / RYR / anything else (ignored). Returns
    the confusion counts plus the disagreeing accessions, which is what
    makes the assignment auditable rather than asserted (D14b's discipline).
    """
    counts = {"agree": 0, "disagree": 0, "profile_unassigned": 0,
              "truth_unassigned": 0, "not_in_truth": 0}
    disagreements: list[dict] = []
    for row in rows:
        t = truth.get(row["accession"])
        if t is None:
            counts["not_in_truth"] += 1
            continue
        if t not in ("ITPR", "RYR"):
            counts["truth_unassigned"] += 1
            continue
        if row["assignment"] == "unassigned":
            counts["profile_unassigned"] += 1
        elif row["assignment"] == t:
            counts["agree"] += 1
        else:
            counts["disagree"] += 1
            disagreements.append({**row, "truth": t, "truth_source": label})
    total = counts["agree"] + counts["disagree"]
    counts["agreement"] = round(counts["agree"] / total, 4) if total else None
    return {"label": label, "counts": counts, "disagreements": disagreements}
