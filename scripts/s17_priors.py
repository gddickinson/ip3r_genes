"""The priors S17's results are judged against, in one module imported by every
half of the report so a prior cannot be stated twice and drift.

The verdict vocabulary is S8's five-valued one, imported from `s15_priors`
unchanged rather than restated. Three of its values do real work here:

* **`orthogonal`** — a different property of the same object. S9 measured a
  *rate* over a vertebrate tree; S17 measures the *dispersion of residues* in a
  column of 250 orthologues. Both are called "constraint" in prose and they are
  not the same quantity: a site can be invariant across sampled vertebrates and
  still sit in a clade whose ω is high, and a paralog can have the lowest ω and
  not the highest column conservation. Reading a mismatch between them as a
  contradiction would claim a disagreement between two measurements that were
  never of the same thing.
* **`underpowered`** — the measurement did not reach the question. S17 has two
  places this can happen honestly: a paralog whose ClinVar record is too thin
  to score a classifier on, and a functional site set of ten residues, where a
  null result is a statement about ten numbers and not about the site.
* **`not corroborated`** — a related property measured cleanly that does not
  support the prior. Reserved for the case where a *different* instrument
  (identity between paralogs, FEL's per-site rate) is measured properly and
  simply does not point the prior's way.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s15_priors import VERDICTS                                # noqa: E402,F401

#: key -> (what an earlier task concluded, where it said it)
PRIOR = {
    "omega_ranking": (
        "ITPR1 is the most constrained paralog and ITPR2 the least — "
        "one-ratio omega 0.0238 (ITPR1) < 0.0415 (ITPR3) < 0.0430 (ITPR2), "
        "and the two-ratio test puts each paralog clade's own omega "
        "significantly apart from the rest of the family (q <= 2.1e-18)",
        "S9 §4.1/§4.3"),
    "no_positive_selection": (
        "no site in any paralog is positively selected — all three M8-vs-M7 "
        "tests are significant after BH and in every one the extra class "
        "sits exactly at codeml's omega = 1 boundary carrying under 1 % of "
        "sites, with 0, 0 and 1 sites reaching a 0.95 BEB posterior",
        "S9 §4.5"),
    "itpr1_stem_burst": (
        "the ITPR1 stem is the one branch-site result that stands — omega2 "
        "= 5.53 on 11.0 % of sites, stable across four restarts, 8 sites at "
        "BEB >= 0.95 — while the ITPR2 and ITPR3 stems' omega2 is pinned at "
        "codeml's 999 bound and is not an estimate",
        "S9 §4.4"),
    "within_paralog_identity": (
        "mean covered-only identity within a vertebrate paralog group is "
        "0.910, and between groups 0.791 (ITPR1xITPR2), 0.753 "
        "(ITPR2xITPR3), 0.741 (ITPR1xITPR3)",
        "S6 §4.1/§4.2"),
    "saturation": (
        "synonymous sites are saturated *within* a single paralog, not only "
        "between them — 94.2 % of within-ITPR1 pairs exceed dS = 1.5 — so "
        "every ratio in S9 is a tree-based estimate and no pairwise omega "
        "is quoted as one",
        "S9 §4.2"),
    "arg2524cys": (
        "the recurrent multisystem variant ITPR3 p.Arg2524Cys sits 7 "
        "residues past the measured gate (Phe2513/Ile2517 in 6DQN), on the "
        "C-terminal stretch after the pore",
        "S0 review figures, emergent 2026-08-19"),
    "thin_variant_literature": (
        "the literature variant set is thinner than the review reads — of "
        "the nine disease entries the review catalogues, only two carry a "
        "residue the cited source names, both in ITPR3 (p.Thr1424Met and "
        "p.Arg2524Cys); the rest localise to a domain or to the gene",
        "S0 review figures, emergent 2026-08-19"),
    "itpr3_multisystem": (
        "ITPR3's clinical phenotype is broader than neuropathy — the "
        "recurrent de novo p.Arg2524Cys causes a multisystemic disease with "
        "immunodeficiency, so the ITPR3 variant set is not a CMT1J set",
        "S0, emergent 2026-08-18"),
    "afdb_absent": (
        "the family is effectively absent from AlphaFold DB — 13 of 5,861 "
        "census records at or above 2,000 aa have a usable full-length "
        "model, and the human ITPR2 accession returns a 181-residue isoform",
        "S11 §3, emergent 2026-09-08"),
    "ryr_shares_the_domains": (
        "the ryanodine receptors carry every ITPR-diagnostic Pfam domain "
        "(PF08709, PF02815, PF01365, PF08454) and are inside every search "
        "this project runs",
        "D14 / docs/ip3r_background.md"),
    "itpr3_lesion_excess": (
        "ITPR3 carries a disabling-indel excess against sibling family loci "
        "matched on bait identity within the same genome, localised to Aves "
        "(25:2, q = 4.5e-5) — pointing the opposite way to every constraint "
        "result in the project",
        "S15a D47 / S15b §9"),
}


def prior(key: str) -> tuple[str, str]:
    return PRIOR[key]


def line(key: str, verdict: str, detail: str) -> str:
    """One rendered prior line: what was said, where, what S17 measured.

    The verdict is checked against the committed vocabulary rather than
    formatted freely, so a typo cannot invent a sixth verdict.
    """
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not one of {VERDICTS}")
    what, where = PRIOR[key]
    return f"Prior: **{what}** — {where}. **{verdict}** — {detail}"
