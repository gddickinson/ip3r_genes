"""The priors S22's results are judged against, in one module imported by
every half of the report so a prior cannot be stated twice and drift.

The verdict vocabulary is S8's five-valued one, imported from `s15_priors`
unchanged.  Two of its values carry most of the weight here.

* **`orthogonal`** — a different property of the same object.  S22 measures
  two quantities that prose calls the same thing: a column's dispersion
  across 250 orthologues of one paralogue, and a site's non-synonymous rate
  on a 57-tip vertebrate tree.  They can point different ways at the same
  residue without either being wrong, and S17 already reserved this verdict
  for exactly that pair.  It also covers the comparison a reader will most
  want to make and should not: S20's *absence* of the receptor in land
  plants is about whether a gene exists, and S22's pathway table is about
  what sits beside the ones that do.
* **`underpowered`** — the measurement did not reach the question.  The
  brief said in advance that the lineage test is easy to underpower and
  that saying so is a result, so this verdict is not a consolation here;
  it is one of the outcomes the design was built to be able to report, and
  `lineage_power.tsv` is what licenses it.

Two priors are imported from `s17_priors` rather than restated, because
they are the same sentence about the same tables and a second copy could
drift from the first.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s15_priors import VERDICTS                                # noqa: E402,F401
from s17_priors import PRIOR as S17_PRIOR                      # noqa: E402

#: key -> (what an earlier task concluded, where it said it)
PRIOR = {
    "contacts_vs_own_element": (
        "the ten IP3 contacts are more constrained than the Pfam element "
        "each sits in — mean JSD 0.813 / 0.840 / 0.824 against 0.763 / "
        "0.767 / 0.770 for MIR and RIH_N, p = 0.070 / 0.0021 / 0.025",
        "S17 §6.3, functional_site_constraint.tsv"),
    "pore_most_constrained": (
        "the gate and the selectivity filter are the most constrained "
        "elements of the protein on both metrics and in all three "
        "paralogues, and with the luminal loop separated out the channel "
        "domain itself sits at JSD 0.753 / 0.750 / 0.749 above a linker "
        "mean of 0.734",
        "S17 §5"),
    "luminal_loop_least_conserved": (
        "the luminal loop is by a wide margin the least conserved element "
        "in the receptor — JSD 0.488 / 0.452 / 0.485 — and it lives about "
        "fifty residues from the gate inside the same Pfam domain",
        "S17 §5"),
    "ligand_core_names_the_family": (
        "the signature that names the family is the IP3-binding core "
        "(PF08709, Ins145_P3_rec), and it is the one diagnostic domain the "
        "ryanodine receptors do not carry functionally",
        "docs/ip3r_background.md; D14"),
    "ryr_shares_the_domains": S17_PRIOR["ryr_shares_the_domains"],
    "omega_ranking": S17_PRIOR["omega_ranking"],
    "family_absent_in_plants_and_dikarya": (
        "no land plant and no dikaryan fungal reference proteome carries an "
        "ITPR — 15 of 432 Viridiplantae and 28 of 1,527 fungal proteomes "
        "score a family hit at all, and the plant and fungal records that "
        "survive being chased are chlorophyte algae and early-diverging "
        "fungi, not embryophytes or Dikarya",
        "S20 §4, S23 §5"),
    "family_identity_to_ryr": (
        "ITPR and RyR are 0.249 identical over mutually covered columns "
        "while within-family identity is 0.828 — the two families share "
        "every diagnostic Pfam and almost no sequence",
        "S1 benchmark, confirmed on the S6 alignment (S6 §4.1)"),
}


def prior(key: str) -> tuple[str, str]:
    return PRIOR[key]


def line(key: str, verdict: str, detail: str) -> str:
    """One rendered prior line: what was said, where, what S22 measured."""
    if verdict not in VERDICTS:
        raise ValueError(f"{verdict!r} is not one of {VERDICTS}")
    what, where = PRIOR[key]
    return f"Prior: **{what}** — {where}. **{verdict}** — {detail}"
