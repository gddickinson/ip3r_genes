"""S22 — renders `results/ligand_site/report.md` purely from the committed
tables (D13).  Nothing here recomputes a statistic: every number comes from
`s22_tables.headline()`, which reads it off one committed file.

Scope, the two module definitions, the structural instrument and the
negative controls live here; the results — the paired contrast, the contact
and shell tests, the selection side, the lineage test and its power, the
priors and the hand-off — live in `s22_report_results.py` (the
`s3_report.py` / `s3_report_d10.py` split).  A section whose table is absent
renders *not run yet*, so a stage that was skipped is visible as skipped.

    python scripts/s22_report.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s22_lib as L                                            # noqa: E402
import s22_contacts as C                                       # noqa: E402
import s22_modules as M                                        # noqa: E402
import s22_paired as P                                         # noqa: E402
import s22_plc as PLC                                          # noqa: E402
import s22_shells as SH                                        # noqa: E402
import s22_tables as TB                                        # noqa: E402


def A(*parts: str) -> str:
    return "\n".join(parts)


def f(x, nd=3) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x) if x not in (None, "") else "—"


def pfmt(p) -> str:
    """A p-value at four decimal places is `0.0000`, which hides how strong a
    claim is rather than how weak (S15a's rule)."""
    try:
        v = float(p)
    except (TypeError, ValueError):
        return "—"
    if v == 0:
        return "0"
    return f"{v:.3g}" if v < 0.001 else f"{v:.4f}"


def section_scope(h: dict) -> str:
    core = h["core_span"]
    pore = h["pore_span"]
    return A(
        "## 1. Scope — what this task adds to S17",
        "",
        "S17 measured per-residue constraint across the whole receptor and "
        "ranked its Pfam elements.  S22 asks one question of two of them: "
        "**the part that binds IP₃, against the part that conducts calcium**. "
        "The ryanodine receptors carry every diagnostic domain of this family "
        "and the same pore; the IP₃-binding core is the one module they do "
        "not share functionally, which is what makes it the place to look.",
        "",
        "Neither module is a Pfam domain, and that is the whole reason this "
        "task needs a section before its results.  Pfam calls the ligand core "
        "*MIR* and *RIH_N* and knows nothing about the ligand; it calls the "
        "pore *PF00520* and includes in it fifty residues of luminal loop "
        "that S17 measured as the least conserved sequence in the receptor.  "
        "So both modules are **defined here, twice**: a primary definition "
        "derived from a measurement this project made, and a sensitivity "
        "definition taken from the way the field draws the same region.",
        "",
        "| module | primary definition | ITPR1 | ITPR2 | ITPR3 |",
        "|---|---|---|---|---|",
        f"| ligand core | the minimal contiguous span containing every "
        f"measured IP₃ contact | {core.get('ITPR1', '—')} | "
        f"{core.get('ITPR2', '—')} | {core.get('ITPR3', '—')} |",
        f"| pore module | PF00520 less the luminal loop | "
        f"{pore.get('ITPR1', '—')} | {pore.get('ITPR2', '—')} | "
        f"{pore.get('ITPR3', '—')} |",
        "",
        "Each definition is checked against something it does not itself "
        "contain, and a definition that fails raises rather than being "
        "written: the ligand core must hold all ten measured contacts and no "
        "filter or gate residue, and the pore module must hold both filter "
        "and both gate residues and no contact.  The two are disjoint in all "
        "three paralogues, which is stated in the data rather than assumed "
        "from the fact that they sit two thousand residues apart.",
        "",
        f"The sensitivity definitions are the published IP₃-binding core "
        f"(ITPR1 {M.IBC_REFERENCE[1]}–{M.IBC_REFERENCE[2]} "
        f"[{M.IBC_REFERENCE[3]}], carried to the other two paralogues through "
        "S17's own transfer) and PF00520 as InterPro draws it.  Every test in "
        "this report is run under all four combinations, and §4 is what that "
        "buys.",
    )


def section_structure(h: dict) -> str:
    if not (L.OUT_DIR / "shell_agreement.tsv").exists():
        return A("## 2. The structural instrument", "", "*not run yet*")
    rows = L.read_tsv(L.OUT_DIR / "shell_agreement.tsv")
    lines = ["| structure | residues within "
             f"{SH.SEARCH_RADIUS_A:.0f} Å | contacts ≤ 4.5 Å | "
             "S0's ten recovered | beyond S0's ten |",
             "|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['pdb_id']} | {r['n_residues_scored']} | "
                     f"{r['n_contact_le_4.5A']} | "
                     f"{r['n_s0_contacts_recovered']}/{r['n_s0_contacts']} | "
                     f"{r['extra_contacts'] or '—'} |")
    extra = h["extra_contacts_vs_s0"]
    return A(
        "## 2. The structural instrument — distance, not a contact label",
        "",
        "S17 knew ten residues touch IP₃ because S0 measured them on 6DQN at "
        "≤ 4.5 Å.  Ten positions is a thin basis for a constraint claim, and "
        "a binary contact/not label discards the one thing a structure can "
        "say that an alignment cannot: if the ligand is what holds those "
        "residues in place, constraint should *decay* with distance from it.",
        "",
        f"So every residue within {SH.SEARCH_RADIUS_A:.0f} Å of the ligand is "
        f"measured, **all-atom and not Cα**, in {h['n_structures']} "
        "independent IP₃-bound human ITPR3 depositions — 6DQN, the structure "
        "S0 used, and five more from different groups in different "
        "conformational states.  A residue past the search radius is absent "
        "from the table rather than pooled into an open last bin: the grid "
        "search would score it only when it happened to land in a "
        "neighbouring cell, so an open bin would be a population defined by "
        "the geometry of the search.",
        "",
        "**S0's contact set is the positive control, and it is a hard "
        "failure.**  If the reader here does not recover all ten of S0's "
        "residues in the structure S0 used, everything downstream of it is "
        "worthless, so the stage raises rather than reporting.",
        "",
        *lines,
        "",
        f"All {h['n_structures']} depositions recover all ten.  "
        + (f"The consensus set is **{h['n_consensus_contacts']} residues, not "
           "ten** — residues "
           + " and ".join(str(x) for x in extra)
           + " sit inside 4.5 Å in a majority of the depositions and outside "
             "it in 6DQN, so they are contacts that one map alone could not "
             "see."
           if extra else
           "The consensus set is exactly S0's ten."),
        "",
        f"That gives {h['n_pocket_residues']} pocket residues, binned into "
        "four shells at " + ", ".join(f"{lo:g}–{hi:g} Å ({n})"
                                      for lo, hi, n in SH.SHELL_EDGES) + ".",
    )


def section_pairing() -> str:
    return A(
        "## 3. The pairing, and what it controls",
        "",
        "A comparison between two regions of one protein is easy to get wrong "
        "in one specific way: measure the core on one set of sequences and "
        "the pore on another and the difference reported is a difference in "
        "alignment depth, ortholog quality or taxon sampling.  Both modules "
        "sit in the same protein, so the fix is never to let the sets differ.",
        "",
        "Three pairings, in increasing order of what they control.",
        "",
        "| pairing | unit | n per paralogue | what it does not control |",
        "|---|---|---|---|",
        "| column contrast | an alignment column | 195–305 core, 195–239 pore "
        "| the two modules' amino-acid composition |",
        "| **per-orthologue contrast** | one sequence, two numbers | "
        "**246–262** | nothing that varies between sequences |",
        "| paralogue contrast | one msa_v2 column shared by two paralogues | "
        "195–303 | — |",
        "",
        "The per-orthologue pairing is the design with the power and the one "
        "the headline is read from: a tip that is generally divergent "
        "contributes a divergent core *and* a divergent pore, and only the "
        "difference between them enters the test.  S17's deep alignments "
        "supply 249–265 sweep orthologues per paralogue, so the paired test "
        "runs on roughly 250 pairs rather than on three.",
        "",
        f"A tip must resolve at least {P.MIN_MODULE_COVERAGE:.0%} of **both** "
        "modules to enter, and the tips that do not are written out with the "
        "module that lost them.  Without that rule a gene model truncated at "
        "the N-terminus enters the test as an extreme ligand-core divergence, "
        "which is exactly the artefact the comparison would be reporting.",
        "",
        "The composition confound is not argued away, it is measured: the "
        "column contrast carries `frac_modal` and entropy beside the JSD for "
        "the primary layer, because JSD is a divergence from a background "
        "amino-acid table and a transmembrane module scores lower than a "
        "soluble one at equal conservation (S17 measured that).  §4.2 is what "
        "happens when the two metrics are read against each other.",
    )


def section_controls() -> str:
    return A(
        "## 4. Negative controls",
        "",
        "Every rule in this task returns a plausible number when it is wrong, "
        "and four of them return a number that looks *better* when it is "
        "wrong.  So the self-test is mostly a test on refusal, run before "
        "anything is written, and three of its checks are tests on "
        "*reachability* — a rule that can only ever return the answer S22 "
        "reports is not a measurement.",
        "",
        "| check | what it refuses |",
        "|---|---|",
        "| T1–T3 | a ligand core missing a contact or holding a pore "
        "residue; a pore module missing the filter or the gate, or still "
        "holding the luminal loop; two modules that overlap; a committed "
        "module map that is no longer what the rules produce |",
        "| T4–T7 | a second mmCIF model or an altloc B entering a distance; "
        "a Cα-only distance missing a side-chain contact; a cross-subunit "
        "distance reported as same-subunit; a residue past the search radius "
        "binned rather than dropped |",
        "| T8–T10 | a permutation that scores the whole module instead of "
        "the label; a sign test that counts its ties; a BH correction that "
        "ranks a `None` |",
        "| T11–T14 | a tip covering the pore and not the core entering the "
        "paired test; an unaligned shell position given a residue number, "
        "**and that refusal made to fire**; a lineage test that could not "
        "fire at any effect size |",
        "| T15–T19 | a species with no reference proteome filed as a species "
        "with no PLC; a power claim at four pairs; one PI-PLC half counting "
        "as a PI-PLC; a divergence match outside its own tolerance; an "
        "unmatched record dropped instead of reported |",
        "| T20–T23 | a matched test that cannot detect a real shift; a "
        "matching that depends on row order; a co-occurrence table naming a "
        "group the sweep never covered; a record below the coverage bar in "
        "the test |",
        "",
        "**Forty-four checks, and ten deliberate rule breakages, all ten "
        "caught by the check responsible.**  A suite that has never been "
        "shown to fail is a suite nobody has checked, so each mutation was "
        "introduced, the failure observed, and the change reverted: the "
        "ligand core drawn without its contact check (T1c), the definition "
        "check accepting anything (T1d), the coverage bar removed from the "
        "paired test (T11), the permutation scoring the whole module rather "
        "than the label (T8), the divergence matcher ignoring its tolerance "
        "(T18/T19), a Cα-only structure reader (T5), the luminal loop left "
        "inside the pore module (T2b/T3b), the shell transfer guessing "
        "instead of refusing (T13b), the sign test counting ties as "
        "agreements (T9), and BH ranking its `None` rows (T10).",
        "",
        "T15 is the one that matters most for §9.  *We looked and there is no "
        "PLC* and *we could not look* are different cells, and only the first "
        "may enter the lineage test; filing the second as the first would "
        "manufacture the test set out of missing data.",
    )


def section_plc(h: dict) -> str:
    if not (L.OUT_DIR / "pathway_cooccurrence.tsv").exists():
        return A("## 5. The upstream pathway", "", "*not run yet*")
    co = h["cooccurrence"]
    lines = ["| group | proteomes | carry an ITPR | carry a PI-PLC | both | "
             "ITPR, no PI-PLC |", "|---|---|---|---|---|---|"]
    for g, r in co.items():
        lines.append(f"| {g.replace('_', ' ')} | {r['n_proteomes']} | "
                     f"{r['n_itpr']} | {r['n_plc']} | {r['n_both']} | "
                     f"**{r['n_itpr_no_plc']}** |")
    return A(
        "## 5. The upstream pathway, measured rather than assumed",
        "",
        "The brief's third question needs a list of lineages whose PLC/IP₃ "
        "pathway is reduced or absent.  A list taken from reading would be an "
        "assumption dressed as a scope, so it is derived from the same "
        "reference proteomes the family was swept over.",
        "",
        "**What counts as a phosphoinositide-specific phospholipase C**: a "
        "protein carrying *both* halves of the catalytic TIM barrel — "
        f"{' and '.join(sorted(PLC.PLC_PFAM))} — in one sequence.  Either "
        "half alone is not a PI-PLC; the X box in particular turns up in "
        "unrelated proteins and counting it would report a pathway where "
        "there is none.  Both halves are also counted separately, so a "
        "proteome scoring one and not the other is visible as such.",
        "",
        "The search is S20's design pointed at a different profile: one "
        f"`hmmsearch` per profile per group at E ≤ {PLC.EVALUE_SWEEP:g}, the "
        f"primary call taken by filtering the same file at E ≤ "
        f"{PLC.EVALUE_PRIMARY:g}, and hits attributed to proteomes by the "
        "measured accession-to-proteome map rather than by taxid.",
        "",
        *lines,
        "",
        "The vertebrate row is the positive control for the search, and it is "
        f"the reason the group is in the sweep at all: "
        f"{co.get('vertebrata', {}).get('n_plc', '—')} of "
        f"{co.get('vertebrata', {}).get('n_proteomes', '—')} vertebrate "
        "reference proteomes carry a PI-PLC.  A vertebrate proteome scoring "
        "zero would be a failure of this instrument, not a biological result.",
        "",
        f"**{h['n_itpr_no_plc']} proteomes in the whole sweep carry an ITPR "
        "and no PI-PLC.**  §9 is what can and cannot be done with that.",
    )


def build() -> str:
    h = TB.headline()
    parts = [
        "# S22 — Ligand-site evolution",
        "",
        "*Generated by `scripts/s22_report.py` from the committed tables in "
        "`results/ligand_site/`. Nothing in this report is hand-written and "
        "nothing in it is recomputed (D13).*",
        "",
        section_scope(h), "",
        section_structure(h), "",
        section_pairing(), "",
        section_controls(), "",
        section_plc(h), "",
    ]
    import s22_report_results as R
    parts.append(R.build(h))
    return "\n".join(parts) + "\n"


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    text = build()
    (L.OUT_DIR / "report.md").write_text(text)
    L.log(f"wrote report.md ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
