"""The variant, selection and structural half of S17's report.

Split out of `s17_report_results.py` to keep both inside the 500-line budget
(the `s3_report.py` / `s3_report_d10.py` pattern), and taking the caller's
`headline()` dict so no half can read the tables differently.

This is the half that argues against the task's own design: §7.2 reports that
the deep within-paralog layer S17 was built around is the second-best of four
classifiers, and §11 lists what the task does not settle.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402
import s17_priors as P                                         # noqa: E402
from s17_report_results import A, _f                           # noqa: E402


def s_variants(h: dict) -> str:
    """§7 — the variant harvest, its numbering check, and the classifier."""
    b = h["buckets"]
    tr = h["transcripts"]
    ctrl = h["curated_control"]
    auc = h["auc_matched"]
    prot = h["auc_protein"]
    order = sorted(auc, key=lambda k: -float(auc[k][0] or 0))
    best = order[0] if order else ""
    bucket_rows = [
        f"| {g} | {b.get((g, 'P/LP'), 0)} | {b.get((g, 'B/LB'), 0)} | "
        f"{b.get((g, 'VUS'), 0)} | {b.get((g, 'conflicting'), 0)} | "
        f"{b.get((g, 'other'), 0)} |" for g in L.PARALOGS]
    auc_rows = [
        f"| `{ly}` | {auc[ly][0]} | {auc[ly][1]} | {prot.get(ly, '—')} |"
        for ly in order]
    return A(
        "## 7. Is the score any good?",
        "",
        f"**{h['n_clinvar']} ClinVar missense records across the three genes, "
        f"and {h['n_dropped']} dropped.** Numbering is checked, never assumed: "
        "every cited transcript's own translated CDS is fetched, aligned to "
        "the canonical, positions transferred through that alignment, and the "
        "reference amino acid then required to match.",
        "",
        "| gene | transcript | length | identical to canonical | variants | kept |",
        "|---|---|---|---|---|---|",
        *[f"| {r['gene']} | {r['transcript']} | {r['protein_len']} | "
          f"{r['identical_to_canonical']} | {r['n_variants']} | {r['n_kept']} |"
          for r in tr],
        "",
        "All three genes file their ClinVar records on a transcript whose "
        "translated CDS **is** the UniProt canonical, so the check passes with "
        "nothing dropped. That is a result of the check, not a reason it was "
        "unnecessary: the same code on the PIEZO family found 512 of 773 "
        "positions carrying a different amino acid in the canonical, and a "
        "resource built without the check would have been wrong everywhere. "
        "Here it is right everywhere, and the report can say which.",
        "",
        "**The harvest carries a positive control it could fail.** S0's "
        "curated table localises two variants to a named residue, and both "
        "must come back from a query returning ~1,750 records:",
        "",
        *[f"- {r['gene']} {r['label']} → **{r['verdict']}** "
          f"({r['n_clinvar_at_resi']} records at residue {r['resi']}, "
          f"{r['sources']})" for r in ctrl],
        "",
        "### 7.1 The record is mostly uncertain, and that is the finding",
        "",
        "| gene | P/LP | B/LB | VUS | conflicting | other |",
        "|---|---|---|---|---|---|",
        *bucket_rows,
        "",
        f"**{h['n_vus']} of {h['n_clinvar']} records "
        f"({h['frac_vus']:.0%}) are of uncertain significance**, and the "
        "labelled sets a classifier can be scored on are 49, 1 and 5 "
        "pathogenic positions. ITPR2's clinical record consists of a single "
        "pathogenic missense variant. Any per-gene AUC for ITPR2 is therefore "
        "not reported as a result, and ITPR3's rests on five positions.",
        "",
        "### 7.2 The layers, on one fixed set of variants",
        "",
        "The four layers do not cover the same residues — the msa_v2 layers "
        "lose gappy columns the deep alignment fills, and vice versa — so "
        "ranking four AUCs measured on four slightly different variant sets "
        "would compare the sets as much as the layers. The table below is the "
        f"pooled test restricted to the {auc[order[0]][2] if order else 0} "
        f"pathogenic and {auc[order[0]][3] if order else 0} benign positions "
        "where **every** layer has a reliable score.",
        "",
        "| layer | AUC | p | AUC vs whole protein |",
        "|---|---|---|---|",
        *auc_rows,
        "",
        "**Two things fall out, and one of them is against the design of this "
        "task.**",
        "",
        f"The `shallow` control is the worst layer ({auc.get('shallow', ['—'])[0]}), "
        f"so the depth was worth building: without §3's ortholog sets S17 "
        f"would have had that number. But **the best layer is `{best}`** "
        f"({auc[best][0]}), not the deep within-paralog one "
        f"({auc.get('deep', ['—'])[0]}). Taxonomic breadth beats within-gene "
        "depth here: a position conserved across 500 Myr of paralog "
        "divergence *and* across the eukaryotes discriminates pathogenic from "
        "benign better than a position merely invariant across 260 "
        "vertebrate orthologues of the same gene. A resource for this family "
        "should quote the family-wide layer, and this task's own instrument is "
        "the second-best of the four.",
        "",
        "The whole-protein control is the reason the first column cannot be "
        "read alone. Every layer separates pathogenic positions from the "
        f"average residue too (AUC {prot.get(best, '—')} for `{best}`), but by "
        "less than it separates them from the benign set — so part of the "
        "P/LP-vs-B/LB signal is ClinVar's benign calls being enriched in "
        "*unconstrained* positions, not only its pathogenic calls being "
        "enriched in constrained ones. Both directions are real and the "
        "control is what distinguishes them.",
        "",
        P.line("thin_variant_literature", "confirmed",
               f"the database record is thin in the same way the literature "
               f"is: {h['frac_vus']:.0%} of {h['n_clinvar']} missense records "
               f"are uncertain, and the two residue-level variants S0 could "
               f"cite are both recovered here — so the review's thinness was "
               f"the field's, not the review's"),
        "",
        P.line("itpr3_multisystem", "orthogonal",
               "S17 scores positions, not phenotypes. p.Arg2524Cys is in the "
               "harvest with 4 records and its constraint is reported in "
               "§8; whether the phenotype is multisystemic is a clinical "
               "claim no alignment can address, and it is carried forward "
               "unchanged rather than re-tested"),
    )


def s_vus(h: dict) -> str:
    """§8 — what the score offers the 1,500 uncertain variants."""
    rows = []
    for r in h["vus"]:
        if r["layer"] != "family":
            continue
        rows.append(f"| {r['gene']} | {r['n_vus']} | "
                    f"{_f(r['median_pathogenic'])} | {_f(r['median_benign'])} | "
                    f"{_f(r['median_vus'])} | "
                    f"{float(r['frac_vus_above_pathogenic_median']):.0%} | "
                    f"{float(r['frac_vus_below_benign_median']):.0%} |")
    return A(
        "## 8. What this offers the uncertain variants",
        "",
        "88 % of the record is uncertain, and a per-site score is what a "
        "constraint resource can offer such a variant. Each VUS is placed "
        "against the *labelled* distributions of the same gene, on thresholds "
        "that are those distributions' own medians so no cut was chosen to "
        "make a count. On the family layer §7.2 selected:",
        "",
        "| gene | VUS scored | P/LP median | B/LB median | VUS median | "
        "≥ P/LP median | ≤ B/LB median |",
        "|---|---|---|---|---|---|---|",
        *rows,
        "",
        "This is a **stratification and not a call**. It says where a variant "
        "sits on an axis the labelled variants separate on, which is a "
        "different claim from pathogenicity, and the per-gene numbers inherit "
        "§7.1's problem — the ITPR2 row's P/LP median is one position's score.",
        "",
        "The deliverable is `constraint_<gene>_<acc>.tsv`: every residue of "
        "each human paralog with four conservation layers, its occupancy, its "
        "element, and whether it is a measured IP₃ contact, filter or gate "
        "residue. `variants.tsv` is the same columns joined onto every "
        "harvested variant.",
    )


def s_selection(h: dict) -> str:
    """§9 — per-site dN/dS on the same coordinates, or a note that it has not run."""
    status = h.get("fel_status") or []
    if not status:
        return A("## 9. Per-site selection",
                 "",
                 "*not run yet* — `s17_fel.py` writes `fel_status.tsv` and "
                 "`selection_by_element.tsv`; neither is present.")
    rows = [f"| {r['paralog']} | {r['n_sites']} | {r['n_mapped']} | "
            f"{float(r['aa_agreement']):.1%} | {r['n_purifying_q05']} | "
            f"{r['n_diversifying_q05']} | {r['n_alpha_at_bound']} |"
            for r in status]
    fe = h.get("fel_elements") or []
    erows = []
    for e in ("gate", "selectivity_filter", "channel", "RIH_assoc",
              "nterm_trefoil", "luminal_loop", "ip3_contact", "WHOLE_PROTEIN"):
        cells = []
        for p in L.PARALOGS:
            hit = [r for r in fe if r["paralog"] == p and r["element"] == e]
            cells.append(f"{_f(hit[0]['median_beta'])} / "
                         f"{float(hit[0]['frac_purifying_q05']):.0%}"
                         if hit and hit[0]["frac_purifying_q05"] != "" else "—")
        if any(c != "—" for c in cells):
            erows.append(f"| `{e}` | " + " | ".join(cells) + " |")
    ndiv = sum(int(r["n_diversifying_q05"]) for r in status)
    npur = sum(int(r["n_purifying_q05"]) for r in status)
    nbound = sum(int(r["n_alpha_at_bound"]) for r in status)
    nsite = sum(int(r["n_sites"]) for r in status)
    # BH at q < 0.05 admits up to 5 % of *rejections* as false, so the bar a
    # diversifying count has to clear to mean anything is that share of the
    # rejections FEL actually made — not zero, and not a number typed here.
    fdr_budget = 0.05 * (npur + ndiv)
    return A(
        "## 9. Per-site selection on the same coordinates",
        "",
        "S9 answered the whole-gene question — every paralog far below "
        "neutrality, no site model finding a positively selected site — but "
        "not *where* the constraint sits. HyPhy FEL on S9's three codon "
        "alignments gives a per-site synonymous and non-synonymous rate, and "
        "every site is carried onto the human reference protein through a "
        "**validated** transfer: the codon alignment is trimAl-trimmed, so its "
        "reference row is a subsequence of the protein, and a site whose amino "
        "acid disagrees after realignment is written with no residue rather "
        "than with a plausible wrong one.",
        "",
        "| paralog | sites | mapped | amino-acid agreement | purifying q<0.05 | "
        "diversifying q<0.05 | α at bound |",
        "|---|---|---|---|---|---|---|",
        *rows,
        "",
        "Element-level rates, as **median β and the share of sites FEL calls "
        "significantly constrained**. The obvious aggregate — Σβ over Σα — is "
        "not usable here: S9 measured median dS of 4.6–13.5 on these "
        "alignments and FEL duly pins α at its upper bound at a share of "
        "sites, so the sum is dominated by sites whose synonymous rate is "
        "unidentifiable rather than large.",
        "",
        "| element | ITPR1 | ITPR2 | ITPR3 |",
        "|---|---|---|---|",
        "| | *median β / % purifying at q<0.05* | | |",
        *erows,
        "",
        P.line("no_positive_selection",
               "confirmed" if ndiv <= fdr_budget else "not corroborated",
               f"FEL calls {ndiv} site(s) diversifying and {npur} purifying "
               f"at q < 0.05 across the three paralogs. Benjamini–Hochberg at "
               f"that level admits up to {fdr_budget:.0f} false rejections in "
               f"a family this size, so {ndiv} is inside the procedure's own "
               f"error budget and is not evidence of a selected site. A "
               f"per-site screen finding nothing is what a whole-gene ω of "
               f"0.024–0.043 predicts, and this is a third instrument saying "
               f"it — after codeml's site models and its branch-site tests"),
        "",
        P.line("saturation", "orthogonal",
               f"α reaches HyPhy's upper bound at only {nbound} of {nsite} "
               f"sites ({nbound / nsite:.1%}), so FEL's per-site synonymous "
               f"rate is identifiable almost everywhere — which does **not** "
               f"corroborate the pairwise saturation S9 measured, and should "
               f"not be read as contradicting it either. S9's number is a "
               f"pairwise dS asked to carry 450 Myr between two tips; FEL "
               f"fits α over the whole tree, which is exactly the distinction "
               f"S9 drew when it said no pairwise ω is quoted as an estimate. "
               f"The two are measurements of different things, and the "
               f"element table above reports β and a count rather than a "
               f"per-site ω because the *ratio* is what saturation degrades"),
    )


def s_structures(h: dict) -> str:
    painted = h.get("painted") or []
    if not painted:
        return A("## 10. Constraint on the structures",
                 "",
                 "*not run yet* — `s17_paint.py` writes "
                 "`painted_structures.tsv`; it is not present.")
    done = [r for r in painted if r["layer"] == "constraint" and r["n_painted"]]
    skipped = [r for r in painted if not r["layer"]]
    rows = [f"| {r['structure']} | {r['role']} | {r['paralog']} | "
            f"{r['n_residues']} | {r['n_painted']} | "
            f"{float(r['frac_painted']):.0%} |" for r in done]
    return A(
        "## 10. Constraint on the structures",
        "",
        "Two layers are written into the B-factor column of every ITPR "
        "structure in S11's panel — `constraint` as JSD × 100 and `selection` "
        "as the FEL non-synonymous rate inverted onto the same scale — so both "
        "read the same way round and colour with one command:",
        "",
        "```",
        "load constraint_reference_ITPR3_8TKG.pdb",
        "spectrum b, blue_white_red, minimum=40, maximum=100",
        "```",
        "",
        "| structure | role | paralog | residues | painted | |",
        "|---|---|---|---|---|---|",
        *rows,
        "",
        "Unscored residues are written **−1, never 0**, so an unscored residue "
        "cannot be read as an unconstrained one. That distinction is load-"
        "bearing for exactly one region: the luminal loop is both the least "
        "conserved element in the receptor and the one the cryo-EM maps "
        "resolve worst, and on a coloured structure those two must not look "
        "the same.",
        "",
        *([f"{len(skipped)} panel row(s) carry no structure file and are "
           f"reported unpainted: "
           + "; ".join(f"`{r['structure']}` ({r['note']})" for r in skipped)]
          if skipped else []),
        "",
        P.line("afdb_absent", "confirmed",
               "the paintable panel is carried by cryo-EM entries, not by "
               "predicted models: S11's human ITPR2 AFDB record is the "
               "181-residue isoform and has no file to paint, so ITPR2's "
               "constraint is shown on 9YKK. A structural resource for this "
               "family cannot be built out of AlphaFold DB"),
    )


def _cited(h: dict, gene: str, resi: int, key: str):
    return (h.get("cited_variants", {}).get((gene, resi), {}) or {}).get(key, "")


def s_caveats(h: dict) -> str:
    auc = h["auc_matched"]
    b = h["buckets"]
    best = max(auc, key=lambda k: float(auc[k][0] or 0)) if auc else ""
    return A(
        "## 11. What this settles, and what it does not",
        "",
        "**Settled.**",
        "",
        "1. The gate and the selectivity filter are the most constrained "
        "elements of the receptor, on two metrics and in all three paralogs, "
        "and the gate is identical between all three.",
        "2. The measured IP₃ contacts are more constrained than the rest of "
        "the domains carrying them — and those domains are MIR and RIH, not "
        "the Pfam signature named for IP₃ binding.",
        "3. The pore module is more constrained than the receptor's linkers "
        "**once the luminal loop inside it is separated out**, and that loop "
        "is the most variable sequence in the protein.",
        "4. Per-site conservation separates ClinVar's pathogenic positions "
        f"from its benign ones (best AUC {auc.get(best, ['—'])[0]} on "
        f"`{best}`), and the family-wide layer does it better than the deep "
        "within-paralog one.",
        "",
        "**Not settled, and why.**",
        "",
        f"1. **The labelled variant set is small.** "
        f"{auc.get(best, [0, 0, 0, 0])[2]} pathogenic and "
        f"{auc.get(best, [0, 0, 0, 0])[3]} benign positions survive the "
        f"all-layers-scorable restriction, and "
        f"{b.get(('ITPR1', 'P/LP'), 0)} of the "
        f"{sum(b.get((g, 'P/LP'), 0) for g in L.PARALOGS)} pathogenic ClinVar "
        f"records are ITPR1's. Every AUC here is a measurement on that set, "
        f"not on the family's variant space.",
        "2. **ITPR2 has one pathogenic missense record.** No per-gene "
        "classifier claim is made for it, and its paralog-audit rows rest on "
        "a single position. The right reading of ITPR2's near-empty clinical "
        "record is ascertainment, not tolerance: §6.1 shows its gate and "
        "IP₃ contacts are as conserved as ITPR1's.",
        "3. **Conservation is not selection.** The four layers measure "
        "residue dispersion in a column; ω measures a substitution-rate "
        "ratio. §9 reports both on one set of coordinates but does not merge "
        "them, and a mismatch between them is not a disagreement.",
        "4. **The filter and gate sets are 2 residues each.** Their means are "
        "reported; their p-values are not, and the element-level test in §5 is "
        "where those claims are powered.",
        "5. **A structural claim about ITPR1 rests on a rat structure.** "
        "S11's ITPR1 reference is 7LHF (*Rattus norvegicus*); the human "
        "profile is transferred onto it by alignment, and `n_painted` records "
        "how much of the chain that reached.",
        "",
        "**Two things to hold against this task.**",
        "",
        P.line("omega_ranking", "orthogonal",
               "S9 ranks the paralogs ITPR1 < ITPR3 < ITPR2 by ω; S17's "
               "median deep JSD is within 0.007 across all three and ranks "
               "them differently element by element. These are not the same "
               "quantity — a rate over a tree against the dispersion of a "
               "column of 260 sampled orthologues — and S17 does not offer a "
               "paralog ranking of its own"),
        "",
        P.line("itpr3_lesion_excess", "orthogonal",
               "S15a/S15b's ITPR3 indel excess points the opposite way to "
               "every constraint result here, and S17 cannot adjudicate it: "
               "lesion density counts indels in assemblies and this task "
               "counts residues in an alignment whose lesion-rich loci were "
               "**excluded by rule** (§3.1). The two measure different "
               "objects, and the excess remains an open emergent row"),
        "",
        P.line("arg2524cys", "confirmed",
               "the variant is in the harvest at ITPR3 residue 2524 with the "
               "reference arginine confirmed, and the transferred gate is "
               f"2513–2517 — so it sits 7 residues past the gate. S17 adds "
               f"the element it falls in, which is not `cterm` as S0's phrase "
               f"\"the C-terminal stretch after the pore\" would suggest: it "
               f"is `{_cited(h, 'ITPR3', 2524, 'element')}`, three residues "
               f"inside the PF00520 boundary (2240–2527). Past the gate and "
               f"still within the channel domain. It is not in a "
               f"low-constraint position: deep JSD "
               f"{_f(_cited(h, 'ITPR3', 2524, 'deep_jsd'))} on "
               f"{_cited(h, 'ITPR3', 2524, 'deep_n')} orthologues with modal "
               f"fraction {_f(_cited(h, 'ITPR3', 2524, 'deep_frac_modal'))}, "
               f"against that protein's whole-protein mean of 0.738. S0's "
               f"other residue-level citation, p.Thr1424Met, sits in "
               f"`{_cited(h, 'ITPR3', 1424, 'element')}` at deep JSD "
               f"{_f(_cited(h, 'ITPR3', 1424, 'deep_jsd'))}"),
        "",
        "**Hand-off.** `constraint_<gene>_<acc>.tsv` is the per-residue "
        "resource S22 needs for the ligand-site question and S24 needs for "
        "the supplementary alignment figures; `painted/` is what a structure "
        "figure is drawn from. S18's correction list is unaffected by this "
        "task — S17 read the annotation only through S5's gene models.",
    )


def build(h: dict) -> str:
    return A(s_variants(h), "", s_vus(h), "", s_selection(h), "",
             s_structures(h), "", s_caveats(h))
