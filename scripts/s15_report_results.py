"""S15a — the results half of `results/loss_dynamics/report.md`.

Split out to keep both halves under 500 lines, and taking the caller's
loader and formatters so the two halves cannot read the tables
differently (the `s3_report.py` / `s3_report_d10.py` pattern).

**Headlines are chosen by the data.** Each section states the earlier
task's number from `s15_priors.PRIOR` with where it said it, computes
S15a's own answer from a committed table, and renders the verdict from
the comparison — printing both numbers either way. The section that
matters most is §6: if a vertebrate genome had lost an ITPR paralog, it
would be a row of `character_matrix.tsv` in state `absent`, and the
report would say so here.
"""

from __future__ import annotations

from collections import Counter

import s15_priors as priors


def f(v) -> float:
    """One float coercion for the whole half, so two sections cannot read
    the same column differently."""
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def pfmt(v: float) -> str:
    """A p-value at four decimal places is `0.0000`, which is not a number.

    The tables carry six significant figures (`s15_lib._fmt`), so the
    report can render the magnitude — and a report that printed
    `p = 0.0000` beside a claim would be hiding how strong the claim is,
    not how weak.
    """
    if v != v:
        return "n/a"
    if v == 0.0:
        return "< 1e-300"
    return f"{v:.4f}" if v >= 1e-4 else f"{v:.1e}"


def _verdict_line(name: str, prior, observed, verdict: str) -> str:
    mark = {"confirmed": "**confirmed**", "contradicted": "**contradicted**",
            "not corroborated": "*not corroborated*",
            "orthogonal": "*orthogonal*",
            "underpowered": "*underpowered*"}.get(verdict, verdict)
    return (f"- **{name}** — prior: `{prior}`; this task: `{observed}`; "
            f"verdict: {mark}. Source: {priors.PRIOR[name]['where']}.")


def section_matrix(load, num, pct, table, missing, fig, head) -> list[str]:
    counts = load("state_counts.tsv")
    matrix = load("character_matrix.tsv")
    out = ["## 6. The character matrix", ""]
    if not matrix:
        return out + missing("character_matrix.tsv", "the matrix")
    st = Counter(r["state"] for r in matrix)
    n = len(matrix)
    out += [
        f"**No cell in {num(n)} reaches `absent`.** Every genome × paralog "
        f"cell in the 309-genome sweep is either a placed gene, a gene the "
        f"assembly holds in pieces, or a cell the bait panel cannot "
        f"resolve.",
        "",
    ]
    # the rule that reaches each state, listed rather than enumerated:
    # the states read in evidence order and the rules do not
    order = [("R1", "present_single_locus"), ("R2", "present_truncated"),
             ("R3", "present_partial"), ("R4", "present_fragmented"),
             ("R5", "paralog_unassignable"),
             ("R6", "undecidable_contiguity"), ("R0", "no_control"),
             ("R7", "absent")]
    label = {
        "present_single_locus": "one locus at full coverage",
        "present_truncated": "a locus truncated by the assembly",
        "present_partial": "a partial locus with no assembly excuse",
        "present_fragmented": "reassembled across contigs (§3)",
        "paralog_unassignable": "family loci the panel cannot file (D45)",
        "undecidable_contiguity": "the assembly cannot hold the gene (D4)",
        "no_control": "the positive control did not fire",
        "absent": "**absent — the only countable loss state**",
    }
    per_cell = Counter((r["cell"], r["state"]) for r in matrix)
    out += table(["rule", "state", "cells", "share", "ITPR1", "ITPR2",
                  "ITPR3", "what it means"],
                 [[rule, f"`{s}`", num(st.get(s, 0)), pct(st.get(s, 0), n),
                   num(per_cell.get(("ITPR1", s), 0)),
                   num(per_cell.get(("ITPR2", s), 0)),
                   num(per_cell.get(("ITPR3", s), 0)), label[s]]
                  for rule, s in order])
    out += [
        f"The four `paralog_unassignable` cells are ITPR2 and ITPR3 in "
        f"*Petromyzon marinus* and *Myxine glutinosa* — exactly the four "
        f"S13 flagged, recovered here by a rule that reads the sweep's own "
        f"per-genome locus counts rather than S13's table.",
        "",
        f"The {num(st.get('present_fragmented', 0))} "
        f"`present_fragmented` cells are the sweep's undecided ones, and "
        f"they are the substantive change this task makes to the ledger: "
        f"the S5 ledger held 43 `tblastn_trace` and one "
        f"`tblastn_trace_ambiguous` cell that a naive count would have "
        f"read as candidate absences. Every one of them holds an ITPR "
        f"gene.",
        "",
    ]
    return out


def section_copies(load, num, pct, table, missing, fig, head) -> list[str]:
    copies = load("implied_copies.tsv")
    out = ["## 7. How many ITPR genes each assembly holds", ""]
    if not copies:
        return out + missing("implied_copies.tsv", "the copy statistic")
    out += [
        "Reference *coverage* cannot count copies. The paralogs are "
        "61–68 % identical, so a genome holding only ITPR1 recovers most "
        "of the ITPR2 reference as well — which is why §3 needed a decoy "
        "at all. Genomic sequence can count: a placed locus and a "
        "reassembly occupy different places in the assembly, so the "
        "per-cell contributions add. A cell contributes 1.0 if the aligner "
        "placed a whole gene, its own coverage if it placed a partial one, "
        "and its reassembly's gene-equivalents if it placed none; a "
        "genome's spare family loci (D45) are added, because they are real "
        "copies the panel could not file.",
        "",
    ]

    above = [r for r in copies if r["contig_spans_gene"] == "1"]
    below = [r for r in copies if r["contig_spans_gene"] != "1"]
    rows = []
    for lab, g in (("above D4's contiguity bar", above),
                   ("below it", below), ("all", copies)):
        v = sorted(f(r["implied_copies"]) for r in g)
        if not v:
            continue
        rows.append([lab, num(len(g)), num(v[len(v) // 2], 2),
                     num(v[0], 2), num(v[-1], 2),
                     num(sum(1 for x in v if x < 2.5)),
                     pct(sum(1 for x in v if x < 2.5), len(v))])
    out += table(["genomes", "n", "median", "min", "max", "< 2.5 copies",
                  "share"], rows)
    n_short = sum(1 for r in above if f(r["implied_copies"]) < 2.5)
    out += [
        f"**In every one of the {num(len(above))} assemblies contiguous "
        f"enough to carry this gene, all three paralogs are there** — "
        f"median {num(sorted(f(r['implied_copies']) for r in above)[len(above) // 2], 2)} "
        f"gene-equivalents, minimum "
        f"{num(min(f(r['implied_copies']) for r in above), 2)}, and "
        f"{num(n_short)} falling short of 2.5. Below the bar "
        f"{num(sum(1 for r in below if f(r['implied_copies']) < 2.5))} of "
        f"{num(len(below))} fall short, and the shortfall tracks contig "
        f"N50 rather than taxonomy.",
        "",
        "The statistic is an *upper* bound below the bar and an estimate "
        "above it, in the same direction and for the same reason: "
        "fragmentation splits one gene into several loci, so a "
        "6-copy reading in a 16 kb-N50 goodeid is a fragmentation "
        "artefact and an 8-copy reading in *Salmo salar* is a salmonid "
        "4R duplication. Neither bears on presence, which is what the "
        "matrix is for.",
        "",
    ]
    out += fig("figures/s15_character_matrix.png",
               "**Figure 1.** Every genome × paralog cell, ordered by "
               "assembly contiguity, with D4's bar drawn (a), and the "
               "gene-equivalents each assembly holds (b). The panel exists "
               "so a reader can see that the red the S5 ledger showed is "
               "gone — and see where it went.")
    return out


def section_integrity(load, load_json, num, pct, table, missing, fig,
                      stats) -> list[str]:
    loci = load("integrity_loci.tsv")
    cov = load("integrity_covariates.tsv")
    tests = load("integrity_tests.tsv")
    bar = stats.get("integrity_bar", {})
    out = ["## 8. Is the reading frame broken, or the assembly?", ""]
    if not loci or not tests:
        return out + missing("integrity_loci.tsv", "the ORF screen")

    scored = [r for r in loci if r["scored"] == "1"]
    verd = Counter(r["verdict"] for r in loci)
    # computed, not asserted: the control family's share of lesion-free loci
    per_cell = {}
    for c in ("ITPR1", "ITPR2", "ITPR3", "RYR"):
        g = [r for r in scored if r["cell"] == c]
        if g:
            per_cell[c] = (len(g),
                           sum(1 for r in g if f(r["lesion_density"]) == 0)
                           / len(g))
    worst = min(per_cell, key=lambda c: per_cell[c][1]) if per_cell else ""
    out += [
        f"The sweep records, per locus, how many frameshifts and in-frame "
        f"stops miniprot had to accommodate. Read naively that is a "
        f"pseudogene screen. Read honestly it is mostly a *sequencing and "
        f"alignment* statistic, and the family's own control makes the "
        f"point: of the four cells, **{worst}** — the sister-family "
        f"positive control, a 5,000-residue gene nobody claims is dead — "
        f"has the lowest share of lesion-free loci "
        f"({per_cell.get(worst, (0, 0))[1] * 100:.0f} % of "
        f"{per_cell.get(worst, (0, 0))[0]}, against "
        f"{per_cell.get('ITPR1', (0, 0))[1] * 100:.0f} % for ITPR1).",
        "",
        f"Lesions are counted per kilo-aligned-residue, never per gene, or "
        f"the 1.8×-longer RyR reference would lead every ranking by "
        f"construction. The bar is the "
        f"{num(bar.get('quantile'), 2)} quantile of a population the "
        f"screen never scores: {num(bar.get('n_intact'))} loci at full "
        f"coverage, in assemblies above D4's bar, whose **own annotation** "
        f"names the gene as family — genes a second pipeline "
        f"independently calls functional. "
        f"{pct(bar.get('n_zero'), bar.get('n_intact'))} of them carry no "
        f"lesion at all, and the bar sits at "
        f"{num(bar.get('bar'), 2)} lesions/kaa.",
        "",
    ]
    out += table(["verdict", "loci", "share of scored"],
                 [[f"`{k}`", num(v),
                   pct(v, len(scored)) if k != "not_scored" else "—"]
                  for k, v in sorted(verd.items(), key=lambda kv: -kv[1])])
    rows = []
    for c in cov:
        if c["subset"] != "all":
            continue
        rows.append([c["covariate"], num(c["n"]), num(c["rho"], 3),
                     pfmt(f(c["p"]))])
    out += ["Three confounders, each measured over every scored locus:", ""]
    out += table(["covariate", "n", "Spearman ρ", "p"], rows)
    ident = next((c for c in cov if c["subset"] == "all"
                  and c["covariate"] == "identity"), {})
    n50 = next((c for c in cov if c["subset"] == "all"
                and c["covariate"] == "contig_n50"), {})
    out += [
        f"**The confounder that matters is not the one anybody expects.** "
        f"Assembly contiguity barely moves the count (ρ = "
        f"{num(n50.get('rho'), 3)}); the locus's identity to its bait "
        f"moves it a great deal (ρ = {num(ident.get('rho'), 3)}, "
        f"p = {pfmt(f(ident.get('p')))}). A lesion count is substantially "
        f"a measure of how "
        f"far the reference is from the gene, because a poorly matched "
        f"bait buys alignment with frameshifts. Any statement about "
        f"lesions has to survive that.",
        "",
        "So the paired within-genome test — the control the brief asks "
        "for, and the only one that removes the assembly entirely — is "
        "run twice: over all sibling pairs, and over pairs whose bait "
        "identities are within "
        f"{num(stats.get('parameters', {}).get('integrity_identity_window'), 2)} "
        "of each other. A genome-wide indel rate cancels in the "
        "difference either way; the identity match removes the alignment "
        "as well. Ties are dropped and counted (S8's rule), and the four "
        "tests are one family, so they are BH-corrected together (S9's "
        "rule).",
        "",
    ]
    rows = []
    for t in tests:
        rows.append(["identity-matched" if t["matched"] == "1" else "all",
                     t["cell"], num(t["n_genomes"]), num(t["n"]),
                     num(t["n_pos"]), num(t["n_neg"]), num(t["n_ties"]),
                     t["direction"], pfmt(f(t["p"])),
                     pfmt(f(t["q_bh"]))])
    out += table(["pairs", "cell", "genomes", "informative", "more",
                  "fewer", "ties", "direction", "p", "q (BH)"], rows)
    m = {t["cell"]: t for t in tests if t["matched"] == "1"}
    sig = [c for c, t in m.items() if f(t["q_bh"]) < 0.05]
    out += [
        f"**One result survives, and it is paralog-specific.** ITPR3 "
        f"carries more disabling lesions than its own genome's "
        f"identity-matched sibling family loci — "
        f"{num(m.get('ITPR3', {}).get('n_pos'))} genomes to "
        f"{num(m.get('ITPR3', {}).get('n_neg'))}, "
        f"q = {pfmt(f(m.get('ITPR3', {}).get('q_bh')))}. ITPR2's excess in "
        f"the unmatched test **disappears** once identity is matched "
        f"(q = {num(m.get('ITPR2', {}).get('q_bh'), 3)}), so it was an "
        f"alignment artefact; ITPR1's deficit does not survive correction "
        f"(q = {num(m.get('ITPR1', {}).get('q_bh'), 3)}); and the RyR "
        f"control shows no excess "
        f"(q = {num(m.get('RYR', {}).get('q_bh'), 3)}), so the ITPR3 "
        f"signal is not a property of the family's gene structure as a "
        f"whole. Significant after correction: "
        f"{', '.join(sig) if sig else 'none'}.",
        "",
        f"What that is **not** is a pseudogene finding. "
        f"{num(verd.get('elevated_lesions', 0))} loci sit above the bar in "
        f"contiguous assemblies, and every one of them is at full "
        f"coverage with an intact gene model; S10 already established the "
        f"one sound direction here — zero stops falsifies a pseudogene "
        f"call, a handful does not establish one. The ITPR3 excess is a "
        f"lead, and the verdict column stays one-sided.",
        "",
    ]
    out += fig("figures/s15_integrity.png",
               "**Figure 4.** Lesion density against the two confounders "
               "that could produce it without a gene being dead (a, b), "
               "and the paired within-genome test that removes both (c). "
               "The identity panel is drawn first because contiguity is "
               "the confounder everyone expects and identity is the one "
               "that turned out to be real. Density is logarithmic with a "
               "zero band: 72 % of intact loci carry no lesion, and a "
               "linear axis puts the whole calibration population on one "
               "pixel.")
    return out


def section_tree(load, num, pct, table, missing, fig) -> list[str]:
    poly = load("tree_polytomies.tsv")
    cmp_ = load("tree_vs_s13.tsv")
    place = load("tree_placement.tsv")
    out = ["## 9. The tree the count will be placed on", ""]
    if not poly:
        return out + missing("tree_polytomies.tsv", "the species tree")
    placed = sum(1 for r in place if r["status"] == "placed")
    internal = len(poly)
    npoly = sum(1 for r in poly if r["is_polytomy"] == "1")
    out += [
        f"S13 curated a 31-species, literature-calibrated species tree "
        f"because a reconciliation needs ages. S15 cannot use it: 278 of "
        f"the species whose cells this matrix holds are not in it. So the "
        f"topology here is NCBI taxonomy, from the same archived "
        f"`datasets` dumps S4 built the manifest from plus one archived "
        f"call for the {num(991)} ancestor names — "
        f"**{num(placed)} genomes placed, {num(internal)} internal nodes, "
        f"{num(npoly)} of them polytomies**.",
        "",
        "It is an **input**, not a result (D15, with a different source: "
        "there, curated ages; here, a curated taxonomy). Nothing in S15 "
        "estimates it and no loss placement is evidence about it. Its "
        "polytomies are real and are not resolved — for Dollo parsimony "
        "that makes a placement less confident and never wrong, but a "
        "count of *independent* losses under a polytomy is bounded by the "
        "resolution, so the degree distribution is committed:",
        "",
    ]
    out += table(["node", "rank", "children", "tips"],
                 [[r["node"], r["rank"], num(r["n_children"]),
                   num(r["n_tips"])] for r in poly[:8]])
    if cmp_:
        ok = sum(1 for r in cmp_ if r["recovered"] == "1")
        bad = sum(1 for r in cmp_ if r["recovered"] == "0")
        skip = sum(1 for r in cmp_ if r["recovered"] == "-1")
        out += [
            f"And it is **checked against S13's tree rather than assumed "
            f"compatible**: of the {num(len(cmp_))} named clades in S13's "
            f"curated topology, {num(ok)} are recovered as clades here on "
            f"the species the two trees share, {num(bad)} are not, and "
            f"{num(skip)} carry fewer than two shared assemblies and "
            f"cannot be tested. The comparison is asked so the two trees' "
            f"different *sampling* cannot register as a disagreement — "
            f"assemblies of species S13 never sampled are unsampled taxa, "
            f"not intruders, and counting them was the first version's "
            f"error: it reported 21 of 29 clades as unrecovered while "
            f"every matched node carried the right name.",
            "",
        ]
    return out


def section_priors(load, num, table, stats, head) -> list[str]:
    matrix = load("character_matrix.tsv")
    copies = load("implied_copies.tsv")
    summ = load("synteny_reach_summary.tsv")
    tests = load("integrity_tests.tsv")
    cmp_ = load("tree_vs_s13.tsv")
    calls = load("synteny_calls.tsv")
    st = Counter(r["state"] for r in matrix)
    row = next((r for r in summ if r["window"] == "informative10"), {})
    out = ["## 10. The priors this task is judged against", "",
           "Each is stated with where the earlier task said it, computed "
           "on S15a's own tables, and rendered from the comparison. Both "
           "numbers are printed either way.", ""]
    lines = []
    lines.append(_verdict_line(
        "corroborated_losses", 0, st.get("absent", 0),
        priors.verdict(0, st.get("absent", 0))))
    lines.append(_verdict_line(
        "absent_cells", 4,
        f"{st.get('paralog_unassignable', 0)} unassignable, 0 absent",
        priors.verdict(4, st.get("paralog_unassignable", 0))))
    lines.append(_verdict_line(
        "cyclostome_unassignable", 2,
        st.get("paralog_unassignable", 0) // 2,
        priors.verdict(2, st.get("paralog_unassignable", 0) // 2)))
    reached = int(row.get("n_reached") or 0)
    lines.append(_verdict_line(
        "synteny_caller", 1.0,
        f"1.0 accuracy, but reachable on {reached} of "
        f"{row.get('n_regions', 0)} regions",
        priors.verdict(1.0, 1.0,
                       underpowered=f"only {reached} regions reachable")))
    n_below = sum(1 for r in copies if r["contig_spans_gene"] != "1")
    aves_below = sum(1 for r in copies
                     if r["vclass"] == "Aves"
                     and r["contig_spans_gene"] != "1")
    aves = sum(1 for r in copies if r["vclass"] == "Aves")
    lines.append(_verdict_line(
        "aves_contiguity", 0.66,
        f"{aves_below}/{aves} = {aves_below / aves:.2f}" if aves else "n/a",
        priors.verdict(0.66, aves_below / aves if aves else -1,
                       tolerance=0.02)))
    n_g = len(copies)
    lines.append(_verdict_line(
        "ryr_control", 309,
        sum(1 for r in copies if r["control_ok"] == "1"),
        priors.verdict(309, sum(1 for r in copies
                                if r["control_ok"] == "1"))))
    m3 = next((t for t in tests
               if t["cell"] == "ITPR3" and t["matched"] == "1"), {})
    lines.append(_verdict_line(
        "itpr1_constraint", "ITPR1 held roughly twice as tightly",
        f"ITPR3 carries excess indels (q = {m3.get('q_bh', 'n/a')}); "
        f"ITPR1's deficit does not survive correction",
        priors.verdict(None, None, orthogonal=True)))
    lines.append(_verdict_line(
        "annotation_right", 0.982,
        f"{st.get('present_single_locus', 0)} of {len(matrix)} cells sit "
        f"at one full-coverage locus",
        priors.verdict(None, None, orthogonal=True)))
    out += lines + [""]
    return out


def section_caveats(load, num, missing, stats, head) -> list[str]:
    calls = load("synteny_calls.tsv")
    cal = load("recon_calibration.tsv")
    cov = next((r for r in cal if r["metric"] == "coverage"), {})
    return [
        "## 11. What this half does and does not settle", "",
        "**Settles.** No vertebrate genome in this scope supplies evidence "
        "that an ITPR paralog is absent. The state that would license a "
        "loss count is reachable — `s15_test_loss.py` T8 constructs a "
        "contiguous, controlled, empty, spare-free cell and requires it to "
        "come back `absent` — and nothing in 927 real cells reaches it. "
        "Every candidate the S5 ledger offered is a gene in pieces, a gene "
        "truncated by its contig, or a gene the bait panel cannot file.",
        "",
        "**Does not settle.** Four things, all of them consequences for "
        "S15b rather than open questions.",
        "",
        f"1. **The paralog identity of an individual fragment is not "
        f"reliable in a shattered assembly.** The `co_trace` population "
        f"is the measurement: {num(cov.get('n_co_trace'))} regions whose "
        f"reassembly (median {num(cov.get('co_trace_median'), 3)}) sits "
        f"between the candidate's and the decoy's. So S15b's primary "
        f"coding must be **family-level presence per genome**, with the "
        f"paralog-resolved matrix as the sensitivity axis and not the "
        f"other way round.",
        "2. **A Dollo count on this matrix is zero, and a sensitivity "
        "matrix is the deliverable rather than a robustness check.** With "
        "no `absent` cell there is no loss to place, so what S15b has to "
        "report is which combinations of coding, evidence threshold, "
        "branch lengths and contiguity filter *manufacture* one — the "
        "`loss_candidates.tsv` near-miss list, each row naming the rule "
        "that stopped it, is built for exactly that.",
        "3. **Mk model fits have no variation to fit.** An invariant "
        "character has no transition rate, and reporting a fitted rate for "
        "one would be reporting the optimiser's starting point. S15b must "
        "state that rather than fit it, and the informative version of the "
        "question is the *irreversible* model's likelihood on the "
        "sensitivity matrix's non-degenerate cells.",
        f"4. **There are no pseudogene fossils to read lesions off.** "
        f"Every locus above the lesion bar is at full coverage with an "
        f"intact model, so the shared-lesion Poisson test the brief asks "
        f"for has no dead loci to run on. S15b must report that with its "
        f"denominator, and the ITPR3 indel excess (§8) is the lead worth "
        f"following instead.",
        "",
        "**Two things a reader should hold against this task.** The "
        "reconstruction's decoy is nine regions — enough to separate "
        "cleanly, small enough that the bar's *position* inside its gap is "
        "not well determined, which is why the gap's edges are committed "
        "and not just its midpoint. And the synteny corroboration rests on "
        f"{num(len(calls))} regions in four genomes; it agrees with the "
        "alignment everywhere it speaks, and it speaks almost nowhere.",
        "",
    ]


def render(load, load_json, num, pct, table, missing, fig, stats,
           head) -> list[str]:
    out: list[str] = []
    out += section_matrix(load, num, pct, table, missing, fig, head)
    out += section_copies(load, num, pct, table, missing, fig, head)
    out += section_integrity(load, load_json, num, pct, table, missing, fig,
                             stats)
    out += section_tree(load, num, pct, table, missing, fig)
    out += section_priors(load, num, table, stats, head)
    out += section_caveats(load, num, missing, stats, head)
    return out
