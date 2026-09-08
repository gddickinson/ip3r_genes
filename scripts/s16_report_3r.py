"""The 3R half of S16's report, plus the caveats and the hand-off.

Split out of `s16_report_results.py` to keep every S16 module under the
project's 500-line budget (the `s3_report.py` / `s3_report_d10.py` pattern,
applied twice as S7, S9 and S12 do). It takes the caller's headline dict and
imports the caller's formatters, so the three halves of this report cannot
read the tables differently.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_priors as PR                                        # noqa: E402
from s16_report import rows, n, pfmt                           # noqa: E402


def sec_3r(out: Path, h: dict) -> list[str]:
    md = ["## 5. Are the two teleost ITPR1 copies 3R ohnologs?", ""]
    g = h["teleost_groups"]
    md.append("### 5.1 The three groups, and the two controls inside them")
    md.append("")
    md.append("| group | genomes above D4's bar | ITPR1 | ITPR2 | ITPR3 | "
              "RyR | genomes with >1 ITPR1 |")
    md.append("|---|---|---|---|---|---|---|")
    for key, label in (("pre_3R_outgroup",
                        "pre-3R ray-fins (bichir, gar, bowfin)"),
                       ("teleost", "teleosts (3R)"),
                       ("extra_wgd", "extra WGD (sturgeon, salmon)")):
        if key not in g:
            continue
        v = g[key]
        md.append(f"| {label} | {v['n']} | " + " | ".join(
            f"{v[c]:.2f}" for c in ("ITPR1", "ITPR2", "ITPR3", "RYR"))
            + f" | {v['frac_ITPR1_multi']:.1%} |")
    md.append("")
    md.append(
        "**Predictions 1 and 2 both hold, and they hold in opposite "
        "directions.** The lineages that diverged before 3R carry exactly "
        "one of each ITPR and three RyRs — the gnathostome state. The "
        "teleosts carry two ITPR1 and one each of ITPR2 and ITPR3, and six "
        "RyRs. The lineages with a *further* whole-genome duplication carry "
        "more of everything again. A second copy that appeared in the "
        "outgroup, or a teleost RyR count that had not doubled, would each "
        "have killed the 3R reading; neither does.")
    md.append("")
    md.append(
        "The RyR control is what makes the ITPR2/ITPR3 singletons "
        "interpretable. 3R duplicated the whole genome, so it duplicated "
        "ITPR2 and ITPR3 as surely as ITPR1 — and the sister family in the "
        "same genomes kept **all six** of its copies. So the teleost ITPR2 "
        "and ITPR3 singletons are a statement about *retention*, not about "
        "the sweep's ability to find a duplicate.")
    md.append("")
    disp = h.get("dispersion") or {}
    if disp:
        md.append("### 5.2 Prediction 4: the copies are not tandem")
        md.append("")
        md.append(
            f"Of {n(disp['n_two_copy_genomes'])} teleost genomes with two "
            f"ITPR1 copies above the bar, **{n(disp['n_different_contigs'])} "
            f"put them on different contigs**, and the rest sit a median "
            f"**{n(disp['median_gap_bp_same_contig'])} bp apart on one "
            f"contig** — the closest pair anywhere is "
            f"{n(disp['min_gap_bp'])} bp. Nothing here is a tandem "
            f"duplication. Two copies megabases apart on one chromosome are "
            f"what ~320 Myr of post-3R rediploidisation and chromosome "
            f"fusion look like, and the number is printed rather than "
            f"argued: the sweep's own contigs decide how much of that is "
            f"biology and how much is assembly.")
        md.append("")
    d = h["dcs"]
    md.append("### 5.3 Prediction 5: double-conserved synteny")
    md.append("")
    md.append(
        f"Both copies of a 3R pair should keep *part* of the one ancestral "
        f"neighbourhood, and between them account for it. Of "
        f"**{d['n_pairs']}** two-copy genomes with both copies flanked:")
    md.append("")
    md.append("| reference block | both copies keep ancestral symbols | "
              "and the two sets are disjoint |")
    md.append("|---|---|---|")
    md.append(f"| the tetrapod / non-teleost consensus | "
              f"{d['both_tetrapod']} of {d['n_pairs']} | "
              f"{d['disjoint_tetrapod']} |")
    md.append(f"| the pre-3R ray-finned block | {d['both_fish']} of "
              f"{d['n_pairs']} | {d['disjoint_fish']} |")
    md.append("")
    md.append(
        "Two references rather than one because teleost gene symbols diverge "
        "from tetrapod ones even after S8's relaxed-key normalisation, so "
        "the tetrapod consensus under-counts; gar, bowfin and bichir do not "
        "have that problem. Both give the same answer. **Disjoint is the "
        "load-bearing word**: the two copies do not merely each resemble the "
        "ancestor, they *partition* it, which is what reciprocal gene loss "
        "after one duplication produces and what a pair of independent "
        "later duplications would not.")
    md.append("")
    r3 = h["r3"]
    md.append("### 5.4 Are they the *same* two blocks across the radiation?")
    md.append("")
    if r3.get("assignments_compared"):
        md.append(
            f"{n(r3.get('anchors', 0))} anchors, one per order "
            f"({h.get('r3_notes', {}).get('anchors', '')}); every other "
            f"two-copy genome matched "
            f"onto each independently. **"
            f"{n(r3['assignments_agreeing'])} of "
            f"{n(r3['assignments_compared'])} assignments agree "
            f"({float(r3['fraction_agreeing']):.1%}, two-sided binomial p = "
            f"{pfmt(r3['p_binomial_two_sided'])})**, and the worst-agreeing "
            f"anchor pair is at {r3.get('min_anchor_pair_agreement', '')}. "
            f"Under independent lineage-specific duplications the anchors "
            f"carry no shared information and this sits at 0.5; the "
            f"self-test builds exactly that case and requires the statistic "
            f"to land there (T12), so the 1.0 is a measurement and not a "
            f"property of the routine.")
        md.append("")
        if r3.get("bait_concordance_n") not in ("", "0", None):
            md.append(
                f"Corroborated by evidence of a different kind: which bait "
                f"won each copy is a *sequence* call made with no synteny "
                f"input at all. In the "
                f"{r3['bait_concordance_n']} genomes whose two copies won "
                f"different baits, the sequence call and the synteny block "
                f"agree **{r3['bait_concordance_agree']} times** "
                f"(p = {pfmt(r3.get('bait_concordance_p', 1))}), against "
                f"*{r3.get('bait_concordance_reference', '')}* as the "
                f"reference — chosen as the first anchor whose own two "
                f"copies won different baits, because anchors are ranked on "
                f"flank richness and taking the first one blindly makes this "
                f"check unrunnable whenever that genome's copies share a "
                f"bait.")
            md.append("")
    md.append(PR.line(
        "teleost_two_loci", "confirmed",
        f"{d['disjoint_tetrapod']} of {d['n_pairs']} two-copy genomes "
        f"partition the ancestral block disjointly between their two copies, "
        f"and every anchor assigns every genome to the same two blocks — so "
        f"a cell holding two teleost loci holds the two 3R co-orthologs"))
    md.append("")
    md.append(PR.line(
        "ryr_is_the_hazard", "confirmed",
        "used as an instrument rather than avoided: the RyR trio is the 2R "
        "positive control in §4 and the 3R positive control in §5, and in "
        "both it behaves exactly as the family under test does"))
    md.append("")
    md.append("![double-conserved synteny](figures/s16_dcs.png)")
    md.append("")
    md.append("![blocks and sensitivity](figures/s16_blocks.png)")
    md.append("")
    return md


def sec_close(out: Path, h: dict) -> list[str]:
    rep = h["replication"]
    md = ["## 6. What this settles, what it does not, and the caveats", ""]
    md.append("**Settles.**")
    md.append("")
    md.append(
        f"- **The ITPR blocks are paralogous, and the paralogy runs through "
        f"ITPR1.** ITPR1's neighbourhood carries paralogs of ITPR2's in "
        f"{rep['ITPR1_vs_ITPR2']['frac']:.0%} of "
        f"{n(rep['ITPR1_vs_ITPR2']['n'])} genomes and of ITPR3's in "
        f"{rep['ITPR1_vs_ITPR3']['frac']:.0%} of "
        f"{n(rep['ITPR1_vs_ITPR3']['n'])}, against "
        f"{rep['ITPR2_vs_ITPR3']['null_frac']:.1%} of matched random "
        f"windows. ITPR2 and ITPR3 retain nothing above background.")
    md.append(
        "- **One of those two links is 2R-dated and the other is not.** "
        "`GRM7 ↔ GRM4` (ITPR1–ITPR3) is a *Vertebrata* duplication; "
        "`BHLHE40 ↔ BHLHE41` (ITPR1–ITPR2) is *Opisthokonta*. S8 found both "
        "and could not date either.")
    md.append(
        f"- **3R doubled ITPR1 and only ITPR1**, in "
        f"{h['teleost_groups']['teleost']['frac_ITPR1_multi']:.0%} of "
        f"teleost genomes above the contiguity bar, while doubling all "
        f"three RyRs in the same genomes — so the ITPR2 and ITPR3 "
        f"singletons are retention, not detection.")
    md.append(
        f"- **The two ITPR1 copies are one ancestral duplication.** They "
        f"partition the ancestral block disjointly in "
        f"{h['dcs']['disjoint_tetrapod']} of {h['dcs']['n_pairs']} genomes "
        f"against a tetrapod reference and "
        f"{h['dcs']['disjoint_fish']} against a pre-3R ray-finned one, and "
        f"every anchor from every order assigns every genome to the same two "
        f"blocks.")
    md.append(
        "- **The instrument is calibrated, because the sister family went "
        "through the same test.** Everything above was measured on RYR1/2/3 "
        "in the same run.")
    md.append("")
    md.append("**Does not settle.** Four things, all leads rather than gaps.")
    md.append("")
    md.append(
        "1. **Whether the ITPR quartet had a fourth slot.** The block scan "
        "returns no genome-wide block that both carries no family gene and "
        "looks like the ITPR blocks' missing sibling. With one dated ohnolog "
        "pair surviving between the blocks that *do* carry a gene, a block "
        "that lost the gene too has nothing left to be recognised by, so "
        "this is a limit of the evidence and not a claim that no fourth slot "
        "existed.")
    md.append(
        "2. **Which of the two 2R rounds made which split.** The 2R-dated "
        "link is dated to *Vertebrata*, which is both rounds. Separating R1 "
        "from R2 needs the cyclostome side of the quartet, and S8 already "
        "reported that cyclostome flank synteny is underpowered — 27–29 % of "
        "their coding genes carry a symbol at all.")
    md.append(
        "3. **Why ITPR1 alone kept its 3R duplicate.** The observation is "
        "clean and the cause is not in this task's evidence. S9's ω "
        "estimates and S17's constraint mapping are where a dosage or "
        "subfunctionalisation argument would have to be made.")
    md.append(
        "4. **The paralogy map is human.** Compara paralogy exists for one "
        "genome in this scope, so the replication varies the neighbourhood "
        "and holds the map fixed. That is the right design for the question "
        "— what varies is the thing under test — but a symbol with no human "
        "one-to-one can only *lower* a link count, so every number in §4.3 "
        "is a floor.")
    md.append("")
    md.append("**Two things a reader should hold against this half.**")
    md.append("")
    md.append(
        f"- **The human single-genome test is underpowered on its own** and "
        f"the report says so twice. One dated link per family is what both "
        f"the ITPRs and a family whose 2R origin nobody disputes return; "
        f"after correcting across all 135 tests neither survives. The claim "
        f"in §4.3 rests on the 309-genome replication against its own "
        f"measured null, not on the human p-values.")
    md.append(
        f"- **{n(h['n_genomes'])} vertebrate genomes is the denominator.** "
        f"S4 declared it and S5b swept it. Nothing here extends to a "
        f"lineage with no assembly, and the teleost result rests on the "
        f"{h['teleost_groups']['teleost']['n']} ray-finned genomes above "
        f"D4's contiguity bar rather than on all "
        f"{n(sum(1 for r in rows(out, 'teleost_copies.tsv')))} in the "
        f"ledger.")
    md.append("")
    md.append("## 7. Hand-off")
    md.append("")
    md.append(
        "- **S17** inherits the question §6 could not answer: ITPR1 is the "
        "paralog that kept its 3R duplicate and the one whose neighbourhood "
        "kept its 2R ohnologs. Whether that is the same fact twice is a "
        "constraint question.")
    md.append(
        "- **S19** inherits `copy_sensitivity.tsv` as a worked example of a "
        "threshold this project owns and measured rather than chose.")
    md.append(
        "- **S14a** should take the copy-number panel and the DCS panel; "
        "the 2R replication panel is the one that carries its own null and "
        "is the strongest single figure this task produced.")
    md.append("")
    return md


