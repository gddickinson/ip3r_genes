"""S16 — renders `results/duplication/report.md` purely from the committed
tables (D13). Scope, the instrument and the negative controls live here; the
results live in `s16_report_results.py` (the `s3_report.py` /
`s3_report_d10.py` split).

Nothing in this report recomputes anything. Every number comes back through
`s16_tables.headline()`, which reads the TSVs, so the report and the data
cannot drift.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_paralogon as P                                      # noqa: E402
import s16_tables as TB                                        # noqa: E402
import s16_teleost as T                                        # noqa: E402


def rows(out: Path, name: str) -> list[dict]:
    p = out / name
    return L.read_tsv(p) if p.exists() else []


def n(x, digits: int = 0) -> str:
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{v:,.{digits}f}" if digits else f"{int(round(v)):,}"


def pfmt(x) -> str:
    """A p-value at four decimal places is 0.0000, which hides how strong a
    claim is rather than how weak it is (S15a's rule)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{v:.3g}" if v >= 1e-4 else f"{v:.2e}"


def scope(out: Path, h: dict) -> list[str]:
    md = ["## 1. What is in the analysis", ""]
    md.append(
        "S16 asks two questions about where ITPR1/2/3 came from: are they a "
        "**2R** quartet, and are the two teleost ITPR1 copies **3R** "
        "ohnologs. Both are questions about *blocks* of genome rather than "
        "about the gene, so both are asked of neighbourhoods.")
    md.append("")
    md.append(
        f"The unit of observation is the **locus**, not the ledger cell. The "
        f"S5 ledger holds one row per genome x cell and therefore only each "
        f"cell's best locus, which is exactly the information a duplication "
        f"question needs — S8 and S15a both hit this. So every alignment is "
        f"read back out of the sweep's own `summary.json`: "
        f"**{n(h['n_loci'])} loci across {n(h['n_genomes'])} genomes**, of "
        f"which **{n(h['n_copies'])}** clear the copy rule below.")
    md.append("")
    md.append("| cell | loci | copies | genomes carrying >1 copy, above D4's bar |")
    md.append("|---|---|---|---|")
    loci = rows(out, "loci.tsv")
    for cell in list(L.PARALOGS) + [L.CONTROL_CELL]:
        nl = sum(1 for r in loci if r["cell"] == cell)
        md.append(f"| {cell} | {n(nl)} | {n(h[f'copies_{cell}'])} | "
                  f"{n(h[f'multi_{cell}_above_bar'])} of "
                  f"{n(h['n_above_bar'])} |")
    md.append("")
    md.append(
        "**The RyR cell is this task's positive control and it is not a "
        "fourth paralog.** RYR1/2/3 are a three-member vertebrate family of "
        "the same age and the same 2R candidacy as ITPR1/2/3, and D14 puts "
        "them inside every search this project runs anyway. Every "
        "measurement below is made on both families with the same "
        "instrument, in the same genomes, in the same run. A rule that finds "
        "one family's duplicates and not the other's is a rule about the "
        "instrument, and there is no other way to see that.")
    md.append("")
    return md


def instrument(out: Path, h: dict) -> list[str]:
    md = ["## 2. The instrument", ""]

    md.append("### 2.1 What counts as a copy, and what counts as one gene")
    md.append("")
    md.append(
        f"Two rules, pulling opposite ways. A **copy** is a locus whose own "
        f"model covers at least {L.COV_FULL:.0%} of its bait, over "
        f"{n(L.MIN_COPY_ALIGNED_AA)} aligned residues, at or above the "
        f"sweep's own recording floor of {L.MIN_COPY_IDENTITY:.0%} identity "
        f"— read back from the sweep rather than retyped. And two models are "
        f"**one gene** when they sit on one strand of one contig within "
        f"{n(L.MERGE_GAP_BP)} bp *and* their bait spans are complementary "
        f"rather than repeated: two alignments each covering residues 1–2748 "
        f"of the same bait are two genes however close they are, and two "
        f"covering 1–1300 and 1310–2748 are one gene however far apart.")
    md.append("")
    md.append(
        f"The merge is the brief's requirement and it is load-bearing in one "
        f"direction only — a split model manufactures a duplication, which "
        f"is the claim this task is testing. **It fired on "
        f"{n(h['n_merges'])} of {n(h['n_loci'])} loci**: S5's own clustering "
        f"already chained miniprot's alignments into loci, so the artefact "
        f"this rule exists to catch is not present in this sweep. That zero "
        f"is only worth printing because the rule that returns it is one "
        f"that folds a constructed split — `s16_test_dup.py` T1 builds two "
        f"halves that each clear the copy bar on their own and requires the "
        f"merge to return one copy, which is the only configuration in which "
        f"the merge changes a count.")
    md.append("")
    sens = rows(out, "copy_sensitivity.tsv")
    md.append(
        "The coverage bar is a threshold this task owns, so every count is "
        "recomputed at seven of them in `copy_sensitivity.tsv` "
        f"({sens[0]['cov_bar'] if sens else '?'}–"
        f"{sens[-1]['cov_bar'] if sens else '?'}). Figure "
        "`s16_blocks` panel b draws the result: across the whole range the "
        "number of genomes carrying a second ITPR1 moves by one.")
    md.append("")

    md.append("### 2.2 The 2R test: paralogy, not synteny")
    md.append("")
    md.append(
        "S8 measured microsynteny between the three neighbourhoods and got a "
        "flat **zero** — cross-paralog symbol Jaccard 0.000 in every "
        "direction. That is what a 2R signal looks like to a test that is "
        "looking for the same word twice: after ~500 Myr the neighbours are "
        "no longer the same genes, they are *paralogs* of each other. So the "
        "test needs a paralogy map.")
    md.append("")
    notes = out / "paralogy_map_notes.md"
    release = ""
    if notes.exists():
        for line in notes.read_text().splitlines():
            if line.startswith("- release served by that host"):
                release = line.split("**")[1] if "**" in line else ""
    pmap = rows(out, "paralogy_map.tsv")
    md.append(
        f"The map is Ensembl Compara, pulled from **BioMart against a pinned "
        f"dated archive host** (`{P.M.ARCHIVE_HOST.split('/')[2]}`, serving "
        f"**{release or 'an unrecorded release'}**). Pinning is D24's "
        f"discipline applied to a database release: the rolling host follows "
        f"the release cycle, so a rerun a month later would score the same "
        f"windows against a different Compara tree. The archive's own "
        f"registry is committed as `biomart_registry.xml`, so the release "
        f"the numbers were made on is evidence in this directory and not a "
        f"sentence in a report. **{n(len(pmap))} undirected paralog pairs** "
        f"touch a real neighbourhood gene.")
    md.append("")
    md.append(
        "Three rules make the test a test rather than a restatement:")
    md.append("")
    md.append(
        f"1. **Every ITPR *and* RyR gene is removed from every window** "
        f"(`{', '.join(P.FAMILY_PREFIXES)}`). Counting the ITPR1–ITPR2 "
        f"paralogy itself as evidence that their blocks are paralogous is "
        f"circular — that pair is the thing being explained — and leaving a "
        f"RyR in an ITPR window would import the control's answer into the "
        f"test.")
    md.append(
        f"2. **The null is drawn from real genomic windows** (D17): "
        f"{n(P.N_PERMUTATIONS)} draws per direction from the actual gene "
        f"order at the matched gene count, so the clustering of gene "
        f"families that a shuffled gene set throws away is kept. That is the "
        f"conservative choice — a tandem array inflates the real window and "
        f"the null alike. A Poisson cross-check is committed beside it.")
    md.append(
        f"3. **Links are dated.** Compara puts a duplication node on every "
        f"pair, and that column is the difference between a paralogon test "
        f"and a paralogy test: an older duplication whose two copies happen "
        f"to sit in these blocks is not a 2R ohnolog pair. Two nested "
        f"vocabularies are reported so the answer does not rest on one line "
        f"— `2R_core` = {', '.join(P.LEVELS_2R_CORE)}; `2R_window` adds "
        f"{', '.join(P.LEVELS_2R_WINDOW[len(P.LEVELS_2R_CORE):])}, where "
        f"Compara's LCA reconstruction puts a real 2R pair when the deep "
        f"outgroup genes are missing from its tree. **This changed the "
        f"answer** (§4.1).")
    md.append("")
    md.append(
        f"The primary setting is **±{TB.PRIMARY_WINDOW} genes, "
        f"`{TB.PRIMARY_LEVEL_SET}`** — ±10 is S8's own window and the "
        f"highest-power one available, because the null mean rises with "
        f"window size, so a wider window buys genes and loses signal. ±20 "
        f"and ±30 are committed beside it as the sensitivity.")
    md.append("")
    md.append(
        "Multiple testing is corrected **twice**, because the two families "
        "of tests answer different questions and reporting only one of them "
        "would be a choice made after seeing the numbers (S9's rule, stated "
        "at both scopes). `q_stratum` corrects across the 15 pairs of one "
        "(level set, window) — the family a reader actually reads a row "
        "against, since the three level sets are nested and the three "
        "windows are nested. `q_global` corrects across every test the stage "
        "ran, which is the conservative bound. And `pooled_test.tsv` asks "
        "the question a reader actually has — *do these three "
        "neighbourhoods retain more dated ohnologs than random windows do* "
        "— as **one** hypothesis per family instead of three, which is where "
        "the power is.")
    md.append("")

    md.append("### 2.3 The 3R test: five predictions, all checked")
    md.append("")
    md.append(
        "| # | prediction | how it is checked |")
    md.append("|---|---|---|")
    md.append("| 1 | the pre-3R ray-fins carry one copy | bichir, gar and "
              "bowfin in the sweep, above D4's bar |")
    md.append("| 2 | a lineage with a *further* WGD carries more | "
              "Acipenseriformes and Salmoniformes, reported separately and "
              "never folded into the 3R count |")
    md.append("| 3 | the second copy is spread across the radiation | orders "
              "and classes carrying it |")
    md.append("| 4 | the copies are not tandem | measured separation, not "
              "asserted |")
    md.append("| 5 | double-conserved synteny | both copies keep part of the "
              "*same* ancestral block, against a tetrapod **and** a pre-3R "
              "ray-finned reference |")
    md.append("")
    md.append(
        f"Flank sets come from S8's committed `flanks.tsv`, which flanked "
        f"**every locus** and not one per cell, so the brief's *flank sets "
        f"for every copy* is a join on (accession, cell, contig, start) "
        f"rather than a second extraction (D13) — and S16 therefore cannot "
        f"disagree with S8 about what a flank is. A copy enters the DCS test "
        f"with at least {T.MIN_INFORMATIVE} informative symbols; the "
        f"ancestral block is the set carried by at least {T.REF_FRAC:.0%} of "
        f"the single-copy non-teleost loci.")
    md.append("")
    md.append(
        "The sixth check is the one that separates *one* ancestral "
        "duplication from a series of independent lineage-specific ones: "
        "**cross-anchor block identity**. The obvious statistic — does the "
        "matched pairing of two genomes beat the crossed one — has no power, "
        "because each genome's own copy labels are arbitrary and the winner "
        "is a coin flip either way. So several well-flanked genomes from "
        "different orders are taken as anchors, every other genome's two "
        "copies are matched onto each anchor independently, and the anchors "
        "are asked whether they agree. Agreement between two anchors is "
        "defined only up to one global flip, so the statistic is "
        "max(f, 1−f) against a two-sided binomial at p = 0.5.")
    md.append("")
    return md


def controls(out: Path) -> list[str]:
    md = ["### 2.4 The negative controls, run before anything is written", ""]
    md.append(
        "`scripts/s16_test_dup.py` (the pattern of "
        "`s5_bait_screen.self_test()` and `s15_test_loss.py`) constructs "
        "cases for every rule above and aborts the run on a failure. Two of "
        "S16's failure modes are silent, which is why the tests are mostly "
        "tests on *refusal* and on *reachability*:")
    md.append("")
    md.append("| test | what it refuses, or requires to be reachable |")
    md.append("|---|---|")
    for row in [
        ("T1", "a split model whose halves *each* clear the copy bar must "
               "merge to one copy — the only configuration in which the "
               "merge changes a count"),
        ("T2", "two models each covering the whole bait stay two copies "
               "however close, or the rule would delete real duplications"),
        ("T3", "complementary spans on opposite strands, or on two contigs, "
               "are not one split gene"),
        ("T4", "the merge is order-invariant, or a copy number is a fact "
               "about a JSON file"),
        ("T5", "each copy floor fires on its own violation and nothing else"),
        ("T6", "no ITPR or RyR symbol survives in a window — the "
               "circularity guard"),
        ("T7", "the duplication-node filter discriminates, **and the null is "
               "built from the same link set as the test**"),
        ("T8", "the permutation null is non-degenerate, and a window with no "
               "paralogs scores 0"),
        ("T9", "nine adjacent zinc fingers are one family, not nine"),
        ("T10", "BH is monotone, bounded, and never below its own p"),
        ("T11", "every 3R group is assigned by rule, including the "
                "species-level cyprinid 4R this scope contains none of"),
        ("T12", "**random blocks must agree at about chance** — a "
                "cross-anchor statistic that returned 1.0 on noise would "
                "make the 3R result unfalsifiable"),
        ("T13", "and two genuinely shared blocks must be recovered even when "
                "each genome's copy order is scrambled"),
        ("T14", "the binomial tail is exact"),
    ]:
        md.append(f"| {row[0]} | {row[1]} |")
    md.append("")
    return md


def render(out_dir: Path | None = None) -> Path:
    import s16_report_results as R
    out = out_dir or L.out_dir()
    h = TB.headline(out)
    md = ["# S16 — duplication history of the IP3 receptor family", ""]
    md.append(f"*Generated {time.strftime('%Y-%m-%d %H:%M')}; "
              f"{n(h['n_loci'])} loci in {n(h['n_genomes'])} genomes; "
              f"{n(P.N_PERMUTATIONS)} permutations per null, seed "
              f"{P.SEED}.*")
    md.append("")
    md.append(
        "**What this task asks.** S7 gave the family's tree and S13 placed "
        "its two duplications on the species tree. Neither says what kind of "
        "event made them. This task asks the genome instead: a whole-genome "
        "duplication copies a *block*, so the three ITPR neighbourhoods "
        "should still be paralogous to one another, and the teleost "
        "duplicate should still sit in the block its single ancestor "
        "occupied.")
    md.append("")
    md += scope(out, h)
    md += instrument(out, h)
    md += controls(out)
    md += R.render(out, h)
    path = out / "report.md"
    path.write_text("\n".join(md) + "\n")
    return path


def main() -> None:
    print(render())


if __name__ == "__main__":
    main()
