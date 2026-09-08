"""The results half of S16's report, split to keep both halves under the
project's 500-line budget and taking the caller's headline dict so the two
halves cannot read the tables differently (the `s3_report.py` /
`s3_report_d10.py` pattern).

This half carries the copy-number landscape and the 2R question; the 3R
question, the caveats and the hand-off are in `s16_report_3r.py` (the split
applied twice, as S7, S9 and S12 do).

**Headlines chosen by the data.** Every prior is stated in `s16_priors.PRIOR`
with where the earlier task said it, computed on S16's own tables, and
rendered with both numbers printed either way. The comparison that has to be
sayable is the one that goes badly: if the ITPR neighbourhoods had scored at
their null while the RyR control fired, the table it would be said from is
`paralogon_species_summary.tsv`, and §4.3 is where it would appear.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_priors as PR                                        # noqa: E402
import s16_tables as TB                                        # noqa: E402
import s16_teleost as T                                        # noqa: E402
from s16_report import rows, n, pfmt                           # noqa: E402


def sec_copies(out: Path, h: dict) -> list[str]:
    md = ["## 3. The copy-number landscape", ""]
    clade = rows(out, "copy_number_by_clade.tsv")
    md.append(
        "Above D4's contiguity bar, **every non-teleost gnathostome genome "
        "in the sweep carries exactly one of each ITPR and three RyRs**. "
        "The whole of this task's copy-number variation is in the ray-fins.")
    md.append("")
    md.append("| class | genomes above the bar | ITPR1 | ITPR2 | ITPR3 | "
              "RyR | genomes with >1 ITPR1 |")
    md.append("|---|---|---|---|---|---|---|")
    order = sorted({r["vclass"] for r in clade})
    for vclass in order:
        sub = {r["cell"]: r for r in clade if r["vclass"] == vclass}
        if "ITPR1" not in sub or int(sub["ITPR1"]["n_above_bar"]) < 1:
            continue
        c = sub["ITPR1"]
        md.append(f"| {vclass} | {n(c['n_above_bar'])} | "
                  + " | ".join(n(sub[k]["mean_copies_above_bar"], 2)
                               for k in ("ITPR1", "ITPR2", "ITPR3", "RYR"))
                  + f" | {float(c['frac_multi_above_bar']):.1%} |")
    md.append("")
    md.append(
        f"{n(h['n_lesion_rich_copies'])} of {n(h['n_copies'])} copies are "
        f"`lesion_rich` on S15a's own bar, read back from its committed "
        f"verdicts rather than recomputed. The call is one-sided by "
        f"decision: S10 established that zero lesions falsifies a "
        f"pseudogene call and a handful does not establish one.")
    md.append("")
    md.append(PR.line(
        "no_losses", "confirmed",
        "every genome above the contiguity bar carries all three, and the "
        "only cells at zero copies are below it or cyclostome — S16 counts "
        "copies where S15b counted presence, and they agree"))
    md.append("")
    md.append(PR.line(
        "extra_copies_same_paralog", "confirmed",
        f"{n(h['multi_ITPR1_above_bar'])} of {n(h['n_above_bar'])} genomes "
        f"above the bar carry a second ITPR1 and §5 shows both copies sit in "
        f"the same ancestral block, so a cell holding two loci holds two "
        f"copies of one gene"))
    md.append("")
    md.append("![copy number](figures/s16_copy_number.png)")
    md.append("")
    return md


def sec_2r(out: Path, h: dict) -> list[str]:
    md = ["## 4. Are ITPR1/2/3 a 2R quartet?", ""]

    md.append("### 4.1 What the human windows retain, and how old it is")
    md.append("")
    md.append(
        "Two paralog links survive between the three ITPR neighbourhoods at "
        "every window size, and they are **exactly the two flanking families "
        "S8 found by a completely different instrument** — S8 matched "
        "normalised gene-symbol roots across 309 genomes; this matches "
        "Compara paralogy in one. Two instruments, one answer:")
    md.append("")
    md.append("| pair | link | Compara duplication node | 2R-dated? |")
    md.append("|---|---|---|---|")
    links = rows(out, "paralogy_links.tsv")
    seen = set()
    for r in links:
        if r["pair_class"] != "ITPR_vs_ITPR":
            continue
        key = (r["group_a"], r["group_b"], r["symbol_a"], r["symbol_b"])
        if key in seen:
            continue
        seen.add(key)
        md.append(f"| {r['group_a']}–{r['group_b']} | {r['symbol_a']} ↔ "
                  f"{r['symbol_b']} | {r['duplication_level']} | "
                  f"{'**yes**' if r['in_2R_window'] == '1' else 'no'} |")
    md.append("")
    md.append(
        "**Dating the links changed the answer, and it is the one thing here "
        "S8 could not do.** `BHLHE40 ↔ BHLHE41` is an *Opisthokonta* "
        "duplication — a pair far older than the vertebrates whose two "
        "copies happen to sit beside ITPR1 and ITPR2. It is real paralogy "
        "and it is not a 2R ohnolog pair. `GRM7 ↔ GRM4` is dated to "
        "*Vertebrata*, which is what 2R means. So the single 2R-dated "
        "retained ohnolog pair among these blocks links **ITPR1 and ITPR3**, "
        "and ITPR2 retains none with anybody.")
    md.append("")

    md.append("### 4.2 Against the null, and against the RyR control")
    md.append("")
    prim = h["paralogon_primary"]
    md.append(f"At the primary setting (±{TB.PRIMARY_WINDOW} genes, "
              f"`{TB.PRIMARY_LEVEL_SET}`):")
    md.append("")
    md.append("| pair | dated links | permutation null (mean) | enrichment | "
              "p | q within the 15-pair stratum |")
    md.append("|---|---|---|---|---|---|")
    for key, v in sorted(prim.items(),
                         key=lambda kv: (kv[1]["class"], kv[0])):
        if v["class"].startswith("ITPR_vs_RYR"):
            continue
        tag = " *(control)*" if v["class"].startswith("RYR") else ""
        enr = (v["links"] / v["null_mean"]) if v["null_mean"] else 0.0
        md.append(f"| {key.replace('_vs_', '–')}{tag} | {v['links']} | "
                  f"{v['null_mean']:.4f} | {enr:.0f}× | {pfmt(v['p'])} | "
                  f"{pfmt(v['q_stratum'])} |")
    md.append("")
    pool = h["pooled_primary"]
    if pool:
        md.append(
            "Asked once per family rather than three times per pair — which "
            "is the question, and is where the power is:")
        md.append("")
        md.append("| family | dated links over its 3 pairs | null (mean) | "
                  "p |")
        md.append("|---|---|---|---|")
        for fam in ("ITPR", "RYR"):
            if fam not in pool:
                continue
            v = pool[fam]
            tag = " *(control)*" if fam == "RYR" else ""
            md.append(f"| {fam}{tag} | {v['links']} | {v['null_mean']:.4f} | "
                      f"{pfmt(v['p'])} |")
        md.append("")
    md.append(
        "**The two families behave the same way, and that is the result.** "
        "Each retains exactly one vertebrate-dated ohnolog pair between two "
        "of its three neighbourhoods; each clears its own permutation null "
        "pooled; and neither survives correction across all 135 tests the "
        "stage ran. The RyR trio's 2R origin is not in question, so what "
        "this measures is the *instrument's* ceiling on a single human "
        "genome rather than a difference between the families — and it is "
        "why the replication in §4.3 exists.")
    md.append("")

    md.append("### 4.3 Replicated across the sweep, against its own null")
    md.append("")
    rep = h["replication"]
    nul = h["replication_null"]
    md.append(
        f"The same human paralogy map, asked of S8's flank sets in every "
        f"swept genome, against **matched random neighbourhoods in the same "
        f"genomes** — {n(nul['pairs'])} control window pairs drawn by "
        f"`s8_control.sample_windows` unchanged, of which "
        f"{n(nul['with_link'])} carry a link.")
    md.append("")
    md.append("| pair | genomes compared | with any link | with a 2R-dated "
              "link | vertebrate classes | matched random windows | p |")
    md.append("|---|---|---|---|---|---|---|")
    for pair, v in sorted(rep.items()):
        md.append(f"| {pair.replace('_vs_', '–')} | {n(v['n'])} | "
                  f"{n(v['hit'])} ({v['frac']:.1%}) | {n(v['hit_2R'])} "
                  f"({v['frac_2R']:.1%}) | {v['classes']} | "
                  f"{v['null_frac']:.1%} | {pfmt(v['p'])} |")
    md.append("")
    md.append(
        "**This is the 2R answer.** The ITPR1 neighbourhood is paralogous to "
        "both the ITPR2 and the ITPR3 neighbourhood in most vertebrate "
        "genomes and across most vertebrate classes, ~20–30× above a null "
        "measured in the same genomes with the same map and the same "
        "symbol normalisation. The ITPR2 and ITPR3 neighbourhoods are "
        "paralogous to each other **at exactly the background rate** "
        f"({rep['ITPR2_vs_ITPR3']['frac']:.1%} against "
        f"{rep['ITPR2_vs_ITPR3']['null_frac']:.1%}, p = "
        f"{pfmt(rep['ITPR2_vs_ITPR3']['p'])}). And the dated column splits "
        "the two surviving links cleanly: the ITPR1–ITPR3 link is "
        "vertebrate-dated in "
        f"{n(rep['ITPR1_vs_ITPR3']['hit_2R'])} genomes, the ITPR1–ITPR2 link "
        f"in {n(rep['ITPR1_vs_ITPR2']['hit_2R'])}.")
    md.append("")
    md.append(PR.line(
        "paralogon_through_itpr1", "confirmed",
        f"a second instrument, on different evidence, recovers the same "
        f"shape and the same two families: ITPR1 with ITPR2 "
        f"{rep['ITPR1_vs_ITPR2']['frac']:.1%} of genomes, ITPR1 with ITPR3 "
        f"{rep['ITPR1_vs_ITPR3']['frac']:.1%}, ITPR2 with ITPR3 "
        f"{rep['ITPR2_vs_ITPR3']['frac']:.1%} against a "
        f"{rep['ITPR2_vs_ITPR3']['null_frac']:.1%} background. What S16 adds "
        f"is the date, and it removes one of the two links from the 2R "
        f"account"))
    md.append("")
    md.append(PR.line(
        "sister_pair", "not corroborated",
        "the tree's sister pair is still the one pair whose neighbourhoods "
        "share nothing — ITPR2 and ITPR3 sit at the background rate at both "
        "dating levels. These remain different measurements and neither "
        "overturns the other: a tree estimates the order of duplication, a "
        "retained flanking ohnolog records which copies survived deletion "
        "beside each gene, and 2R quartets lose flank copies independently "
        "of the duplication order. What can now be added is that the same "
        "asymmetry survives being dated"))
    md.append("")
    md.append(PR.line(
        "duplication_placement", "orthogonal",
        "S13 places the ITPR1 split on the vertebrate stem and the "
        "ITPR2/ITPR3 split on the gnathostome stem; S16 measures which "
        "neighbours each block kept. A tree node and a retained neighbour "
        "are not the same quantity. They are consistent — a 2R-dated "
        "ohnolog beside ITPR1 and ITPR3 is compatible with both splits — but "
        "the flank evidence cannot separate them and is not offered as "
        "support"))
    md.append("")
    blocks = rows(out, "paralogon_blocks.tsv")
    quart = rows(out, "quartet_test.tsv")
    if blocks and quart:
        honest = [r for r in quart if r["selection_circular"] == "0"]
        sig = sum(1 for r in honest if float(r["p_permutation"]) < 0.05)
        top = [r for r in blocks if r["rank"] == "1"]
        md.append("### 4.4 The fourth slot: the block scan and the quartet")
        md.append("")
        md.append(
            f"Two rounds of duplication make four copies of an ancestral "
            f"block; three carry a family gene, and a fourth surviving with "
            f"its gene deleted would still be paralogous. The scan ranks "
            f"every genome-wide block by the number of distinct gene "
            f"*families* it shares with each window — families rather than "
            f"raw hits, because a tandem array is one duplication and nine "
            f"adjacent zinc fingers must score 1. **No top block carries a "
            f"family gene** ({sum(1 for r in top if r['carries_family_gene'] == 'yes')}"
            f" of {len(top)}).")
        md.append("")
        md.append(
            f"The quartet test then measures paralogy between every pair of "
            f"the six windows and the six top blocks. **{len(quart) - len(honest)} "
            f"of its {len(quart)} pairs are circular by construction and are "
            f"flagged as such in the table**: a window against *its own* top "
            f"block tests the selection — that block was chosen out of "
            f"~23,000 as the one most paralogous to this window — and not "
            f"the quartet. All six of them are significant, which is what "
            f"the selection guarantees. Of the {len(honest)} pairs that are "
            f"not circular, **{sig} are enriched at p < 0.05**.")
        md.append("")
        md.append("| window | top block | shared families | carries a family "
                  "gene |")
        md.append("|---|---|---|---|")
        for r in top:
            md.append(f"| {r['against']} | {r['contig']}:{n(r['start'])}–"
                      f"{n(r['end'])} | {r['n_families']} | "
                      f"{r['carries_family_gene']} |")
        md.append("")
        md.append(
            "So the scan does not hand back a clean fourth ITPR slot, and "
            "the flag is the reason the table cannot be read as though it "
            "did. That is the expected outcome rather than a negative "
            "result: with one "
            "dated ohnolog pair surviving between the blocks that *do* carry "
            "a gene, a block that lost its gene as well as most of its "
            "neighbours has nothing left to be recognised by. It is reported "
            "because a scan that could only ever confirm is not a scan.")
        md.append("")
    md.append("![the 2R test](figures/s16_paralogon.png)")
    md.append("")
    return md


def render(out: Path, h: dict) -> list[str]:
    """The whole results half, both modules, so the report has one entry
    point and the two halves cannot be assembled in different orders."""
    import s16_report_3r as R3
    md: list[str] = []
    md += sec_copies(out, h)
    md += sec_2r(out, h)
    md += R3.sec_3r(out, h)
    md += R3.sec_close(out, h)
    return md
