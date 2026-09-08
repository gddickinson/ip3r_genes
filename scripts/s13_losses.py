"""S13 step 3 — the loss map, and what it is actually evidence of.

The reconciliation reports 53 implied losses. Taken at face value that is a
spectacular claim, and it is wrong: `msa_v2` is a **representative**
alignment — one sequence per clade per species by S6's eight rules — so a
species with no ITPR2 tip is almost always a species whose ITPR2 was never
sampled, not one that lost it. A reconciliation cannot tell those apart; only
a genome can.

So every species x paralog cell is checked against the S5 genome ledger,
which *is* an exhaustive per-genome presence/absence call with per-locus
evidence:

  in the gene tree                          -> sampled_in_gene_tree
  absent from it + ledger `found_*`         -> sampling artefact (the gene is
                                               there; S6 did not pick it)
  absent + ledger `absent` + the genome
    carries unassigned family loci          -> paralog_unassignable
  absent + ledger `absent` otherwise        -> corroborated loss
  absent + ledger `tblastn_trace`           -> loss with a remnant
  absent + `fragment`/`assembly_gap`/
           `tblastn_trace_ambiguous`        -> undecidable
  absent + no genome in the manifest        -> out of scope

**`paralog_unassignable` is the rule this audit could not do without.** Run
without it, the only four "corroborated losses" in the whole table were ITPR2
and ITPR3 in *Myxine glutinosa* and *Petromyzon marinus* — and both genomes
carry **three** ITPR loci apiece, all filed by the sweep in the ITPR1 cell
because S5 has no cyclostome-labelled bait to offer the other two (its six
unfilled slots include all three paralogs in cyclostomes). The ledger says
`absent` there about a *cell*, not about a gene, and S7 says the same thing
from the other side: every cyclostome tip sits outside all three paralog
clades. A loss claim built on that cell would report a bait-panel limit as
biology. The rule is a positive test on two independently committed numbers —
the genome's total ITPR locus count against its occupied cells (S5), and the
species' tree-unplaced tips (S7) — and both are written into the row.

Two further things this project measured make the audit sharper than a status
lookup.

**S10 and S12 settled what a database absence means here.** A gene can be
annotated and still reach no protein database — two *Nibea albiflora* loci are
filed as pseudogenes — and all 7 of S10's annotation failures are transcribed,
spliced genes (S12). So `found_unannotated` and `found_no_annotation` are
*presence*, not partial evidence, and they are counted as such.

**A `fragment` is not an absence.** S5's ledger separates `fragment` from
`absent`, and D25's margin calibration is why: a partial alignment can be a
real gene in a broken assembly. It is undecidable here rather than a loss.

    python scripts/s13_losses.py
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s6_lib import binomial                                        # noqa: E402
from s13_lib import (PARALOGS, PHYLO_DIR, RECON_DIR, group_resolver,  # noqa: E402
                     load_species_tree, load_tree, read_tsv, tip_metadata,
                     vertebrate_labels, write_tsv)

LEDGER = ROOT / "results" / "genome_ledger" / "genome_ledger.tsv"

# S5 ledger statuses grouped by what they mean for a loss claim.
PRESENT = {"found_annotated", "found_unannotated", "found_no_annotation"}
ABSENT = {"absent"}
# A decisive tblastn trace is not an undecidable cell — it is a loss that left
# a fossil, which is the strongest kind of loss evidence there is. The
# *ambiguous* traces stay undecidable, as does a fragment (D25).
REMNANT = {"tblastn_trace"}
UNDECIDABLE = {"assembly_gap", "tblastn_trace_ambiguous", "fragment"}

ORDER = ["sampled_in_gene_tree", "sampling_artefact", "paralog_unassignable",
         "corroborated_loss", "loss_with_remnant", "undecidable",
         "no_genome_in_manifest"]


def ledger_by_species() -> tuple:
    """(status per species x paralog, spare-locus count per species).

    The second is the number of ITPR loci the sweep found in that genome
    beyond the number of paralog cells it could fill — the count that says an
    `absent` cell may be a locus filed under another name.
    """
    status: dict = defaultdict(dict)
    loci: dict = Counter()
    filled: dict = Counter()
    for r in read_tsv(LEDGER):
        if r["class"] not in PARALOGS:
            continue                              # the RyR control is not ours
        key = binomial(r["organism"])
        status[key][r["class"]] = r["status"]
        loci[key] += int(r["n_loci"] or 0)
        if r["status"] in PRESENT:
            filled[key] += 1
    spare = {k: loci[k] - filled[k] for k in loci}
    return status, spare, loci


def verdict(in_tree: bool, status, spare: int, unplaced: int) -> str:
    if in_tree:
        return "sampled_in_gene_tree"
    if status is None:
        return "no_genome_in_manifest"
    if status in PRESENT:
        return "sampling_artefact"          # present in the genome, unsampled
    if status in ABSENT:
        if spare > 0 or unplaced > 0:
            return "paralog_unassignable"
        return "corroborated_loss"
    if status in REMNANT:
        return "loss_with_remnant"
    if status in UNDECIDABLE:
        return "undecidable"
    return f"unclassified:{status}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=RECON_DIR)
    args = ap.parse_args()

    meta = tip_metadata()
    sp = load_species_tree()
    led, spare, n_loci = ledger_by_species()
    ml = load_tree(PHYLO_DIR / "rooted.nwk")
    group_of, _ = group_resolver(ml)
    species = sorted(n.name for n in sp.by_name.values() if n.is_leaf)

    sampled: dict = defaultdict(list)
    unplaced: dict = defaultdict(list)
    for lab in vertebrate_labels(meta):
        g = group_of(lab)
        key = binomial(meta[lab]["species"]).replace(" ", "_")
        if g in PARALOGS:
            sampled[(key, g)].append(lab)
        else:
            unplaced[key].append(lab)

    rows, tally = [], {p: Counter() for p in PARALOGS}
    for sp_name in species:
        plain = sp_name.replace("_", " ")
        for p in PARALOGS:
            tips = sampled.get((sp_name, p), [])
            status = led.get(plain, {}).get(p)
            n_unpl = len(unplaced.get(sp_name, []))
            v = verdict(bool(tips), status, spare.get(plain, 0), n_unpl)
            tally[p][v] += 1
            rows.append([plain, p, len(tips), ";".join(sorted(tips)),
                         status or "", v, n_unpl,
                         n_loci.get(plain, ""), spare.get(plain, "")])
    header = ["species", "paralog", "n_gene_tree_tips", "gene_tree_tips",
              "s5_ledger_status", "verdict", "n_tree_unplaced_tips",
              "n_ledger_itpr_loci", "n_spare_loci"]
    write_tsv(args.out / "paralog_presence.tsv", header, rows)

    summary = []
    for p in PARALOGS:
        for v in ORDER + sorted(set(tally[p]) - set(ORDER)):
            if tally[p][v]:
                summary.append([p, v, tally[p][v]])
    total = Counter()
    for p in PARALOGS:
        total.update(tally[p])
    for v in ORDER + sorted(set(total) - set(ORDER)):
        if total[v]:
            summary.append(["all", v, total[v]])
    write_tsv(args.out / "loss_verdicts.tsv",
              ["paralog", "verdict", "n_species"], summary)

    # Every corroborated or remnant loss, named, with its evidence — the rows
    # a reader is entitled to check one by one.
    detail = [r for r in rows
              if r[5] in ("corroborated_loss", "loss_with_remnant",
                          "paralog_unassignable")]
    write_tsv(args.out / "corroborated_losses.tsv", header, detail)

    for p in PARALOGS:
        line = "  ".join(f"{v}={tally[p][v]}" for v in ORDER if tally[p][v])
        print(f"{p:8s} {line}")
    print(f"\nall      " + "  ".join(f"{v}={total[v]}" for v in ORDER if total[v]))
    print(f"\n{len(detail)} loss-claim cells written to corroborated_losses.tsv")
    print(f"wrote paralog_presence.tsv + loss_verdicts.tsv + "
          f"corroborated_losses.tsv to {args.out}")


if __name__ == "__main__":
    main()
