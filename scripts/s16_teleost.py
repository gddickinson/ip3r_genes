"""S16 stage `teleost` — are the two teleost ITPR1 copies the 3R ohnologs?

The copy stage found the asymmetry. Above D4's contiguity bar, 93.6 % of
ray-finned genomes carry two ITPR1 copies against 6.4 % for ITPR2 and 5.1 %
for ITPR3 — while the **RyR control**, a three-member family in the same
genomes under the same duplication, is doubled in 100 % of them.

The 3R hypothesis makes five predictions and every one is checked here off
data that already exists:

1. **Outgroup.** The pre-3R ray-fins — bichir (Cladistia), gar and bowfin
   (Holostei) — should carry one copy. If they carry two, the duplication is
   older than 3R and 3R is the wrong label for it.
2. **Extra-WGD control.** Lineages with a *further* duplication —
   Acipenseriformes and Salmoniformes — should carry more, and are reported
   separately rather than folded into the 3R count. They are the positive
   control for "more WGD, more copies".
3. **Breadth.** The second copy should be spread across the teleost
   radiation, not concentrated in one order.
4. **Dispersion.** WGD copies are not tandem. Measured, not assumed: this
   family answers it in a way the PIEZO project's did not.
5. **Double-conserved synteny.** Each teleost copy should keep part of the
   *same* ancestral neighbourhood, so both share flanks with the single
   pre-3R block while being mutually distinguishable — the standard 3R
   signature and the one test a series of independent lineage-specific
   duplications cannot fake. Run against a tetrapod reference **and** a
   pre-3R ray-finned one, because teleost symbols diverge from tetrapod ones
   even after S8's relaxed-key normalisation and the tetrapod consensus
   therefore under-counts.

Flank sets come from S8's committed `flanks.tsv`, which flanked **every
locus** and not one per cell — so the brief's "flank sets for every copy"
is a join, not a re-extraction (D13).

Outputs (`results/duplication/`): `teleost_copies.tsv`, `dcs_flanks.tsv`,
`dcs_test.tsv`, `r3_summary.tsv`, `block_assignments.tsv`.
"""

from __future__ import annotations

import collections
import math
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s8_flank_lib as F                                       # noqa: E402

#: Ray-finned lineages that diverged **before** the teleost-specific genome
#: duplication. Cladistia (bichir, reedfish) and Holostei (gar, bowfin) are
#: the internal control for prediction 1.
PRE_3R_CLASSES = {"Cladistia"}
PRE_3R_ORDERS = {"Semionotiformes", "Lepisosteiformes", "Amiiformes",
                 "Polypteriformes"}
#: Lineages carrying a WGD *other than* 3R: Acipenseriformes (sturgeon,
#: paddlefish) sit outside Teleostei with their own duplication, and
#: salmonids carry Ss4R on top of 3R. Extra copies here are expected.
EXTRA_WGD_ORDERS = {"Salmoniformes", "Acipenseriformes"}
#: The cyprinid 4R is species-level, not order-level: *Danio* is a plain
#: teleost and these are not. The list is kept even when this scope contains
#: none of them, so a later scope that does cannot silently mis-group them.
EXTRA_WGD_SPECIES = {"Cyprinus carpio", "Carassius auratus",
                     "Carassius gibelio", "Sinocyclocheilus grahami",
                     "Sinocyclocheilus anshuiensis",
                     "Sinocyclocheilus rhinocerous"}

#: The cell whose duplication this stage is about, and the control cell.
FOCAL = "ITPR1"
MIN_INFORMATIVE = 5
#: Consensus bar for the ancestral block: a symbol carried by this fraction
#: of the single-copy non-teleost loci.
REF_FRAC = 0.40


def group_of(row: dict) -> str:
    if (row.get("vclass") in PRE_3R_CLASSES
            or row.get("vorder") in PRE_3R_ORDERS):
        return "pre_3R_outgroup"
    if (row.get("vorder") in EXTRA_WGD_ORDERS
            or row.get("organism") in EXTRA_WGD_SPECIES):
        return "extra_wgd"
    if row.get("vclass") == "Actinopteri":
        return "teleost"
    return "non_actinopterygian"


# ------------------------------------------------------------ copy census

CENSUS_HEADER = (["accession", "organism", "vclass", "vorder", "group",
                  "contig_spans_gene"]
                 + [f"{p}_copies" for p in L.PARALOGS]
                 + ["RYR_copies", f"{FOCAL}_contigs", f"{FOCAL}_same_contig"])


def census(cn: dict[tuple[str, str], dict], log=print
           ) -> tuple[list[list], list[dict]]:
    by_acc: dict[str, dict] = collections.defaultdict(dict)
    for (acc, cell), c in cn.items():
        by_acc[acc][cell] = c
    rows, recs = [], []
    for acc, cells in by_acc.items():
        m = next(iter(cells.values()))
        if m["vclass"] not in ("Actinopteri", "Cladistia"):
            continue
        grp = group_of(m)
        focal = cells.get(FOCAL, {})
        contigs = sorted(set(focal.get("copy_contigs") or []))
        rec = {"accession": acc, "organism": m["organism"],
               "vclass": m["vclass"], "vorder": m["vorder"], "group": grp,
               "spans": m["contig_spans_gene"],
               "contigs": contigs}
        for p in L.PARALOGS + (L.CONTROL_CELL,):
            rec[p] = cells.get(p, {}).get("n_copies", 0)
        recs.append(rec)
        rows.append([acc, rec["organism"], rec["vclass"], rec["vorder"], grp,
                     int(rec["spans"])]
                    + [rec[p] for p in L.PARALOGS] + [rec[L.CONTROL_CELL],
                       ";".join(contigs),
                       int(rec[FOCAL] > 1 and len(contigs) < rec[FOCAL])])
    rows.sort(key=lambda r: (r[4], r[1]))
    for grp in ("pre_3R_outgroup", "teleost", "extra_wgd"):
        sub = [r for r in recs if r["group"] == grp and r["spans"]]
        if not sub:
            log(f"[teleost] {grp}: no genome above the contiguity bar")
            continue
        log(f"[teleost] {grp} (above D4's bar, n={len(sub)}): "
            + "  ".join(f"{p} mean {statistics.mean(r[p] for r in sub):.2f}"
                        for p in L.PARALOGS + (L.CONTROL_CELL,)))
    return [CENSUS_HEADER] + rows, recs


# ------------------------------------------------------------- dispersion

DISP_HEADER = ["group", "cell", "n_two_copy_genomes", "n_different_contigs",
               "n_same_contig", "frac_different", "median_gap_bp_same_contig",
               "min_gap_bp", "note"]


def dispersion(recs_loci: list[dict], census_recs: list[dict],
               log=print) -> list[list]:
    """Are the two copies on two chromosomes, or side by side?

    A WGD copies a whole chromosome, so the two copies start on different
    ones; a tandem duplication puts them adjacent. Between those extremes
    sits the case this family actually shows, and reporting the median
    separation is the only way to say which one it is rather than assert it.
    """
    grp_of = {r["accession"]: r["group"] for r in census_recs}
    spans_of = {r["accession"]: r["spans"] for r in census_recs}
    by: dict[tuple, list[dict]] = collections.defaultdict(list)
    for r in recs_loci:
        if r["accession"] in grp_of and L.is_copy(r):
            by[(r["accession"], r["cell"])].append(r)
    stats: dict[tuple, list] = collections.defaultdict(list)
    for (acc, cell), rs in by.items():
        if len(rs) != 2 or not spans_of.get(acc):
            continue
        a, b = sorted(rs, key=lambda r: (r["contig"], r["start"]))
        same = a["contig"] == b["contig"]
        gap = (b["start"] - a["end"]) if same else None
        stats[(grp_of[acc], cell)].append((same, gap))
    rows = []
    for (grp, cell), vals in sorted(stats.items()):
        same = [g for s, g in vals if s]
        n_diff = sum(1 for s, _ in vals if not s)
        note = ("tandem" if same and statistics.median(same) < 1_000_000
                else "dispersed on one chromosome" if same else "")
        rows.append([grp, cell, len(vals), n_diff, len(same),
                     round(n_diff / len(vals), 4) if vals else 0.0,
                     int(statistics.median(same)) if same else "",
                     min(same) if same else "", note])
    for r in rows:
        if r[1] != FOCAL:
            continue
        tail = (f"; the rest sit a median {r[6]:,} bp apart on one contig"
                if r[6] != "" else "")
        log(f"[teleost] dispersion {r[0]} {r[1]}: {r[3]}/{r[2]} genomes put "
            f"the two copies on different contigs{tail}")
    return rows


# ------------------------------------------------- double-conserved synteny

FLANK_HEADER = ["accession", "organism", "vclass", "vorder", "group",
                "cell", "copy_id", "label", "contig", "start", "end",
                "n_informative", "keys"]


def flank_sets(recs_loci: list[dict], log=print) -> list[dict]:
    """Every copy's flank set, joined out of S8's committed tables.

    S8 flanked every locus (its labels carry `CELL.idx`), so the brief's
    "flank sets for every copy" is a join on (accession, cell, contig,
    start) rather than a second extraction — which also means S16 cannot
    disagree with S8 about what a flank is.
    """
    s8_loci = {}
    for r in L.read_tsv(L.PROJECT_ROOT / "results" / "synteny" / "loci.tsv"):
        s8_loci[(r["accession"], r["cell"], r["contig"],
                 int(r["start"]))] = r["label"]
    keys: dict[str, set[str]] = collections.defaultdict(set)
    for r in L.read_tsv(L.PROJECT_ROOT / "results" / "synteny" / "flanks.tsv"):
        if r.get("window") == "fixed10" and r.get("informative") == "1":
            k = (r.get("relaxed_key") or "").upper()
            if k:
                keys[r["label"]].add(k)
    out = []
    for r in recs_loci:
        if not L.is_copy(r):
            continue
        label = s8_loci.get((r["accession"], r["cell"], r["contig"],
                             r["start"]))
        if not label:
            continue
        ks = keys.get(label, set())
        if len(ks) < MIN_INFORMATIVE:
            continue
        out.append({**r, "label": label, "group": group_of(r),
                    "keys": frozenset(ks)})
    log(f"[teleost] joined flanks for {len(out)} copies with "
        f">={MIN_INFORMATIVE} informative symbols")
    return out


DCS_HEADER = ["accession", "organism", "vorder", "cell", "copy_a", "copy_b",
              "n_keys_a", "n_keys_b", "j_a_vs_b", "shared_tetrapod_a",
              "shared_tetrapod_b", "both_share_tetrapod",
              "tetrapod_partition_disjoint", "shared_fish_a",
              "shared_fish_b", "both_share_fish", "fish_partition_disjoint",
              "tetrapod_symbols_a", "tetrapod_symbols_b"]


def dcs(sets: list[dict], cell: str = FOCAL, log=print
        ) -> tuple[list[list], dict]:
    focal = [s for s in sets if s["cell"] == cell]
    n_copies: dict[str, int] = collections.Counter(s["accession"]
                                                   for s in focal)
    singles = [s for s in focal
               if n_copies[s["accession"]] == 1
               and s["group"] == "non_actinopterygian"]
    counts: collections.Counter = collections.Counter()
    for s in singles:
        counts.update(s["keys"])
    ref = {k for k, c in counts.items() if c >= REF_FRAC * max(1, len(singles))}
    fish = set().union(*[s["keys"] for s in focal
                         if s["group"] == "pre_3R_outgroup"]) \
        if any(s["group"] == "pre_3R_outgroup" for s in focal) else set()
    log(f"[teleost] ancestral {cell} block: {len(ref)} symbols from "
        f"{len(singles)} single-copy non-teleost loci; pre-3R ray-finned "
        f"block {len(fish)} symbols")

    by_acc: dict[str, list[dict]] = collections.defaultdict(list)
    for s in focal:
        if s["group"] in ("teleost", "extra_wgd"):
            by_acc[s["accession"]].append(s)
    rows = []
    for acc, v in sorted(by_acc.items()):
        if len(v) != 2:
            continue
        a, b = sorted(v, key=lambda s: (s["contig"], s["start"]))
        ra, rb = a["keys"] & ref, b["keys"] & ref
        fa, fb = a["keys"] & fish, b["keys"] & fish
        rows.append([acc, a["organism"], a["vorder"], cell, a["label"],
                     b["label"], len(a["keys"]), len(b["keys"]),
                     round(F.jaccard(a["keys"], b["keys"]), 4),
                     len(ra), len(rb), int(bool(ra and rb)),
                     int(bool(ra and rb and not (ra & rb))),
                     len(fa), len(fb), int(bool(fa and fb)),
                     int(bool(fa and fb and not (fa & fb))),
                     ";".join(sorted(ra)), ";".join(sorted(rb))])
    stat = {
        "n_pairs": len(rows),
        "both_tetrapod": sum(r[11] for r in rows),
        "disjoint_tetrapod": sum(r[12] for r in rows),
        "both_fish": sum(r[15] for r in rows),
        "disjoint_fish": sum(r[16] for r in rows),
        "ref_size": len(ref), "fish_size": len(fish),
        "n_singles": len(singles),
    }
    log(f"[teleost] DCS vs tetrapod block: {stat['both_tetrapod']}/"
        f"{stat['n_pairs']} genomes have BOTH copies in it "
        f"({stat['disjoint_tetrapod']} partition it disjointly); vs pre-3R "
        f"fish block {stat['both_fish']}/{stat['n_pairs']} "
        f"({stat['disjoint_fish']} disjoint)")
    return [DCS_HEADER] + rows, stat


# ------------------------------------------------- cross-anchor block identity

def binom_p(k: int, n: int, p: float = 0.5) -> float:
    if n == 0:
        return 1.0
    return max(0.0, min(1.0, sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i)
                                 for i in range(k, n + 1))))


def _assign(target: tuple[frozenset, frozenset],
            anchor: tuple[frozenset, frozenset]) -> tuple[int, float]:
    t1, t2 = target
    a1, a2 = anchor
    s0 = F.jaccard(t1, a1) + F.jaccard(t2, a2)
    s1 = F.jaccard(t1, a2) + F.jaccard(t2, a1)
    return (0 if s0 >= s1 else 1), abs(s0 - s1)


R3_HEADER = ["metric", "value", "note"]
ASSIGN_HEADER = ["accession", "organism", "vorder", "anchor",
                 "anchor_organism", "orientation", "margin"]


def block_consistency(sets: list[dict], cell: str = FOCAL, log=print,
                      n_anchors: int = 6) -> tuple[list[list], list[list]]:
    """Are the two teleost blocks the *same* two blocks across species?

    One ancestral duplication predicts two neighbourhoods that stay mutually
    recognisable across the radiation: copy-A here matches copy-A there.
    Independent lineage-specific duplications predict no such correspondence.

    The obvious statistic — does the matched pairing of two genomes beat the
    crossed one — has **no power**, because each genome's own copy labels are
    arbitrary and the winner is a coin flip either way. The test used instead
    is *cross-anchor agreement*: several well-flanked genomes from different
    orders are anchors, every other genome's two copies are matched onto each
    anchor independently, and the anchors are asked whether they agree.
    Agreement between two anchors is defined only up to one global flip, so
    the statistic is max(f, 1-f) against a two-sided binomial at p = 0.5.
    """
    focal = [s for s in sets
             if s["cell"] == cell and s["group"] in ("teleost", "extra_wgd")]
    by_acc: dict[str, list[dict]] = collections.defaultdict(list)
    for s in focal:
        by_acc[s["accession"]].append(s)
    genomes = {a: sorted(v, key=lambda s: (s["contig"], s["start"]))
               for a, v in by_acc.items() if len(v) == 2}
    if len(genomes) < 4:
        return [R3_HEADER], [ASSIGN_HEADER]

    # One anchor per order, falling back to the genus where the manifest has
    # no order: two order-less margin species would otherwise each count as
    # their own order and could take two of six anchor slots between them.
    best_by_order: dict[str, str] = {}
    for acc, v in genomes.items():
        order = v[0]["vorder"] or v[0]["organism"].split()[0]
        score = min(len(v[0]["keys"]), len(v[1]["keys"]))
        cur = best_by_order.get(order)
        if cur is None or score > min(len(genomes[cur][0]["keys"]),
                                      len(genomes[cur][1]["keys"])):
            best_by_order[order] = acc
    anchors = sorted(best_by_order.values(),
                     key=lambda a: (-min(len(genomes[a][0]["keys"]),
                                         len(genomes[a][1]["keys"])), a))
    anchors = anchors[:n_anchors]

    detail, calls = [], {}
    for acc, v in sorted(genomes.items()):
        target = (v[0]["keys"], v[1]["keys"])
        for anc in anchors:
            if anc == acc:
                continue
            av = genomes[anc]
            orient, margin = _assign(target, (av[0]["keys"], av[1]["keys"]))
            if margin <= 0:
                continue
            calls[(anc, acc)] = orient
            detail.append([acc, v[0]["organism"], v[0]["vorder"], anc,
                           genomes[anc][0]["organism"], orient,
                           round(margin, 4)])

    agree_rows, tot_agree, tot_n = [], 0, 0
    for i in range(len(anchors)):
        for j in range(i + 1, len(anchors)):
            ai, aj = anchors[i], anchors[j]
            shared = [g for g in genomes
                      if (ai, g) in calls and (aj, g) in calls
                      and g not in (ai, aj)]
            if len(shared) < 5:
                continue
            same = sum(1 for g in shared if calls[(ai, g)] == calls[(aj, g)])
            agree_rows.append((ai, aj, len(shared),
                               round(max(same, len(shared) - same)
                                     / len(shared), 4)))
            tot_agree += max(same, len(shared) - same)
            tot_n += len(shared)
    frac = tot_agree / tot_n if tot_n else 0.0
    p = 2 * binom_p(tot_agree, tot_n) if tot_n else 1.0

    # Corroboration from a different kind of evidence: which bait won each
    # copy is a *sequence* call made with no synteny input at all.
    # The reference for this check is the first anchor whose own two copies
    # won *different* baits: anchors are ranked on flank richness, which has
    # nothing to do with whether the sequence evidence can distinguish the
    # copies, so taking anchors[0] blindly makes the check unrunnable
    # whenever that genome's copies happen to share a bait.
    bait_agree = bait_n = 0
    ref = next((a for a in anchors
                if genomes[a][0]["bait"].split("|")[0]
                != genomes[a][1]["bait"].split("|")[0]), None)
    ref_baits = ([s["bait"].split("|")[0] for s in genomes[ref]]
                 if ref else ["", ""])
    if ref is not None:
        for acc, v in genomes.items():
            if acc == ref or (ref, acc) not in calls:
                continue
            names = [s["bait"].split("|")[0] for s in v]
            if names[0] == names[1]:
                continue
            seq = (0 if names[0] == ref_baits[0]
                   else 1 if names[0] == ref_baits[1] else None)
            if seq is None:
                continue
            bait_n += 1
            bait_agree += int(calls[(ref, acc)] == seq)

    rows = [
        ["cell", cell, "the duplicated paralog cell under test"],
        ["two_copy_genomes_flanked", len(genomes),
         f"genomes with both {cell} copies flanked at >={MIN_INFORMATIVE} "
         f"symbols"],
        ["anchors", len(anchors),
         "; ".join(f"{genomes[a][0]['organism']} "
                   f"({genomes[a][0]['vorder'] or 'order unassigned'})"
                   for a in anchors)],
        ["anchor_pairs_compared", len(agree_rows), ""],
        ["assignments_compared", tot_n, "genome x anchor-pair calls"],
        ["assignments_agreeing", tot_agree, ""],
        ["fraction_agreeing", round(frac, 4),
         "0.5 under independent lineage-specific duplications"],
        ["p_binomial_two_sided", L.fmt(p), ""],
        ["min_anchor_pair_agreement",
         min((r[3] for r in agree_rows), default=0.0), ""],
        ["bait_concordance_reference",
         genomes[ref][0]["organism"] if ref else "",
         "first anchor whose own two copies won different baits" if ref
         else "no anchor's two copies won different baits — check not run"],
        ["bait_concordance_n", bait_n,
         "genomes where the two copies won different baits"],
        ["bait_concordance_agree", bait_agree,
         "synteny block agrees with the sequence (bait) call"],
        ["bait_concordance_p",
         L.fmt(2 * binom_p(max(bait_agree, bait_n - bait_agree), bait_n))
         if bait_n else "", "independent evidence, same partition"],
    ]
    log(f"[teleost] cross-anchor agreement {tot_agree}/{tot_n} = {frac:.1%} "
        f"(p = {p:.2e}); bait concordance {bait_agree}/{bait_n}")
    return [R3_HEADER] + rows, [ASSIGN_HEADER] + detail


def run(out_dir: Path | None = None, log=print) -> dict:
    import s16_copy_number as CN
    out = out_dir or L.out_dir()
    recs, _ = CN.build_loci(log=lambda *a: None)
    cn = CN.copy_number(recs)
    cen, cen_recs = census(cn, log)
    L.write_tsv(out / "teleost_copies.tsv", cen[0], cen[1:])
    L.write_tsv(out / "dispersion.tsv", DISP_HEADER,
                dispersion(recs, cen_recs, log))
    sets = flank_sets(recs, log)
    L.write_tsv(out / "dcs_flanks.tsv", FLANK_HEADER,
                [[s["accession"], s["organism"], s["vclass"], s["vorder"],
                  s["group"], s["cell"], s["label"], s["label"], s["contig"],
                  s["start"], s["end"], len(s["keys"]),
                  ";".join(sorted(s["keys"]))] for s in sets])
    test, stat = dcs(sets, log=log)
    L.write_tsv(out / "dcs_test.tsv", test[0], test[1:])
    r3, assign = block_consistency(sets, log=log)
    L.write_tsv(out / "r3_summary.tsv", r3[0], r3[1:])
    L.write_tsv(out / "block_assignments.tsv", assign[0], assign[1:])
    return {"census": cen, "dcs": test, "dcs_stat": stat, "r3": r3}


def main() -> None:
    run()


if __name__ == "__main__":
    main()
