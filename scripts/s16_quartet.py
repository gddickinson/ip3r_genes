"""S16 stage `quartet` — the fourth 2R slot, and the same test in 309 genomes.

Split out of `s16_paralogon.py` to keep both under the project's 500-line
budget (the `s3_report.py` / `s3_report_d10.py` pattern). Two tests:

**The quartet.** Two rounds of whole-genome duplication make four copies of
an ancestral block. Three of them carry an ITPR gene. If the fourth survives
with its ITPR deleted, it is still a paralogous block — so the test takes
the three ITPR windows plus the best non-family block the scan found for
each, and measures paralogy between every pair against its own permutation
null (D17). The RyR windows go through the same test as the positive
control, because a quartet statistic with nothing to be compared against is
a number, not a result.

**Cross-species replication.** The human answer is one genome. S8 already
extracted flanking symbols for every ITPR and RyR locus in 309 genomes, so
the same paralogy map is asked of each of them: does *this* species'
ITPR1 neighbourhood carry paralogs of *this* species' ITPR2 neighbourhood?
The map stays human — it is the only Compara paralogy this project has —
so what varies across the replication is the neighbourhood, which is the
thing under test.

Outputs (`results/duplication/`): `quartet_test.tsv`, `paralogon_species.tsv`,
`paralogon_species_summary.tsv`.
"""

from __future__ import annotations

import collections
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_ensembl_map as M                                    # noqa: E402
import s16_paralogon as P                                      # noqa: E402
import s8_control as C                                         # noqa: E402
import s8_flank_lib as F                                       # noqa: E402

#: A genome contributes a replication pair only when both neighbourhoods
#: resolve at least this many symbols to a human gene. Below that a zero is
#: a statement about the annotation, not about the neighbourhood — S8's own
#: informativeness bar, reused rather than re-chosen.
MIN_MAPPED = 5

POOLED_HEADER = ["level_set", "window_n", "family", "n_pairs",
                 "observed_links", "n_families", "null_mean", "enrichment",
                 "p_permutation", "null_max", "n_permutations"]


def pooled_test(by_contig: dict[str, list[tuple]],
                windows: dict[tuple, list[tuple]], sym2ensg: dict[str, str],
                ensg2sym: dict[str, str], pairs: dict[str, list[tuple]],
                groups: list[str], rng: random.Random, log=print
                ) -> list[list]:
    """One question per family instead of three, because that is the question.

    "Do these three neighbourhoods retain more dated ohnologs than random
    windows do" is a single hypothesis about a family, and asking it three
    times as three pairs spends power on a multiple-comparison correction
    for tests nobody wanted separately. The pooled statistic is the total
    link count over a family's three within-family pairs, against a null
    that draws three matched random windows and pools them the same way.
    """
    rows = []
    for n in P.WINDOW_SIZES:
        sliding = P.sliding_windows(by_contig, 2 * n + 1, sym2ensg, step=1)
        fam_pairs = {"ITPR": [], "RYR": []}
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a, b = groups[i], groups[j]
                a_it, b_it = a in P.TEST_CELLS, b in P.TEST_CELLS
                if a_it and b_it:
                    fam_pairs["ITPR"].append((a, b))
                elif not a_it and not b_it:
                    fam_pairs["RYR"].append((a, b))
        for set_name, levels in P.LEVEL_SETS.items():
            for fam, plist in fam_pairs.items():
                obs, fams, draws = 0, set(), []
                for a, b in plist:
                    wa = [x for x in (sym2ensg.get(g[5].upper(), "")
                                      for g in windows[(a, n)]) if x]
                    wb = [x for x in (sym2ensg.get(g[5].upper(), "")
                                      for g in windows[(b, n)]) if x]
                    links = P.links_between(wa, set(wb), pairs, levels)
                    obs += len(links)
                    fams |= {P.family_root(ensg2sym.get(t, t))
                             for _, t, _ in links}
                    draws.append(P.permutation_null(
                        P.paralogs_of(wa, pairs, levels), len(wb), sliding,
                        rng))
                if not draws:
                    continue
                null = [sum(d[k] for d in draws) for k in range(len(draws[0]))]
                mean_null = sum(null) / len(null)
                p = (sum(1 for x in null if x >= obs) + 1) / (len(null) + 1)
                rows.append([set_name, n, fam, len(plist), obs, len(fams),
                             round(mean_null, 4),
                             round(obs / mean_null, 3) if mean_null else "inf",
                             L.fmt(p), max(null), len(null)])
                if set_name == "2R_window":
                    log(f"[quartet] pooled +/-{n:>2} {fam:5s} dated links="
                        f"{obs} null={mean_null:.3f} p={p:.1e}")
    return [POOLED_HEADER] + rows


QUARTET_HEADER = ["block_a", "block_b", "class", "selection_circular",
                  "region_a", "region_b", "n_genes_a", "n_genes_b",
                  "n_families", "n_links", "null_mean", "enrichment",
                  "p_permutation", "carries_family_a", "carries_family_b",
                  "link_symbols"]


def quartet_test(by_contig: dict[str, list[tuple]],
                 windows: dict[tuple, list[tuple]],
                 sym2ensg: dict[str, str], blocks: list[list],
                 loci: dict[str, tuple], groups: list[str],
                 rng: random.Random, log=print, n: int = 20) -> list[list]:
    cand: list[tuple] = []
    for g in groups:
        wg = windows[(g, n)]
        cand.append((f"{g}_window", wg[0][0], wg[0][1], wg[-1][2], "yes",
                     "ITPR" if g in P.TEST_CELLS else "RYR"))
    for g in groups:
        for r in blocks[1:]:
            if r[1] == g and r[9] == "no":
                cand.append((f"top_vs_{g}", r[2], int(r[3]), int(r[4]), "no",
                             "ITPR" if g in P.TEST_CELLS else "RYR"))
                break

    ids_of: dict[str, list[str]] = {}
    for name, contig, s, e, _, _ in cand:
        genes = [g for g in by_contig.get(contig, [])
                 if g[2] >= s and g[1] <= e
                 and not g[5].upper().startswith(P.FAMILY_PREFIXES)]
        ids_of[name] = [x for x in (sym2ensg.get(g[5].upper(), "")
                                    for g in genes) if x]
    universe = {x for genes in by_contig.values()
                for x in (sym2ensg.get(g[5].upper(), "") for g in genes) if x}
    universe |= {g for v in ids_of.values() for g in v}
    adj = M.full_pairs(universe)
    sym_of = {}
    for genes in by_contig.values():
        for g in genes:
            gid = sym2ensg.get(g[5].upper(), "")
            if gid:
                sym_of[gid] = g[5]

    sliding = P.sliding_windows(by_contig, 2 * n + 1, sym2ensg, step=1)
    rows = []
    for i in range(len(cand)):
        for j in range(i + 1, len(cand)):
            na, ca, sa, ea, pa, fa = cand[i]
            nb, cb, sb, eb, pb, fb = cand[j]
            if ca == cb and not (ea < sb or sa > eb):
                continue
            a_ids, b_ids = ids_of[na], set(ids_of[nb])
            links = P.dedupe_links([(g, t, lev) for g in a_ids
                                    for t, lev in adj.get(g, [])
                                    if t in b_ids])
            fams = {P.family_root(sym_of.get(t, t)) for _, t, _ in links}
            a_par = {t for g in a_ids for t, _ in adj.get(g, [])}
            null = P.permutation_null(a_par, len(b_ids), sliding, rng)
            obs = len(links)
            p_perm = (sum(1 for x in null if x >= obs) + 1) / (len(null) + 1)
            mean_null = sum(null) / len(null)
            # A window against *its own* top block is circular: that block
            # was selected, out of ~23,000, as the one most paralogous to
            # this window, so its p-value tests the selection and not the
            # quartet. Flagged in the data rather than left for the prose,
            # because a reader scanning the table would otherwise read six
            # significant rows as six results.
            circular = int(nb == f"top_vs_{na.replace('_window', '')}"
                           or na == f"top_vs_{nb.replace('_window', '')}")
            rows.append([
                na, nb, f"{fa}/{fb}", circular,
                f"{ca}:{sa}-{ea}", f"{cb}:{sb}-{eb}",
                len(a_ids), len(b_ids), len(fams), obs, round(mean_null, 4),
                round(obs / mean_null, 3) if mean_null > 0 else "inf",
                L.fmt(p_perm), pa, pb,
                ";".join(sorted({f"{sym_of.get(g, g)}-{sym_of.get(t, t)}"
                                 for g, t, _ in links}))[:600]])
    honest = [r for r in rows if not r[3]]
    sig = sum(1 for r in honest if float(r[12]) < 0.05)
    log(f"[quartet] {sig}/{len(honest)} non-circular block pairs enriched at "
        f"p < 0.05 ({len(rows) - len(honest)} window-vs-own-top-block pairs "
        f"excluded as circular)")
    return [QUARTET_HEADER] + rows


# --------------------------------------------------------- cross-species

#: Control windows drawn per genome for the replication's own null. Four
#: windows give six control pairs per genome, which is twice the three real
#: pairs a genome can contribute — deliberately, so the null is the better
#: measured of the two.
N_CONTROL_WINDOWS = 4

SPECIES_HEADER = ["accession", "organism", "vclass", "vorder", "pair",
                  "pair_class", "n_mapped_a", "n_mapped_b", "n_links",
                  "n_links_2R", "n_families", "link_symbols"]
SPECIES_SUM_HEADER = ["pair", "pair_class", "n_genomes", "n_with_link",
                      "frac_with_link", "n_with_2R_link", "frac_with_2R_link",
                      "mean_links", "max_links", "n_classes", "classes"]


def cross_species(sym2ensg: dict[str, str], pairs: dict[str, list[tuple]],
                  ensg2sym: dict[str, str], log=print
                  ) -> tuple[list[list], list[list]]:
    """Replicate the window test on S8's flank sets in every swept genome.

    The RyR cell is one cell holding up to three genes, so it is **not** a
    neighbourhood and cannot be one side of a pair: pooling RYR1, RYR2 and
    RYR3's flanks would hand the control three windows' worth of symbols
    and make it win by construction. The replication therefore runs on the
    three ITPR pairs only, and the RyR control lives in the human test
    where the three loci can be told apart by the annotation.
    """
    flanks_path = L.PROJECT_ROOT / "results" / "synteny" / "flanks.tsv"
    if not flanks_path.exists():
        return [SPECIES_HEADER], [SPECIES_SUM_HEADER]
    by: dict[tuple, set[str]] = collections.defaultdict(set)
    meta: dict[str, tuple] = {}
    gidx = L.genome_index()
    for r in L.read_tsv(flanks_path):
        if r.get("informative") != "1" or r.get("window") != "fixed10":
            continue
        sym = (r.get("relaxed_key") or r.get("symbol") or "").upper()
        if sym and not sym.startswith(P.FAMILY_PREFIXES):
            by[(r["accession"], r["cell"])].add(sym)
            g = gidx.get(r["accession"], {})
            meta[r["accession"]] = (r["organism"], r["vclass"],
                                    g.get("vorder", ""))
    rows = []
    for acc, (organism, vclass, vorder) in sorted(meta.items(),
                                                  key=lambda kv: kv[1][0]):
        for a, b in (("ITPR1", "ITPR2"), ("ITPR1", "ITPR3"),
                     ("ITPR2", "ITPR3")):
            sa = {sym2ensg[s] for s in by.get((acc, a), ()) if s in sym2ensg}
            sb = {sym2ensg[s] for s in by.get((acc, b), ()) if s in sym2ensg}
            if len(sa) < MIN_MAPPED or len(sb) < MIN_MAPPED:
                continue
            links = P.links_between(sorted(sa), sb, pairs)
            dated = [x for x in links if x[2] in P.LEVELS_2R_WINDOW]
            fams = {P.family_root(ensg2sym.get(t, t)) for _, t, _ in links}
            rows.append([acc, organism, vclass, vorder, f"{a}_vs_{b}",
                         "ITPR_vs_ITPR", len(sa), len(sb), len(links),
                         len(dated), len(fams),
                         ";".join(sorted(f"{ensg2sym.get(x, x)}-"
                                         f"{ensg2sym.get(y, y)}"
                                         for x, y, _ in links))[:400]])
    summary = []
    by_pair: dict[str, list[list]] = collections.defaultdict(list)
    for r in rows:
        by_pair[r[4]].append(r)
    for pair, rs in sorted(by_pair.items()):
        with_link = [r for r in rs if r[8] > 0]
        with_2r = [r for r in rs if r[9] > 0]
        classes = sorted({r[2] for r in with_link})
        summary.append([pair, rs[0][5], len(rs), len(with_link),
                        round(len(with_link) / len(rs), 4) if rs else 0.0,
                        len(with_2r),
                        round(len(with_2r) / len(rs), 4) if rs else 0.0,
                        round(sum(r[8] for r in rs) / len(rs), 4) if rs else 0.0,
                        max((r[8] for r in rs), default=0), len(classes),
                        ";".join(classes)])
        log(f"[quartet] cross-species {pair}: {len(with_link)}/{len(rs)} "
            f"genomes show >=1 paralog link, in {len(classes)} classes")
    return [SPECIES_HEADER] + rows, [SPECIES_SUM_HEADER] + summary


CTRL_HEADER = ["accession", "organism", "vclass", "window_a", "window_b",
               "n_mapped_a", "n_mapped_b", "n_links"]


def replication_null(accessions: list[str], sym2ensg: dict[str, str],
                     pairs: dict[str, list[tuple]], log=print
                     ) -> tuple[list[list], dict]:
    """The replication's own null: random neighbourhoods in the same genomes.

    A genome where 80 % of ITPR1's flanks have a paralog somewhere in ITPR2's
    flanks proves nothing until the same genome's *random* windows have been
    asked the same question, with the same map, the same symbol
    normalisation and the same mapping loss. Windows come from
    `s8_control.sample_windows` unchanged, so the null S16 reads its
    replication against is the null S8 built and calibrated.
    """
    gidx = L.genome_index()
    rows, hit, total = [], 0, 0
    for acc in accessions:
        g = gidx.get(acc, {})
        genes = F.load_genes(acc)
        if not genes:
            continue
        index = F.index_by_contig(genes)
        wins = C.sample_windows(index, acc, g.get("organism", ""),
                                g.get("vclass", ""), n_rep=N_CONTROL_WINDOWS)
        ids = []
        for w in wins:
            keys = {k for k in w.keys_relaxed
                    if not k.startswith(P.FAMILY_PREFIXES)}
            ids.append({sym2ensg[k] for k in keys if k in sym2ensg})
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                if len(ids[i]) < MIN_MAPPED or len(ids[j]) < MIN_MAPPED:
                    continue
                n = len(P.links_between(sorted(ids[i]), ids[j], pairs))
                total += 1
                hit += int(n > 0)
                rows.append([acc, g.get("organism", ""), g.get("vclass", ""),
                             f"ctrl{i}", f"ctrl{j}", len(ids[i]),
                             len(ids[j]), n])
    frac = hit / total if total else 0.0
    log(f"[quartet] replication null: {hit}/{total} random window pairs "
        f"({frac:.1%}) carry >=1 paralog link")
    return [CTRL_HEADER] + rows, {"n_pairs": total, "n_with_link": hit,
                                  "frac": frac}


def _fisher_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact: is a/(a+b) greater than c/(c+d)?

    Stdlib, because the replication's whole claim is a 2x2 table and adding
    a third-party version number inside it buys nothing.
    """
    import math
    n = a + b + c + d
    if n == 0 or (a + c) == 0:
        return 1.0

    def hyper(k):
        return (math.comb(a + b, k) * math.comb(c + d, a + c - k)
                / math.comb(n, a + c))
    lo = max(0, a + c - (c + d))
    hi = min(a + b, a + c)
    return max(0.0, min(1.0, sum(hyper(k) for k in range(a, hi + 1))
                        if a >= lo else 1.0))


def run(out_dir: Path | None = None, par: dict | None = None,
        log=print) -> dict:
    out = out_dir or L.out_dir()
    rng = random.Random(P.SEED + 1)
    mapping = M.load(out)
    sym2ensg = {s: g for s, g in mapping["sym2ensg"].items() if g}
    pairs, ensg2sym = mapping["pairs"], mapping["ensg2sym"]
    if par is None:
        loci = P.human_loci()
        groups = list(P.TEST_CELLS) + P.control_names(loci)
        by_contig = P.contig_order(P.human_genes())
        windows = {(g, n): P.window_genes(by_contig, *loci[g][:3], n)
                   for g in groups for n in P.WINDOW_SIZES}
        blocks = [P.BLOCK_HEADER] + [[r[c] for c in P.BLOCK_HEADER]
                                     for r in L.read_tsv(
                                         out / "paralogon_blocks.tsv")]
    else:
        loci, groups = par["loci"], par["groups"]
        by_contig, windows, blocks = (par["by_contig"], par["windows"],
                                      par["blocks"])
    pooled = pooled_test(by_contig, windows, sym2ensg, ensg2sym, pairs,
                         groups, rng, log)
    L.write_tsv(out / "pooled_test.tsv", pooled[0], pooled[1:])
    q = quartet_test(by_contig, windows, sym2ensg, blocks, loci, groups,
                     rng, log)
    L.write_tsv(out / "quartet_test.tsv", q[0], q[1:])
    sp, spsum = cross_species(sym2ensg, pairs, ensg2sym, log)
    L.write_tsv(out / "paralogon_species.tsv", sp[0], sp[1:])
    accs = sorted({r[0] for r in sp[1:]})
    ctrl, cstat = replication_null(accs, sym2ensg, pairs, log)
    L.write_tsv(out / "replication_null.tsv", ctrl[0], ctrl[1:])
    for row in spsum[1:]:
        obs_hit, obs_n = row[3], row[2]
        row.extend([cstat["n_pairs"], cstat["n_with_link"],
                    round(cstat["frac"], 4),
                    L.fmt(_fisher_greater(obs_hit, obs_n - obs_hit,
                                          cstat["n_with_link"],
                                          cstat["n_pairs"]
                                          - cstat["n_with_link"]))])
    spsum[0].extend(["null_pairs", "null_with_link", "null_frac",
                     "p_fisher_vs_null"])
    L.write_tsv(out / "paralogon_species_summary.tsv", spsum[0], spsum[1:])
    return {"quartet": q, "pooled": pooled, "species": sp,
            "species_summary": spsum}


def main() -> None:
    run()


if __name__ == "__main__":
    main()
