"""S16 stage `paralogon` — are ITPR1/2/3 a 2R quartet?

S8 measured microsynteny between the three neighbourhoods and got **zero**
shared symbols in every direction, which is what a 2R signal looks like to a
test that is looking for the same word twice. The 2R test asks the next
question: are the neighbours *paralogs of each other*? A whole-genome
duplication copies a chromosomal block, so the block around ITPR1 should be
studded with paralogs of the genes around ITPR2 and ITPR3, far above what
two unrelated windows would share.

Everything is anchored on the human genome, where Compara paralogy is best,
and the windows are cut out of the sweep's own `genes_slim.tsv` — the file
S8 cut its flanks from — so the two tasks cannot disagree about what a
neighbourhood is.

**The RyR trio is the positive control, and it is what makes the answer
readable.** RYR1/2/3 are a three-member vertebrate family of the same age
and the same 2R candidacy; D14 puts them inside every search anyway. Running
the identical test on the identical instrument in the same genome gives a
measured yardstick for what a 2R quartet scores here — so an ITPR result can
be compared against something other than zero.

Three tests:

1. **Cross-window paralogy** at +/-10, 20 and 30 genes, against a
   permutation null of **real genomic windows** of the same gene count
   (D17) with the family's own genes removed from every window (or the
   test is circular), plus a Poisson cross-check.
2. **Block scan.** Slide the same window genome-wide and rank every block by
   how many distinct gene *families* it shares with each ITPR window. This
   is where a fourth 2R slot — a paralogous block that no longer carries an
   ITPR gene — would show up, without assuming where it is.
3. **Quartet + cross-species replication** — `s16_quartet.py`.

Outputs (`results/duplication/`): `windows.tsv`, `paralogy_links.tsv`,
`paralogon_test.tsv`, `paralogon_blocks.tsv`.
"""

from __future__ import annotations

import collections
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s16_lib as L                                            # noqa: E402
import s16_ensembl_map as M                                    # noqa: E402
import s8_flank_lib as F                                       # noqa: E402

HUMAN = "GCF_000001405.40"
WINDOW_SIZES = (10, 20, 30)
N_PERMUTATIONS = 5000
SEED = 20260908

#: Every gene of both families is dropped from every window, not just the
#: three under test. Counting the ITPR1-ITPR2 paralogy itself as evidence
#: that their blocks are paralogous is circular — that pair is the thing
#: being explained — and leaving a RyR in an ITPR window would import the
#: control's answer into the test.
FAMILY_PREFIXES = ("ITPR", "RYR")

#: Compara puts a **duplication node** on every paralog pair — the last
#: common ancestor at which the two copies split. That column is what tells a
#: 2R ohnolog from an older duplication whose two copies happen to sit in
#: these blocks, and it is the difference between a paralogon test and a
#: paralogy test. Two nested vocabularies, both reported, so the answer does
#: not rest on where one line is drawn: `2R_core` is the vertebrate stem
#: itself, `2R_window` adds the two nodes immediately below it, where
#: Compara's LCA reconstruction puts a real 2R pair whenever the deep
#: outgroup genes are missing from its tree.
LEVELS_2R_CORE = ("Chordata", "Vertebrata")
LEVELS_2R_WINDOW = LEVELS_2R_CORE + ("Gnathostomata", "Euteleostomi")
LEVEL_SETS = {
    "all_levels": None,
    "2R_window": set(LEVELS_2R_WINDOW),
    "2R_core": set(LEVELS_2R_CORE),
}

#: Human loci, straight out of the S5 sweep's own `summary.json` rather than
#: typed from memory. Filled by `human_loci()`; the constant records which
#: cells are the test and which are the control.
TEST_CELLS = ("ITPR1", "ITPR2", "ITPR3")
CONTROL_GROUP = "RYR"


def human_loci() -> dict[str, tuple]:
    """name -> (contig, start, end, note), from the sweep's human summary.

    The three ITPR cells hold one locus each; the RyR *cell* holds three
    loci, which the sweep's own annotation names RYR1/RYR2/RYR3 — so the
    control's three windows are read off the assembly's annotation and are
    not a list this module carries.
    """
    recs = [r for r in L.iter_loci([HUMAN]) if L.is_copy(r)]
    out: dict[str, tuple] = {}
    for r in recs:
        if r["cell"] in TEST_CELLS:
            out[r["cell"]] = (r["contig"], r["start"], r["end"],
                              r["annot_gene"] or r["cell"])
        elif r["cell"] == CONTROL_GROUP:
            name = r["annot_gene"] or f"{CONTROL_GROUP}?"
            out[name] = (r["contig"], r["start"], r["end"], name)
    if not all(c in out for c in TEST_CELLS):
        raise SystemExit(f"human sweep summary is missing an ITPR cell: "
                         f"{sorted(out)}")
    return out


def control_names(loci: dict[str, tuple]) -> list[str]:
    return sorted(k for k in loci if k not in TEST_CELLS)


# --------------------------------------------------------------- windows

def human_genes() -> list[tuple]:
    genes = F.load_genes(HUMAN)
    if not genes:
        raise SystemExit(f"no genes_slim.tsv for {HUMAN}")
    return [g for g in genes if g[6] == "protein_coding"]


def contig_order(genes: list[tuple]) -> dict[str, list[tuple]]:
    by: dict[str, list[tuple]] = collections.defaultdict(list)
    for g in genes:
        by[g[0]].append(g)
    for v in by.values():
        v.sort(key=lambda g: (g[1], g[2]))
    return dict(by)


def window_genes(by_contig: dict[str, list[tuple]], contig: str,
                 start: int, end: int, n: int) -> list[tuple]:
    """The n nearest coding genes each side, plus anything overlapping,
    with every ITPR and RyR gene removed (see FAMILY_PREFIXES)."""
    genes = [g for g in by_contig.get(contig, [])
             if not g[5].upper().startswith(FAMILY_PREFIXES)]
    up = [g for g in genes if g[2] < start]
    down = [g for g in genes if g[1] > end]
    over = [g for g in genes if g[2] >= start and g[1] <= end]
    return up[-n:] + over + down[:n]


def human_window_symbols(n: int = max(WINDOW_SIZES)) -> set[str]:
    by = contig_order(human_genes())
    syms: set[str] = set()
    for contig, start, end, _ in human_loci().values():
        for g in window_genes(by, contig, start, end, n):
            syms.add(g[5].upper())
    return syms


# ------------------------------------------------------------ link counts

def links_between(a_ids: list[str], b_ids: set[str],
                  pairs: dict[str, list[tuple]],
                  levels: set[str] | None = None) -> list[tuple]:
    out = []
    for gid in a_ids:
        for tid, level in pairs.get(gid, []):
            if tid in b_ids and tid != gid:
                if levels is not None and level not in levels:
                    continue
                out.append((gid, tid, level))
    return dedupe_links(out)


def dedupe_links(links: list[tuple]) -> list[tuple]:
    seen, out = set(), []
    for a, b, lev in links:
        key = tuple(sorted((a, b)))
        if key in seen:
            continue
        seen.add(key)
        out.append((a, b, lev))
    return out


def poisson_p(observed: int, expected: float) -> float:
    """P(X >= observed) for X ~ Poisson(expected)."""
    if expected <= 0:
        return 1.0 if observed == 0 else 0.0
    cum, term = 0.0, math.exp(-expected)
    for k in range(observed):
        cum += term
        term *= expected / (k + 1)
    return max(0.0, min(1.0, 1.0 - cum))


def paralogs_of(ids: list[str], pairs: dict[str, list[tuple]],
                levels: set[str] | None = None) -> set[str]:
    """Every gene paralogous to one of `ids`, at the stated duplication
    nodes. The null must be built from the same link set as the test, or
    the dated test is scored against an undated background."""
    return {t for g in ids for t, lev in pairs.get(g, [])
            if levels is None or lev in levels}


def permutation_null(a_paralogs: set[str], b_size: int,
                     all_windows: list[list[str]], rng: random.Random,
                     n: int = N_PERMUTATIONS) -> list[int]:
    """How many of window A's paralogs land in a *real* random window of |B|?

    D17: the windows are drawn from the real gene order, so the null keeps
    the clustering of gene families that a shuffled gene set throws away.
    That is the conservative choice — a tandem array inflates the link count
    in the real window and in the null alike.
    """
    counts = []
    pool = [w for w in all_windows if len(w) == b_size] or all_windows
    for _ in range(n):
        w = rng.choice(pool)
        counts.append(sum(1 for g in w if g in a_paralogs))
    return counts


def sliding_windows(by_contig: dict[str, list[tuple]], size: int,
                    sym2ensg: dict[str, str], step: int = 1) -> list[list[str]]:
    """Every sliding window of `size` genes genome-wide, as Ensembl ids."""
    out = []
    for genes in by_contig.values():
        ids = [sym2ensg.get(g[5].upper(), "") for g in genes]
        for i in range(0, max(0, len(ids) - size + 1), step):
            out.append([x for x in ids[i:i + size] if x])
    return out


def family_root(symbol: str) -> str:
    """The alphabetic prefix before the first digit: ZNF763 -> ZNF,
    SLC22A10 -> SLC, GNAS -> GNAS. Crude, and only ever used to collapse a
    tandem array into one family so nine adjacent zinc fingers score 1."""
    out = []
    for ch in symbol.upper():
        if ch.isdigit():
            break
        out.append(ch)
    return "".join(out) or symbol.upper()


# ------------------------------------------------------------------- run

WIN_HEADER = ["group", "window_n", "locus_note", "contig", "start", "end",
              "strand", "symbol", "ensembl_gene", "map_status",
              "index_in_window"]
LINK_HEADER = ["window_n", "group_a", "group_b", "pair_class", "symbol_a",
               "symbol_b", "gene_a", "gene_b", "duplication_level",
               "in_2R_window", "in_2R_core"]
TEST_HEADER = ["level_set", "window_n", "group_a", "group_b", "pair_class",
               "n_genes_a", "n_genes_b", "observed_links", "n_families",
               "null_mean", "enrichment", "p_permutation", "p_poisson",
               "null_max", "n_permutations", "q_stratum", "q_global"]


def pair_class(a: str, b: str) -> str:
    a_itpr, b_itpr = a in TEST_CELLS, b in TEST_CELLS
    if a_itpr and b_itpr:
        return "ITPR_vs_ITPR"
    if not a_itpr and not b_itpr:
        return "RYR_vs_RYR (positive control)"
    return "ITPR_vs_RYR (cross-family)"


def run(out_dir: Path | None = None, log=print) -> dict:
    out = out_dir or L.out_dir()
    rng = random.Random(SEED)
    mapping = M.load(out)
    sym2ensg = {s: g for s, g in mapping["sym2ensg"].items() if g}
    ensg2sym = mapping["ensg2sym"]
    pairs = mapping["pairs"]
    if not pairs:
        raise SystemExit("run the `map` stage first (paralogy_map.tsv empty)")

    loci = human_loci()
    groups = list(TEST_CELLS) + control_names(loci)
    genes = human_genes()
    by_contig = contig_order(genes)
    log(f"[paralogon] human background: {len(genes)} coding genes on "
        f"{len(by_contig)} contigs; groups {', '.join(groups)}")

    win_rows, test_rows, link_rows = [], [], []
    windows: dict[tuple, list[tuple]] = {}
    for n in WINDOW_SIZES:
        for g in groups:
            contig, start, end, note = loci[g]
            wg = window_genes(by_contig, contig, start, end, n)
            windows[(g, n)] = wg
            for rank, gene in enumerate(wg):
                gid = sym2ensg.get(gene[5].upper(), "")
                win_rows.append([g, n, note, gene[0], gene[1], gene[2],
                                 gene[3], gene[5], gid,
                                 "mapped" if gid else "unmapped", rank])

    for n in WINDOW_SIZES:
        sliding = sliding_windows(by_contig, 2 * n + 1, sym2ensg, step=1)
        for i in range(len(groups)):
            for j in range(i + 1, len(groups)):
                a, b = groups[i], groups[j]
                wa = [x for x in (sym2ensg.get(g[5].upper(), "")
                                  for g in windows[(a, n)]) if x]
                wb = [x for x in (sym2ensg.get(g[5].upper(), "")
                                  for g in windows[(b, n)]) if x]
                cls = pair_class(a, b)
                for ga, gb, lev in links_between(wa, set(wb), pairs):
                    link_rows.append([n, a, b, cls, ensg2sym.get(ga, ga),
                                      ensg2sym.get(gb, gb), ga, gb, lev,
                                      int(lev in LEVELS_2R_WINDOW),
                                      int(lev in LEVELS_2R_CORE)])
                for set_name, levels in LEVEL_SETS.items():
                    links = links_between(wa, set(wb), pairs, levels)
                    fams = {family_root(ensg2sym.get(t, t))
                            for _, t, _ in links}
                    null = (permutation_null(paralogs_of(wa, pairs, levels),
                                             len(wb), sliding, rng)
                            + permutation_null(paralogs_of(wb, pairs, levels),
                                               len(wa), sliding, rng))
                    obs = len(links)
                    p_perm = ((sum(1 for x in null if x >= obs) + 1)
                              / (len(null) + 1))
                    mean_null = sum(null) / len(null)
                    test_rows.append([
                        set_name, n, a, b, cls, len(wa), len(wb), obs,
                        len(fams), round(mean_null, 4),
                        round(obs / mean_null, 3) if mean_null > 0 else "inf",
                        L.fmt(p_perm), L.fmt(poisson_p(obs, mean_null)),
                        max(null), len(null)])
                    if set_name == "2R_window":
                        log(f"[paralogon] +/-{n:>2} {a:>5} vs {b:<5} "
                            f"2R-dated links={obs:<3} fams={len(fams):<3} "
                            f"null={mean_null:.2f} p={p_perm:.1e}")

    # BH twice, because the two families of tests answer different
    # questions and reporting only one of them would be a choice made after
    # seeing the numbers (S9's rule, stated at both scopes). `q_stratum`
    # corrects across the 15 pairs of one (level set, window) — the family a
    # reader actually reads a row against, since the three level sets are
    # nested and the three windows are nested. `q_global` corrects across
    # every test this stage ran, which is the conservative bound.
    strata: dict[tuple, list[int]] = {}
    for i, r in enumerate(test_rows):
        strata.setdefault((r[0], r[1]), []).append(i)
    q_str = [1.0] * len(test_rows)
    for idxs in strata.values():
        for i, q in zip(idxs, L.benjamini_hochberg(
                [float(test_rows[i][11]) for i in idxs])):
            q_str[i] = q
    q_glob = L.benjamini_hochberg([float(r[11]) for r in test_rows])
    for r, a, b in zip(test_rows, q_str, q_glob):
        r.extend([L.fmt(a), L.fmt(b)])

    L.write_tsv(out / "windows.tsv", WIN_HEADER, win_rows)
    L.write_tsv(out / "paralogy_links.tsv", LINK_HEADER, link_rows)
    L.write_tsv(out / "paralogon_test.tsv", TEST_HEADER, test_rows)

    blocks = block_scan(by_contig, windows, sym2ensg, ensg2sym, pairs,
                        loci, groups, log)
    L.write_tsv(out / "paralogon_blocks.tsv", blocks[0], blocks[1:])
    return {"tests": test_rows, "links": link_rows, "blocks": blocks,
            "windows": windows, "by_contig": by_contig, "groups": groups,
            "loci": loci}


# -------------------------------------------------------------- block scan

BLOCK_HEADER = ["rank", "against", "contig", "start", "end", "n_genes",
                "n_families", "n_paralog_hits", "n_query_genes_hit",
                "carries_family_gene", "hit_symbols"]


def block_scan(by_contig: dict[str, list[tuple]],
               windows: dict[tuple, list[tuple]],
               sym2ensg: dict[str, str], ensg2sym: dict[str, str],
               pairs: dict[str, list[tuple]], loci: dict[str, tuple],
               groups: list[str], log=print,
               n: int = 20, top: int = 25) -> list[list]:
    """Rank every genomic block by paralogy to each window.

    Two rounds of whole-genome duplication make **four** copies of an
    ancestral block. Three carry an ITPR gene; if the fourth survives with
    its ITPR deleted it is still a paralogous block, and this is where it
    would appear. The scan is agnostic — it ranks every block genome-wide
    and lets a quartet fall out or not.

    Blocks are ranked by distinct gene **families**, not by raw hit count:
    a tandem array is one duplication, not a paralogon, and nine adjacent
    zinc fingers must score 1.
    """
    size = 2 * n + 1
    spans = {(c, s, e) for c, s, e, _ in loci.values()}
    rows = []
    for p in groups:
        wg = [x for x in (sym2ensg.get(g[5].upper(), "")
                          for g in windows[(p, n)]) if x]
        wset = set(wg)
        par: dict[str, set[str]] = collections.defaultdict(set)
        for g in wg:
            for t, _ in pairs.get(g, []):
                if t not in wset:
                    par[t].add(g)
        scored = []
        for contig, genes in by_contig.items():
            ids = [(g, sym2ensg.get(g[5].upper(), "")) for g in genes]
            for i in range(0, max(0, len(ids) - size + 1)):
                chunk = ids[i:i + size]
                hits = [(g[5], gid) for g, gid in chunk if gid in par]
                if len(hits) < 2:
                    continue
                fams = {family_root(sym) for sym, _ in hits}
                queries = set().union(*(par[gid] for _, gid in hits))
                s, e = chunk[0][0][1], chunk[-1][0][2]
                carries = any(c == contig and not (e < st or s > en)
                              for c, st, en in spans)
                scored.append((len(fams), len(hits), len(queries), contig,
                               s, e, len(chunk), [h[0] for h in hits],
                               carries))
        scored.sort(key=lambda x: (-x[0], -x[2], -x[1], x[3], x[4]))
        kept: list[tuple] = []
        for item in scored:
            if any(item[3] == k[3] and not (item[5] < k[4] or item[4] > k[5])
                   for k in kept):
                continue
            kept.append(item)
            if len(kept) >= top:
                break
        for rank, it in enumerate(kept, 1):
            nf, nh, nq, contig, s, e, ng, hits, carries = it
            rows.append([rank, p, contig, s, e, ng, nf, nh, nq,
                         "yes" if carries else "no", ";".join(sorted(hits))])
        if kept:
            log(f"[paralogon] block scan vs {p}: top block {kept[0][0]} "
                f"families on {kept[0][3]} "
                f"(carries a family gene: {'yes' if kept[0][8] else 'no'})")
    return [BLOCK_HEADER] + rows


def main() -> None:
    run()


if __name__ == "__main__":
    main()
