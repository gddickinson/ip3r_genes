"""s15_integrity.py — is a recovered gene's reading frame broken, or is its
assembly broken?

The sweep records, per locus, how many frameshifts and in-frame stop
codons miniprot had to accommodate to align a reference through it.  Read
naively that is a pseudogene screen: a gene carrying 43 frameshifts is
dead.  Read honestly it is mostly a *sequencing* statistic — a frameshift
is one indel, and an assembly built from short reads at modest coverage
puts indels in coding sequence at a rate that has nothing to do with
whether the gene works.  The family's own longest member makes the point:
the RyR control cell, a 5,000-residue gene nobody claims is dead, is 51
of the 230 loci in this sweep carrying any lesion at all.

So the density is measured and then three confounders are controlled, all
three on the sweep's own committed numbers.

**Length.** Lesions are counted per kilo-aligned-residue, never per gene,
because the RyR reference is 1.8x the ITPR references and would otherwise
lead every ranking by construction.

**The assembly, and the alignment.**  `covariates()` is Spearman of
density against the genome's own contig N50, the locus's aligned length
and — the one that turned out to matter — the locus's percent identity to
its bait.  Contiguity is nearly irrelevant here (rho = -0.08) and
identity is not (rho = -0.40): a lesion count is substantially a measure
of how far the reference is from the gene, because a poorly matched bait
buys alignment with frameshifts.  Any statement about lesions has to
survive that, which is why the paired test below is also run
identity-matched.

**The genome.** `paired_within_genome()` is the control the brief asks
for and the only one that removes the assembly entirely: each ITPR cell's
density is compared against the *same genome's* other family loci, and
the paired differences go to a sign test (ties dropped and counted,
S8's rule).  A genome-wide indel rate cancels in the difference, so what
survives is a lesion excess specific to that paralog in that species —
which is what a pseudogene would look like and a bad assembly would not.

The bar is calibrated on a population the screen never scores: loci at
full coverage, in assemblies above D4's contiguity bar, whose own
annotation names the gene as a member of this family.  Those are genes a
second pipeline independently calls functional, so their lesion density
is what an intact gene looks like on this instrument, and the bar is
their upper quantile.  A locus above it is `elevated_lesions` — a
candidate, never a verdict, because `s10_orf.py` already established
that the one-sided direction is the only sound one: zero stops falsifies
a pseudogene call, a handful does not establish one.
"""

from __future__ import annotations

import collections

import s15_lib as lib
import s5_calibration as s5cal

#: the quantile of the intact population that sets the bar
BAR_QUANTILE = 0.99
#: coverage at or above which a locus is scored (the sweep's own COV_FOUND)
COV_FULL = 0.70
#: a locus must align at least this many residues to be scored at all
MIN_ALIGNED_AA = 500


def _q(xs: list[float], p: float) -> float:
    s = sorted(xs)
    if not s:
        return float("nan")
    i = min(len(s) - 1, max(0, int(round(p * (len(s) - 1)))))
    return s[i]


def loci_table(summaries: dict[str, dict], manifest: dict[str, dict]
               ) -> list[dict]:
    """One row per locus the sweep placed, with its lesion density.

    Contiguity comes from `s5_calibration.spans_a_gene` — D4's own bar,
    computed from the committed span measurement — and not from the
    ledger's per-cell `contig_spans_gene` column, which is a property of
    that cell's best locus and cannot be keyed by genome.
    """
    rows: list[dict] = []
    for acc, s in summaries.items():
        n50 = int(manifest.get(acc, {}).get("contig_n50") or 0)
        for cell, c in (s.get("cells") or {}).items():
            for i, loc in enumerate(c.get("loci", []) or []):
                aa = int(loc.get("aligned_aa") or 0)
                if aa < MIN_ALIGNED_AA:
                    continue
                fs = int(loc.get("frameshifts") or 0)
                st = int(loc.get("stop_codons") or 0)
                dens = (fs + st) / (aa / 1000.0)
                rows.append(dict(
                    accession=acc, organism=s.get("organism", ""),
                    vclass=s.get("vclass", ""), cell=cell, locus_idx=i,
                    contig=loc.get("contig", ""),
                    coverage=round(float(loc.get("coverage") or 0.0), 4),
                    identity=round(float(loc.get("identity") or 0.0), 4),
                    aligned_aa=aa, frameshifts=fs, stop_codons=st,
                    lesions=fs + st, lesion_density=round(dens, 4),
                    contig_edge=int(bool(loc.get("contig_edge"))),
                    n_gap=int(bool(loc.get("n_gap"))),
                    longest_n_run=int(loc.get("longest_n_run") or 0),
                    annot_paralog_matches=int(
                        bool(loc.get("annot_paralog_matches"))),
                    annot_gene=(loc.get("annot_gene") or {}).get("name", "")
                    if isinstance(loc.get("annot_gene"), dict)
                    else (loc.get("annot_gene") or ""),
                    contig_n50=n50,
                    contig_spans_gene=int(s5cal.spans_a_gene(n50)),
                    scored=int(float(loc.get("coverage") or 0.0) >= COV_FULL),
                ))
    return rows


def intact_population(rows: list[dict]) -> list[dict]:
    """The calibration set: full coverage, contiguous assembly, and the
    assembly's *own* annotation names the gene as family."""
    return [r for r in rows if r["scored"] and r["contig_spans_gene"]
            and r["annot_paralog_matches"]]


def calibrate_bar(rows: list[dict]) -> dict:
    pop = intact_population(rows)
    dens = [r["lesion_density"] for r in pop]
    bar = _q(dens, BAR_QUANTILE)
    return dict(n_intact=len(pop), quantile=BAR_QUANTILE, bar=bar,
                intact_median=lib.median(dens),
                intact_max=max(dens) if dens else float("nan"),
                n_zero=sum(1 for d in dens if d == 0),
                frac_zero=(sum(1 for d in dens if d == 0) / len(dens))
                if dens else float("nan"))


#: how close two loci's bait identities must be to be an identity-matched pair
IDENTITY_WINDOW = 0.02

_COVARIATES = (("contig_n50", lambda r: float(r["contig_n50"])),
               ("identity", lambda r: r["identity"]),
               ("aligned_aa", lambda r: float(r["aligned_aa"])))


def covariates(rows: list[dict]) -> list[dict]:
    """Spearman of lesion density against each confounder, overall and per
    cell — the three things that could produce a lesion count without a
    gene being dead."""
    out = []
    scored = [r for r in rows if r["scored"] and r["contig_n50"] > 0]
    for name, fn in _COVARIATES:
        for label, g in [("all", scored)] + [
                (c, [r for r in scored if r["cell"] == c])
                for c in lib.ALL_CELLS]:
            rho, p, n = lib.spearman([fn(r) for r in g],
                                     [r["lesion_density"] for r in g])
            out.append(dict(covariate=name, subset=label, n=n, rho=rho, p=p,
                            median_density=lib.median(
                                [r["lesion_density"] for r in g])))
    return out


def paired_within_genome(rows: list[dict], match_identity: bool = False
                         ) -> tuple[list[dict], list[dict]]:
    """Each ITPR cell against its own genome's other family loci.

    The comparison is within one assembly, so a genome-wide indel rate
    cancels: what a positive difference means is that *this* paralog in
    *this* species carries more lesions than its siblings do, which is
    the only form the pseudogene hypothesis can take here.
    """
    by_acc: dict[str, dict[str, tuple]] = collections.defaultdict(dict)
    for r in rows:
        if not r["scored"]:
            continue
        # the genome's best-covered locus per cell, as the ledger does
        cur = by_acc[r["accession"]].get(r["cell"])
        if cur is None or r["coverage"] > cur[0]:
            by_acc[r["accession"]][r["cell"]] = (r["coverage"],
                                                 r["lesion_density"],
                                                 r["identity"])
    pairs: list[dict] = []
    for acc, cells in by_acc.items():
        if len(cells) < 2:
            continue
        for cell, (_, dens, ident) in cells.items():
            others = [(d, i) for c, (_, d, i) in cells.items() if c != cell]
            if match_identity:
                others = [(d, i) for d, i in others
                          if abs(i - ident) <= IDENTITY_WINDOW]
            if not others:
                continue
            pairs.append(dict(accession=acc, cell=cell, density=dens,
                              identity=ident,
                              others_median=lib.median([d for d, _ in others]),
                              others_identity_median=lib.median(
                                  [i for _, i in others]),
                              n_others=len(others),
                              matched=int(match_identity),
                              diff=dens - lib.median([d for d, _ in others])))
    tests = []
    for cell in lib.ALL_CELLS:
        g = [p for p in pairs if p["cell"] == cell]
        t = lib.sign_test([p["diff"] for p in g])
        tests.append(dict(cell=cell, matched=int(match_identity),
                          n_genomes=len(g),
                          median_diff=lib.median([p["diff"] for p in g]),
                          direction=("excess" if t["n_pos"] > t["n_neg"]
                                     else "deficit" if t["n_neg"] > t["n_pos"]
                                     else "none"),
                          **t))
    # one family of four tests of the same hypothesis, so the smallest of
    # them is not a p-value (S9's rule, D-corrected the same way)
    for r, q in zip(tests, _bh([t["p"] for t in tests])):
        r["q_bh"] = q
    return pairs, tests


def _bh(ps: list[float]) -> list[float]:
    """Benjamini-Hochberg, stdlib, order preserved."""
    idx = [i for i, p in enumerate(ps) if p == p]      # drop NaN
    m = len(idx)
    out = [float("nan")] * len(ps)
    if not m:
        return out
    order = sorted(idx, key=lambda i: ps[i])
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = m - rank + 1
        prev = min(prev, ps[i] * m / k)
        out[i] = min(1.0, prev)
    return out


def verdicts(rows: list[dict], bar: dict) -> list[dict]:
    """A one-sided candidate flag, with the confounder named when it fires."""
    b = bar["bar"]
    out = []
    for r in rows:
        if not r["scored"]:
            v, why = "not_scored", f"coverage {r['coverage']:.2f} < {COV_FULL}"
        elif r["lesion_density"] <= b:
            v, why = "intact", (f"density {r['lesion_density']:.2f} <= bar "
                                f"{b:.2f}")
        elif not r["contig_spans_gene"]:
            v, why = "assembly_explained", (
                f"density {r['lesion_density']:.2f} > bar {b:.2f} but the "
                f"assembly's contig N50 ({r['contig_n50']:,} bp) is below "
                f"D4's bar for this gene")
        else:
            v, why = "elevated_lesions", (
                f"density {r['lesion_density']:.2f} > bar {b:.2f} in a "
                f"contiguous assembly")
        out.append(dict(r, verdict=v, verdict_reason=why))
    return out


LOCUS_COLS = ["accession", "organism", "vclass", "cell", "locus_idx",
              "contig", "coverage", "identity", "aligned_aa", "frameshifts",
              "stop_codons", "lesions", "lesion_density", "contig_edge",
              "n_gap", "longest_n_run", "annot_gene",
              "annot_paralog_matches", "contig_n50", "contig_spans_gene",
              "scored", "verdict", "verdict_reason"]
PAIR_COLS = ["accession", "cell", "matched", "density", "identity",
             "others_median", "others_identity_median", "n_others", "diff"]
TEST_COLS = ["cell", "matched", "n_genomes", "median_diff", "direction",
             "n", "n_pos", "n_neg", "n_ties", "p", "q_bh"]
CORR_COLS = ["covariate", "subset", "n", "rho", "p", "median_density"]
