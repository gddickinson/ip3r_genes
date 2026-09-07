#!/usr/bin/env python3
"""S8's null: matched random neighbourhoods in the same genomes.

A within-paralog Jaccard of 0.4 means nothing on its own. Two neighbourhoods
in two well-annotated mammals share vocabulary for reasons that have nothing
to do with orthology, and two neighbourhoods in a lamprey and a hagfish
share almost none whatever their history, because 27 % of their coding genes
carry a symbol at all.

So every real pair (locus in genome A, locus in genome B) is scored against
matched control pairs (random coding gene in A, random coding gene in B),
drawn with the same window rule and the same key rule. The control therefore
holds constant: the two genomes, their annotation depth, their naming
conventions, the window size and the key normalisation. What is left is the
locus.

Sampling is seeded per accession, so the same genome yields the same control
windows in every run and in every replicate ordering (D24's reproducibility
discipline applied to a random control).
"""
from __future__ import annotations

import random
from collections import defaultdict
from math import sqrt
from statistics import mean, median

import s8_flank_lib as lib

N_REPLICATES = 3        # control windows drawn per genome
CONTROL_SEED = 20260907


def sample_windows(index: dict, accession: str, organism: str, vclass: str,
                   n_rep: int = N_REPLICATES, n: int = lib.DEFAULT_N,
                   window: str = "fixed10") -> list[lib.FlankSet]:
    """n_rep random coding genes in this genome, each with its own flank set.

    The drawn gene plays the part of the locus: it is excluded from its own
    flanks exactly as an ITPR gene is. Contigs carrying fewer than 2n+1
    coding genes are skipped, because a window truncated by a contig end is
    not a fair match for a locus on a chromosome -- the real loci that sit
    on short contigs are handled by reporting their own n_flanks, not by
    degrading the null.
    """
    usable = [(c, gs) for c, gs in index.items() if len(gs) >= 2 * n + 1]
    if not usable:
        return []
    rng = random.Random(f"{CONTROL_SEED}:{accession}")
    total = sum(len(gs) - 2 * n for c, gs in usable)
    out = []
    for r in range(n_rep):
        # weight contigs by how many interior positions they offer
        pick = rng.randrange(total)
        for c, gs in usable:
            room = len(gs) - 2 * n
            if pick < room:
                g = gs[n + pick]
                break
            pick -= room
        loc = lib.Locus(
            accession=accession, organism=organism, vclass=vclass, vorder="",
            cell=f"control{r}", status="control", contig=g[0],
            start=g[1], end=g[2], idx=r, bait="", bait_paralog="",
            identity=0.0, coverage=0.0, annot_gene=g[5], annot_paralog="")
        out.append(lib.extract_flanks(index, loc, n=n, window=window))
    return out


def matched_null(pairs: list[tuple], controls: dict[str, list[lib.FlankSet]],
                 kind: str = "relaxed") -> tuple[list[float], list[float]]:
    """For each real pair, the Jaccard of the matched control pairs.

    Returns (observed J per pair, mean control J per pair) over the pairs for
    which both genomes have control windows.
    """
    obs, null = [], []
    for _, _, j, _, acc_a, acc_b in pairs:
        ca, cb = controls.get(acc_a), controls.get(acc_b)
        if not ca or not cb:
            continue
        js = [lib.jaccard(lib.keys_of(x, kind), lib.keys_of(y, kind))
              for x, y in zip(ca, cb)]
        if not js:
            continue
        obs.append(j)
        null.append(mean(js))
    return obs, null


def sign_test(obs: list[float], null: list[float]) -> dict:
    """How often does the real pair beat its own matched control pair?

    A sign test rather than a t-test: Jaccard is bounded, zero-inflated and
    nowhere near normal, and the paired design already removes the genome
    effect the difference would otherwise be dominated by. Ties (both zero,
    which is common in sparse annotations) are dropped and counted, because
    counting them as failures would make an underpowered comparison look
    like a negative result.
    """
    up = sum(1 for a, b in zip(obs, null) if a > b)
    down = sum(1 for a, b in zip(obs, null) if a < b)
    ties = len(obs) - up - down
    n = up + down
    if n == 0:
        return dict(n_pairs=len(obs), n_used=0, n_ties=ties, frac_up=0.0, z=0.0)
    frac = up / n
    z = (abs(up - n / 2) - 0.5) / sqrt(n / 4)
    return dict(n_pairs=len(obs), n_used=n, n_ties=ties, frac_up=frac,
                z=z if up >= down else -z)


def describe(js: list[float]) -> dict:
    if not js:
        return dict(n=0, mean=0.0, median=0.0, max=0.0, frac_pos=0.0)
    return dict(n=len(js), mean=mean(js), median=median(js), max=max(js),
                frac_pos=sum(1 for j in js if j > 0) / len(js))


def pair_class_row(name: str, pairs: list[tuple],
                   controls: dict[str, list[lib.FlankSet]],
                   kind: str = "relaxed") -> dict:
    """One row of pair_stats.tsv: the class, its matched null, and the test."""
    obs, null = matched_null(pairs, controls, kind=kind)
    o, nl = describe(obs), describe(null)
    st = sign_test(obs, null)
    return dict(
        pair_class=name, key=kind, n_pairs=o["n"],
        mean_j=o["mean"], median_j=o["median"], max_j=o["max"],
        frac_positive=o["frac_pos"],
        control_mean_j=nl["mean"], control_median_j=nl["median"],
        control_frac_positive=nl["frac_pos"],
        excess_mean=o["mean"] - nl["mean"],
        ratio=(o["mean"] / nl["mean"]) if nl["mean"] else float("inf"),
        n_used=st["n_used"], n_ties=st["n_ties"],
        frac_beats_control=st["frac_up"], z=st["z"])


def background_prevalence(controls: dict[str, list[lib.FlankSet]],
                          kind: str = "root") -> tuple[dict[str, float], int]:
    """How often each key appears in a random neighbourhood.

    This is what tells a paralogon signal from a vocabulary artefact: a root
    shared between two ITPR neighbourhoods is only evidence if it is not in
    every third random window as well.
    """
    counts: defaultdict = defaultdict(int)
    n = 0
    for sets_ in controls.values():
        for fs in sets_:
            keys = lib.keys_of(fs, kind)
            if not keys:
                continue
            n += 1
            for k in keys:
                counts[k] += 1
    return ({k: c / n for k, c in counts.items()} if n else {}), n


def stratify_pairs(pairs: list[tuple], vclass_of: dict[str, str]
                   ) -> dict[str, list[tuple]]:
    """Split a pair list into same-class and cross-class comparisons.

    A within-paralog mean Jaccard pools two very different questions: do two
    mammals share the neighbourhood (they do, trivially), and does a mammal
    share it with a teleost (the question about the locus). Reporting only
    the pooled number lets a clade-restricted signal read as a vertebrate-
    wide one.
    """
    same, cross = [], []
    for p in pairs:
        a, b = vclass_of.get(p[0], ""), vclass_of.get(p[1], "")
        (same if (a and a == b) else cross).append(p)
    return {"same_vclass": same, "cross_vclass": cross}
