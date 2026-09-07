#!/usr/bin/env python3
"""S8's two family-specific questions.

**1. Is there a paralogon?** ITPR1/2/3 are a 2R product, so their
neighbourhoods should be *paralogous*, not identical: the flanking genes of
ITPR1 and ITPR2 should be the two surviving copies of the same ancestral
gene families, carrying different symbols. A symbol Jaccard cannot see that
at all -- it is looking for the same word twice. The root key
(`s8_flank_lib.root_key`, BHLHE40 and BHLHE41 -> BHLHE) is what makes an
ohnolog pair visible to a set comparison, and the random-window background
is what stops a promiscuous root family (ZNF, SLC) from manufacturing one.

**2. Which neighbourhood does an unplaced locus sit in?** The cyclostome
loci S5 filed into the ITPR1 cell -- three in the sea lamprey, three in the
hagfish -- and every fragment or assembly-gap locus. This is scored against
a per-paralog flank consensus built with the query's own genome left out,
and the caller is calibrated on the loci whose paralog identity their own
assembly's annotation already establishes. An assignment rule that cannot
recover the answers already known is not evidence about the ones that are
not.
"""
from __future__ import annotations

from collections import Counter, defaultdict

import s8_flank_lib as lib

# a key is in a paralog's consensus if this fraction of its loci carry it
CONSENSUS_FRAC = 0.50
# a root has to reach this prevalence on both sides to be a paralogon
# candidate. Deliberately low, and the summary counts at three bars, because
# a shared-ohnolog claim should not depend on where one line is drawn.
PARALOGON_MIN_FRAC = 0.10
PARALOGON_BARS = (0.10, 0.25, 0.50)
# a locus needs this many informative keys before it is scored at all
MIN_KEYS = 4


# ------------------------------------------------------------- paralogon

def side_prevalence(sets_: list[lib.FlankSet], kind: str
                    ) -> tuple[dict, int, dict]:
    """Fraction of these loci carrying each key, counting one vote per
    species so a densely-sampled clade cannot carry a root on its own.

    Also returns, per key, the vertebrate classes it was seen in: 62 % of
    species is a different claim if those species are all mammals.
    """
    seen: defaultdict = defaultdict(set)
    classes: defaultdict = defaultdict(set)
    species = set()
    for fs in sets_:
        keys = lib.keys_of(fs, kind)
        if not keys:
            continue
        sp = fs.locus.organism
        species.add(sp)
        for k in keys:
            seen[k].add(sp)
            classes[k].add(fs.locus.vclass)
    n = len(species)
    return (({k: len(v) / n for k, v in seen.items()} if n else {}), n,
            {k: sorted(v) for k, v in classes.items()})


def shared_roots(sets_a: list[lib.FlankSet], sets_b: list[lib.FlankSet],
                 background: dict[str, float], kind: str = "root",
                 min_frac: float = PARALOGON_MIN_FRAC) -> list[dict]:
    """Root families carried by both neighbourhoods in most species, with
    the random-window background they have to beat."""
    pa, na, ca = side_prevalence(sets_a, kind)
    pb, nb, cb = side_prevalence(sets_b, kind)
    rows = []
    for k in set(pa) & set(pb):
        if pa[k] < min_frac or pb[k] < min_frac:
            continue
        bg = background.get(k, 0.0)
        rows.append(dict(
            key=k, frac_a=pa[k], frac_b=pb[k], n_species_a=na, n_species_b=nb,
            n_vclass_a=len(ca.get(k, [])), n_vclass_b=len(cb.get(k, [])),
            vclass_a=";".join(ca.get(k, [])), vclass_b=";".join(cb.get(k, [])),
            background=bg,
            enrichment=(min(pa[k], pb[k]) / bg) if bg else float("inf")))
    # the key is the final tiebreak, not the dict's iteration order: set
    # iteration is hash-seeded per process, so without it two runs of the
    # same data write different files
    rows.sort(key=lambda r: (-r["enrichment"], -min(r["frac_a"], r["frac_b"]),
                             r["key"]))
    return rows


def paralogon_summary(name: str, rows: list[dict], bg_n: int) -> dict:
    """Counted at three prevalence bars and two background bars, so a reader
    can see whether the answer moves when the line does."""
    def n_at(bar):
        return sum(1 for r in rows
                   if min(r["frac_a"], r["frac_b"]) >= bar)
    return dict(
        pair=name, n_shared=len(rows),
        n_at_10=n_at(0.10), n_at_25=n_at(0.25), n_at_50=n_at(0.50),
        n_shared_bg_lt_05=sum(1 for r in rows if r["background"] < 0.05),
        n_shared_bg_lt_01=sum(1 for r in rows if r["background"] < 0.01),
        n_control_windows=bg_n,
        top="; ".join(f"{r['key']}(bg {r['background']:.3f}, "
                      f"{min(r['frac_a'], r['frac_b']):.2f})" for r in rows[:8]))


# ------------------------------------------------- consensus-based calling

class ConsensusCaller:
    """Per-paralog flank consensus, with leave-one-genome-out support.

    Counts are held per (class, key) as the set of species carrying the key
    and per class as the set of species present, so removing one genome is a
    set difference rather than a rebuild.
    """

    def __init__(self, sets_by_class: dict[str, list[lib.FlankSet]], kind: str):
        self.kind = kind
        self.species: dict[str, set] = {}
        self.carriers: dict[str, dict[str, set]] = {}
        self.acc_species: dict[str, str] = {}
        for cls, sets_ in sets_by_class.items():
            sp, car = set(), defaultdict(set)
            for fs in sets_:
                keys = lib.keys_of(fs, kind)
                if not keys:
                    continue
                s = fs.locus.organism
                self.acc_species[fs.locus.accession] = s
                sp.add(s)
                for k in keys:
                    car[k].add(s)
            self.species[cls] = sp
            self.carriers[cls] = dict(car)

    def consensus(self, cls: str, drop_species: str | None = None,
                  frac: float | None = None) -> frozenset:
        # resolved at call time, not bound as a default: frac_sweep() varies
        # the module constant and a default argument would ignore it
        frac = CONSENSUS_FRAC if frac is None else frac
        sp = self.species[cls] - ({drop_species} if drop_species else set())
        n = len(sp)
        if not n:
            return frozenset()
        out = {k for k, car in self.carriers[cls].items()
               if len(car - ({drop_species} if drop_species else set())) / n >= frac}
        return frozenset(out)

    def score(self, fs: lib.FlankSet, classes: tuple[str, ...],
              leave_out: bool = True) -> dict:
        """Overlap of this locus's keys with each class consensus."""
        keys = lib.keys_of(fs, self.kind)
        drop = fs.locus.organism if leave_out else None
        scores, sizes = {}, {}
        for cls in classes:
            cons = self.consensus(cls, drop_species=drop)
            sizes[cls] = len(cons)
            scores[cls] = len(keys & cons)
        ranked = sorted(classes, key=lambda c: (-scores[c], c))
        best, second = ranked[0], ranked[1] if len(ranked) > 1 else None
        margin = scores[best] - (scores[second] if second else 0)
        return dict(
            n_keys=len(keys), scores=scores, consensus_sizes=sizes,
            best=best if scores[best] > 0 else "",
            best_score=scores[best], margin=margin,
            call=(best if (scores[best] > 0 and margin > 0
                           and len(keys) >= MIN_KEYS) else "no_call"))


def calibrate(caller: ConsensusCaller, labelled: list[tuple[lib.FlankSet, str]],
              classes: tuple[str, ...]) -> tuple[list[dict], dict]:
    """Leave-one-genome-out accuracy on loci whose paralog their own
    assembly's annotation already establishes.

    The annotation is evidence this caller never sees: it is scored on the
    symbols of the *neighbouring* genes, and the truth is the symbol of the
    gene itself.
    """
    rows = []
    for fs, truth in labelled:
        s = caller.score(fs, classes, leave_out=True)
        rows.append(dict(
            label=fs.locus.label, organism=fs.locus.organism,
            vclass=fs.locus.vclass, cell=fs.locus.cell, truth=truth,
            call=s["call"], best=s["best"], best_score=s["best_score"],
            margin=s["margin"], n_keys=s["n_keys"],
            **{f"score_{c}": s["scores"][c] for c in classes}))
    called = [r for r in rows if r["call"] != "no_call"]
    correct = [r for r in called if r["call"] == r["truth"]]
    per_class = {}
    for c in classes:
        sub = [r for r in rows if r["truth"] == c]
        subc = [r for r in sub if r["call"] != "no_call"]
        per_class[c] = dict(
            n=len(sub), n_called=len(subc),
            n_correct=sum(1 for r in subc if r["call"] == c),
            accuracy=(sum(1 for r in subc if r["call"] == c) / len(subc))
            if subc else 0.0)
    summary = dict(
        n=len(rows), n_called=len(called),
        call_rate=len(called) / len(rows) if rows else 0.0,
        n_correct=len(correct),
        accuracy=len(correct) / len(called) if called else 0.0,
        per_class=per_class)
    return rows, summary


def top_consensus(caller: ConsensusCaller, cls: str, limit: int = 25):
    """The consensus keys of one class, ranked by how many species carry them."""
    n = len(caller.species[cls]) or 1
    rows = [(k, len(car), len(car) / n)
            for k, car in caller.carriers[cls].items()]
    rows.sort(key=lambda r: (-r[1], r[0]))
    return rows[:limit], n


def counter_of(sets_: list[lib.FlankSet], kind: str) -> Counter:
    c: Counter = Counter()
    for fs in sets_:
        c.update(lib.keys_of(fs, kind))
    return c


# ------------------------------------------------------- caller controls

def caller_null(caller: "ConsensusCaller", controls: dict,
                classes: tuple[str, ...]) -> tuple[list[dict], dict]:
    """The caller run on the random control windows.

    A consensus overlap of 1 or 2 keys means nothing until it is known how
    often a random neighbourhood in the same genomes reaches it. This is
    what turns the cyclostome scores from a weak answer into no answer.
    """
    rows = []
    for acc, sets_ in sorted(controls.items()):
        for fs in sets_:
            s = caller.score(fs, classes, leave_out=False)
            rows.append(dict(
                label=f"{fs.locus.organism.replace(' ', '_')}|{acc}|"
                      f"{fs.locus.cell}",
                organism=fs.locus.organism, vclass=fs.locus.vclass,
                accession=acc, anchor_gene=fs.locus.annot_gene,
                n_keys=s["n_keys"], call=s["call"], best=s["best"],
                best_score=s["best_score"], margin=s["margin"],
                **{f"score_{c}": s["scores"][c] for c in classes}))
    n = len(rows)
    called = [r for r in rows if r["call"] != "no_call"]
    dist: dict[int, int] = {}
    for r in rows:
        dist[r["best_score"]] = dist.get(r["best_score"], 0) + 1
    summary = dict(
        n=n, n_called=len(called),
        false_call_rate=len(called) / n if n else 0.0,
        max_best_score=max((r["best_score"] for r in rows), default=0),
        frac_best_ge_1=sum(1 for r in rows if r["best_score"] >= 1) / n if n else 0.0,
        frac_best_ge_2=sum(1 for r in rows if r["best_score"] >= 2) / n if n else 0.0,
        frac_best_ge_3=sum(1 for r in rows if r["best_score"] >= 3) / n if n else 0.0,
        score_distribution={str(k): dist[k] for k in sorted(dist)})
    return rows, summary


def select_frac(sweep: list[dict]) -> tuple[float, str]:
    """The consensus threshold, chosen from the sweep rather than typed.

    A loose threshold admits more keys into each consensus, which buys calls
    on loci whose neighbourhood is only partly conserved -- and admits
    calls on neighbourhoods that are not ITPR neighbourhoods at all. Call
    rate alone therefore selects the loosest setting on offer, which is how
    an instrument gets tuned into agreeing with itself.

    So the sweep measures both, on data the choice cannot borrow from: call
    rate on the annotation-confirmed loci, and the false-call rate on the
    random control windows. The rule is to maximise their difference -- one
    call gained is worth exactly one false call avoided -- with ties broken
    toward the larger frac, which is the stricter consensus.

    Accuracy on the labelled loci is *not* what is optimised (it is 1.000
    across the whole sweep, so it separates nothing); it is reported.
    """
    if not sweep:
        return CONSENSUS_FRAC, "empty sweep; kept the declared default"
    scored = [(r["call_rate"] - r.get("false_call_rate", 0.0),
               r["consensus_frac"], r) for r in sweep]
    best = max(s for s, _, _ in scored)
    pick = max(f for s, f, _ in scored if s >= best - 1e-12)
    row = [r for s, f, r in scored if f == pick][0]
    return pick, (f"maximises call rate - random-window false-call rate "
                  f"({row['call_rate']:.3f} - "
                  f"{row.get('false_call_rate', 0.0):.3f} = {best:.3f}); "
                  f"accuracy on labelled loci {row['accuracy']:.3f}")


def frac_sweep(sets_by_class: dict, labelled: list, classes: tuple[str, ...],
               controls: dict | None = None,
               fracs=(0.20, 0.30, 0.40, 0.50, 0.60, 0.70)) -> list[dict]:
    """What the consensus threshold buys and what it costs, measured.

    Both sides at every setting: the call rate on the annotation-confirmed
    loci, and the rate at which random neighbourhoods in the same genomes
    are called an ITPR paralog. Committing the sweep means the chosen value
    is visible as a choice with both its numbers attached.
    """
    global CONSENSUS_FRAC
    keep = CONSENSUS_FRAC
    out = []
    try:
        for f in fracs:
            CONSENSUS_FRAC = f
            caller = ConsensusCaller(sets_by_class, "relaxed")
            _, summ = calibrate(caller, labelled, classes)
            nullsumm = ({} if controls is None
                        else caller_null(caller, controls, classes)[1])
            row = dict(consensus_frac=f, n=summ["n"], n_called=summ["n_called"],
                       call_rate=summ["call_rate"], accuracy=summ["accuracy"],
                       n_control=nullsumm.get("n", 0),
                       n_control_called=nullsumm.get("n_called", 0),
                       false_call_rate=nullsumm.get("false_call_rate", 0.0),
                       margin=summ["call_rate"]
                       - nullsumm.get("false_call_rate", 0.0))
            for c in classes:
                row[f"consensus_size_{c}"] = len(caller.consensus(c))
                row[f"call_rate_{c}"] = (
                    summ["per_class"][c]["n_called"] / summ["per_class"][c]["n"]
                    if summ["per_class"][c]["n"] else 0.0)
                row[f"accuracy_{c}"] = summ["per_class"][c]["accuracy"]
            out.append(row)
    finally:
        CONSENSUS_FRAC = keep
    return out


def null_tail(null_summary: dict) -> tuple[dict[int, float], int]:
    """P(a random neighbourhood reaches at least this consensus overlap).

    Built from the committed null score distribution, so a call's own tail
    probability travels with it in the table rather than living in a
    sentence of the report.
    """
    dist = {int(k): v for k, v in null_summary.get("score_distribution", {}).items()}
    n = sum(dist.values())
    tail, run = {}, 0
    for s in sorted(dist, reverse=True):
        run += dist[s]
        tail[s] = run / n if n else 0.0
    return tail, null_summary.get("max_best_score", 0)


def null_verdict(best_score: int, call: str, max_null: int) -> str:
    """A call is only evidence if it beats every random neighbourhood.

    The bar is the null's own maximum rather than a probability cut, because
    with 726 control windows a 1-in-726 tail is not a rate anyone should
    build a claim on.
    """
    if call == "no_call":
        return "no_call"
    return "supported" if best_score > max_null else "within_null"
