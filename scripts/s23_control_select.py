"""s23_control_select.py — which profile controls which clade, decided by data.

`s23_control_profiles` says what the candidate control profiles are and why
each is a candidate. This module measures them and picks one per control
clade.

**The rule.** For every candidate profile and every control clade, count how
many of that clade's swept reference proteomes carry a hit. The winner is the
profile with the highest clade coverage; ties go to the one whose proteins are
longer, because a control on a longer multi-exon gene exercises more of the
spliced aligner. The whole table is written out, not just the winner, so the
choice can be read as a measurement rather than taken on trust — and so a
clade where *every* candidate is thin is visible as such instead of being
represented by the least-bad one without comment.

**Why the winner is not simply the broadest profile.** PF00004 (AAA) is in
every proteome of every clade and would win everywhere on coverage alone. It
is kept last in `CANDIDATES` and used as the fallback it is: a control that
proves the assembly is searchable and little else. The ranking therefore
prefers, among profiles that clear a coverage bar, the one carrying the
longest proteins — which is the one whose recovery says the most about
whether a 2,700 aa receptor could have been recovered.

**What this does not decide.** Family membership. Every control bait still
goes through `s5_bait_screen` and must be assigned to *neither* ITPR nor RyR;
a control protein the profiles call a family member is a finding.
"""

from __future__ import annotations

import statistics
import sys
from collections import defaultdict
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import s23_bait_spec as spec                               # noqa: E402
import s23_control_profiles as prof                        # noqa: E402
import s23_controls as ctl                                 # noqa: E402
import s23_scope as scope                                  # noqa: E402

OUT_DIR = PROJECT_ROOT / "results" / "s23_baits"

#: A control bait's length band. Wider than `spec.MIR_BAND_AA` because the
#: candidates are different protein families with different natural sizes; the
#: floor is the same principle (a bait shorter than the domain it is supposed
#: to carry is a fragment) and the ceiling keeps a control from being the size
#: of the receptor itself.
CONTROL_BAND_AA = (200, 3_000)

#: How much of the candidate profile's own model a control bait must carry.
#: Measured from the domtblout's model coordinates, not from the protein's
#: length: a 900 aa protein matching 15 % of the model is a fragment of the
#: control family however long it is.
MIN_MODEL_COV = 0.60

#: The clade-coverage bar a profile has to clear before length is allowed to
#: break the tie. Below it, coverage alone decides — a profile present in 90 %
#: of the clade beats a longer one present in 20 %, whatever its proteins
#: weigh.
GOOD_COVERAGE = 0.75

#: `control_strength` — a strong control is present across a quarter of the
#: clade's swept proteomes, carries most of its model, and is a substantial
#: multi-exon protein rather than a lone domain.
STRONG_MIN_CLADE_FRAC = spec.CONTROL_STRONG_MIN_CLADE_FRAC
STRONG_MIN_AA = spec.CONTROL_STRONG_MIN_AA
STRONG_MIN_MODEL_COV = 0.80


# ------------------------------------------------------------------ measuring

def _usable(row: dict) -> bool:
    lo, hi = CONTROL_BAND_AA
    return (lo <= row["length"] <= hi
            and row.get("model_cov", 0.0) >= MIN_MODEL_COV)


def hits_for(pfam: str, threads: int, force: bool) -> dict[str, dict[str, dict]]:
    """group -> accession -> best hit row, for one candidate profile."""
    return {g: ctl.parse_hits(
                ctl.run_search(g, threads, force, pfam=pfam,
                               hmm=prof.profile_path(pfam)), g)
            for g in ctl.GROUP_DBS}


def measure(presence: list[dict], taxonomy: dict[str, dict],
            threads: int = 6, force: bool = False,
            candidates: tuple | None = None
            ) -> tuple[dict[str, dict], dict[str, dict]]:
    """(profile -> (rank, clade) -> stats, profile -> group -> hits)."""
    candidates = candidates if candidates is not None else prof.primary()
    swept: dict[tuple, set] = defaultdict(set)
    for r in presence:
        tax = taxonomy.get(r["taxid"], {})
        for rank in ("class", "phylum"):
            if tax.get(rank):
                swept[(rank, tax[rank])].add(r["taxid"])

    stats: dict[str, dict] = {}
    all_hits: dict[str, dict] = {}
    for cand in candidates:
        pfam = cand["pfam"]
        hits = hits_for(pfam, threads, force)
        all_hits[pfam] = hits
        taxa_hit: dict[tuple, set] = defaultdict(set)
        lens: dict[tuple, list] = defaultdict(list)
        for group_hits in hits.values():
            for row in group_hits.values():
                if not _usable(row):
                    continue
                tax = taxonomy.get(row["taxid"], {})
                for rank in ("class", "phylum"):
                    key = (rank, tax.get(rank) or "")
                    if key[1]:
                        taxa_hit[key].add(row["taxid"])
                        lens[key].append(row["length"])
        stats[pfam] = {
            key: {"swept": len(taxa), "with_hit": len(taxa_hit.get(key, ())),
                  "frac": (len(taxa_hit.get(key, ())) / len(taxa)
                           if taxa else 0.0),
                  "median_len": (int(statistics.median(lens[key]))
                                 if lens.get(key) else 0)}
            for key, taxa in swept.items()}
    return stats, all_hits


# ------------------------------------------------------------------ choosing

def _rank_key(st: dict) -> tuple:
    """Coverage first, then protein length once coverage is good enough."""
    good = st["frac"] >= GOOD_COVERAGE
    return (-int(good), -round(st["frac"], 3) if not good else 0,
            -st["median_len"], -st["frac"])


def choose(clades: list[dict], stats: dict[str, dict],
           candidates: tuple | None = None) -> list[dict]:
    """One profile per control clade, with every candidate's number beside it."""
    candidates = candidates if candidates is not None else prof.CANDIDATES
    out = []
    for c in clades:
        key = (c["rank"], c["clade"])
        scored = []
        for cand in candidates:
            st = stats.get(cand["pfam"], {}).get(key)
            if st and st["with_hit"]:
                scored.append((cand, st))
        if not scored:
            out.append(dict(c, pfam="", why="no candidate control profile has "
                                            "a usable hit in any swept "
                                            "proteome of this clade"))
            continue
        scored.sort(key=lambda t: _rank_key(t[1]))
        cand, st = scored[0]
        others = ", ".join(f"{c2['pfam']} {s2['frac']:.0%}"
                           for c2, s2 in scored[1:])
        out.append(dict(
            c, pfam=cand["pfam"], label=cand["label"],
            frac=round(st["frac"], 4), with_hit=st["with_hit"],
            swept_proteomes=st["swept"], median_len=st["median_len"],
            n_candidates=len(scored),
            why=(f"{cand['label']} ({cand['pfam']}) is in "
                 f"{st['with_hit']}/{st['swept']} of the clade's swept "
                 f"proteomes ({st['frac']:.0%}), median {st['median_len']} aa"
                 + (f"; runners-up {others}" if others else ""))))
    return out


def control_strength(length: int, model_cov: float, clade_frac: float
                     ) -> tuple[str, str]:
    """Is this control strong, and if not, what is missing from it?

    Graded rather than counted, for the reason S23a gave: an absence claim
    standing on a 233 aa control found in one proteome of 36 is a weaker claim
    than one standing on a 950 aa protein found across a whole class, and
    presenting both as "controlled" hides the difference the reader needs.
    """
    problems = []
    if clade_frac < STRONG_MIN_CLADE_FRAC:
        problems.append(f"present in only {clade_frac:.0%} of the clade's "
                        "swept proteomes")
    if length < STRONG_MIN_AA:
        problems.append(f"{length} aa — a lone domain rather than a "
                        f"substantial multi-exon gene (< {STRONG_MIN_AA})")
    if model_cov < STRONG_MIN_MODEL_COV:
        problems.append(f"carries {model_cov:.0%} of the control profile's "
                        "model")
    if not problems:
        return "strong", (f"{length} aa carrying {model_cov:.0%} of its "
                          f"profile, present across {clade_frac:.0%} of the "
                          "clade")
    return "weak", "; ".join(problems)


def select(chosen: list[dict], taxonomy: dict[str, dict],
           all_hits: dict[str, dict], stats: dict[str, dict]
           ) -> tuple[list[dict], list[dict]]:
    """One control bait per clade, drawn from that clade's chosen profile."""
    baits, unfilled = [], []
    for c in chosen:
        pfam = c.get("pfam")
        if not pfam:
            unfilled.append(dict(c, why=c["why"]))
            continue
        key = (c["rank"], c["clade"])
        pool = []
        for group, hits in all_hits[pfam].items():
            for row in hits.values():
                if not _usable(row):
                    continue
                tax = taxonomy.get(row["taxid"])
                if not tax or (tax.get(c["rank"]) or "") != c["clade"]:
                    continue
                pool.append(dict(row, rank=c["rank"], clade=c["clade"],
                                 lineage_kingdom=tax.get("kingdom", ""),
                                 lineage_phylum=tax.get("phylum", "")))
        if not pool:
            unfilled.append(dict(c, why=f"{pfam} covers the clade but no hit "
                                        "sits in the control bait band"))
            continue
        # Longest first: within one profile and one clade, the longest
        # full-model protein is the hardest target and therefore the most
        # informative control. Score breaks ties.
        pool.sort(key=lambda r: (-r["length"], -r["score"], r["accession"]))
        pick = pool[0]
        frac = c.get("frac", 0.0)
        strength, why = control_strength(pick["length"], pick["model_cov"],
                                         frac)
        baits.append(dict(
            pick, pfam=pfam, profile_label=c.get("label", ""),
            control_for=c["clade"], control_rank=c["rank"],
            swept=c["swept"], candidates=len(pool),
            clade_frac=round(frac, 4), strength=strength, strength_note=why,
            reason=(f"{c.get('label', pfam)} control for {c['clade']}: "
                    f"{pick['length']} aa, {pick['score']:.0f} bits, "
                    f"{pick['model_cov']:.0%} of the model; "
                    f"{c['why']}; the clade S20 swept {c['swept']} proteomes "
                    f"of for 0 ITPR; control is {strength} — {why}")))
    return baits, unfilled


# ------------------------------------------------------------------ driver

def build(threads: int = 6, force: bool = False):
    """(control bait rows, unfilled clades, sequences, coverage stats, choice)."""
    presence = scope.load_presence()
    taxonomy = scope.load_taxonomy()
    clades = scope.control_clades(presence, taxonomy)
    searched = prof.primary()
    stats, all_hits = measure(presence, taxonomy, threads, force, searched)
    chosen = choose(clades, stats, searched)
    # The fallback is a fallback: searched only if the ranked candidates leave
    # a control clade with nothing at all.
    if any(not c.get("pfam") for c in chosen) and prof.fallbacks():
        searched = searched + prof.fallbacks()
        stats, all_hits = measure(presence, taxonomy, threads, force, searched)
        chosen = choose(clades, stats, searched)
    baits, unfilled = select(chosen, taxonomy, all_hits, stats)
    seqs = ctl.fetch_sequences(baits)
    missing = [r["accession"] for r in baits if r["accession"] not in seqs]
    if missing:
        raise SystemExit(
            f"{len(missing)} control baits have no sequence in the archived "
            f"proteome DBs: {missing[:5]}\n"
            "  the domtblout and the FASTA are out of step — rerun with "
            "--force")
    return baits, unfilled, seqs, stats, chosen


def write_tables(stats: dict, chosen: list[dict],
                 clades: list[dict] | None = None) -> None:
    """Every candidate x clade, not only the winners.

    A candidate that was never searched is written as `not_searched` rather
    than left out: the fallback is only consulted when the ranked candidates
    fail, and a blank row would be indistinguishable from a profile that was
    searched and found nothing.
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    keys = {(c["rank"], c["clade"]) for c in chosen}
    with open(OUT_DIR / "control_profile_coverage.tsv", "w") as fh:
        fh.write("rank\tclade\tpfam\tlabel\tsearched\tswept_proteomes\t"
                 "with_hit\tfrac\tmedian_len\tchosen\n")
        pick = {(c["rank"], c["clade"]): c.get("pfam", "") for c in chosen}
        for rank, clade in sorted(keys):
            for cand in prof.CANDIDATES:
                searched = cand["pfam"] in stats
                st = stats.get(cand["pfam"], {}).get((rank, clade))
                if searched and not st:
                    continue
                fh.write(f"{rank}\t{clade}\t{cand['pfam']}\t{cand['label']}\t"
                         f"{int(searched)}\t"
                         + (f"{st['swept']}\t{st['with_hit']}\t"
                            f"{st['frac']:.4f}\t{st['median_len']}\t"
                            if st else "\t\t\t\t")
                         + f"{int(pick.get((rank, clade)) == cand['pfam'])}\n")
    with open(OUT_DIR / "control_profile_choice.tsv", "w") as fh:
        fh.write("rank\tclade\tswept\tpfam\tlabel\tfrac\twith_hit\t"
                 "median_len\tn_candidates\twhy\n")
        for c in chosen:
            fh.write(f"{c['rank']}\t{c['clade']}\t{c.get('swept', '')}\t"
                     f"{c.get('pfam', '')}\t{c.get('label', '')}\t"
                     f"{c.get('frac', 0):.4f}\t{c.get('with_hit', 0)}\t"
                     f"{c.get('median_len', 0)}\t{c.get('n_candidates', 0)}\t"
                     f"{c['why']}\n")


def build_and_write(threads: int = 6, force: bool = False):
    baits, unfilled, seqs, stats, chosen = build(threads, force)
    write_tables(stats, chosen)
    n_strong = sum(1 for b in baits if b["strength"] == "strong")
    print(f"{len(baits)} control baits over {len(chosen)} clades "
          f"({n_strong} strong / {len(baits) - n_strong} weak), "
          f"{len(unfilled)} unfilled")
    from collections import Counter
    for pfam, n in Counter(b["pfam"] for b in baits).most_common():
        label = next(c["label"] for c in prof.CANDIDATES if c["pfam"] == pfam)
        print(f"  {pfam} {label:14s} {n:3d} clade(s)")
    for u in unfilled:
        print(f"  ! {u['clade']:22s} {u['why']}")
    return baits, unfilled, seqs, stats, chosen


if __name__ == "__main__":
    build_and_write()
