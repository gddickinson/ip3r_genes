"""s23_calibration.py — S23's measured thresholds, read back rather than typed.

The discipline D13 applies to reports, applied to constants: every number here
is read out of `results/s23_baits/span_calibration.json`, which
`s23_span_calibration.py` measured from NCBI's own gene annotations. Nothing
in the sweep may hard-code a span or a max-intron, because both are results
and a retyped result silently stops tracking the table it came from.

Two thresholds, and the reason each one is a measurement:

  `max_intron_for(group)`   miniprot's `-G`. Too small does not lose a gene,
                            it **splits** one, and a split ITPR reads out of
                            the copy-number ledger as `fragment` — worse here
                            than in S5, because in this task the count *is*
                            the result. Floored at miniprot's own default so
                            an unmeasured group can never get a smaller `-G`
                            than the tool ships with.

  `contiguity_bar(group)`   D4's bar. An assembly whose contigs are shorter
                            than the gene cannot carry it, so its silence is
                            not evidence. **Per group**, because the family's
                            span varies ~100x outside the vertebrates against
                            6.5x inside: metazoan genes have a median span of
                            83 kb and protist and fungal ones 7-9 kb, and one
                            bar across that range is either far too strict for
                            the protists or far too lax for the metazoans.

The per-group split is also why S5b's headline caveat mostly does not apply
here. There, 120 of 309 genomes could not hold the gene on one contig and
recovery fell from 98-99 % above the bar to 57-70 % below it. The genomes
carrying this task's negative claims — land plants, Dikarya, Apicomplexa — are
measured against a ~9 kb bar, which almost any modern assembly clears.
"""

from __future__ import annotations

import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CALIBRATION = PROJECT_ROOT / "results" / "s23_baits" / "span_calibration.json"


def _load() -> dict:
    if not CALIBRATION.exists():
        raise SystemExit(
            f"S23 span calibration missing: {CALIBRATION}\n"
            "  run: python3 scripts/s23_span_calibration.py\n"
            "  the sweep will not run on typed thresholds — both of these are "
            "measurements, and one of them (-G) decides the ledger's call.")
    return json.loads(CALIBRATION.read_text())


def calibration() -> dict:
    return _load()


def max_intron_for(group: str) -> int:
    """miniprot `-G` for one kingdom-level group."""
    cal = _load()
    g = cal.get("by_group", {}).get(group)
    return int(g["max_intron_bp"] if g else cal["derived_max_intron_bp"])


def contiguity_bar(group: str) -> tuple[int, str]:
    """(bar in bp, how it was derived) for one group.

    The fallback is stated in the return value rather than hidden, because a
    group with no measured span is a group whose bar is borrowed — S23a
    measured no Viridiplantae gene span at all (the chlorophyte records carry
    locus tags NCBI has no gene record for), so every plant genome is judged
    against the global median and the report has to say so.
    """
    cal = _load()
    g = cal.get("by_group", {}).get(group)
    if g:
        return int(g["contiguity_bar_bp"]), (
            f"median span of {g['n']} measured {group} gene(s)")
    return int(cal["contiguity_bar_bp"]), (
        f"no {group} gene span was measured; the global median of "
        f"{cal['n_spans']} genes is used instead")


def spans_a_gene(contig_n50: int, group: str) -> tuple[bool, int, str]:
    """D4 — can this assembly's typical contig hold this group's gene?"""
    bar, why = contiguity_bar(group)
    return int(contig_n50 or 0) >= bar, bar, why


def measured_groups() -> list[str]:
    return sorted(_load().get("by_group", {}))


def summary() -> str:
    cal = _load()
    lines = [f"S23 span calibration: {cal['n_spans']} genes, "
             f"{cal['n_species']} species, {cal['n_bands']} bands",
             f"  global median span {cal['median_span_bp']:,} bp "
             f"(S5 vertebrate bar {cal['s5_vertebrate_bar_bp']:,} bp)"]
    for g, s in cal.get("by_group", {}).items():
        lines.append(f"  {g:14s} n={s['n']:2d}  bar {s['contiguity_bar_bp']:>9,} bp"
                     f"   -G {s['max_intron_bp']:>9,} bp")
    return "\n".join(lines)


if __name__ == "__main__":
    print(summary())
