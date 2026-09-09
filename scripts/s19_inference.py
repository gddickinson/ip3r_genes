"""S19 — where the *inference* methods ran out, not the search.

Four limits this project measured while doing something else. Each of them is
a general result about the method rather than about the IP3 receptors, and
each is recomputed here from the committed tables of the session that met it,
so the methods report cannot drift from those sessions' own reports.

* **A reconciliation's loss count is mostly sampling** (S13). Losses implied
  by a gene tree over a 134-tip representative alignment, checked cell by
  cell against a 309-genome ledger: a "loss" whose gene is present and intact
  in the assembly is a fact about who was sampled.
* **Synteny disambiguation is accurate and unavailable** (S15a). S8's caller
  is right wherever it acts, and it almost never acts on the cells that need
  it, because a fragment's contig carries no neighbours. Accuracy and
  reach are different numbers and only reporting both is honest.
* **Per-site codon models are past their power at these divergences** (S9,
  S17). The synonymous rate is unidentifiable at a share of sites, so any
  per-site or per-element omega built as a ratio of sums collapses.
* **A likelihood tree does not always resolve the question asked of it**
  (S7). The AU test's 95 % confidence set is reported as a count, because
  "the sister is unresolved" and "the sister is X" are different claims and
  the table says which one the data supports.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402

RECON = S.RESULTS / "reconciliation"
LOSS = S.RESULTS / "loss_dynamics"
CONSTRAINT = S.RESULTS / "constraint"
SELECTION = S.RESULTS / "selection"
PHYLO = S.RESULTS / "phylogeny"

#: HyPhy FEL reports alpha at a hard upper bound when a site carries no
#: usable synonymous signal. S9 measured median dS well above 1 for these
#: paralogs, so an alpha at or above 100 is equally uninformative whether or
#: not it is pinned.
ALPHA_CEILING = 100.0
ALPHA_HARD_BOUND = 9_999.0

#: The verdicts in S13's audit that mean "the gene is in the genome" — an
#: implied loss that meets one of these is an artefact of representative
#: sampling, not an event.
ARTEFACT_VERDICTS = ("sampling_artefact", "in_genome_unsampled",
                     "present_unsampled")
REAL_VERDICTS = ("corroborated_loss", "loss_with_remnant")


def reconciliation_losses(log=S.log) -> list[list]:
    path = RECON / "loss_verdicts.tsv"
    if not path.exists():
        return []
    # S13 committed its own `all` aggregate row in this table. Summing every
    # row would count each cell twice, so the aggregate is dropped here and
    # recomputed — and then checked against S13's, because two aggregates
    # that disagree are a bug in exactly one of them.
    by_paralog: dict[str, dict[str, int]] = {}
    committed_all: dict[str, int] = {}
    for r in S.read_tsv(path):
        n = int(S.fnum(r["n_species"], float, 0))
        if r["paralog"].lower() == "all":
            committed_all[r["verdict"]] = n
            continue
        by_paralog.setdefault(r["paralog"], {})[r["verdict"]] = n
    # `sampled_in_gene_tree` is not an implied loss — it is the paralog
    # being present at a tip — so it never enters the denominator.
    out = []
    seen_verdicts: set[str] = set()
    for paralog, counts in sorted(by_paralog.items()):
        counts = {k: v for k, v in counts.items()
                  if k != "sampled_in_gene_tree"}
        seen_verdicts |= set(counts)
        implied = sum(counts.values())
        artefact = sum(v for k, v in counts.items() if k in ARTEFACT_VERDICTS)
        real = sum(v for k, v in counts.items() if k in REAL_VERDICTS)
        other = implied - artefact - real
        out.append([paralog, implied, artefact, real, other,
                    round(artefact / implied, 4) if implied else "",
                    ";".join(f"{k}={v}" for k, v in sorted(counts.items()))])
    tot = [sum(r[i] for r in out) for i in (1, 2, 3, 4)]
    check = sum(v for k, v in committed_all.items()
                if k != "sampled_in_gene_tree")
    out.append(["ALL", *tot, round(tot[1] / tot[0], 4) if tot[0] else "",
                f"recomputed from the per-paralog rows; S13's own `all` row "
                f"totals {check} — {'agrees' if check == tot[0] else 'DISAGREES'}"])
    S.write_tsv(S.out_dir() / "reconciliation_loss_audit.tsv",
                ["paralog", "n_implied_loss_cells", "n_sampling_artefact",
                 "n_corroborated", "n_other_verdict", "frac_artefact",
                 "verdict_breakdown"], out)
    log(f"reconciliation_loss_audit: {tot[1]}/{tot[0]} implied-loss cells are "
        f"sampling artefacts, {tot[2]} corroborated")
    return out


def synteny_power(log=S.log) -> list[list]:
    """Accuracy where it acts, against how often it acts. Both, never one."""
    rows = []
    for r in S.read_tsv(LOSS / "synteny_reach_summary.tsv"):
        n = S.fnum(r["n_regions"], float, 0)
        rows.append(["reach", r["window"], int(n),
                     int(S.fnum(r["n_reached"], float, 0)),
                     round(S.fnum(r["n_reached"], float, 0) / n, 4) if n else "",
                     f"no neighbourhood {r['n_no_neighbourhood']}; too few "
                     f"keys {r['n_too_few_keys']}; no gene table "
                     f"{r['n_no_gene_table']}"])
    for r in S.read_tsv(LOSS / "caller_accuracy_by_keys.tsv"):
        n = S.fnum(r["n"], float, 0)
        # A key bin the caller never acted in has no accuracy, and S15a
        # writes `nan` there. Carrying that through as a number would put a
        # NaN in a report; the honest cell is empty with the call rate beside
        # it saying why.
        called = S.fnum(r["n_called"], float, 0)
        rows.append(["accuracy", f"{r['keys_lo']}-{r['keys_hi']} keys", int(n),
                     int(S.fnum(r["n_correct"], float, 0)),
                     S.fnum(r["accuracy"], float, 0.0) if called else "",
                     f"call rate {round(S.fnum(r['call_rate'], float, 0), 4)}"
                     + ("" if called else " — the caller never acted here")])
    S.write_tsv(S.out_dir() / "synteny_power.tsv",
                ["measure", "scope", "n", "n_positive", "frac", "detail"],
                rows)
    log(f"synteny_power: {len(rows)} rows")
    return rows


def codon_model_power(log=S.log) -> list[list]:
    """Per-site: how often the synonymous rate is identifiable at all."""
    path = CONSTRAINT / "fel_sites.tsv"
    if not path.exists():
        return []
    sites = [r for r in S.read_tsv(path) if (r.get("resi") or "").strip()]
    rows = []
    for paralog in list(S.PARALOGS) + ["ALL"]:
        sub = (sites if paralog == "ALL"
               else [r for r in sites if r["paralog"] == paralog])
        if not sub:
            continue
        alphas = [S.fnum(r.get("alpha"), float, 0.0) for r in sub]
        zero = sum(1 for a in alphas if a <= 0)
        high = sum(1 for a in alphas if a >= ALPHA_CEILING)
        pinned = sum(1 for a in alphas if a >= ALPHA_HARD_BOUND)
        usable = [a for a in alphas if 0 < a < ALPHA_CEILING]
        rows.append([paralog, len(sub), zero, high, pinned, zero + high,
                     round((zero + high) / len(sub), 4),
                     round(S.median(usable), 3) if usable else ""])
    S.write_tsv(S.out_dir() / "codon_model_power.tsv",
                ["paralog", "n_sites", "n_alpha_zero", "n_alpha_ge_100",
                 "n_alpha_at_hyphy_bound", "n_unidentifiable",
                 "frac_unidentifiable", "median_usable_alpha"], rows)
    tot = rows[-1]
    log(f"codon_model_power: alpha unidentifiable at {tot[5]}/{tot[1]} sites "
        f"({tot[6]:.1%})")
    return rows


def saturation(log=S.log) -> list[list]:
    """Pairwise dS at these divergences, and how much of it is saturated.

    The gene-scale companion to the per-site table: a pair whose synonymous
    sites are saturated carries no usable rate however the model is fitted,
    and S9 flagged them in its own committed table rather than leaving the
    reader to infer it from a large dS.
    """
    path = SELECTION / "pairwise_dnds.tsv"
    if not path.exists():
        return []
    rows = []
    by_set: dict[str, list[dict]] = {}
    for r in S.read_tsv(path):
        by_set.setdefault(r["set"], []).append(r)
    for name, sub in sorted(by_set.items()) + [("ALL", S.read_tsv(path))]:
        sat = sum(1 for r in sub
                  if str(r.get("saturated")).lower() in ("1", "true", "yes"))
        ds = [S.fnum(r.get("dS"), float, 0.0) for r in sub]
        rows.append([name, len(sub), sat, round(sat / max(1, len(sub)), 4),
                     round(S.median(ds), 3),
                     round(max(ds), 3) if ds else ""])
    S.write_tsv(S.out_dir() / "saturation.tsv",
                ["set", "n_pairs", "n_saturated", "frac_saturated",
                 "median_dS", "max_dS"], rows)
    log(f"saturation: {rows[-1][2]}/{rows[-1][1]} pairs flagged saturated")
    return rows


def tree_resolution(log=S.log) -> list[list]:
    """How many topologies the AU test leaves in the 95 % confidence set."""
    path = PHYLO / "au_test.tsv"
    if not path.exists():
        return []
    rows = []
    kept = 0
    for r in S.read_tsv(path):
        p_au = S.fnum(r.get("p_AU"), float, 0.0)
        in_set = p_au >= 0.05
        kept += int(in_set)
        rows.append([r.get("hypothesis", r.get("tree", "")),
                     S.fnum(r.get("logL"), float, 0.0),
                     S.fnum(r.get("deltaL"), float, 0.0), p_au, int(in_set)])
    rows.append(["95 % confidence set size", "", "", "", kept])
    S.write_tsv(S.out_dir() / "tree_resolution.tsv",
                ["hypothesis", "logL", "delta_logL", "p_AU",
                 "in_95pct_set"], rows)
    log(f"tree_resolution: {kept} of {len(rows) - 1} topologies in the 95 % "
        f"confidence set")
    return rows


def run(log=S.log) -> dict:
    recon = reconciliation_losses(log)
    syn = synteny_power(log)
    fel = codon_model_power(log)
    sat = saturation(log)
    tree = tree_resolution(log)
    summary = {
        "reconciliation": ({"implied": recon[-1][1], "artefact": recon[-1][2],
                            "corroborated": recon[-1][3],
                            "other": recon[-1][4]} if recon else {}),
        "synteny_rows": len(syn),
        "fel": ({"n_sites": fel[-1][1], "n_unidentifiable": fel[-1][5],
                 "frac": fel[-1][6]} if fel else {}),
        "saturation": ({"n_pairs": sat[-1][1], "n_saturated": sat[-1][2],
                        "frac": sat[-1][3]} if sat else {}),
        "tree_confidence_set": tree[-1][4] if tree else None,
    }
    S.write_json(S.out_dir() / "inference_summary.json", summary)
    return summary


if __name__ == "__main__":
    run()
