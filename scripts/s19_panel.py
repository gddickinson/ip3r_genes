"""S19 — how much of the sweep's recall depended on the bait panel.

The S5 sweep ran miniprot with a 38-protein panel. miniprot aligns each query
independently, so every genome's retained `miniprot.gff` holds a **separable**
alignment record per bait: dropping baits from the file and re-clustering
reproduces exactly what the sweep would have reported had those baits never
been in the panel. Nothing is re-downloaded and nothing is re-aligned — the
ablation is exact, not a model of one, and `panel_validation.tsv` scores the
full-panel reproduction cell for cell against the committed ledger before any
ablation is read.

Two experiments.

* **Named panels.** Drop a clade band, drop a paralog's baits, drop the
  additions S5a made to the S3 seed set, keep one bait per paralog, keep only
  human. And the family-specific one: drop the three **unlabelled**
  `vertebrate_basal` baits. Those exist because no labelled ITPR1/ITPR2
  record exists for cyclostomes or (for two paralogs) chondrichthyans, and
  S5a found they were competing as a fourth paralog until rescue attribution
  was fixed. What they are worth has never been measured.
* **Single-bait recall against bait-target identity.** Every (genome, cell,
  bait) triple the sweep produced is a measurement of what one bait alone
  recovers at a known sequence distance. Pooled, that is the design rule: how
  close a panel must reach into a clade before a search over it is
  trustworthy.

The simulated call is the coverage half of `s5_classify.classify_class` —
"some locus this cell owns carries an alignment covering at least COV_FOUND
of its bait". Annotation and gap grading need the genome FASTA, which the
sweep deleted after searching; the found/not-found decision does not.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s5_classify as CL                                        # noqa: E402
import s5_sweep_lib as SW                                       # noqa: E402
import s19_lib as S                                             # noqa: E402

COV_FOUND = SW.COV_FOUND

#: Clade bands grouped for the taxon-drop ablations. The band vocabulary is
#: the panel's own (`s5_bait_spec.band_of`), so an ablation cannot be run on
#: a grouping the panel was not built under.
BAND_GROUP = {
    "mammalia": "mammal", "aves": "bird", "reptilia": "reptile",
    "amphibia": "amphibian", "actinopteri": "ray-finned fish",
    "chondrichthyes": "cartilaginous fish",
    "sarcopterygian_fish": "lobe-finned fish", "cyclostomata": "cyclostome",
}


def panels(meta: dict[str, dict]) -> dict[str, dict]:
    """name -> {baits, description}. Every panel is a subset of the 38."""
    ids = list(meta)

    def where(fn):
        return {b for b in ids if fn(meta[b])}

    itpr = where(lambda m: m["family"] == "ITPR")
    out: dict[str, tuple[set[str], str]] = {
        "full38": (set(ids), "the panel as run (S5b)"),
        "s3_seeds_only": (
            where(lambda m: m["origin"] != "s5_addition"),
            "the S3 seed set alone, before S5a's additions"),
        "labelled_only": (
            where(lambda m: m["clade"] != "vertebrate_basal"),
            "no unlabelled `vertebrate_basal` bait — the gar, chimaera and "
            "lamprey seeds that carry no paralog assignment"),
        "basal_only": (
            where(lambda m: m["clade"] == "vertebrate_basal"),
            "only the three unlabelled baits"),
        "no_ryr_control": (itpr, "the ITPR baits alone, no sister control"),
        "amniote_only": (
            where(lambda m: BAND_GROUP.get(m["band"]) in
                  ("mammal", "bird", "reptile")),
            "mammal + bird + reptile baits only — a tetrapod-centric panel"),
        "human_only": (
            where(lambda m: m["species"].startswith("Homo_sapiens")),
            "the human baits alone"),
    }
    for band, group in sorted(set(BAND_GROUP.items())):
        name = f"drop_{group.replace(' ', '_').replace('-', '_')}"
        if name in out:
            continue
        out[name] = (where(lambda m, b=band: m["band"] != b),
                     f"no {group} bait")
    for p in S.PARALOGS:
        out[f"drop_{p}_baits"] = (
            where(lambda m, p=p: m["clade"] != p),
            f"no bait labelled {p} — can the other paralogs' baits and the "
            f"unlabelled ones still place its gene?")
    # One bait per cell, deliberately **not** human: `human_only` already
    # measures the human quartet, and a "one per cell" panel that resolved to
    # the same four baits would be the same experiment reported twice. The
    # pick is the longest non-human bait of each cell, which is the choice a
    # panel builder with no reference species would make.
    one_each = set()
    for clade in S.PARALOGS + (S.CONTROL_CELL,):
        cands = sorted((b for b in ids if meta[b]["clade"] == clade
                        and "Homo_sapiens" not in meta[b]["species"]),
                       key=lambda b: -meta[b]["length"])
        if cands:
            one_each.add(cands[0])
    out["one_per_cell_nonhuman"] = (
        one_each, "one non-human bait per cell — the longest of each")
    return {k: {"baits": v[0], "description": v[1]} for k, v in out.items()}


# ---------------------------------------------------------------- simulation

def _cell_calls(alns: list, keep: set[str]) -> dict[str, dict]:
    """Per-cell evidence from a bait subset, by the sweep's own rules.

    The full chain is reproduced, not approximated: cluster, then apply the
    identity floor (`filter_loci` — MIN_LOCUS_IDENTITY, which S5b added after
    a 2 Mbp `-G` chained shared-module hits into apparent loci), then
    `s5_classify.cell_loci`, which is where D14 lives — a locus won by the
    RyR baits is never offered to an ITPR cell, and only inside the winning
    family is the paralog question asked. Scoring the best coverage over
    every alignment instead would let an ITPR3 bait's alignment at the ITPR1
    gene stand in as ITPR3 evidence, which makes every panel look alike.
    """
    subset = [a for a in alns if a.bait in keep]
    calls = {c: {"coverage": 0.0, "identity": 0.0, "bait": "", "n_loci": 0,
                 "locus": None} for c in S.CELLS}
    if not subset:
        return calls
    loci = SW.filter_loci(SW.cluster_loci(subset))
    for cell in S.CELLS:
        primary, secondary = CL.cell_loci(loci, cell)
        mine = primary + secondary
        if not mine:
            continue
        mine.sort(key=lambda L: (L.best_of(cell) or L.best).score, reverse=True)
        best = mine[0].best_of(cell) or mine[0].best
        calls[cell] = {"coverage": best.coverage, "identity": best.identity,
                       "bait": best.bait, "n_loci": len(mine),
                       "locus": mine[0]}
    return calls


def _single_bait_at(alns: list, locus, bait: str) -> tuple[float, float]:
    """(coverage, identity) this bait alone reaches at an already-found gene.

    Measured **at the locus the full panel assigned to the cell**, not
    wherever the bait aligns best in the genome: a lone ITPR3 bait has no
    competition, so scored genome-wide it would report its hit at the ITPR1
    gene as an ITPR3 recovery. Restricting to the true locus is what makes
    this a measurement of reach into a clade rather than of mis-assignment.
    """
    if locus is None:
        return (0.0, 0.0)
    best = None
    for a in alns:
        if a.bait != bait or a.contig != locus.contig:
            continue
        if a.end < locus.start or a.start > locus.end:
            continue
        if best is None or a.score > best.score:
            best = a
    return (best.coverage, best.identity) if best else (0.0, 0.0)


ALN_FIELDS = ("contig", "start", "end", "strand", "score", "identity", "bait",
              "clade", "family", "q_start", "q_end", "bait_len", "aligned_aa",
              "mp_id")


def genome_alignments(acc: str, meta: dict, use_cache: bool = True) -> list:
    """Every miniprot alignment for one genome, cached in compact form.

    The retained GFFs total 939 MB and re-parsing them once per panel is the
    slow step; the simulation needs a small subset of the fields (no CDS
    blocks, no translations), so they are cached per accession under the data
    root. `cds_blocks` is deliberately not cached and not needed: coverage is
    already computed into `aligned_aa` by the parser.
    """
    cache = S.cache_dir() / "alns" / f"{acc}.tsv"
    if use_cache and cache.exists() and cache.stat().st_size:
        out = []
        for r in S.read_tsv(cache):
            out.append(SW.Aln(
                contig=r["contig"], start=int(r["start"]), end=int(r["end"]),
                strand=r["strand"], score=float(r["score"]),
                identity=float(r["identity"]), bait=r["bait"],
                clade=r["clade"], family=r["family"],
                q_start=int(r["q_start"]), q_end=int(r["q_end"]),
                bait_len=int(r["bait_len"]),
                aligned_aa=int(r["aligned_aa"]), mp_id=r["mp_id"]))
        return out
    gff = S.sweep_dir(acc) / "miniprot.gff"
    if not gff.exists():
        return []
    alns = SW.parse_miniprot_gff(gff, meta)
    cache.parent.mkdir(parents=True, exist_ok=True)
    S.write_tsv(cache, list(ALN_FIELDS),
                [[getattr(a, f) for f in ALN_FIELDS] for a in alns])
    return alns


def run_panels(log=S.log) -> tuple[list[list], dict]:
    meta = S.bait_meta()
    plans = panels(meta)
    ledger = {(c["accession"], c["cell"]): c for c in S.ledger_cells()}
    genomes: dict[str, dict] = {}
    for (acc, _c), row in ledger.items():
        genomes.setdefault(acc, row)
    log(f"panel simulation: {len(genomes)} genomes x {len(plans)} panels")

    rows: list[list] = []
    single: list[list] = []
    per_panel = {n: {"found": 0, "cells": 0} for n in plans}
    mismatches: list[list] = []
    for i, (acc, g) in enumerate(sorted(genomes.items()), 1):
        alns = genome_alignments(acc, meta)
        if not alns:
            log(f"  !! no alignments retained for {acc}")
            continue
        for name, plan in plans.items():
            calls = _cell_calls(alns, plan["baits"])
            for cell in S.CELLS:
                c = calls[cell]
                found = c["coverage"] >= COV_FOUND
                per_panel[name]["cells"] += 1
                per_panel[name]["found"] += int(found)
                rows.append([acc, g["organism"], g["vclass"], cell, name,
                             len(plan["baits"]), int(found),
                             round(c["coverage"], 4), round(c["identity"], 4),
                             c["bait"]])
                if name == "full38":
                    truth = S.is_found(ledger[(acc, cell)]["status"])
                    if truth != found:
                        mismatches.append([acc, g["organism"], cell,
                                           ledger[(acc, cell)]["status"],
                                           int(found), round(c["coverage"], 4),
                                           ledger[(acc, cell)]["best_coverage"]])
        full = _cell_calls(alns, plans["full38"]["baits"])
        for bait, m in meta.items():
            cell = m["clade"]
            if cell not in S.CELLS:
                continue                # an unlabelled bait owns no cell
            if full[cell]["coverage"] < COV_FOUND:
                continue                # no gene here to recover
            cov, ident = _single_bait_at(alns, full[cell]["locus"], bait)
            single.append([acc, g["organism"], g["vclass"], cell, bait,
                           m["species"], m["band"],
                           BAND_GROUP.get(m["band"], "other"), m["origin"],
                           round(ident, 4), round(cov, 4),
                           int(cov > 0), int(cov >= COV_FOUND)])
        if i % 50 == 0:
            log(f"  {i}/{len(genomes)} genomes")

    S.write_tsv(S.out_dir() / "panel_cells.tsv",
                ["accession", "organism", "vclass", "cell", "panel",
                 "n_baits", "found", "best_coverage", "best_identity",
                 "best_bait"], rows)
    S.write_tsv(S.out_dir() / "panel_single_bait.tsv",
                ["accession", "organism", "vclass", "cell", "bait",
                 "bait_species", "bait_band", "bait_group", "bait_origin",
                 "identity", "coverage", "aligned", "found"], single)
    S.write_tsv(S.out_dir() / "panel_validation.tsv",
                ["accession", "organism", "cell", "ledger_status",
                 "simulated_found", "simulated_coverage",
                 "ledger_best_coverage"], mismatches)
    log(f"panel_validation: {len(mismatches)} of {len(genomes) * len(S.CELLS)} "
        f"cells disagree with the ledger under the full panel")
    return rows, {"per_panel": per_panel, "mismatch": len(mismatches),
                  "n_genomes": len(genomes), "plans": plans}


# ------------------------------------------------------------------ summaries

def panel_summary(state: dict, cells: list[list], log=S.log) -> list[list]:
    plans = state["plans"]
    by_panel: dict[str, dict] = {}
    for acc, _org, _vc, cell, panel, _nb, found, *_rest in cells:
        by_panel.setdefault(panel, {})[(acc, cell)] = found
    full = by_panel["full38"]
    base_itpr = sum(1 for (a, c), v in full.items() if c in S.PARALOGS and v)
    rows = []
    for name, plan in plans.items():
        mine = by_panel[name]
        itpr = {k: v for k, v in mine.items() if k[1] in S.PARALOGS}
        n_found = sum(itpr.values())
        lost = sum(1 for k, v in full.items()
                   if k[1] in S.PARALOGS and v and not itpr.get(k, 0))
        gained = sum(1 for k, v in itpr.items()
                     if v and not full.get(k, 0))
        ryr = {k: v for k, v in mine.items() if k[1] == S.CONTROL_CELL}
        rows.append([name, plan["description"], len(plan["baits"]),
                     len(itpr), n_found,
                     round(n_found / max(1, len(itpr)), 4),
                     n_found - base_itpr, lost, gained,
                     sum(ryr.values()), len(ryr)])
    rows.sort(key=lambda r: -r[4])
    S.write_tsv(S.out_dir() / "panel_recall.tsv",
                ["panel", "description", "n_baits", "n_itpr_cells",
                 "n_itpr_found", "itpr_recall", "delta_vs_full", "n_lost",
                 "n_gained", "n_ryr_found", "n_ryr_cells"], rows)
    log("panel_recall: " + ", ".join(f"{r[0]}={r[5]}" for r in rows[:5]))
    return rows


def panel_changes(cells: list[list], log=S.log) -> list[list]:
    """Every cell whose call differs from the full panel, both sides shown.

    Two mechanisms live in this table and they are opposite. A **loss** is a
    cell the removed bait was the only one close enough to place. A **gain**
    is the rule's own edge: the cell's coverage is read off the cell's
    *top-scoring* bait, not its best-covering one, so removing a competitor
    can promote a bait whose coverage is fractionally higher and push a
    marginal cell across `COV_FOUND`. Both are reported with the two
    coverages, so a gain cannot be read as an ablation improving recall.
    """
    full = {(r[0], r[3]): r for r in cells if r[4] == "full38"}
    rows = []
    for r in cells:
        if r[4] == "full38":
            continue
        f = full.get((r[0], r[3]))
        if f is None or f[6] == r[6]:
            continue
        rows.append([r[4], r[0], r[1], r[2], r[3],
                     f[6], f[7], f[9], r[6], r[7], r[9],
                     "gain" if r[6] > f[6] else "loss",
                     round(abs(f[7] - COV_FOUND), 4)])
    S.write_tsv(S.out_dir() / "panel_changes.tsv",
                ["panel", "accession", "organism", "vclass", "cell",
                 "full_found", "full_coverage", "full_bait", "panel_found",
                 "panel_coverage", "panel_bait", "direction",
                 "full_coverage_distance_from_bar"], rows)
    gains = sum(1 for r in rows if r[11] == "gain")
    log(f"panel_changes: {len(rows)} call changes across all panels "
        f"({gains} gains, {len(rows) - gains} losses)")
    return rows


def panel_by_group(cells: list[list], log=S.log) -> list[list]:
    """Where each ablation costs recall — by vertebrate class and cell."""
    idx: dict[tuple, dict] = {}
    for _acc, _org, vclass, cell, panel, _nb, found, *_rest in cells:
        d = idx.setdefault((panel, vclass, cell), {"cells": 0, "found": 0})
        d["cells"] += 1
        d["found"] += found
    base = {k[1:]: v for k, v in idx.items() if k[0] == "full38"}
    rows = []
    for (panel, vclass, cell), d in sorted(idx.items()):
        b = base.get((vclass, cell), {"found": 0})
        rows.append([panel, vclass, cell, d["cells"], d["found"],
                     round(d["found"] / max(1, d["cells"]), 4),
                     d["found"] - b["found"]])
    S.write_tsv(S.out_dir() / "panel_by_group.tsv",
                ["panel", "vclass", "cell", "n_cells", "n_found", "recall",
                 "delta_vs_full"], rows)
    log(f"panel_by_group: {len(rows)} rows")
    return rows


def unfilled_slot_effect(cells: list[list], log=S.log) -> list[list]:
    """What the six unfilled panel slots cost, where they would have acted.

    S5's panel has no labelled ITPR1/ITPR2 bait for chondrichthyans, no
    ITPR2 for the coelacanth grade and none of the three for cyclostomes,
    because no labelled full-length record exists there. The counterfactual
    cannot be run — there is no bait to add — but its consequence can be
    measured from the other side: in exactly those bands, how much of the
    recall the full panel achieves is delivered by an unlabelled bait?
    """
    # The unfilled slots are recorded in S5's own build statistics as
    # "ITPR1@chondrichthyes (0/1)", so they are parsed from there rather
    # than re-derived — a second derivation could disagree about which
    # slots the panel actually failed to fill.
    blob = S.read_json(S.RESULTS / "s5_baits" / "bait_build_stats.json")
    want = set()
    for entry in blob.get("unfilled_slots") or []:
        slot = entry.split(" ")[0]
        if "@" in slot:
            cell, band = slot.split("@", 1)
            want.add((cell, band))
    band_class = {"chondrichthyes": ("Chondrichthyes",),
                  "cyclostomata": ("Hyperoartia", "Myxini"),
                  "sarcopterygian_fish": ("Coelacanthimorpha", "Dipnoi"),
                  "actinopteri": ("Actinopteri",), "aves": ("Aves",),
                  "mammalia": ("Mammalia",), "reptilia": ("Lepidosauria",
                                                          "Testudines",
                                                          "Crocodylia"),
                  "amphibia": ("Amphibia",)}
    full = [r for r in cells if r[4] == "full38"]
    labelled = {(r[0], r[3]): r[6] for r in cells if r[4] == "labelled_only"}
    rows = []
    for cell, band in sorted(want):
        classes = band_class.get(band, ())
        sub = [r for r in full if r[3] == cell and r[2] in classes]
        if not sub:
            rows.append([cell, band, ";".join(classes), 0, 0, 0, "",
                         "no genome of this band in scope"])
            continue
        n_found = sum(r[6] for r in sub)
        n_lab = sum(labelled.get((r[0], r[3]), 0) for r in sub)
        rows.append([cell, band, ";".join(classes), len(sub), n_found, n_lab,
                     n_found - n_lab,
                     "cells the unlabelled baits carry in a band with no "
                     "labelled bait for this paralog"])
    S.write_tsv(S.out_dir() / "panel_unfilled_slots.tsv",
                ["cell", "band", "vclasses", "n_cells", "n_found_full_panel",
                 "n_found_labelled_only", "n_carried_by_unlabelled", "note"],
                rows)
    log(f"panel_unfilled_slots: {len(rows)} slots")
    return rows


IDENTITY_BINS = [(0.0, 0.50), (0.50, 0.60), (0.60, 0.70), (0.70, 0.80),
                 (0.80, 0.90), (0.90, 0.95), (0.95, 1.01)]


def design_rule(log=S.log) -> tuple[list[list], dict]:
    """Single-bait recall as a function of bait-target identity."""
    single = S.read_tsv(S.out_dir() / "panel_single_bait.tsv")
    no_aln = [r for r in single if r["aligned"] == "0"]
    aligned = [r for r in single if r["aligned"] == "1"]
    rows = [["no alignment", len(no_aln), 0, 0.0, 0.0, 0.0,
             len({r["accession"] for r in no_aln}),
             len({r["bait"] for r in no_aln})]]
    for lo, hi in IDENTITY_BINS:
        sub = [r for r in aligned if lo <= S.fnum(r["identity"], float, 0) < hi]
        if not sub:
            continue
        found = sum(1 for r in sub if r["found"] == "1")
        wlo, whi = S.wilson(found, len(sub))
        rows.append([f"{lo:.2f}-{hi:.2f}", len(sub), found,
                     round(found / len(sub), 4), round(wlo, 4), round(whi, 4),
                     len({r["accession"] for r in sub}),
                     len({r["bait"] for r in sub})])
    S.write_tsv(S.out_dir() / "panel_identity_recall.tsv",
                ["identity_bin", "n_measurements", "n_found", "recall",
                 "wilson_lo", "wilson_hi", "n_genomes", "n_baits"], rows)

    # Per cell, the *lowest* identity at which some single bait still
    # sufficed. This is the form the design rule is quoted in: not "recall at
    # identity x", but "one bait at identity >= x recovered the gene in N% of
    # the cells that have one".
    per_cell: dict[tuple, float] = {}
    for r in single:
        if r["found"] != "1":
            continue
        key = (r["accession"], r["cell"])
        ident = S.fnum(r["identity"], float, 0.0)
        per_cell[key] = min(per_cell.get(key, 1.0), ident)
    vals = sorted(per_cell.values())
    pct = {}
    for q in (50, 75, 90, 95, 99, 100):
        if vals:
            pct[f"p{q}"] = round(
                vals[min(len(vals) - 1, int(round(q / 100 * (len(vals) - 1))))],
                4)
    S.write_tsv(S.out_dir() / "panel_min_identity.tsv",
                ["accession", "cell", "min_sufficient_identity"],
                [[a, c, round(v, 4)] for (a, c), v in sorted(per_cell.items())])
    rule = next((r[0] for r in rows[1:] if r[4] >= 0.95), None)
    log(f"design_rule: single-bait recall reaches 95 % (lower CI) at "
        f"identity {rule or 'not reached'}; percentiles {pct}")
    return rows, {"identity_bin_for_95pct_recall": rule,
                  "n_measurements": len(single),
                  "n_no_alignment": len(no_aln), "n_aligned": len(aligned),
                  "min_sufficient_identity_percentiles": pct,
                  "n_cells_with_single_bait_rescue": len(per_cell)}


def run(log=S.log) -> dict:
    cells, state = run_panels(log)
    recall = panel_summary(state, cells, log)
    panel_changes(cells, log)
    panel_by_group(cells, log)
    unfilled_slot_effect(cells, log)
    _, rule = design_rule(log)
    summary = {
        "n_genomes": state["n_genomes"],
        "full_panel_vs_ledger_mismatches": state["mismatch"],
        "cov_found": COV_FOUND,
        "panels": {r[0]: {"n_baits": r[2], "itpr_recall": r[5],
                          "delta_vs_full": r[6], "n_lost": r[7],
                          "n_gained": r[8]} for r in recall},
        "design_rule": rule,
    }
    S.write_json(S.out_dir() / "panel_summary.json", summary)
    return summary


if __name__ == "__main__":
    run()
