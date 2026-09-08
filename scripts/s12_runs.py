"""S12 step 2 — choose the public RNA-seq runs to quantify.

Selection is keyword-driven but **attribute-verified**: a tissue keyword
queries SRA, and a returned run is kept only if its own sample attributes
name that same tissue.  A run that merely mentions "brain" in a study title
is recorded with `title` as its tissue evidence, so the report can say how
many of its tissue calls rest on a curated attribute and how many on prose.

Filters, all recorded in `runs_considered.tsv` so the selection is auditable:
  * LIBRARY_STRATEGY = RNA-Seq, LIBRARY_SOURCE = TRANSCRIPTOMIC
  * Illumina platform, read length >= MIN_READ_LEN
  * >= MIN_SPOTS spots — a shallow run cannot support an absence claim
  * at most `--per-tissue` runs per tissue, preferring distinct studies so
    one lab's library prep cannot carry a tissue on its own

Two things differ from the PIEZO port, and both come from the species.
These are non-model fish, so a **tissue-blind sweep runs first**: the
tissue-keyword queries are a way of *spreading* the selection, not the only
route in, and a species whose libraries carry no tissue metadata at all
would otherwise select nothing.  And the per-species run counts are
measured and written out, because P2 in `s12_panel` excludes a species for
want of libraries and that exclusion has to rest on a number in a table.

Outputs (committed):
    results/expression/runs_considered.tsv
    results/expression/runs_selected.tsv
    results/expression/run_availability.tsv

Run:  python scripts/s12_runs.py [--species ...] [--per-tissue 2]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s12_panel  # noqa: E402
from s12_lib import (  # noqa: E402
    OUT_DIR, TISSUES, SraRun, assign_tissue, efetch_sra, esearch_count,
    esearch_sra, live, write_tsv,
)

MIN_SPOTS = 3_000_000
MIN_READ_LEN = 50
MAX_BLIND = 24      # tissue-blind runs kept per species

BASE_TERM = ('"{sp}"[Organism] AND "rna seq"[Strategy] '
             'AND "transcriptomic"[Source] AND "illumina"[Platform]')

CONSIDERED_COLS = ["species", "organism", "query", "run", "study", "sample",
                   "strategy", "source", "selection", "layout", "platform",
                   "model", "spots", "bases", "read_len", "tissue_called",
                   "tissue_evidence", "verdict", "study_title"]
SELECTED_COLS = ["species", "organism", "run", "tissue", "tissue_evidence",
                 "tissue_from", "study", "study_title", "sample", "layout",
                 "model", "spots", "read_len", "selection"]
AVAIL_COLS = ["species", "organism", "term", "n_runs"]


def _verdict(r: SraRun, want: str | None) -> str:
    if r.strategy.upper() != "RNA-SEQ":
        return "reject:strategy"
    if r.platform.upper() != "ILLUMINA":
        return "reject:platform"
    if r.spots < MIN_SPOTS:
        return "reject:shallow"
    if r.read_len < MIN_READ_LEN:
        return "reject:read_len"
    if not assign_tissue(r, want=want):
        return "reject:tissue_unconfirmed" if want else "keep:tissue_unknown"
    return "keep"


def _tissue_from(r: SraRun) -> str:
    """Whether the tissue came from a curated attribute or from free text."""
    key = r.tissue_source.split("=")[0]
    return "title" if key in ("sample_title", "sample_description",
                              "experiment_title") else "attribute"


def search_species(sp, per_tissue: int, retmax: int, log=print
                   ) -> tuple[list[dict], list[SraRun], list[dict]]:
    considered: list[dict] = []
    selected: list[SraRun] = []
    seen: set[str] = set()
    avail: list[dict] = []

    def record(r: SraRun, query: str, verdict: str) -> None:
        considered.append({
            "species": sp.short, "organism": sp.name, "query": query,
            "run": r.run, "study": r.study, "sample": r.sample,
            "strategy": r.strategy, "source": r.source,
            "selection": r.selection, "layout": r.layout,
            "platform": r.platform, "model": r.model, "spots": r.spots,
            "bases": r.bases, "read_len": r.read_len,
            "tissue_called": r.tissue or "",
            "tissue_evidence": r.tissue_source.replace("\t", " "),
            "verdict": verdict,
            "study_title": r.study_title.replace("\t", " ")[:160]})

    base = BASE_TERM.format(sp=sp.name)
    avail.append({"species": sp.short, "organism": sp.name, "term": "base",
                  "n_runs": esearch_count(base)})

    # Tissue-directed queries first: they spread the selection across the
    # organs the species has libraries for.
    for tissue in TISSUES:
        term = f"{base} AND {tissue}[All Fields]"
        uids = esearch_sra(term, retmax=retmax)
        if not uids:
            continue
        kept: list[SraRun] = []
        studies: set[str] = set()
        for r in efetch_sra(uids):
            v = _verdict(r, tissue)
            record(r, tissue, v)
            if v != "keep" or r.run in seen:
                continue
            if len(kept) < per_tissue and r.study not in studies:
                kept.append(r)
                studies.add(r.study)
        for r in kept:
            seen.add(r.run)
            selected.append(r)
        if kept:
            log(f"    {tissue:13s} {len(kept)} run(s) "
                f"[{', '.join(r.run for r in kept)}]")

    # Then a tissue-blind pass, so a species whose submitters recorded no
    # organ still gets libraries.  These runs carry tissue "unknown" and the
    # report counts them separately rather than folding them into an organ.
    #
    # Capped, and the cap is spread over studies before depth.  *Nibea
    # albiflora* alone returns 114 qualifying runs and streaming them all
    # would spend hours re-asking one question of one lab's libraries;
    # taking the deepest run of each distinct study first buys independence
    # instead, which is what an absence claim actually needs.
    uids = esearch_sra(base, retmax=retmax * 4)
    pool: list[SraRun] = []
    for r in efetch_sra(uids):
        if r.run in seen:
            continue
        v = _verdict(r, None)
        record(r, "base", v)
        if v.startswith("keep"):
            pool.append(r)
    pool.sort(key=lambda x: -x.spots)
    blind: list[SraRun] = []
    for want_new_study in (True, False):
        studies = {x.study for x in blind}
        for r in pool:
            if len(blind) >= MAX_BLIND:
                break
            if r.run in seen:
                continue
            if want_new_study and r.study in studies:
                continue
            blind.append(r)
            studies.add(r.study)
            seen.add(r.run)
    for r in blind:
        if r.tissue is None:
            r.tissue, r.tissue_source = "unknown", ""
        selected.append(r)
    if blind:
        log(f"    {'(tissue-blind)':13s} {len(blind)} further run(s) from "
            f"{len({x.study for x in blind})} stud"
            f"{'y' if len({x.study for x in blind}) == 1 else 'ies'}")
    return considered, selected, avail


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--species", default="")
    ap.add_argument("--per-tissue", type=int, default=2)
    ap.add_argument("--retmax", type=int, default=40)
    args = ap.parse_args()

    panel, _ = s12_panel.build_panel()
    if args.species:
        want = set(args.species.split(","))
        panel = [s for s in panel if s.short in want]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    considered: list[dict] = []
    selected: list[dict] = []
    avail: list[dict] = []

    for i, sp in enumerate(panel, 1):
        print(f"[runs] {sp.name}", flush=True)
        live("runs", i - 1, len(panel), sp.name)
        con, sel, av = search_species(sp, args.per_tissue, args.retmax)
        considered += con
        avail += av
        for r in sel:
            selected.append({
                "species": sp.short, "organism": sp.name, "run": r.run,
                "tissue": r.tissue or "unknown",
                "tissue_evidence": r.tissue_source.replace("\t", " "),
                "tissue_from": _tissue_from(r) if r.tissue_source else "none",
                "study": r.study,
                "study_title": r.study_title.replace("\t", " ")[:160],
                "sample": r.sample, "layout": r.layout, "model": r.model,
                "spots": r.spots, "read_len": r.read_len,
                "selection": r.selection})
        tissues = sorted({r["tissue"] for r in selected
                          if r["species"] == sp.short})
        print(f"  -> {len([r for r in selected if r['species'] == sp.short])}"
              f" runs across {len(tissues)} tissue labels: "
              f"{', '.join(tissues)}")

    write_tsv(OUT_DIR / "runs_considered.tsv", CONSIDERED_COLS, considered)
    write_tsv(OUT_DIR / "runs_selected.tsv", SELECTED_COLS, selected)
    write_tsv(OUT_DIR / "run_availability.tsv", AVAIL_COLS, avail)
    live("runs", len(panel), len(panel), "done")
    print(f"\nwrote runs_selected.tsv ({len(selected)} runs), "
          f"runs_considered.tsv ({len(considered)} rows)")


if __name__ == "__main__":
    main()
