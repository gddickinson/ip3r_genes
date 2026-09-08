"""The S10 driver — ordered stages, each resumable, each writing tables only.

    select -> evidence -> orf -> tile -> probes -> figures -> report

`--only` / `--from` / `--list` as in `s9_run.py`, and for the same reason: a
stage that fails does not stop the ones after it that do not depend on it. The
report renders *not run yet* for a section whose table is absent, so a partial
S10 says which parts are partial.

`--no-remote` skips the network half of `probes` (the transcript search) and
keeps the offline half (exon-boundary concordance), which is what makes the
whole task re-runnable on a machine with no internet.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_case_spec as spec                                   # noqa: E402
import s10_context as ctx                                      # noqa: E402
import s10_evidence as ev                                      # noqa: E402
import s10_orf                                                 # noqa: E402
import s10_probes as pr                                        # noqa: E402
import s10_tables as tb                                        # noqa: E402
import s10_tile as tile                                        # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "annotation_bugs"
STAGES = ["select", "evidence", "orf", "tile", "probes", "figures", "report"]


def _data_root() -> Path:
    sys.path.insert(0, str(ROOT))
    from src.utils.data_root import require_data_root
    return require_data_root()


def _live(stage: str, done: int, total: int) -> None:
    """Feed the dashboard's live panel (session protocol, during-session)."""
    try:
        (ROOT / "results" / "session_live.json").write_text(json.dumps({
            "task": "S10", "workers": 1,
            "steps": [{"label": s, "done": (STAGES.index(s) < done)}
                      for s in STAGES],
            "current": stage, "progress": f"{done}/{total}"}) + "\n")
    except OSError:
        pass


def load_cases() -> list[dict]:
    path = OUT / "cases.tsv"
    if not path.exists():
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))


def load_ranking() -> list[dict]:
    path = OUT / "case_ranking.tsv"
    if not path.exists():
        return []
    return list(csv.DictReader(open(path), delimiter="\t"))


# --------------------------------------------------------------------------
def stage_select() -> None:
    subprocess.run([sys.executable, str(Path(__file__).parent / "s10_select.py")],
                   check=True)


def stage_evidence(dr: Path) -> None:
    exons, introns, models, blocks, neighbours = [], [], [], [], []
    audits, family, dbrecs, flankchk = [], [], [], []
    archive = dr / "raw_api" / "ncbi_datasets" / "s10"
    for case in load_cases():
        acc = case["accession"]
        gdir, sdir = dr / "genomes" / acc, dr / "genome_sweep" / acc
        g = ev.gather(case, gdir, sdir)
        exons += tb.stamp(g["exons"], case)
        introns += tb.stamp(g["introns"], case)
        models += tb.stamp(g["models"], case)
        blocks += tb.stamp([g["blocks"]], case)
        neighbours += tb.stamp(g["neighbourhood"], case)
        flankchk += ctx.flank_consensus_check(
            case, g["neighbourhood"],
            ROOT / "results" / "synteny" / "flank_consensus.tsv")
        audits.append(ctx.assembly_audit(case, gdir, archive))
        family += [{"case_id": case["case_id"], "accession": acc, **r}
                   for r in ctx.family_named_models(sdir / "genes_slim.tsv")]
        dbrecs += [{"case_id": case["case_id"], "species": case["organism"],
                    **r} for r in ctx.database_records(
                        case["organism"],
                        ROOT / "results" / "census_v3" / "census_v3.tsv")]
    tb.write_tsv(OUT / "case_exons.tsv", exons, tb.EXON_COLS)
    tb.write_tsv(OUT / "case_introns.tsv", introns, tb.INTRON_COLS)
    tb.write_tsv(OUT / "annotated_models.tsv", models, tb.MODEL_COLS)
    tb.write_tsv(OUT / "block_accounting.tsv", blocks, tb.BLOCK_COLS)
    tb.write_tsv(OUT / "locus_neighbourhood.tsv", neighbours, tb.NEIGHBOUR_COLS)
    tb.write_tsv(OUT / "assembly_audit.tsv", audits, ctx.AUDIT_COLS)
    tb.write_tsv(OUT / "family_named_models.tsv", family, ctx.FAMILY_MODEL_COLS)
    tb.write_tsv(OUT / "database_records.tsv", dbrecs, ctx.DB_RECORD_COLS)
    tb.write_tsv(OUT / "flank_consensus_check.tsv", flankchk,
                 ctx.FLANK_CHECK_COLS)
    print(f"[s10] evidence: {len(exons)} exons, {len(introns)} introns, "
          f"{len(models)} annotated models over {len(load_cases())} cases")


def stage_orf(dr: Path, threads: int) -> None:
    """The reading-frame test, run on the case locus *and* on every other
    family locus of the same genome.

    The comparison set is the point: Case A's genome files two complete
    receptor gene models as pseudogenes, and a table holding only the case
    locus could not show that the demotion is systematic rather than a judgment
    about one gene.
    """
    rows = []
    for case in load_cases():
        acc = case["accession"]
        gdir, sdir = dr / "genomes" / acc, dr / "genome_sweep" / acc
        summary = json.loads((sdir / "summary.json").read_text())
        for cell_name, cell in summary["cells"].items():
            loci = cell.get("loci") or []
            if not loci:
                continue
            best = loci[0]
            sub = {**case, "cell": cell_name, "contig": best["contig"],
                   "mp_id": best["mp_id"],
                   "start": best["start"], "end": best["end"],
                   "strand": best["strand"], "identity": best["identity"],
                   "frameshifts": best.get("frameshifts", 0),
                   "stop_codons": best.get("stop_codons", 0)}
            row = s10_orf.analyse(sub, gdir, sdir, threads=threads)
            row["is_case_locus"] = int(cell_name == case["cell"])
            row["cell_status"] = cell["status"]
            rows.append(row)
            print(f"[s10]   orf {acc} {cell_name}: {row.get('verdict')} "
                  f"stops={row.get('internal_stops')}")
    tb.write_tsv(OUT / "reading_frame.tsv", rows,
                 s10_orf.ORF_COLS + ["is_case_locus", "cell_status"])


def stage_tile(dr: Path) -> None:
    """Tile every family-relevant annotated protein onto the genome's loci.

    The query set is the annotated models inside the case locus *plus* every
    family-named model anywhere in the genome. The second half is what
    adjudicates a naming disagreement: a model called ITPR2 five megabases away
    on another chromosome is tiled against all of this genome's recovered loci
    and the sequence decides which one it is a piece of.
    """
    rows = []
    for case in load_cases():
        acc = case["accession"]
        gdir, sdir = dr / "genomes" / acc, dr / "genome_sweep" / acc
        fna = next(gdir.glob("*_genomic.fna"))
        gff = gdir / "genomic.gff.gz"
        subjects = tile.read_fasta(sdir / "novel_models.faa")
        # Add every locus the sweep placed in a family cell, not only the
        # census-v4 non-redundant ones — see `tile.model_protein`.
        from s5_genome_io import build_fai, read_fai
        idx = read_fai(build_fai(fna))
        summary = json.loads((sdir / "summary.json").read_text())
        for cell_name, cell in summary["cells"].items():
            for loc in cell.get("loci") or []:
                key = (f"{acc}|{cell_name}|{loc['contig']}:{loc['start']}-"
                       f"{loc['end']}{loc['strand']}|sweep")
                if any(k.startswith(key.rsplit("|", 1)[0]) for k in subjects):
                    continue
                m = ev.read_miniprot_model(sdir / "miniprot.gff", loc["mp_id"])
                if m.get("cds"):
                    subjects[key] = tile.model_protein(fna, idx, m)
        g = ev.gather(case, gdir, sdir)
        queries = tile.annotated_proteins(g["genes_inside"], fna,
                                          g["model"]["contig"])
        seen = {q["gene_id"] for q in queries}
        for fam in ctx.family_named_models(sdir / "genes_slim.tsv"):
            gid = fam["gene_id"]
            if gid in seen:
                continue
            win = ev.read_annotation_window(gff, fam["contig"],
                                            fam["start"], fam["end"])
            hits = [x for x in win["genes"] if x["gene_id"] == gid]
            queries += tile.annotated_proteins(hits, fna, fam["contig"])
            seen.add(gid)
        for r in tile.tile(queries, subjects):
            r["at_own_locus"] = tile.at_own_locus(r)
            r["naming"] = tile.naming_verdict(r, tile.NAME_HINTS)
            rows.append({"case_id": case["case_id"], "accession": acc, **r})
        print(f"[s10]   tile {acc}: {len(queries)} annotated proteins")
    tb.write_tsv(OUT / "fragment_tiling.tsv", rows, tile.TILE_COLS)


def stage_probes(dr: Path, remote: bool) -> None:
    """Junction probes against the species' own records, plus the offline
    boundary control.

    `remote` gates only the *fetch*: once a species' records are archived the
    search is local and always runs, which is what makes the whole task
    re-derive without a network.
    """
    probes_rows, conc_rows, summaries = [], [], []
    ranking = load_ranking()
    for case in load_cases():
        acc = case["accession"]
        gdir, sdir = dr / "genomes" / acc, dr / "genome_sweep" / acc
        audit = next((r for r in csv.DictReader(
            open(OUT / "assembly_audit.tsv"), delimiter="\t")
            if r["case_id"] == case["case_id"]), {})
        taxid = int(audit.get("taxid") or 0)
        cache = dr / "annotation_bugs" / case["case_id"]
        g = ev.gather(case, gdir, sdir)
        probes = pr.build_probes(g["region"], g["model"], g["introns"])

        sets = (pr.fetch_species_nucleotides(taxid, cache) if remote
                else {"mrna": {"status": "cached",
                               "path": cache / "nuccore_mrna.fasta",
                               "n_records": (cache / "nuccore_mrna.fasta"
                                             ).read_text().count(">")
                               if (cache / "nuccore_mrna.fasta").exists() else 0,
                               "n_total": 0}})
        # The spanning criterion's own negative control, built from the locus
        # this case is about: many hits, none of which may span.
        sets["genomic_control"] = {
            "status": "local", "path": pr.genomic_control(g["region"], cache),
            "n_records": 1, "n_total": 1}
        res = {}
        for tag, info in sets.items():
            hit = pr.search_probes(probes, info.get("path"), cache, tag)
            res[f"n_records_{tag}"] = info.get("n_records", 0)
            res[f"n_total_{tag}"] = info.get("n_total", 0)
            res[f"status_{tag}"] = hit["status"]
            res[f"spanning_{tag}"] = sum(
                v["n_spanning"] for v in hit["by_probe"].values())
            res[f"hits_{tag}"] = sum(
                v["n_hits"] for v in hit["by_probe"].values())
            for p in probes:
                r = hit["by_probe"].get(
                    p["junction_index"],
                    {"status": hit["status"], "n_hits": 0, "n_spanning": 0,
                     "best_identity": "", "best_subject": ""})
                probes_rows.append({
                    "case_id": case["case_id"], "accession": acc, "db": tag,
                    **{k: v for k, v in p.items() if k != "seq"}, **r})

        refs = pr.reference_models(case, ranking, dr / "genome_sweep")
        rows, summ = pr.concordance(g["model"], refs)
        conc_rows += [{"case_id": case["case_id"], "accession": acc, **r}
                      for r in rows]
        summaries.append({"case_id": case["case_id"], "accession": acc,
                          "cell": case["cell"], "taxid": taxid,
                          "n_probes": len(probes),
                          "n_reference_refseq": sum(
                              1 for r in refs
                              if r["accession"].startswith("GCF_")),
                          "reference_genomes": ";".join(
                              r["accession"] for r in refs),
                          **summ, **res})
        print(f"[s10]   probes {acc}: {len(probes)} junctions; "
              f"mrna {res.get('spanning_mrna', 0)}/{res.get('hits_mrna', 0)} "
              f"spanning/hits over {res.get('n_records_mrna', 0)} records; "
              f"control {res.get('spanning_other', 0)}/"
              f"{res.get('hits_other', 0)} over "
              f"{res.get('n_records_other', 0)}; "
              f"{summ['frac_shared_with_majority']:.1%} boundaries shared")
    tb.write_tsv(OUT / "junction_probes.tsv", probes_rows, pr.PROBE_COLS)
    tb.write_tsv(OUT / "boundary_concordance.tsv", conc_rows,
                 pr.CONCORDANCE_COLS)
    tb.write_tsv(OUT / "probe_summary.tsv", summaries,
                 sorted({k for r in summaries for k in r},
                        key=lambda k: (k != "case_id", k)))


def _write_stats(args) -> None:
    """`annotation_bugs_stats.json` — parameters, self-tests, table hashes."""
    import s10_case_spec as sp
    cases = load_cases()
    tables = sorted(OUT.glob("*.tsv"))
    tb.write_stats(OUT, {
        "task": "S10",
        "rules": {code: text for code, text in sp.ELIGIBILITY},
        "modes": {m: sp.MODE_RULE[m] for m in sp.FAILURE_MODES},
        "parameters": {
            "COV_FOUND": sp.COV_FOUND, "HEADROOM_X": sp.HEADROOM_X,
            "MIN_LONGER": sp.MIN_LONGER,
            "ANNOT_CDS_FRAC": sp.ANNOT_CDS_FRAC,
            "FRAG_MIN_FRAC_CDS": sp.FRAG_MIN_FRAC_CDS,
            "FLANK_BP": ev.FLANK_BP,
            "PROBE_FLANK_NT": pr.PROBE_FLANK_NT,
            "MIN_ANCHOR": pr.MIN_ANCHOR,
            "MIN_PROBE_COVERAGE": pr.MIN_PROBE_COVERAGE,
            "BOUNDARY_TOL_AA": pr.BOUNDARY_TOL_AA,
            "remote_search": not args.no_remote,
        },
        "self_tests": "s10_test_evidence.py passed before this run",
        "cases": [{k: c[k] for k in ("case_id", "selected_as", "accession",
                                     "organism", "cell", "loss")}
                  for c in cases],
    }, tables)


def _run_module(name: str) -> None:
    subprocess.run([sys.executable, str(Path(__file__).parent / name)],
                   check=True)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", choices=STAGES)
    ap.add_argument("--from", dest="from_stage", choices=STAGES)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--no-remote", action="store_true",
                    help="skip the network transcript search; keep the "
                         "offline exon-boundary concordance")
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    if args.list:
        for s in STAGES:
            print(s)
        return 0
    todo = args.only or (STAGES[STAGES.index(args.from_stage):]
                         if args.from_stage else STAGES)

    # The negative controls run before anything is written, as
    # `s8_run_synteny.py` does. Every rule in S10 produces a plausible number
    # when it is wrong, so a build that cannot fail on purpose does not run.
    from s10_test_evidence import main as run_tests
    if run_tests() != 0:
        print("[s10] negative controls failed — refusing to write tables")
        return 2

    dr = _data_root()
    failed = []
    for i, stage in enumerate(todo):
        _live(stage, STAGES.index(stage), len(STAGES))
        t0 = time.time()
        print(f"[s10] === {stage} ===")
        try:
            if stage == "select":
                stage_select()
            elif stage == "evidence":
                stage_evidence(dr)
            elif stage == "orf":
                stage_orf(dr, args.threads)
            elif stage == "tile":
                stage_tile(dr)
            elif stage == "probes":
                stage_probes(dr, not args.no_remote)
            elif stage == "figures":
                _run_module("s10_figures.py")
            elif stage == "report":
                _run_module("s10_report.py")
        except Exception as exc:                      # noqa: BLE001
            print(f"[s10] !! {stage} failed: {type(exc).__name__}: {exc}")
            failed.append(stage)
        print(f"[s10] {stage} in {time.time() - t0:.1f}s")
    _write_stats(args)
    _live("done", len(STAGES), len(STAGES))
    if failed:
        print(f"[s10] failed stages: {', '.join(failed)}")
        return 1
    print(f"[s10] complete — {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
