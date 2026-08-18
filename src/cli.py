"""Headless CLI for running a preset and saving a results bundle.

Used by `python run.py --headless --preset NAME [--save-results]`. Keeps the
GUI module out of the import path so this works on machines without a
display.
"""

from __future__ import annotations

import json
import queue
import time
from pathlib import Path
from typing import Optional

from .analysis import analyse
from .analysis.pipeline import write_analysis
from .core.cache import DiskCache
from .core.models import ProteinVariant, SearchQuery, SearchResult, SearchStatus
from .core.search import SearchOrchestrator
from .databases.interpro import batch_fetch_family_signatures
from .discovery import DiscoveryConfig, discover_novel_paralogs, write_discovery
from .investigation import (
    Investigator,
    InvestigationOptions,
    write_investigation,
)
from .utils.report import write_report
from .utils.results_writer import make_bundle_dir, write_bundle


def load_preset(presets_dir: Path, name: str) -> dict:
    path = presets_dir / f"{name}.json"
    if not path.exists():
        raise FileNotFoundError(f"preset not found: {path}")
    return json.loads(path.read_text())


def query_from_preset(preset: dict, species_override: str = "", max_override: int = 0) -> SearchQuery:
    return SearchQuery(
        gene_symbols=list(preset.get("genes", [])),
        species=species_override or preset.get("species", "") or "",
        sources=list(preset.get("sources", ["NCBI", "Ensembl", "UniProt"])),
        max_results_per_source=max_override or int(preset.get("max_per_source", 50)),
        include_sequence=bool(preset.get("include_sequence", False)),
    )


def run_headless(
    project_root: Path,
    preset: str,
    email: str = "",
    species_override: str = "",
    max_override: int = 0,
    save_results: bool = False,
    timeout_s: int = 300,
    label_override: str = "",
    analyze: bool = False,
    discover: bool = False,
    known_paralogs: Optional[list[str]] = None,
    interpro_lookup: bool = False,
) -> dict:
    presets_dir_local = project_root / "presets"
    preset_data_peek = load_preset(presets_dir_local, preset)
    if preset_data_peek.get("mode") == "domain_scan":
        from .discovery.domain_scan import run_domain_scan
        # In domain-scan mode the whole point is to score candidates —
        # enable analyse + discover by default. The user can still pass
        # --no-analyze / --no-discover semantics via the preset itself.
        return run_domain_scan(
            project_root=project_root,
            preset_data=preset_data_peek,
            preset_name=preset,
            email=email,
            save_results=save_results,
            label_override=label_override,
            analyze=True,
            discover=True,
        )
    if preset_data_peek.get("mode") == "exhaustive":
        from .discovery.exhaustive import run_exhaustive_hunt
        return run_exhaustive_hunt(
            project_root=project_root,
            preset_data=preset_data_peek,
            preset_name=preset,
            email=email,
            save_results=save_results or True,
            label_override=label_override,
        )
    presets_dir = project_root / "presets"
    preset_data = load_preset(presets_dir, preset)
    query = query_from_preset(preset_data, species_override, max_override)
    print(f"[headless] preset={preset}  genes={query.gene_symbols}  "
          f"species={query.species or '(any)'}  sources={query.sources}  "
          f"max={query.max_results_per_source}")

    cache = DiskCache(project_root / "cache")
    # 90 s: Ensembl REST answers a symbol lookup in under a second most
    # days and in ~55 s when loaded (measured 2026-08-18). 45 s cost a
    # whole source on the first smoke test.
    orch = SearchOrchestrator(cache=cache, email=email, timeout_s=90)
    started = orch.start(query)
    if started == 0:
        print("[headless] no workers started (no enabled sources)")
        return {}

    results: list[SearchResult] = []
    deadline = time.time() + timeout_s
    while len(results) < started and time.time() < deadline:
        try:
            r = orch.queue.get(timeout=5)
            results.append(r)
            print(f"[headless]   [{r.status.value:5s}] {r.source:9s} "
                  f"{len(r.variants):4d} hits  ({r.elapsed_s:.1f}s)  {r.message}")
        except queue.Empty:
            print(f"[headless]   …waiting (got {len(results)}/{started})")

    all_variants = [v for r in results for v in r.variants]
    print(f"[headless] total: {len(all_variants)} variants across {len(results)} source(s)")

    analysis_result = None
    discovery_report = None
    discovery_report_text = ""
    if analyze:
        print("[headless] running sequence analysis…")
        analysis_result = analyse(all_variants, identity_threshold=0.5)
        print(f"[headless]   analysed {analysis_result.n_analyzed}/{analysis_result.n_input} variants with sequences")

    if discover:
        print("[headless] running novel-paralog discovery…")
        known = known_paralogs or list(preset_data.get("known_paralogs", [])) \
            or list(query.gene_symbols)
        domain_hits: dict = {}
        if interpro_lookup:
            print("[headless]   querying InterPro for Pfam signatures…")
            uniprot_accs = sorted({v.accession for v in all_variants if v.source == "UniProt"})
            domain_hits = batch_fetch_family_signatures(uniprot_accs)
            n_with = sum(1 for hits in domain_hits.values() if hits)
            print(f"[headless]   {n_with}/{len(uniprot_accs)} UniProt entries carry the family signature")

        # Build the label-mapping the discovery scorer needs.
        from .analysis.pipeline import _label_for
        from .discovery import build_signature_set
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v) for v in all_variants}
        foldseek_hits = [v for v in all_variants if v.source == "Foldseek"]
        sset = None
        if analysis_result and analysis_result.msa:
            sset = build_signature_set(all_variants, analysis_result.msa, label_map, known)
            if sset:
                print(f"[headless]   derived {len(sset.signatures)} family-signature blocks")
        # Presets may tune scorer thresholds via a "discovery" block, e.g.
        # {"discovery": {"novelty_ceiling": 0.55}}.
        discovery_report = discover_novel_paralogs(
            all_variants,
            DiscoveryConfig(known_paralogs=known, **preset_data.get("discovery", {})),
            analysis_label_for=label_map,
            distances=(analysis_result.distances if analysis_result else None),
            domain_hits=domain_hits,
            foldseek_hits=foldseek_hits,
            signature_set=sset,
        )
        discovery_report_text = discovery_report.text_summary()
        print(f"[headless]   discovery: "
              f"{sum(1 for c in discovery_report.candidates if c.score >= 60)} promising, "
              f"{sum(1 for c in discovery_report.candidates if 40 <= c.score < 60)} worth review")

    if save_results:
        results_root = project_root / "results"
        results_root.mkdir(exist_ok=True)
        label = label_override or "_".join(query.gene_symbols[:3]) or preset
        bundle_dir = make_bundle_dir(results_root, label)
        summary = write_bundle(
            bundle_dir, all_variants, query, results,
            notes=f"Headless run of preset '{preset}'.",
        )
        if analysis_result:
            analysis_dir = bundle_dir / "analysis"
            write_analysis(analysis_dir, analysis_result)
        if discovery_report is not None:
            write_discovery(bundle_dir / "discovery", discovery_report)
        report_paths = write_report(
            bundle_dir, all_variants, query,
            analysis=analysis_result,
            discovery_summary=discovery_report_text,
            notes=f"Headless run of preset '{preset}'.",
        )
        print(f"[headless] bundle saved → {summary['dir']}")
        print(f"[headless]   report.md  : {report_paths['report_md']}")
        print(f"[headless]   report.html: {report_paths['report_html']}")
        summary.update(report_paths)
        return summary
    return {"n_variants": len(all_variants)}


def run_investigate(
    project_root: Path,
    accession: str,
    email: str = "",
    panel_csv: str = "",
    run_foldseek: bool = False,
    save_results: bool = False,
    label_override: str = "",
) -> dict:
    """Deep-dive investigation of a single candidate accession.

    Result lands in results/<timestamp>_investigate_<acc>/investigations/<acc>/
    when save_results is on.
    """
    panel: list[tuple[str, str]] = []
    if panel_csv:
        for acc in panel_csv.split(","):
            acc = acc.strip()
            if acc:
                panel.append((f"panel_{acc}", acc))
    options = InvestigationOptions(
        panel=panel or list(InvestigationOptions().panel),
        run_foldseek=run_foldseek,
        run_literature=True,
        run_synteny=True,
        run_phylo=True,
        email=email,
    )
    print(f"[investigate] accession={accession}  "
          f"panel={[l for l,_ in options.panel]}  "
          f"foldseek={'on' if run_foldseek else 'off'}")

    def progress(label: str) -> None:
        print(f"[investigate]   {label}")

    result = Investigator(accession, options).run(on_progress=progress)

    from .investigation.synthesis import assess_signals, verdict_from_signals
    signals = assess_signals(result)
    verdict, pct = verdict_from_signals(signals)
    print(f"[investigate] VERDICT: {verdict}  ({pct}/100)")
    for s in signals:
        print(f"[investigate]   • {s.name:32s} {s.strength:9s} {s.rationale[:70]}")

    if save_results:
        results_root = project_root / "results"
        results_root.mkdir(exist_ok=True)
        label = label_override or f"investigate_{accession}"
        bundle_dir = make_bundle_dir(results_root, label)
        inv_dir = bundle_dir / "investigations" / accession
        paths = write_investigation(inv_dir, result)
        print(f"[investigate] bundle saved → {bundle_dir}")
        print(f"[investigate]   case file: {paths.get('case_file_md','')}")
        return {"bundle": str(bundle_dir), **paths}
    return {"verdict": verdict, "pct": pct}
