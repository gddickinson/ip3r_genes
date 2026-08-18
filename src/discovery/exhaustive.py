"""Exhaustive hunt — every ITPR-like sequence, every analysis, one run.

The all-in-one pipeline behind the `"mode": "exhaustive"` preset. It
combines the three harvest routes that each found something the others
missed, then runs the full analysis battery:

    Harvest
    1. Compara paralog mine       — every within-species paralog of the
                                    known genes across the vertebrate panel
                                    (name-agnostic; finds unnamed models).
    2. Ortholog expansion         — for every NON-known gene the mine
                                    surfaced, pull its own Compara
                                    orthologs across all vertebrates (this
                                    is how the opossum/tortoise/cichlid
                                    family set was found), keep the best
                                    per species, fetch sequences.
    3. InterPro domain scan       — every UniProt protein carrying the
                                    family Pfam signature, minus known
                                    names (catches proteins Ensembl lacks).

    Merge → dedupe → save the census (FASTA + CSV) → analyse (MSA,
    identity matrix, NJ tree) → discovery scoring (with MSA-signature
    fallback) → presence/absence retention matrix → highlighted tree
    figure → bundle + report.

Returns the same shape as `run_domain_scan` (variants/query/results/
analysis/discovery_report…) so the GUI adopts the run as its result set,
plus `presence_report`, `sequence_fasta`, `sequence_csv`, `tree_png`.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Optional

import requests

from ..analysis import analyse
from ..analysis.figures import HAVE_MPL, tree_figure
from ..analysis.pipeline import _label_for, write_analysis
from ..analysis.presence import build_presence_matrix
from ..analysis.tree import build_nj_tree
from ..core.models import ProteinVariant, SearchQuery, SearchResult, SearchStatus
from ..databases.compara import ComparaClient, ENSEMBL_REST, _PAUSE_S
from ..databases.interpro import fetch_uniprot_sequence, list_proteins_with_pfam
from ..utils.family import KNOWN_PARALOGS
from ..discovery.candidates import (
    DiscoveryConfig,
    _is_known,
    build_signature_set,
    discover_novel_paralogs,
    write_discovery,
)
from ..utils.report import write_report
from ..utils.results_writer import make_bundle_dir, write_bundle


def _get_json(path: str, timeout_s: int = 30, **params):
    params.setdefault("content-type", "application/json")
    r = requests.get(f"{ENSEMBL_REST}{path}", params=params, timeout=timeout_s)
    time.sleep(_PAUSE_S)
    if r.status_code in (400, 404):
        return None
    r.raise_for_status()
    return r.json()


def expand_orthologs(
    seed: ProteinVariant,
    max_per_gene: int = 12,
    fetch_sequences: bool = True,
    log=print,
) -> list[ProteinVariant]:
    """All-vertebrate Compara orthologs of one seed gene, best per species."""
    gene_id, species = seed.gene_id, seed.species.lower().replace(" ", "_")
    data = _get_json(f"/homology/id/{species}/{gene_id}",
                     type="orthologues", sequence="none")
    if not data or not data.get("data"):
        return []
    best_per_species: dict[str, dict] = {}
    for entry in data["data"]:
        for h in entry.get("homologies", []):
            t = h.get("target", {})
            sp, pid = t.get("species", ""), t.get("protein_id", "")
            if not sp or not pid:
                continue
            one2one = h.get("type") == "ortholog_one2one"
            score = (int(one2one), t.get("perc_id") or 0.0)
            prev = best_per_species.get(sp)
            if prev is None or score > prev["_score"]:
                best_per_species[sp] = {**t, "_score": score, "_type": h.get("type", "")}
    ranked = sorted(best_per_species.values(),
                    key=lambda t: t["_score"], reverse=True)[:max_per_gene]
    out: list[ProteinVariant] = []
    for t in ranked:
        sp_slug = t["species"]
        v = ProteinVariant(
            source="Compara",
            accession=t["protein_id"],
            gene_symbol=f"{seed.gene_symbol}_ortholog",
            species=sp_slug.replace("_", " ").capitalize(),
            length_aa=None,
            description=(f"ortholog of {seed.gene_symbol} ({seed.species}) — "
                         f"{t.get('_type','')}, {t.get('perc_id') or 0:.0f}% id"),
            gene_id=t.get("id", ""),
            url=f"https://www.ensembl.org/{sp_slug.capitalize()}/Gene/Summary?g={t.get('id','')}",
            raw={"protein_id": t["protein_id"], "homology_type": t.get("_type", ""),
                 "anchor": seed.gene_symbol, "expanded_from": seed.accession},
        )
        if fetch_sequences:
            seq_data = _get_json(f"/sequence/id/{t['protein_id']}", type="protein")
            if seq_data and seq_data.get("seq"):
                v.sequence = seq_data["seq"]
                v.length_aa = len(v.sequence)
        out.append(v)
    log(f"[exhaustive]   {seed.gene_symbol}: +{len(out)} ortholog(s) across vertebrates")
    return out


def save_census(out_dir: Path, variants: list[ProteinVariant]) -> dict[str, str]:
    """The deliverable list: every ITPR-like sequence found, FASTA + CSV."""
    out_dir.mkdir(parents=True, exist_ok=True)
    fasta = out_dir / "itpr_like_sequences.fasta"
    with fasta.open("w") as f:
        for v in variants:
            if v.sequence:
                f.write(f">{v.source}|{v.accession}|{v.gene_symbol}|{v.species.replace(' ', '_')}\n")
                for i in range(0, len(v.sequence), 60):
                    f.write(v.sequence[i:i + 60] + "\n")
    csv = out_dir / "itpr_like_census.csv"
    with csv.open("w") as f:
        f.write("source,accession,gene_symbol,species,length_aa,has_sequence,description\n")
        for v in variants:
            desc = v.description.replace('"', "'")
            f.write(f'{v.source},{v.accession},{v.gene_symbol},"{v.species}",'
                    f'{v.length_aa or ""},{bool(v.sequence)},"{desc}"\n')
    return {"sequence_fasta": str(fasta), "sequence_csv": str(csv)}


def run_exhaustive_hunt(
    project_root: Path,
    preset_data: dict,
    preset_name: str = "ip3r_all",
    email: str = "",
    save_results: bool = True,
    label_override: str = "",
    log=print,
) -> dict:
    known = list(preset_data.get("known_paralogs", list(KNOWN_PARALOGS)))
    genes = list(preset_data.get("genes", known))
    pfam_ids = list(preset_data.get("pfam_ids", []))
    max_per_source = int(preset_data.get("max_per_source", 60))
    max_orthologs = int(preset_data.get("max_orthologs_per_candidate", 12))
    max_pfam = int(preset_data.get("max_proteins_per_pfam", 200))
    min_len = int(preset_data.get("min_length_aa", 1500))
    max_len = int(preset_data.get("max_length_aa", 3500))
    max_analyse = int(preset_data.get("max_analysis_variants", 80))
    known_upper = {g.upper() for g in known}

    # ---- 1. Compara paralog mine ----------------------------------------
    log(f"[exhaustive] 1/5  Compara paralog mine: {genes} across vertebrate panel…")
    query = SearchQuery(gene_symbols=genes, species="", sources=["Compara"],
                        max_results_per_source=max_per_source, include_sequence=True)
    mined = ComparaClient(email=email, timeout_s=45).search(query)
    log(f"[exhaustive]   {len(mined)} paralog record(s)")

    # ---- 2. Ortholog expansion of every non-known gene -------------------
    seeds = [v for v in mined if not _is_known(v.gene_symbol, known_upper)]
    log(f"[exhaustive] 2/5  expanding orthologs of {len(seeds)} novel/unnamed gene(s)…")
    expanded: list[ProteinVariant] = []
    for seed in seeds:
        try:
            expanded.extend(expand_orthologs(seed, max_per_gene=max_orthologs, log=log))
        except requests.RequestException as e:
            log(f"[exhaustive]   {seed.gene_symbol}: expansion failed ({e})")

    # ---- 3. InterPro domain enumeration ----------------------------------
    interpro_variants: list[ProteinVariant] = []
    if pfam_ids:
        log(f"[exhaustive] 3/5  InterPro enumeration for {pfam_ids}…")
        seen_acc: set[str] = set()
        rows = []
        for pfam in pfam_ids:
            batch = list_proteins_with_pfam(pfam, max_results=max_pfam)
            log(f"[exhaustive]   {pfam}: {len(batch)} entries")
            rows.extend(batch)
        exclude = [k.lower() for k in known]
        kept = 0
        for row in rows:
            if row.accession in seen_acc:
                continue
            seen_acc.add(row.accession)
            blob = f"{row.name} {row.gene}".lower()
            if any(k in blob for k in exclude):
                continue
            if row.length and not (min_len <= row.length <= max_len):
                continue
            seq = fetch_uniprot_sequence(row.accession) if kept < 25 else ""
            interpro_variants.append(ProteinVariant(
                source="InterPro",
                accession=row.accession,
                gene_symbol=row.gene or f"unnamed_{row.pfam_id}",
                species=row.species,
                taxon_id=row.taxon_id,
                length_aa=row.length or None,
                description=f"[{row.pfam_id}] {row.name}"
                           + (" [reviewed]" if row.reviewed else " [unreviewed]"),
                gene_id=row.accession,
                sequence=seq,
                url=f"https://www.uniprot.org/uniprotkb/{row.accession}",
                raw={"pfam_id": row.pfam_id, "reviewed": str(row.reviewed)},
            ))
            kept += 1
        log(f"[exhaustive]   kept {kept} non-known InterPro candidate(s)")

    # ---- merge + dedupe --------------------------------------------------
    merged: dict[str, ProteinVariant] = {}
    for v in mined + expanded + interpro_variants:
        key = v.accession.split(".")[0]
        prev = merged.get(key)
        if prev is None or (v.sequence and not prev.sequence):
            merged[key] = v
    variants = list(merged.values())
    n_seq = sum(1 for v in variants if v.sequence)
    log(f"[exhaustive] merged census: {len(variants)} unique protein(s), {n_seq} with sequences")

    # ---- 4. analyses -----------------------------------------------------
    log(f"[exhaustive] 4/5  MSA + tree + discovery + presence matrix…")
    analysis = analyse(variants, identity_threshold=0.4, max_variants=max_analyse)
    label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v) for v in variants}
    sset = build_signature_set(variants, analysis.msa, label_map, known)
    if sset:
        log(f"[exhaustive]   {len(sset.signatures)} family-signature blocks derived")
    discovery = discover_novel_paralogs(
        variants,
        DiscoveryConfig(known_paralogs=known, **preset_data.get("discovery", {})),
        analysis_label_for=label_map,
        distances=analysis.distances,
        signature_set=sset,
    )
    bins = discovery.by_verdict()
    log(f"[exhaustive]   discovery: {len(bins['promising'])} promising, "
        f"{len(bins['worth manual review'])} worth review")
    presence = build_presence_matrix(variants, distances=analysis.distances,
                                     analysis_label_for=label_map)

    # Labels to highlight in the tree = promising + worth-review candidates.
    novel_labels = {
        label_map.get(c.label, "") for c in discovery.candidates if c.score >= 40
    } - {""}

    out: dict = {
        "n_variants": len(variants), "variants": variants,
        "query": query, "analysis": analysis, "discovery_report": discovery,
        "discovery_text": discovery.text_summary(),
        "presence_report": presence, "novel_labels": novel_labels,
        "results": [SearchResult(
            source="Exhaustive", query=query, status=SearchStatus.OK,
            variants=variants,
            message=f"census: {len(mined)} mined + {len(expanded)} expanded "
                    f"+ {len(interpro_variants)} InterPro",
            elapsed_s=0.0)],
    }

    # ---- 5. persist ------------------------------------------------------
    if save_results:
        log(f"[exhaustive] 5/5  saving bundle…")
        results_root = project_root / "results"
        results_root.mkdir(exist_ok=True)
        bundle_dir = make_bundle_dir(results_root, label_override or "exhaustive_hunt")
        summary = write_bundle(bundle_dir, variants, query, out["results"],
                               notes=f"Exhaustive hunt from preset '{preset_name}'.")
        out.update(summary)
        out.update(save_census(bundle_dir, variants))
        write_analysis(bundle_dir / "analysis", analysis)
        write_discovery(bundle_dir / "discovery", discovery)
        (bundle_dir / "presence_matrix.txt").write_text(presence.summary())
        if HAVE_MPL and analysis.distances:
            tree = build_nj_tree(analysis.distances)
            png = tree_figure(tree, bundle_dir / "tree_highlighted.png",
                              title="ITPR family — NJ tree",
                              highlight_labels=novel_labels)
            if png:
                out["tree_png"] = str(png)
        report_paths = write_report(
            bundle_dir, variants, query, analysis=analysis,
            discovery_summary=discovery.text_summary() + "\n\n" + presence.summary(),
            notes=f"Exhaustive hunt (preset '{preset_name}').")
        out.update(report_paths)
        log(f"[exhaustive] bundle → {summary['dir']}")
    return out
