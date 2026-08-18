"""Domain-bait scan — enumerate every protein with a family Pfam signature
that is *not* already named as a known paralog. The remainder is the
candidate pool for the discovery scorer.

This is the v1.4 "fourth-paralog hunt" path. The standard preset search
(NCBI/Ensembl/UniProt by gene symbol) only returns proteins already named
ITPR1/2/3 — by definition no novel paralogs. The domain-bait route flips
the query: ask InterPro "give me every protein with PF15917" and filter
out the things we already know.

Public entry point: `run_domain_scan(...)`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from ..analysis import analyse
from ..analysis.pipeline import write_analysis, _label_for
from ..core.models import ProteinVariant, SearchQuery, SearchResult, SearchStatus
from ..databases.interpro import (
    ProteinWithDomain,
    fetch_uniprot_sequence,
    list_proteins_with_pfam,
)
from ..discovery.candidates import DiscoveryConfig, discover_novel_paralogs, write_discovery
from ..utils.family import KNOWN_PARALOGS
from ..utils.report import write_report
from ..utils.results_writer import make_bundle_dir, write_bundle


def run_domain_scan(
    project_root: Path,
    preset_data: dict,
    preset_name: str,
    email: str = "",
    save_results: bool = False,
    label_override: str = "",
    analyze: bool = True,
    discover: bool = True,
    log=print,
) -> dict:
    pfam_ids = list(preset_data.get("pfam_ids", []))
    exclude_name = [s.lower() for s in preset_data.get("exclude_name_substrings", [])]
    exclude_gene = [s.lower() for s in preset_data.get("exclude_gene_substrings", [])]
    known = list(preset_data.get("known_paralogs", []))
    max_per_pfam = int(preset_data.get("max_proteins_per_pfam", 400))
    max_with_seq = int(preset_data.get("max_with_sequences", 40))
    reviewed_only = bool(preset_data.get("reviewed_only", False))
    min_len = int(preset_data.get("min_length_aa", 0))
    max_len = int(preset_data.get("max_length_aa", 10**9))

    if not pfam_ids:
        log("[domain-scan] preset is missing 'pfam_ids' — nothing to do")
        return {}

    log(f"[domain-scan] preset={preset_name}  pfam_ids={pfam_ids}  "
          f"exclude_name={exclude_name}  reviewed_only={reviewed_only}")

    # 1. Enumerate via InterPro
    rows: list[ProteinWithDomain] = []
    for pfam in pfam_ids:
        log(f"[domain-scan] listing proteins with {pfam} (cap {max_per_pfam})…")
        batch = list_proteins_with_pfam(pfam, max_results=max_per_pfam)
        log(f"[domain-scan]   {len(batch)} entries")
        rows.extend(batch)

    # 2. Filter
    filtered: list[ProteinWithDomain] = []
    skipped_known = skipped_size = skipped_reviewed = 0
    seen_accessions: set[str] = set()
    for row in rows:
        if row.accession in seen_accessions:
            continue
        seen_accessions.add(row.accession)
        name_lc = row.name.lower()
        gene_lc = row.gene.lower()
        if any(s in name_lc for s in exclude_name) or any(s in gene_lc for s in exclude_gene):
            skipped_known += 1
            continue
        if reviewed_only and not row.reviewed:
            skipped_reviewed += 1
            continue
        if row.length and (row.length < min_len or row.length > max_len):
            skipped_size += 1
            continue
        filtered.append(row)
    log(f"[domain-scan] after filtering: {len(filtered)} candidate(s)  "
          f"(excluded known={skipped_known}, size={skipped_size}, reviewed-only={skipped_reviewed})")

    if not filtered:
        log("[domain-scan] no candidates — try relaxing size band, reviewed-only, or excludes")
        return {}

    # 3. Score-rank the candidate pool by simple heuristics so we pick the
    # right N for sequence fetching: prefer reviewed entries, then larger
    # proteins, then ones with AlphaFold predictions.
    filtered.sort(
        key=lambda r: (-int(r.reviewed), -int(r.in_alphafold), -(r.length or 0)),
    )
    to_seqfetch = filtered[:max_with_seq]

    # 4. Build ProteinVariants. Fetch sequences for the top N so the
    # downstream MSA / tree have something to chew on. Other candidates
    # still appear in the result set, just without sequences.
    variants: list[ProteinVariant] = []
    log(f"[domain-scan] fetching sequences for top {len(to_seqfetch)} candidate(s)…")
    fetched = 0
    for i, row in enumerate(filtered):
        seq = ""
        if i < max_with_seq:
            seq = fetch_uniprot_sequence(row.accession)
            if seq:
                fetched += 1
        # Use a synthetic gene_symbol so the discovery scorer sees these as
        # "unknown family members" — pick something stable and meaningful.
        candidate_label = row.gene if row.gene else f"unnamed_{row.pfam_id}"
        variants.append(ProteinVariant(
            source="InterPro",
            accession=row.accession,
            gene_symbol=candidate_label,
            species=row.species,
            taxon_id=row.taxon_id,
            length_aa=row.length or None,
            description=f"[{row.pfam_id}] {row.name}"
                       + (" [reviewed]" if row.reviewed else " [unreviewed]"),
            gene_id=row.accession,
            sequence=seq,
            url=f"https://www.uniprot.org/uniprotkb/{row.accession}",
            raw={
                "pfam_id": row.pfam_id,
                "reviewed": str(row.reviewed),
                "in_alphafold": str(row.in_alphafold),
                "interpro_name": row.name,
                "interpro_gene": row.gene,
            },
        ))
    log(f"[domain-scan]   fetched {fetched} sequence(s)")

    # 5. Run analysis on the sequence-bearing subset
    analysis_result = None
    if analyze:
        log("[domain-scan] running sequence analysis on the candidate pool…")
        analysis_result = analyse(variants, identity_threshold=0.4, max_variants=60)
        log(f"[domain-scan]   analysed {analysis_result.n_analyzed}/{analysis_result.n_input}")

    # 6. Score each candidate with the discovery pipeline. domain_hits is
    # pre-populated from the InterPro list itself (we know every row carries
    # at least one Pfam ID from the bait set).
    discovery_report = None
    discovery_text = ""
    if discover:
        log("[domain-scan] scoring candidates…")
        # Build domain_hits dict so the discovery scorer gets credit for the
        # Pfam signature without a second InterPro round-trip.
        from ..databases.interpro import DomainHit
        domain_hits = {
            v.accession: [DomainHit(
                accession=v.accession,
                pfam_id=v.raw.get("pfam_id", ""),
                pfam_name=v.raw.get("interpro_name", ""),
                n_hits=1,
            )]
            for v in variants
        }
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v) for v in variants}
        discovery_report = discover_novel_paralogs(
            variants,
            DiscoveryConfig(known_paralogs=known or list(KNOWN_PARALOGS)),
            analysis_label_for=label_map,
            distances=(analysis_result.distances if analysis_result else None),
            domain_hits=domain_hits,
            foldseek_hits=None,
        )
        discovery_text = discovery_report.text_summary()
        log(f"[domain-scan]   promising:    "
              f"{sum(1 for c in discovery_report.candidates if c.score >= 60)}")
        log(f"[domain-scan]   worth review: "
              f"{sum(1 for c in discovery_report.candidates if 40 <= c.score < 60)}")
        log(f"[domain-scan]   weak:         "
              f"{sum(1 for c in discovery_report.candidates if c.score < 40)}")

    # 7. Package the outcome. The pseudo query/results make bundle metadata
    # sensible and let the GUI adopt the scan as its current result set.
    pseudo_query = SearchQuery(
        gene_symbols=pfam_ids,
        species="",
        sources=["InterPro"],
        max_results_per_source=max_per_pfam,
        include_sequence=True,
    )
    pseudo_results = [SearchResult(
        source="InterPro",
        query=pseudo_query,
        status=SearchStatus.OK,
        variants=variants,
        message=f"domain-scan: kept {len(filtered)} of {len(rows)} after filtering",
        elapsed_s=0.0,
    )]
    out = {
        "n_variants": len(variants),
        "variants": variants,
        "query": pseudo_query,
        "results": pseudo_results,
        "analysis": analysis_result,
        "discovery_report": discovery_report,
        "discovery_text": discovery_text,
    }

    # 8. Persist the bundle if requested
    if save_results:
        results_root = project_root / "results"
        results_root.mkdir(exist_ok=True)
        label = label_override or "domain_scan_" + "_".join(pfam_ids[:2]).lower()
        bundle_dir = make_bundle_dir(results_root, label)
        summary = write_bundle(
            bundle_dir, variants, pseudo_query, pseudo_results,
            notes=f"Domain-bait scan from preset '{preset_name}'.",
        )
        if analysis_result:
            write_analysis(bundle_dir / "analysis", analysis_result)
        if discovery_report is not None:
            write_discovery(bundle_dir / "discovery", discovery_report)
        report_paths = write_report(
            bundle_dir, variants, pseudo_query,
            analysis=analysis_result,
            discovery_summary=discovery_text,
            notes=f"Domain-bait scan (fourth-paralog hunt) — preset '{preset_name}'.",
        )
        log(f"[domain-scan] bundle saved → {summary['dir']}")
        log(f"[domain-scan]   report.html: {report_paths['report_html']}")
        out.update(summary)
        out.update(report_paths)
    return out
