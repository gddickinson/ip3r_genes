"""Synthesize an investigation into a verdict + case file.

The verdict combines six independent signals (one per evidence line) into a
weighted score. Each signal returns a strength in {strong, moderate, weak,
absent} along with a one-line rationale. The renderer outputs a single
markdown "case file" suitable for sharing with collaborators.

Designed so a human can immediately see *which* evidence is driving the
verdict, not just the bottom-line score.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from ..utils.family import LITERATURE_KEYWORD
from .pipeline import InvestigationResult


Strength = Literal["strong", "moderate", "weak", "absent"]


@dataclass
class Signal:
    name: str
    strength: Strength
    rationale: str


def _score(s: Strength) -> int:
    return {"strong": 3, "moderate": 2, "weak": 1, "absent": 0}[s]


def assess_signals(result: InvestigationResult) -> list[Signal]:
    signals: list[Signal] = []
    family_kw = result.options.family_keywords

    # 1. Domain architecture
    arch = result.domain_arch
    if arch and arch.protein_length > 0:
        coverage = arch.find_family_coverage(family_kw)
        family_signs = [a for a in arch.annotations
                        if any(k in a.name.lower() for k in family_kw)]
        if coverage >= 0.6 and len(family_signs) >= 2:
            signals.append(Signal(
                "Domain architecture", "strong",
                f"{len(family_signs)} family signatures covering "
                f"{coverage:.0%} of the {arch.protein_length} aa protein "
                f"(e.g. {family_signs[0].name})",
            ))
        elif coverage >= 0.3:
            signals.append(Signal(
                "Domain architecture", "moderate",
                f"{len(family_signs)} family signatures covering "
                f"{coverage:.0%}",
            ))
        elif family_signs:
            signals.append(Signal(
                "Domain architecture", "weak",
                f"Only {len(family_signs)} family signature(s); coverage {coverage:.0%}",
            ))
        else:
            signals.append(Signal("Domain architecture", "absent",
                                   "No family domain signatures hit"))
    else:
        signals.append(Signal("Domain architecture", "absent",
                               "InterPro lookup failed or returned no entries"))

    # 2. UniProt evidence
    u = result.uniprot
    if u:
        channel_terms = u.channel_related_go_terms()
        if u.n_transmembrane >= 20 and channel_terms:
            signals.append(Signal(
                "UniProt evidence", "strong",
                f"{u.n_transmembrane} transmembrane regions annotated; "
                f"GO includes channel/ion-transport terms",
            ))
        elif u.n_transmembrane >= 10:
            signals.append(Signal(
                "UniProt evidence", "moderate",
                f"{u.n_transmembrane} transmembrane regions; "
                f"{len(channel_terms)} channel-related GO term(s)",
            ))
        elif u.n_transmembrane > 0:
            signals.append(Signal(
                "UniProt evidence", "weak",
                f"Only {u.n_transmembrane} TM region(s) annotated",
            ))
        else:
            signals.append(Signal(
                "UniProt evidence", "absent",
                "No transmembrane features in UniProt",
            ))
    else:
        signals.append(Signal("UniProt evidence", "absent", "UniProt fetch failed"))

    # 3. Predicted structure
    s = result.structure
    if s and s.alphafold_available:
        plddt_global = s.global_plddt or 0
        plddt_tm = s.plddt_in_tm
        if plddt_global >= 70 and (plddt_tm is None or plddt_tm >= 70):
            signals.append(Signal(
                "AlphaFold structure", "strong",
                f"Global pLDDT={plddt_global:.0f}"
                + (f"; mean pLDDT in TM regions={plddt_tm:.0f}" if plddt_tm else ""),
            ))
        elif plddt_global >= 60:
            signals.append(Signal(
                "AlphaFold structure", "moderate",
                f"Global pLDDT={plddt_global:.0f}",
            ))
        else:
            signals.append(Signal(
                "AlphaFold structure", "weak",
                f"Low confidence (pLDDT={plddt_global:.0f})",
            ))
    else:
        signals.append(Signal("AlphaFold structure", "absent",
                               "No AlphaFold model available"))

    # 4. Foldseek (only if it was run)
    if s and s.foldseek_hits:
        signals.append(Signal(
            "Structural homology (Foldseek)", "strong",
            f"{len(s.foldseek_hits)} structural hit(s) returned",
        ))
    elif result.options.run_foldseek:
        signals.append(Signal(
            "Structural homology (Foldseek)", "absent",
            "Foldseek run requested but no hits — could indicate sequence "
            "issues or genuinely novel fold",
        ))

    # 5. Phylogenetic placement
    p = result.phylo
    if p and p.nearest_panel_identity:
        ident = p.nearest_panel_identity
        if 0.20 <= ident <= 0.45:
            signals.append(Signal(
                "Phylogenetic placement", "strong",
                f"Nearest panel member {p.nearest_panel_label} at "
                f"{ident:.0%} identity — sits in the paralog twilight zone "
                f"({p.topology_class})",
            ))
        elif ident > 0.45:
            signals.append(Signal(
                "Phylogenetic placement", "moderate",
                f"{ident:.0%} identity to {p.nearest_panel_label} — possibly "
                f"an existing paralog or close ortholog ({p.topology_class})",
            ))
        else:
            signals.append(Signal(
                "Phylogenetic placement", "weak",
                f"Only {ident:.0%} identity to nearest panel member — "
                f"likely outside the family or extreme outlier",
            ))
    else:
        signals.append(Signal("Phylogenetic placement", "absent",
                               "Phylogenetic context not built"))

    # 6. Literature evidence
    lit = result.literature
    if lit:
        if lit.direct_pmids:
            signals.append(Signal(
                "Literature", "strong",
                f"{len(lit.direct_pmids)} PubMed result(s) cite this accession",
            ))
        elif lit.organism_family_pmids:
            signals.append(Signal(
                "Literature", "moderate",
                f"No direct hits, but {len(lit.organism_family_pmids)} paper(s) "
                f"on '{result.uniprot.organism if result.uniprot else 'organism'} "
                f"+ {LITERATURE_KEYWORD}'",
            ))
        else:
            signals.append(Signal("Literature", "absent",
                                   "No relevant PubMed results — uncharted protein"))

    # 7. Synteny — informational only (don't down-weight if missing for
    # organisms where EBI doesn't have coordinates)
    syn = result.synteny
    if syn and syn.flanking_genes:
        signals.append(Signal(
            "Genomic context", "moderate",
            f"Flanking gene set: {', '.join(syn.flanking_genes[:5])}",
        ))
    elif syn and syn.chromosome:
        signals.append(Signal(
            "Genomic context", "weak",
            f"Located on {syn.chromosome} but no flanking genes returned",
        ))

    return signals


def verdict_from_signals(signals: list[Signal]) -> tuple[str, int]:
    total = sum(_score(s.strength) for s in signals)
    max_total = 3 * len(signals)
    pct = 100 * total / max(1, max_total)
    if pct >= 70:
        return ("STRONG novel-paralog candidate — pursue with wet-lab follow-up", int(pct))
    if pct >= 50:
        return ("MODERATE — convergent evidence supports family membership, "
                "but at least one line is weak or absent", int(pct))
    if pct >= 30:
        return ("WEAK — some evidence consistent with the family, but several "
                "key lines fail to fire", int(pct))
    return ("UNLIKELY — does not look like a real family paralog", int(pct))


def render_case_file(result: InvestigationResult) -> str:
    signals = assess_signals(result)
    verdict, pct = verdict_from_signals(signals)

    L: list[str] = []
    u = result.uniprot
    arch = result.domain_arch
    L.append(f"# Investigation case file — `{result.accession}`")
    L.append("")
    if u:
        L.append(f"**Protein name:** {u.protein_name}")
        L.append(f"**Organism:** *{u.organism}* (taxon {u.taxon_id})")
        L.append(f"**Length:** {u.sequence_length} aa")
        L.append(f"**Review status:** {'reviewed (SwissProt)' if u.reviewed else 'unreviewed (TrEMBL)'}")
        L.append(f"**Evidence level:** {u.evidence_level}")
    L.append("")
    L.append(f"## Verdict — **{verdict}** ({pct}/100)")
    L.append("")
    L.append("| Evidence line | Strength | Rationale |")
    L.append("|---|---|---|")
    for s in signals:
        L.append(f"| {s.name} | **{s.strength}** | {s.rationale} |")
    L.append("")

    # Domain architecture detail
    if arch and arch.annotations:
        L.append("## 1. Domain architecture")
        L.append("")
        L.append(f"InterPro reports {len(arch.annotations)} entries across the {arch.protein_length} aa protein.")
        L.append("")
        L.append("| Source | Accession | Type | Name | Locations |")
        L.append("|---|---|---|---|---|")
        for a in arch.annotations[:25]:
            locs = ", ".join(f"{s}–{e}" for s, e in a.locations[:3])
            L.append(f"| {a.source} | {a.accession} | {a.type} | {a.name[:60]} | {locs} |")
        L.append("")

    # UniProt detail
    if u:
        L.append("## 2. UniProt detail")
        L.append("")
        L.append(f"- **Transmembrane regions annotated:** {u.n_transmembrane}")
        if u.transmembrane_regions:
            tm_summary = ", ".join(f"{s}–{e}" for s, e in u.transmembrane_regions[:8])
            more = f", … and {len(u.transmembrane_regions) - 8} more" if len(u.transmembrane_regions) > 8 else ""
            L.append(f"  - First TM regions: {tm_summary}{more}")
        if u.channel_related_go_terms():
            L.append("- **Channel-related GO terms:**")
            for gid, name, asp in u.channel_related_go_terms()[:10]:
                L.append(f"  - {gid} *(\"{name}\", {asp})*")
        if u.go_terms:
            L.append(f"- **All GO terms ({len(u.go_terms)}):** "
                     + ", ".join(g[0] for g in u.go_terms[:10])
                     + (" …" if len(u.go_terms) > 10 else ""))
        if u.pdb_ids:
            L.append(f"- **PDB entries:** {', '.join(u.pdb_ids[:5])}")
        if u.xref_counts:
            xr = ", ".join(f"{k}({v})" for k, v in sorted(u.xref_counts.items(), key=lambda x: -x[1])[:10])
            L.append(f"- **Cross-references:** {xr}")
        L.append("")

    # Structure
    s = result.structure
    if s:
        L.append("## 3. Predicted 3D structure")
        L.append("")
        if s.alphafold_available:
            L.append(f"- AlphaFold model available: `{s.alphafold_url}`")
            if s.global_plddt is not None:
                L.append(f"- Global mean pLDDT: **{s.global_plddt:.1f}** "
                         f"({'high confidence' if s.global_plddt >= 70 else 'moderate' if s.global_plddt >= 60 else 'low'})")
            if s.plddt_in_tm is not None:
                L.append(f"- Mean pLDDT across TM regions: **{s.plddt_in_tm:.1f}**")
            if s.pae_image_url:
                L.append(f"- PAE plot: `{s.pae_image_url}`")
        else:
            L.append("- No AlphaFold model in the public database — "
                     "structure not yet predicted at scale for this entry.")
            L.append("  → Recommended: submit the sequence to ColabFold / AlphaFold Server.")
        if s.foldseek_hits:
            L.append("")
            L.append(f"- Foldseek structural hits: {len(s.foldseek_hits)}")
            for h in s.foldseek_hits[:5]:
                L.append(f"  - {h.get('target','?')} — {h.get('description','')[:80]}")
        L.append("")

    # Phylogeny
    p = result.phylo
    if p:
        L.append("## 4. Phylogenetic placement")
        L.append("")
        if p.nearest_panel_label:
            L.append(f"- Nearest panel member: **{p.nearest_panel_label}** "
                     f"({p.nearest_panel_identity:.0%} identity)")
            L.append(f"- Topology class: `{p.topology_class}`")
        if p.pairwise_identities:
            L.append("")
            L.append("| Panel member | Identity |")
            L.append("|---|---|")
            for k, v in sorted(p.pairwise_identities.items(), key=lambda x: -x[1]):
                L.append(f"| {k} | {v:.0%} |")
        if p.ascii_tree:
            L.append("")
            L.append("```")
            L.append(p.ascii_tree.rstrip())
            L.append("```")
        L.append("")

    # Synteny
    syn = result.synteny
    if syn:
        L.append("## 5. Genomic context (synteny lite)")
        L.append("")
        if syn.genome_assembly or syn.chromosome:
            L.append(f"- Genome assembly: {syn.genome_assembly or '(unknown)'}")
            if syn.chromosome:
                L.append(f"- Chromosome / contig: {syn.chromosome}")
            if syn.location:
                L.append(f"- Coordinates: {syn.location}")
        if syn.flanking_genes:
            L.append(f"- Flanking genes: {', '.join(syn.flanking_genes)}")
        if syn.notes:
            L.append(f"- _Note:_ {syn.notes}")
        L.append("")

    # Literature
    lit = result.literature
    if lit:
        L.append("## 6. Literature")
        L.append("")
        if lit.direct_titles:
            L.append("**Direct accession hits:**")
            for pmid, title in lit.direct_titles:
                L.append(f"- PubMed [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) — {title}")
        elif lit.direct_pmids:
            L.append(f"**{len(lit.direct_pmids)} direct hit(s)** — see `investigation.json` for the full list.")
        else:
            L.append("- No PubMed entries cite this accession.")
        L.append("")
        if lit.organism_family_titles:
            L.append(f"**Organism + '{LITERATURE_KEYWORD}' hits:**")
            for pmid, title in lit.organism_family_titles:
                L.append(f"- PubMed [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/) — {title}")
        L.append("")

    # Suggested follow-up
    L.append("## 7. Suggested follow-up experiments")
    L.append("")
    follow_up = _suggest_followup(result, signals)
    for item in follow_up:
        L.append(f"- {item}")
    L.append("")
    L.append("---")
    L.append("_Generated by the Protein Variant Finder investigation pipeline._")
    return "\n".join(L) + "\n"


def _suggest_followup(result: InvestigationResult, signals: list[Signal]) -> list[str]:
    out: list[str] = []
    s = result.structure
    arch = result.domain_arch
    p = result.phylo

    if s and not s.alphafold_available:
        out.append("Predict the 3D structure with ColabFold / AlphaFold Server "
                   "and verify the canonical IP3R fold: a beta-trefoil IP3-binding\n"
                   "core over a six-TM pore, tetrameric.")
    if s and s.alphafold_available and not result.options.run_foldseek:
        out.append("Run Foldseek with this protein's structure as bait against "
                   "the AlphaFold DB to find structural homologs that BLAST misses.")
    if arch and arch.has_family_signature(result.options.family_keywords) and arch.find_family_coverage(result.options.family_keywords) >= 0.5:
        out.append("Re-annotate this entry — propose a name change to UniProt "
                   "since the domain architecture is dominated by the family "
                   "signature, not the existing 'recommended name'.")
    if p and 0.20 <= p.nearest_panel_identity <= 0.45:
        out.append(f"Heterologously express in HEK293T and assay calcium influx "
                   f"after IP3 uncaging or agonist stimulation — the standard "
                   f"functional test for an IP3-gated release channel.")
    if any(s.strength in ("strong", "moderate") for s in signals if s.name == "Domain architecture"):
        out.append("Search the genome assembly for additional family-member paralogs — "
                   "if this is a real paralog, others may exist in the same lineage.")
    out.append("Cross-check with OrthoDB and OMA browser entries (see UniProt xrefs) "
               "to see whether existing orthology databases group it correctly.")
    return out
