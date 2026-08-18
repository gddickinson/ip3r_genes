"""Investigator — orchestrates the seven lines of evidence for one candidate.

Each line is independent and degrades gracefully when its source isn't
available (e.g., AlphaFold returns 404 for many TrEMBL entries — the
pipeline records the gap and moves on; the verdict logic accounts for it).

Public API:
    inv = Investigator(accession, options).run()
    write_investigation(out_dir, inv)
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional

from .domain_arch import DomainArchitecture, fetch_domain_architecture
from .literature import LiteratureEvidence, search_literature
from .phylo_context import PhyloContext, build_phylo_context
from .structure import StructuralEvidence, gather_structural_evidence
from .synteny_lite import SyntenyContext, fetch_synteny_context
from ..utils.family import LITERATURE_KEYWORD, REFERENCE_PANEL
from .uniprot_detail import UniProtDetail, fetch_uniprot_detail


# Default phylogenetic comparison panel: the three human paralogs, declared
# once in `src/utils/family.py`. Override with the `panel` option — e.g. add
# a ryanodine receptor when the question is which side of the superfamily
# split a candidate falls on.
DEFAULT_PHYLO_PANEL: list[tuple[str, str]] = list(REFERENCE_PANEL)

DEFAULT_FAMILY_KEYWORDS = [LITERATURE_KEYWORD]


@dataclass
class InvestigationOptions:
    panel: list[tuple[str, str]] = field(default_factory=lambda: list(DEFAULT_PHYLO_PANEL))
    family_keywords: list[str] = field(default_factory=lambda: list(DEFAULT_FAMILY_KEYWORDS))
    run_foldseek: bool = False
    run_literature: bool = True
    run_synteny: bool = True
    run_phylo: bool = True
    email: str = ""
    timeout_s: int = 30


@dataclass
class InvestigationResult:
    accession: str
    options: InvestigationOptions
    domain_arch: Optional[DomainArchitecture] = None
    uniprot: Optional[UniProtDetail] = None
    structure: Optional[StructuralEvidence] = None
    phylo: Optional[PhyloContext] = None
    synteny: Optional[SyntenyContext] = None
    literature: Optional[LiteratureEvidence] = None
    notes: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "accession": self.accession,
            "options": {
                "panel": self.options.panel,
                "family_keywords": self.options.family_keywords,
                "run_foldseek": self.options.run_foldseek,
            },
            "domain_arch": (asdict(self.domain_arch) if self.domain_arch else None),
            "uniprot": (asdict(self.uniprot) if self.uniprot else None),
            "structure": (asdict(self.structure) if self.structure else None),
            "phylo": (asdict(self.phylo) if self.phylo else None),
            "synteny": (asdict(self.synteny) if self.synteny else None),
            "literature": (asdict(self.literature) if self.literature else None),
            "notes": self.notes,
        }


class Investigator:
    def __init__(
        self,
        accession: str,
        options: Optional[InvestigationOptions] = None,
        candidate_sequence: str = "",
    ) -> None:
        self.accession = accession
        self.options = options or InvestigationOptions()
        self.candidate_sequence = candidate_sequence

    def run(self, on_progress=None) -> InvestigationResult:
        result = InvestigationResult(accession=self.accession, options=self.options)
        def tick(label: str) -> None:
            if on_progress:
                on_progress(label)

        tick("Fetching domain architecture (InterPro)…")
        result.domain_arch = fetch_domain_architecture(self.accession, self.options.timeout_s)

        tick("Fetching UniProt detail…")
        result.uniprot = fetch_uniprot_detail(self.accession, self.options.timeout_s)

        if result.domain_arch and result.uniprot:
            result.domain_arch.n_transmembrane = result.uniprot.n_transmembrane
            result.domain_arch.family_signatures = [
                a.accession for a in result.domain_arch.annotations
                if any(k in a.name.lower() for k in self.options.family_keywords)
            ]

        tick("Fetching AlphaFold prediction…")
        tm_regions = result.uniprot.transmembrane_regions if result.uniprot else None
        if not self.candidate_sequence and result.uniprot and result.uniprot.sequence_length:
            # Sequence not yet known — phylogenetic line will fetch it.
            pass
        result.structure = gather_structural_evidence(
            self.accession,
            transmembrane_regions=tm_regions,
            foldseek_bait_sequence=self.candidate_sequence,
            run_foldseek=self.options.run_foldseek,
            timeout_s=self.options.timeout_s,
        )

        if self.options.run_phylo:
            tick("Building phylogenetic context (panel MSA + NJ tree)…")
            result.phylo = build_phylo_context(
                self.accession,
                self.candidate_sequence,
                self.options.panel,
                timeout_s=self.options.timeout_s,
            )

        if self.options.run_synteny:
            tick("Fetching genomic coordinates (synteny lite)…")
            result.synteny = fetch_synteny_context(self.accession, self.options.timeout_s)

        if self.options.run_literature:
            tick("Searching PubMed…")
            organism = result.uniprot.organism if result.uniprot else ""
            uniprot_id = result.uniprot.uniprot_id if result.uniprot else ""
            result.literature = search_literature(
                self.accession,
                uniprot_id=uniprot_id,
                organism=organism,
                email=self.options.email,
                timeout_s=self.options.timeout_s,
            )

        tick("Done.")
        return result


def write_investigation(out_dir: Path, result: InvestigationResult) -> dict:
    """Persist the investigation. Returns paths written."""
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, str] = {}

    raw_path = out_dir / "investigation.json"
    raw_path.write_text(json.dumps(result.to_dict(), indent=2, default=str))
    paths["json"] = str(raw_path)

    # Phylo Newick if available
    if result.phylo and result.phylo.newick:
        nwk = out_dir / "phylo_context.newick"
        nwk.write_text(result.phylo.newick + "\n")
        paths["phylo_newick"] = str(nwk)
        ascii_path = out_dir / "phylo_context_ascii.txt"
        ascii_path.write_text(result.phylo.ascii_tree)
        paths["phylo_ascii"] = str(ascii_path)

    # Case file
    from .synthesis import render_case_file
    md_path = out_dir / "case_file.md"
    md_text = render_case_file(result)
    md_path.write_text(md_text)
    paths["case_file_md"] = str(md_path)

    return paths
