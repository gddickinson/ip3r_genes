"""Score and rank candidate novel paralogs (e.g. a fourth ITPR).

Finding an unnamed paralog is a manual integration of several
weak signals. This module automates that integration: each candidate gets
a transparent additive score in [0, 100] from up to six lines of evidence,
so the user can see *why* it's flagged and where the evidence is strong
or weak.

Inputs
------
- The full set of ProteinVariants returned by a search.
- A list of *known paralog* gene symbols (e.g. ["ITPR1","ITPR2","ITPR3"]).
- Optional InterPro Pfam domain hits per UniProt accession.
- Optional Foldseek structural hits already in `variants` (source=='Foldseek').

Scoring (component points; total capped at 100)
-----------------------------------------------
+20  Sequence outlier:    identity to the nearest KNOWN paralog sits in
                          [floor, ceiling] — the paralog twilight zone.
                          (Measured against *known* members so a clade of
                          mutually-similar novel orthologs can't mask
                          itself.)
+20  Structural fold:     a Foldseek hit links the candidate to a known
                          paralog at TM-score >= 0.5
+20  Domain signature:    candidate carries one of the family's defining
                          Pfam IDs (`family.FAMILY_PFAM_IDS`); OR (+15,
                          fallback for candidates absent from UniProt/
                          InterPro) it hits ≥ `signature_min_coverage` of
                          the MSA-derived family signature blocks
+15  Size plausibility:   length_aa within configured min/max (IP3Rs:
                          2000-3000 aa). Filters out fragments + huge fusions
+15  Taxonomic breadth:   present in >=2 species — same symbol OR corroborated
                          by *other candidates* in different species at
                          >= `breadth_identity_min` (unnamed ortholog sets
                          have a different id per species; mutual similarity
                          between candidates is the retention signal)
+10  Cluster exclusion:   does not co-cluster with any known paralog at
                          the family threshold
+10  Homology provenance: surfaced by a family-linked homology search
                          (Compara paralogy / BLAST bait) rather than name
+10  Split-annotation:    sits suspiciously close to another family locus —
                          the fragmented-gene-model pattern that hides a
                          real paralog behind two partial models

A candidate that scores >=60 is reported as "promising"; >=40 "worth
manual review"; <40 is shown for completeness but not promoted.

Promotion gate (S2 decision, benchmark-driven): a candidate can only be
promoted (score >= 40) if at least ONE family-specific evidence component
fired — twilight-zone outlier, domain signature (Pfam or MSA-signature
fallback), structural fold, or split-annotation. The generic components
(size, breadth, cluster, homology provenance) describe *any* large widely
conserved protein: the S1 negative controls showed a cross-species decoy
pair (TLN1) reaching exactly 40 on size+breadth+cluster alone. Candidates
failing the gate are capped at 39 with an explanatory evidence note
(`require_evidence_gate=False` restores the raw additive score).

Outputs
-------
- A `DiscoveryReport` with the ranked candidate list and per-candidate
  evidence breakdown.
- Written to `discovery/candidates.tsv` and `discovery/report.md` in the
  results bundle.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..analysis.distance import DistanceTable
from ..analysis.motifs import FamilySignatureSet, derive_family_signatures
from ..core.models import ProteinVariant
from ..utils.family import (FAMILY_PFAM_IDS, MAX_LENGTH_AA,
                             MIN_LENGTH_AA)


@dataclass
class DiscoveryConfig:
    known_paralogs: list[str]
    family_pfam_ids: list[str] = field(
        default_factory=lambda: list(FAMILY_PFAM_IDS))
    novelty_floor: float = 0.20          # below this is noise
    novelty_ceiling: float = 0.40        # above this is just an ortholog of a known
    cluster_threshold: float = 0.50      # ≥ this with a known paralog = same cluster
    min_length_aa: int = MIN_LENGTH_AA   # family size band (utils/family.py)
    max_length_aa: int = MAX_LENGTH_AA
    structural_tm_min: float = 0.5       # Foldseek TM-score for "same fold"
    min_taxa_breadth: int = 2            # candidate must appear in ≥N species
    breadth_identity_min: float = 0.35   # cross-species candidate corroboration
    signature_min_coverage: float = 0.5  # MSA-signature fallback threshold
    require_evidence_gate: bool = True   # promotion needs ≥1 family-specific
                                         # component (outlier/domain/fold/split)


@dataclass
class Candidate:
    label: str               # source|accession or accession
    variant: ProteinVariant
    score: int = 0
    evidence: dict[str, str] = field(default_factory=dict)

    def verdict(self) -> str:
        if self.score >= 60: return "promising"
        if self.score >= 40: return "worth manual review"
        return "weak"


@dataclass
class DiscoveryReport:
    config: DiscoveryConfig
    candidates: list[Candidate] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def by_verdict(self) -> dict[str, list[Candidate]]:
        bins: dict[str, list[Candidate]] = {"promising": [], "worth manual review": [], "weak": []}
        for c in self.candidates:
            bins[c.verdict()].append(c)
        return bins

    def text_summary(self) -> str:
        bins = self.by_verdict()
        lines = [
            f"Discovery report — searching for novel paralogs of "
            f"{', '.join(self.config.known_paralogs)}",
            f"  Promising:           {len(bins['promising'])}",
            f"  Worth manual review: {len(bins['worth manual review'])}",
            f"  Weak:                {len(bins['weak'])}",
        ]
        for verdict in ("promising", "worth manual review"):
            for c in bins[verdict]:
                lines.append("")
                lines.append(f"[{c.score:3d}] {verdict.upper()} — {c.label}")
                lines.append(f"      species: {c.variant.species or 'unknown'}")
                lines.append(f"      length:  {c.variant.length_aa or '?'} aa")
                for k, v in c.evidence.items():
                    lines.append(f"      {k}: {v}")
        if self.notes:
            lines.append("")
            lines.append("Notes:")
            for n in self.notes:
                lines.append(f"  • {n}")
        return "\n".join(lines)


def _label(v: ProteinVariant) -> str:
    return f"{v.source}|{v.accession}|{v.gene_symbol}"


def _is_known(symbol: str, known_upper: set[str]) -> bool:
    """Symbol counts as a known paralog on exact or prefix match, so
    lineage-specific duplicates (itpr1a / itpr1b) don't pollute the
    candidate list."""
    u = symbol.upper()
    return any(u == k or u.startswith(k) for k in known_upper)


def build_signature_set(
    variants: list[ProteinVariant],
    msa_rows,                                  # analysis MSA (AlignmentRow list)
    analysis_label_for: dict[str, str],
    known_paralogs: list[str],
) -> Optional[FamilySignatureSet]:
    """Derive the family-signature PSSM set from the KNOWN-paralog rows of
    an existing analysis MSA. Returns None when too few known rows exist."""
    known = {g.upper() for g in known_paralogs}
    known_labels = {
        analysis_label_for.get(_label(v))
        for v in variants
        if _is_known(v.gene_symbol, known) and analysis_label_for.get(_label(v))
    }
    rows = [(r.label, r.aligned) for r in (msa_rows or []) if r.label in known_labels]
    return derive_family_signatures(rows)


def _nearest_identity(
    label: str, distances: DistanceTable, restrict: Optional[set[str]] = None,
) -> tuple[Optional[str], float]:
    """Nearest neighbour by identity; `restrict` limits the comparison set
    (e.g. to known-paralog labels)."""
    if label not in distances.labels:
        return None, 0.0
    i = distances.labels.index(label)
    best = (-1, 0.0)
    for j, other in enumerate(distances.labels):
        if j == i or (restrict is not None and other not in restrict):
            continue
        if distances.identity[i][j] > best[1]:
            best = (j, distances.identity[i][j])
    if best[0] < 0:
        return None, 0.0
    return distances.labels[best[0]], best[1]


def discover_novel_paralogs(
    variants: Iterable[ProteinVariant],
    config: DiscoveryConfig,
    analysis_label_for: Optional[dict[str, str]] = None,
    distances: Optional[DistanceTable] = None,
    domain_hits: Optional[dict[str, list]] = None,
    foldseek_hits: Optional[list[ProteinVariant]] = None,
    signature_set=None,   # analysis.motifs.FamilySignatureSet, optional
) -> DiscoveryReport:
    """Score every variant against the known-paralog set.

    `analysis_label_for` maps `_label(v)` → the short label used in the
    analysis MSA/distance matrix (`AnalysisResult._label_for`). If None,
    distances are skipped (cluster + sequence-outlier components score 0).
    """
    variants = list(variants)
    known = {g.upper() for g in config.known_paralogs}
    candidates: list[Candidate] = []
    notes: list[str] = []

    # Group by gene symbol to compute taxonomic breadth per unknown gene.
    by_gene_species: dict[str, set[str]] = defaultdict(set)
    for v in variants:
        if v.species:
            by_gene_species[v.gene_symbol.upper()].add(v.species)

    domain_hits = domain_hits or {}
    foldseek_by_acc = {fv.accession: fv for fv in (foldseek_hits or [])}

    # Maps needed for the distance-based components: analysis label → species
    # and the set of labels belonging to known paralogs.
    label_species: dict[str, str] = {}
    known_labels: set[str] = set()
    if analysis_label_for:
        for v2 in variants:
            albl = analysis_label_for.get(_label(v2))
            if albl:
                label_species[albl] = v2.species
                if _is_known(v2.gene_symbol, known):
                    known_labels.add(albl)

    for v in variants:
        # Skip variants that already represent known paralogs.
        if _is_known(v.gene_symbol, known):
            continue
        # Skip Foldseek raw hits — they're evidence, not candidates per se.
        if v.source == "Foldseek":
            continue

        c = Candidate(label=_label(v), variant=v)
        # Family-specific evidence (outlier/domain/fold/split) — required
        # for promotion when the evidence gate is on.
        gate_evidence = False

        # --- Component 1: size plausibility (cheap filter) ---
        if v.length_aa and config.min_length_aa <= v.length_aa <= config.max_length_aa:
            c.score += 15
            c.evidence["size"] = f"{v.length_aa} aa within family band {config.min_length_aa}-{config.max_length_aa}"
        else:
            c.evidence["size"] = f"{v.length_aa or '?'} aa — outside family band"

        # --- Component 2: domain signature ---
        # Primary: a real Pfam hit — either a domain_hits lookup, or the
        # candidate was *enumerated from InterPro by family Pfam ID* (the
        # pfam_id rides in v.raw; not crediting it would deny the census's
        # own selection criterion to every enumerated row). Fallback for
        # candidates that aren't in UniProt/InterPro at all (Compara gene
        # models, fresh BLAST hits): coverage of the MSA-derived family
        # signature blocks.
        sig_hits = domain_hits.get(v.accession, [])
        hit_ids = {getattr(h, 'pfam_id', '') for h in sig_hits}
        if v.source == "InterPro":
            hit_ids |= set(str(v.raw.get("pfam_id", "")).split("+"))
        hit_ids &= set(config.family_pfam_ids)
        if hit_ids:
            c.score += 20
            gate_evidence = True
            c.evidence["domain"] = f"Pfam signature present ({','.join(sorted(hit_ids))})"
        elif signature_set is not None and v.sequence:
            n_hit, n_total, _ = signature_set.scan(v.sequence)
            cov = n_hit / n_total if n_total else 0.0
            if cov >= config.signature_min_coverage:
                c.score += 15
                gate_evidence = True
                c.evidence["domain"] = (
                    f"MSA-signature fallback: {n_hit}/{n_total} conserved "
                    f"family blocks present (no InterPro record needed)")
            else:
                c.evidence["domain"] = (
                    f"weak family-signature coverage ({n_hit}/{n_total} blocks)")
        else:
            c.evidence["domain"] = "no family Pfam hit (or InterPro lookup skipped)"

        # --- Component 3: structural fold via Foldseek ---
        fs = foldseek_by_acc.get(v.accession)
        if fs:
            tm_raw = fs.raw.get("tmScore") or fs.raw.get("tm_score") or "0"
            try:
                tm = float(tm_raw)
            except (TypeError, ValueError):
                tm = 0.0
            if tm >= config.structural_tm_min:
                c.score += 20
                gate_evidence = True
                c.evidence["fold"] = f"Foldseek TM={tm:.2f} ≥ {config.structural_tm_min}"
            else:
                c.evidence["fold"] = f"Foldseek hit but TM={tm:.2f} below threshold"
        else:
            c.evidence["fold"] = "no Foldseek hit (skipped or not run)"

        # --- Component 4 + 6: sequence outlier + cluster exclusion ---
        # Both are measured against the nearest KNOWN paralog: novelty means
        # "twilight-zone distance from the named family members", regardless
        # of how similar the candidate is to sibling candidates (a clade of
        # unnamed orthologs must not mask itself).
        mapped = (analysis_label_for or {}).get(_label(v))
        cluster_taxa: set[str] = set()
        if distances and mapped and mapped in distances.labels:
            nearest_any, id_any = _nearest_identity(mapped, distances)
            if nearest_any is not None:
                c.evidence["nearest_in_set"] = f"{nearest_any} at {id_any:.0%} identity"
            nearest_known, id_known = _nearest_identity(mapped, distances,
                                                        restrict=known_labels or None)
            if nearest_known is not None:
                if config.novelty_floor <= id_known <= config.novelty_ceiling:
                    c.score += 20
                    gate_evidence = True
                    c.evidence["outlier"] = (
                        f"identity to nearest known ({nearest_known}) = {id_known:.0%} — "
                        f"paralog twilight zone "
                        f"[{config.novelty_floor:.0%}-{config.novelty_ceiling:.0%}]"
                    )
                else:
                    c.evidence["outlier"] = (
                        f"identity to nearest known ({nearest_known}) = {id_known:.0%} "
                        f"— outside twilight zone")
                if id_known < config.cluster_threshold:
                    c.score += 10
                    c.evidence["cluster"] = (
                        f"outside known paralog clusters (best link "
                        f"{nearest_known} = {id_known:.0%} < {config.cluster_threshold:.0%})")
                else:
                    c.evidence["cluster"] = f"links to {nearest_known} at {id_known:.0%}"
            # Cross-species corroboration for the breadth component: other
            # NON-known variants (sibling candidates) in any species that
            # are ≥ breadth_identity_min identical. Known paralogs are
            # excluded — every family member matches those.
            i = distances.labels.index(mapped)
            cluster_taxa.add(v.species)
            for j, other in enumerate(distances.labels):
                if j == i or other in known_labels:
                    continue
                if distances.identity[i][j] >= config.breadth_identity_min:
                    sp = label_species.get(other, "")
                    if sp:
                        cluster_taxa.add(sp)
        else:
            c.evidence["outlier"] = "distance matrix unavailable"

        # --- Component 5: taxonomic breadth ---
        # Same-symbol occurrences plus cross-species sibling candidates at
        # >= breadth_identity_min (unnamed ortholog sets have a different
        # id per species — mutual candidate similarity is the signal).
        taxa = set(by_gene_species.get(v.gene_symbol.upper(), set())) | cluster_taxa
        if len(taxa) >= config.min_taxa_breadth:
            c.score += 15
            c.evidence["breadth"] = f"present in {len(taxa)} species: {', '.join(sorted(taxa)[:6])}"
        else:
            c.evidence["breadth"] = f"only {len(taxa)} species — could be one-off"

        # --- Component 7: homology provenance ---
        if v.source in ("Compara", "BLAST"):
            c.score += 10
            c.evidence["homology"] = (
                f"surfaced by family-linked homology search ({v.source}), "
                "not by gene name")

        # --- Component 8: split-annotation suspect ---
        if v.raw.get("split_suspect"):
            c.score += 10
            gate_evidence = True
            c.evidence["annotation"] = (
                f"possible split gene model — {v.raw['split_suspect']} "
                "(verify the locus against NCBI's gene span)")

        c.score = min(c.score, 100)

        # --- Promotion gate: generic components alone can't promote ---
        if config.require_evidence_gate and c.score >= 40 and not gate_evidence:
            c.evidence["gate"] = (
                f"capped {c.score} → 39: no family-specific evidence "
                "(twilight-zone outlier, domain signature, fold, or "
                "split-annotation)")
            c.score = 39
        candidates.append(c)

    candidates.sort(key=lambda c: -c.score)
    # Dedup by (gene_symbol, species, length) — same biological hit from
    # multiple sources collapses into the top-scored representative.
    seen: set[tuple[str, str, Optional[int]]] = set()
    deduped: list[Candidate] = []
    for c in candidates:
        v = c.variant
        key = (v.gene_symbol.upper(), v.species, v.length_aa)
        if key in seen:
            continue
        seen.add(key)
        deduped.append(c)

    if not deduped:
        notes.append("No candidate genes found outside the known paralog set. "
                     "Broaden the search (more species, lower max-per-source) "
                     "or enable Foldseek / InterPro for structural / domain hits.")
    return DiscoveryReport(config=config, candidates=deduped, notes=notes)


def write_discovery(out_dir: Path, report: DiscoveryReport) -> dict:
    """Persist a discovery report as `report.txt` + `candidates.tsv`.

    Shared by the CLI, the domain scan, and the GUI bundle writer so the
    on-disk schema stays identical everywhere. Returns the paths written.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    text_path = out_dir / "report.txt"
    text_path.write_text(report.text_summary())
    tsv_path = out_dir / "candidates.tsv"
    with tsv_path.open("w") as f:
        f.write("score\tverdict\tlabel\taccession\tspecies\tlength_aa\treviewed\t"
                "size\tdomain\tfold\toutlier\tbreadth\tcluster\tnearest_in_set\n")
        for c in report.candidates:
            ev = c.evidence
            f.write("\t".join([
                str(c.score), c.verdict(), c.label, c.variant.accession,
                c.variant.species or "", str(c.variant.length_aa or ""),
                c.variant.raw.get("reviewed", ""),
                ev.get("size", ""), ev.get("domain", ""), ev.get("fold", ""),
                ev.get("outlier", ""), ev.get("breadth", ""), ev.get("cluster", ""),
                ev.get("nearest_in_set", ""),
            ]) + "\n")
    return {"report_txt": str(text_path), "candidates_tsv": str(tsv_path)}
