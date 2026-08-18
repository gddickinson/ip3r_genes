"""Ensembl Compara homology client — phylogenomic paralog mining.

The gene-symbol clients only return proteins already *named* like the query
gene. Compara flips that: Ensembl's gene trees already cluster every gene in
every genome, so asking "what are the within-species paralogues of ITPR1?"
surfaces family members regardless of what they're called — including
clone-named placeholders (`si:dkey-…`, `zgc:…`) and gene models with *no*
name at all — which is how a real paralog stays uncounted.

Strategy per (gene symbol × species):

    GET /homology/symbol/{species}/{symbol}?type=paralogues;sequence=none
        → anchor gene + paralogue targets with %identity + taxonomy level
    POST /lookup/id  (bulk)
        → display names, descriptions, genomic coordinates, protein lengths
    GET /sequence/id/{protein_id}  (optional, when include_sequence)

If the query has no species, a vertebrate panel is swept (see
lineages). As a bonus, paralogue loci that sit within `SPLIT_SUSPECT_BP` of
each other on the same chromosome are flagged as possible split gene models
— the annotation-error pattern that hides a real gene behind two partial
models, which is common in ~2,700-residue, ~58-exon genes like these.
"""

from __future__ import annotations

import time
from typing import Any, Optional

import requests

from ..core.models import ProteinVariant, SearchQuery
from ..utils.species import ensembl_species_slug
from .base import DatabaseClient

ENSEMBL_REST = "https://rest.ensembl.org"

# Vertebrate sweep used when the query has no species: the classic models
# plus the lineages that decide the ITPR paralog questions — a teleost with
# 3R duplicates, a monotreme, and the chondrichthyan outgroup that sits
# below the teleost duplication.
COMPARA_PANEL = [
    "homo_sapiens",
    "mus_musculus",
    "danio_rerio",
    "gallus_gallus",
    "xenopus_tropicalis",
    "ornithorhynchus_anatinus",
    "callorhinchus_milii",
]

SPLIT_SUSPECT_BP = 300_000   # paralogue loci closer than this = suspect split
_PAUSE_S = 0.12              # stay well under Ensembl's 15 req/s limit


class ComparaClient(DatabaseClient):
    name = "Compara"

    def search(self, query: SearchQuery) -> list[ProteinVariant]:
        species_slugs = (
            [ensembl_species_slug(query.species)] if query.species else list(COMPARA_PANEL)
        )
        # gene_id -> collected info
        genes: dict[str, dict[str, Any]] = {}
        for symbol in query.gene_symbols:
            for slug in species_slugs:
                for entry in self._paralogues(slug, symbol):
                    self._collect(genes, entry, slug, symbol)

        if not genes:
            return []

        lookup = self._bulk_lookup(
            list(genes.keys())
            + [g["protein_id"] for g in genes.values() if g.get("protein_id")]
        )
        variants = [self._to_variant(gene_id, info, lookup)
                    for gene_id, info in genes.items()]
        variants = [v for v in variants if v is not None][: query.max_results_per_source]
        self._flag_split_suspects(variants)
        if query.include_sequence:
            self._fetch_sequences(variants)
        return variants

    # ---- REST calls ------------------------------------------------------

    def _get_json(self, path: str, **params: Any) -> Any:
        params.setdefault("content-type", "application/json")
        r = requests.get(f"{ENSEMBL_REST}{path}", params=params, timeout=self.timeout_s)
        time.sleep(_PAUSE_S)
        if r.status_code == 404 or r.status_code == 400:
            return None
        r.raise_for_status()
        return r.json()

    def _paralogues(self, slug: str, symbol: str) -> list[dict]:
        """Homology entries for the symbol; tries case variants (human genes
        are upper-case, fish lower-case)."""
        for sym in dict.fromkeys([symbol, symbol.lower(), symbol.upper()]):
            try:
                data = self._get_json(
                    f"/homology/symbol/{slug}/{sym}",
                    type="paralogues", sequence="none",
                )
            except requests.RequestException:
                continue
            if data and data.get("data"):
                return data["data"]
        return []

    def _collect(self, genes: dict, entry: dict, slug: str, symbol: str) -> None:
        anchor_id = entry.get("id", "")
        homologies = entry.get("homologies", [])
        if anchor_id and anchor_id not in genes:
            src = (homologies[0].get("source", {}) if homologies else {})
            genes[anchor_id] = {
                "protein_id": src.get("protein_id", ""),
                "species_slug": slug, "anchor": symbol,
                "perc_id": None, "taxonomy_level": "", "type": "anchor",
            }
        for h in homologies:
            t = h.get("target", {})
            gid = t.get("id", "")
            if not gid:
                continue
            # Keep the best-identity record if the same gene shows up as a
            # paralogue of several anchors.
            prev = genes.get(gid)
            perc = t.get("perc_id")
            if prev is not None and prev.get("perc_id") is not None \
                    and (perc is None or perc <= prev["perc_id"]):
                continue
            genes[gid] = {
                "protein_id": t.get("protein_id", ""),
                "species_slug": t.get("species", slug), "anchor": symbol,
                "perc_id": perc, "taxonomy_level": h.get("taxonomy_level", ""),
                "type": h.get("type", "paralog"),
            }

    def _bulk_lookup(self, ids: list[str]) -> dict[str, dict]:
        out: dict[str, dict] = {}
        for i in range(0, len(ids), 500):
            chunk = ids[i : i + 500]
            try:
                r = requests.post(
                    f"{ENSEMBL_REST}/lookup/id",
                    json={"ids": chunk},
                    params={"content-type": "application/json"},
                    timeout=self.timeout_s,
                )
                time.sleep(_PAUSE_S)
                r.raise_for_status()
                data = r.json()
            except (requests.RequestException, ValueError):
                continue
            for k, v in data.items():
                if isinstance(v, dict):
                    out[k] = v
        return out

    def _fetch_sequences(self, variants: list[ProteinVariant]) -> None:
        for v in variants:
            pid = v.raw.get("protein_id", "")
            if not pid:
                continue
            try:
                data = self._get_json(f"/sequence/id/{pid}", type="protein")
            except requests.RequestException:
                continue
            if data and data.get("seq"):
                v.sequence = data["seq"]
                if not v.length_aa:
                    v.length_aa = len(v.sequence)

    # ---- variant assembly ------------------------------------------------

    def _to_variant(self, gene_id: str, info: dict, lookup: dict) -> Optional[ProteinVariant]:
        gene = lookup.get(gene_id, {})
        protein = lookup.get(info.get("protein_id", ""), {})
        slug = info["species_slug"]
        species = slug.replace("_", " ").capitalize()
        display = gene.get("display_name", "") or ""
        desc = (gene.get("description", "") or "").split(" [Source:")[0]
        if info["type"] == "anchor":
            rel = f"Compara anchor gene ({info['anchor']})"
        else:
            pid = info.get("perc_id")
            pid_txt = f"{pid:.0f}% id" if pid is not None else "id n/a"
            rel = (f"Compara {info['type']} of {info['anchor']} "
                   f"({pid_txt} @ {info.get('taxonomy_level') or '?'})")
        if not display:
            rel += " — UNNAMED gene model"
        raw = {
            "protein_id": info.get("protein_id", ""),
            "homology_type": info["type"],
            "anchor": info["anchor"],
            "perc_id": "" if info.get("perc_id") is None else f"{info['perc_id']:.1f}",
            "taxonomy_level": info.get("taxonomy_level", ""),
            "chrom": str(gene.get("seq_region_name", "")),
            "start": str(gene.get("start", "")),
            "end": str(gene.get("end", "")),
        }
        return ProteinVariant(
            source=self.name,
            accession=info.get("protein_id") or gene_id,
            gene_symbol=display or gene_id,
            species=species,
            length_aa=protein.get("length") or None,
            description=(desc + " · " if desc else "") + rel,
            gene_id=gene_id,
            url=f"https://www.ensembl.org/{slug.capitalize()}/Gene/Summary?g={gene_id}",
            raw=raw,
        )

    def _flag_split_suspects(self, variants: list[ProteinVariant]) -> None:
        """Two family loci within SPLIT_SUSPECT_BP on one chromosome usually
        mean one real gene split across two models (the split-model pattern)."""
        by_chrom: dict[tuple[str, str], list[ProteinVariant]] = {}
        for v in variants:
            chrom, start = v.raw.get("chrom"), v.raw.get("start")
            if chrom and start:
                by_chrom.setdefault((v.species, chrom), []).append(v)
        for group in by_chrom.values():
            group.sort(key=lambda v: int(v.raw["start"]))
            for a, b in zip(group, group[1:]):
                gap = int(b.raw["start"]) - int(a.raw["end"] or a.raw["start"])
                if gap < SPLIT_SUSPECT_BP:
                    for v, other in ((a, b), (b, a)):
                        v.raw["split_suspect"] = f"within {gap:,} bp of {other.gene_symbol}"
                        v.description += (f" ⚠ possible split gene model — "
                                          f"{gap:,} bp from {other.gene_symbol}")
