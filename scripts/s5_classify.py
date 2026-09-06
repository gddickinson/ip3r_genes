"""s5_classify.py — per-locus annotation matching and per-cell status.

Turns clustered miniprot loci into the ledger's genome x class cells: the
CDS-footprint test that decides whether a locus is already an annotated gene,
and the family-then-paralog assignment behind each cell's status.

**D14 lives here.** The family call comes first and is a positive test: a
locus belongs to whichever family's baits win it on alignment score, and a
locus won by the RyR baits is never offered to an ITPR paralog cell — not
capped, not down-weighted, not admitted. Only inside the winning family does
the paralog question get asked. Without that ordering the sweep would do
exactly what S1 measured the discovery scorer doing before the sister test
existed: promote every ryanodine receptor, because they carry the shared
architecture and there are three of them in every vertebrate genome.

The RYR cell is classified by the same machinery but read as a control:
presence, not which RyR.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.utils import family                    # noqa: E402
from s5_genome_io import longest_n_run          # noqa: E402
from s5_sweep_lib import (CLASSES, CONTROL_CLASS, COV_FOUND,  # noqa: E402
                          EDGE_BP, NGAP_RUN, UNRESOLVED_ITPR_CLADES, Locus)

#: Substrings that make an annotated gene name a family name. `ip3r` and
#: `insp3r` appear in RefSeq gene models that predate the ITPR symbols;
#: `iplA` is the *Dictyostelium* name S0 found the family signature does not
#: even reach, kept here because the same annotation habit reaches other
#: non-model assemblies.
#: Sourced from the family definition rather than retyped (CLAUDE.md's rule:
#: the family lives in one place). The local copy had drifted — it was missing
#: `itr-1`, the *C. elegans* receptor's own symbol, so S23's pilot read a
#: correctly-recovered nematode gene as an annotation that names something
#: else. No vertebrate symbol changes, so S5's calls are unaffected.
ITPR_NAME_HINTS = tuple(
    sorted({h.lower() for h in family.KNOWN_NAME_SUBSTRINGS}
           | {"inositol 1,4,5-trisphosphate receptor"}))
RYR_NAME_HINTS = ("ryr", "ryanodine receptor")

ANNOT_CDS_FRAC = 0.50     # gene must cover this much of the alignment's CDS bp


def _locus_flags(locus: Locus, seqlens: dict[str, int], region_seq: str) -> dict:
    contig_len = seqlens.get(locus.contig, 0)
    edge = locus.start <= EDGE_BP or (contig_len and locus.end >= contig_len - EDGE_BP)
    ngap = longest_n_run(region_seq) if region_seq else 0
    return {"contig_edge": bool(edge), "longest_n_run": ngap,
            "n_gap": ngap >= NGAP_RUN}


def locus_dict(locus: Locus, clade: str) -> dict:
    a = locus.best_of(clade) or locus.best
    return {"contig": locus.contig, "strand": locus.strand,
            "start": locus.start, "end": locus.end,
            "bait": a.bait, "identity": round(a.identity, 4),
            "coverage": round(a.coverage, 4),
            "span_coverage": round(a.span_coverage, 4),
            "aligned_aa": a.aligned_aa, "score": a.score,
            "q_span": [a.q_start, a.q_end], "bait_len": a.bait_len,
            "frameshifts": a.frameshifts, "stop_codons": a.stop_codons,
            "mp_id": a.mp_id, "locus_top_clade": locus.clade,
            "locus_family": locus.family,
            "family_margin": locus.family_margin(),
            "paralog_margin": locus.paralog_margin(clade),
            "n_baits": len(locus.alns)}


def name_family(name: str) -> str:
    """Which family an annotated gene *name* claims, or "" if neither.

    Used only to describe an annotation, never to make a call — the call is
    the alignment's. A locus whose alignment says ITPR and whose overlapping
    gene is named RYR is precisely the correction S18 is looking for, so both
    have to be recorded side by side rather than reconciled here.
    """
    low = (name or "").lower()
    if any(h in low for h in RYR_NAME_HINTS):
        return "RYR"
    if any(h in low for h in ITPR_NAME_HINTS):
        return "ITPR"
    return ""


def name_paralog(name: str) -> str:
    """The paralog an annotated gene name claims (ITPR2, RYR3), or ""."""
    fam = name_family(name)
    if not fam:
        return ""
    low = (name or "").lower()
    for digit in "123":
        if f"{fam.lower()}{digit}" in low or f"type {digit}" in low:
            return f"{fam}{digit}"
        if fam == "ITPR" and (f"ip3r{digit}" in low or f"insp3r{digit}" in low):
            return f"ITPR{digit}"
    return ""


def _overlap_bp(blocks: list[tuple[int, int]], s: int, e: int) -> int:
    return sum(max(0, min(e, b2) - max(s, b1) + 1) for b1, b2 in blocks)


def _annotate_locus(d: dict, gene_index, cds_blocks: list, cell_class: str) -> None:
    """Attach the annotated gene that actually covers the alignment's exons.

    Scoring against the CDS footprint (not the locus span) is what keeps a
    small gene sitting inside an ITPR intron from being mistaken for the ITPR
    annotation — an intronic passenger covers ~0 % of the CDS blocks. ITPR
    introns run to 152 kb (`intron_calibration.tsv`), so there is room for
    several.
    """
    d["overlapping_genes"] = []
    d["annot_gene"] = None
    d["annot_family"] = ""
    d["annot_paralog"] = ""
    d["annot_paralog_matches"] = False
    if gene_index is None:
        return
    genes = gene_index.overlapping(d["contig"], d["start"], d["end"])
    span = max(1, d["end"] - d["start"] + 1)
    cds_total = max(1, sum(b2 - b1 + 1 for b1, b2 in cds_blocks))
    scored = []
    for g in genes:
        ov = min(d["end"], g[2]) - max(d["start"], g[1]) + 1
        if ov <= 0:
            continue
        cds_ov = _overlap_bp(cds_blocks, g[1], g[2])
        scored.append({"name": g[5] or g[4], "biotype": g[6], "overlap_bp": ov,
                       "frac_locus": round(ov / span, 3),
                       "frac_gene": round(ov / max(1, g[2] - g[1] + 1), 3),
                       "frac_cds": round(cds_ov / cds_total, 3),
                       "start": g[1], "end": g[2], "strand": g[3]})
    scored.sort(key=lambda x: x["frac_cds"], reverse=True)
    d["overlapping_genes"] = scored[:8]
    d["cds_bp"] = cds_total
    good = [g for g in scored if g["frac_cds"] >= ANNOT_CDS_FRAC]
    want = CONTROL_CLASS if cell_class == CONTROL_CLASS else "ITPR"
    named = [g for g in good if name_family(g["name"]) == want]
    pick = named[0] if named else (good[0] if good else None)
    if pick:
        d["annot_gene"] = pick
        d["annot_family"] = name_family(pick["name"])
        d["annot_paralog"] = name_paralog(pick["name"])
        d["annot_paralog_matches"] = (
            d["annot_paralog"] == cell_class if cell_class in CLASSES
            else d["annot_family"] == CONTROL_CLASS)


def best_canonical_paralog(locus: Locus) -> str | None:
    """Which of ITPR1/2/3 scores highest at this locus. ITPR loci only."""
    hits = [a for a in locus.alns if a.clade in CLASSES]
    return max(hits, key=lambda a: a.score).clade if hits else None


def cell_loci(loci: list[Locus], cell_class: str) -> tuple[list, list]:
    """(primary, secondary) loci for one cell — D14 ordered family then paralog.

    Primary: the locus's top-scoring bait is this cell's own clade.

    Secondary (ITPR cells only): the locus is won by the **ITPR** baits, its
    top bait carries no paralog assignment (`vertebrate_basal` — the
    chimaera, lamprey and gar seeds, which exist precisely because no
    labelled paralog record exists in those bands), and among the paralog
    baits that do align there this cell's scores highest. Admitting these
    keeps a real gene in a shark or a lamprey from reading as absent merely
    because the species' own nearest bait outranks every labelled one.

    No coverage bar on the offer. The PIEZO port gated an unresolved-clade
    locus *before* it entered the cell while gating a primary locus *after*,
    so identical partial evidence produced `fragment` for one and a flat
    `absent` for the other; admitting the locus and letting the coverage
    logic grade it is the consistent rule.

    A locus whose family is RYR never enters an ITPR cell by either route.
    """
    if cell_class == CONTROL_CLASS:
        return [loc for loc in loci if loc.family == CONTROL_CLASS], []
    primary = [loc for loc in loci
               if loc.family == "ITPR" and loc.clade == cell_class]
    secondary = [loc for loc in loci
                 if loc.family == "ITPR"
                 and loc.clade in UNRESOLVED_ITPR_CLADES
                 and best_canonical_paralog(loc) == cell_class]
    return primary, secondary


def classify_class(loci: list[Locus], cell_class: str, gene_index,
                   seqlens: dict, fetch_fn) -> dict:
    """Classify one cell. fetch_fn(contig, start, end) -> str."""
    primary, secondary = cell_loci(loci, cell_class)
    mine = primary + secondary
    cell: dict = {"class": cell_class, "n_loci": len(mine),
                  "n_primary": len(primary), "n_secondary": len(secondary),
                  "is_control": cell_class == CONTROL_CLASS, "loci": []}
    if not mine:
        cell["status"] = "no_locus"        # rescue decides absent vs trace
        return cell
    mine.sort(key=lambda L: (L.best_of(cell_class) or L.best).score, reverse=True)
    for loc in mine:
        d = locus_dict(loc, cell_class)
        d["assigned_via"] = ("primary" if loc in primary
                             else f"unresolved:{loc.clade}")
        region = fetch_fn(loc.contig, loc.start - 2000, loc.end + 2000)
        d.update(_locus_flags(loc, seqlens, region))
        _annotate_locus(d, gene_index,
                        (loc.best_of(cell_class) or loc.best).cds_blocks,
                        cell_class)
        cell["loci"].append(d)
    best = cell["loci"][0]
    if best["coverage"] >= COV_FOUND:
        if gene_index is None:
            cell["status"] = "found_no_annotation"   # assembly has no GFF
        else:
            cell["status"] = ("found_annotated" if best["annot_gene"]
                              else "found_unannotated")
    elif best["contig_edge"] or best["n_gap"]:
        cell["status"] = "assembly_gap"
    else:
        cell["status"] = "fragment"
    cell["best_coverage"] = best["coverage"]
    cell["best_identity"] = best["identity"]
    cell["best_family_margin"] = best["family_margin"]
    cell["best_paralog_margin"] = best["paralog_margin"]
    cell["family_named_annotation"] = bool(best.get("annot_family"))
    cell["paralog_named_annotation"] = bool(best.get("annot_paralog_matches"))
    cell["annot_paralog"] = best.get("annot_paralog", "")
    cell["best_assigned_via"] = best.get("assigned_via", "primary")
    return cell


def other_loci(loci: list[Locus], gene_index=None, seqlens: dict | None = None,
               fetch_fn=None) -> list[dict]:
    """Loci no cell claimed — the review list.

    An ITPR-family locus topped by an unresolved bait that no paralog cell
    took, and anything the panel scores but cannot place. Reported rather
    than folded into a cell: in the chondrichthyan, coelacanth and cyclostome
    bands the panel carries no labelled ITPR1/ITPR2 bait at all (six unfilled
    slots in `bait_manifest.tsv`), so an unplaced locus there is expected and
    is S7's question, not this task's.
    """
    claimed = set()
    for cell_class in (*CLASSES, CONTROL_CLASS):
        p, s = cell_loci(loci, cell_class)
        claimed.update(id(loc) for loc in p + s)
    out = []
    for loc in loci:
        if id(loc) in claimed:
            continue
        d = locus_dict(loc, loc.clade)
        d["clade"] = loc.clade
        d["best_canonical_paralog"] = best_canonical_paralog(loc)
        if seqlens is not None and fetch_fn is not None:
            region = fetch_fn(loc.contig, loc.start - 2000, loc.end + 2000)
            d.update(_locus_flags(loc, seqlens, region))
        _annotate_locus(d, gene_index, loc.best.cds_blocks, loc.family)
        out.append(d)
    out.sort(key=lambda d: d["score"], reverse=True)
    return out
