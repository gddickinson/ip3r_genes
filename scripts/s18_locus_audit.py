"""S18 stage `loci` — every gene-scale locus in the genome scope, scored.

One pass per genome: the sweep's archived `miniprot.gff` for the loci's own
coding blocks, the assembly's `genomic.gff.gz` for what the annotation puts
over them. Both readers are `s10_gff`'s multi-window ones, so a 274-genome
audit costs one stream of each file rather than one per locus.

The bar is applied **after** the sweep, not during it: every locus's
measurements and both name verdicts are bar-independent, so the audit measures
first, calibrates on what it measured, and only then assigns a state. That is
also what makes `sensitivity()` cheap — the whole grid is a re-labelling of one
measurement pass, and every count in it is the same evidence read at a
different bar.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_gff                                                    # noqa: E402
import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402
from s5_sweep_lib import parse_miniprot_gff                       # noqa: E402

WINDOW_PAD = 20_000     # a model may start outside the alignment's first exon


def _locus_rows(acc: str, summary: dict, genome: dict) -> list[dict]:
    """Measure every locus of one genome. No state, no bar."""
    loci = L.loci_of(acc, summary)
    if not loci:
        return []
    sw = L.sweep_dir() / acc
    mp = sw / "miniprot.gff"
    # The sweep's own parser, not a second one. miniprot numbers alignments
    # from MP000001 *per run*, so a chunked genome's concatenated GFF repeats
    # every ID and `s5_sweep_lib.parse_miniprot_gff` disambiguates the
    # duplicates as `MP000001#<n>` — which is the id the summary records.
    # Reading the raw `ID=` attribute instead missed 21 loci in 4 giant
    # genomes, and the fallback then measured the annotation against the whole
    # locus *span* (2.4 Mb for one ITPR2) rather than its coding blocks, which
    # forces `unannotated` whatever the annotation holds.
    models = {a.mp_id: a for a in parse_miniprot_gff(mp, {})} \
        if mp.exists() else {}

    gff = L.gff_path(acc)
    gene_set = ("present" if gff else
                ("unavailable" if genome.get("annotated") == "Y" else "none"))
    windows = []
    for l in loci:
        windows.append((l["contig"], max(1, int(l["start"]) - WINDOW_PAD),
                        int(l["end"]) + WINDOW_PAD))
    ann = s10_gff.read_annotation_windows(gff, windows) if gff else {}

    out = []
    for l, win in zip(loci, windows):
        model = models.get(l["mp_id"])
        # No fallback to the locus span: a span includes every intron, and
        # substituting one would put a 2.4 Mb denominator under a 8 kb gene.
        # A locus whose model cannot be recovered is reported unscorable.
        cds = s10_gff.merge(model.cds_blocks) if model and model.cds_blocks \
            else []
        over = R.models_over(ann.get(win, {"genes": []}), cds,
                             l["strand"]) if cds else []
        cell = l["cell"]
        out.append({
            "accession": acc, "organism": l["organism"], "vclass": l["vclass"],
            "vorder": l["vorder"], "cell": cell, "locus_idx": l["locus_idx"],
            "is_control": int(bool(l["is_control"])),
            "source": L.annotation_source(acc),
            "assembly_level": l["assembly_level"],
            "contig": l["contig"], "start": int(l["start"]),
            "end": int(l["end"]), "strand": l["strand"],
            "bait": L.bait_id(l.get("bait", "")),
            "coverage": l.get("coverage", 0.0),
            "identity": l.get("identity", 0.0),
            "frameshifts": l.get("frameshifts", 0),
            "stop_codons": l.get("stop_codons", 0),
            "cell_status": l["cell_status"],
            "gene_set": gene_set if cds else "cds_unavailable",
            "_models": over, "_cds": cds,
        })
    return out


def apply_verdicts(rows: list[dict]) -> None:
    """Name and symbol verdicts, derived from `_models`, in place.

    Kept out of `measure()` on purpose. The measurement pass is four minutes
    of gzipped GFF parsing and is cached, and the first version cached the
    verdicts with it — so when `name_verdict` gained `paralog_unspecified`
    the committed table still carried the old labels and mouse *Itpr1* read as
    the annotation naming a different paralog. Everything a rule decides is
    recomputed on load; only what the parser found is cached.
    """
    for r in rows:
        over = r["_models"]
        top = over[0] if over else {}
        best_coding = next((m for m in over if m["coding"]), {})
        namer = R.naming_model(over)
        family = "RYR" if r["cell"] == "RYR" else "ITPR"
        expected = "" if r["cell"] == "RYR" else r["cell"]
        sym_v, sym_claim = R.name_verdict(namer.get("symbol", ""),
                                          expected, family)
        prod_v, prod_claim = R.name_verdict(namer.get("product", ""),
                                            expected, family)
        r.update({
            "best_gene_id": best_coding.get("gene_id", ""),
            "best_symbol": best_coding.get("symbol", ""),
            "best_product": best_coding.get("product", ""),
            "best_biotype": best_coding.get("biotype", ""),
            "namer_gene_id": namer.get("gene_id", ""),
            "namer_symbol": namer.get("symbol", ""),
            "namer_product": namer.get("product", ""),
            "namer_biotype": namer.get("biotype", ""),
            "namer_coding": int(bool(namer.get("coding"))),
            "namer_frac_cds": namer.get("frac_cds", 0.0),
            "namer_own_frac": namer.get("own_frac", 0.0),
            "namer_model_aa": (namer.get("model_bp", 0) or 0) // 3,
            "top_overlapper_biotype": top.get("biotype", ""),
            "top_overlapper_pseudo": int(bool(top.get("pseudo"))),
            "symbol_verdict": sym_v, "symbol_claim": sym_claim,
            "product_verdict": prod_v, "product_claim": prod_claim,
        })


def measure(summaries: dict[str, dict], genomes: dict[str, dict],
            on_progress=None) -> list[dict]:
    rows = []
    for i, (acc, s) in enumerate(sorted(summaries.items()), start=1):
        rows.extend(_locus_rows(acc, s, genomes.get(acc, {})))
        if on_progress and i % 25 == 0:
            on_progress(i, len(summaries))
    return rows


def apply_bar(rows: list[dict], bar: float,
              min_piece: float = R.MIN_PIECE_FRAC) -> None:
    """Assign each locus its state at `bar`, in place."""
    for r in rows:
        st = R.locus_state(r["_models"], r["_cds"], r["gene_set"], bar,
                           min_piece)
        r.update(st)


def contiguity_join(rows: list[dict], ledger: dict, integrity: dict) -> None:
    """D4's bar and S15's ORF verdict beside every locus (the D6 veto).

    Both are joined rather than recomputed: `contig_spans_gene` is S5's own
    measurement and the integrity verdict is S15's, and a second derivation of
    either would be a fork in what this project means by 'the assembly can
    carry the gene' and 'the reading frame is broken'.
    """
    for r in rows:
        led = ledger.get((r["accession"], r["cell"]), {})
        r["contig_spans_gene"] = led.get("contig_spans_gene", "")
        r["contig_n50"] = led.get("contig_n50", "")
        ig = integrity.get((r["accession"], r["cell"], r["locus_idx"]), {})
        r["integrity_verdict"] = ig.get("verdict", "not_scored")
        r["lesion_density"] = ig.get("lesion_density", "")


def sensitivity(rows: list[dict], bars: list[float],
                pieces: list[float]) -> list[dict]:
    """Every state count over the bar × piece grid.

    A state boundary is a decision, and the honest form of a decision is the
    grid of what it changes (S15b's sensitivity matrix, applied to the one
    threshold S18 owns). `unannotated` and `noncoding` do not move with the
    bar by construction, which is worth being able to see rather than assert.
    """
    scorable = [r for r in rows
                if r["gene_set"] == "present" and not r["is_control"]]
    out = []
    for bar in bars:
        for mp in pieces:
            counts = {s: 0 for s in R.STATES}
            for r in scorable:
                st = R.locus_state(r["_models"], r["_cds"], "present", bar, mp)
                counts[st["state"]] += 1
            out.append({"bar": bar, "min_piece": mp, "n_loci": len(scorable),
                        **{f"n_{k}": v for k, v in counts.items()}})
    return out


AUDIT_COLS = [
    "accession", "organism", "vclass", "vorder", "source", "assembly_level",
    "cell", "locus_idx", "is_control", "contig", "start", "end", "strand",
    "bait", "coverage", "identity", "frameshifts", "stop_codons",
    "cell_status", "contig_spans_gene", "contig_n50", "gene_set", "state",
    "state_reason", "locus_cds_bp", "n_models", "n_coding_models",
    "n_coding_pieces", "n_noncoding_models", "best_model_frac",
    "union_coding_frac", "noncoding_biotype", "top_overlapper_biotype",
    "top_overlapper_pseudo", "best_gene_id", "best_symbol", "best_product",
    "best_biotype", "namer_gene_id", "namer_symbol", "namer_product",
    "namer_biotype", "namer_coding", "namer_frac_cds", "namer_own_frac",
    "namer_model_aa",
    "symbol_verdict", "symbol_claim", "product_verdict",
    "product_claim", "integrity_verdict", "lesion_density",
]
