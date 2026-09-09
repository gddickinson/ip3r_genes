"""S18's deliverable — the corrections, and the veto that keeps them honest.

A correction list is the one artefact in this project addressed to somebody
else, so every row has to survive being read by the curator who owns the
record. Three rules make that possible.

**D6 — the integrity veto.** If S15's ORF screen says a locus is lesion-rich,
the audit does not tell RefSeq to resurrect it. The veto is applied as a
*column*, not a filter: the row is written with `vetoed = 1` and the reason,
because a locus this audit declined to correct is evidence about the audit and
a reader should be able to count them.

**Every correction carries the evidence that would falsify it.** `evidence`
names the archived file a curator can open — the sweep's own miniprot GFF for
a genome correction, the blastp hit table for a protein one — and the row
carries the numbers the class fired on rather than a verdict alone.

**Priority is a rule, not an impression.** `high` needs the gene demonstrably
present *and* the assembly demonstrably able to carry it *and* the reading
frame intact; anything resting on a partial recovery or an unscored ORF is
`medium`; anything below D4's contiguity bar is `low`, because there the
annotation's silence may be the assembly's fault and not the annotator's.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402

#: A genome correction needs this much of the gene recovered before the
#: annotation's silence is the annotation's fault.
MIN_RECOVERY = 0.50
#: …and this much for the correction to be `high` priority.
STRONG_RECOVERY = 0.90
#: A protein-record family correction needs an absolute score as well as D7's
#: relative margin. `s3_assign.MIN_SCORE` (30 bits) is the floor for *calling*
#: a record; telling a database its name is wrong is a stronger claim, and the
#: 5 candidates here score 90-191 bits over ~2,700 residues, which is weak
#: evidence for a rename however large the relative margin looks.
MIN_RENAME_BITS = 200.0

CLASSES = {
    "C1_unannotated": "a recovered gene with no same-strand annotated feature",
    "C2_noncoding_only": "the annotation names the gene and files it as "
                         "non-coding, so no protein is served",
    "C3_incomplete": "the annotated coding models deliver less than the bar",
    "C4_wrong_family_name": "the covering model is named for the sister "
                            "family (D14's failure, at the annotation)",
    "C5_wrong_paralog_name": "the covering model names a different paralog "
                             "from the one the alignment assigns",
    "C6_protein_wrong_family": "a protein record named for the sister family "
                               "that the bait panel calls this one",
}


def _priority(row: dict, cov: float) -> tuple[str, str]:
    spans = str(row.get("contig_spans_gene")) in ("1", "True", "true")
    intact = row.get("integrity_verdict") == "intact"
    if not spans:
        return "low", ("the contig cannot carry the whole gene (D4) — the "
                       "annotation's silence may be the assembly's")
    if cov >= STRONG_RECOVERY and intact:
        return "high", (f"gene recovered at coverage {cov:.2f}, contig spans "
                        f"it, reading frame intact (S15)")
    if cov >= MIN_RECOVERY:
        return "medium", (f"gene recovered at coverage {cov:.2f}; integrity "
                          f"{row.get('integrity_verdict')}")
    return "low", f"only {cov:.2f} of the gene recovered"


def from_loci(rows: list[dict], bar: float = R.COMPLETE_FRAC) -> list[dict]:
    """Genome-side corrections, one row per locus that needs one."""
    out = []
    for r in rows:
        if r["gene_set"] != "present":
            continue
        cov = float(r.get("coverage") or 0)
        if cov < MIN_RECOVERY:
            continue
        state, best = r["state"], float(r.get("best_model_frac") or 0)
        cls = proposal = ""
        if state == "unannotated":
            cls = "C1_unannotated"
            # A locus can reach this state with a coding model on it that
            # covers under 1 % of the footprint — a passenger in an intron.
            # Saying "no model" there would be false to the curator reading it.
            near = (f" (the nearest same-strand coding model, "
                    f"'{r['namer_symbol'] or r['namer_product']}', reaches "
                    f"{best:.3f} of it)" if int(r["n_coding_models"] or 0)
                    else "")
            proposal = (f"add a protein-coding gene model over the "
                        f"{r['locus_cds_bp']} bp of aligned coding sequence "
                        f"on {r['contig']}:{r['start']}-{r['end']} "
                        f"({r['strand']}){near}")
        elif state == "noncoding":
            cls = "C2_noncoding_only"
            proposal = (f"the locus is held only by {r['n_noncoding_models']} "
                        f"non-coding feature(s) ({r['noncoding_biotype']}); the "
                        f"alignment gives a coding model over the same exons")
        elif state in ("fragmentary", "split"):
            cls = "C3_incomplete"
            proposal = (f"{r['n_coding_models']} coding model(s) deliver "
                        f"{r['union_coding_frac']:.2f} of the gene, best alone "
                        f"{best:.2f}; extend or merge to the aligned exons")
        if r["symbol_verdict"] == "wrong_family" or \
                r["product_verdict"] == "wrong_family":
            cls = "C4_wrong_family_name"
            proposal = (f"the covering model is named "
                        f"'{r['namer_symbol'] or r['namer_product']}' — the "
                        f"alignment assigns it to {r['cell']} with family "
                        f"margin from the sweep")
        elif r["symbol_verdict"] == "correct_family_wrong_paralog" or \
                r["product_verdict"] == "correct_family_wrong_paralog":
            cls = "C5_wrong_paralog_name"
            proposal = (f"the covering model claims "
                        f"{r['symbol_claim'] or r['product_claim']}; the "
                        f"alignment assigns {r['cell']}")
        if not cls:
            continue
        prio, prio_reason = _priority(r, cov)
        vetoed = int(r.get("integrity_verdict") == "elevated_lesions")
        out.append({
            "correction_id": f"{r['accession']}:{r['cell']}:{r['locus_idx']}",
            "side": "genome", "cls": cls, "class_rule": CLASSES[cls],
            "assembly": r["accession"], "source": r["source"],
            "organism": r["organism"], "vclass": r["vclass"],
            "cell": r["cell"], "contig": r["contig"], "start": r["start"],
            "end": r["end"], "strand": r["strand"],
            "current_state": r["state"],
            "current_name": r["namer_symbol"] or r["namer_product"] or "(none)",
            "current_biotype": r["namer_biotype"] or r["top_overlapper_biotype"],
            "coverage": cov, "best_model_frac": best,
            "union_coding_frac": r["union_coding_frac"],
            "contig_spans_gene": r.get("contig_spans_gene", ""),
            "integrity_verdict": r.get("integrity_verdict", ""),
            "proposal": proposal, "priority": prio,
            "priority_reason": prio_reason,
            "vetoed": vetoed,
            "veto_reason": ("S15 scored the reading frame lesion-rich (D6) — "
                            "reported, not proposed") if vetoed else "",
            "evidence": f"<data_root>/genome_sweep/{r['accession']}/miniprot.gff"
                        f" (model over {r['contig']}:{r['start']}-{r['end']})",
        })
    return out


def from_proteins(rows: list[dict]) -> list[dict]:
    """Protein-side corrections — naming only; a database record has no exons.

    Only the family question is proposed. A wrong *paralog* on a protein record
    is reported in the audit table but not raised here: the panel's paralog
    margins on the disputed records run 0.10-0.36, which clears D7 and is still
    thin evidence for asking a curator to rename a reviewed entry, and the
    three ITPR cases sit at 0.10-0.12 — on the threshold itself.
    """
    out = []
    for r in rows:
        if r["name_verdict"] != "wrong_family" and \
                r["symbol_verdict"] != "wrong_family":
            continue
        bits = max(float(r["itpr_bits"]), float(r["ryr_bits"]))
        strong = (r["seq_family"] in ("ITPR", "RYR")
                  and bits >= MIN_RENAME_BITS)
        out.append({
            "correction_id": r["accession"], "side": "protein",
            "cls": "C6_protein_wrong_family",
            "class_rule": CLASSES["C6_protein_wrong_family"],
            "assembly": "", "source": "UniProtKB", "organism": r["species"],
            "vclass": r.get("tax_class", ""), "cell": r["seq_family"],
            "contig": "", "start": "", "end": "", "strand": "",
            "current_state": "protein_record",
            "current_name": r["protein_name"] or r["gene"] or "(none)",
            "current_biotype": "", "coverage": "",
            "best_model_frac": "", "union_coding_frac": "",
            "contig_spans_gene": "", "integrity_verdict": "",
            "proposal": (f"named for the sister family; best bait score "
                         f"{r['itpr_bits']} ITPR vs {r['ryr_bits']} RyR "
                         f"(relative margin {r['family_rel_margin']})"),
            "priority": "high" if strong else "low",
            "priority_reason": (
                f"best family score {bits:.0f} bits "
                f"{'>=' if strong else '<'} {MIN_RENAME_BITS:g}; a rename needs "
                f"an absolute score as well as D7's relative margin"),
            "vetoed": 0, "veto_reason": "",
            "evidence": "<data_root>/s18/protein_audit_hits.tsv",
        })
    return out


def summarise(rows: list[dict]) -> list[dict]:
    """Counts by class × priority, with the vetoed rows kept visible."""
    keys: dict[tuple[str, str], dict] = {}
    for r in rows:
        k = (r["cls"], r["priority"])
        d = keys.setdefault(k, {"cls": r["cls"], "class_rule": r["class_rule"],
                                "priority": r["priority"], "n": 0,
                                "n_vetoed": 0, "n_refseq": 0, "n_genbank": 0})
        d["n"] += 1
        d["n_vetoed"] += int(r["vetoed"])
        d["n_refseq"] += int(r["source"] == "RefSeq")
        d["n_genbank"] += int(r["source"] == "GenBank")
    order = list(CLASSES) + ["(other)"]
    prio = {"high": 0, "medium": 1, "low": 2}
    return sorted(keys.values(),
                  key=lambda d: (order.index(d["cls"]) if d["cls"] in order
                                 else 99, prio.get(d["priority"], 9)))


COLS = ["correction_id", "side", "cls", "class_rule", "priority",
        "priority_reason", "vetoed", "veto_reason", "assembly", "source",
        "organism", "vclass", "cell", "contig", "start", "end", "strand",
        "current_state", "current_name", "current_biotype", "coverage",
        "best_model_frac", "union_coding_frac", "contig_spans_gene",
        "integrity_verdict", "proposal", "evidence"]
