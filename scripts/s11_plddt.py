"""S11 — per-domain pLDDT, with the domain boundaries transferred by
pairwise alignment.

The one-number-per-model version of this ("mean pLDDT 73") is close to
useless for a 2,700-residue multi-domain channel: it averages a
well-predicted β-trefoil with the long low-complexity linkers that make up
much of the subunit, and a model can post a respectable mean while the
part a claim rests on is unmodelled. So confidence is reported **per
domain**, and the domains come from the coordinates S0 measured on the
human paralogs (`review_figures/domain_coords.tsv`) rather than from a
residue list typed out of a paper.

Transfer is by pairwise alignment against the human paralog, and which
paralog is chosen is decided by measurement — the model is aligned to all
three and the best-scoring one wins — because most of this panel is not
vertebrate and has no paralog identity to inherit.

Two things the code refuses to do:

  * **Report a domain it could not place.** A reference domain whose span
    aligns to fewer than `MIN_DOMAIN_COVER` of its residues in the model is
    written out with `placed = 0` and no pLDDT. An unplaced domain is
    evidence about the model; a domain silently averaged over whatever
    residues happened to align is not.
  * **Average a B-factor as a pLDDT.** Only AFDB models are scored here.
"""

from __future__ import annotations

import sys
from pathlib import Path
from statistics import median

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s11_lib import PROJECT_ROOT, load_tsv, uniprot_record      # noqa: E402
from s11_struct_io import read_structure, largest_chain         # noqa: E402

sys.path.insert(0, str(PROJECT_ROOT))
from src.analysis.alignment import pairwise_align               # noqa: E402

DOMAIN_COORDS = (PROJECT_ROOT / "results" / "s0_baseline" / "review_figures"
                 / "domain_coords.tsv")

#: A domain must land on at least this fraction of its reference span in
#: the model before a confidence number is reported for it.
MIN_DOMAIN_COVER = 0.60

#: AlphaFold's own published bands. Not thresholds this project chose.
PLDDT_CONFIDENT = 70.0
PLDDT_VERY_HIGH = 90.0

COLUMNS = ["id", "accession", "paralog", "group", "reference", "ref_symbol",
           "align_identity", "pfam", "label", "shared_class", "ref_start",
           "ref_end", "ref_span", "model_start", "model_end", "n_residues",
           "placed", "mean_plddt", "median_plddt", "frac_ge_70",
           "frac_ge_90", "note"]


def reference_domains(call: str = "ITPR") -> dict[str, list[dict]]:
    """S0's measured Pfam coordinates, per human reference accession."""
    out: dict[str, list[dict]] = {}
    for row in load_tsv(DOMAIN_COORDS):
        if row.get("family") != call:
            continue
        out.setdefault(row["accession"], []).append(row)
    for rows in out.values():
        rows.sort(key=lambda r: int(r["start"]))
    return out


def reference_sequences(accessions: list[str]) -> dict[str, str]:
    seqs = {}
    for acc in accessions:
        rec = uniprot_record(acc)
        if rec:
            seqs[acc] = (rec.get("sequence") or {}).get("value", "")
    return seqs


def _chain_of(path: Path):
    return largest_chain(read_structure(path))


def covered_identity(ref_aln: str, qry_aln: str) -> float:
    """Identity over mutually covered columns.

    `pairwise_align` returns a raw BLOSUM62 *score*, which scales with
    length and sign — the first version of this module ranked the three
    human references on it and printed "identity = 10181.00" into a table.
    The metric here is the project's own covered-only identity
    (`src/analysis/distance.py`), which is what makes a 1,000-residue
    fungal model and a 2,700-residue vertebrate one comparable.
    """
    matches = covered = 0
    for a, b in zip(ref_aln, qry_aln):
        if a == "-" or b == "-":
            continue
        covered += 1
        matches += (a == b)
    return round(matches / covered, 4) if covered else 0.0


def _index_map(ref_aln: str, qry_aln: str) -> dict[int, int]:
    """Reference residue number → model residue index (1-based, ungapped)."""
    mapping: dict[int, int] = {}
    ri = qi = 0
    for a, b in zip(ref_aln, qry_aln):
        if a != "-":
            ri += 1
        if b != "-":
            qi += 1
        if a != "-" and b != "-":
            mapping[ri] = qi
    return mapping


def score_model(row: dict, refs: dict[str, list[dict]],
                ref_seqs: dict[str, str]) -> list[dict]:
    """Per-domain pLDDT for one AFDB model."""
    path = Path(row["path"])
    chain = _chain_of(path)
    plddt = [r.bfactor for r in chain.residues]
    seq = chain.sequence

    best = None
    for acc, dom_rows in sorted(refs.items()):
        ref_seq = ref_seqs.get(acc, "")
        if not ref_seq:
            continue
        a, b, _score = pairwise_align(ref_seq, seq)
        ident = covered_identity(a, b)
        if best is None or ident > best[2]:
            best = (acc, dom_rows, ident, a, b)
    if best is None:
        return [{"id": row["id"], "note": "no reference sequence available"}]
    acc, dom_rows, ident, ref_aln, qry_aln = best
    mapping = _index_map(ref_aln, qry_aln)

    out: list[dict] = []
    covered: set[int] = set()
    for dom in dom_rows:
        rs, re_ = int(dom["start"]), int(dom["end"])
        hits = [mapping[i] for i in range(rs, re_ + 1) if i in mapping]
        span = re_ - rs + 1
        base = {
            "id": row["id"], "accession": row.get("source_id", ""),
            "paralog": row.get("paralog", ""), "group": row.get("group", ""),
            "reference": acc, "ref_symbol": dom.get("symbol", ""),
            "align_identity": ident,
            "pfam": dom["pfam"], "label": dom["label"],
            "shared_class": dom.get("shared_class", ""),
            "ref_start": rs, "ref_end": re_, "ref_span": span,
        }
        if len(hits) < MIN_DOMAIN_COVER * span:
            out.append({**base, "model_start": "", "model_end": "",
                        "n_residues": len(hits), "placed": 0,
                        "mean_plddt": "", "median_plddt": "",
                        "frac_ge_70": "", "frac_ge_90": "",
                        "note": f"only {len(hits)}/{span} reference residues "
                                "aligned"})
            continue
        vals = [plddt[i - 1] for i in hits if 1 <= i <= len(plddt)]
        covered.update(hits)
        out.append({**base, "model_start": min(hits), "model_end": max(hits),
                    "n_residues": len(vals), "placed": 1,
                    "mean_plddt": round(sum(vals) / len(vals), 2),
                    "median_plddt": round(median(vals), 2),
                    "frac_ge_70": round(sum(v >= PLDDT_CONFIDENT
                                            for v in vals) / len(vals), 4),
                    "frac_ge_90": round(sum(v >= PLDDT_VERY_HIGH
                                            for v in vals) / len(vals), 4),
                    "note": ""})

    # The contrast row: everything the reference architecture does not
    # cover. Without it "the pore is at 85" has nothing to be high against.
    rest = [plddt[i - 1] for i in range(1, len(plddt) + 1) if i not in covered]
    if rest:
        out.append({
            "id": row["id"], "accession": row.get("source_id", ""),
            "paralog": row.get("paralog", ""), "group": row.get("group", ""),
            "reference": acc, "ref_symbol": "", "align_identity": ident,
            "pfam": "-", "label": "outside annotated domains",
            "shared_class": "none", "ref_start": "", "ref_end": "",
            "ref_span": "", "model_start": "", "model_end": "",
            "n_residues": len(rest), "placed": 1,
            "mean_plddt": round(sum(rest) / len(rest), 2),
            "median_plddt": round(median(rest), 2),
            "frac_ge_70": round(sum(v >= PLDDT_CONFIDENT
                                    for v in rest) / len(rest), 4),
            "frac_ge_90": round(sum(v >= PLDDT_VERY_HIGH
                                    for v in rest) / len(rest), 4),
            "note": ""})
    return out


def run(manifest: list[dict], on_progress=None) -> list[dict]:
    models = [r for r in manifest
              if r.get("role") == "model" and r.get("status") == "ok"]
    refs = reference_domains("ITPR")
    ref_seqs = reference_sequences(sorted(refs))
    rows: list[dict] = []
    for i, row in enumerate(models, start=1):
        rows += score_model(row, refs, ref_seqs)
        if on_progress:
            on_progress(i, len(models))
    return rows


def domain_summary(rows: list[dict]) -> list[list]:
    """Per domain across the panel — the table the claim is read off."""
    by_label: dict[str, list[dict]] = {}
    for r in rows:
        if r.get("placed") != 1:
            continue
        by_label.setdefault(r["label"], []).append(r)
    out = []
    order = {d["label"]: i for i, d in
             enumerate(sum(reference_domains("ITPR").values(), []))}
    for label in sorted(by_label, key=lambda l: order.get(l, 99)):
        rs = by_label[label]
        means = sorted(float(r["mean_plddt"]) for r in rs)
        out.append([
            label, rs[0].get("pfam", ""), rs[0].get("shared_class", ""),
            len(rs), round(sum(means) / len(means), 2),
            round(median(means), 2), means[0], means[-1],
            sum(1 for m in means if m >= PLDDT_CONFIDENT),
        ])
    return out


SUMMARY_COLUMNS = ["label", "pfam", "shared_class", "n_models", "mean_plddt",
                   "median_plddt", "min_plddt", "max_plddt", "n_ge_70"]
