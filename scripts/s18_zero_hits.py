"""S18 — the S3 proteomes that returned nothing, resolved against their genomes.

S3 swept 764 vertebrate reference proteomes with both profiles and 15 came back
with no family hit at all. On its own that is uninterpretable: an ITPR-shaped
hole in a proteome is either a gene the species does not have, or a gene its
*gene caller* did not find. The two are told apart by asking the genome, and
this project has already asked — S5 swept an assembly of every one of those 15
species with an alignment-based instrument that owes the gene caller nothing.

The verdict vocabulary is deliberately four-valued, and the fourth value is the
one that keeps the answer honest: `undecidable_no_genome` exists because a
species with no assembly in the S4 scope cannot be resolved either way, and
reporting it as an absence would report the genome scope as biology. It happens
to fire on none of the 15 here — every zero-hit species has a genome in scope —
which is worth being able to say rather than assume.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_lib as L                                               # noqa: E402
from s6_lib import binomial                                       # noqa: E402

#: A cell counts as recovered from the genome at this alignment coverage. The
#: sweep's own `found` statuses already encode a bar; this one is stated here
#: because the question asked of it is different — not "did the sweep place the
#: gene" but "is there enough of it in the assembly to prove the proteome is
#: missing something".
GENOME_FOUND_COVERAGE = 0.50

VERDICTS = ("gene_caller_missed_it", "genome_also_empty",
            "assembly_cannot_carry_it", "undecidable_no_genome")


def resolve(zero_rows: list[dict], genomes: dict[str, dict],
            summaries: dict[str, dict]) -> list[dict]:
    """One row per zero-hit proteome, with the number its verdict fired on."""
    by_species: dict[str, list[str]] = {}
    for acc, g in genomes.items():
        by_species.setdefault(binomial(g["organism"]), []).append(acc)

    out = []
    for z in zero_rows:
        sp = binomial(z["Organism"])
        accs = sorted(by_species.get(sp, []))
        row = {
            "proteome_id": z["Proteome Id"], "organism": z["Organism"],
            "species": sp, "taxon_id": z.get("Organism Id", ""),
            "proteome_proteins": int(z.get("Protein count") or 0),
            "genome_accession": ";".join(accs),
            "n_genomes_in_scope": len(accs),
        }
        if not accs:
            out.append({**row, "verdict": "undecidable_no_genome",
                        "verdict_reason": "no assembly of this species is in "
                                          "the S4 genome scope",
                        "n_itpr_loci": "", "n_recovered": "",
                        "best_coverage": "", "vclass": "",
                        "contig_spans_any": ""})
            continue
        acc = accs[0]
        s = summaries.get(acc, {})
        cells = (s.get("cells") or {})
        itpr = {k: v for k, v in cells.items() if not v.get("is_control")}
        n_loci = sum(len(c.get("loci") or []) for c in itpr.values())
        covs = [float(c.get("best_coverage") or 0) for c in itpr.values()]
        n_rec = sum(1 for c in covs if c >= GENOME_FOUND_COVERAGE)
        spans = any(str(loc.get("contig_edge")) in ("0", "False")
                    for c in itpr.values() for loc in (c.get("loci") or []))
        row.update(vclass=s.get("vclass", ""), n_itpr_loci=n_loci,
                   n_recovered=n_rec,
                   best_coverage=round(max(covs), 4) if covs else 0.0,
                   contig_spans_any=int(bool(spans)))
        if n_rec:
            v = ("gene_caller_missed_it",
                 f"{n_rec} of 3 ITPR cells recovered from the assembly at "
                 f"coverage >= {GENOME_FOUND_COVERAGE:g} while the proteome "
                 f"holds none among {row['proteome_proteins']} proteins")
        elif n_loci:
            v = ("assembly_cannot_carry_it",
                 f"{n_loci} partial ITPR loci in the assembly, none reaching "
                 f"coverage {GENOME_FOUND_COVERAGE:g} — the gene is there but "
                 f"neither record delivers it")
        else:
            v = ("genome_also_empty",
                 "the genomic sweep placed no ITPR locus in this assembly "
                 "either")
        out.append({**row, "verdict": v[0], "verdict_reason": v[1]})
    return out


def summarise(rows: list[dict]) -> list[dict]:
    counts: dict[str, int] = {v: 0 for v in VERDICTS}
    for r in rows:
        counts[r["verdict"]] = counts.get(r["verdict"], 0) + 1
    n = len(rows)
    return [{"verdict": v, "n_proteomes": counts.get(v, 0),
             "frac": L.frac(counts.get(v, 0), n)} for v in VERDICTS]


COLS = ["proteome_id", "organism", "species", "taxon_id", "vclass",
        "proteome_proteins", "genome_accession", "n_genomes_in_scope",
        "n_itpr_loci", "n_recovered", "best_coverage", "contig_spans_any",
        "verdict", "verdict_reason"]
