"""S19 — per gene: what a protein-database search would have missed.

Split out of `s19_contribution` to keep both inside the 500-line budget, and
because it asks a different question with a different denominator. The record
counts in `s19_contribution` are about the *archive*; this module is about
the *genes*, asked per genome x cell, with the genome sweep as the ground
truth for where they are.

Three ways to miss a gene, kept apart because they call for different fixes:
the species has no reference proteome at all, it has one and no family record
exists for the species, or records exist and none resolves to that paralog.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402

_SCOPE: dict[str, str] = {}


def _scope_of(accession: str) -> str:
    """Why S4 put this genome in scope: order representative, margin, or both.

    Read from the manifest's own `reasons` field rather than re-derived — S4
    stated four margin rules with thresholds and 61 genomes carry more than
    one, so a second derivation here could disagree with the denominator the
    whole project is scoped against.
    """
    if not _SCOPE:
        for r in S.read_tsv(S.MANIFEST):
            reasons = {x for x in (r.get("reasons") or "").split(";") if x}
            # S4 writes the order-representative reason as `order_rep:<order>`
            # and the four margin rules by name, so the test is on the prefix.
            is_rep = any(x.startswith("order_rep") for x in reasons)
            is_margin = any(not x.startswith("order_rep") for x in reasons)
            _SCOPE[r["accession"]] = ("both" if is_rep and is_margin else
                                      "margin species" if is_margin else
                                      "order representative")
    return _SCOPE.get(accession, "order representative")


def _binomial(name: str) -> str:
    return " ".join((name or "").replace("_", " ").split()[:2])


def _proteome_species() -> tuple[dict[str, str], dict[str, str]]:
    """(binomial -> proteome id, taxid -> proteome id) over every swept DB.

    Both keys are kept because they disagree: the proteome sets carry
    subspecies and sister taxids to the assemblies in S4's manifest. Matching
    on the binomial is the more generous reading of "a proteome search could
    have seen this gene" and therefore the conservative one for the claim
    made here; both counts are reported.
    """
    out, by_taxid = {}, {}
    paths = [S.HMM_DIR / "proteome_manifest.tsv"]
    paths += [S.S20_DIR / f"proteome_manifest_{g}.tsv"
              for g in S.NONVERT_GROUPS]
    for path in paths:
        if not path.exists():
            continue
        for r in S.read_tsv(path):
            name = (r.get("Organism") or "").strip()
            pid = (r.get("Proteome Id") or "").strip()
            taxid = (r.get("Organism Id") or "").strip()
            if name:
                out.setdefault(_binomial(name), pid)
            if taxid:
                by_taxid.setdefault(taxid, pid)
    return out, by_taxid


#: Census sources that are protein-database records rather than gene models
#: this project produced. A genome-sweep row in the census cannot answer
#: "would a protein search have found this gene", because it *is* the genome
#: search's answer.
PROTEIN_SOURCES = ("InterPro", "HMM", "S20 sweep")


def _protein_records() -> tuple[dict[tuple, int], dict[str, int]]:
    """(binomial, cell) -> n resolvable records, and binomial -> n records.

    Two tiers, and they answer different questions.

    The **species** tier is every family record in census v6 that came from a
    protein database, at any length — the generous reading of "a protein
    search could have seen something here", and therefore the conservative
    one for the claim made from it.

    The **cell** tier needs the record's paralog, and takes it from S18's
    **sequence** call against the 38-bait panel rather than from a gene
    symbol: a record filed under the wrong paralog still holds the gene, and
    a record with no name at all still holds it. S18 scores only full-length
    records, so a species can have family records and still resolve to no
    cell — which is a real state of the archive and gets its own channel
    rather than being folded into either neighbour.
    """
    by_cell: dict[tuple, int] = {}
    by_species: dict[str, int] = {}
    for r in S.census_rows("v6"):
        if r.get("source") not in PROTEIN_SOURCES:
            continue
        if r.get("call") not in ("ITPR", "RYR"):
            continue
        sp = _binomial(r.get("species", ""))
        if sp:
            by_species[sp] = by_species.get(sp, 0) + 1
    for r in S.read_tsv(S.RESULTS / "annotation_audit" / "protein_audit.tsv"):
        sp = _binomial(r.get("species", ""))
        cell = r.get("seq_paralog", "")
        if not sp:
            continue
        if cell in S.PARALOGS or cell.startswith("RYR"):
            key = (sp, S.CONTROL_CELL if cell.startswith("RYR") else cell)
            by_cell[key] = by_cell.get(key, 0) + 1
    return by_cell, by_species


def gene_recovery(log=S.log) -> tuple[list[list], dict]:
    """Per genome x cell: which channels could have found this gene.

    The headline needs its scope composition printed beside it. S4 chose 169
    of the 309 genomes as **margin species** precisely because their
    proteomes are empty, missing a paralog or fragment-only, so a high
    protein-database-invisible rate over the whole scope is partly the scope
    working as designed. The rate is therefore also reported split by whether
    a genome is in scope as an order representative or as a margin species.
    """
    prot_species, prot_taxid = _proteome_species()
    by_cell, by_species = _protein_records()
    rows = []
    taxid_matched = set()
    for c in S.ledger_cells():
        sp = _binomial(c["organism"])
        present = (c["cell"] == S.CONTROL_CELL
                   or c["s15_state"] in S.PRESENT_STATES)
        has_proteome = sp in prot_species
        if str(c.get("taxid", "")) in prot_taxid:
            taxid_matched.add(c["accession"])
        n_cell = by_cell.get((sp, c["cell"]), 0)
        n_any = by_species.get(sp, 0)
        n_resolved = sum(v for (s, _cc), v in by_cell.items() if s == sp)
        if not present:
            channel = "n/a (no gene demonstrated in this cell)"
        elif not has_proteome:
            channel = "genome_only:no_reference_proteome"
        elif n_any == 0:
            channel = "genome_only:no_family_record_for_species"
        elif n_cell == 0 and n_resolved == 0:
            channel = "genome_only:records_exist_but_none_full_length"
        elif n_cell == 0:
            channel = "genome_only:no_record_resolves_to_this_paralog"
        else:
            channel = "protein_database_and_genome"
        rows.append([c["accession"], c["organism"], c["vclass"], c["cell"],
                     c["status"], int(present), int(has_proteome),
                     prot_species.get(sp, ""), n_any, n_resolved, n_cell,
                     channel])
    S.write_tsv(S.out_dir() / "gene_recovery.tsv",
                ["accession", "organism", "vclass", "cell", "ledger_status",
                 "gene_present", "has_reference_proteome", "proteome_id",
                 "n_family_protein_records_species",
                 "n_records_resolving_to_any_cell",
                 "n_records_resolving_to_cell", "recovery_channel"], rows)

    genes = [r for r in rows if r[5]]
    by_channel: dict[str, int] = {}
    for r in genes:
        by_channel[r[11]] = by_channel.get(r[11], 0) + 1
    invisible = sum(v for k, v in by_channel.items()
                    if k.startswith("genome_only"))
    summary = {
        "genomes": len({r[0] for r in rows}),
        "genomes_with_reference_proteome": len({r[0] for r in rows if r[7]}),
        "genomes_with_reference_proteome_taxid_match": len(taxid_matched),
        "genes_present": len(genes), "by_channel": by_channel,
        "protein_db_invisible": invisible,
        "protein_db_invisible_frac": round(invisible / max(1, len(genes)), 4),
    }
    scope_rows = []
    for scope in ("order representative", "margin species", "both"):
        sub = [r for r in genes if _scope_of(r[0]) == scope]
        if not sub:
            continue
        miss = [r for r in sub if r[11].startswith("genome_only")]
        lo, hi = S.wilson(len(miss), len(sub))
        scope_rows.append([scope, len({r[0] for r in sub}), len(sub),
                           len(miss), round(len(miss) / len(sub), 4),
                           round(lo, 4), round(hi, 4)])
    S.write_tsv(S.out_dir() / "gene_recovery_by_scope.tsv",
                ["scope", "n_genomes", "n_genes_present",
                 "n_protein_db_invisible", "frac", "wilson_lo", "wilson_hi"],
                scope_rows)
    summary["by_scope"] = {r[0]: {"n_genomes": r[1], "n_genes": r[2],
                                  "n_invisible": r[3], "frac": r[4]}
                           for r in scope_rows}

    per_cell = []
    for cell in S.CELLS:
        sub = [r for r in genes if r[3] == cell]
        miss = [r for r in sub if r[11].startswith("genome_only")]
        lo, hi = S.wilson(len(miss), len(sub))
        per_cell.append([cell, len(sub), len(miss),
                         round(len(miss) / max(1, len(sub)), 4),
                         round(lo, 4), round(hi, 4)])
    S.write_tsv(S.out_dir() / "gene_recovery_by_cell.tsv",
                ["cell", "n_genes_present", "n_protein_db_invisible", "frac",
                 "wilson_lo", "wilson_hi"], per_cell)
    log(f"gene_recovery: {invisible}/{len(genes)} demonstrated genes are "
        f"invisible to any protein-database search")
    return rows, summary


def db_status(log=S.log) -> list[list]:
    """How the databases hold each gene model the sweeps produced.

    S5b and S23c already recorded this per locus; S19 puts the two sweeps'
    counts on one table, because "how many ITPR genes exist only as DNA" is a
    methods number and it is currently split across two census versions.
    """
    rows = []
    for label, path in (("vertebrate (S5b)",
                         S.RESULTS / "census_v4" / "db_status_counts.tsv"),
                        ("non-vertebrate (S23c)",
                         S.RESULTS / "census_v6" / "db_status_counts.tsv")):
        for r in S.read_tsv(path):
            # The two sweeps wrote the same table under different column
            # names (S5b `cell`/`db_status`/`n_loci`, S23c `group`/`status`/
            # `loci`), so the keys are resolved rather than assumed.
            rows.append([label, r.get("cell") or r.get("group", ""),
                         r.get("db_status") or r.get("status", ""),
                         int(S.fnum(r.get("n_loci") or r.get("loci"),
                                    float, 0))])
    S.write_tsv(S.out_dir() / "db_status_combined.tsv",
                ["sweep", "cell_or_group", "db_status", "n_loci"], rows)
    log(f"db_status_combined: {len(rows)} rows")
    return rows


