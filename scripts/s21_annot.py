"""The assemblies' own annotations against the sweep's exon boundaries.

One pass per annotated genome, two questions, and they are different questions.

**Does the instrument agree with an independent pipeline?** Every annotated CDS
block edge over an in-scope locus is compared with the sweep's own exon
boundaries. The assemblies' gene sets were built by NCBI's and the submitters'
pipelines from evidence miniprot never saw, so agreement is corroboration of
the exon structure S21 measures — S10 could only ask this of two cases against
35 and 40 other genomes; here it is asked of every locus in every annotated
genome in the scope.

**Where does a fragmentary annotation stop?** S18 called 264 loci
`fragmentary` and 27 `split`: coding sequence is there but no single model
covers the gene. That says nothing about *where* the annotation stopped, and the
difference is the whole claim — a model that stops at a genuine junction has
produced a plausible short gene, while one that stops in the middle of an exon
has produced a boundary no splicing machinery could make.

So the unit for the second question is the **annotated model's own terminus**,
not its internal exon edges: a model's internal boundaries are its own splice
sites and are canonical by construction, which is why scoring them would answer
the first question a second time and call it the second. Each terminus is scored
on two axes that know nothing about each other — against the alignment's exon
boundaries, and against the two genomic bases immediately outside it read in
gene orientation — and a terminus coinciding with the gene's own end is excluded
by rule, because a real gene legitimately starts and stops inside an exon.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_gff as GFF                                             # noqa: E402
import s21_blocks as B                                            # noqa: E402
import s21_lib as L                                               # noqa: E402

#: S18's states for a locus whose annotation does not deliver the gene.
TARGET_STATES = ("split", "fragmentary")

#: How far an annotated edge may sit from an alignment exon edge and still count
#: as the same boundary. 2 bp allows the terminal-codon ambiguity a spliced
#: aligner and an annotation pipeline can genuinely differ on; 0 is reported
#: beside it so the tolerance is never load-bearing on its own.
EDGE_TOL = 2

CONC_COLS = ("accession", "organism", "vclass", "source", "cell", "locus_idx",
             "state", "n_annot_edges", "n_exact", "n_within_tol", "n_no_match",
             "frac_exact", "frac_within_tol", "n_model_exons",
             "n_annot_blocks", "median_offset")

EDGE_COLS = ("accession", "organism", "cell", "locus_idx", "state", "strand",
             "gene_id", "symbol", "biotype", "model_index", "n_models",
             "terminus", "position", "vs_alignment", "distance_to_exon_edge",
             "dinucleotide", "expected_dinucleotide", "splice_ok",
             "is_gene_terminus", "both_axes_agree")

LOCUS_COLS = ("accession", "organism", "vclass", "cell", "locus_idx", "state",
              "in_architecture_scope", "architecture_excluded_by",
              "n_models", "n_internal_termini", "n_at_exon_boundary",
              "n_mid_exon", "n_in_intron", "n_splice_ok", "n_both_agree",
              "frac_at_exon_boundary", "frac_splice_ok", "verdict",
              "verdict_reason")


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


def _exon_edges(exons: list[dict]) -> tuple[list[int], list[tuple[int, int]]]:
    edges = sorted({e["start"] for e in exons} | {e["end"] for e in exons})
    spans = sorted((e["start"], e["end"]) for e in exons)
    return edges, spans


def _vs_alignment(pos: int, edges: list[int],
                  spans: list[tuple[int, int]]) -> tuple[str, int]:
    near = min((abs(pos - x) for x in edges), default=10 ** 9)
    if near <= EDGE_TOL:
        return "at_exon_boundary", near
    for s, e in spans:
        if s <= pos <= e:
            return "mid_exon", near
    return "in_intron", near


# --------------------------------------------------------------------------
# question 1 — boundary concordance
# --------------------------------------------------------------------------
def annot_blocks(window: dict, contig: str, start: int, end: int,
                 strand: str) -> list[tuple[int, int]]:
    """The annotation's merged coding blocks over the locus, same strand only.

    Same strand because an antisense gene's CDS is not this gene's coding
    sequence, and scored on **CDS blocks and never gene spans** (S18's rule): a
    span covers its own introns, and this family's reach 152 kb.
    """
    blocks = []
    for g in window.get("genes", []):
        if g["strand"] != strand:
            continue
        for s, e in g["cds"]:
            if e >= start and s <= end:
                blocks.append((max(s, start), min(e, end)))
    return GFF.merge(blocks)


def concordance_row(arch: dict, blocks: list[tuple[int, int]],
                    exons: list[dict], source: str) -> dict:
    """How many annotated block edges land on an alignment exon boundary.

    Terminal edges of the annotated set are excluded on the same reasoning as
    the fragment test: the first and last coding positions of a model are its
    start and stop codons, which are inside an exon by definition.
    """
    edges, spans = _exon_edges(exons)
    lo = min(s for s, _ in blocks)
    hi = max(e for _, e in blocks)
    offsets, exact, within, none = [], 0, 0, 0
    n = 0
    for s, e in blocks:
        for pos in (s, e):
            if pos in (lo, hi):
                continue
            n += 1
            verdict, dist = _vs_alignment(pos, edges, spans)
            offsets.append(dist)
            if dist == 0:
                exact += 1
            if verdict == "at_exon_boundary":
                within += 1
            else:
                none += 1
    return {"accession": arch["accession"], "organism": arch["organism"],
            "vclass": arch["vclass"], "source": source, "cell": arch["cell"],
            "locus_idx": arch["locus_idx"], "state": arch.get("state", ""),
            "n_annot_edges": n, "n_exact": exact, "n_within_tol": within,
            "n_no_match": none,
            "frac_exact": round(exact / n, 4) if n else 0.0,
            "frac_within_tol": round(within / n, 4) if n else 0.0,
            "n_model_exons": len(exons), "n_annot_blocks": len(blocks),
            "median_offset": round(L.quantile([float(x) for x in offsets], 0.5), 1)
            if offsets else float("nan")}


# --------------------------------------------------------------------------
# question 2 — where a fragmentary annotation stops
# --------------------------------------------------------------------------
def coding_models(window: dict, contig: str, start: int, end: int,
                  strand: str) -> list[dict]:
    """The annotated coding models over the locus, each with its own termini.

    One entry per gene (its transcripts' CDS blocks merged), because two
    transcripts of one gene are one piece of evidence about where the
    annotation stops.
    """
    out = []
    for g in window.get("genes", []):
        if g["strand"] != strand or not g["cds"]:
            continue
        blocks = [(s, e) for s, e in g["cds"] if e >= start and s <= end]
        if not blocks:
            continue
        out.append({"gene_id": g["gene_id"], "symbol": g.get("name", ""),
                    "biotype": g.get("biotype", ""),
                    "low": min(s for s, _ in blocks),
                    "high": max(e for _, e in blocks),
                    "n_blocks": len(blocks)})
    return sorted(out, key=lambda d: d["low"])


def terminus_rows(arch: dict, region: L.LocusRegion, models: list[dict],
                  exons: list[dict], strand: str) -> list[dict]:
    """One row per annotated model terminus, scored on both axes."""
    edges, spans = _exon_edges(exons)
    gene_lo = min(e["start"] for e in exons)
    gene_hi = max(e["end"] for e in exons)
    rows = []
    for i, m in enumerate(models):
        for side, pos in (("low", m["low"]), ("high", m["high"])):
            # In gene order a plus-strand model's high edge is where it stops
            # (a donor position) and its low edge where it starts (an acceptor).
            if strand == "+":
                role, gene_end = ("donor", gene_hi) if side == "high" \
                    else ("acceptor", gene_lo)
            else:
                role, gene_end = ("acceptor", gene_hi) if side == "high" \
                    else ("donor", gene_lo)
            is_gene_end = int(abs(pos - gene_end) <= EDGE_TOL)
            if role == "donor":
                raw = (region.at(pos + 1, pos + 2) if strand == "+"
                       else _revcomp(region.at(pos - 2, pos - 1)))
                expect = "GT"
            else:
                raw = (region.at(pos - 2, pos - 1) if strand == "+"
                       else _revcomp(region.at(pos + 1, pos + 2)))
                expect = "AG"
            verdict, dist = _vs_alignment(pos, edges, spans)
            ok = int(raw == expect)
            rows.append({
                "accession": arch["accession"], "organism": arch["organism"],
                "cell": arch["cell"], "locus_idx": arch["locus_idx"],
                "state": arch.get("state", ""), "strand": strand,
                "gene_id": m["gene_id"], "symbol": m["symbol"],
                "biotype": m["biotype"], "model_index": i,
                "n_models": len(models),
                "terminus": f"{role}_{side}", "position": pos,
                "vs_alignment": verdict, "distance_to_exon_edge": dist,
                "dinucleotide": raw, "expected_dinucleotide": expect,
                "splice_ok": ok, "is_gene_terminus": is_gene_end,
                "both_axes_agree": int(ok and verdict == "at_exon_boundary")})
    return rows


def locus_verdict(arch: dict, rows: list[dict]) -> dict:
    """The per-locus verdict, one-sided in the way S10's ORF screen is.

    A mid-exon terminus *falsifies* the reading that the annotation stopped at a
    real gene boundary — nothing splices inside an exon. A terminus inside one of
    the model's introns is a different observation and gets its own verdict: the
    annotation is calling coding sequence the model calls intronic, which is the
    two pipelines disagreeing about the structure rather than one of them
    inventing a boundary. And every terminus landing on a boundary does not prove
    the pieces are separate genes; it says the annotation broke the gene at
    junctions the gene has, which is a weaker statement, so that verdict names
    what it is.
    """
    internal = [r for r in rows if not r["is_gene_terminus"]]
    n = len(internal)
    at_edge = sum(1 for r in internal if r["vs_alignment"] == "at_exon_boundary")
    mid = sum(1 for r in internal if r["vs_alignment"] == "mid_exon")
    intr = sum(1 for r in internal if r["vs_alignment"] == "in_intron")
    ok = sum(r["splice_ok"] for r in internal)
    both = sum(r["both_axes_agree"] for r in internal)
    if not n:
        verdict = "no_internal_terminus"
        why = "every annotated model terminus is the gene's own end"
    elif mid:
        verdict = "annotation_failure"
        why = (f"{mid} of {n} internal model termini sit inside an exon of the "
               "gene model, where nothing splices")
    elif intr:
        verdict = "structure_disagreement"
        why = (f"{intr} of {n} internal model termini sit inside an intron of "
               "the gene model — the annotation calls coding what the model "
               "calls intronic")
    elif at_edge == n and ok == n:
        verdict = "broken_at_real_junctions"
        why = (f"all {n} internal termini are exon boundaries the genome reads "
               "as splice sites — the annotation broke the gene between real "
               "junctions, not at invented ones")
    else:
        verdict = "mixed"
        why = (f"{at_edge} of {n} internal termini are exon boundaries and "
               f"{ok} of {n} read as splice sites")
    return {"accession": arch["accession"], "organism": arch["organism"],
            "vclass": arch["vclass"], "cell": arch["cell"],
            "locus_idx": arch["locus_idx"], "state": arch.get("state", ""),
            "in_architecture_scope": arch.get("in_scope", 0),
            "architecture_excluded_by": arch.get("excluded_by", ""),
            "n_models": rows[0]["n_models"] if rows else 0,
            "n_internal_termini": n, "n_at_exon_boundary": at_edge,
            "n_mid_exon": mid, "n_in_intron": intr, "n_splice_ok": ok,
            "n_both_agree": both,
            "frac_at_exon_boundary": round(at_edge / n, 4) if n else 0.0,
            "frac_splice_ok": round(ok / n, 4) if n else 0.0,
            "verdict": verdict, "verdict_reason": why}


# --------------------------------------------------------------------------
# the pass
# --------------------------------------------------------------------------
def run(arch: list[dict], audit: dict, min_intron: int, on_progress=None
        ) -> tuple[list[dict], list[dict], list[dict]]:
    """Concordance rows, terminus rows and fragment verdicts, in one pass.

    Two scopes, and they are deliberately different. **Concordance** is asked of
    every in-scope locus, because corroborating an exon structure needs a whole
    gene model to corroborate. **The fragment question** is asked of every locus
    S18 called `split` or `fragmentary` whatever the architecture rules say
    about it: those loci are the ones the annotation failed, they sit in the
    poorer assemblies by construction, and restricting the question to the
    architecture scope would answer it on 50 of 291 loci and quietly drop the
    hardest ones. Each fragment row carries the architecture rule it fails, so
    a reader can see which loci would not have been asked.
    """
    for r in arch:
        key = (r["accession"], r["cell"], r["locus_idx"])
        r["state"] = (audit.get(key) or {}).get("state", "")
    scoped = {id(r) for r in arch if r["in_scope"]}
    wanted = [r for r in arch
              if r["in_scope"] or r["state"] in TARGET_STATES]
    by_acc: dict[str, list[dict]] = {}
    for r in wanted:
        by_acc.setdefault(r["accession"], []).append(r)
    crows: list[dict] = []
    erows: list[dict] = []
    lrows: list[dict] = []
    accs = sorted(by_acc)
    for i, acc in enumerate(accs, 1):
        gff = L.gff_path(acc)
        fna = L.genome_fna(acc)
        rows = by_acc[acc]
        if gff is None or fna is None:
            continue
        source = "RefSeq" if acc.startswith("GCF_") else "GenBank"
        windows = [(r["contig"], r["start"], r["end"]) for r in rows]
        wins = GFF.read_annotation_windows(gff, windows)
        models = B.read_models(L.sweep_dir(acc) / "miniprot.gff",
                               {r["mp_id"] for r in rows})
        for r in rows:
            m = models.get(r["mp_id"])
            if m is None:
                continue
            ex, _ = B.exons(m, min_intron)
            win = wins.get((r["contig"], r["start"], r["end"])) or {}
            blocks = annot_blocks(win, r["contig"], r["start"], r["end"],
                                  r["strand"])
            if blocks and id(r) in scoped:
                crows.append(concordance_row(r, blocks, ex, source))
            if r["state"] not in TARGET_STATES:
                continue
            cmods = coding_models(win, r["contig"], r["start"], r["end"],
                                  r["strand"])
            if not cmods:
                continue
            region = L.LocusRegion(acc, fna, r["contig"], r["start"] - 10,
                                   r["end"] + 10)
            trows = terminus_rows(r, region, cmods, ex, r["strand"])
            erows += trows
            lrows.append(locus_verdict(r, trows))
        if on_progress:
            on_progress(i, len(accs))
    return crows, erows, lrows


def concordance_summary(crows: list[dict]) -> list[dict]:
    """Pooled concordance by paralogue, by annotation source and by state."""
    out = []

    def block(label: str, key: str, sub: list[dict]) -> dict:
        edges = sum(r["n_annot_edges"] for r in sub)
        return {"grouping": label, "group": key, "n_loci": len(sub),
                "n_genomes": len({r["accession"] for r in sub}),
                "n_annot_edges": edges,
                "n_exact": sum(r["n_exact"] for r in sub),
                "n_within_tol": sum(r["n_within_tol"] for r in sub),
                "frac_exact": round(sum(r["n_exact"] for r in sub) / edges, 4)
                if edges else 0.0,
                "frac_within_tol": round(
                    sum(r["n_within_tol"] for r in sub) / edges, 4)
                if edges else 0.0}
    out.append(block("all", "all", crows))
    for grouping in ("cell", "source", "vclass", "state"):
        for key in sorted({r[grouping] for r in crows}):
            out.append(block(grouping, key or "(none)",
                             [r for r in crows if r[grouping] == key]))
    return out


def fragment_summary(lrows: list[dict]) -> list[dict]:
    """Verdict counts by state — the brief's question, answered."""
    verdicts = sorted({r["verdict"] for r in lrows})
    out = []
    for state in list(TARGET_STATES) + ["all"]:
        sub = [r for r in lrows if state == "all" or r["state"] == state]
        if not sub:
            continue
        d = {"state": state, "n_loci": len(sub),
             "n_genomes": len({r["accession"] for r in sub})}
        for v in verdicts:
            d[v] = sum(1 for r in sub if r["verdict"] == v)
        term = sum(r["n_internal_termini"] for r in sub)
        d["internal_termini"] = term
        d["termini_at_exon_boundary"] = sum(r["n_at_exon_boundary"] for r in sub)
        d["termini_mid_exon"] = sum(r["n_mid_exon"] for r in sub)
        d["termini_in_intron"] = sum(r["n_in_intron"] for r in sub)
        d["termini_splice_ok"] = sum(r["n_splice_ok"] for r in sub)
        d["frac_termini_at_boundary"] = round(
            d["termini_at_exon_boundary"] / term, 4) if term else 0.0
        d["frac_termini_splice_ok"] = round(
            d["termini_splice_ok"] / term, 4) if term else 0.0
        out.append(d)
    return out


def verdicts_from_termini(erows: list[dict], arch: list[dict]) -> list[dict]:
    """Re-derive every fragment verdict from the cached terminus rows.

    S18's D55 applies here: the *measurement* is expensive and cached, and the
    verdicts are not. A rule change that survived into a committed table because
    the verdicts were cached with the measurement has happened once in this
    project already, so nothing here reads a stored verdict.
    """
    idx = {(r["accession"], r["cell"], r["locus_idx"]): r for r in arch}
    groups: dict[tuple, list[dict]] = {}
    for r in erows:
        groups.setdefault((r["accession"], r["cell"], r["locus_idx"]),
                          []).append(r)
    out = []
    for key, rows in sorted(groups.items()):
        a = idx.get(key)
        if a is None:
            continue
        rows.sort(key=lambda x: (x["model_index"], x["terminus"]))
        out.append(locus_verdict(a, rows))
    return out
