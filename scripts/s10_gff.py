"""S10's data layer — annotated features, aligned exons, and interval algebra.

Two GFFs meet in this task and they say different things, so they are parsed
separately and never merged into one feature list:

* the **assembly's own annotation** (`genomic.gff.gz`, downloaded with the
  genome and still on the drive), which is the thing being audited; and
* the **sweep's miniprot alignment** (`miniprot.gff`, archived per genome),
  which is the evidence it is being audited against.

Both are read out of the archive rather than regenerated, so every S10 number
re-derives offline from what S5 actually ran (the discipline `s2_interpro.py`
established and `s3_run_sweep.py --parse-only` applied to HMMER).

**Pseudogenes are kept.** A GFF3 `pseudogene` may still carry a complete set of
`CDS` features — *Nibea albiflora*'s ITPR locus carries 51 — and dropping them
would make S10 report "no coding model here" about a locus that has a coding
model the submitter demoted. The biotype is recorded per feature instead, and
whether the demotion is justified is a question the evidence answers rather
than one the parser decides.

Coordinates throughout are **1-based inclusive**, matching GFF3 and S5's locus
records, and blocks are always returned sorted and non-overlapping.
"""

from __future__ import annotations

import gzip
from pathlib import Path
from urllib.parse import unquote

Block = tuple[int, int]

CODING_TYPES = {"CDS"}
EXON_TYPES = {"exon"}
GENE_TYPES = {"gene", "pseudogene"}
RNA_SUFFIX = ("RNA", "transcript", "_segment")


# --------------------------------------------------------------------------
# interval algebra
# --------------------------------------------------------------------------
def merge(blocks: list[Block]) -> list[Block]:
    """Sorted, non-overlapping union. Touching blocks (b, b+1) merge."""
    out: list[Block] = []
    for s, e in sorted(blocks):
        if out and s <= out[-1][1] + 1:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def overlap_bp(a: list[Block], b: list[Block]) -> int:
    """Total bp shared by two block lists. Both are merged first."""
    a, b = merge(a), merge(b)
    i = j = total = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if hi >= lo:
            total += hi - lo + 1
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def subtract(a: list[Block], b: list[Block]) -> list[Block]:
    """`a` minus `b`, as merged blocks."""
    out: list[Block] = []
    b = merge(b)
    for s, e in merge(a):
        cur = s
        for bs, be in b:
            if be < cur or bs > e:
                continue
            if bs > cur:
                out.append((cur, min(e, bs - 1)))
            cur = max(cur, be + 1)
            if cur > e:
                break
        if cur <= e:
            out.append((cur, e))
    return out


def span_bp(blocks: list[Block]) -> int:
    return sum(e - s + 1 for s, e in merge(blocks))


def gaps_between(blocks: list[Block]) -> list[Block]:
    """The introns implied by a sorted exon list."""
    m = merge(blocks)
    return [(m[i][1] + 1, m[i + 1][0] - 1)
            for i in range(len(m) - 1) if m[i + 1][0] > m[i][1] + 1]


def n_disjoint(blocks: list[Block]) -> int:
    """How many disjoint pieces a block list collapses to.

    The brief asks for *disjoint coding model blocks, not genes*, and this is
    why: two annotated genes 359 bp apart are two genes but, once their CDS
    blocks are merged, the coding evidence is still two pieces — while one
    gene whose CDS is interrupted by a 20 kb unmodelled stretch is one gene
    and two pieces. Counting genes measures the annotation's bookkeeping;
    counting disjoint blocks measures what a reader can actually recover.
    """
    return len(merge(blocks))


# --------------------------------------------------------------------------
# the assembly's annotation
# --------------------------------------------------------------------------
def _attrs(field: str) -> dict[str, str]:
    out = {}
    for part in field.rstrip(";").split(";"):
        if "=" in part:
            k, _, v = part.partition("=")
            out[k.strip()] = unquote(v)
    return out


def read_annotation_window(gff_gz: Path, contig: str, start: int, end: int
                           ) -> dict:
    """Every annotated feature touching `contig:start-end`.

    Returns ``{"genes": [...], "by_gene": {gene_id: {...}}}`` where each gene
    carries its span, name, biotype, and the merged CDS and exon blocks of all
    its transcripts. Transcript-level detail is kept in `mrnas` because the
    tiling step needs one protein per transcript, not per gene.

    The whole file is streamed once and filtered on the contig, because a
    vertebrate GFF3 is 5–8 MB gzipped and random access into it would need an
    index this project does not build. Parent links are resolved in a second
    pass over the retained rows, since GFF3 does not guarantee that a parent
    is written before its child.
    """
    rows = []
    op = gzip.open if str(gff_gz).endswith(".gz") else open
    with op(gff_gz, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != contig:
                continue
            s, e = int(f[3]), int(f[4])
            if e < start or s > end:
                continue
            rows.append((f[2], s, e, f[6], _attrs(f[8]),
                         int(f[7]) if f[7] not in (".", "") else 0))

    rows_phase = rows
    rows = [r[:5] for r in rows_phase]
    genes, mrna_parent, by_gene = {}, {}, {}
    for kind, s, e, strand, at in rows:
        if kind in GENE_TYPES:
            gid = at.get("ID", "")
            genes[gid] = {
                "gene_id": gid, "start": s, "end": e, "strand": strand,
                "name": at.get("Name") or at.get("gene") or
                        at.get("locus_tag") or gid,
                "locus_tag": at.get("locus_tag", ""),
                "biotype": at.get("gene_biotype", kind),
                "pseudo": at.get("pseudo", "") == "true" or kind == "pseudogene",
                "description": at.get("description", ""),
                "cds": [], "exons": [], "mrnas": {},
            }
    for kind, s, e, strand, at in rows:
        if kind.endswith(RNA_SUFFIX) or kind == "mRNA":
            rid, par = at.get("ID", ""), at.get("Parent", "")
            if rid:
                mrna_parent[rid] = par
                if par in genes:
                    genes[par]["mrnas"].setdefault(
                        rid, {"cds": [], "exons": [], "cds_phase": {},
                              "protein_id": "",
                              "product": at.get("product", "")})
    for kind, s, e, strand, at, phase in rows_phase:
        if kind not in CODING_TYPES and kind not in EXON_TYPES:
            continue
        par = at.get("Parent", "")
        gid = mrna_parent.get(par, par if par in genes else "")
        if gid not in genes:
            continue
        key = "cds" if kind in CODING_TYPES else "exons"
        genes[gid][key].append((s, e))
        tx = genes[gid]["mrnas"].get(par)
        if tx is not None:
            tx[key].append((s, e))
            if kind in CODING_TYPES:
                # Phase is kept per block, keyed on the block, because a
                # partial model's first coding block is often not the start of
                # a codon and translating it at phase 0 frameshifts exactly the
                # fragments this task is about.
                tx.setdefault("cds_phase", {})[(s, e)] = phase
                if at.get("protein_id"):
                    tx["protein_id"] = at["protein_id"]
    for g in genes.values():
        g["cds"], g["exons"] = merge(g["cds"]), merge(g["exons"])
        for tx in g["mrnas"].values():
            tx["cds"], tx["exons"] = merge(tx["cds"]), merge(tx["exons"])
        by_gene[g["gene_id"]] = g
    return {"genes": sorted(genes.values(), key=lambda g: g["start"]),
            "by_gene": by_gene}


def annotation_meta(gff_gz: Path) -> dict:
    """The GFF3 header: which build the annotation was written against.

    This is the assembly-version audit's primary evidence, and it is the one
    place a version mismatch would be visible without guessing — NCBI's
    annotwriter records `genome-build` and `genome-build-accession` in the
    header, so an annotation transferred from an earlier build says so.
    """
    out: dict[str, str] = {}
    op = gzip.open if str(gff_gz).endswith(".gz") else open
    with op(gff_gz, "rt") as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            # NCBI writes these as `#!genome-build ...`, so both the comment
            # marker and the pragma bang have to come off. Stripping only `#`
            # leaves `!genome-build` and every field reads as absent, which
            # looks exactly like an annotation with no build stamp.
            t = line.strip().lstrip("#!").strip()
            for key, field in (("genome-build-accession", "build_accession"),
                               ("genome-build", "build_name"),
                               ("processor", "processor"),
                               ("gff-spec-version", "gff_spec")):
                if t.startswith(key) and field not in out:
                    out[field] = t[len(key):].strip(" :")
    return out


# --------------------------------------------------------------------------
# the sweep's alignment
# --------------------------------------------------------------------------
def read_miniprot_model(gff: Path, mp_id: str) -> dict:
    """One miniprot model's mRNA record and its CDS blocks, in order.

    Blocks are returned in **target order** (the order the bait's residues are
    aligned in), not coordinate order, so a minus-strand model's first block is
    the one carrying residue 1. Every downstream step — exon numbering, the
    junction probes, the reconstructed CDS — depends on that, and sorting by
    coordinate instead silently reverses every minus-strand gene.
    """
    model: dict = {"mp_id": mp_id, "cds": []}
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            at = _attrs(f[8].replace("Target=", "Target="))
            if f[2] == "mRNA" and at.get("ID") == mp_id:
                tgt = at.get("Target", "").split()
                model.update(
                    contig=f[0], start=int(f[3]), end=int(f[4]),
                    strand=f[6], score=float(f[5]) if f[5] != "." else 0.0,
                    bait=tgt[0] if tgt else "",
                    identity=float(at.get("Identity", 0)),
                    frameshifts=int(at.get("Frameshift", 0)),
                    stop_codons=int(at.get("StopCodon", 0)))
            elif f[2] == "CDS" and at.get("Parent") == mp_id:
                tgt = at.get("Target", "").split()
                model["cds"].append({
                    "start": int(f[3]), "end": int(f[4]), "strand": f[6],
                    "phase": int(f[7]) if f[7] != "." else 0,
                    "identity": float(at.get("Identity", 0)),
                    "q_start": int(tgt[1]) if len(tgt) > 2 else 0,
                    "q_end": int(tgt[2]) if len(tgt) > 2 else 0})
    model["cds"].sort(key=lambda b: b["q_start"])
    return model


def model_blocks(model: dict) -> list[Block]:
    return merge([(b["start"], b["end"]) for b in model["cds"]])
