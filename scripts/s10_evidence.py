"""S10 step 2 — what the annotation actually put on the gene, exon by exon.

The claim S5b makes about a case locus is a claim about *coverage*: some
fraction of a real gene's coding footprint has no gene model on it. That is a
summary statistic, and a summary statistic is not molecular validation. This
module turns it into per-feature evidence:

* **every aligned exon** classified by what the annotation has at that exact
  interval — annotated coding sequence, annotated non-coding sequence (a UTR or
  an intron of some gene), or nothing at all;
* **every intron** classified by *which annotated feature its two flanking
  exons fall in*, which is what distinguishes the three shapes a broken gene
  can take: an intron inside one model (the annotation is intact there), an
  intron *between two models* (the annotation split the gene at that point),
  and an intron running into unannotated sequence;
* **disjoint coding blocks**, counted rather than genes, on both sides;
* **splice-site dinucleotides** read off the genome at every intron, because
  the exon structure being audited is one an aligner proposed and a reader is
  entitled to ask whether it is spliceable.

The last of those is a check on S10's own evidence rather than on the
annotation, and it is the reason this module reads the genome at all. An
alignment whose introns do not begin GT and end AG is an alignment, not a gene,
and reporting it as a missed gene would be the same error the task exists to
catch — one instrument's output taken as fact because nothing tested it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s10_gff import (Block, gaps_between, merge, n_disjoint,  # noqa: E402
                     overlap_bp, read_annotation_window, read_miniprot_model,
                     subtract)
from s5_genome_io import build_fai, fetch_region, read_fai  # noqa: E402

#: How far either side of the locus the annotation is read. The neighbourhood
#: has to be wide enough to name the genes that bracket the gap — for an
#: omission those flanking genes are the evidence that the annotation was
#: working right up to the missing gene — but not so wide that a distant gene
#: is offered as an explanation. 200 kb is ~1.5x the widest ITPR intron
#: measured (152 kb, `intron_calibration.tsv`).
FLANK_BP = 200_000

#: Canonical and the common non-canonical splice pair. U12 introns are rare
#: (~0.3 % of human introns) and real, so an AT..AC intron is recorded as
#: `minor` rather than as a failure.
SPLICE_CANONICAL = {("GT", "AG")}
SPLICE_MINOR = {("AT", "AC"), ("GC", "AG")}

EXON_CLASSES = ("in_annotated_cds", "in_annotated_noncoding", "intergenic")
INTRON_CLASSES = ("within_one_model", "between_models",
                  "model_to_gap", "unannotated")


def _revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))[::-1]


class Region:
    """The genome under a locus plus its flanks, fetched once.

    Coordinates handed in are absolute and 1-based; the class does the offset
    arithmetic, because doing it at each call site is how an off-by-one gets
    into a splice-site table and is never noticed.
    """

    def __init__(self, fna: Path, contig: str, start: int, end: int):
        idx = read_fai(build_fai(fna))
        self.contig, self.start = contig, max(1, start)
        self.seq = fetch_region(fna, idx, contig, self.start, end).upper()
        self.end = self.start + len(self.seq) - 1

    def at(self, start: int, end: int) -> str:
        if start < self.start or end > self.end or end < start:
            return ""
        return self.seq[start - self.start:end - self.start + 1]


def splice_pair(region: Region, intron: Block, strand: str) -> tuple[str, str]:
    """The donor and acceptor dinucleotides of one intron, read in gene order.

    On the minus strand the genomic first two bases of the interval are the
    *acceptor's* reverse complement, so the pair has to be built from the
    revcomp of both ends and swapped. Getting this wrong turns every canonical
    minus-strand intron into `CT..AC` — which looks like a systematic
    non-canonical signal and would falsely discredit the alignment.
    """
    s, e = intron
    left, right = region.at(s, s + 1), region.at(e - 1, e)
    if len(left) < 2 or len(right) < 2:
        return "", ""
    if strand == "-":
        return _revcomp(right), _revcomp(left)
    return left, right


def splice_class(donor: str, acceptor: str) -> str:
    if not donor or not acceptor:
        return "unreadable"
    if (donor, acceptor) in SPLICE_CANONICAL:
        return "canonical"
    if (donor, acceptor) in SPLICE_MINOR:
        return "minor"
    return "non_canonical"


def _gene_at(genes: list[dict], block: Block) -> list[dict]:
    s, e = block
    return [g for g in genes if g["start"] <= e and g["end"] >= s]


def classify_exons(model: dict, genes: list[dict]) -> list[dict]:
    """One row per aligned exon, with what the annotation holds there.

    An exon is scored against annotated **CDS blocks**, not gene spans. A gene
    span covers its own introns, so scoring against spans would credit the
    annotation with coding sequence it never called — which is exactly the
    over-credit `frac_cds` carries in S5's summary and the reason S10 goes back
    to the GFF instead of reusing that number here.
    """
    rows = []
    for i, blk in enumerate(model["cds"], start=1):
        b: Block = (blk["start"], blk["end"])
        length = b[1] - b[0] + 1
        hits = _gene_at(genes, b)
        cds_ov = overlap_bp([b], [x for g in hits for x in g["cds"]])
        gene_ov = overlap_bp([b], [(g["start"], g["end"]) for g in hits])
        if cds_ov > 0:
            klass = "in_annotated_cds"
        elif gene_ov > 0:
            klass = "in_annotated_noncoding"
        else:
            klass = "intergenic"
        best = max(hits, key=lambda g: overlap_bp([b], g["cds"] or
                                                  [(g["start"], g["end"])]),
                   default=None)
        rows.append({
            "exon_index": i, "start": b[0], "end": b[1], "length": length,
            "q_start": blk["q_start"], "q_end": blk["q_end"],
            "aln_identity": blk["identity"], "phase": blk["phase"],
            "class": klass,
            "annot_gene": best["name"] if best else "",
            "annot_gene_id": best["gene_id"] if best else "",
            "annot_biotype": best["biotype"] if best else "",
            "annot_pseudo": int(bool(best and best["pseudo"])),
            "annot_cds_overlap_bp": cds_ov,
            "annot_cds_overlap_frac": round(cds_ov / length, 4),
        })
    return rows


def classify_introns(model: dict, genes: list[dict], region: Region
                     ) -> list[dict]:
    """One row per implied intron, classified by where its flanks land.

    The four classes answer different questions, and only the second is an
    annotation bug in its own right:

    * `within_one_model` — the annotation has both flanking exons in the same
      gene, so the gene model is intact across this intron;
    * `between_models` — the flanking exons are in two *different* annotated
      genes: this is a splice junction the annotation broke the gene at, and
      it is the junction the transcript probes are aimed at;
    * `model_to_gap` — one flank is annotated and the other is not, so this is
      where a model ends and the unannotated part of the gene begins;
    * `unannotated` — neither flank has a model.
    """
    blocks = sorted([(b["start"], b["end"]) for b in model["cds"]])
    rows = []
    for i, intr in enumerate(gaps_between(blocks), start=1):
        left = (blocks[i - 1][0], blocks[i - 1][1])
        right = (blocks[i][0], blocks[i][1])
        lg = {g["gene_id"] for g in _gene_at(genes, left)
              if overlap_bp([left], g["cds"]) > 0}
        rg = {g["gene_id"] for g in _gene_at(genes, right)
              if overlap_bp([right], g["cds"]) > 0}
        if lg and rg:
            klass = "within_one_model" if lg & rg else "between_models"
        elif lg or rg:
            klass = "model_to_gap"
        else:
            klass = "unannotated"
        donor, acceptor = splice_pair(region, intr, model["strand"])
        rows.append({
            "intron_index": i, "start": intr[0], "end": intr[1],
            "length": intr[1] - intr[0] + 1, "class": klass,
            "left_models": ";".join(sorted(lg)), "right_models": ";".join(sorted(rg)),
            "donor": donor, "acceptor": acceptor,
            "splice_class": splice_class(donor, acceptor),
        })
    return rows


def annotated_models(model: dict, genes: list[dict]) -> list[dict]:
    """The annotated genes that carry any of the aligned coding sequence.

    Reported with `pseudo` and `biotype` on every row. A GFF3 `pseudogene` may
    carry a full set of CDS features and still emit no protein, so "how much
    of the gene is annotated" and "how much of it reaches a protein database"
    are different numbers and both are needed.
    """
    aligned: list[Block] = [(b["start"], b["end"]) for b in model["cds"]]
    total = max(1, sum(e - s + 1 for s, e in merge(aligned)))
    rows = []
    for g in genes:
        ov = overlap_bp(aligned, g["cds"])
        if ov <= 0:
            continue
        covered = [b for b in model["cds"]
                   if overlap_bp([(b["start"], b["end"])], g["cds"]) > 0]
        rows.append({
            "gene_id": g["gene_id"], "name": g["name"],
            "locus_tag": g["locus_tag"], "biotype": g["biotype"],
            "pseudo": int(g["pseudo"]), "start": g["start"], "end": g["end"],
            "strand": g["strand"], "span": g["end"] - g["start"] + 1,
            "n_cds_blocks": len(g["cds"]), "cds_bp": sum(e - s + 1 for s, e in g["cds"]),
            "n_transcripts": len(g["mrnas"]),
            "protein_ids": ";".join(sorted(
                {t["protein_id"] for t in g["mrnas"].values() if t["protein_id"]})),
            "overlap_aligned_cds_bp": ov,
            "frac_aligned_cds": round(ov / total, 4),
            "n_aligned_exons_covered": len(covered),
            "q_start": min([b["q_start"] for b in covered], default=0),
            "q_end": max([b["q_end"] for b in covered], default=0),
            "description": g["description"][:200],
        })
    rows.sort(key=lambda r: -r["frac_aligned_cds"])
    return rows


def block_accounting(model: dict, genes: list[dict]) -> dict:
    """Disjoint coding blocks on each side, and what is missing from which.

    `n_annotated_gene_models` and `n_annotated_cds_blocks` differ whenever the
    annotation splits a gene without leaving a gap between the pieces — and
    they differ the other way whenever one model's CDS is interrupted. Both
    numbers are reported because the brief asks for blocks and a reader will
    ask for genes.
    """
    aligned = merge([(b["start"], b["end"]) for b in model["cds"]])
    annot = merge([blk for g in genes for blk in g["cds"]])
    annot_coding = merge([blk for g in genes if not g["pseudo"]
                          for blk in g["cds"]])
    uncovered = subtract(aligned, annot)
    untranslated = subtract(aligned, annot_coding)
    total = max(1, sum(e - s + 1 for s, e in aligned))
    return {
        "n_aligned_exons": len(model["cds"]),
        "n_aligned_cds_blocks": n_disjoint(aligned),
        "aligned_cds_bp": total,
        "n_annotated_gene_models": sum(1 for g in genes if overlap_bp(aligned, g["cds"]) > 0),
        "n_annotated_cds_blocks": n_disjoint(annot),
        "annotated_cds_bp_on_gene": overlap_bp(aligned, annot),
        "frac_aligned_cds_annotated": round(overlap_bp(aligned, annot) / total, 4),
        "frac_aligned_cds_translated": round(
            overlap_bp(aligned, annot_coding) / total, 4),
        "n_uncovered_blocks": len(uncovered),
        "uncovered_bp": sum(e - s + 1 for s, e in uncovered),
        "n_untranslated_blocks": len(untranslated),
        "untranslated_bp": sum(e - s + 1 for s, e in untranslated),
    }


def neighbourhood(genes: list[dict], locus: Block, n_each: int = 5
                  ) -> list[dict]:
    """The nearest annotated genes on each side that do not touch the locus.

    For an omission these are the load-bearing rows: the annotation naming
    genes correctly up to both edges of a gap is what makes the gap a hole in
    the annotation rather than a hole in the assembly.
    """
    s, e = locus
    left = [g for g in genes if g["end"] < s][-n_each:]
    right = [g for g in genes if g["start"] > e][:n_each]
    rows = []
    for side, gs in (("upstream", left), ("downstream", right)):
        for g in gs:
            rows.append({
                "side": side, "gene_id": g["gene_id"], "name": g["name"],
                "biotype": g["biotype"], "pseudo": int(g["pseudo"]),
                "start": g["start"], "end": g["end"], "strand": g["strand"],
                "distance_bp": (s - g["end"]) if side == "upstream"
                               else (g["start"] - e),
                "description": g["description"][:200],
            })
    return rows


def gather(case: dict, genome_dir: Path, sweep_dir: Path) -> dict:
    """All of the above for one case, from the archive plus the genome."""
    gff_gz = genome_dir / "genomic.gff.gz"
    fna = next(genome_dir.glob("*_genomic.fna"))
    model = read_miniprot_model(sweep_dir / "miniprot.gff", case["mp_id"])
    if not model.get("cds"):
        raise SystemExit(f"[s10] no CDS blocks for {case['mp_id']}")
    lo = min(b["start"] for b in model["cds"])
    hi = max(b["end"] for b in model["cds"])
    ann = read_annotation_window(gff_gz, model["contig"],
                                 lo - FLANK_BP, hi + FLANK_BP)
    region = Region(fna, model["contig"], lo - 100, hi + 100)
    inside = [g for g in ann["genes"] if g["start"] <= hi and g["end"] >= lo]
    return {
        "model": model, "genes_all": ann["genes"], "genes_inside": inside,
        "exons": classify_exons(model, inside),
        "introns": classify_introns(model, inside, region),
        "models": annotated_models(model, inside),
        "blocks": block_accounting(model, inside),
        "neighbourhood": neighbourhood(ann["genes"], (lo, hi)),
        "locus": (lo, hi), "region": region,
    }
