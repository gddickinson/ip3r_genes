"""Turn a recovered locus into a reference sequence with tagged junctions.

This is the module every S12 number rests on, so it is worth saying exactly
what the reference sequence *is* and what it is not.

**It is the spliced genomic coding sequence**: the miniprot model's CDS
blocks, taken from the archived sweep GFF, fetched out of the assembly by
coordinate, reverse-complemented on the minus strand, and concatenated in
**target order** (the order the bait's residues align in, which for a minus-
strand gene is the reverse of coordinate order).

It is deliberately **not** S9's frameshift-corrected CDS.  S9 needed a clean
reading frame because codeml counts codons; S12 needs the sequence a read
actually comes from.  Correcting a frameshift inserts or deletes 1-2 bp
relative to the genome, and every read crossing that point would then carry
an indel against the reference — which hisat2 end-to-end penalises, so the
correction would *cost* reads exactly at the sites the sweep already flagged
as difficult.  The uncorrected splice is the honest target.

The reference is validated against the same object S9 validated against —
the sweep's own ``##STA`` protein — but **not** by an identity threshold.
Translating and scoring identity looked obvious and is the wrong
instrument: a frameshift costs the reading frame from where it sits, so a
correct reference scores 1.00 with no frameshifts and 0.92 with nine, and
any floor drawn across that range rejects correct references for carrying
frameshifts the sweep already recorded.

The test used instead is **colinear block placement**, and it has no tuned
threshold.  Each block is translated in its own recorded phase and its
first ``MIN_PLACEMENT_AA`` residues are searched for verbatim in the
expected protein, from the previous block's match onwards.  A block places
only if it is really this protein's sequence *and* it comes after the block
before it.  The pass rule is then a statement the model makes about itself:
**a block may fail to place only if the model's own frameshift count can
explain it.**  Measured against constructed failures, a correct reference
places 100 % of blocks, the same model with its blocks reversed places
1.7 %, and the same model read off the wrong strand places 0 %.

Block count is checked separately, against the ``cds_bp`` the sweep
recorded for this locus, because placement is a fraction of the blocks
present and a dropped block would simply be one fewer test.

**Junctions are derived in target order too**, and each is tagged with what
the annotation holds there — S10's four intron classes, reused unchanged.
`within_one_model` is a junction the annotation models; the other three are
junctions it does not, and those are the ones S10 handed to this task.  The
tag is matched to the junction by **genomic coordinate**, never by index:
a coordinate-sorted list and a target-sorted list have the same length and
opposite order on a minus-strand gene, so an index-to-index mapping is right
half the time and silently wrong the other half.
"""

from __future__ import annotations

import difflib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s10_evidence import Region, classify_introns  # noqa: E402
from s10_gff import read_annotation_window, read_miniprot_model  # noqa: E402
from s12_lib import revcomp, translate  # noqa: E402

MIN_PLACEMENT_AA = 10   # verbatim residues a block must contribute to place

# The four S10 intron classes, and which of them the annotation models.
ANNOTATED_CLASSES = {"within_one_model"}
UNANNOTATED_CLASSES = {"between_models", "model_to_gap", "unannotated"}


class LocusError(RuntimeError):
    """The reference for this locus could not be built or did not validate.

    Raised rather than returning a partial sequence: a reference that is
    silently short by one exon still aligns reads, and every count taken
    against it would be wrong in a way no downstream check could see.
    """


def block_sequence(region, start: int, end: int, strand: str) -> str:
    """The genomic sequence of one CDS block, in translation orientation.

    `region.at()` returns "" for anything outside the fetched window rather
    than raising, so a block the region does not cover would otherwise
    contribute nothing and shorten the reference silently.
    """
    seq = region.at(start, end).upper()
    if len(seq) != end - start + 1:
        raise LocusError(f"block {start}-{end} outside fetched region")
    return revcomp(seq) if strand == "-" else seq


def splice(model: dict, region) -> tuple[str, list[dict]]:
    """(spliced sequence, per-block rows) in target order."""
    blocks = model["cds"]
    if not blocks:
        raise LocusError(f"{model.get('mp_id')}: model has no CDS blocks")
    parts: list[str] = []
    rows: list[dict] = []
    offset = 0
    for i, b in enumerate(blocks, start=1):
        seq = block_sequence(region, b["start"], b["end"], model["strand"])
        parts.append(seq)
        rows.append({"exon_index": i, "start": b["start"], "end": b["end"],
                     "length": len(seq), "phase": b["phase"],
                     "q_start": b["q_start"], "q_end": b["q_end"],
                     "cds_start": offset + 1, "cds_end": offset + len(seq)})
        offset += len(seq)
    return "".join(parts), rows


def frame_aware_translation(model: dict, region) -> tuple[str, list[int]]:
    """Translate the spliced blocks using the GFF `phase` column as the frame.

    Plain concatenation is right only while the reading frame runs
    continuously through the gene, and miniprot reports these loci with
    frameshifts: the sweep's own ``##STA`` protein skips them, so from the
    first frameshift onwards a frame-0 translation of the concatenated
    blocks is noise.  Measured on *Nibea albiflora* ITPR1 that is 1,213
    correct residues followed by 1,538 wrong ones, at an overall identity of
    0.467 — a number low enough to look like the wrong gene and high enough
    that a floor set by eye might have passed it.

    GFF3 `phase` is the authority: it is the number of bases at the start of
    a feature that finish the codon begun in the previous one.  So a codon
    is completed from the carry only when the carry and the phase add to 3;
    when they do not, the frame has jumped and the carry is dropped.  That
    costs one residue per frameshift and keeps the rest of the gene in
    frame.  Returns the protein and the block indices where the frame broke.
    """
    prot: list[str] = []
    carry = ""
    breaks: list[int] = []
    for i, b in enumerate(model["cds"]):
        seq = block_sequence(region, b["start"], b["end"], model["strand"])
        phase = b["phase"] % 3
        if carry:
            if len(carry) + phase == 3:
                prot.append(translate(carry + seq[:phase]))
            else:
                breaks.append(i)
        elif phase and i:
            breaks.append(i)
        body = seq[phase:]
        end = len(body) - len(body) % 3
        prot.append(translate(body[:end]))
        carry = body[end:]
    return "".join(prot), breaks


def translation_identity(model: dict, region, expected_prot: str
                         ) -> tuple[float, int, list[int]]:
    """Frame-aware translation vs the sweep's own protein, as a diagnostic.

    Reported, never gated on: it falls with the model's frameshift count
    (1.00 at none, 0.92 at nine) and a floor across that range would reject
    correct references.  `autojunk` is off deliberately — it is on by
    default and treats any element in more than 1 % of a sequence over 200
    long as junk, which in a protein is most of the amino acid alphabet.
    """
    obs, breaks = frame_aware_translation(model, region)
    obs_n = obs.replace("*", "X")
    exp = expected_prot.replace("*", "X")
    if not obs_n or not exp:
        return 0.0, 0, breaks
    sm = difflib.SequenceMatcher(None, obs_n, exp, autojunk=False)
    matched = sum(b.size for b in sm.get_matching_blocks())
    return 2 * matched / (len(obs_n) + len(exp)), matched, breaks


def colinear_placement(model: dict, region, expected_prot: str,
                       minlen: int = MIN_PLACEMENT_AA) -> tuple[int, int]:
    """(blocks placed, blocks tested) — the reference's validation.

    A block is translated from its own phase and placed by searching the
    expected protein for `minlen` of its residues **from the previous
    block's match onwards**, so order is part of the test rather than a
    separate check.  Blocks too short to carry `minlen` residues are not
    tested; a gene of nothing but micro-exons would therefore be untested
    and the caller refuses on `tested == 0`.
    """
    exp = expected_prot.replace("*", "X")
    cursor = 0
    placed = tested = 0
    for b in model["cds"]:
        seq = block_sequence(region, b["start"], b["end"], model["strand"])
        body = seq[b["phase"] % 3:]
        tr = translate(body[:len(body) - len(body) % 3]).replace("*", "X")
        if len(tr) < minlen:
            continue
        tested += 1
        for i in range(0, len(tr) - minlen + 1):
            j = exp.find(tr[i:i + minlen], max(0, cursor - minlen))
            if j >= 0:
                placed += 1
                cursor = max(cursor, j - i + len(tr))
                break
    return placed, tested


def junctions(model: dict, block_rows: list[dict], intron_rows: list[dict]
              ) -> list[dict]:
    """One row per junction, in CDS coordinates, tagged by annotation class.

    The junction between target exons k and k+1 sits at CDS offset
    `block_rows[k]["cds_end"]`.  Its genomic gap is looked up among S10's
    classified introns **by coordinate**; a junction whose gap matches none
    of them (a model whose blocks are not colinear on the genome) is kept
    with class `non_colinear` rather than dropped, because a dropped
    junction is indistinguishable from a spliced one that nothing crossed.
    """
    by_coord = {(r["start"], r["end"]): r for r in intron_rows}
    strand = model["strand"]
    out: list[dict] = []
    for k in range(len(block_rows) - 1):
        a, b = block_rows[k], block_rows[k + 1]
        if strand == "-":
            lo, hi = b["end"] + 1, a["start"] - 1
        else:
            lo, hi = a["end"] + 1, b["start"] - 1
        info = by_coord.get((lo, hi))
        if info is None or hi < lo:
            klass, donor, acceptor, splice_class = "non_colinear", "", "", ""
            left = right = ""
        else:
            klass = info["class"]
            donor, acceptor = info["donor"], info["acceptor"]
            splice_class = info["splice_class"]
            left, right = info["left_models"], info["right_models"]
        out.append({
            "junction_index": k + 1,
            "cds_offset": a["cds_end"],
            "left_exon": a["exon_index"], "right_exon": b["exon_index"],
            "intron_start": lo, "intron_end": hi,
            "intron_length": max(0, hi - lo + 1),
            "class": klass,
            "annotated": int(klass in ANNOTATED_CLASSES),
            "donor": donor, "acceptor": acceptor,
            "splice_class": splice_class,
            "left_models": left, "right_models": right,
        })
    return out


def validation_controls(model: dict, region, expected_prot: str) -> dict:
    """Placement under the model as built and under two broken versions.

    The pass rule for a reference is that its blocks place colinearly in
    the sweep's own protein, and a rule is only worth reporting beside the
    failures it rejects.  Both controls are constructed from **this** model
    rather than from a synthetic gene, so the numbers in the report are
    measurements on the actual references: blocks reversed (target order
    destroyed, coordinates and strand intact) and the strand flipped
    (coordinates and order intact, sequence read the wrong way).
    """
    import copy
    out = {}
    for tag, mutate in (
            ("as_built", lambda m: m),
            ("reversed_order", lambda m: {**m, "cds": list(reversed(m["cds"]))}),
            ("wrong_strand", lambda m: {**m,
                                        "strand": "+" if m["strand"] == "-"
                                        else "-"})):
        placed, tested = colinear_placement(mutate(copy.deepcopy(model)),
                                            region, expected_prot)
        out[tag] = (placed, tested)
    return out


def build(genome_fna: Path, gff_gz: Path | None, sweep_gff: Path,
          mp_id: str, expected_prot: str, expected_cds_bp: int = 0,
          pad: int = 200) -> dict:
    """Everything S12 needs from one locus, or a LocusError saying why not."""
    model = read_miniprot_model(sweep_gff, mp_id)
    if not model.get("cds"):
        raise LocusError(f"{mp_id}: not found in {sweep_gff}")
    lo = min(b["start"] for b in model["cds"]) - pad
    hi = max(b["end"] for b in model["cds"]) + pad
    region = Region(genome_fna, model["contig"], lo, hi)
    if not region.seq:
        raise LocusError(f"{mp_id}: could not fetch "
                         f"{model['contig']}:{lo}-{hi}")

    seq, block_rows = splice(model, region)
    ident, n_aa, frame_breaks = translation_identity(
        model, region, expected_prot)
    placed, tested = colinear_placement(model, region, expected_prot)
    unplaced = tested - placed
    if tested == 0:
        raise LocusError(f"{mp_id}: no block long enough to place")
    if unplaced > model.get("frameshifts", 0):
        raise LocusError(
            f"{mp_id}: {unplaced} of {tested} blocks do not place colinearly "
            f"in the sweep's own protein, and the model records only "
            f"{model.get('frameshifts', 0)} frameshift(s) to explain them")
    if expected_cds_bp and len(seq) != expected_cds_bp:
        raise LocusError(
            f"{mp_id}: spliced {len(seq)} nt against the {expected_cds_bp} nt "
            f"the sweep recorded for this locus")

    controls = validation_controls(model, region, expected_prot)

    ann = read_annotation_window(gff_gz, model["contig"], lo, hi) \
        if gff_gz and Path(gff_gz).exists() else {"genes": []}
    genes = ann.get("genes", [])
    intron_rows = classify_introns(model, genes, region)
    junc = junctions(model, block_rows, intron_rows)

    return {
        "mp_id": mp_id, "contig": model["contig"], "strand": model["strand"],
        "start": model["start"], "end": model["end"],
        "bait": model.get("bait", ""),
        "frameshifts": model.get("frameshifts", 0),
        "stop_codons": model.get("stop_codons", 0),
        "sequence": seq, "length": len(seq),
        "n_exons": len(block_rows), "exons": block_rows,
        "junctions": junc,
        "n_junctions": len(junc),
        "n_junctions_annotated": sum(j["annotated"] for j in junc),
        "n_junctions_unannotated": sum(
            1 for j in junc if j["class"] in UNANNOTATED_CLASSES),
        "translation_identity": round(ident, 4),
        "translation_aa": n_aa,
        "frame_breaks": len(frame_breaks),
        "blocks_placed": placed, "blocks_tested": tested,
        "placement": round(placed / tested, 4),
        "controls": controls,
        "n_annotated_genes": len(genes),
    }
