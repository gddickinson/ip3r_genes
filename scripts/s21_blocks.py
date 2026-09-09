"""The exon blocks, read out of the sweep's own GFFs — S21's instrument.

A `miniprot.gff` holds one mRNA record per alignment and one CDS record per
aligned block, each carrying its genomic interval, its phase and its span in
the **bait's** residue coordinates. That is a spliced gene model, and the
whole of S21 is read off it. Four things this module gets right and a naive
reader does not.

**The keys.** miniprot numbers alignments from `MP000001` per run, so the two
giant assemblies' chunked GFFs repeat every id. `s5_sweep_lib.parse_miniprot_gff`
disambiguates a repeat as `MP000001#<n>` where n is the number of alignments
seen so far, and the sweep's `summary.json` records *that* key. So this reader
reproduces the same keying exactly rather than inventing its own, and
`s21_test_arch` requires the two readers to return identical genomic blocks for
every model in a genome — a mismatch would silently move every exon.

**Gene order.** Blocks are returned in target order, so a minus-strand gene's
first block is the one carrying bait residue 1. An intron is the gap between
consecutive blocks *in that order*, which on the minus strand runs downwards
in genomic coordinate; computing it as `next.start - prev.end` reverses every
minus-strand intron and reports a negative length or, worse, a plausible one
from the wrong pair.

**Frameshifts are inside blocks, not between them.** The brief expected
miniprot to emit two CDS records either side of a frameshift, which would make
a broken locus look exon-rich. Measured over all 309 genomes it does not:
query spans are contiguous across every consecutive pair (no residue is
emitted twice, no residue skipped) and an indel shows up as a block whose
genomic length differs from three times its residue count. So the merge S21
needs is not a frameshift merge but a **sub-intron gap** merge, and its
threshold is measured against splice dinucleotides in
`s21_calibrate_intron.py` rather than chosen.

**A block is not an exon.** An exon is a maximal run of blocks not separated
by a real intron. `exons()` applies the calibrated gap rule and records how
many merges it made, so a locus whose exon count depended on the rule is
visible in the data.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_lib as L                                               # noqa: E402

#: Sub-intron gap floor, in bp. Two consecutive blocks closer than this are one
#: exon interrupted by an indel, not two exons separated by an intron. The
#: value is the calibration's, read back from the committed table by
#: `min_intron_bp()`; this is only the fallback used by the calibration itself
#: and by the negative controls, and it is miniprot's own minimum intron.
MIN_INTRON_FALLBACK = 30


def _attr(field: str, key: str) -> str:
    for part in field.split(";"):
        if part.startswith(key + "="):
            return part[len(key) + 1:]
    return ""


def read_models(gff: Path, wanted: set[str] | None = None) -> dict[str, dict]:
    """Every alignment in one GFF as `{key: model}`, blocks in target order.

    `wanted` restricts the returned models but **not** the pass: the keying
    depends on how many alignments have been seen, so every mRNA record has to
    be counted whether or not it is kept.
    """
    models: dict[str, dict] = {}
    seen: set[str] = set()
    n_aln = 0
    cur_raw = cur_key = ""
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            if f[2] == "mRNA":
                raw = _attr(f[8], "ID")
                key = raw if raw not in seen else f"{raw}#{n_aln}"
                seen.add(raw)
                cur_raw, cur_key = raw, key
                n_aln += 1
                if wanted is not None and key not in wanted:
                    continue
                tgt = _attr(f[8], "Target").split()
                models[key] = {
                    "mp_id": key, "contig": f[0], "start": int(f[3]),
                    "end": int(f[4]), "strand": f[6],
                    "score": float(f[5]) if f[5] not in (".", "") else 0.0,
                    "bait": tgt[0] if tgt else "",
                    "identity": float(_attr(f[8], "Identity") or 0.0),
                    "frameshifts": int(_attr(f[8], "Frameshift") or 0),
                    "stop_codons": int(_attr(f[8], "StopCodon") or 0),
                    "blocks": []}
            elif f[2] == "CDS":
                parent = _attr(f[8], "Parent")
                key = cur_key if parent == cur_raw else parent
                m = models.get(key)
                if m is None:
                    continue
                tgt = _attr(f[8], "Target").split()
                m["blocks"].append({
                    "start": int(f[3]), "end": int(f[4]),
                    "phase": int(f[7]) if f[7] != "." else 0,
                    "identity": float(_attr(f[8], "Identity") or 0.0),
                    "q_start": int(tgt[1]) if len(tgt) > 2 else 0,
                    "q_end": int(tgt[2]) if len(tgt) > 2 else 0})
    for m in models.values():
        m["blocks"].sort(key=lambda b: b["q_start"])
    return models


# --------------------------------------------------------------------------
# gene order, introns, exons
# --------------------------------------------------------------------------
def gap_bp(prev: dict, nxt: dict, strand: str) -> int:
    """The genomic gap between two consecutive blocks, in **gene order**.

    Negative when the two blocks overlap, which miniprot produces at an
    insertion. Doing this without the strand branch reverses every
    minus-strand intron.
    """
    if strand == "-":
        return prev["start"] - nxt["end"] - 1
    return nxt["start"] - prev["end"] - 1


def intron_interval(prev: dict, nxt: dict, strand: str) -> tuple[int, int]:
    """The intron's genomic interval (low, high), 1-based inclusive."""
    if strand == "-":
        return (nxt["end"] + 1, prev["start"] - 1)
    return (prev["end"] + 1, nxt["start"] - 1)


def block_frame_deviation(b: dict) -> int:
    """Genomic length minus three times the residues the block carries.

    Zero to within the phase slop (-2..+2) for a clean block: a block's two
    ends usually share a codon with its neighbours. The terminal block of a
    complete model is +3, because miniprot writes the stop codon into it.
    """
    return (b["end"] - b["start"] + 1) - 3 * (b["q_end"] - b["q_start"] + 1)


def next_phase(b: dict) -> int:
    """The phase the block after this one must carry if the frame is intact.

    GFF phase is the number of bases to drop from the start of the feature to
    reach the first base of a codon, so a block of length L opening at phase p
    leaves `(L - p) % 3` bases of a codon hanging over its end, and the next
    block must open at the complement of that. A block pair that disagrees is
    a frame step — which is how *this* output records a frameshift, since the
    query spans are contiguous and no second CDS record is emitted.
    """
    length = b["end"] - b["start"] + 1
    return (3 - ((length - b["phase"]) % 3)) % 3


def frame_steps(blocks: list[dict]) -> int:
    """Junctions at which the reading frame does not carry over.

    A lower bound on miniprot's own `Frameshift` count, not a substitute for
    it: two frameshifts inside one block can cancel at the junction, and one
    in the terminal block has no junction after it. Measured over the whole
    sweep the two agree on 96 % of models and this count is never the larger,
    which is what makes it usable as the per-**boundary** trust flag S21 needs
    — an exon boundary adjacent to a frame step is one the alignment is not
    entitled to call.
    """
    return sum(1 for a, b in zip(blocks, blocks[1:])
               if next_phase(a) != b["phase"])


def exons(model: dict, min_intron: int) -> tuple[list[dict], int]:
    """Blocks merged into exons, plus the number of merges made.

    A gap below `min_intron` is not an intron, so the two blocks either side of
    it are one exon. The merged exon keeps the union of the genomic interval,
    the union of the residue span and the **first** block's phase (the phase of
    an exon is the phase at which it starts), and records the smallest gap it
    swallowed so a locus whose exon count depended on the rule is visible in
    the data.
    """
    bl = model["blocks"]
    if not bl:
        return [], 0
    strand = model["strand"]
    out: list[dict] = []
    merges = 0
    cur = dict(bl[0], n_blocks=1, min_gap_bp=-1)
    for b in bl[1:]:
        g = gap_bp({"start": cur["start"], "end": cur["end"]}, b, strand)
        if g < min_intron:
            merges += 1
            cur["start"] = min(cur["start"], b["start"])
            cur["end"] = max(cur["end"], b["end"])
            cur["q_end"] = max(cur["q_end"], b["q_end"])
            cur["n_blocks"] += 1
            cur["min_gap_bp"] = g if cur["min_gap_bp"] < 0 \
                else min(cur["min_gap_bp"], g)
        else:
            out.append(cur)
            cur = dict(b, n_blocks=1, min_gap_bp=-1)
    out.append(cur)
    return out, merges


def introns(exon_list: list[dict], strand: str) -> list[dict]:
    """One row per intron: gene-order index, genomic interval, length."""
    out = []
    for i, (p, n) in enumerate(zip(exon_list, exon_list[1:])):
        lo, hi = intron_interval(p, n, strand)
        out.append({"index": i, "start": lo, "end": hi,
                    "length": hi - lo + 1,
                    "q_after": p["q_end"], "q_before": n["q_start"],
                    "phase_next": n["phase"]})
    return out


def architecture(model: dict, min_intron: int) -> dict:
    """The per-locus architecture numbers, all from the blocks.

    `cds_bp` is the union of the exon intervals rather than their sum, because
    an insertion makes two blocks overlap and summing would count the overlap
    twice. `span_bp` is the model's own outer span, which is what a database
    reports as the gene's length. `frame_balance_bp` is the coding length the
    model carries beyond three bases per aligned residue, with the terminal
    stop codon removed — the net indel load, and the number that says whether
    the exon structure is being read off a model in frame with its bait.
    """
    ex, merges = exons(model, min_intron)
    intr = introns(ex, model["strand"])
    spans = sorted((e["start"], e["end"]) for e in ex)
    cds_bp, cs, ce = 0, *spans[0]
    for s_, e_ in spans[1:]:
        if s_ <= ce + 1:
            ce = max(ce, e_)
        else:
            cds_bp += ce - cs + 1
            cs, ce = s_, e_
    cds_bp += ce - cs + 1
    q_lo = min(e["q_start"] for e in ex)
    q_hi = max(e["q_end"] for e in ex)
    residues = q_hi - q_lo + 1
    ilen = [i["length"] for i in intr]
    stop_bp = 3 if block_frame_deviation(model["blocks"][-1]) == 3 else 0
    return {
        "n_blocks": len(model["blocks"]), "n_exons": len(ex),
        "n_merges": merges, "n_introns": len(intr),
        "frame_steps": frame_steps(model["blocks"]),
        "cds_bp": cds_bp, "span_bp": model["end"] - model["start"] + 1,
        "q_start": q_lo, "q_end": q_hi, "residues": residues,
        "total_intron_bp": sum(ilen),
        "max_intron_bp": max(ilen) if ilen else 0,
        "median_intron_bp": L.quantile([float(x) for x in ilen], 0.5)
        if ilen else 0.0,
        "min_intron_bp_seen": min(ilen) if ilen else 0,
        "stop_codon_bp": stop_bp,
        "frame_balance_bp": cds_bp - stop_bp - 3 * residues,
        "mean_exon_bp": round(cds_bp / len(ex), 2) if ex else 0.0,
        "exons": ex, "introns": intr}


def min_intron_bp() -> int:
    """The calibrated sub-intron gap floor, read back from its own table."""
    path = L.OUT / "intron_floor_calibration.tsv"
    if path.exists():
        rows = L.read_tsv(path)
        for r in rows:
            if r.get("chosen") == "1":
                return int(float(r["floor_bp"]))
    return MIN_INTRON_FALLBACK
