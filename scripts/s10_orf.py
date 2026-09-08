"""S10 step 3 — is the gene an open reading frame, or is it really a pseudogene?

Case A's annotation makes a specific, falsifiable claim: the model at the locus
is a `pseudogene`, `pseudo=true`, and therefore emits no protein. That claim is
testable and this module tests it, because "the annotation calls it a
pseudogene and we disagree" is an opinion until the reading frame is counted.

**The test.** Reconstruct the spliced coding sequence at the locus, translate
it, and count what a pseudogene is defined by: premature stop codons and
frame-disrupting indels. A coding sequence 2,700 codons long that has been
evolving without selection carries nonsense mutations — the expected number
under neutrality is computed here rather than asserted, from the locus's own
base composition and the observed synonymous divergence, so the comparison is
against this gene rather than against a rule of thumb.

**Where the CDS comes from.** `s9_miniprot_cds.cds_for_model`, unchanged: it
realigns the locus with `miniprot --aln`, whose `##ATA` row keeps the
frameshift residues that the sweep's `##STA` protein skips, and it **refuses**
a reconstruction that does not reproduce the sweep's own protein. S9 wrapped
twelve negative controls around that function precisely because a
frame-shifted CDS is the one artefact that stays invisible downstream. Reusing
it means S10's reading-frame claim rests on machinery that has already been
made to fail on purpose.

A frameshift reported by miniprot is *not* by itself evidence of pseudogeny at
these divergences: the bait is 13–14 % divergent from the target in both cases,
and a single-base alignment artefact at a short exon boundary produces the same
signal. What separates them is that a real pseudogene accumulates stops as well.
So the verdict rule is stated on stops, with frameshifts reported beside it.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import (GENETIC_CODE, cds_from_miniprot_model,  # noqa: E402
                        translate)

#: The three terminators, read out of the genetic code rather than listed, so
#: this module and `s9_cds_lib.translate` cannot disagree about what a stop is.
STOPS = frozenset(c for c, aa in GENETIC_CODE.items() if aa == "*")

#: A pseudogene of this length is expected to carry stops. The verdict rule is
#: deliberately one-sided: zero internal stops falsifies the pseudogene call,
#: while a handful of stops would not on its own establish it (they could be
#: sequencing error in a single-pass assembly). S10 only needs the direction it
#: has evidence for.
MAX_STOPS_FOR_ORF = 0


def _read_fasta(path: Path) -> dict[str, str]:
    out, name, buf = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name:
                out[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        elif name:
            buf.append(line.strip())
    if name:
        out[name] = "".join(buf)
    return out


def _model_protein(sweep_dir: Path, mp_id: str, cell: str, contig: str,
                   start: int, end: int, strand: str) -> tuple[str, str]:
    """The sweep's own protein for this locus, and where it was read from.

    The sweep GFF's `##STA` row is the authority: it is the protein
    `s9_miniprot_cds.find_model` matches on, and it exists for **every** model.
    `novel_models.faa` is only the fallback, because it holds the
    *non-redundant* census-v4 set — a genome with two loci of one paralog
    contributes one of them, and the other's reading frame would then read as
    `no_model_protein` when the evidence for it is sitting in the archive.
    (Measured: 3 of the 8 family loci across the two case genomes.)

    The fallback matches on the **locus**, not on the cell name: a genome can
    hold more than one model per cell and taking the first that names the right
    paralog would silently analyse a different gene.
    """
    gff = sweep_dir / "miniprot.gff"
    if mp_id and gff.exists():
        from s9_miniprot_cds import index_models
        for m in index_models(gff):
            if m["id"] == mp_id and m["protein"]:
                return m["protein"].replace("*", "X"), "sweep_gff"
    want = f"{contig}:{start}-{end}{strand}"
    novel = sweep_dir / "novel_models.faa"
    if novel.exists():
        for header, seq in _read_fasta(novel).items():
            parts = header.split("|")
            if len(parts) > 2 and parts[1] == cell and parts[2] == want:
                return seq, "novel_models"
    return "", "none"


def _model_span(sweep_dir: Path, mp_id: str) -> tuple[int, int]:
    gff = sweep_dir / "miniprot.gff"
    if not (mp_id and gff.exists()):
        return 0, 0
    from s9_miniprot_cds import index_models
    for m in index_models(gff):
        if m["id"] == mp_id:
            return m["start"], m["end"]
    return 0, 0


def expected_stops(cds: str, protein: str, identity: float) -> dict:
    """Stops expected if this CDS had been evolving neutrally since duplication.

    A crude but honest calculation, and it is written out rather than quoted:
    take the observed protein divergence `1 - identity` as a floor on the
    number of substitutions per codon, and multiply by the per-substitution
    probability of creating a stop, computed from the locus's *own* codon
    usage by enumerating every single-base change of every codon in it. The
    result is a lower bound on the expected count, because protein identity
    counts only the changes that were not lethal to the reading frame.
    """
    codons = [cds[i:i + 3] for i in range(0, len(cds) - 2, 3)]
    codons = [c for c in codons if len(c) == 3 and "N" not in c.upper()]
    if not codons:
        return {"expected_stops": None, "p_stop_per_sub": None, "n_codons": 0}
    n_paths = n_to_stop = 0
    for c in codons:
        c = c.upper()
        for pos in range(3):
            for base in "ACGT":
                if base == c[pos]:
                    continue
                n_paths += 1
                if c[:pos] + base + c[pos + 1:] in STOPS:
                    n_to_stop += 1
    p = n_to_stop / max(1, n_paths)
    subs = max(0.0, 1.0 - identity) * len(codons)
    return {"expected_stops": round(p * subs, 1),
            "p_stop_per_sub": round(p, 5), "n_codons": len(codons),
            "implied_substitutions": round(subs, 1)}


def analyse(case: dict, genome_dir: Path, sweep_dir: Path,
            threads: int = 8) -> dict:
    """Reconstruct and grade the reading frame at one case locus."""
    prot, src = _model_protein(sweep_dir, case.get("mp_id", ""), case["cell"],
                               case["contig"], int(case["start"]),
                               int(case["end"]), case["strand"])
    out: dict = {"case_id": case.get("case_id", ""),
                 "accession": case["accession"], "cell": case["cell"],
                 "mp_id": case.get("mp_id", ""), "protein_source": src,
                 "model_protein_aa": len(prot)}
    if not prot:
        out.update(route="no_model_protein", verdict="not_tested")
        return out
    # The locus handed to `find_model` is the miniprot mRNA's own span, not
    # the summary's: S5 pads and clusters, so the two differ by a few bases and
    # the coordinate pre-filter then returns nothing.
    lo, hi = _model_span(sweep_dir, case.get("mp_id", ""))
    cds, note = cds_from_miniprot_model(
        case["accession"], prot, threads=threads,
        locus=(case["contig"], lo or int(case["start"]),
               hi or int(case["end"]), case["strand"]))
    out["route"] = note
    if cds is None:
        out.update(verdict="not_tested")
        return out
    aa = translate(cds)
    internal = aa[:-1] if aa.endswith("*") else aa
    n_stop = internal.count("*")
    n_x = internal.count("X")
    out.update({
        "cds_nt": len(cds), "cds_codons": len(cds) // 3,
        "translated_aa": len(aa), "internal_stops": n_stop,
        "masked_codons": n_x,
        "terminal_stop": int(aa.endswith("*")),
        "miniprot_frameshifts": int(case.get("frameshifts") or 0),
        "miniprot_stop_codons": int(case.get("stop_codons") or 0),
        "bait_identity": float(case.get("identity") or 0),
    })
    out.update(expected_stops(cds, aa, float(case.get("identity") or 0)))
    out["verdict"] = ("open_reading_frame" if n_stop <= MAX_STOPS_FOR_ORF
                      else "disrupted")
    return out


ORF_COLS = ["case_id", "accession", "cell", "mp_id", "protein_source",
            "route", "model_protein_aa",
            "cds_nt", "cds_codons", "translated_aa", "internal_stops",
            "masked_codons", "terminal_stop", "miniprot_frameshifts",
            "miniprot_stop_codons", "bait_identity", "n_codons",
            "implied_substitutions", "p_stop_per_sub", "expected_stops",
            "verdict"]
