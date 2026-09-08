"""Negative controls for the S10 rules, run on every build.

`s5_bait_screen.self_test()`'s pattern: constructed cases, each of which must
be rejected by the rule responsible. These are checks on **refusal**. Every
step in S10 produces a plausible-looking number when it is wrong — a
coordinate parser off by a strand gives every intron a non-canonical splice
site, a tiling step with an incomplete subject set gives a confident naming
disagreement, an exon classifier scoring against gene spans instead of CDS
blocks credits the annotation with coding sequence it never called — and none
of those failures is visible in the output.

Run with `python scripts/s10_test_evidence.py`; non-zero exit on any failure.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s10_case_spec as spec                                   # noqa: E402
import s10_evidence as ev                                      # noqa: E402
import s10_gff as gff                                          # noqa: E402
import s10_probes as pr                                        # noqa: E402
import s10_tile as tile                                        # noqa: E402

FAILURES: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    print(f"  {'ok  ' if ok else 'FAIL'}  {name}" + (f" — {detail}" if detail
                                                     else ""))
    if not ok:
        FAILURES.append(name)


class FakeRegion:
    """A genome region with known sequence, for the splice-site tests."""

    def __init__(self, seq: str, start: int = 1):
        self.seq, self.start = seq.upper(), start
        self.end = start + len(seq) - 1

    def at(self, s: int, e: int) -> str:
        if s < self.start or e > self.end or e < s:
            return ""
        return self.seq[s - self.start:e - self.start + 1]


# --------------------------------------------------------------------------
def t1_interval_algebra() -> None:
    """Merge, subtract and disjoint-count on cases that break naive code."""
    check("T1a touching blocks merge",
          gff.merge([(1, 10), (11, 20)]) == [(1, 20)])
    check("T1b nested blocks merge",
          gff.merge([(1, 100), (20, 30)]) == [(1, 100)])
    check("T1c subtract splits a block",
          gff.subtract([(1, 100)], [(40, 60)]) == [(1, 39), (61, 100)])
    check("T1d subtract to nothing",
          gff.subtract([(10, 20)], [(1, 100)]) == [])
    check("T1e disjoint count is of merged pieces, not of inputs",
          gff.n_disjoint([(1, 10), (11, 20), (100, 110)]) == 2)
    check("T1f overlap of disjoint sets is zero",
          gff.overlap_bp([(1, 10)], [(20, 30)]) == 0)


def t2_splice_strand() -> None:
    """A canonical minus-strand intron must read GT..AG, not CT..AC.

    The single most likely silent error in this module. Genomic GT..AG on the
    plus strand appears as CT..AC when the gene is on the minus strand, so
    reading the dinucleotides without the strand correction makes every
    minus-strand gene look non-canonical — a systematic signal that would
    discredit exactly the alignments this task is defending.
    """
    plus = FakeRegion("AAAA" + "GT" + "N" * 20 + "AG" + "TTTT", start=1)
    donor, acc = ev.splice_pair(plus, (5, 28), "+")
    check("T2a plus-strand canonical intron reads GT..AG",
          (donor, acc) == ("GT", "AG"), f"{donor}..{acc}")
    check("T2b and is classed canonical",
          ev.splice_class(donor, acc) == "canonical")
    # the same intron on a gene transcribed the other way: genomic CT..AC
    minus = FakeRegion("AAAA" + "CT" + "N" * 20 + "AC" + "TTTT", start=1)
    donor, acc = ev.splice_pair(minus, (5, 28), "-")
    check("T2c minus-strand canonical intron also reads GT..AG",
          (donor, acc) == ("GT", "AG"), f"{donor}..{acc}")
    donor, acc = ev.splice_pair(minus, (5, 28), "+")
    check("T2d reading it on the wrong strand is *not* canonical",
          ev.splice_class(donor, acc) == "non_canonical", f"{donor}..{acc}")


def t3_exon_class_uses_cds_not_span() -> None:
    """An exon inside an annotated gene's intron is not annotated coding.

    ITPR introns run to 152 kb and carry passenger genes; scoring an exon
    against gene *spans* would report a gene as annotated because something
    else overlaps it.
    """
    model = {"strand": "+", "cds": [
        {"start": 1000, "end": 1100, "q_start": 1, "q_end": 34,
         "identity": 1.0, "phase": 0}]}
    covering = [{"gene_id": "g1", "name": "PASSENGER", "biotype":
                 "protein_coding", "pseudo": False, "start": 500, "end": 2000,
                 "strand": "+", "cds": [(500, 600), (1800, 2000)],
                 "exons": [], "mrnas": {}, "description": ""}]
    rows = ev.classify_exons(model, covering)
    check("T3a exon in another gene's intron is not 'in_annotated_cds'",
          rows[0]["class"] == "in_annotated_noncoding", rows[0]["class"])
    covering[0]["cds"] = [(990, 1110)]
    rows = ev.classify_exons(model, covering)
    check("T3b exon inside an annotated CDS block is",
          rows[0]["class"] == "in_annotated_cds", rows[0]["class"])
    rows = ev.classify_exons(model, [])
    check("T3c exon with no gene at all is intergenic",
          rows[0]["class"] == "intergenic", rows[0]["class"])


def t4_intron_between_models() -> None:
    """The junction the annotation split a gene at must be found, and only it."""
    model = {"strand": "+", "cds": [
        {"start": 100, "end": 200, "q_start": 1, "q_end": 34, "identity": 1.0,
         "phase": 0},
        {"start": 400, "end": 500, "q_start": 35, "q_end": 68, "identity": 1.0,
         "phase": 0},
        {"start": 800, "end": 900, "q_start": 69, "q_end": 102,
         "identity": 1.0, "phase": 0}]}
    genes = [
        {"gene_id": "A", "name": "A", "biotype": "protein_coding",
         "pseudo": False, "start": 90, "end": 510, "strand": "+",
         "cds": [(100, 200), (400, 500)], "exons": [], "mrnas": {},
         "description": ""},
        {"gene_id": "B", "name": "B", "biotype": "protein_coding",
         "pseudo": False, "start": 790, "end": 910, "strand": "+",
         "cds": [(800, 900)], "exons": [], "mrnas": {}, "description": ""}]
    region = FakeRegion("N" * 1200, start=1)
    rows = ev.classify_introns(model, genes, region)
    classes = [r["class"] for r in rows]
    check("T4a intron inside one model is 'within_one_model'",
          classes[0] == "within_one_model", classes[0])
    check("T4b intron between two models is 'between_models'",
          classes[1] == "between_models", classes[1])
    rows = ev.classify_introns(model, [genes[0]], region)
    check("T4c model to nothing is 'model_to_gap'",
          [r["class"] for r in rows][1] == "model_to_gap")
    rows = ev.classify_introns(model, [], region)
    check("T4d no annotation at all is 'unannotated'",
          {r["class"] for r in rows} == {"unannotated"})


def t5_probe_spans_the_junction() -> None:
    """Only an HSP crossing the junction *and* covering the probe counts.

    T5e is the case the genomic control caught on its first run: blastn runs a
    high-scoring alignment a dozen bases past the true junction into the
    intron, which clears an 8 nt anchor. Six of 55 junctions "spanned" against
    genomic DNA, where nothing can.
    """
    probe = {"junction_index": 1, "junction_offset_in_probe": 90,
             "probe_len": 180}
    def row(qs, qe):
        return ("j1\tsub\t99.0\t%d\t%d\t%d\t1\t100\t1e-40\t150\ttitle"
                % (qe - qs + 1, qs, qe))
    check("T5a an HSP inside one exon does not count",
          pr.count_spanning(row(1, 88), probe)["n_spanning"] == 0)
    check("T5b an HSP crossing by 1 nt does not count",
          pr.count_spanning(row(60, 91), probe)["n_spanning"] == 0)
    # The anchor is necessary and not sufficient: a 16 nt HSP straddling the
    # junction clears it and is still not evidence of splicing.
    anchor_only = pr.count_spanning(
        row(90 - pr.MIN_ANCHOR + 1, 90 + pr.MIN_ANCHOR), probe)
    check("T5c an HSP crossing by exactly the anchor clears the anchor rule",
          anchor_only["n_spanning_anchor_only"] == 1)
    check("T5c2 ...but does not count as spanning on its own",
          anchor_only["n_spanning"] == 0, str(anchor_only))
    check("T5d hits are counted even when none spans",
          pr.count_spanning(row(1, 88), probe)["n_hits"] == 1)
    # The genomic artefact the control caught: an HSP over one exon's flank
    # that has run a dozen bases past the junction into the intron. It clears
    # the anchor and must still be rejected.
    genomic = pr.count_spanning(row(80, 180), probe)
    check("T5e an HSP covering half the probe does not span",
          genomic["n_spanning"] == 0, str(genomic))
    check("T5f but it is still visible as anchor-only",
          genomic["n_spanning_anchor_only"] == 1)
    check("T5g a hit over the whole probe does span",
          pr.count_spanning(row(1, 180), probe)["n_spanning"] == 1)


def t6_spliced_cds_respects_strand() -> None:
    """A minus-strand model splices in target order, not coordinate order.

    Splicing by coordinate reverses the transcript, and the probe for junction
    1 would then be built at the far end of the gene. Nothing downstream would
    notice.
    """
    region = FakeRegion("A" * 100 + "CCC" + "T" * 100 + "GGG" + "A" * 100,
                        start=1)
    model = {"strand": "-", "cds": [
        {"start": 204, "end": 206, "q_start": 1, "q_end": 1, "identity": 1.0,
         "phase": 0},                     # GGG, first in target order
        {"start": 101, "end": 103, "q_start": 2, "q_end": 2, "identity": 1.0,
         "phase": 0}]}                    # CCC, second
    cds, offsets = pr.spliced_cds(region, model)
    check("T6a minus-strand splice starts at the first target block",
          cds.startswith("CCC"), cds)
    check("T6b junction offset is after the first block",
          offsets == [3], str(offsets))


def t7_mode_and_eligibility() -> None:
    """The selection rules must classify and exclude the right things."""
    check("T7a a locus a single model covers is not a failure",
          spec.mode_of(1, 0.20) == "")
    check("T7b no model at all is an omission",
          spec.mode_of(0, 1.0) == "omission")
    check("T7c one partial model is a truncation, and is still a failure",
          spec.mode_of(1, 0.88) == "truncation")
    check("T7d two partial models are fragmentation",
          spec.mode_of(3, 0.58) == "fragmentation")
    base = {"annotated": True, "has_annotation_index": True, "coverage": 1.0,
            "contig_edge": False, "n_gap": False, "contig_spans_gene": True,
            "span": 50_000, "contig_n50": 5_000_000, "n_genes_longer": 400}
    check("T7e a clean locus is eligible", spec.eligible(base)[0])
    check("T7f an unannotated assembly fails E1",
          spec.eligible({**base, "annotated": False})[1] == "E1")
    check("T7g a partial recovery fails E2",
          spec.eligible({**base, "coverage": 0.4})[1] == "E2")
    check("T7h a fragmented assembly fails E3",
          spec.eligible({**base, "contig_spans_gene": False})[1] == "E3")
    check("T7i a marginal contig fails E4",
          spec.eligible({**base, "contig_n50": 100_000})[1] == "E4")
    check("T7j an annotation with a length ceiling fails E5",
          spec.eligible({**base, "n_genes_longer": 1})[1] == "E5")


def t8_selection_is_deterministic_and_one_per_genome() -> None:
    """Two cases, two genomes, one per mode, stable across a reshuffle."""
    recs = [
        {"accession": "G1", "cell": "ITPR1", "eligible": True,
         "mode": "omission", "loss": 1.0, "control_strength": 2,
         "identity": 0.9},
        {"accession": "G1", "cell": "ITPR2", "eligible": True,
         "mode": "fragmentation", "loss": 0.9, "control_strength": 2,
         "identity": 0.9},
        {"accession": "G2", "cell": "ITPR3", "eligible": True,
         "mode": "fragmentation", "loss": 0.6, "control_strength": 1,
         "identity": 0.8},
        {"accession": "G3", "cell": "ITPR1", "eligible": True,
         "mode": "omission", "loss": 1.0, "control_strength": 0,
         "identity": 0.99},
    ]
    got = spec.select_cases(recs)
    check("T8a one case per mode", [c["selected_as"] for c in got]
          == list(spec.MODES))
    check("T8b the two cases are in different genomes",
          len({c["accession"] for c in got}) == 2)
    check("T8c control strength outranks identity at equal loss",
          got[0]["accession"] == "G1",
          f"picked {got[0]['accession']}")
    check("T8d selection is order-invariant",
          [(c["accession"], c["cell"]) for c in spec.select_cases(recs[::-1])]
          == [(c["accession"], c["cell"]) for c in got])
    check("T8e an ineligible locus is never selected",
          not spec.select_cases([{**r, "eligible": False} for r in recs]))


def t9_tiling_needs_the_right_locus() -> None:
    """A model whose best hit is not at its own coordinates is not evidence.

    The failure this guards is real and was measured during the build: with an
    incomplete subject set, *Nibea albiflora*'s ITPR1 pseudogene matched an
    ITPR1 locus 19 Mb away and the row read as a confident naming call.
    """
    row = {"start": 1000, "end": 2000,
           "best_locus": "GCA_x|ITPR1|CTG1:900-2100-|sweep", "name": "ITPR2"}
    check("T9a a model over its own locus is placed",
          tile.at_own_locus(row) == 1)
    far = {**row, "best_locus": "GCA_x|ITPR1|CTG1:900000-910000-|sweep"}
    check("T9b a model matching a distant locus is not placed",
          tile.at_own_locus(far) == 0)
    check("T9c an unplaced model gets no naming verdict",
          tile.naming_verdict(far, tile.NAME_HINTS) == "undetermined_locus")
    check("T9d a placed model with the wrong name is a mismatch",
          tile.naming_verdict(row, tile.NAME_HINTS) == "name_mismatch")
    check("T9e a placed model with the right name matches",
          tile.naming_verdict({**row, "name": "ITPR1"}, tile.NAME_HINTS)
          == "name_matches_sequence")
    check("T9f a locus tag claims no paralog",
          tile.naming_verdict({**row, "name": "KUDE01_018466"},
                              tile.NAME_HINTS) == "unnamed")


def t10_concordance_can_fail() -> None:
    """The boundary control must be able to report a boundary nobody shares."""
    def model(qs):
        return {"cds": [{"q_start": 1, "q_end": q, "start": q, "end": q,
                         "identity": 1.0, "phase": 0} for q in qs]}
    case = model([100, 200, 300, 400])
    refs = [{"accession": "A", "organism": "a", "model": model([100, 200, 999,
                                                                400])},
            {"accession": "B", "organism": "b", "model": model([100, 200, 998,
                                                                400])}]
    rows, summ = pr.concordance(case, refs)
    shared = {r["q_residue"]: r["n_refs_with_boundary"] for r in rows}
    check("T10a a shared boundary is counted", shared[100] == 2)
    check("T10b a boundary no reference has scores zero", shared[300] == 0)
    check("T10c the summary reports the shortfall",
          summ["n_shared_with_majority"] == 2 and summ["n_boundaries"] == 3,
          str(summ))
    within = model([102, 200, 300, 400])
    rows, _ = pr.concordance(within, refs[:1])
    check("T10d a boundary within tolerance still matches",
          rows[0]["n_refs_with_boundary"] == 1)


def t11_pseudogene_cds_is_kept() -> None:
    """A GFF3 pseudogene carrying CDS features must not be dropped.

    Dropping it makes S10 report 'no coding model here' about a locus that has
    one the submitter demoted — which is the opposite of Case A's finding.
    """
    text = "\n".join([
        "##gff-version 3",
        "#!genome-build TESTv1",
        "#!genome-build-accession NCBI_Assembly:GCA_000000000.1",
        "\t".join(["ctg", "Genbank", "pseudogene", "100", "500", ".", "+", ".",
                   "ID=gene-P1;Name=ITPR2;gene_biotype=pseudogene;pseudo=true"]),
        "\t".join(["ctg", "Genbank", "mRNA", "100", "500", ".", "+", ".",
                   "ID=rna-P1;Parent=gene-P1;pseudo=true"]),
        "\t".join(["ctg", "Genbank", "CDS", "100", "200", ".", "+", "0",
                   "ID=cds-P1;Parent=rna-P1"]),
        "\t".join(["ctg", "Genbank", "CDS", "400", "500", ".", "+", "2",
                   "ID=cds-P1;Parent=rna-P1"]),
    ]) + "\n"
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.gff"
        p.write_text(text)
        win = gff.read_annotation_window(p, "ctg", 1, 1000)
        meta = gff.annotation_meta(p)
    genes = win["genes"]
    check("T11a the pseudogene is present", len(genes) == 1)
    check("T11b it keeps its CDS blocks",
          genes and genes[0]["cds"] == [(100, 200), (400, 500)])
    check("T11c it is flagged pseudo", genes and genes[0]["pseudo"])
    check("T11d and it keeps its family name",
          genes and genes[0]["name"] == "ITPR2")
    tx = list(genes[0]["mrnas"].values())[0] if genes else {}
    check("T11e per-block phase is recorded",
          tx.get("cds_phase", {}).get((400, 500)) == 2)
    check("T11f the build stamp parses past the '#!' pragma",
          meta.get("build_accession") == "NCBI_Assembly:GCA_000000000.1",
          str(meta))


def main() -> int:
    print("[s10] negative controls")
    for fn in (t1_interval_algebra, t2_splice_strand,
               t3_exon_class_uses_cds_not_span, t4_intron_between_models,
               t5_probe_spans_the_junction, t6_spliced_cds_respects_strand,
               t7_mode_and_eligibility,
               t8_selection_is_deterministic_and_one_per_genome,
               t9_tiling_needs_the_right_locus, t10_concordance_can_fail,
               t11_pseudogene_cds_is_kept):
        print(f" {fn.__name__}")
        fn()
    if FAILURES:
        print(f"[s10] {len(FAILURES)} FAILED: {', '.join(FAILURES)}")
        return 1
    print("[s10] all negative controls pass")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
