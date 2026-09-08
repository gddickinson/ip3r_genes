"""S10 — which annotation failures get validated, and why those two.

The brief asks for *the two worst annotation failures S5b surfaced*. "Worst"
has to be a rule, not a preference, because the whole point of the task is to
prove that a handful of hand-picked loci are not what the audit rests on. So
this module writes the rule out and `s10_select.py` applies it to every locus
the sweep found.

**The measurement.** For each locus the sweep recovered, `annotation_loss` is
``1 - max(frac_cds)`` over the annotated genes overlapping it: the fraction of
the gene's real coding footprint that *no single annotated gene model*
delivers. It is read off S5's own per-locus record, so S10 does not re-derive
what S5 already measured (D13). A locus at loss 0 is a gene the annotation got
right; a locus at loss 1 has no coding model on it at all.

**The five eligibility rules**, each a positive test, in the order they are
applied. E1–E4 establish that a failure *is* one; E5 establishes that it is a
failure of the annotation rather than of something else.

* **E1 — there is an annotation to audit.** The assembly carries a gene set
  and the sweep built an index from it. Without this every locus is trivially
  "unannotated" and the ranking measures nothing.
* **E2 — the sweep recovered the gene.** Bait coverage ≥ `COV_FOUND`, no
  contig edge, no N-gap inside the locus. A partial alignment cannot show that
  an annotation is missing anything.
* **E3/E4 — the assembly can carry the gene.** D4's contiguity bar, in its
  per-locus form: the contig holds the whole locus, and the assembly's contig
  N50 is at least `HEADROOM_X` times the locus span. A fragmented gene model
  on a fragmented assembly is an assembly result.
* **E5 — the annotation can build a gene this long.** At least `MIN_LONGER`
  annotated protein-coding genes elsewhere in the *same* genome are longer
  than this locus. This is the control that does the work, and it is not
  optional: an annotation whose longest gene anywhere is 42 kb has not failed
  *at the IP3 receptor*, it has failed everywhere, and validating one of its
  loci would report a genome-wide property as a locus-specific bug. Measured
  over the 309-genome sweep it removes six loci in two genomes and nothing
  else — see `annotation_depth.tsv`.

**Two modes, one case each.** S5b surfaced two distinct kinds of failure and
they need different evidence:

* `omission` — no annotated coding model overlaps the gene's CDS at all.
* `fragmentation` — two or more annotated models overlap it and none of them
  covers half.

Taking the top two of one ranking would have given two omissions and left the
fragmentation claim unvalidated, so the rule takes **the worst of each mode**,
at most one case per genome. The full ranking is committed either way
(`case_ranking.tsv`), so a reader can see what a single-ranking selection
would have chosen.

A third quantity, `control_strength`, breaks ties: how many *other* loci of
this family in the same genome the annotation gets right. It is the same idea
as `s23_controls.control_strength()` — grade the control, do not merely check
that one exists. A missing gene in a genome where the annotation nails the
other two paralogs is a far sharper result than the same gene missing from a
genome that got nothing right, and at equal loss that is what should decide.
"""

from __future__ import annotations

# --- E2 -------------------------------------------------------------------
#: Bait coverage a locus needs before its annotation can be said to be missing
#: anything. Imported rather than retyped so S10's notion of "the sweep found
#: this gene" is S5's (see `s5_sweep_lib.COV_FOUND`).
try:                                                  # pragma: no cover
    from s5_sweep_lib import COV_FOUND
except Exception:                                     # pragma: no cover
    COV_FOUND = 0.90

# --- E4 -------------------------------------------------------------------
#: Contig N50 must exceed the locus span by this factor. E3 (the contig holds
#: the locus) is the bar D4 sets; this is the stronger version, and it is
#: deliberately strong: at 10x an assembly is nowhere near the regime where
#: contiguity could explain a broken gene model.
HEADROOM_X = 10.0

# --- E5 -------------------------------------------------------------------
#: How many annotated protein-coding genes elsewhere in the genome must be
#: longer than the locus before the annotation is credited with being able to
#: build a gene this long. Ten rather than one, so a single chimeric model
#: cannot certify an annotation that otherwise tops out well below the locus.
MIN_LONGER = 10

# --- modes ----------------------------------------------------------------
#: A model has to reach this fraction of the locus's CDS footprint to count as
#: one of the fragments, rather than a passenger gene sitting in an intron.
#: ITPR introns run to 152 kb (`intron_calibration.tsv`), so there is room for
#: several, and an intronic passenger scores ~0.
FRAG_MIN_FRAC_CDS = 0.01

#: `annot_gene` in S5 is set when one gene covers this much of the CDS, so a
#: locus at or below it is one no single model delivers. Imported for the same
#: reason as COV_FOUND.
try:                                                  # pragma: no cover
    from s5_classify import ANNOT_CDS_FRAC
except Exception:                                     # pragma: no cover
    ANNOT_CDS_FRAC = 0.50

#: The two modes a case may be selected from.
MODES = ("omission", "fragmentation")

#: `truncation` is a failure too, and it is *reported* — it is simply not
#: selectable, because the brief's method (tile the fragments, classify the
#: introns between them) needs more than one fragment to work on. Leaving it
#: unlabelled would make it read in the ranking as "the annotation got this
#: right", which is the opposite of what it is.
FAILURE_MODES = MODES + ("truncation",)

MODE_RULE = {
    "omission": (
        "no annotated coding model overlaps the recovered gene's CDS "
        f"footprint (loss = 1.0, 0 models above {FRAG_MIN_FRAC_CDS:g} of it)"),
    "fragmentation": (
        f"two or more annotated models overlap the CDS footprint and none "
        f"covers {ANNOT_CDS_FRAC:.0%} of it"),
    "truncation": (
        "exactly one annotated model overlaps the CDS footprint and it "
        f"covers less than {ANNOT_CDS_FRAC:.0%} of it"),
}

ELIGIBILITY = [
    ("E1", "assembly carries a gene set and the sweep indexed it"),
    ("E2", f"sweep recovered the gene: coverage >= {COV_FOUND}, "
           "no contig edge, no N-gap"),
    ("E3", "the contig holds the whole locus (D4)"),
    ("E4", f"contig N50 >= {HEADROOM_X:g}x the locus span"),
    ("E5", f">= {MIN_LONGER} annotated protein-coding genes elsewhere in the "
           "same genome are longer than the locus"),
]


def mode_of(n_fragments: int, loss: float) -> str:
    """Which failure mode a locus is in, or '' if it is not a failure.

    The three are exclusive by construction, on the number of annotated
    models overlapping the CDS footprint: zero, one, or more than one. Only
    the first and third are selectable (see `MODES`); the middle one is
    reported so that a locus whose single annotated model delivers a third of
    the gene is not filed alongside the loci that are annotated correctly.
    """
    if loss <= 1.0 - ANNOT_CDS_FRAC:
        return ""                       # a model covers half: not a failure
    if n_fragments == 0:
        return "omission"
    if n_fragments == 1:
        return "truncation"
    return "fragmentation"


def eligible(rec: dict) -> tuple[bool, str]:
    """Apply E1–E5 to one locus record. Returns (ok, first failing rule)."""
    if not rec.get("annotated") or not rec.get("has_annotation_index"):
        return False, "E1"
    if (rec.get("coverage") or 0) < COV_FOUND:
        return False, "E2"
    if rec.get("contig_edge") or rec.get("n_gap"):
        return False, "E2"
    if not rec.get("contig_spans_gene"):
        return False, "E3"
    span = max(1, rec.get("span") or 1)
    if (rec.get("contig_n50") or 0) < HEADROOM_X * span:
        return False, "E4"
    if (rec.get("n_genes_longer") or 0) < MIN_LONGER:
        return False, "E5"
    return True, ""


def severity_key(rec: dict) -> tuple:
    """Sort key for the ranking: worst first.

    Loss is the measurement. `control_strength` breaks ties, because at equal
    loss the sharper case is the one whose own genome proves the annotation
    could have done better. Identity breaks what that leaves, since it is the
    evidence that the recovered model really is this gene; and the accession
    makes the order total, so a re-run cannot reorder a tie.
    """
    return (-(rec.get("loss") or 0.0),
            -(rec.get("control_strength") or 0),
            -(rec.get("identity") or 0.0),
            rec.get("accession", ""), rec.get("cell", ""))


def select_cases(records: list[dict]) -> list[dict]:
    """The worst eligible locus of each mode, at most one case per genome.

    Order of the returned cases follows `MODES`, not severity, so the case
    labels (`case_a`, `case_b`) are stable across re-runs even if a future
    sweep changes which mode holds the worst locus overall.
    """
    ok = [r for r in records if r.get("eligible") and r.get("mode")]
    ok.sort(key=severity_key)
    cases, used_genomes = [], set()
    for mode in MODES:
        for r in ok:
            if r["mode"] != mode or r["accession"] in used_genomes:
                continue
            case = dict(r)
            case["case_id"] = f"case_{chr(ord('a') + len(cases))}"
            case["selected_as"] = mode
            case["selection_rank_in_mode"] = 1 + sum(
                1 for x in ok if x["mode"] == mode
                and severity_key(x) < severity_key(r))
            cases.append(case)
            used_genomes.add(r["accession"])
            break
    return cases


def control_strength(cells: dict, this_cell: str) -> int:
    """How many *other* family loci in this genome the annotation gets right.

    The family's own loci are the right control here: same genome, same
    annotation run, same class of gene (very long, many-exon, deeply
    conserved). A housekeeping gene would not test the same thing.
    """
    n = 0
    for name, cell in cells.items():
        if name == this_cell:
            continue
        if cell.get("status") == "found_annotated":
            n += 1
    return n
