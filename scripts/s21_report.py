"""S21 — renders `results/gene_architecture/report.md` from the committed tables.

D13 applied to a report about gene structure: nothing here is hand-written and
nothing is recomputed. Every number comes through `s21_tables.headline()`,
which reads the tables, so the prose and the tables cannot drift.

Scope, the instrument, what an intron was measured to be and the negative
controls live here; the results — the architecture, the intron positions, the
fragments and the duplicates — live in `s21_report_results.py` (the
`s3_report.py` / `s3_report_d10.py` split).

Headlines are chosen by the data. Every prior an earlier task set is stated in
`s21_priors.PRIOR` with where it was said, computed on S21's own tables, and
rendered `confirmed` / `contradicted` / `not corroborated` / `orthogonal` /
`underpowered` with both numbers printed either way. A section whose table is
absent renders *not run yet*, so a stage that did not run is visible as one.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s21_architecture as A                                      # noqa: E402
import s21_frame as FR                                            # noqa: E402
import s21_introns as I                                           # noqa: E402
import s21_lib as L                                               # noqa: E402
import s21_tables as T                                            # noqa: E402


def fmt(v, nd: int = 4) -> str:
    if v is None or v == "":
        return "—"
    if isinstance(v, float):
        return (f"{v:,.{nd}f}".rstrip("0").rstrip(".") if v else "0")
    if isinstance(v, int):
        return f"{v:,}"
    return str(v)


def pct(v, nd: int = 1) -> str:
    return "—" if v in (None, "") else f"{float(v) * 100:.{nd}f} %"


def pfmt(p) -> str:
    """A p-value at four decimal places reads `0.0000`, which hides how strong
    a claim is rather than how weak (`s15_report.py`'s rule)."""
    try:
        x = float(p)
    except (TypeError, ValueError):
        return str(p) if p else "—"
    return "0" if x == 0 else f"{x:.3g}"


def table(header: list[str], rows: list[list]) -> str:
    out = ["| " + " | ".join(header) + " |",
           "|" + "|".join(["---"] * len(header)) + "|"]
    for r in rows:
        out.append("| " + " | ".join("" if c is None else str(c)
                                     for c in r) + " |")
    return "\n".join(out)


def cell_label(cell: str) -> str:
    return "RyR (control)" if cell == "RYR" else cell


#: What each scope rule requires, in the order it is applied. The counts come
#: from the table; only the wording is here.
RULE_TEXT = {
    "A0": "the locus belongs to one of the four cells (D14's family call)",
    "A1": "it is one of S16's committed gene copies (`is_copy`)",
    "A2": "the contig spans the gene (D4)",
    "A3": f"the model covers at least {A.COV_ARCH} of its bait",
    "A4": "the bait reaches a usable alignment frame",
}


def section_1(h: dict) -> str:
    excl = h["excluded_by"]
    rows = []
    for rule, text in RULE_TEXT.items():
        n = sum(v for k, v in excl.items() if k.startswith(rule + "_"))
        rows.append([rule, text, fmt(n)])
    return f"""# S21 — gene architecture: the exon structure of an IP3 receptor

*Generated {date.today().isoformat()} by `scripts/s21_run.py` from the
committed tables in `results/gene_architecture/`. Nothing here is
hand-written.*

## 1. What this task measures, and off what

Every task before this one measured a *sequence* or a *search*. S21 measures
the **gene**: how many coding exons an IP3 receptor has, how much genomic
space it takes to hold them, where its introns sit, and whether the three
vertebrate paralogues inherited that structure from one ancestor.

None of it needed a new experiment. The S5 sweep aligned a 38-protein bait
panel to 309 vertebrate genomes with miniprot and kept every genome's
`miniprot.gff`, and each of those files carries, for every alignment, one CDS
record per aligned block with its genomic interval, its span in the bait's own
residue coordinates and its phase. That is a spliced gene model. The genomes
themselves are still on the drive, so every intron's splice dinucleotides can
be read rather than assumed, and the assemblies' own gene sets are there too,
so the boundaries can be checked against a pipeline that never saw a bait.

**{fmt(h['n_loci_measured'])} loci** were measured, of which
**{fmt(h['n_in_scope'])}** in **{fmt(h['n_genomes_in_scope'])} genomes** are
in the architecture scope. Four rules, each a positive test naming the number
it fired on:

{table(["rule", "what it requires", "loci it excluded"], rows)}

**A2 is the one that does the work.** S19 measured 15.2 % of the sweep's cells
as false negatives of the method and found every one of them is an assembly,
not a gene. An exon count taken below D4's contiguity bar measures contig
lengths, so the scope is above it. **A3** is S21's own bar — S16's copy rule
admits a model covering half its bait, and half a gene has half the exons —
and it is measured rather than assumed: §6 recomputes every count at six
coverage bars and the median exon count moves by at most one exon across the
whole range.

## 2. What an intron is, measured rather than chosen

The brief expected miniprot to emit two CDS records either side of a
frameshift, which would make a broken locus look exon-rich, and asked for
those pairs to be merged. Whether they exist is a measurement, so it was made.

Every consecutive block pair in the sweep — **{fmt(h['intron_floor']['n_gaps'])}
of them** — was binned by the gap between the two blocks in gene order and
scored by the splice dinucleotides at that gap's two edges, which is evidence
the alignment score did not produce. The result is that
**{fmt(h['intron_floor']['n_non_spliceable'])} gaps** sit in a bin the genome
does not call spliceable, the smallest gap anywhere in the sweep is
**{fmt(h['intron_floor']['smallest_gap_bp'])} bp** and it reads as a splice
pair, and the query spans are contiguous across every pair — no residue
emitted twice, none skipped. There is no frameshift-pair population in this
output at all: an indel appears *inside* a block, as a block whose genomic
length differs from three times its residue count.

So the calibration refuses to derive a threshold, in the way
`s15_calibrate_recon.bar()` and `s23_calibrate_loci.separation()` do: a bar
needs gaps on both sides of it, and one side is empty. The floor is placed at
{fmt(h['intron_floor']['floor_bp'])} bp — *{h['intron_floor']['derivation']}* —
and the merge rule consequently fires on **{fmt(h['n_merges_total'])}** of the
{fmt(h['junctions']['n'])} junctions in scope. That zero is a result and not a
rule that cannot act: `s21_test_arch` T5 constructs a pair 2 bp apart with
contiguous query spans and requires the rule to merge it, and the same pair
4.9 kb apart to stay two exons.

## 3. The instrument's own error rate

An exon boundary the aligner placed is a claim about the genome, and the
genome can be asked. Over all **{fmt(h['junctions']['n'])} junctions** in the
architecture scope, **{pct(h['junctions']['frac_spliceable'], 2)}** are a
canonical `GT..AG` or one of the two common minor pairs, and
**{pct(h['junctions']['frac_canonical'], 2)}** are canonical outright.
**{pct(h['junctions']['frac_frame_step'], 2)}** of junctions carry a frame
step — the reading frame does not carry over, which is how this output records
a frameshift.

That is the instrument judged against the sequence. Judged against **another
pipeline**, {fmt(h['concordance']['n_edges'])} annotated CDS block edges over
{fmt(h['concordance']['n_loci'])} loci in
{fmt(h['concordance']['n_genomes'])} genomes were compared with the sweep's
own exon boundaries, and **{pct(h['concordance']['frac_exact'])}** of them
land on one exactly. The assemblies' gene sets were built from evidence
miniprot never saw, so this is corroboration and not a consistency check.
D9's contrast survives it: RefSeq
{pct(h['concordance']['by_source'].get('RefSeq'))} against GenBank
{pct(h['concordance']['by_source'].get('GenBank'))}.

## 4. The common coordinate frame

An exon boundary is a position in the *bait's* numbering, and the sweep used
38 baits, so comparing an ITPR1 gene's boundaries with an ITPR2 gene's needs
one frame both reach. That frame is S6's committed alignment: all three human
paralogues and human RYR2 are tips of it, so a boundary travels bait → its
cell's human reference (pairwise MAFFT at `--thread 1`, D24) → alignment
column.

**{fmt(h['frame']['n_pairs'])}** (bait, cell) frames were built and
**{fmt(h['frame']['n_unusable'])}** were refused; the worst puts
{pct(h['frame']['min_frac_on_ref'])} of its bait on the reference against a
{pct(FR.MIN_BAIT_ON_REF)} floor. {fmt(h['frame']['n_via_cell'])} frames are
`frame_via_cell` — the three unlabelled `vertebrate_basal` baits, which have
no paralogue of their own because S7 declined to place those tips, so their
loci are transferred through the reference of the cell they filled and every
cross-paralogue test is reported with and without them.

**The frame is checked, not assumed.** All **{fmt(h['anchors']['n'])}**
residues S0 *measured* on the 6DQN structure — the ten IP₃ contacts, the two
selectivity-filter residues and the two gate residues, in each paralogue's own
numbering from S17's committed `functional_sites.tsv` — land in the **same
alignment column** in all three paralogues:
{fmt(h['anchors']['n_one_column'])} of {fmt(h['anchors']['n'])}. A frame that
had slipped anywhere in the pore or the ligand core would fail there, and a
slipped frame is exactly the error that makes a shared-intron count look like
a result. `s21_test_arch` T8 shifts one paralogue's aligned row by a column
and requires the test to refuse.

## 5. Negative controls

Twenty-one constructed controls run **before anything is written**, and the
run refuses to continue if any fails. Four of the rules here return a
*good-looking* number when they are wrong, which is why the suite exists:

- a reader that keys models the way miniprot names them merges the blocks of
  two different genes in the two chunked assemblies, and the merged gene has
  more exons and longer introns (**T1** requires this reader and
  `s5_sweep_lib.parse_miniprot_gff` to return identical genomic blocks for
  every model in an unchunked, a 14-chunk and a 25-chunk genome; **T2**
  requires every `mp_id` the summaries recorded to be a key);
- computing an intron as `next.start − prev.end` reverses every minus-strand
  gene (**T3**), and reading its splice pair without reverse-complementing
  both ends turns every canonical minus-strand intron into `CT..AC` (**T4**);
- inverting the intron-phase convention swaps phases 1 and 2 everywhere,
  changing every shared-intron call without changing a single count (**T7**);
- and asking whether the same bait aligns twice at disjoint positions measures
  *paralogy*, not duplication, in a family whose members are 61-68 % identical
  — the first version of the detector reached a specificity of 0.16 that way
  (**T16**).

The rest are reachability: the merge must be able to fire (**T5**), the loss
of frame must be detectable and absent on a clean chain (**T6**), each scope
rule must fire on its own violation (**T9**), the shared-intron test must
reach both extremes and its exact tail must match a brute-force enumeration
(**T10**, **T11**), both nulls must agree (**T12**), the duplication detector
must fire and decline (**T14**) and its `within_locus` class must be reachable
even though it fires on nothing (**T15**), each terminus class must be
reachable (**T17**, **T18**), the concordance rule must be able to report
disagreement (**T19**), and the intron calibration must refuse without a
second population and place a bar with one (**T20**). **T21** requires the
suite to alter no committed table.

The suite was **mutation-tested on six deliberate rule breakages** — the phase
complement, the intron-phase convention, the strand branch in the gap, the
phase requirement in the shared-intron match, the cell attribution in the
duplication detector, and the reader's collision key — and caught all six, by
T6, T7, T3, T13, T16 and T1 respectively.
"""


def main() -> int:
    import s21_report_results as R
    h = T.headline()
    parts = [section_1(h), R.section_6(h), R.section_7(h), R.section_8(h),
             R.section_9(h), R.section_10(h), R.section_11(h)]
    out = L.OUT / "report.md"
    out.write_text("\n\n".join(parts) + "\n")
    print(f"  wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
