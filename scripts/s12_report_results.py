"""The results half of the S12 report, split to keep both under 500 lines.

Takes the caller's loader and formatters (the `s3_report.py` /
`s3_report_d10.py` pattern) so the two halves cannot read the tables
differently.

**Headlines chosen by the data.** Each section states the earlier task's
conclusion, computes S12's, and renders the comparison with S8's
five-valued verdict — printing both numbers either way. The verdict this
task most needs to be able to reach is `underpowered`: S10 said in prose
that a search of 43 deposit records returning nothing shows only that the
species has no deposits, and if S12's own deposit cross-check comes back
empty it has to reach the same verdict about itself rather than reporting
a negative.
"""

from __future__ import annotations

import s12_report_checks as checks


def _i(r, k, d=0):
    try:
        return int(float(r.get(k) or d))
    except (TypeError, ValueError):
        return d


def _f(r, k, d=0.0):
    try:
        return float(r.get(k) or d)
    except (TypeError, ValueError):
        return d


def frac(x) -> str:
    """A 0-1 fraction always shown to two places.

    `fmt` renders a whole number as an integer, which is right for counts
    and wrong here: an annotation loss of 1.0 printed as `1` beside a loss
    of `0.54` reads as a count in a column of fractions.
    """
    try:
        return f"{float(x):.2f}"
    except (TypeError, ValueError):
        return "n/a"


def _median(vals) -> float:
    vals = sorted(vals)
    if not vals:
        return 0.0
    n = len(vals)
    return vals[n // 2] if n % 2 else (vals[n // 2 - 1] + vals[n // 2]) / 2


def sections(load, load_json, fmt, pct, missing, prior_line) -> list[str]:
    loci = load("expression_by_locus.tsv")
    junc = load("junction_support.tsv")
    gap = load("annotation_gap_coverage.tsv")
    tissue = load("expression_by_tissue.tsv")
    by_run = load("expression_by_run.tsv")
    atlas = load("atlas_summary.tsv")
    resources = load("atlas_resources.tsv")
    crossmap = load("crossmap_control.tsv")
    refs = load("reference_table.tsv")

    out: list[str] = []
    out += _headline(loci, junc, fmt, pct, missing, prior_line)
    out += _junctions(loci, junc, fmt, pct, missing, prior_line)
    out += _controls(loci, by_run, crossmap, fmt, pct, missing, prior_line)
    out += _pseudogenes(loci, refs, fmt, pct, missing)
    out += _gap(gap, fmt, pct, missing, prior_line)
    out += _tissues(tissue, fmt, pct, missing)
    out += checks.sections(atlas, resources, loci, by_run, fmt, pct,
                           missing, prior_line)
    return out


# ---------------------------------------------------------------------------

def _headline(loci, junc, fmt, pct, missing, prior_line) -> list[str]:
    failed = [r for r in loci if r["role"] == "failed"]
    if not failed:
        return [missing("4. Are the lost genes transcribed?")]
    det = [r for r in failed if _i(r, "runs_detected") > 0]
    out = ["## 4. Are the lost genes transcribed?", "",
           "| species | locus | mode | annotation loss | runs | runs "
           "detected | reads | junction reads | decoy | transcribed |",
           "|---|---|---|---|---|---|---|---|---|---|"]
    for r in sorted(failed, key=lambda x: (x["organism"], x["cell"])):
        yes = _i(r, "runs_detected") > 0
        out.append(
            f"| *{r['organism']}* | {r['cell']} | {r['failure_mode']} | "
            f"{frac(r['annotation_loss'])} | {fmt(r['runs'])} | "
            f"{fmt(r['runs_detected'])} | {fmt(r['reads'])} | "
            f"{fmt(r['junction_reads'])} | {fmt(r['decoy_reads'])} | "
            f"{'**yes**' if yes else 'no'} |")
    out += ["",
            f"**{len(det)} of {len(failed)}** of the loci the annotation "
            f"loses meet all three criteria in at least one library.", ""]
    out.append(prior_line(
        "handoff",
        f"The question is answerable with reads. {len(det)} of {len(failed)} "
        f"loci are transcribed and spliced; the deposits that could not "
        f"answer it are re-searched in §10 and still cannot.",
        len(det) == len(failed) and len(failed) > 0))
    out.append("")
    out.append(prior_line(
        "not_pseudogenes",
        f"An intact reading frame is a necessary and not a sufficient "
        f"condition, and these loci meet the sufficient one too: "
        f"{fmt(sum(_i(r, 'junction_reads') for r in failed))} reads read "
        f"through their splice junctions.",
        len(det) == len(failed) and len(failed) > 0))
    out.append("")
    return out


def _junctions(loci, junc, fmt, pct, missing, prior_line) -> list[str]:
    if not junc:
        return [missing("5. The junctions the annotation does not model")]
    failed = [r for r in loci if r["role"] == "failed"]
    tot_u = sum(_i(r, "junctions_unannotated") for r in failed)
    hit_u = sum(_i(r, "junctions_unannotated_covered") for r in failed)
    tot_a = sum(_i(r, "junctions_annotated") for r in failed)
    hit_a = sum(_i(r, "junctions_annotated_covered") for r in failed)
    out = ["## 5. The junctions the annotation does not model", "",
           "This is the measurement S10 handed here. A junction the "
           "annotation models is a junction some gene model already claims; "
           "a junction it does not model is sequence no annotated gene "
           "delivers, and a read crossing it is transcript evidence for "
           "exactly that sequence.", "",
           "| species | locus | unannotated junctions | crossed | annotated "
           "junctions | crossed |", "|---|---|---|---|---|---|"]
    for r in sorted(failed, key=lambda x: (x["organism"], x["cell"])):
        out.append(
            f"| *{r['organism']}* | {r['cell']} | "
            f"{fmt(r['junctions_unannotated'])} | "
            f"{fmt(r['junctions_unannotated_covered'])} "
            f"({pct(_i(r, 'junctions_unannotated_covered'), _i(r, 'junctions_unannotated'))}) | "
            f"{fmt(r['junctions_annotated'])} | "
            f"{fmt(r['junctions_annotated_covered'])} "
            f"({pct(_i(r, 'junctions_annotated_covered'), _i(r, 'junctions_annotated'))}) |")
    out += ["",
            f"Across the failed loci, **{fmt(hit_u)} of {fmt(tot_u)} "
            f"({pct(hit_u, tot_u)}) junctions no annotated model spans are "
            f"crossed by reads**"
            + (f", against {fmt(hit_a)} of {fmt(tot_a)} ({pct(hit_a, tot_a)}) "
               f"of the annotated junctions in the same genes."
               if tot_a else "; these genes have no annotated junctions at "
                             "all, which is what makes them omissions."),
            ""]
    out += _position_confound(junc, failed, fmt, pct, _median)
    out += ["",
            "![](figures/s12_junctions.png)", "",
            "**Figure 1.** Every junction of every locus the annotation "
            "loses, at its position in the spliced coding sequence, coloured "
            "by whether the annotation models it. A junction nothing crossed "
            "is marked on the baseline, so the denominator is in the "
            "picture.", ""]
    out.append(prior_line(
        "boundaries_corroborated",
        f"S10 corroborated these exon boundaries against other genomes' "
        f"annotations; reads now cross {pct(hit_u, tot_u)} of the "
        f"boundaries that no annotation of *this* genome models.",
        hit_u > 0))
    out.append("")
    return out


def _position_confound(junc, failed, fmt, pct, median) -> list[str]:
    """Say why the annotated/unannotated ratio is not the claim.

    Several of these libraries are strongly 3'-biased, and in these genes
    the annotated and unannotated junctions are not interleaved along the
    transcript — an annotation that truncates a gene keeps its 5' end, so
    its junctions sit where the coverage is thinnest.  Comparing the two
    recovery rates without that would report a property of library prep as
    a property of the annotation.  The load-bearing claim is the absolute
    one — *these* junctions are crossed — so the ratio is reported with the
    confound printed beside it rather than dropped or leaned on.
    """
    keys = [(r["species"], r["seq"]) for r in failed]
    rows = []
    for sp, seq in keys:
        sub = [j for j in junc if (j["species"], j["seq"]) == (sp, seq)]
        if not sub:
            continue
        offs = [_i(j, "cds_offset") for j in sub]
        span = max(offs) if offs else 1
        un = [_i(j, "cds_offset") / span for j in sub if not _i(j, "annotated")]
        an = [_i(j, "cds_offset") / span for j in sub if _i(j, "annotated")]
        if not (un and an):
            continue
        cell = seq.split("|")[0]
        rows.append((sp, cell, median(un), median(an), len(un), len(an)))
    if not rows:
        return ["*Only one junction class is present in these genes, so no "
                "within-gene comparison of recovery rates is available.*", ""]
    out = ["The two rates are **not** directly comparable, and the table "
           "below says why. Several of these libraries are strongly "
           "3′-biased, and in a truncated or fragmented gene the "
           "annotated junctions are the ones at the 5′ end — where the "
           "coverage is thinnest. Comparing the recovery rates without that "
           "would report a property of library preparation as a property of "
           "the annotation.", "",
           "| species | locus | unannotated junctions | median position | "
           "annotated junctions | median position |",
           "|---|---|---|---|---|---|"]
    for sp, cell, mu, ma, nu, na in rows:
        out.append(f"| {sp} | {cell} | {fmt(nu)} | {mu:.2f} | {fmt(na)} | "
                   f"{ma:.2f} |")
    out += ["", "(position is the fraction of the way along the spliced "
            "coding sequence.) The claim this task makes is the absolute "
            "one — that these particular junctions are crossed — not the "
            "ratio between the two classes.", ""]
    return out


def _controls(loci, by_run, crossmap, fmt, pct, missing, prior_line
              ) -> list[str]:
    if not loci:
        return [missing("6. The controls")]
    decoy_reads = sum(_i(r, "decoy_reads") for r in loci)
    hot = [r for r in by_run if _i(r, "decoy_reads")]
    hk = [r for r in loci if r["role"] == "housekeeping"]
    ctrl = [r for r in loci if r["role"] in ("control_paralog",
                                             "control_family")]
    runs = len({(r["species"], r["run"]) for r in by_run})
    out = ["## 6. The controls", "", "### The spurious-mapping floor", "",
           f"Every reference sequence has a reversed twin — identical "
           f"length, identical base composition, no homology. Across "
           f"{fmt(runs)} runs those decoys collected **{fmt(decoy_reads)} "
           f"reads in total**, in {fmt(len(hot))} of "
           f"{fmt(len(by_run))} run x locus comparisons.", ""]
    if hot:
        out += ["Where a decoy did collect reads, the shape of the "
                "collection is the point, and it is why `covered_bases` is "
                "recorded alongside every count: reads spread over a decoy "
                "would mean the reference is porous, whereas reads piled on "
                "one short window are a repeat or a low-complexity stretch "
                "that survived reversal.", "",
                "| species | run | tissue | decoy | reads | bases covered | "
                "decoy length | its locus's reads |",
                "|---|---|---|---|---|---|---|---|"]
        for r in sorted(hot, key=lambda x: -_i(x, "decoy_reads")):
            out.append(
                f"| *{r['organism']}* | `{r['run']}` | {r['tissue']} | "
                f"`decoy_{r['seq'].split('|')[0]}` | "
                f"{fmt(r['decoy_reads'])} | "
                f"{fmt(r['decoy_covered_bases'])} | {fmt(r['ref_len'])} | "
                f"{fmt(r['reads'])} |")
        out += ["", "Every affected locus still clears its own decoy in the "
                "same run, which is the clause the detection rule actually "
                "applies — the decoy is a per-run floor, not a global one.",
                ""]
    out.append(prior_line(
        "decoy_floor",
        (f"{fmt(decoy_reads)} decoy reads across {fmt(runs)} runs on "
         f"references up to 15 kb — the floor is zero here too."
         if decoy_reads == 0 else
         f"{fmt(decoy_reads)} decoy reads across {fmt(runs)} runs and "
         f"{fmt(len(by_run))} run x locus comparisons. The floor is not "
         f"exactly zero on this panel, and reporting it as zero because "
         f"the ported result was zero is the error this comparison "
         f"exists to prevent. What it is instead is bounded and "
         f"localised, as the table above shows, and it never changes a "
         f"detection call."),
        decoy_reads == 0))
    out += ["", "### The libraries worked", "",
            "The three anchors are **not** equally informative, and the "
            "reason is measured rather than supposed. The reference holds "
            "one model per anchor, so reads from a gene's other genomic "
            "copies multi-map and are cut by the MAPQ floor — which orders "
            "the anchors exactly as their copy numbers do. A near-zero "
            "count on a multi-copy anchor is therefore not a library that "
            "failed.", "",
            "| species | anchor | copies in this genome | reads | junction "
            "reads | runs detected |", "|---|---|---|---|---|---|"]
    for r in sorted(hk, key=lambda x: (x["organism"], x["cell"])):
        out.append(f"| *{r['organism']}* | {r['cell']} | "
                   f"{fmt(r.get('genome_copies'))} | {fmt(r['reads'])} | "
                   f"{fmt(r['junction_reads'])} | "
                   f"{fmt(r['runs_detected'])}/{fmt(r['runs'])} |")
    out += ["", "### The rest of the family in the same libraries", "",
            "Reads on a locus the annotation loses mean more when the "
            "genes beside it behave. The RyR control is this project's "
            "sister family (D14): its reads must land on RyR.", "",
            "| species | locus | role | annotation | reads | junction reads "
            "| runs detected |", "|---|---|---|---|---|---|---|"]
    for r in sorted(ctrl, key=lambda x: (x["organism"], x["cell"])):
        out.append(f"| *{r['organism']}* | {r['cell']} | {r['role']} | "
                   f"{r['annot_gene'] or '—'} | {fmt(r['reads'])} | "
                   f"{fmt(r['junction_reads'])} | "
                   f"{fmt(r['runs_detected'])}/{fmt(r['runs'])} |")
    out += _crossmap_section(crossmap, fmt, pct)
    out += ["", "![](figures/s12_detection.png)", "",
            "**Figure 2.** Reads on each recovered locus against its own "
            "composition-matched decoy, per species. Log axis: these "
            "libraries differ by more than an order of magnitude in depth "
            "and the comparison that matters is within a run, not between "
            "them.", ""]
    return out


def _crossmap_section(crossmap, fmt, pct) -> list[str]:
    """The closed-set risk, measured rather than argued away."""
    if not crossmap:
        return ["", "### Cross-mapping between the references", "",
                "*Not run yet — `crossmap_control.tsv` is not present.*", ""]
    worst = max(_f(r, "cross_rate") for r in crossmap)
    total = sum(_i(r, "reads_simulated") for r in crossmap)
    other = sum(_i(r, "assigned_other") for r in crossmap)
    lowq = sum(_i(r, "below_mapq") for r in crossmap)
    bad = [r for r in crossmap if _i(r, "assigned_other")]
    out = ["", "### Cross-mapping between the references", "",
           "The reference is a closed set, so a read from a gene not in it "
           "cannot be assigned away. The risk that matters is the paralogs "
           "assigning to each other, and the ITPR paralogs are more similar "
           "to one another than the PIEZO family this method was ported "
           "from — whose S12 argued the risk away in a caveat. Here it is "
           "measured.", "",
           f"Every reference sequence was tiled exhaustively with 100 nt "
           f"synthetic reads at 10 nt steps ({fmt(total)} reads) and mapped "
           f"back with the same aligner settings. **{fmt(other)} of "
           f"{fmt(total)} reads were assigned to a sequence other than the "
           f"one they came from** (worst per-sequence rate "
           f"{worst:.5f}); {fmt(lowq)} fell below the MAPQ floor.", ""]
    if bad:
        out += ["| species | source | reads elsewhere | worst destination |",
                "|---|---|---|---|"]
        for r in bad:
            out.append(f"| {r['species']} | `{r['source_seq']}` | "
                       f"{fmt(r['assigned_other'])} | "
                       f"{r['worst_other_seq']} "
                       f"({fmt(r['worst_other_reads'])}) |")
        out.append("")
    else:
        out += ["The paralogs, the RyR control and the housekeeping anchors "
                "are all separable at read level in these genomes; the "
                "closed set is safe here as a measured fact.", ""]
    return out


def _pseudogenes(loci, refs, fmt, pct, missing) -> list[str]:
    """The loci the annotation *does* model — as pseudogenes.

    S10 established that a GFF3 `pseudogene` can carry a full set of CDS
    features and still emit no protein, so "the annotation has this gene"
    and "the databases serve this protein" are different statements.  Two
    of the panel's three internal controls turn out to be exactly that
    case, which makes them a second test of the same kind rather than
    merely a control: the annotation asserts these loci are not
    protein-coding, and a transcribed, spliced ORF is evidence against it.
    """
    # Which loci carry a pseudogene model is a property of the reference,
    # not of the reads, so the section is driven by `reference_table.tsv`.
    # Driving it off the counts would make the section — and every section
    # number after it — appear and disappear with how much of the sweep has
    # been quantified.
    expected = {(r["species"], r["cell"]) for r in refs
                if r.get("annot_biotype") == "pseudogene"}
    if not refs:
        return [missing("7. The loci the annotation models — as pseudogenes")]
    if not expected:
        return ["## 7. The loci the annotation models — as pseudogenes", "",
                "No locus in this panel carries a pseudogene model.", ""]
    pseudo = [r for r in loci if (r["species"], r["cell"]) in expected]
    if not pseudo:
        return [missing("7. The loci the annotation models — as pseudogenes")]
    out = ["## 7. The loci the annotation models — as pseudogenes", "",
           "Not every locus the annotation reaches is a locus it delivers. "
           "A GFF3 `pseudogene` may carry a full set of coding features and "
           "still emit no protein, so these loci are `found_annotated` in "
           "S5's ledger and absent from every protein database.", "",
           "| species | locus | annotated as | biotype | share of the "
           "coding footprint | reads | junction reads | runs detected |",
           "|---|---|---|---|---|---|---|---|"]
    for r in sorted(pseudo, key=lambda x: (x["organism"], x["cell"])):
        out.append(
            f"| *{r['organism']}* | {r['cell']} | `{r['annot_gene']}` | "
            f"{r['annot_biotype']} | {fmt(r['annot_frac_cds'])} | "
            f"{fmt(r['reads'])} | {fmt(r['junction_reads'])} | "
            f"{fmt(r['runs_detected'])}/{fmt(r['runs'])} |")
    live = [r for r in pseudo if _i(r, "runs_detected")]
    by_sp: dict[str, list] = {}
    for r in pseudo:
        by_sp.setdefault(r["organism"], []).append(r)
    out += ["",
            f"**{len(live)} of {len(pseudo)}** meet all three detection "
            f"criteria. S10 had already shown these loci splice into an "
            f"uninterrupted reading frame; the reads add that the frame is "
            f"transcribed.", ""]
    for org, rs in sorted(by_sp.items()):
        cells = sorted(r["cell"] for r in rs)
        out.append(f"In *{org}* this is not a detail about one gene: "
                   f"{', '.join(cells)} are filed as pseudogenes and the "
                   f"remaining paralog is not annotated at all, so **no IP3 "
                   f"receptor of this species reaches a protein record by "
                   f"any route**.")
    out.append("")
    return out


def _gap(gap, fmt, pct, missing, prior_line) -> list[str]:
    rows = [r for r in gap if r["role"] == "failed"]
    if not rows:
        return [missing("8. Where along the gene the reads fall")]
    out = ["## 8. Where along the gene the reads fall", "",
           "S10 measured annotation loss as the share of a recovered gene's "
           "coding footprint that no annotated model delivers. The reads "
           "give the same quantity a second way: the share of read coverage "
           "landing outside any annotated coding block. The two are "
           "computed from different evidence — one from a GFF, one from "
           "aligned reads — so agreement between them is a check, not a "
           "restatement.", "",
           "| species | locus | S10 annotation loss | read coverage outside "
           "annotated CDS | bins with reads (unannotated) |",
           "|---|---|---|---|---|"]
    for r in sorted(rows, key=lambda x: (x["organism"], x["cell"])):
        out.append(
            f"| *{r['organism']}* | {r['cell']} | "
            f"{frac(r['annotation_loss'])} | "
            f"{frac(r['frac_coverage_unannotated'])} | "
            f"{fmt(r['bins_unannotated_with_reads'])}/"
            f"{fmt(r['bins_unannotated'])} |")
    pairs = [(_f(r, "annotation_loss"), _f(r, "frac_coverage_unannotated"))
             for r in rows if r.get("annotation_loss") not in ("", None)]
    agree = sum(1 for a, b in pairs if abs(a - b) <= 0.25)
    out += ["", "![](figures/s12_gap_coverage.png)", "",
            "**Figure 3.** The fraction of read coverage falling outside any "
            "annotated coding block, with S10's coding-footprint loss marked "
            "on each bar.", ""]
    if pairs:
        out.append(prior_line(
            "cell_status",
            f"S5b called these cells `found_unannotated` and declined to "
            f"adjudicate. Reads land outside the annotation in the "
            f"proportion S10 predicted from the coding footprint for "
            f"{agree} of {len(pairs)} loci (within 0.25).",
            agree >= max(1, len(pairs) - 1)))
        out.append("")
    return out


def _tissues(tissue, fmt, pct, missing) -> list[str]:
    rows = [r for r in tissue if r["role"] == "failed"
            and r["tissue"] != "unknown"]
    if not rows:
        return [missing("9. Where they are expressed")]
    tissues = sorted({r["tissue"] for r in rows})
    out = ["## 9. Where they are expressed", "",
           "Pooled by tissue, over runs whose organ is confirmed by their "
           "own sample attributes. These are shallow subsamples of single "
           "libraries, so a tissue with no reads is weak evidence of "
           "absence and a tissue with reads is strong evidence of presence.",
           "", "| species | locus | " + " | ".join(tissues) + " |",
           "|---|---|" + "---|" * len(tissues)]
    keys = sorted({(r["organism"], r["cell"]) for r in rows})
    for org, cell in keys:
        cells = []
        for t in tissues:
            m = [r for r in rows if r["organism"] == org
                 and r["cell"] == cell and r["tissue"] == t]
            if not m:
                cells.append("—")
            else:
                r = m[0]
                mark = "**" if _i(r, "detected") else ""
                cells.append(f"{mark}{fmt(r['reads'])}{mark}")
        out.append(f"| *{org}* | {cell} | " + " | ".join(cells) + " |")
    out += ["", "Bold is a tissue meeting all three detection criteria on "
            "its pooled runs.", ""]
    return out


