"""s23_report_results.py — the results half of the S23 report.

Split from `s23_report.py` to keep both under 500 lines, following the
`s3_report.py` / `s3_report_d10.py` pattern: this module takes the caller's
table loader and formatter rather than re-importing its own, so the two halves
cannot read the same tables differently.

**Every headline here is chosen by the data.** S23a's pilot answered 14
genomes and S20a's absences were measured on *proteomes* — what a gene-caller
found, not what is in the DNA. Both are hypotheses at 194 genomes, and a
generator written to narrate them would print them whatever the sweep said.
So each section computes its statistic, compares it against the recorded prior
in `PRIOR`, and renders `confirmed` / `contradicted` / `underpowered` from the
comparison, printing both numbers either way.

The comparison that matters most is the one that could go badly: S20a's
"Streptophyta 0/384" and "Dikarya 0/1,353" are absences in *annotation*. If
this sweep finds a receptor in an assembly whose proteome had none, the
absence was a database fact and the headline changes. The report has to be
able to say that.
"""

from __future__ import annotations

from collections import Counter, defaultdict

#: What was believed before this sweep ran, recorded so it can be compared
#: rather than silently replaced.
#:   `s20_absences`  proteome-level zero-hit counts from S20a's presence table
#:                   (`results/s20_sweep/report.md`, 2026-09-05)
#:   `pilot_*`       the 14-genome S23a pilot (`results/s23_scope/report.md`)
PRIOR = {
    "s20_absences": {
        "Streptophyta": 384, "Ascomycota": 1034, "Basidiomycota": 319,
        "Apicomplexa": 60, "Microsporidia": 29, "Glomeromycota": 27,
        "Mortierellomycota": 18, "Kickxellomycota": 35,
        "Bacillariophyta": 16, "Rhodophyta": 12, "Cestoda": 11,
    },
    "pilot_genomes": 14,
    "pilot_found": 6,
    "pilot_no_locus": 8,
    "pilot_uncontrolled": 1,
    "pilot_uncontrolled_who": "Toxoplasma gondii",
    "pilot_max_copy_number": 1,
    "s5_min_locus_identity": 0.40,
    "s5_confirmed_floor": 0.759,
}

#: Below this many genomes in a clade, an absence claim about it is reported
#: with its n rather than given a verdict.
MIN_GENOMES_FOR_CLADE_VERDICT = 1


def _rate(hit: int, tot: int) -> str:
    return f"{hit}/{tot} ({hit / tot:.0%})" if tot else "n/a"


def _i(v, d: int = 0) -> int:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return d


def _f(v, d: float = 0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


# ------------------------------------------------------ the locus thresholds

def section_locus_calibration(A, table, cal: dict, spans: list[dict]) -> None:
    """The two thresholds S23b measured, and the one it could not."""
    A("\n## The locus thresholds, measured for this scope\n")
    if not cal:
        A("*Not yet measured — run `scripts/s23_calibrate_loci.py` after the "
          "sweep.*\n")
        return
    counts = cal.get("evidence_counts", {})
    conf = counts.get("confirmed", 0)
    A(f"S5b's `MIN_LOCUS_IDENTITY` = {PRIOR['s5_min_locus_identity']:.2f} came "
      f"from a wide empty gap in the vertebrates: confirmed loci from "
      f"{PRIOR['s5_confirmed_floor']:.3f} up, junk from 0.23 to 0.34, and "
      "every genome with a bait from its own class. Neither holds here — the "
      "bands are whole phyla and 13 bait slots are unfilled — so this sweep "
      f"**recorded** every cluster down to {cal.get('record_min_identity')} "
      "and the threshold was measured afterwards from what it recorded. A "
      "threshold measured from the population it has already filtered is "
      "circular; this one is not (D27).\n")
    A(f"Over {cal.get('genomes')} genomes the sweep recorded "
      f"{cal.get('loci_recorded')} clusters. Two evidence axes were used, and "
      "the first is thin for a reason worth stating: **only "
      f"{conf + counts.get('contradicted', 0) + counts.get('sister', 0)} of "
      f"{cal.get('loci_recorded')} clusters sit on a gene whose name says "
      "anything at all.** Outside the vertebrates most gene models carry "
      "locus tags — *Chlamydomonas* files its receptor as "
      "`CHLRE_16g665450v5`, *Strongylocentrotus* as `LOC594527` — so the "
      "evidence S5b calibrated 571 loci against does not exist at that depth "
      "in this scope.\n")
    A(table(["evidence axis", "confirmed", "contradicted"],
            [["the assembly's own annotation", conf,
              counts.get("contradicted", 0) + counts.get("sister", 0)],
             ["the S3 profiles (itpr.hmm / ryr.hmm)",
              counts.get("profile_confirmed", 0),
              counts.get("profile_contradicted", 0)],
             ["**pooled**", counts.get("confirmed_pooled", 0),
              counts.get("contradicted_pooled", 0)]]))

    seps = cal.get("separation") or {}
    if seps:
        A("\n**No threshold on either statistic separates them.**\n")
        A(table(["statistic", "confirmed reach down to", "contradicted up to",
                 "separates?"],
                [[k, f"{s.get('confirmed_min', 0):.3f}",
                  f"{s.get('contradicted_max', 0):.3f}",
                  "yes" if s.get("separates") else "**no — they overlap**"]
                 for k, s in seps.items()]))
    A(f"\n**`call_min_identity` = {cal.get('call_min_identity'):.2f}.** "
      f"{cal.get('call_min_identity_why')}\n")
    gate = cal.get("profile_gate_vs_annotation") or {}
    if gate:
        A(f"\n**The gate that replaces it, scored against the one axis it does "
          f"not share.** {gate.get('why')}.\n")
        for a in gate.get("admitted", []):
            A(f"The one exception is *{a['organism']}* gene `{a['gene']}` — "
              f"{a['coverage']} of its bait at {a['itpr_score']} bits against "
              f"{a['ryr_score']} for `ryr.hmm`. `IPR1` is not in the family "
              "name list, so the annotation axis scored it as naming a "
              "different gene; on the evidence it is more likely the name "
              "list is short than that the gate is wrong. It is reported as a "
              "conflict rather than reconciled — adding `ipr1` as a name "
              "substring would collide with every InterPro accession.\n")
    scan = cal.get("floor_scan") or []
    if scan:
        A("\nWhat each candidate floor would have cost, for the record:\n")
        A(table(["floor", "confirmed kept", "confirmed lost",
                 "contradicted kept", "no evidence either way"],
                [[f"{s['floor']:.2f}", s["confirmed_kept"],
                  s["confirmed_lost"], s["contradicted_kept"],
                  s.get("no_evidence_kept", s.get("unannotated_kept", ""))]
                 for s in scan]))
    misses = cal.get("near_miss_no_locus") or []
    A(f"\n**{len(misses)} `no_locus` genome(s) carry a cluster within 0.05 of "
      "the floor.** The brief asked for this number explicitly: a genome "
      "called absent because its best evidence scored just under the bar is a "
      "different claim from one whose best evidence was nowhere near it.\n")
    if misses:
        A(table(["identity", "organism", "group", "best bait"],
                [[f"{m['identity']:.3f}", f"*{m['organism']}*",
                  m.get("group", ""), m.get("bait", "")[:44]]
                 for m in misses[:12]]))
    inflate = cal.get("span_inflation_by_group") or {}
    if inflate:
        A("\n### A locus is much bigger than the gene inside it\n")
        A("`-G` is 650 kb for the metazoa, so a cluster is a large object. "
          "In the S23a pilot *Drosophila*'s 22 kb *Itpr* sat inside a "
          "297 kb cluster — harmless for a status call, fatal for a copy "
          "count, because two genes inside one such chain would be counted "
          "once. Copy number is therefore counted on **non-overlapping "
          "complete alignments**, not on clusters.\n")
        A(table(["group", "loci", "median span ÷ CDS", "max", "≥ 10×"],
                [[g, s["n"], f"{s['median']:.1f}×", f"{s['max']:.0f}×",
                  s["over_10x"]] for g, s in sorted(inflate.items())]))
        A(f"\n`copy_max_overlap` = {cal.get('copy_max_overlap'):.2f} — "
          f"{cal.get('copy_max_overlap_why')}\n")


# ------------------------------------------------------------- the control

def section_control(A, table, ledger: list[dict], controls: list[dict],
                    choice: list[dict], stats: dict) -> None:
    A("\n## The positive control\n")
    verdicts = Counter(r.get("verdict", "") for r in controls)
    admissible = sum(v for k, v in verdicts.items()
                     if k.startswith("controlled"))
    unc = verdicts.get("uncontrolled", 0)
    A("S5b could write \"the RyR control fired in all 309\" because every "
      "vertebrate has three RyRs. Outside Metazoa that sentence is not "
      "available (D26), and S23a's replacement — the MIR-domain sharer — has "
      "its own hole: both apicomplexan classes carry **one** PF02815 protein "
      "each across 60 swept proteomes, and the pilot's *Toxoplasma gondii* "
      "came back with neither a receptor nor a control.\n")
    A("**The fix is not a different profile, it is not fixing the profile in "
      "advance.** A control's job is to fire in the clade whose absence is "
      "the claim, so which family makes the best control is a property of the "
      "clade and is measurable. Six candidate profiles — all large, deeply "
      "conserved, multi-exon eukaryotic families — were run over each clade's "
      "own swept reference proteomes, and each clade takes the one that is "
      "actually there.\n")
    if choice:
        by_pfam = Counter(c.get("pfam", "") for c in choice if c.get("pfam"))
        A(table(["profile", "clades it controls", "example clade"],
                [[f"{p} ({next((c['label'] for c in choice if c.get('pfam') == p), '')})",
                  n_, next((c["clade"] for c in choice
                            if c.get("pfam") == p), "")]
                 for p, n_ in by_pfam.most_common()]))
        unfilled = [c for c in choice if not c.get("pfam")]
        if unfilled:
            A(f"\n{len(unfilled)} clade(s) have no candidate profile at all: "
              + ", ".join(f"*{c['clade']}*" for c in unfilled[:8]) + ".\n")
    A(f"\n**{admissible} of {len(controls)} genomes are controlled** "
      f"({_rate(admissible, len(controls))}); {unc} are not, and no absence "
      "claim rests on those.\n")
    A(table(["verdict", "genomes", "what it means"],
            [[f"`{k}`", v, _verdict_gloss(k)]
             for k, v in verdicts.most_common() if k]))
    xk = verdicts.get("controlled_cross_kingdom", 0)
    A(f"\n**{xk} genome(s) are controlled across a kingdom boundary** — a "
      "control bait from a different kingdom-level group aligns across the "
      "same locus. That is the tier that matters for an absence claim: it "
      "shows the search reaches across the divergence any receptor here "
      "would have to be found across, not merely that the assembly is "
      "readable.\n")
    prior_unc = PRIOR["pilot_uncontrolled"]
    A(f"*Against the pilot:* {prior_unc} of {PRIOR['pilot_genomes']} genomes "
      f"were uncontrolled there ({PRIOR['pilot_uncontrolled_who']}); here "
      f"{unc} of {len(controls)}.\n")


def _verdict_gloss(v: str) -> str:
    return {
        "controlled_by_target": "the ITPR baits found a full locus — the "
                                "search reached this assembly",
        "controlled_cross_kingdom": "a control locus, and a bait from another "
                                    "kingdom aligns across it",
        "controlled": "a control locus, from this genome's own clade only",
        "controlled_ryr_only": "no control locus, but a RyR where RyR is "
                               "expected",
        "no_control_bait": "the panel carries no control bait for this "
                           "genome's clade — a silent control, not an absent "
                           "result",
        "uncontrolled": "a control was available and none fired — supports no "
                        "absence claim",
    }.get(v, "")


# ------------------------------------------------------------- copy number

def section_copy_number(A, table, ledger: list[dict], stats: dict) -> None:
    A("\n## Copy number\n")
    ok = [r for r in ledger
          if (r.get("control_verdict") or "").startswith("controlled")]
    counts = Counter(_i(r.get("n_full")) for r in ok)
    A("Outside the vertebrates the paralog cells do not exist — ITPR1/2/3 are "
      "a 2R product and S5b found the trio absent below the cyclostomes — so "
      "the question is how many receptors a genome has.\n")
    A(f"Across the {len(ok)} controlled genomes:\n")
    A(table(["complete ITPR gene models", "genomes"],
            [[str(k) , v] for k, v in sorted(counts.items())]))
    multi = [r for r in ok if _i(r.get("n_full")) > 1]
    top = sorted(multi, key=lambda r: -_i(r.get("n_full")))[:12]
    if top:
        A(f"\n**{len(multi)} genome(s) carry more than one complete model.**\n")
        A(table(["organism", "group", "phylum", "copies", "clusters"],
                [[f"*{r['organism']}*", r.get("group", ""),
                  r.get("phylum", ""), r.get("n_full", ""),
                  r.get("n_full_loci", "")] for r in top]))
    diff = stats.get("genomes_copy_differs_from_loci", 0)
    A(f"\n**The copy rule changed the count in {diff} genome(s)** relative to "
      "counting full-graded clusters — the two columns above. Where they "
      "differ, a cluster either chained two genes or several baits hit one.\n")
    prior = PRIOR["pilot_max_copy_number"]
    hi = max(counts) if counts else 0
    A(f"*Against the pilot:* the 14 anchors topped out at {prior} copy; here "
      f"the maximum is {hi}.\n")


# ----------------------------------------------------------- the absences

def section_absences(A, table, absences: list[dict], ledger: list[dict],
                     stats: dict) -> None:
    A("\n## The absence claims, taken to assembly level\n")
    A("A proteome absence is a fact about what a gene-caller found. Only an "
      "assembly search makes it a fact about the genome. Each clade below "
      "counts **controlled genomes only** in both columns — an uncontrolled "
      "genome contributes to neither, which is the point of carrying a "
      "control.\n")
    rows = [a for a in absences if _i(a.get("genomes_swept"))]
    rows.sort(key=lambda a: -_i(a.get("swept_proteomes")))
    A(table(["clade", "proteomes (S20)", "genomes", "controlled",
             "with full ITPR", "trace only", "verdict"],
            [[a["clade"], a["swept_proteomes"], a["genomes_swept"],
              a["genomes_controlled"], a["genomes_with_full_itpr"],
              a["genomes_with_trace_only"], a["verdict"]] for a in rows]))
    found = [a for a in rows if _i(a.get("genomes_with_full_itpr"))]
    held = [a for a in rows if _i(a.get("genomes_controlled"))
            and not _i(a.get("genomes_with_full_itpr"))]
    A(f"\n**{len(held)} absence(s) hold at assembly level; {len(found)} do "
      "not.**\n")
    for a in found:
        prior = PRIOR["s20_absences"].get(a["clade"])
        A(f"- **{a['clade']} — contradicted.** S20 found 0 ITPR across "
          f"{prior if prior else a['swept_proteomes']} reference proteomes; "
          f"this sweep finds a complete gene model in "
          f"{a['genomes_with_full_itpr']} of {a['genomes_controlled']} "
          "controlled genomes. The proteome absence was an annotation fact.\n")
    for a in held:
        prior = PRIOR["s20_absences"].get(a["clade"])
        if prior:
            A(f"- **{a['clade']} — confirmed.** 0 ITPR across {prior} "
              f"reference proteomes (S20) and none in "
              f"{a['genomes_controlled']} controlled assembly search(es).\n")
    uncovered = [a for a in absences if not _i(a.get("genomes_controlled"))]
    if uncovered:
        A(f"\n{len(uncovered)} absence clade(s) have no controlled genome and "
          "therefore no assembly-level claim: "
          + ", ".join(f"*{a['clade']}*" for a in uncovered[:10]) + ".\n")


# ------------------------------------------------------------------- D14

def section_d14(A, table, loci: list[dict], controls: list[dict]) -> None:
    A("\n## D14 outside the vertebrates\n")
    itpr = [r for r in loci if r.get("grade") in ("full", "fragment")]
    margins = [_f(r.get("family_margin")) for r in itpr
               if r.get("family_margin") not in ("", None)]
    ctl_called = [r for r in loci if (r.get("role") or "") == "CONTROL"]
    A("The family call is a positive test on alignment score: a locus won by "
      "the RyR or the control baits is never offered to the family, and the "
      "margin says by how much the winner won.\n")
    if margins:
        margins.sort()
        A(f"**{len(margins)} graded ITPR loci**, family margin median "
          f"{margins[len(margins) // 2]:.3f}, minimum {margins[0]:.3f}. "
          f"{sum(1 for m in margins if m >= 0.99)} locus/loci have a margin "
          "of 1.0 — the RyR and control baits put no alignment there at "
          "all.\n")
    A(f"**{len(ctl_called)} control locus/loci were called ITPR.** A control "
      "protein the family call claims would be the sharpest possible failure "
      "of D14, and the panel carries one in every genome precisely so that "
      "failure has somewhere to show up.\n")
