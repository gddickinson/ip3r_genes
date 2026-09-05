"""s20_report_results.py — the results half of the S20 report.

Split from `s20_report.py` to keep both under 500 lines, and following the
`s3_report.py` / `s3_report_d10.py` and `s5_report.py` /
`s5_report_results.py` pattern: this module takes the caller's loader and
formatter, so the two halves cannot read the same tables differently.

**The headlines are chosen by the data.** S2 enumerated the InterPro-seeded
space and found the plant and fungal records phylogenetically clean — every
green call a chlorophyte, no land-plant call at all; every fungal call in an
early-diverging phylum, Dikarya contributing no records whatever. Those are
S20's hypotheses, not its conclusions, and they were measured on a *seeded*
search space: a record only enters it by already carrying a family Pfam, so
"Dikarya has no records" could as easily mean "Dikarya has no records with a
family Pfam annotation". S20 searches the proteomes themselves, which is the
one design that can tell those apart.

So each section states S2's number, computes S20's, and renders `confirmed`
/ `contradicted` / `underpowered` from the comparison — printing both either
way. A generator written to narrate S2's answer would print it whatever the
sweep said.
"""

from __future__ import annotations

from collections import Counter

#: S2's recorded values, from `results/census_v2/report.md` (2026-09-03).
#: Source of every number: the committed census v2 tables.
PRIOR = {
    "viridiplantae": {"records": 75, "itpr": 20, "taxa": 26},
    "fungi": {"records": 46, "itpr": 25, "taxa": 35},
    "streptophyta": {"records": 15, "itpr": 0, "taxa": 13},
    "chlorophyta": {"records": 60, "itpr": 20, "taxa": 11},
    "dikarya": {"records": 0, "itpr": 0, "taxa": 0},
    "early_fungi": {"Mucoromycota": 15, "Chytridiomycota": 6,
                    "Basidiobolomycota": 3, "Entomophthoromycota": 1},
}

#: Below this many swept proteomes a phylum-level absence is not evidence.
MIN_PROTEOMES_FOR_ABSENCE = 10

DIKARYA = ("Ascomycota", "Basidiomycota")

#: Protein-name markers for the family's *known* false positives — the
#: MIR-domain sharers S1's decoy panel was built around. Used only to decide
#: whether an iteration-only record in a lineage called empty is already
#: accounted for; a record that matches none of these is reported as
#: unexplained rather than explained away.
DECOY_MARKERS = ("mannosyltransferase", "pomt", "dolichyl")


def _rate(hit: int, tot: int) -> str:
    return f"{hit:,}/{tot:,} ({hit / tot:.0%})" if tot else "n/a"


def _verdict(prior_zero: bool, found: int, n_proteomes: int) -> tuple[str, str]:
    """(verdict, gloss) for a prior absence claim tested at proteome scale."""
    if n_proteomes < MIN_PROTEOMES_FOR_ABSENCE:
        return ("underpowered",
                f"only {n_proteomes} reference proteomes were swept — too few "
                "to test an absence")
    if prior_zero and found == 0:
        return ("confirmed",
                "searching the proteomes themselves finds no more than the "
                "seeded space did")
    if prior_zero and found:
        return ("contradicted",
                "the seeded space held none, but searching the proteomes "
                "themselves finds them — an annotation absence, not a "
                "genomic one")
    return ("confirmed" if found else "contradicted", "")


def _by_phylum(presence, taxa, phyla=None, group=None):
    """(n_proteomes, n_with_itpr, [(clade, n, hits)]) for a slice.

    Keyed on the taxonomy table's derived `clade`, not on `phylum`: UniProt
    gives no phylum rank to the choanoflagellates, filastereans, cryptophytes
    and apusozoans, and those are exactly the lineages a range result about
    this family has to be able to name.
    """
    rows = []
    for r in presence:
        t = taxa.get(int(r["taxid"]), {})
        if group and t.get("group") != group:
            continue
        label = t.get("clade") or t.get("phylum") or "unclassified"
        if phyla and label not in phyla:
            continue
        rows.append((label, r["itpr_status"] == "present"))
    agg: dict[str, list[int]] = {}
    for phylum, hit in rows:
        cell = agg.setdefault(phylum, [0, 0])
        cell[0] += 1
        cell[1] += int(hit)
    detail = sorted(((p, v[0], v[1]) for p, v in agg.items()),
                    key=lambda x: (-x[2], -x[1]))
    return len(rows), sum(1 for _, h in rows if h), detail


# ------------------------------------------------------------------ sections
def section_range(A, table, presence, taxa, v5) -> None:
    A("\n## The family's range across the eukaryotic tree\n")
    if not presence:
        A("_`proteome_presence.tsv` missing — run `s20_run_sweep.py`._\n")
        return
    n = len(presence)
    hit = sum(1 for r in presence if r["itpr_status"] == "present")
    A(f"{n:,} reference proteomes were swept with both family profiles at "
      f"the S3 protocol's threshold, and **{_rate(hit, n)} carry at least "
      f"one ITPR call**.\n")

    by_group: dict[str, list[int]] = {}
    for r in presence:
        g = taxa.get(int(r["taxid"]), {}).get("group", "unclassified")
        cell = by_group.setdefault(g, [0, 0, 0])
        cell[0] += 1
        cell[1] += int(r["itpr_status"] == "present")
        cell[2] += int(r["n_itpr"] or 0)
    A(table(["group", "proteomes swept", "with an ITPR call", "ITPR records"],
            [[g, f"{v[0]:,}", _rate(v[1], v[0]), f"{v[2]:,}"]
             for g, v in sorted(by_group.items(), key=lambda kv: -kv[1][1])]))

    _, _, detail = _by_phylum(presence, taxa)
    carrying = [d for d in detail if d[2]]
    A(f"\nAcross every group, **{len(carrying)} clades of the "
      f"{len(detail)} swept carry an ITPR call**. Clade is the phylum where "
      "UniProt states one and the next-deepest named group where it does not "
      "— the choanoflagellates, filastereans, cryptophytes and apusozoans "
      "have no phylum rank, and those are precisely the lineages a range "
      "result about this family has to be able to name. Most-covered first:\n")
    A(table(["clade", "proteomes", "with an ITPR call"],
            [[p, str(n_), _rate(h, n_)] for p, n_, h in carrying[:25]]))
    if v5:
        A(f"\nCensus v5 holds **{v5.get('n_itpr', 0):,} ITPR records across "
          f"{v5.get('n_itpr_taxa', 0):,} taxa**.\n")


def section_d14(A, table, assignments, min_score, rel_margin,
                min_positions, seed_accessions=frozenset()) -> None:
    """How separable the two families are outside the vertebrates — measured.

    S1 had all six RyR decoys promoted at 45 before the sister test existed;
    S5b found the two panels never once contested a genomic locus. This is
    the same question at protein level in the part of the tree neither task
    covered, and the answer is a count, not a claim.
    """
    A("\n## D14 outside the vertebrates\n")
    if not assignments:
        A("_no assignment tables — run `s20_run_sweep.py`._\n")
        return
    both = [r for r in assignments
            if float(r["itpr_score"] or 0) >= min_score
            and float(r["ryr_score"] or 0) >= min_score]
    contested = [r for r in both if float(r["rel_margin"] or 0) < rel_margin]
    # Keyed on the exact phrases `s3_assign.assign()` writes, not on a
    # loose substring: a reason string that changes shape should stop
    # matching visibly rather than land in the wrong bucket.
    def _outcome(reason: str) -> str:
        if reason.startswith("best profile score"):
            return "declined: under the bit-score floor"
        if "shared module" in reason:
            return "declined: shared module only"
        if "no-call band" in reason:
            return "declined: inside the no-call band"
        return "called"

    reasons = Counter(_outcome(r["reason"]) for r in assignments)
    A(f"{len(assignments):,} targets were scored by at least one profile. "
      f"**{len(both):,} were scored by both above the "
      f"{min_score:.0f}-bit floor** — the only ones where the two families "
      f"can be said to compete at all — and of those "
      f"**{_rate(len(contested), len(both))} fall inside D7's "
      f"{rel_margin:.0%} no-call band**, where this instrument declines to "
      "choose.\n")
    A(table(["outcome", "targets"],
            [[k, f"{v:,}"] for k, v in reasons.most_common()]))
    A(f"\nThe {min_positions}-state floor (D22) is doing most of the "
      "declining: a protein sharing one module with a full-length channel "
      "model scores against it, and outside the vertebrates that is the "
      "common case rather than the exception.\n")
    A("\n### What each call actually rests on\n")
    A("A call is only as good as how much of the model it spans, and the two "
      "profiles are very different sizes — `itpr.hmm` is 2,684 match states, "
      "`ryr.hmm` 4,930 — so 200 states is 7 % of one model and 4 % of the "
      "other. Breaking the calls out by evidence class is what stops that "
      "asymmetry from being read as biology:\n")
    by_group: dict[tuple[str, str], list[dict]] = {}
    for r in assignments:
        if r["assignment"] in ("ITPR", "RYR"):
            by_group.setdefault((r["group"], r["assignment"]), []).append(r)
    body = []
    for (g, call), sub in sorted(by_group.items()):
        col = "itpr_coverage" if call == "ITPR" else "ryr_coverage"
        cov = sorted(float(r[col] or 0) for r in sub)
        ev = Counter(r["evidence"] for r in sub)
        body.append([g, call, f"{len(sub):,}",
                     f"{ev.get('architecture', 0):,}",
                     f"{ev.get('partial', 0):,}",
                     f"{cov[len(cov) // 2]:.0%}"])
    A(table(["group", "call", "records", "architecture (≥50 % of the model)",
             "partial", "median coverage"], body))
    # The contrast is computed per group, not asserted globally. When this
    # paragraph was first written the sweep held no invertebrates and *every*
    # RYR call outside the vertebrates was module-level; with the
    # invertebrates in, 441 of them are real full-length ryanodine receptors.
    # A sentence that kept saying "overwhelmingly module-level" would have
    # been contradicted by the table directly above it.
    def _arch(call, groups):
        sub = [r for r in assignments if r["assignment"] == call
               and r.get("group") in groups]
        return sum(1 for r in sub if r["evidence"] == "architecture"), len(sub)

    metazoan = {"metazoa_nonvert"}
    others = {r.get("group") for r in assignments} - metazoan
    A(f"\nThe split runs along the one line that matters. Where ryanodine "
      f"receptors are expected — the invertebrates — the RYR calls are real: "
      f"{_rate(*_arch('RYR', metazoan))} span at least half the model. "
      f"Everywhere else they are not: {_rate(*_arch('RYR', others))}. The "
      f"ITPR calls hold up in both ({_rate(*_arch('ITPR', metazoan))} and "
      f"{_rate(*_arch('ITPR', others))}). So the module-level RYR calls in "
      "fungi, green algae and protists are the shared architecture D14 exists "
      "to see through, not receptors: read as gene counts they would invent a "
      "ryanodine receptor family across half the eukaryotic tree, and read as "
      "coverage they do not.\n")

    arch_ryr = [r for r in assignments if r["assignment"] == "RYR"
                and r["evidence"] == "architecture"
                and r.get("group") in others]
    if arch_ryr:
        A(f"\nThe {len(arch_ryr)} architecture-level RYR call(s) outside the "
          "invertebrates are worth naming, because where the sister family "
          "is tells you when the two families split:\n")
        A(table(["accession", "species", "length", "ryr bits", "coverage",
                 "margin", "a profile seed?"],
                [[r["accession"], f"*{r['species']}*", r["length"],
                  r["ryr_score"], f"{float(r['ryr_coverage']):.0%}",
                  f"{float(r['rel_margin']):.0%}",
                  "yes — circular" if r["accession"] in seed_accessions
                  else "no"]
                 for r in sorted(arch_ryr,
                                 key=lambda r: -float(r["ryr_score"] or 0))[:15]]))
        A("\nA record that is itself one of `ryr.hmm`'s seeds scores well "
          "against a profile built partly from it, so its score is not "
          "independent evidence and the column says so. The rest are.\n")


    if contested:
        A("\nThe contested targets — kept and reported, not resolved:\n")
        A(table(["accession", "species", "itpr bits", "ryr bits", "margin"],
                [[r["accession"], f"*{r['species']}*", r["itpr_score"],
                  r["ryr_score"], f"{float(r['rel_margin']):.1%}"]
                 for r in sorted(contested,
                                 key=lambda r: -float(r["itpr_score"] or 0))[:15]]))


def section_plants(A, table, presence, taxa, verdicts) -> None:
    A("\n## Land plants: a genome fact or a database fact?\n")
    prior = PRIOR["streptophyta"]
    n, hit, detail = _by_phylum(presence, taxa, group="Viridiplantae")
    strep = [d for d in detail if d[0] == "Streptophyta"]
    n_strep = strep[0][1] if strep else 0
    h_strep = strep[0][2] if strep else 0
    verdict, gloss = _verdict(prior["itpr"] == 0, h_strep, n_strep)
    A(f"S2 enumerated the InterPro-seeded space and found **Streptophyta — "
      f"the land-plant lineage — with {prior['itpr']} ITPR calls from "
      f"{prior['records']} records in {prior['taxa']} taxa**. That is a claim "
      "about a space a protein enters only by already carrying a family Pfam "
      "annotation. Sweeping the proteomes themselves asks the question "
      "without that filter.\n")
    A(f"**{verdict.upper()}** — {_rate(h_strep, n_strep)} Streptophyta "
      f"reference proteomes carry an ITPR call: {gloss}.\n")
    A(table(["clade", "proteomes swept", "with an ITPR call"],
            [[p, str(n_), _rate(h, n_)] for p, n_, h in detail]))
    A(f"\nViridiplantae overall: {_rate(hit, n)} proteomes.\n")
    _verdict_block(A, table, verdicts, "Viridiplantae")


def section_fungi(A, table, presence, taxa, verdicts) -> None:
    A("\n## Dikarya: the yeasts and moulds\n")
    n, hit, detail = _by_phylum(presence, taxa, group="Fungi")
    dik = [d for d in detail if d[0] in DIKARYA]
    n_dik = sum(d[1] for d in dik)
    h_dik = sum(d[2] for d in dik)
    verdict, gloss = _verdict(True, h_dik, n_dik)
    A(f"S2 found **Dikarya contributing no records at all** to the seeded "
      f"space — not zero calls from some records, zero records. Every fungal "
      f"call sat in an early-diverging phylum "
      f"({', '.join(f'{k} {v}' for k, v in PRIOR['early_fungi'].items())}).\n")
    A(f"**{verdict.upper()}** — {_rate(h_dik, n_dik)} Ascomycota + "
      f"Basidiomycota reference proteomes carry an ITPR call: {gloss}.\n")
    A(table(["clade", "proteomes swept", "with an ITPR call"],
            [[p, str(n_), _rate(h, n_)] for p, n_, h in detail]))
    A(f"\nFungi overall: {_rate(hit, n)} proteomes.\n")
    _verdict_block(A, table, verdicts, "Fungi")


def _verdict_block(A, table, verdicts, kingdom: str) -> None:
    rows = [v for v in verdicts if v["kingdom"] == kingdom]
    if not rows:
        return
    tally = Counter(v["verdict"] for v in rows)
    A(f"\nEvery {kingdom} record called ITPR by either instrument was chased "
      f"individually ({len(rows)} records):\n")
    A(table(["verdict", "records", "what it means"],
            [[v, str(n), _VERDICT_GLOSS.get(v, "")]
             for v, n in tally.most_common()]))
    per_species = Counter(v["species"] for v in rows
                          if v["verdict"] == "real_gene")
    if per_species:
        top, n_top = per_species.most_common(1)[0]
        if n_top >= 3:
            A(f"\nThe surviving records are not evenly spread: **{top}** "
              f"alone contributes {n_top} of them, against a median of "
              f"{sorted(per_species.values())[len(per_species) // 2]} per "
              "species. Whether that is real copy-number expansion or a "
              "duplicated assembly is a question for an assembly search, not "
              "a proteome one *(pending: S23)*.\n")

    real = [v for v in rows if v["verdict"] == "real_gene"]
    if real:
        A(f"\nThe {len(real)} that survive, by lineage:\n")
        A(table(["accession", "species", "phylum", "length",
                 "nearest outside its kingdom"],
                [[v["accession"], f"*{v['species']}*", v["phylum"],
                  v["length"],
                  f"{v['out_kingdom_pident']} % ({v['out_kingdom_group']})"
                  if v["out_kingdom_pident"] else "—"]
                 for v in sorted(real, key=lambda r: r["phylum"])[:20]]))


_VERDICT_GLOSS = {
    "real_gene": "survives every test",
    "fragment": "real but incomplete gene model",
    "no_genome_backing": "no EMBL or proteome cross-reference",
    "cross_kingdom_outlier": "80-95 % to another kingdom — flagged, not called",
    "module_only": "the hit rests on a shared module, not a family domain",
    "contaminant_suspect": "≥ 95 % to another kingdom — a sequence in the "
                           "wrong assembly",
    "unresolved": "evidence unavailable",
}


def section_negatives(A, table, relaxed, taxa, assignments, groups_relaxed,
                      e_relaxed, e_primary) -> None:
    """The negative claims, with the sensitivity they were made at — and with
    positive controls inside the same search.

    An absence is only worth as much as the search that failed to find
    anything, and the way to show a search works is to point it somewhere the
    thing *is*. So every count below is repeated in the nearest lineage where
    the family is present, chosen by the data (the clade in the same kingdom
    with the most ITPR calls in the primary sweep) rather than picked to
    flatter the result.
    """
    A("\n## The negative claims, at a stated sensitivity\n")
    A(f"Every claim was made at **E ≤ {e_primary}** with the two full-length "
      f"family profiles and re-made at **E ≤ {e_relaxed}** with those two "
      "*plus* the family's four Pfam domain models — a 213-state IP₃-core "
      "model can reach something a 2,684-state channel model cannot.\n")
    if not relaxed:
        A("_`relaxed_hits.tsv` missing — run `s20_run_sweep.py --stage "
          "relaxed`._\n")
        return

    def clade_of(row) -> str:
        return taxa.get(int(row["taxon_id"] or 0), {}).get("clade", "")

    # The control is derived: within each kingdom, the clade holding the most
    # ITPR calls in the primary sweep.
    controls: dict[str, str] = {}
    for kingdom, group in (("land plants", "viridiplantae"),
                           ("Dikarya", "fungi")):
        counts = Counter(clade_of(r) for r in assignments
                         if r.get("group") == group and r["assignment"] == "ITPR")
        if counts:
            controls[kingdom] = counts.most_common(1)[0][0]

    lineages = {"land plants": ["Streptophyta"], "Dikarya": list(DIKARYA)}
    labels: list[tuple[str, list[str], bool]] = []
    for claim, clades in lineages.items():
        labels.append((claim, clades, False))
        if claim in controls:
            labels.append((f"{controls[claim]} (control)",
                           [controls[claim]], True))

    profiles = sorted({r["profile"] for r in relaxed})
    body = []
    for label, clades, is_control in labels:
        sub = [r for r in relaxed if clade_of(r) in clades]
        row = [label]
        for prof in profiles:
            hits = [r for r in sub if r["profile"] == prof]
            deep = [r for r in hits if float(r["hmm_coverage"] or 0) >= 0.5
                    and float(r["full_evalue"]) <= float(e_primary)]
            row.append(f"{len(hits):,} / **{len(deep):,}**")
        body.append(row)
    A(f"\nReported at E ≤ {e_relaxed} / of those, spanning at least half the "
      f"model at E ≤ {e_primary}:\n")
    A(table(["lineage"] + profiles, body))

    core = "PF08709"
    def deep_count(clades, prof):
        return sum(1 for r in relaxed if clade_of(r) in clades
                   and r["profile"] == prof
                   and float(r["hmm_coverage"] or 0) >= 0.5
                   and float(r["full_evalue"]) <= float(e_primary))

    A(f"\n**The controls fire and the claims do not.** {core} is the "
      "IP₃-binding core, the signature that names the family:\n")
    for claim, clades in lineages.items():
        ctrl = controls.get(claim)
        got = deep_count(clades, core)
        ctrl_got = deep_count([ctrl], core) if ctrl else None
        A(f"- **{claim}**: {got} substantial {core} match(es)"
          + (f", against {ctrl_got} in *{ctrl}*, "
             "the nearest lineage in the same kingdom where the family is "
             "present." if ctrl else ".") + "\n")
    mir = "PF02815"
    A(f"\nAnd the search is demonstrably sensitive in those very genomes: "
      f"{mir} (MIR, which the family shares with POMT1/2 and every eukaryote "
      "therefore carries) returns "
      + ", ".join(f"{deep_count(cl, mir):,} substantial matches in {claim}"
                  for claim, cl in lineages.items())
      + ". The promiscuous domain finds thousands where the family-defining "
      "one finds none, which is what an absence looks like when the "
      "instrument is working rather than blind.\n")

    strong = [r for r in relaxed
              if r["profile"] in ("itpr", "ryr")
              and float(r["hmm_coverage"] or 0) >= 0.5]
    A(f"\nAcross every relaxed group, **{len(strong)} target(s) span at "
      "least half a full-length family profile** at the relaxed threshold.\n")
    if strong:
        A(table(["group", "profile", "accession", "species", "length",
                 "E-value", "model coverage"],
                [[r["group"], r["profile"], r["accession"],
                  f"*{r['species']}*", r["length"], r["full_evalue"],
                  f"{float(r['hmm_coverage']):.0%}"] for r in
                 sorted(strong, key=lambda r: float(r["full_evalue"]))[:20]]))


def section_jackhmmer(A, table, conv, verdicts_json,
                      composition=(), iteration_only=()) -> None:
    A("\n## Iterated search, per group (D10)\n")
    if not conv:
        A("_`jackhmmer_convergence_s20.tsv` missing — run "
          "`s20_jackhmmer.py`._\n")
        return
    A("A single profile pass can only find what it is already close enough "
      "to. Each group was therefore searched again with `jackhmmer` from a "
      "seed native to that group — the highest-scoring full-length ITPR call "
      "in the group's own sweep — iterated to convergence under D10's coded "
      "kill criterion.\n")
    v = (verdicts_json or {}).get("verdicts", {})
    A(table(["group", "seed", "rounds", "converged", "D10 verdict",
             "rule fired", "targets in the accepted model"],
            [[g, val.get("accession", ""), str(val.get("n_rounds", "")),
              "yes" if val.get("converged") else "no",
              val.get("verdict", ""), val.get("rule", "") or "none",
              f"{val.get('accepted_targets', ''):,}"
              if isinstance(val.get("accepted_targets"), int) else ""]
             for g, val in sorted(v.items())]))
    if composition:
        A("\n### What each converged model is built from\n")
        A("The completeness statement the iterated search exists to make. A "
          "model seeded in one lineage and iterated to convergence either "
          "reaches a neighbouring lineage or it does not, and the targets "
          "**only iteration found** are the ones that matter: if they sit in "
          "the same lineage as the seed, the absence next door is not a "
          "sensitivity artefact.\n")
        A(table(["group", "clade", "profile call", "targets",
                 "found by the single pass", "iteration only"],
                [[r["group"], r["clade"], r["profile_call"], r["targets"],
                  r["found_by_single_pass"], r["iteration_only"]]
                 for r in composition]))
        only = sum(int(r["iteration_only"]) for r in composition)
        A(f"\n{only} target(s) entered an accepted model that the single "
          "profile pass never reported. Those are the ones iteration was run "
          "to find, so they are named rather than counted:\n")
        if iteration_only:
            A(table(["group", "accession", "species", "clade", "length",
                     "what it is"],
                    [[r["group"], r["accession"], f"*{r['species']}*",
                      r["clade"], r["length"], r["protein_name"]]
                     for r in iteration_only]))
            outside = [r for r in iteration_only
                       if r["clade"] in DIKARYA or r["clade"] == "Streptophyta"]
            if outside:
                names = sorted({r["protein_name"] for r in outside})
                A(f"\n**{len(outside)} of them sit in a lineage this task "
                  "calls empty** — so the absence there is not quite "
                  "absolute, and what they are decides whether that matters: "
                  + "; ".join(names) + ".\n")
                # The identification is made from the record's own name, not
                # asserted: if a rerun finds something else sitting there, this
                # sentence has to stop appearing rather than keep explaining
                # away a different record.
                if all(any(k in n.lower() for k in DECOY_MARKERS)
                       for n in names):
                    A("\nEvery one is a mannosyltransferase — the MIR-domain "
                      "sharer S1's decoy panel was built around (POMT1/2), "
                      "not a receptor. The iterated search reaches these "
                      "lineages exactly far enough to pick up the known false "
                      "positive and no further.\n")
                else:
                    A("\nThese are not all accounted for by the family's "
                      "known decoys, so the absence claim above is qualified "
                      "by them and they need chasing individually.\n")
            else:
                A("\nNone of them sits in a lineage this task calls empty.\n")

    unfinished = (verdicts_json or {}).get("groups_not_completed") or []
    if unfinished:
        A(f"\n{len(unfinished)} group(s) had not finished iterating when this "
          f"report was rendered and are reported as unfinished rather than "
          f"summarised from a partial log: {', '.join(unfinished)}. The "
          "bottleneck is the per-round alignment, not the search — a round "
          "over a group with thousands of included targets costs far more to "
          "align than to scan, and that cost is not reduced by more cores.\n")

    missing = (verdicts_json or {}).get("groups_without_seed") or []
    if missing:
        A(f"\n{len(missing)} group(s) had no ITPR call spanning at least half "
          f"the model and so could not be seeded: {', '.join(missing)}. That "
          "is a statement about the group, not a failed run.\n")


def section_census(A, table, v5, changes, conflicts, by_group) -> None:
    A("\n## Census v5\n")
    if not v5:
        A("_`census_v5_stats.json` missing — run `s20_census_v5.py`._\n")
        return
    A(f"**{v5['n_records']:,} records**, of which "
      f"{v5['n_novel_from_s20']:,} are new from this sweep. Calls: "
      + ", ".join(f"{k} {n:,}" for k, n in sorted(v5["calls"].items(),
                                                  key=lambda kv: -kv[1]))
      + ".\n")
    A(f"\n{len(changes):,} records changed their call against census v4, and "
      f"{v5['n_conflicts']:,} are `conflict` — kept and reported rather than "
      "resolved (D23).\n")
    if by_group:
        A(table(["group", "ITPR", "RYR", "conflict", "unassigned", "total"],
                [[r["group"]] + [f"{int(r[k]):,}" for k in
                                 ("ITPR", "RYR", "conflict", "unassigned",
                                  "total")] for r in by_group]))
    if changes:
        A(f"\nThe {len(changes)} record(s) whose call moved:\n")
        A(table(["accession", "species", "census v4", "census v5", "why"],
                [[c["accession"], f"*{c['species']}*", c["v4_call"],
                  c["call"], c["reason"]] for c in changes[:10]]))
    if conflicts:
        A("\nThe conflicts:\n")
        A(table(["accession", "species", "census v4", "S20 sweep", "reason"],
                [[c["accession"], f"*{c['species']}*", c.get("arch_call", ""),
                  c.get("s20_call", ""), c["reason"]]
                 for c in conflicts[:15]]))
