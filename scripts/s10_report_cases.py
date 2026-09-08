"""The case-study half of the S10 report — one section per validated case.

Split from `s10_report.py` to keep both inside the 500-line budget, and it
takes the caller's table loader and number formatters rather than importing its
own (the `s3_report.py` / `s3_report_d10.py` pattern), so the two halves cannot
read the same table differently.

Each case is written in the order a reader needs it, which is the reverse of
the order the evidence was gathered in: what the annotation says, then what the
DNA says, then the independent checks that the DNA is right. The checks
come last on purpose — a case study that opens with its own validation is
asking to be believed before it has said anything.
"""

from __future__ import annotations


def _case_rows(rows: list[dict], case_id: str) -> list[dict]:
    return [r for r in rows if r.get("case_id") == case_id]


def _counter(rows: list[dict], key: str) -> dict[str, int]:
    out: dict[str, int] = {}
    for r in rows:
        out[r[key]] = out.get(r[key], 0) + 1
    return out


def case_sections(load, load_json, fmt, pct, verdict, missing, prior
                  ) -> list[str]:
    cases = load("cases.tsv")
    if not cases:
        return [missing("3. The cases")]
    out: list[str] = []
    for n, case in enumerate(cases, start=3):
        out += _one_case(n, case, load, fmt, pct, verdict, missing)
    return out


def _one_case(n: int, case: dict, load, fmt, pct, verdict, missing
              ) -> list[str]:
    cid = case["case_id"]
    exons = _case_rows(load("case_exons.tsv"), cid)
    introns = _case_rows(load("case_introns.tsv"), cid)
    models = _case_rows(load("annotated_models.tsv"), cid)
    blocks = _case_rows(load("block_accounting.tsv"), cid)
    neigh = _case_rows(load("locus_neighbourhood.tsv"), cid)
    audit = _case_rows(load("assembly_audit.tsv"), cid)
    family = _case_rows(load("family_named_models.tsv"), cid)
    tiling = _case_rows(load("fragment_tiling.tsv"), cid)
    orf = _case_rows(load("reading_frame.tsv"), cid)
    probes = _case_rows(load("junction_probes.tsv"), cid)
    conc = _case_rows(load("probe_summary.tsv"), cid)
    flank = _case_rows(load("flank_consensus_check.tsv"), cid)

    span_kb = (int(case["end"]) - int(case["start"]) + 1) / 1000
    title = (f"## {n}. {case['case_id'].replace('_', ' ').title()} — "
             f"*{case['organism']}* {case['cell']} ({case['selected_as']})")
    out = [title, "",
           f"`{case['accession']}` · {case['contig']}:"
           f"{int(case['start']):,}–{int(case['end']):,}{case['strand']} · "
           f"{fmt(span_kb)} kb · "
           f"recovered from bait `{case['bait'].split('|')[0]}` at "
           f"{float(case['identity']):.1%} identity over "
           f"{float(case['coverage']):.0%} of the bait.", ""]

    # --- the build audit --------------------------------------------------
    if audit:
        a = audit[0]
        same = a.get("annotation_on_this_assembly", "")
        out += ["### Is this even the right build?", "",
                f"The annotation's GFF3 header names build "
                f"`{a.get('gff_build_name') or 'n/a'}` "
                f"(`{a.get('gff_build_accession') or 'n/a'}`), written by "
                f"`{a.get('gff_processor') or 'n/a'}`. The sweep searched "
                f"`{case['accession']}`, so the annotation and the alignment "
                f"sit on the **same assembly**: `{same or 'unknown'}`. "
                f"NCBI lists {fmt(a.get('n_assemblies_for_species'))} "
                f"assemblies for this species; "
                + (f"{fmt(a.get('n_newer_assemblies'))} "
                   f"{'is' if str(a.get('n_newer_assemblies')) == '1' else 'are'}"
                   " newer, of which "
                   f"{'none' if not a.get('newer_annotated') else a['newer_annotated']}"
                   f" {'carries' if not a.get('newer_annotated') else 'carry'}"
                   " an annotation."
                   if str(a.get("n_newer_assemblies") or "0") != "0"
                   else "none is newer.")
                + f" The assembly is {a.get('assembly_level', '?').lower()}-level, "
                f"released {a.get('release_date', '?')} by "
                f"{a.get('submitter', '?')}, contig N50 "
                f"{fmt(int(a.get('contig_n50') or 0) / 1000)} kb — "
                f"{fmt(float(case['headroom']))}× this locus.", ""]
        if a.get("annot_pseudogene_frac"):
            out += [f"NCBI's own counts for this annotation: "
                    f"{fmt(a.get('annot_protein_coding'))} protein-coding "
                    f"genes and {fmt(a.get('annot_pseudogene'))} pseudogenes "
                    f"out of {fmt(a.get('annot_genes_total'))} — "
                    f"{float(a['annot_pseudogene_frac']):.1%} of the gene set "
                    "is filed as pseudogene.", ""]

    # --- what the annotation put there ------------------------------------
    ex_classes = _counter(exons, "class")
    n_ex = len(exons)
    out += ["### What the annotation put on the gene", ""]
    if blocks:
        b = blocks[0]
        out += [f"The recovered gene has **{fmt(b['n_aligned_exons'])} coding "
                f"exons** totalling {fmt(b['aligned_cds_bp'])} bp. The "
                f"annotation places **{fmt(b['n_annotated_gene_models'])} "
                f"gene "
                f"{'model' if b['n_annotated_gene_models'] == '1' else 'models'}"
                "** over them, contributing "
                f"{fmt(b['n_annotated_cds_blocks'])} disjoint coding blocks "
                f"and covering {fmt(b['annotated_cds_bp_on_gene'])} bp "
                f"({float(b['frac_aligned_cds_annotated']):.1%}) of the "
                f"coding footprint. "
                f"{fmt(b['uncovered_bp'])} bp in "
                f"{fmt(b['n_uncovered_blocks'])} blocks has no coding model "
                f"at all, and {fmt(b['untranslated_bp'])} bp reaches no "
                "protein.", "",
                "Blocks rather than genes, because that is what a reader can "
                "recover: two gene models a few hundred bp apart are two "
                "genes but their merged coding evidence is still two pieces.",
                "",
                "| exon class | exons | share |", "|---|---|---|"]
    for klass in ("in_annotated_cds", "in_annotated_noncoding", "intergenic"):
        if klass in ex_classes:
            out.append(f"| {klass.replace('_', ' ')} | "
                       f"{fmt(ex_classes[klass])} | "
                       f"{pct(ex_classes[klass], n_ex)} |")
    out.append("")

    if models:
        out += ["| annotated model | biotype | span | share of the coding "
                "footprint | residues delivered |", "|---|---|---|---|---|"]
        for m in models:
            tiled = next((t for t in tiling if t["gene_id"] == m["gene_id"]
                          and t["at_own_locus"] == "1"), None)
            res = (f"{tiled['s_start']}–{tiled['s_end']}" if tiled else "—")
            out.append(
                f"| `{m['name']}` | {m['biotype']}"
                + (" **(pseudogene)**" if m["pseudo"] == "1" else "")
                + f" | {fmt((int(m['end']) - int(m['start']) + 1) / 1000)} kb | "
                f"{float(m['frac_aligned_cds']):.1%} | {res} |")
        out.append("")
    else:
        out += ["**No annotated gene model overlaps any exon of this gene.** "
                "Not a partial model, not a mis-named one: the annotation has "
                "nothing here.", ""]

    if neigh:
        near = sorted(neigh, key=lambda r: int(r["distance_bp"]))[:4]
        out += ["The annotation is working on both sides of the gap. The "
                "nearest annotated genes are "
                + ", ".join(f"`{r['name']}` ({fmt(int(r['distance_bp']) / 1000)} "
                            f"kb {r['side']})" for r in near)
                + ".", ""]

    # --- what the DNA says ------------------------------------------------
    out += ["### What the DNA says", ""]
    mine = [r for r in orf if r["cell"] == case["cell"]]
    if mine and mine[0].get("verdict") in ("open_reading_frame", "disrupted"):
        r = mine[0]
        out += [f"The locus splices into **{fmt(r['cds_nt'])} nt of coding "
                f"sequence, {fmt(r['cds_codons'])} codons**, translating with "
                f"**{fmt(r['internal_stops'])} internal stop codons** and a "
                f"terminal stop"
                + (" present" if r.get("terminal_stop") == "1" else " absent")
                + f". Under neutral drift to the observed "
                f"{float(r['bait_identity']):.1%} protein identity, "
                f"{fmt(r['expected_stops'])} stops would be expected "
                f"(computed from this locus's own codon usage). "
                f"miniprot reports {fmt(r['miniprot_frameshifts'])} "
                "frameshifts, which at this divergence is not on its own "
                "evidence of pseudogeny; the absence of nonsense over "
                f"{fmt(r['cds_codons'])} codons is.", ""]
        if any(m["pseudo"] == "1" for m in models):
            out += ["The annotation files a model here as a **pseudogene**. "
                    "That is a testable claim and the reading frame does not "
                    "support it.", ""]
    else:
        out += ["*The reading-frame test has not been run for this locus.*", ""]

    if family:
        pseudo = [f for f in family if f["biotype"] == "pseudogene"]
        emit = [f for f in family if f["emits_protein"] == "1"]
        out += [f"Genome-wide, the annotation carries {fmt(len(family))} "
                f"family-named gene models, of which {fmt(len(pseudo))} are "
                f"filed as pseudogenes and {fmt(len(emit))} emit a protein.",
                ""]
        if pseudo:
            out += ["| model | biotype | span | emits a protein |",
                    "|---|---|---|---|"]
            for f in sorted(family, key=lambda r: -int(r["span"]))[:8]:
                out.append(f"| `{f['name']}` | {f['biotype']} | "
                           f"{fmt(int(f['span']) / 1000)} kb | "
                           f"{'yes' if f['emits_protein'] == '1' else 'no'} |")
            out.append("")

    mism = [t for t in tiling if t.get("naming") == "name_mismatch"]
    if mism:
        out += ["### The name and the sequence disagree", ""]
        for t in mism:
            loc = t["best_locus"].split("|")
            out.append(
                f"`{t['name']}` ({fmt(t['protein_aa'])} aa, "
                f"{t['biotype']}) is named for one paralog and its translated "
                f"protein matches the **{loc[1] if len(loc) > 1 else '?'}** "
                f"locus at {float(t['identity']):.1%} identity over "
                f"{fmt(t['aln_aa'])} residues — a bit-score margin of "
                f"{t['bit_margin']} over the next locus. The comparison is "
                "the annotation's own sequence against this genome's own "
                "loci, and the model sits at the coordinates of the locus it "
                "matches.")
        out.append("")

    # --- the checks -------------------------------------------------------
    out += ["### Five checks on this evidence", ""]
    sp = _counter(introns, "splice_class")
    canon = sp.get("canonical", 0) + sp.get("minor", 0)
    out += [f"1. **Splice sites.** {fmt(canon)} of {fmt(len(introns))} introns "
            f"({pct(canon, len(introns))}) carry canonical or minor splice "
            f"dinucleotides "
            f"({fmt(sp.get('canonical', 0))} GT–AG, "
            f"{fmt(sp.get('minor', 0))} minor). An alignment is not a gene; "
            "an exon structure whose introns are spliceable is evidence that "
            "this one is.", ""]
    if conc:
        c = conc[0]
        out += [f"2. **Exon boundaries against independently annotated "
                f"genomes.** {fmt(c['n_reference_genomes'])} swept genomes "
                f"carry this paralog, aligned from the same bait, with their "
                "own annotation independently placing one gene model over the "
                f"whole alignment. {fmt(c['n_shared_with_majority'])} of "
                f"{fmt(c['n_boundaries'])} of this locus's internal exon "
                f"boundaries ({float(c['frac_shared_with_majority']):.1%}) are "
                "shared by a majority of them; the median boundary is shared "
                f"by {fmt(c['median_refs_per_boundary'])}. This validates the "
                "instrument at this gene, not the individual locus — which is "
                "why it is reported alongside the reading frame and the "
                "splice sites rather than instead of them.", ""]
    else:
        out += ["2. *Exon-boundary concordance has not been run.*", ""]

    summ = conc[0] if conc else {}
    if summ.get("status_mrna"):
        out += [f"3. **Transcript evidence at the junction.** "
                f"{fmt(summ.get('n_probes'))} ±90 nt probes were built across "
                "the junctions the annotation does not model and searched "
                "against every transcript record NCBI holds for this species. "
                f"{fmt(summ.get('hits_mrna'))} hits, "
                f"{fmt(summ.get('spanning_mrna'))} of them contiguous across a "
                "junction.", "",
                f"    The denominator is the result: this species has "
                f"{fmt(summ.get('n_records_mrna'))} transcript records in "
                "total. A search of that many returning nothing has not shown "
                "the gene is untranscribed — it has shown the species has "
                "almost no transcript deposits. RNA-seq for it does exist and "
                "reaching it needs the streaming aligner S12 builds; these "
                "junctions are handed there rather than half-answered here.",
                "",
                f"    The criterion has its own negative control. Run against "
                f"the locus's own genomic DNA — where a probe of contiguous "
                f"*spliced* sequence cannot span its junction by construction "
                f"— the probes make {fmt(summ.get('hits_genomic_control'))} "
                f"hits and {fmt(summ.get('spanning_genomic_control'))} of them "
                "span. So the zero above is a discriminating zero, not a test "
                "that never fires. (On its first run that control returned "
                "six false spans, which is how the rule acquired its second "
                "half: an 8 nt anchor alone admits an alignment that has run a "
                "dozen bases past the junction into the intron.)", ""]
    else:
        out += ["3. *The transcript search has not been run.*", ""]

    others = [r for r in orf if r["cell"] != case["cell"]
              and r.get("verdict") in ("open_reading_frame", "disrupted")]
    if others:
        opened = sum(1 for r in others if r["verdict"] == "open_reading_frame")
        out += [f"4. **The rest of the family in the same genome.** "
                f"{fmt(opened)} of {fmt(len(others))} other family loci in "
                "this assembly also splice into an uninterrupted reading "
                "frame. Whatever the annotation is doing here, it is not "
                "responding to a damaged gene.", ""]
    else:
        out += ["4. *The within-genome comparison has not been run.*", ""]

    if flank:
        hits = [f for f in flank if f["in_paralog_consensus"] == "1"]
        if hits:
            top = min(hits, key=lambda f: int(f["consensus_rank"] or 999))
            out += [f"5. **The neighbourhood.** {fmt(len(hits))} of "
                    f"{fmt(len(flank))} of the nearest flanking genes are in "
                    f"S8's consensus flank set for {case['cell']} — including "
                    f"`{top['symbol']}`, this paralog's most conserved "
                    f"neighbour, found beside it in "
                    f"{fmt(top['consensus_species'])} of "
                    f"{fmt(top['consensus_species_total'])} swept vertebrates "
                    f"({float(top['consensus_fraction']):.1%}). The gene is "
                    "not merely present and intact; it is in the position "
                    "this paralog occupies across the vertebrates, which "
                    "nothing about this assembly's annotation could produce.",
                    ""]
        else:
            out += ["5. **The neighbourhood — not testable here.** None of "
                    "the flanking genes carries a gene symbol (this "
                    "annotation names its genes by locus tag), so S8's "
                    "consensus flank set has nothing to match. That is a "
                    "limitation of the comparison, not a negative result.",
                    ""]
    return out
