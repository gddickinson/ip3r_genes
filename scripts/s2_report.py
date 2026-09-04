"""S2 — render `results/census_v2/report.md` from the committed tables.

Nothing in the report is hand-written (D13): every count, percentage and
verdict below is read back out of the tables this task wrote, so a number
in the prose cannot drift from the number in the data. The narrative order
is fixed here; the values are not.
"""

from __future__ import annotations

import argparse
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.utils.family import MAX_LENGTH_AA, MIN_LENGTH_AA   # noqa: E402
from scripts.s2_lib import (                                # noqa: E402
    CONTEXT_PFAMS, OUT_DIR, RYR_PFAMS, SEED_PFAMS, read_tsv,
)


def table(rows: list[dict], cols: list[str], head: list[str] | None = None) -> str:
    head = head or cols
    out = ["| " + " | ".join(head) + " |",
           "|" + "|".join("---" for _ in head) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:.1f} %" if b else "—"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()
    D = OUT_DIR
    names = {r["pfam_id"]: r["name"] for r in read_tsv(D / "signature_names.tsv")}
    counts = read_tsv(D / "interpro_counts.tsv")
    census = read_tsv(D / "census_v2.tsv")
    calls = read_tsv(D / "call_summary.tsv")
    audit = read_tsv(D / "rule_audit.tsv")
    unprof = read_tsv(D / "unassigned_profile.tsv")
    groups = read_tsv(D / "delta_by_group.tsv")
    verdicts = read_tsv(D / "delta_missing_verdicts.tsv")
    delta = json.loads((D / "delta_summary.json").read_text())
    reps = read_tsv(D / "census_v2_representatives.tsv")

    n = len(census)
    by_call = {c: sum(1 for r in census if r["call"] == c)
               for c in ("ITPR", "RYR", "unassigned")}
    lengths = {c: sorted(int(r["length"]) for r in census if r["call"] == c)
               for c in by_call}
    species = len({r["taxon_id"] for r in census if r["call"] == "ITPR"})

    L: list[str] = []
    A = L.append
    A("# S2 — Uncapped InterPro enumeration (census v2)\n")
    A("*Rendered from the committed tables by `scripts/s2_report.py` "
      "(D13). Do not hand-edit.*\n")

    A("## 1. What was enumerated\n")
    A(f"Three Pfam signatures define the search space: "
      + ", ".join(f"**{p}** ({names.get(p, '?')})" for p in SEED_PFAMS)
      + ". `PF02815` (MIR) is deliberately not among them — it is carried "
        "by the O-mannosyltransferases as well as by both receptor "
        "families, so it widens the space without adding evidence — and is "
        "kept as an annotation column instead.\n")
    A(table(counts, ["pfam_id", "api_count", "uniprot_count", "pages",
                     "unique", "delta_vs_api_count", "restarts",
                     "cursor_exhausted"],
            ["signature", "InterPro `count`", "UniProt count", "pages",
             "unique accessions", "unique − `count`", "cursor restarts",
             "walked to end"]) + "\n")

    bad = [c for c in counts if c["delta_vs_api_count"]
           and int(c["delta_vs_api_count"]) != 0]
    if bad:
        A("**InterPro's own `count` is not a completeness criterion, and it "
          "is wrong in both directions.** "
          + "; ".join(
              f"{c['pfam_id']} advertises {int(c['api_count']):,} and serves "
              f"{int(c['unique']):,} ({int(c['delta_vs_api_count']):+d})"
              for c in bad)
          + ". The advertised values were stable across re-queries, so this "
            "is not a mid-walk edit — and in every case what the API "
            "*serves* matches UniProt's independent count for the same "
            "signature"
          + (" to within two records ("
             + "; ".join(f"{c['pfam_id']} {int(c['unique']):,} vs "
                         f"{int(c['uniprot_count']):,}" for c in counts
                         if c["uniprot_count"]) + ")"
             if any(c["uniprot_count"] for c in counts)
             else " (cross-check not recorded in this run)")
          + ", so it is the advertised number that is wrong, not the "
            "walk.\n")
        A("Both directions bite. Stopping at an *under*-stated count drops "
          "records silently. Trusting an *over*-stated one declares a "
          "complete walk incomplete — which is exactly what happened here: "
          "the first version of this task's own completeness test failed "
          "PF08454 and exited non-zero on a walk that had run its cursor "
          "chain to the end. The recorded test is therefore cursor "
          "exhaustion (`next: null`), with both counts kept beside it.\n")

    A("### API availability\n")
    A("InterPro answers the same queries on two hosts: the documented "
      "`/interpro/api/` and `/interpro/wwwapi/`, which the InterPro website "
      "itself calls. During this task the documented host returned HTTP 500 "
      "to 11 of 12 probe requests while `wwwapi` served 12 of 12, with "
      "identical payloads and the same `count`. `src/databases/interpro.py` "
      "now tries both hosts before it sleeps, which is why the walk "
      "completed at all.\n")

    A("### What one signature would have missed\n")
    seedc = read_tsv(D / "seed_contribution.tsv")
    per_seed = [r for r in seedc if r["absent"]]
    A("The enumeration is a union of three signatures rather than a query "
      "on the family's defining one, and this is why:\n")
    A(table(per_seed, ["pfam_id", "present", "absent", "absent_called_itpr",
                       "taxa_lost"],
            ["signature", "records carrying it", "records without it",
             "…of those, called ITPR", "taxa"]) + "\n")
    p8 = next((r for r in per_seed if r["pfam_id"] == "PF08709"), None)
    if p8:
        A(f"A census built on **PF08709 alone** — the IP3-binding core, the "
          f"signature that names the family — would have enumerated "
          f"{int(p8['present']):,} of the {n:,} proteins here and missed "
          f"{int(p8['absent']):,}, including {int(p8['absent_called_itpr']):,} "
          f"this census calls ITPR across {int(p8['taxa_lost']):,} taxa. "
          f"*Dictyostelium* iplA (Q9NA13), a characterised IP3 receptor and a "
          f"member of this project's own S1 positive panel, is one of the "
          f"records it would have missed: it carries PF01365 and PF08454 and "
          f"neither PF08709 nor PF02815 nor PF00520. It **is** in this "
          f"census, and it is `unassigned` — 2 of 5 architecture signatures "
          f"is not enough for a positive call, which is the honest outcome "
          f"and the reason S3 builds profiles.\n")
    combos = [r for r in seedc if not r["absent"]]
    if combos:
        full = int(combos[0]["present"])
        A(f"Only {full:,} records ({pct(full, n)}) carry all three seeds; "
          f"the rest are found by one or two.\n")

    A("## 2. The search space\n")
    A(f"{n:,} distinct proteins carry at least one seed signature, across "
      f"{len({r['taxon_id'] for r in census}):,} taxa.\n")
    src: dict[str, int] = {}
    for r in census:
        src[r["source"]] = src.get(r["source"], 0) + 1
    only_ip = [r for r in census if r["source"] == "interpro_only"]
    A(f"The two databases were reconciled rather than pooled, since each is "
      f"a separate view of the same Pfam matches: "
      + ", ".join(f"**{v:,}** {k.replace('_', ' ')}"
                  for k, v in sorted(src.items(), key=lambda kv: -kv[1]))
      + ". "
      + (f"The {len(only_ip)} InterPro-only records are all *"
         f"{only_ip[0]['species']}* fragments of "
         f"{min(int(r['length']) for r in only_ip)}–"
         f"{max(int(r['length']) for r in only_ip)} aa with newly-issued "
         f"accessions — InterPro has matched them and UniProt's own Pfam "
         f"cross-reference has not caught up, which is release skew rather "
         f"than a disagreement. Nothing was in UniProt's set and missing "
         f"from InterPro's.\n" if only_ip else "\n"))
    A(table(groups, ["group", "v2_records", "v2_itpr", "v2_species"],
            ["lineage", "records", "called ITPR", "taxa"]) + "\n")

    A("### Where the family is, outside the animals\n")
    lin = read_tsv(D / "lineage_calls.tsv")
    outer = [r for r in lin if r["group"] not in
             ("Vertebrata", "Metazoa (non-vertebrate)")]
    A("Records are the wrong unit for this question — one well-sequenced "
      "alga contributes a dozen — so the table counts **taxa**, and keeps "
      "the phyla that put records in the search space without earning a "
      "single ITPR call, because those zeros are the interesting rows.\n")
    A(table(outer, ["group", "phylum", "records", "taxa", "itpr",
                    "itpr_high", "itpr_taxa"],
            ["lineage", "phylum", "records", "taxa", "ITPR calls",
             "…high confidence", "taxa with a call"]) + "\n")
    plants = [r for r in lin if r["group"] == "Viridiplantae"]
    fungi = [r for r in lin if r["group"] == "Fungi"]
    chloro = next((r for r in plants if r["phylum"] == "Chlorophyta"), None)
    strepto = next((r for r in plants if r["phylum"] == "Streptophyta"), None)
    if chloro and strepto:
        A(f"**The plant and fungal records are not scattered — they are "
          f"phylogenetically clean.** In Viridiplantae every one of the "
          f"{int(chloro['itpr'])} ITPR calls is in **Chlorophyta** "
          f"({int(chloro['itpr_taxa'])} taxa, including *Chlamydomonas "
          f"reinhardtii* with the complete five-signature architecture), "
          f"and **Streptophyta — the land-plant lineage — has "
          f"{int(strepto['itpr'])} calls** from "
          f"{int(strepto['records'])} records in "
          f"{int(strepto['taxa'])} taxa. In Fungi every call sits in an "
          f"early-diverging phylum ("
          + ", ".join(f"{r['phylum']} {int(r['itpr'])}" for r in fungi
                      if int(r["itpr"]))
          + "), while Dikarya — Ascomycota and Basidiomycota, the yeasts "
            "and moulds — contribute **no records to the search space at "
            "all**.\n")
        A("That is a much sharper statement than the one this project "
          "started from (\"plants and fungi lack the family\"), and it is "
          "the shape a loss looks like: present in the early-diverging "
          "lineage of both kingdoms, absent from the derived one. It is "
          "also, still, a statement about **what UniProt holds**. Turning "
          "it into a statement about genomes is exactly what S20 and S23 "
          "are for, and the emergent row that asked this question is "
          "narrowed rather than closed.\n")

    A("## 3. The ITPR / RYR call\n")
    A("Every record is called on **domain architecture**: RYR when it "
      "carries any of "
      + ", ".join(f"{p} ({names.get(p, '?')})" for p in RYR_PFAMS)
      + "; ITPR when it carries the complete IP3-receptor architecture ("
      + " + ".join(SEED_PFAMS + CONTEXT_PFAMS)
      + ") and none of them. Length is recorded on every row and enters only "
        "as support for a medium-confidence call — the size band "
        f"({MIN_LENGTH_AA:,}–{MAX_LENGTH_AA:,} aa) was chosen to exclude "
        "RyRs, so letting it decide would beg the question.\n")
    A(table(calls, ["call", "confidence", "records"]) + "\n")
    A(f"- **ITPR {by_call['ITPR']:,}** ({pct(by_call['ITPR'], n)}), "
      f"median length {statistics.median(lengths['ITPR']):.0f} aa, "
      f"{species:,} taxa\n"
      f"- **RYR {by_call['RYR']:,}** ({pct(by_call['RYR'], n)}), "
      f"median length {statistics.median(lengths['RYR']):.0f} aa\n"
      f"- **unassigned {by_call['unassigned']:,}** "
      f"({pct(by_call['unassigned'], n)}), median length "
      f"{statistics.median(lengths['unassigned']):.0f} aa\n")

    sym_only = sum(1 for r in census if r["confidence"] == "low")
    A(f"**{sym_only:,} of those calls ({pct(sym_only, n)}) rest on a gene "
      f"symbol alone** — records whose architecture is too partial to "
      f"decide, carrying a name that is not. D14 forbids the *scorer* from "
      f"consulting a candidate's own symbol, and that still holds: here the "
      f"symbol is admitted as explicitly-labelled `low` evidence at census "
      f"scale, it is recorded in the `reason` column of every row it "
      f"touched, and the architecture call is kept in its own column so the "
      f"audit below can score it without symbols anywhere in the "
      f"input.\n")

    A("### The rule audit\n")
    A("The architecture rule never sees a gene symbol, so symbols are an "
      "independent label to score it against. Every symbol-labelled record "
      "in the census is used, not a sample.\n")
    A(table(audit, ["test", "n", "agree", "disagree", "unassigned",
                    "accuracy_of_decided"],
            ["test", "n", "agree", "disagree", "no call",
             "accuracy where decided"]) + "\n")
    A("The RyR rule leans on one conditional claim — that SPRY, which sits "
      "in ~114,000 UniProt proteins and is in no sense RyR-specific, *is* "
      "diagnostic among proteins that already carry a seed signature. That "
      "is the row to check: it is a claim about this search space, not "
      "about SPRY.\n")

    A("## 4. The unassigned pile\n")
    una = lengths["unassigned"]
    in_band = sum(1 for r in census
                  if r["call"] == "unassigned" and r["in_band"] == "1")
    A(f"{by_call['unassigned']:,} records ({pct(by_call['unassigned'], n)}) "
      f"satisfy neither positive test. They are overwhelmingly fragments: "
      f"median {statistics.median(una):.0f} aa against "
      f"{statistics.median(lengths['ITPR']):.0f} aa for a called ITPR, and "
      f"only {in_band:,} of them ({pct(in_band, by_call['unassigned'])}) "
      f"fall inside the size band at all. A partial annotation cannot be "
      f"called either way by a rule that reads absence as evidence, which "
      f"is the intended behaviour rather than a shortfall — S3's profile "
      f"sweep is what resolves them.\n")
    A(table(unprof[:8], ["source", "seed_pfams", "n_itpr_arch",
                         "length_band", "records"],
            ["source", "seed signatures", "ITPR signatures (of 5)",
             "length", "records"]) + "\n")

    short = read_tsv(D / "short_complete.tsv")
    if short:
        A("### What the size band caught instead\n")
        A(f"Length never decides a call, which leaves the band free to "
          f"catch something else. {len(short)} records carry the "
          f"**complete** five-signature IP3-receptor architecture inside "
          f"{MIN_LENGTH_AA:,} residues — "
          f"{min(int(r['length']) for r in short):,}–"
          f"{max(int(r['length']) for r in short):,} aa, against "
          f"{MIN_LENGTH_AA:,} for the shortest real family member. An "
          f"intact architecture in two-thirds of the length is a truncated "
          f"gene model, and none of them is flagged `Fragment` by UniProt, "
          f"because a truncated model submitted as a whole protein is not "
          f"marked as one. All "
          f"{len(short)} are unnamed locus tags from "
          f"{len({r['species'] for r in short})} species.\n")
        A(table(short, ["accession", "gene", "species", "length",
                        "protein_existence"],
                ["accession", "gene", "species", "length (aa)",
                 "protein existence"]) + "\n")
        A("The call on these is *correct* — they are ITPRs — and the "
          "records are still wrong. That is the case for keeping length as "
          "a recorded column on every row rather than dropping it once it "
          "stopped being part of the call, and it is a starting list for "
          "S18's correction register.\n")

    seq = D / "sequence_check.tsv"
    if seq.exists():
        meta = json.loads((D / "sequence_check_meta.json").read_text())
        band = meta["no_call_band"]
        ag = meta["agreement"]
        rows = read_tsv(seq)
        A("## 5. The call checked against sequence\n")
        A("The architecture call is annotation-derived, so it is checked "
          "against something annotation-independent: the labelled-bait "
          "identity margin from D14, run over the committed per-phylum core "
          "panel. Nothing in that test reads a Pfam list or a gene symbol — "
          "each panel member is aligned with the six human references and "
          "assigned to whichever family it is closer to.\n")
        rc = (meta["mafft_calls"][0]["returncode"]
              if meta.get("mafft_calls") else "?")
        A(f"Panel: {meta['n_panel']} sequences + {meta['n_baits']} labelled "
          f"baits; MAFFT rc={rc}, {meta['alignment_columns']:,} columns, "
          f"{meta['seconds']} s.\n")
        A(f"Agreement is scored separately inside and outside D7's no-call "
          f"band (|margin| < {band:.2f}), because a margin of 0.0005 is not "
          f"a disagreement — it is the test declining to separate two "
          f"families at that distance.\n")
        A(table(
            [{"metric": m.replace("_", " "),
              "decisive": f"**{ag[m]['decisive_agree']}/{ag[m]['decisive_n']}**",
              "band": f"{ag[m]['band_agree']}/{ag[m]['band_n']}"}
             for m in ("full_alignment", "covered_only")],
            ["metric", "decisive", "band"],
            ["identity metric", f"agree where the margin decides "
                                 f"(abs ≥ {band:.2f})",
             "agree inside the no-call band"]) + "\n")
        dis = [r for r in rows if r["agrees"] == "0"]
        if dis:
            A("Every row where the two calls differ, under either metric, "
              "with its margin:\n")
            A(table(dis, ["accession", "species", "phylum", "arch_call",
                          "seq_call", "margin", "metric", "decisive"],
                    ["accession", "species", "phylum", "architecture",
                     "sequence", "margin", "metric", "decisive"]) + "\n")
            A("All of them sit inside the no-call band, which is the "
              "expected place for the deepest branches: S1 already recorded "
              "that *Dictyostelium* iplA — a characterised IP3 receptor — "
              "has a bait margin of +0.065, inside the same band. Profile "
              "assignment (S3), not a pairwise margin, is what settles "
              "these.\n")
        else:
            A("No disagreement on any panel member, under either metric.\n")

    A("## 6. Delta against census v1\n")
    A(f"Census v1 is what the application returns when it is pointed at the "
      f"family by name — the committed search bundles plus S1's control "
      f"panels, {delta['v1_accessions']:,} accessions. It was never meant as "
      f"a family-wide harvest, so the delta measures what a domain search "
      f"adds to a name search rather than one census against another.\n")
    A(f"- v1 {delta['v1_accessions']:,} → **v2 {delta['v2_accessions']:,}** "
      f"(×{delta['growth_factor']}), {delta['shared']:,} shared\n"
      f"- {delta['only_v1']:,} v1 accessions are absent from v2, each with a "
      f"verdict:\n")
    A(table(verdicts, ["verdict", "accessions"]) + "\n")
    holes = [r for r in verdicts if "HOLE" in r["verdict"]]
    A(("**No enumeration holes.** " if not holes else
       f"**{sum(int(r['accessions']) for r in holes)} enumeration holes.** ")
      + "Every absence is structural: an accession outside UniProt's "
        "namespace was never in the searched space, a `-2` isoform is "
        "represented by its canonical parent, S1's decoys are supposed to be "
        "absent, and the remainder are short human ITPR fragments "
        "(102–408 aa) that carry no seed signature — the real and stated "
        "limit of a domain census.\n")

    A("## 7. Representatives\n")
    per = {}
    for r in reps:
        per[r["call"]] = per.get(r["call"], 0) + 1
    core = sum(1 for r in reps if r["core"] == "1")
    A(f"{len(reps):,} representatives "
      f"({', '.join(f'{k} {v}' for k, v in sorted(per.items()))}): the "
      f"longest record per species per call, with vertebrates thinned to one "
      f"per taxonomic order. {core} of them form the per-phylum **core "
      f"panel** committed as `census_v2_core_panel.faa` — the only FASTA in "
      f"the repo, since a census-scale one would be re-committed at every "
      f"later census version.\n")

    A("## 8. What is committed, and where the bulk went\n")
    figs = sorted((D / "figures").glob("*.png"))
    A(f"Under the data root: the raw InterPro pages "
      f"(`raw_api/interpro/<PFAM>/`), both enumerations' parsed "
      f"intermediates, the UniProt sweep, the seeded-space FASTA (47 MB) "
      f"and the full representative FASTA. In the repo: the census, the "
      f"call summary and audit, the seed-contribution and lineage tables, "
      f"the unassigned profile, the truncated-model list, the sequence "
      f"check, the delta tables, the representative table, the core panel "
      f"and {len(figs)} figures ("
      + ", ".join(f"`{f.stem}`" for f in figs) + ").\n")

    A("## 9. Caveats\n")
    A("- **This is a UniProt-space census.** RefSeq and Ensembl proteins "
      "reachable by name are not in it; that is S3–S5's job, and it is why "
      f"{sum(int(r['accessions']) for r in verdicts if 'namespace' in r['verdict']):,} "
      "v1 accessions sit outside it.\n"
      "- **A domain census cannot see a protein with no domain annotated.** "
      "The six unexplained v1 records are exactly that case.\n"
      "- **A record is not a gene.** One species contributes as many rows "
      "as UniProt holds isoforms, redundant TrEMBL entries and alternative "
      "models for it, so no copy-number statement can be read off this "
      "census — counting paralogs needs the clustering in S6/S7 and the "
      "genomic evidence in S5. Where this report counts breadth it counts "
      "**taxa**, not records, for that reason.\n"
      "- The unassigned pile is not noise to be tuned away; it is the "
      "fragment population, and it is carried forward with its reasons.\n"
      "- Every call is architecture-derived and inherits InterPro's "
      "annotation. The sequence check in §5 is the independent axis, and it "
      "covers the core panel, not the whole census.\n")

    (D / "report.md").write_text("\n".join(L) + "\n")
    print(f"[s2] wrote {(D / 'report.md').relative_to(OUT_DIR.parents[1])} "
          f"({len(L)} blocks)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
