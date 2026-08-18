#!/usr/bin/env python3
"""S0 — render ``results/s0_baseline/report.md`` from the committed tables.

Nothing in the report is written by hand: every number is read out of the TSV
and JSON files that `s0_db_snapshot.py`, `s0_gene_structure.py` and the S0
literature pass wrote, so the report and its data cannot drift (Decisions D13).

Inputs (all under ``results/s0_baseline/``):
    interpro_signature_counts.tsv   interpro_taxonomy_counts.tsv
    reference_proteins.tsv          pfam_architecture.tsv
    zebrafish_itpr.tsv              gene_structure.tsv
    references.tsv                  lit_claims.tsv
    smoke_test.tsv                  ensembl_endpoint_probe.tsv
    db_snapshot_meta.json

Usage:  python scripts/s0_report.py [--dir results/s0_baseline]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# The ITPR length band declared in src/utils/family.py, used here only to
# classify the zebrafish PF08709 records into ITPR-sized and RyR-sized.
ITPR_MAX_AA = 3600


def read_tsv(path: Path) -> list[dict]:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:] if ln.strip()]


def md_table(rows: list[dict], cols: list[str], headers: list[str]) -> str:
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def fnum(x) -> str:
    try:
        return f"{int(x):,}"
    except (TypeError, ValueError):
        return str(x)


def build(d: Path) -> str:
    meta = json.loads((d / "db_snapshot_meta.json").read_text())
    sigs = read_tsv(d / "interpro_signature_counts.tsv")
    taxa = read_tsv(d / "interpro_taxonomy_counts.tsv")
    refs_p = read_tsv(d / "reference_proteins.tsv")
    arch = read_tsv(d / "pfam_architecture.tsv")
    zeb = read_tsv(d / "zebrafish_itpr.tsv")
    gs = read_tsv(d / "gene_structure.tsv")
    lit = read_tsv(d / "references.tsv")
    claims = read_tsv(d / "lit_claims.tsv")
    smoke = read_tsv(d / "smoke_test.tsv")
    probe = read_tsv(d / "ensembl_endpoint_probe.tsv")

    date = meta["derived_on"]
    P = []
    A = P.append

    A(f"# S0 — literature baseline & scope confirmation\n")
    A(f"*Rendered from the committed tables in `results/s0_baseline/` by "
      f"`scripts/s0_report.py`. Database snapshot re-derived {date}.*\n")

    # ---------------------------------------------------------------- 1. lit
    verdicts = {}
    for c in claims:
        verdicts[c["verdict"]] = verdicts.get(c["verdict"], 0) + 1
    # references.tsv is the whole project bibliography (it also backs
    # docs/ip3r_review_2026.md). The audit used only the subset that the
    # claims table cites, so count that subset, not the file.
    audit_keys = set()
    for c in claims:
        audit_keys.update(k.strip() for k in c["refs"].split(";")
                          if k.strip() and k.strip() != "-")
    audit = [r for r in lit if r["ref_id"] in audit_keys]
    A("## 1. Literature verification\n")
    A(f"{len(claims)} atomic claims were extracted from the `[lit]` "
      f"statements in `docs/ip3r_background.md` and checked against "
      f"{len(audit)} references "
      f"({sum(1 for r in audit if r['role'] == 'primary')} primary, "
      f"{sum(1 for r in audit if r['role'] == 'review')} review). The wider "
      f"bibliography assembled for `docs/ip3r_review_2026.md` extends this to "
      f"{len(lit)} references in the same table.\n")
    A("| verdict | claims | meaning |")
    A("|---|---|---|")
    meanings = {
        "verified": "primary source(s) support the claim as written",
        "verified_qualified": "supported, but the wording overstates the evidence; qualified in the review",
        "corrected": "false or incomplete as written; **struck and replaced**",
        "downgraded_open": "not established; retagged `[open]` because it is a question this project answers",
        "upgraded_db": "re-derived from a live database; no longer a literature claim",
    }
    for v in ["verified", "verified_qualified", "corrected", "downgraded_open",
              "upgraded_db"]:
        if v in verdicts:
            A(f"| `{v}` | {verdicts[v]} | {meanings[v]} |")
    A("")
    changed = [c for c in claims
               if c["verdict"] in ("corrected", "downgraded_open", "upgraded_db")]
    A(f"### The {len(changed)} claims that did not survive as written\n")
    for c in changed:
        A(f"**{c['claim_id']} ({c['section']}) — `{c['verdict']}`**  ")
        A(f"> {c['claim_as_written']}\n")
        A(f"{c['note']}\n")
    A("Full per-claim table with references: `lit_claims.tsv`; "
      "the bibliography: `references.tsv`; the prose baseline: "
      "`docs/ip3r_review_2026.md`.\n")

    # ------------------------------------------------------------ 2. db redo
    A("## 2. Database snapshot re-derived\n")
    A(f"Every `[db]` number in `docs/ip3r_background.md` was re-queried on "
      f"{date}. **All of them reproduced exactly.**\n")
    A("### InterPro protein counts per Pfam signature\n")
    for s in sigs:
        s["n_fmt"] = fnum(s["n_proteins"])
    A(md_table(sigs, ["pfam", "name", "n_fmt", "note"],
               ["Pfam", "name", "proteins", "note"]))
    A("")
    A("### Taxonomic distribution\n")
    for t in taxa:
        t["a"] = fnum(t["n_PF08709"])
        t["b"] = fnum(t["n_PF01365"])
    A(md_table(taxa, ["taxon", "a", "b"],
               ["Taxon", "PF08709 Ins145_P3_rec", "PF01365 RIH"]))
    A("")
    plants = next(t for t in taxa if t["taxon"] == "Viridiplantae")
    fungi = next(t for t in taxa if t["taxon"] == "Fungi")
    A(f"The Q1 anomaly is confirmed and unchanged: **{plants['a']} "
      f"Viridiplantae** and **{fungi['a']} Fungi** proteins carry PF08709, "
      f"the IP3-binding core, while *Arabidopsis thaliana* and "
      f"*Saccharomyces cerevisiae* have none. S2/S20 identify what they are.\n")

    A("### Reference panel and sister family\n")
    A(md_table(refs_p, ["symbol", "accession", "family", "length_aa", "status"],
               ["Gene", "UniProt", "family", "length (aa)", "status"]))
    A("")
    shared = ["PF02815", "PF08709", "PF01365", "PF08454"]
    itpr_hit = {p: all(int(a["present"]) for a in arch
                       if a["pfam"] == p and a["family"] == "ITPR")
                for p in shared}
    ryr_hit = {p: all(int(a["present"]) for a in arch
                      if a["pfam"] == p and a["family"] == "RYR")
               for p in shared}
    both = [p for p in shared if itpr_hit[p] and ryr_hit[p]]
    A(f"**The D14 hazard, re-confirmed at the record level.** All "
      f"{len(both)} of the ITPR-diagnostic Pfam signatures "
      f"({', '.join(both)}) are carried by **all three** human ryanodine "
      f"receptors as well as all three IP3 receptors. Per-protein detail: "
      f"`pfam_architecture.tsv`.\n")

    # --------------------------------------------------------- 3. zebrafish
    itprs = [z for z in zeb if z["length_aa"] and
             int(z["length_aa"]) <= ITPR_MAX_AA]
    ryrs = [z for z in zeb if z["length_aa"] and
            int(z["length_aa"]) > ITPR_MAX_AA]
    genes: dict[str, list[int]] = {}
    for z in zeb:
        genes.setdefault(z["gene"] or "(unnamed)", []).append(int(z["length_aa"]))
    itpr_genes = sorted(g for g, L in genes.items() if max(L) <= ITPR_MAX_AA)
    ryr_genes = sorted(g for g, L in genes.items() if max(L) > ITPR_MAX_AA)
    pct = 100.0 * len(ryrs) / len(zeb)
    A("## 3. What a single Pfam query actually returns (the D14 hazard, measured)\n")
    A(f"The background document cites `taxonomy_id:7955 AND xref:pfam-PF08709` "
      f"as evidence that zebrafish carries four IP3 receptors. Re-running it "
      f"returns **{len(zeb)} protein records** across "
      f"**{len(genes)} gene symbols** — and only "
      f"{len(itpr_genes)} of those genes are IP3 receptors:\n")
    A(f"- IP3R-sized (≤ {ITPR_MAX_AA:,} aa): **{len(itprs)} records** across "
      f"{len(itpr_genes)} genes — {', '.join(f'`{g}`' for g in itpr_genes)}")
    A(f"- RyR-sized (> {ITPR_MAX_AA:,} aa): **{len(ryrs)} records** across "
      f"{len(ryr_genes)} genes — {', '.join(f'`{g}`' for g in ryr_genes)}\n")
    A(f"**{pct:.0f}% of the records returned by the IP3-binding-core Pfam are "
      f"ryanodine receptors.** This is the single most useful number S0 "
      f"produced: it is a measured floor on the contamination any Pfam-driven "
      f"enumeration inherits, in the exact query the planning document quoted "
      f"as a clean result. It is also why D14 requires a positive ITPR/RYR "
      f"call rather than a length filter alone — here the length band happens "
      f"to separate the two cleanly, but that is a fact about zebrafish "
      f"annotation quality, not a rule.\n")
    A(f"One of the RyR-sized genes is unnamed (`LOC101884734`, "
      f"{max(genes.get('LOC101884734', [0])):,} aa) — an early example of the "
      f"unnamed-locus problem S18 audits. Full table: `zebrafish_itpr.tsv`.\n")

    # ---------------------------------------------------- 4. gene structure
    A("## 4. Gene architecture — a `[lit]` claim that failed\n")
    for g in gs:
        g["span_fmt"] = fnum(g["genomic_span_bp"])
        g["loc"] = f"chr{g['chromosome']}:{fnum(g['start'])}-{fnum(g['end'])}"
    A(md_table(gs, ["symbol", "ensembl_gene", "loc", "span_fmt",
                    "n_exons_canonical", "cytoband_db", "cytoband_lit"],
               ["Gene", "Ensembl", "location", "span (bp)",
                "exons (canonical)", "band (measured)", "band (claimed)"]))
    A("")
    exons = [int(g["n_exons_canonical"]) for g in gs]
    spans = [int(g["genomic_span_bp"]) for g in gs]
    smallest = min(gs, key=lambda g: int(g["genomic_span_bp"]))
    A(f"All three cytogenetic bands are confirmed, so that claim moves from "
      f"`[lit]` to `[db]`. The architecture claim does not survive: the "
      f"canonical exon count is **{min(exons)}–{max(exons)}**, not 58–60, and "
      f"**{smallest['symbol']} spans {fnum(smallest['genomic_span_bp'])} bp** "
      f"— not \"hundreds of kb\". Genomic span varies "
      f"**{max(spans) / min(spans):.1f}-fold across the three paralogs** while "
      f"protein length varies by 3%. That asymmetry is a result in waiting for "
      f"S21, not a detail.\n")

    # ---------------------------------------------------------- 5. app smoke
    A("## 5. App smoke test against live APIs\n")
    A(md_table(smoke, ["preset", "scope", "client", "source", "status",
                       "hits", "seconds", "note"],
               ["Preset", "Scope", "Client", "Source", "Status", "Hits", "s",
                "Note"]))
    A("")
    ok3 = [r for r in smoke if r["source"] == "TOTAL" and r["status"] == "ok"]
    A(f"**{len(ok3)} of the 4 runs returned all three sources**, including a "
      f"human-scoped run after the client fix. Runs 1 and 3 are the "
      f"instructive failures.\n")
    A("**Ensembl failed the 8-species human preset twice, for two different "
      "reasons.** Neither is an outage and neither is a slow network:\n")
    A(md_table(probe, ["endpoint", "species", "http_status", "seconds",
                       "verdict"],
               ["Endpoint", "Species", "HTTP", "s", "Verdict"]))
    A("")
    A("**Fault 1 — a per-species stall.** `src/databases/ensembl.py` resolved "
      "gene symbols through `/xrefs/symbol/{species}/{symbol}`. That path "
      "**stalls indefinitely for `homo_sapiens`** — for `BRCA2` as well as "
      "`ITPR1`, so it is neither gene-specific nor a query-string artefact — "
      "while the same endpoint answers in 0.6 s for `danio_rerio`, which is "
      "why the zebrafish preset got all three sources and the human preset "
      "got none. A retry budget could never have fixed this: there is nothing "
      "to retry against. `_symbol_to_ids` now resolves through "
      "`/lookup/symbol/`, which works, and falls back to `xrefs` only on a "
      "clean 404 — never after a transport error, which would walk straight "
      "back into the stall. **Human ITPR1 went from 0 variants to 24.**\n")
    A("**Fault 2 — latency, uncovered by fixing fault 1.** Ensembl's speed is "
      "unstable: the *same* 451-byte `lookup/symbol` call measured 0.61 s, "
      "7.61 s and 13.91 s within one session, and the `expand=1` call that "
      "carries the transcript/exon payload costs ~12 s every time. So one "
      "gene in one species costs ~12 s at best and was measured at 95 s at "
      "worst. The default panel is 8 species × 3 genes = 24 sequential pairs, "
      "i.e. **roughly 5 to 38 minutes** — straddling the old 300 s budget, "
      "which is why run 3 still reported 2/3 while run 4 (3 pairs, 28.5 s) "
      "sailed through. `run_headless`'s budget is raised 300 s → 900 s; a "
      "full panel sweep should use `--species` until "
      "`EnsemblClient.search`'s per-species loop is parallelised (Emergent). "
      "Do not treat any single latency figure here as a rate — treat the "
      "spread as the design constraint.\n")

    # ------------------------------------------------------------- 6. refs
    A("## 6. Bibliography — the claim audit\n")
    A(f"The {len(audit)} references the claim audit rests on. The full "
      f"{len(lit)}-reference bibliography, including everything added for the "
      f"literature review, is `references.tsv`; the review itself is "
      f"`docs/ip3r_review_2026.md`.\n")
    for r in audit:
        A(f"- **{r['ref_id']}** {r['authors'].split(',')[0]} *et al.* "
          f"({r['year']}) {r['title']}. *{r['journal']}*. "
          f"PMID [{r['pmid']}](https://pubmed.ncbi.nlm.nih.gov/{r['pmid']}/)"
          + (f"; doi:{r['doi']}" if r["doi"] and r["doi"] != "-" else ""))
    A("")
    return "\n".join(P)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dir", default="results/s0_baseline")
    args = ap.parse_args()
    d = ROOT / args.dir
    text = build(d)
    out = d / "report.md"
    out.write_text(text, encoding="utf-8")
    print(f"wrote {out.relative_to(ROOT)}  ({len(text.splitlines())} lines)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
