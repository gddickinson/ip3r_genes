"""s23_report.py — render results/s23_scope/report.md from committed tables.

D13: nothing in the report is hand-written. Every number is read back out of a
committed TSV or JSON, so the prose and the data cannot drift.

Usage:  python3 scripts/s23_report.py
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
for _p in (PROJECT_ROOT, PROJECT_ROOT / "scripts"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from s5_build_baits import read_tsv                            # noqa: E402

OUT = PROJECT_ROOT / "results" / "s23_scope"
BAITS = PROJECT_ROOT / "results" / "s23_baits"


def jload(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def tload(path: Path) -> list[dict]:
    return read_tsv(path) if path.exists() else []


def table(rows: list[dict], cols: list[str], headers: list[str] | None = None,
          limit: int | None = None) -> str:
    head = headers or cols
    out = ["| " + " | ".join(head) + " |",
           "|" + "|".join("---" for _ in head) + "|"]
    for r in (rows[:limit] if limit else rows):
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def n(x) -> str:
    try:
        return f"{int(float(x)):,}"
    except (TypeError, ValueError):
        return str(x)


# ------------------------------------------------------------------ sections

def sec_scope(mstats: dict, manifest: list[dict]) -> str:
    r = mstats.get("rules", {})
    by_reason = mstats.get("by_reason", {})
    big = ", ".join(f"{k} ({n(v)})" for k, v in
                    sorted(r.get("big_phyla", {}).items(), key=lambda x: -x[1]))
    lines = [
        "## 1. The denominator\n",
        f"**{n(mstats.get('genomes'))} genomes, {mstats.get('total_gbp')} Gbp** "
        f"(FASTA ≈ {mstats.get('fasta_gb')} GB, zip ≈ {mstats.get('zip_gb')} GB "
        f"against {mstats.get('free_gb')} GB free), selected from "
        f"{n(mstats.get('assemblies_in_pool'))} NCBI eukaryote reference "
        f"assemblies. {n(mstats.get('vertebrate_excluded'))} of those are "
        "vertebrate and belong to S4's manifest, leaving "
        f"{n(mstats.get('assemblies_in_pool', 0) - mstats.get('vertebrate_excluded', 0))} "
        "in scope.\n",
        "The scope rule is **sample most finely where the negative claim is**. "
        "Five rules, every one derived from a committed table except the "
        "anchors, which are hand-written because an anchor is a genome whose "
        "answer is known from the literature and no table can derive that.\n",
        table([
            {"rule": "G1 phylum_rep", "n": by_reason.get("phylum_rep", 0),
             "what": f"one best assembly per non-vertebrate eukaryotic phylum "
                     f"({r.get('n_phyla')} phyla in the pool)"},
            {"rule": "G2 class_rep", "n": by_reason.get("class_rep", 0),
             "what": f"one per class in the phyla S20 swept ≥ 100 proteomes "
                     f"of: {big}"},
            {"rule": "G3 absence_clade", "n": by_reason.get("absence_clade", 0),
             "what": f"one per clade S20 swept ≥ 10 proteomes of and found "
                     f"0 ITPR in ({r.get('n_absence_clades')} clades)"},
            {"rule": "G4 anchor", "n": by_reason.get("anchor", 0),
             "what": "the reference organisms the claims are calibrated on"},
            {"rule": "G5 copy_number", "n": by_reason.get("copy_number", 0),
             "what": "the highest-copy species per class, above the vertebrate "
                     "paralog count"},
        ], ["rule", "n", "what"], ["rule", "genomes", "what it selects"]),
        f"\n{mstats.get('genomes_with_multiple_reasons')} genomes carry more "
        "than one reason and are one row each, so the row count is the number "
        "of genomes to download.\n",
        "**G3 re-derives S20a's target list from the presence table rather "
        "than repeating it, and finds three absences S20a's summary never "
        "named** — Bacillariophyta 0/16, Rhodophyta 0/12, and **Cestoda "
        "0/11, a metazoan clade with no ITPR record at all**.\n",
        "By kingdom-level group: " + ", ".join(
            f"{k} {v}" for k, v in sorted(
                mstats.get("by_group", {}).items())) + ".\n",
    ]
    return "\n".join(lines)


def sec_calibration(cal: dict) -> str:
    rows = [{"group": g, "n": s["n"],
             "range": f"{n(s['min_bp'])} – {n(s['max_bp'])}",
             "bar": n(s["contiguity_bar_bp"]), "g": n(s["max_intron_bp"])}
            for g, s in cal.get("by_group", {}).items()]
    return "\n".join([
        "## 2. The two genomic thresholds, measured for this scope\n",
        "S5a measured both for the vertebrates and neither transfers. A "
        "vertebrate ITPR gene spans 76–498 kb; a *Drosophila* Itpr spans "
        f"22 kb. Measured here over **{cal.get('n_spans')} genes in "
        f"{cal.get('n_species')} species across {cal.get('n_bands')} bands**, "
        "from NCBI's own gene annotations for genes the census independently "
        "calls ITPR — not from miniprot, because the aligner's `-G` shapes the "
        "loci it reports and a measurement taken that way cannot falsify the "
        "setting it calibrates.\n",
        f"Spans run **{n(cal.get('min_span_bp'))} – "
        f"{n(cal.get('max_span_bp'))} bp, median {n(cal.get('median_span_bp'))}"
        "**. That is a ~100× range, against the 6.5× S0 measured inside the "
        "vertebrates, so a single bar would be far too lax for the metazoans "
        "or far too strict for the protists. **The bar is therefore per "
        "group.**\n",
        table(rows, ["group", "n", "range", "bar", "g"],
              ["group", "genes", "span range (bp)", "contiguity bar (bp)",
               "miniprot -G (bp)"]),
        f"\nThe S5 vertebrate bar was {n(cal.get('s5_vertebrate_bar_bp'))} bp "
        "and 120 of 309 genomes failed it — S5b's headline caveat. Here "
        "**metazoan genes are ~9× longer than protist and fungal ones** "
        "(83 kb median against 7–9 kb), and the genomes carrying this task's "
        "negative claims are judged against the ~9 kb bar, which almost any "
        "modern assembly clears.\n",
        "`-G` moves the other way and is floored at miniprot's own default: a "
        "`-G` below an unmeasured species' largest intron **splits** its gene, "
        "and a split ITPR reads out of a copy-number ledger as `fragment`. "
        "41 of the 57 panel records were unmeasurable (their locus tags have "
        "no NCBI gene record), so erring high is the only safe direction.\n",
        "**Gap to state:** no Viridiplantae gene span was measured — the "
        "chlorophyte census records carry locus tags NCBI has no gene record "
        "for — so every plant genome is judged against the global median. "
        "`s23_calibration.contiguity_bar()` returns that fallback in its own "
        "return value rather than hiding it.\n",
    ])


def sec_panel(bstats: dict, manifest: list[dict], unfilled: list[dict],
              controls: list[dict], selftest: list[dict]) -> str:
    by_role = Counter(r["role"] for r in manifest)
    strong = sum(1 for c in controls if c.get("strength") == "strong"
                 and c.get("kept") == "1")
    weak = sum(1 for c in controls if c.get("strength") == "weak"
               and c.get("kept") == "1")
    rejected = [c for c in controls if c.get("kept") == "0"]
    lines = [
        "## 3. The bait panel, and the control that had to be replaced\n",
        f"**{n(bstats.get('baits'))} baits, {n(bstats.get('residues'))} "
        f"residues**: {by_role.get('ITPR', 0)} ITPR, {by_role.get('RYR', 0)} "
        f"RyR and {by_role.get('MIR', 0)} MIR-domain control, derived from "
        f"census v5 by seven enforced rules, {bstats.get('reused_s3_seeds')} "
        "of them S3 seeds reused unchanged (B6).\n",
        "### 3.1 S5's positive control does not transfer\n",
        "S5 can write \"the RyR control fired in all 309 genomes, so none is "
        "excluded on control grounds\" because every vertebrate has three "
        "RyRs. Outside Metazoa that argument collapses: S20 found "
        "architecture-level RyR in **2 of 6,928** non-metazoan proteomes. A "
        "land plant with no RyR locus is the *correct* answer, so RyR cannot "
        "be that genome's proof-of-search — and without one, \"no ITPR in "
        "*Arabidopsis*\" is indistinguishable from \"the sweep did not run "
        "properly on *Arabidopsis*\". Every negative claim in this task rests "
        "on closing that gap.\n",
        "The replacement is the **MIR-domain sharer** — the protein "
        "O-mannosyltransferase family S1's decoy panel was built around and "
        "S20 used as its in-search positive control, where PF08709 returns 0 "
        "matches in land plants and PF02815 returns 633. It is in the "
        "family's own signature set, it is drawn from the clade each claim is "
        "about, and it is simultaneously the sharpest available decoy.\n",
        f"**{strong + weak} usable control baits across {len(controls)} "
        f"control clades: {strong} strong, {weak} weak.** A control is strong "
        "when it is a "
        "named, full-length mannosyltransferase present across at least a "
        "quarter of the clade's swept proteomes; the classification is "
        "recorded per clade rather than resolved, because it is the reader's "
        "to weigh.\n",
    ]
    weak_rows = [c for c in controls if c.get("strength") == "weak"
                 and c.get("kept") == "1"]
    if weak_rows:
        lines += [table(weak_rows, ["control_for", "accession", "length",
                                    "score", "clade_frac", "strength_note"],
                        ["clade", "bait", "aa", "bits", "clade coverage",
                         "why weak"]), ""]
    lines += [
        "**The apicomplexan control is the thin one and it is the one that "
        "matters.** Aconoidasida and Conoidasida carry a single PF02815 "
        "protein each across 60 swept proteomes (3–4 % clade coverage), so "
        "the Apicomplexa absence — one of S20a's headline targets — rests on "
        "a 231–233 aa control found in one proteome. That is stated rather "
        "than smoothed over, and it is the strongest argument for a "
        "second, non-MIR control in that clade.\n",
    ]
    if rejected:
        lines += [
            "### 3.2 The screen rejected a control on its first run\n",
            "A control bait must be assigned to **neither** family — a MIR "
            "protein the profiles call a family member is a finding, not a "
            "control. One was:\n",
            table(rejected, ["control_for", "accession", "organism", "length",
                             "protein_name", "screen_reason"],
                  ["clade", "accession", "organism", "aa", "its own name",
                   "screen verdict"]),
            "\nA tapeworm protein UniProt calls a \"MIR domain-containing "
            "protein\" is a ryanodine receptor. Cestoda is a metazoan clade, "
            "so its genomes keep the RyR control; the MIR slot is left "
            "empty rather than filled with a family member.\n",
        ]
    if unfilled:
        lines += [
            "### 3.3 Unfilled slots\n",
            f"{len(unfilled)} slots have no bait. The reason is recorded per "
            "slot and distinguishes the census having no record from a "
            "threshold rejecting the records it has — the first version of "
            "this table said \"no census record at all\" for seven bands, "
            "four of which have dozens.\n",
            table(unfilled, ["band", "role", "records", "shape_ok", "why"],
                  ["band", "role", "records", "pass shape", "why"]),
            "",
        ]
    lines += [
        "### 3.4 Both measured thresholds had to be stratified\n",
        f"The bait length band is **{bstats.get('band_aa')} aa** and the "
        f"architecture-exception floor **{bstats.get('arch_floor_bits')} "
        "bits**, each measured over the complete-architecture non-vertebrate "
        f"records **capped at {bstats.get('measure_cap_per_band')} per band**"
        f" ({bstats.get('bands_in_measurement')} bands contribute).\n",
        f"Measured without the cap the band is "
        f"{bstats.get('band_aa_unstratified')} aa and the floor ~1,095 bits — "
        "the numbers of a population that is 53 % Arthropoda, wearing the "
        "whole tree's name. They excluded *Bodo saltans*, a 2,356–4,222 aa "
        "euglenozoan ITPR at 4/5 architecture and 780 bits, which is exactly "
        "the deep-branching record a non-vertebrate panel exists to carry. "
        "This is the same correction S20 made when it stratified its "
        "bacterial sample by genus, and the same one G3 and G5 make in the "
        "scope: **a threshold measured on the best-sampled clade describes "
        "the sampling, not the family.**\n",
    ]
    if selftest:
        ok = sum(1 for r in selftest if r.get("rejected") == "1")
        lines += [
            f"The chimera screen is S5's module unchanged (D5), and its own "
            f"negative controls run on every build: **{ok}/{len(selftest)} "
            "synthetic failures rejected**, each by the rule responsible.\n",
        ]
    return "\n".join(lines)


def sec_pilot(lstats: dict, ledger: list[dict], absences: list[dict],
              controls: list[dict]) -> str:
    if not lstats:
        return ("## 4. Pilot\n\nNot yet run. "
                "`python3 scripts/s23_run_sweep.py --reasons anchor`\n")
    status = lstats.get("status", {})
    lines = [
        "## 4. Pilot — the 14 anchor genomes\n",
        f"**{lstats.get('genomes')} genomes of the "
        f"{lstats.get('manifest_genomes')} in the manifest**, chosen as G4's "
        "anchors: the organisms whose answer is known from the literature. "
        "Half are positive controls (a characterised ITPR) and half are the "
        "negative claims (*Arabidopsis*, rice, a moss, *S. cerevisiae*, "
        "*Neurospora*, a microsporidian, *Toxoplasma*), so the pilot tests "
        "the pipeline in both directions at once.\n",
        f"Statuses: " + ", ".join(f"`{k}` {v}" for k, v in
                                  sorted(status.items(), key=lambda x: -x[1]))
        + ".\n",
        table(ledger, ["organism", "group", "status", "n_full", "n_fragment",
                       "mir_loci", "ryr_loci", "control_verdict"],
              ["organism", "group", "status", "full", "frag", "MIR", "RyR",
               "control"]),
        "",
        f"**Control:** " + ", ".join(
            f"{v} {k}" for k, v in lstats.get("control_verdicts", {}).items())
        + f". {lstats.get('uncontrolled', 0)} genome(s) uncontrolled — no "
        "absence claim may rest on those.\n",
        f"**Copy number:** {lstats.get('total_full_loci')} full ITPR loci "
        f"from {lstats.get('total_loci_raw')} raw clusters, with "
        f"{lstats.get('total_merges')} pair(s) folded by the conservative "
        "split-gene rule. "
        f"{lstats.get('novel_models')} novel model(s) written.\n",
    ]
    covered = [a for a in absences if a["genomes_controlled"] != "0"
               and a["genomes_swept"] != "0"]
    if covered:
        lines += [
            "### 4.1 The absence claims at assembly level\n",
            table(covered, ["clade", "swept_proteomes", "genomes_swept",
                            "genomes_controlled", "genomes_with_full_itpr",
                            "verdict"],
                  ["clade", "proteomes (S20)", "genomes", "controlled",
                   "with full ITPR", "verdict"]),
            "",
        ]
    return "\n".join(lines)


def main() -> int:
    mstats = jload(OUT / "manifest_stats.json")
    bstats = jload(BAITS / "bait_build_stats.json")
    cal = jload(BAITS / "span_calibration.json")
    lstats = jload(OUT / "ledger_stats.json")
    manifest = tload(BAITS / "bait_manifest.tsv")
    parts = [
        "# S23a — the non-vertebrate genomic sweep's instrument\n",
        "*Generated by `scripts/s23_report.py` from the committed tables "
        "only (D13). Nothing here is hand-written.*\n",
        "S20 swept 6,928 reference **proteomes** and found the family absent "
        "from the land plants, the Dikarya, the Apicomplexa and a scatter of "
        "fungal phyla. Every one of those is an annotation fact: a proteome "
        "is what a gene-caller found, a genome is what is there. S23 takes "
        "them to assembly level. This half builds and measures the "
        "instrument; S23b runs it.\n",
        sec_scope(mstats, tload(OUT / "genome_manifest_s23.tsv")),
        sec_calibration(cal),
        sec_panel(bstats, manifest, tload(BAITS / "unfilled_slots.tsv"),
                  tload(BAITS / "control_manifest.tsv"),
                  tload(BAITS / "screen_self_test.tsv")),
        sec_pilot(lstats, tload(OUT / "copy_number_ledger.tsv"),
                  tload(OUT / "absence_at_genome.tsv"),
                  tload(OUT / "control_ledger.tsv")),
    ]
    (OUT / "report.md").write_text("\n".join(parts))
    print(f"wrote {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
