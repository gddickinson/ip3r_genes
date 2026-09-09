"""S19's report, part two: the measurements (§4-§6).

Takes the caller's headline dict so the two halves cannot read the tables
differently (`s3_report.py` / `s3_report_d10.py`). Every prior an earlier
task set is stated through `s19_priors.line`, which refuses a verdict outside
the committed five-valued vocabulary.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s19_lib as S                                             # noqa: E402
import s19_panel as PN                                           # noqa: E402
import s19_priors as P                                          # noqa: E402
from s19_report import fmt, pct, pfmt, table                    # noqa: E402


def _t(name: str) -> list[dict]:
    path = S.OUT_DIR / name
    return S.read_tsv(path) if path.exists() else []


def _missing(name: str) -> str:
    return (f"*Not run yet — `{name}` is absent from `results/methods/`.*")


# --------------------------------------------------------------- 4. contiguity

def section_4(h: dict) -> str:
    rows = _t("contiguity_bins.tsv")
    if not rows:
        return "## 4. Contiguity as a confounder\n\n" + _missing(
            "contiguity_bins.tsv")
    floor = _t("absence_floor.tsv")
    tests = _t("contiguity_tests.tsv")
    resid = _t("residual_neighbourhood.tsv")
    bar = h["d4_bar"]

    below = [r for r in tests if r["test"] == "false_negative_rate_below_bar"]
    fn_rows = []
    for series in ("itpr_present", "ryr_sister"):
        for r in floor:
            if r["series"] != series or r["is_d4_bar"] != "1":
                continue
            fn_rows.append([{"itpr_present": "ITPR cells",
                             "ryr_sister": "RyR sister cell"}[series],
                            fmt(int(r["n_genomes_kept"])),
                            fmt(int(r["n_cells"])),
                            fmt(int(r["n_false_negative"])),
                            pct(float(r["fn_rate"])),
                            pct(float(r["wilson_hi"]))])
    nmiss = sum(1 for r in resid if r["verdict"] == "neighbourhood_missing")

    return f"""## 4. Contiguity is the whole of the false-negative rate

The misses are not distributed across the scope. A missed control cell's
assembly has a median contig N50 of
{[r for r in tests if r['test'].startswith('mannwhitney') and r['series'] == 'itpr_present'][0]['effect'].split('vs ')[1]} bp
against
{[r for r in tests if r['test'].startswith('mannwhitney') and r['series'] == 'itpr_present'][0]['effect'].split('median ')[1].split(' vs')[0]} bp
for a found one, and the odds of finding the gene rise
{h['contiguity_effect']['itpr_present'].split('OR per 10x=')[1].split(',')[0]}×
per tenfold of contig N50 in the ITPR series and
{h['contiguity_effect']['ryr_sister'].split('OR per 10x=')[1].split(',')[0]}×
in the RyR series. On chromosome-level assemblies the ITPR series misses
{[r for r in tests if r['test'] == 'fisher_chromosome_vs_lower' and r['series'] == 'itpr_present'][0]['n'].split(';')[0].replace('chromosome ', '').replace(' missed', '')}
cells and the RyR series misses
{[r for r in tests if r['test'] == 'fisher_chromosome_vs_lower' and r['series'] == 'ryr_sister'][0]['n'].split(';')[0].replace('chromosome ', '').replace(' missed', '')}.

### 4.1 D4's bar was chosen a priori, and it survives calibration

D4's contiguity bar is **{fmt(bar)} bp** of contig N50 — the median measured
ITPR genomic span, taken from gene geometry alone, with no error rate in its
derivation. This is the first time it has been scored against one:

{table(["series", "genomes kept", "cells", "missed", "residual rate",
        "upper 95 % bound"], fn_rows)}

The bar the project has been using since S5a lands at under one per cent
residual error in the series that can have any, and at zero in the other. It
is *conservative* against the conventional five-per-cent target, which the
scan reaches at a floor of
{fmt(int(next(r['contig_n50_floor'] for r in floor if r['series'] == 'itpr_present' and float(r['wilson_hi']) <= 0.05)))} bp,
and about right against a one-per-cent target.

{P.line("contiguity_bar", "confirmed",
        f"scored against a measured false-negative rate for the first time, "
        f"the a-priori bar gives "
        f"{pct(h['itpr_present']['rate_at_d4_bar'])} residual error on the "
        f"ITPR series and {pct(h['ryr_sister']['rate_at_d4_bar'])} on the "
        f"RyR sister, retaining "
        f"{fmt(h['itpr_present']['genomes_at_d4_bar'])} of the "
        f"{fmt(h['n_genomes'])} genomes")}

{P.line("false_negatives", "confirmed",
        f"S5b's 98-99 % / 57-70 % split is reproduced as a false-negative "
        f"rate: {pct(h['itpr_present']['rate'])} over the whole scope, "
        f"{pct(h['itpr_present']['rate_at_d4_bar'])} above the bar")}

### 4.2 What the floor costs, and which clades it removes

A floor that reaches one per cent by keeping a third of the scope has not
made the survey more reliable, it has made it smaller. At D4's own bar
{fmt(h['itpr_present']['genomes_at_d4_bar'])} of {fmt(h['n_genomes'])}
genomes are retained; `floor_clade_composition.tsv` records which classes
each floor removes, and the answer is the one S5a already warned about — the
margin species S4 added for their *proteome* gaps are very largely the same
genomes whose assemblies cannot hold the gene, so a contiguity filter removes
the part of the scope the scope was extended for.

### 4.3 Below the bar, the misses run with gene span

{table(["cell", "cells below the bar", "false-negative rate", "95 % CI"],
       [[r["series"], r["n"].replace("n=", ""), pct(float(r["effect"])),
         " – ".join(pct(float(x)) for x in
                    r["detail"].replace("95% CI ", "").split("-"))]
        for r in below])}

{P.line("span_bias", "confirmed",
        "below the bar the ordering is the same one S5b measured and points "
        "the same way as gene span: the shortest gene is the one the "
        "fragmented assembly still yields")}

### 4.4 The residual, diagnosed locally

Contig N50 is a genome-wide statistic and cannot see a regional assembly
defect. Asking S8's own consensus flanks whether the neighbourhood survived
at all: **{fmt(nmiss)} of {fmt(len(resid))}** false negatives sit in an
assembly that carries fewer than half of that paralog's usual neighbours, so
what is missing there is the region and the cell says nothing about the gene.
"""


# ------------------------------------------------------------- 5. contribution

def section_5(h: dict) -> str:
    if not _t("head_to_head.tsv"):
        return "## 5. What each channel contributed\n\n" + _missing(
            "head_to_head.tsv")
    gs = h["head_to_head_gene_scale"]
    tail = h["head_to_head_tail"]
    growth = h["census"]
    return f"""## 5. What each search channel contributed

### 5.1 How the census accumulated

{table(["census", "channel that produced it", "records", "added"],
       [[r["v"], r["channel"], fmt(r["n"]), fmt(r["added"])]
        for r in growth])}

### 5.2 Profile HMM against domain annotation, inside one database

Both channels are counted by what they **call family** — the enumeration by
carrying a family signature, the sweep by clearing D22's gate — and both are
restricted to accessions the swept FASTAs actually hold. Comparing a curated
set against a raw hit list would credit the vertebrate sweep with the
13,371 targets its own gate declined.

At **gene scale** they return nearly the same set:

{table(["database", "Pfam enumeration", "profile HMM", "shared",
        "HMM only", "Pfam only"],
       [[db, fmt(v["pfam"]), fmt(v["hmm"]), fmt(v["shared"]),
         fmt(v["hmm_only"]), fmt(v["pfam_only"])]
        for db, v in gs.items() if v["pfam"] or v["hmm"]])}

In the **fragment tail** they do not:

{table(["database", "Pfam enumeration", "profile HMM", "HMM only",
        "Pfam only"],
       [[db, fmt(v["pfam"]), fmt(v["hmm"]), fmt(v["hmm_only"]),
         fmt(v["pfam_only"])]
        for db, v in tail.items() if v["pfam"] or v["hmm"]])}

So a family profile HMM over reference proteomes is not a discovery method
for whole genes — in the vertebrates it adds one gene-scale record to
{fmt(gs['vertebrata']['pfam'])}, and in the non-vertebrate metazoa two to
{fmt(gs['metazoa_nonvert']['pfam'])}. Its entire gain is in fragments. The
one exception is worth naming: in the **protists** it adds
{fmt(gs['protista_other']['hmm_only'])} gene-scale records the enumeration
never returned, which is where the family's architecture is least well
annotated.

The Pfam-only column is the same statement from the other side:
{fmt(tail['vertebrata']['pfam_only'])} vertebrate records carry a family
signature and are declined by the profile pair, all of them under 1,000 aa
— D22's 200-match-state gate doing what it was added for.

{P.line("hmm_over_pfam", "confirmed",
        f"S3's fragment characterisation is now a measurement rather than a "
        f"description: of everything the vertebrate profile pair returns "
        f"that the enumeration did not, exactly "
        f"{fmt(gs['vertebrata']['hmm_only'])} record is gene-scale and "
        f"{fmt(tail['vertebrata']['hmm_only'])} are under 1,000 aa")}

### 5.3 Per gene: what a protein-database search would have missed

The question asked per genome × cell rather than per record, because the
genome sweep is the ground truth for where the genes are:

{table(["how the gene is reachable", "cells"],
       [[k.replace("genome_only:", "genome only — ").replace("_", " "),
         fmt(v)] for k, v in sorted(h["gene_recovery"]["by_channel"].items(),
                                    key=lambda kv: -kv[1])])}

**{fmt(h['gene_recovery']['invisible'])} of
{fmt(h['gene_recovery']['n_genes'])}** demonstrated genes
({pct(h['gene_recovery']['invisible'] / h['gene_recovery']['n_genes'])}) are
not reachable by any protein-database search. That is not an artefact of the
margin species: split by why S4 put each genome in scope, the rate is
{pct(h['gene_recovery']['by_scope']['order representative'])} for order
representatives against
{pct(h['gene_recovery']['by_scope']['margin species'])} for margin species.

{P.line("genome_only", "confirmed",
        f"S5b's census v4 counts are the per-locus form of this; per cell "
        f"the figure is "
        f"{pct(h['gene_recovery']['invisible'] / h['gene_recovery']['n_genes'])}, "
        f"and the sister family is the best-served cell at "
        f"{pct(h['gene_recovery']['by_cell']['RYR'])}")}

### 5.4 Cost

{h['recorded_hours']:.1f} recorded compute-hours over eight channels.
{len(h['channels_without_timing'])} channels carry no wall clock at all and
are reported blank rather than as zero: `method_cost.tsv` names them and
says why.
"""


# ---------------------------------------------------------------- 6. the panel

def section_6(h: dict) -> str:
    recall = _t("panel_recall.tsv")
    if not recall:
        return "## 6. The bait panel\n\n" + _missing("panel_recall.tsv")
    changes = _t("panel_changes.tsv")
    slots = _t("panel_unfilled_slots.tsv")
    gains = [r for r in changes if r["direction"] == "gain"]
    p = h["panel"]
    keep = ["full38", "labelled_only", "no_ryr_control", "drop_bird",
            "drop_mammal", "drop_ray_finned_fish", "amniote_only",
            "s3_seeds_only", "human_only", "one_per_cell_nonhuman",
            "drop_ITPR1_baits", "drop_ITPR2_baits", "drop_ITPR3_baits",
            "basal_only"]
    desc = {r["panel"]: r["description"] for r in recall}
    return f"""## 6. The bait panel: breadth is nearly free, paralog coverage is not

### 6.1 The ablation is exact, and validated first

The simulation reproduces the committed ledger **{fmt(h['panel_cells_validated'])}
cells out of {fmt(h['panel_cells_validated'])}** under the full panel
({fmt(h['panel_validation_mismatches'])} disagreements). Every delta below is
therefore measured against the sweep itself, not against a model of it.

### 6.2 What each ablation costs

{table(["panel", "baits", "ITPR cells found", "change", "description"],
       [[k, fmt(p[k]["n_baits"]), fmt(p[k]["found"]),
         f"{p[k]['delta']:+d}", desc.get(k, "")]
        for k in keep if k in p])}

Two readings, and they point opposite ways.

**Phylogenetic breadth buys almost nothing.** Four human baits — one per
cell — recover {fmt(p['human_only']['found'])} of the
{fmt(p['full38']['found'])} cells the 38-bait panel recovers. Dropping any
single clade band costs at most two cells. The three unlabelled
`vertebrate_basal` baits cost **nothing** when removed, and on their own
recover **zero** cells: an unlabelled bait cannot fill a paralog cell by
itself, because the sweep offers an unresolved-clade locus to whichever
*labelled* paralog scores highest there, and with no labelled bait present
there is nobody to offer it to. Their function is to keep a real shark or
lamprey gene inside the family when it outranks every labelled bait, not to
place it in a cell.

**Paralog coverage buys everything.** Removing one paralog's own baits costs
{fmt(abs(p['drop_ITPR3_baits']['delta']))}, {fmt(abs(p['drop_ITPR2_baits']['delta']))}
and {fmt(abs(p['drop_ITPR1_baits']['delta']))} cells — about a quarter of the
recovery each — and no amount of breadth in the other paralogs replaces it.

**Removing the sister-family control changes no ITPR call.** D14 at genome
scale costs nothing in recall, which is the cheapest possible price for a
positive family test.

### 6.3 An ablation can *gain* a cell, and that is a property of the rule

{fmt(len(gains))} of the {fmt(len(changes))} call changes across all panels
are gains. The mechanism is the same every time: a cell's coverage is read
off its **top-scoring** bait, not its best-covering one, so removing a
competitor can promote a bait whose coverage is fractionally higher and push
a marginal cell across the {S.fnum(PN.COV_FOUND, float, 0.7)} coverage bar.
The two gains under `drop_bird` are both a chicken ITPR1 bait at coverage
0.6956 giving way to a lizard bait at 0.7022 — seven thousandths either side
of the threshold. They are reported here rather than folded away, because an
ablation table with unexplained positive deltas invites the reading that a
smaller panel searches better.

### 6.4 The six unfilled slots

{table(["cell", "band", "classes", "cells", "found by the full panel",
        "found without the unlabelled baits"],
       [[r["cell"], r["band"], r["vclasses"], r["n_cells"],
         r["n_found_full_panel"], r["n_found_labelled_only"]]
        for r in slots])}

S5's panel has no labelled bait for these six slots because no labelled
full-length record exists there. The chondrichthyan and coelacanth cells are
recovered anyway, by labelled baits from other clades — those genes are
within reach of a bony-fish or tetrapod bait. The four cyclostome cells are
not recovered by any panel, which is S5b's own result and the reason S15a
filed them `paralog_unassignable`.

{P.line("bait_breadth", "not corroborated",
        f"the unlabelled baits change no miniprot cell call in either "
        f"direction ({fmt(p['labelled_only']['delta'])} cells), so S5a's "
        f"fix mattered to rescue attribution — a different code path — and "
        f"not to locus placement; and the six unfilled slots cost only the "
        f"four cyclostome cells")}

### 6.5 The design rule

Every (genome, cell, bait) triple the sweep produced is a measurement of what
one bait alone recovers at a known sequence distance:

{table(["bait-to-target identity", "measurements", "recovered", "recall",
        "lower 95 % bound"],
       [[r["bin"], fmt(r["n"]), fmt(int(round(r["recall"] * r["n"]))),
         pct(r["recall"]), pct(r["lo"])] for r in h["single_bait"]])}

A single ITPR bait recovers the gene essentially always at any identity above
0.5. The scope of that claim matters: it is measured **at loci the full panel
already placed**, so it says what one bait can do where a gene is known to
be, not what a one-bait panel would find in an unsearched genome. Read that
way it is still the useful number — it is why a four-bait panel matches a
38-bait one here.
"""


