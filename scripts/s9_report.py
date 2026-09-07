"""S9 — renders `results/selection/report.md` purely from the committed
tables (D13). Scope and instrument here; the results in
`s9_report_results.py` (the `s3_report.py` / `s3_report_d10.py` split).

Nothing in this file re-parses an `mlc`, re-runs a model or recomputes a
likelihood ratio. Every number it prints comes out of a TSV or a JSON that
`s9_tables.py` wrote, so the report and the data cannot drift.

Run:  python scripts/s9_report.py
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from s9_cds_lib import OUT_DIR, PARALOGS  # noqa: E402

REPORT = OUT_DIR / "report.md"


def load_tsv(name: str) -> list[dict]:
    path = OUT_DIR / name
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    path = OUT_DIR / name
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:  # noqa: BLE001
        return {}


def num(x, fmt: str = "{:.4g}", dash: str = "—") -> str:
    try:
        return fmt.format(float(x))
    except (TypeError, ValueError):
        return dash


def pct(a: float, b: float) -> str:
    return f"{100 * a / b:.1f} %" if b else "—"


def _scope(status: list[dict], codon_stats: str, tips: list[dict]) -> list[str]:
    ok = [r for r in status if r["status"] == "ok"]
    routes: dict[str, int] = {}
    for r in ok:
        routes[r["route"]] = routes.get(r["route"], 0) + 1
    n_multi = sum(1 for r in ok if int(r.get("n_candidates") or 1) > 1)
    n_masked = sum(int(r["n_masked"]) for r in ok
                   if (r.get("n_masked") or "-").lstrip("-").isdigit())
    n_stop = sum(int(r["n_stop_masked"]) for r in ok
                 if (r.get("n_stop_masked") or "-").lstrip("-").isdigit())
    sets: dict[str, int] = {}
    for t in tips:
        sets[t["set"]] = sets.get(t["set"], 0) + 1

    out = ["## 1. Scope — what a codon alignment of this family can be", "",
           "dN/dS is a statement about codons, so every number below rests "
           "on a nucleotide sequence that provably encodes the *exact* "
           "protein S6 aligned and S7 built its tree from. Nothing here is "
           "a translated database record taken on trust.", "",
           f"- vertebrate family tips in S6's representative set: "
           f"**{len(status)}**",
           f"- validated CDS recovered: **{len(ok)} / {len(status)}** "
           f"({pct(len(ok), len(status))})",
           "- by route: " + ", ".join(f"`{k}` {v}" for k, v in
                                      sorted(routes.items())),
           f"- tips where the first CDS a route offered was **not** the "
           f"aligned isoform: **{n_multi}**",
           f"- codons masked to `NNN` across the whole set: **{n_masked}**, "
           f"of which {n_stop} are internal stops in genome gene models",
           ""]
    out += ["The RyR outgroup is deliberately absent. It roots S7's protein "
            "tree, but synonymous sites do not survive that distance — and "
            "as §4.2 shows, they barely survive the vertebrate span *inside* "
            "a paralog. A model fitted across the family/outgroup split "
            "would be estimating alignment error.", ""]
    out += ["| selection set | tips |", "|---|---|"]
    for k in list(PARALOGS) + ["background"]:
        if k in sets:
            out.append(f"| {k} | {sets[k]} |")
    out += ["",
            "The three paralog sets are **S7's extended paralog clades**, "
            "not the census labels. codeml's branch models mark a *node*: a "
            "foreground that is not a clade does not fail, it silently marks "
            "a larger one and returns a well-formed ω for a hypothesis "
            "nobody asked. `s9_sets.py` re-derives each clade from "
            "`rooted.nwk` with S7's own rule and aborts unless it matches "
            "`paralog_clades.tsv` in size and membership.", ""]
    reassigned = [t for t in tips
                  if t["set"] in PARALOGS and t["census_group"] != t["set"]]
    if reassigned:
        out += [f"The tree nests **{len(reassigned)}** `vertebrate_basal` "
                "loci carrying no paralog label of their own inside a "
                "paralog clade, and they are analysed there. "
                "That is D30 read forwards — outside the vertebrates a "
                "paralog label means nothing, but inside a "
                "100/100-supported vertebrate paralog clade the tree has "
                "just supplied one.", ""]
    bg = [t for t in tips if t["set"] == "background"]
    if bg:
        # The species column, not a slice of the label: the S6 label's
        # group prefix itself contains an underscore, so splitting on "_"
        # reported the cyclostomes as "basal Myxine".
        species: dict[str, int] = {}
        for t in bg:
            key = (t.get("species") or t["label"]).split(" (")[0]
            species[key] = species.get(key, 0) + 1
        listed = ", ".join(f"*{k}* ×{v}" for k, v in sorted(species.items()))
        out += [f"**{len(bg)}** tips sit in no paralog clade at all. They "
                "stay in the whole-tree analyses, because dropping them "
                "would change the branch lengths every other estimate is "
                "made on, and they are in no foreground, because a locus "
                "the tree could not place is not evidence about a paralog.",
                "",
                f"They are {listed} — which is the *same six loci* S7 placed "
                "in two cyclostome-only clades and handed to S8, and that "
                "S8 reported as underpowered because after ~550 Myr there "
                "is no shared flank vocabulary left to compare. A third "
                "instrument now declines the same question for a third "
                "reason: a codon model can only ask about a paralog whose "
                "clade the tree defines, and for these six it defines "
                "none.", ""]
    if codon_stats:
        out += ["### 1.1 The codon alignment", "", codon_stats.strip(), ""]
    return out


def _instrument(stats: dict, omega: list[dict]) -> list[str]:
    jobs = {r["job"]: r for r in omega}
    total_min = sum(float(r["minutes"]) for r in omega
                    if (r.get("minutes") or "").replace(".", "", 1).isdigit())
    out = ["## 2. Instrument — the models, and the traps in them", "",
           "PAL2NAL builds the codon alignment and an **independent "
           "in-house protein→codon mapping is computed beside it**; a "
           "single nucleotide of disagreement aborts the build. trimAl's "
           "`-automated1` columns are chosen on the protein and applied "
           "codon-aware, whole triplets only.", "",
           f"- codeml jobs run: **{len(omega)}**"
           + (f", {total_min / 60:.1f} h of CPU time" if total_min else ""),
           f"- likelihood-ratio tests: **{stats.get('n_lrt', 0)}**, "
           "Benjamini–Hochberg corrected across the whole family",
           f"- branch-site model A restarted from initial ω = "
           + ", ".join(str(w) for w in stats.get("bs_init_omegas", []))
           + " on every paralog stem",
           f"- pairwise dS above **{stats.get('saturated_ds_bar', 1.5)}** is "
           "flagged saturated",
           ""]
    out += ["Three of those are decisions, not settings.", "",
            "**Branch-site model A is restarted by construction.** A nested "
            "alternative cannot have a lower optimum than its own null, yet "
            "codeml reaches one routinely on alignments this size — the "
            "PIEZO project hit exactly that and had to add restarts after "
            "the fact. Running four initial ω from the start makes \"the "
            "best of several optima\" the reported number rather than a "
            "repair, and `bs_restarts.tsv` commits the spread, so a stem "
            "that is *still* stuck is visible in the data.", "",
            "**The branch-site LRT is tested against a 50:50 mixture.** Its "
            "null fixes ω₂ = 1 on the boundary of the parameter space, so "
            "2ΔlnL is distributed as ½χ²₀ + ½χ²₁ and a plain χ²₁ p-value is "
            "twice too small. `lrt_table.tsv` carries both; the halved one "
            "is reported.", "",
            "**BH runs across the whole LRT family.** S9 asks the same "
            "question three times, once per paralog. Reporting the smallest "
            "of three p-values uncorrected is the multiple-testing error "
            "this project's self-tests exist to avoid.", ""]
    if "m0_all" in jobs:
        r = jobs["m0_all"]
        out += ["### 2.1 The tree is an input, not a result", "",
                "Every branch model is run on S7's topology pruned to the "
                "tips with a validated CDS, unrooted for PAML. S9 fits "
                "branch lengths to it and never re-estimates it, so a "
                "selection result cannot quietly re-open the sister "
                f"question. Whole-tree one-ratio: lnL {num(r['lnL'], '{:.1f}')}"
                f", tree length {num(r['tree_length'])}, "
                f"κ {num(r['kappa'])}.", ""]
    return out


def _selftest_section() -> list[str]:
    return ["## 3. Self-tests", "",
            "Twelve constructed negative controls run on every build of the "
            "codon alignment (`s9_test_codon.py`), each rejected by the rule "
            "responsible: a frame-shifted CDS, a CDS encoding a different "
            "protein, an internal stop that must be masked rather than "
            "refused, the requirement that an accepted CDS translate to its "
            "aligned protein modulo `X`, the gap-exactness of the in-house "
            "codon map, whole-triplet gap stripping, **a non-monophyletic "
            "foreground being refused**, `$1`-vs-`#1` labelling, unrooting, "
            "the locus-rerun difference bound in all three directions, and "
            "the S7 cross-check being able to fail at all. A codon "
            "alignment is the one artefact here where a silent error is "
            "invisible downstream — codeml will fit a model to a "
            "frame-shifted alignment and return a plausible ω — so these "
            "are checks on *refusal*, not on output appearing.", ""]


def main() -> None:
    status = load_tsv("cds_status.tsv")
    tips = load_tsv("tip_codes.tsv")
    omega = load_tsv("omega_table.tsv")
    stats = load_json("selection_stats.json")
    codon_stats = ""
    p = OUT_DIR / "codon_aln_stats.md"
    if p.exists():
        codon_stats = p.read_text().split("\n", 1)[1]
    if not omega:
        raise SystemExit("no omega_table.tsv — run scripts/s9_tables.py")

    from s9_report_results import results_sections  # noqa: PLC0415

    lines = ["# S9 — ML selection across the ITPR family", "",
             "Rendered from the committed tables by `scripts/s9_report.py` "
             "(D13). Every number below is in a TSV beside this file.", ""]
    lines += _scope(status, codon_stats, tips)
    lines += _instrument(stats, omega)
    lines += _selftest_section()
    lines += results_sections(load_tsv, load_json, num, pct)
    REPORT.write_text("\n".join(lines) + "\n")
    print(f"-> {REPORT} ({len(lines)} lines)")


if __name__ == "__main__":
    main()
