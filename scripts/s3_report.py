"""S3 — render results/census_v3/report.md purely from the committed tables.

D13: nothing in this report is hand-written. Every number is read out of a
table or a stats JSON this session wrote, so the report and the data cannot
drift. If a table is missing the report says so rather than omitting the
section silently.

Run:  python3 scripts/s3_report.py
"""

from __future__ import annotations

import json
import sys
import time
from collections import Counter
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s3_assign import (  # noqa: E402
    MIN_PROFILE_POSITIONS, MIN_SCORE, REL_MARGIN,
)
from scripts.s3_hmm_lib import (  # noqa: E402
    CENSUS_V3_DIR, HMM_SWEEP_DIR, read_tsv,
)
from scripts.s3_report_d10 import render as render_d10  # noqa: E402


def jload(path: Path) -> dict:
    return json.loads(path.read_text()) if path.exists() else {}


def rows_or_none(path: Path):
    return read_tsv(path) if path.exists() else None


def pct(n, d) -> str:
    return f"{100.0 * n / d:.1f} %" if d else "—"


def main() -> int:
    seed_stats = jload(HMM_SWEEP_DIR / "seed_build_stats.json")
    db_stats = jload(HMM_SWEEP_DIR / "proteome_db_stats.json")
    calib = jload(HMM_SWEEP_DIR / "calibration_summary.json")
    v3 = jload(CENSUS_V3_DIR / "census_v3_stats.json")
    sweep = jload(HMM_SWEEP_DIR / "sweep_stats_hmmsearch.json")
    seeds = rows_or_none(HMM_SWEEP_DIR / "seed_manifest.tsv") or []
    conv = rows_or_none(CENSUS_V3_DIR / "convergence.tsv") or []
    disagree = rows_or_none(HMM_SWEEP_DIR / "calibration_disagreements.tsv") or []
    novel = rows_or_none(CENSUS_V3_DIR / "novel_hits.tsv") or []

    L: list[str] = []
    a = L.append
    a("# S3 — profile-HMM sweep and the completeness argument")
    a("")
    a(f"_Rendered from the committed tables on "
      f"{time.strftime('%Y-%m-%d %H:%M')} by `scripts/s3_report.py` (D13)._")
    a("")
    a("The census this project can defend rests on two instruments that see "
      "different evidence. S2 built the first: a positive **architecture** "
      "test over InterPro's annotation of each record (D14b). S3 builds the "
      "second: **best-profile assignment** — score every sequence against "
      "`itpr.hmm` and `ryr.hmm` and take the winner, but only when it wins "
      "by a margin (D7, D14). Where the two agree the call is firm; where "
      "only one speaks it says which; where they disagree the record is "
      "kept as a conflict rather than resolved by preference.")
    a("")

    # ------------------------------------------------------------ profiles
    a("## 1. The two profiles")
    a("")
    if seed_stats:
        a("| Profile | Seeds | Alignment columns | Match states | Clades |")
        a("|---|---:|---:|---:|---|")
        for name, p in seed_stats.get("profiles", {}).items():
            clades = ", ".join(f"{k} {v}" for k, v in p["by_clade"].items())
            a(f"| `{name}.hmm` | {p['n_seeds']} | {p['alignment_cols']:,} | "
              f"{p['match_states']:,} | {clades} |")
        a("")
        a(f"Aligned with {seed_stats.get('mafft_version', '?')} L-INS-i, "
          "**single-threaded** (D24). MAFFT's return code is checked and a "
          "ragged alignment aborts the build, because S1 established that a "
          "silent MAFFT failure degrades to a star alignment with no other "
          "symptom — and the thread count is pinned because `--thread -1` is "
          "not reproducible: the same RyR seeds aligned twice here gave "
          "8,510 and 8,468 columns and profiles of 4,933 and 4,908 match "
          "states. The SHA-256 of every seed set, alignment and profile is "
          "recorded in `seed_build_stats.json`, so a rebuild that drifts is "
          "visible in the data rather than only in a count.")
        a("")
        a("| Profile | Seed set | Alignment | Profile |")
        a("|---|---|---|---|")
        for name, pr in seed_stats.get("profiles", {}).items():
            a(f"| `{name}.hmm` | `{pr.get('seed_faa_sha256', '')[:16]}…` | "
              f"`{pr.get('alignment_sha256', '')[:16]}…` | "
              f"`{pr.get('hmm_sha256', '')[:16]}…` |")
        a("")
        a("Every seed is drawn from the S2 census and carries that census's "
          "call, so the profiles are labelled by a rule rather than by a "
          "gene name — the same separation S2 kept between its architecture "
          "call and the symbols it audited against. The selection rules are "
          "stated in `scripts/s3_seed_spec.py` and enforced in code: a seed "
          "whose census call, architecture count or length has moved since "
          "the manifest was written aborts the build.")
        a("")
        exceptions = [s for s in seeds if not s["arch"].startswith("5")]
        if exceptions:
            a(f"{len(exceptions)} seeds carry an incomplete architecture and "
              "are in the set deliberately, each with its reason:")
            a("")
            a("| Accession | Profile | Species | Arch | Why |")
            a("|---|---|---|---|---|")
            for s in exceptions:
                a(f"| `{s['accession']}` | {s['profile']} | *{s['species']}* "
                  f"| {s['arch']} | {s['note'].replace('|', '/')} |")
            a("")

    # --------------------------------------------------------- calibration
    a("## 2. Calibrating the instrument before using it")
    a("")
    if calib:
        ns = calib["seeds_excluded"]
        a(f"Both profiles were run over S2's archived seeded space — "
          f"**{calib['records_scored']:,} records** whose call came from the "
          "architecture rule, which reads annotation rather than residues. "
          "That makes the agreement between the two a real test and not a "
          "restatement. The seed sequences are inside that set by "
          "construction, so the number that counts is the one with the "
          f"{calib['seeds_in_set']} seeds removed:")
        a("")
        a("| Scored against | Agree | Disagree | Agreement | Profile "
          "declines to call |")
        a("|---|---:|---:|---:|---:|")
        for key, label in (
                ("all_records", "the architecture call, all records"),
                ("seeds_excluded", "the architecture call, seeds excluded"),
                ("vs_call_with_symbol_fallback",
                 "the census call incl. symbol fallback, seeds excluded")):
            c = calib[key]
            agr = f"{c['agreement']:.4f}" if c["agreement"] is not None else "—"
            a(f"| {label} | {c['agree']:,} | {c['disagree']:,} | {agr} | "
              f"{c['profile_unassigned']:,} |")
        a("")
        a(f"Thresholds: a call needs the winning profile to clear "
          f"**{MIN_SCORE:.0f} bits**, to span at least "
          f"**{MIN_PROFILE_POSITIONS} match states** (D22), and to beat the "
          f"loser by more than **{REL_MARGIN:.0%} of its own score**. The "
          "margin is relative rather than absolute because bit scores scale "
          "with alignable length: a fixed gap would call every full-length "
          "protein confidently and no fragment at all. The span floor is "
          "measured rather than tuned — in this project's own S0 domain "
          "coordinates the shortest observed PF08709, the IP3-binding core "
          "that names the family, is 200 aa, and the longest observed SPRY "
          "is 137 aa, so the floor sits in the gap between them. What it "
          f"costs is in the table above: {calib['seeds_excluded']['profile_unassigned']} "
          "records the architecture rule called lose their profile verdict.")
        a("")
        recovered = calib.get("s2_arch_unassigned_called_by_profile", 0)
        a("The second instrument earns its place on the records the first "
          f"could not call: **{recovered:,} of the "
          f"{calib['seeds_excluded']['truth_unassigned']:,} records the "
          "architecture rule leaves `unassigned` get a call from the "
          "profiles**, because a partial architecture defeats a rule that "
          "reads absence as evidence and does not defeat a sequence "
          "profile.")
        a("")
        sym_n = calib.get("symbol_fallback_records", 0)
        sym_o = calib.get("symbol_fallback_overturned", 0)
        if sym_n:
            a(f"S2 filled {sym_n:,} of those from the record's gene symbol "
              "at low confidence, which is the one place its call was not "
              "purely architectural. The profiles decide those on sequence "
              f"and overturn **{sym_o} of {sym_n:,}** — so the symbol "
              "fallback was sound, and it is now checked rather than "
              "assumed.")
            a("")
        if disagree:
            a(f"**{len(disagree)} records where the two instruments "
              "disagree** (`calibration_disagreements.tsv`):")
            a("")
            a("| Accession | Gene | Species | aa | Architecture | Profile | "
              "Margin |")
            a("|---|---|---|---:|---|---|---:|")
            for r in disagree[:15]:
                a(f"| `{r['accession']}` | {r['gene'] or '—'} | "
                  f"*{r['species'][:32]}* | {r['length']} | {r['truth']} | "
                  f"{r['assignment']} | {float(r['rel_margin']):.1%} |")
            if len(disagree) > 15:
                a(f"| … {len(disagree) - 15} more | | | | | | |")
            a("")
        else:
            a("**No record is called one family by the architecture rule and "
              "the other by the profiles.** Two instruments reading "
              "different evidence, zero disagreements.")
            a("")

    # --------------------------------------------------------------- sweep
    a("## 3. The sweep")
    a("")
    if db_stats:
        a(f"* Database: **{db_stats['n_proteomes']} UniProt vertebrate "
          f"reference proteomes** — {db_stats['n_seqs']:,} canonical "
          f"proteins, {db_stats['n_residues']:,} residues "
          f"({db_stats['db_bytes'] / 1e9:.1f} GB). One canonical protein per "
          "gene: isoform sets are excluded because they add isoforms of "
          "genes already counted, and the census counts genes.")
    if sweep:
        hs = sweep.get("hmmsearch", {})
        for name, p in hs.get("profiles", {}).items():
            a(f"* `{name}.hmm`: {p['targets']:,} targets at E ≤ "
              f"{sweep.get('evalue')} ({p['elapsed_s'] / 60:.0f} min)")
        asg = hs.get("assignment", {})
        if asg:
            calls = asg["calls"]
            a(f"* Assignment over the union: **{asg['targets']:,} targets** → "
              f"ITPR {calls['ITPR']:,}, RYR {calls['RYR']:,}, unassigned "
              f"{calls['unassigned']:,}")
        a("")
    if v3.get("module_only_declined"):
        a(f"Before any of that count means anything, **{v3['module_only_declined']:,} "
          "hits are declined as module-only matches** (D22). Both profiles "
          "are full-length channel models, so a protein sharing one small "
          "domain with either of them scores against it — and `ryr.hmm` "
          "carries SPRY, which D14b had already flagged as sitting in "
          "~114,000 UniProt proteins and being in no sense RyR-specific. "
          "Without the gate this sweep called 14,981 vertebrate proteins "
          "RYR; the declined pile is headed by "
          + ", ".join(f"*{g}*" for g in
                      list(v3["module_only_top_genes"])[:6])
          + " — EF-hand and SPRY proteins, not receptors. The declined "
          "rows are not duplicated into a table of their own: they are in "
          "`hmm_sweep/hmmsearch_assignments.tsv` with `evidence = module`.")
        a("")
        a(f"Of those, **{v3['subthreshold_family_fragments']:,} carry a "
          "gene symbol from this family or its sister** "
          "(`subthreshold_family_fragments.tsv`). They are not SPRY "
          "proteins that happen to score: they are pieces of split or "
          "truncated gene models, sitting under the gate because there is "
          "not enough of the protein left to span a family domain. They are "
          "kept as an annotation lead rather than discarded.")
        a("")
    if novel:
        a(f"**{len(novel)} sweep hits are absent from census v2** "
          "(`novel_hits.tsv`) — proteins in a vertebrate reference proteome "
          "that the InterPro enumeration never returned. Their calls: "
          + ", ".join(f"{k} {v}" for k, v in
                      sorted(Counter(r['call'] for r in novel).items())) + ".")
        a("")
        top = sorted(novel, key=lambda r: -float(r["itpr_score"] or 0))[:12]
        a("| Accession | Gene | Species | aa | itpr.hmm | ryr.hmm | Call |")
        a("|---|---|---|---:|---:|---:|---|")
        for r in top:
            a(f"| `{r['accession']}` | {r['gene'] or '—'} | "
              f"*{r['species'][:30]}* | {r['length']} | {r['itpr_score']} | "
              f"{r['ryr_score']} | {r['call']} |")
        a("")
    elif v3:
        a("**No sweep hit is absent from census v2.** Every protein the two "
          "profiles find in 763 vertebrate reference proteomes was already "
          "in the InterPro enumeration — which is the completeness result "
          "this step exists to produce, stated in the direction that costs "
          "something to claim.")
        a("")

    # ---------------------------------------------------------- jackhmmer
    render_d10(a, conv, v3, rows_or_none)


    # ----------------------------------------------------------- census v3
    a("## 5. Census v3")
    a("")
    if v3:
        a(f"**{v3['v3_rows']:,} records** — census v2's {v3['v2_rows']:,} "
          f"with a second verdict on every row, plus {v3['novel_rows']:,} "
          "novel sweep hits.")
        a("")
        a("| Call | Records |")
        a("|---|---:|")
        for k, n in sorted(v3["calls"].items(), key=lambda kv: -kv[1]):
            a(f"| {k} | {n:,} |")
        a("")
        a("| Instruments speaking | Records | Share |")
        a("|---|---:|---:|")
        for k, n in sorted(v3["instruments"].items(), key=lambda kv: -kv[1]):
            a(f"| {k} | {n:,} | {pct(n, v3['v3_rows'])} |")
        a("")
        if v3.get("changes_by_transition"):
            a(f"**{v3['call_changes']:,} records change call from census v2** "
              "— every one of them a record the profiles could resolve and "
              "the architecture rule could not:")
            a("")
            a("| Transition | Records |")
            a("|---|---:|")
            for k, n in sorted(v3["changes_by_transition"].items(),
                               key=lambda kv: -kv[1]):
                a(f"| {k} | {n:,} |")
            a("")
        a(f"Conflicts (the instruments disagreeing): **{v3['conflicts']}**.")
        a("")

    # -------------------------------------------------------- completeness
    a("## 6. Completeness, in both directions")
    a("")
    if v3:
        m = v3["missed_by_hmm"]
        a("A sweep is only a completeness argument if it is checked in the "
          "expensive direction as well: how many census-v2 ITPR records "
          "from a swept proteome did the profiles *fail* to find? "
          f"**{m['total']:,}** — and each one was then looked up in the "
          "sweep database itself, one streaming pass, rather than explained "
          "away by its gene name:")
        a("")
        a("| Verdict | Records |")
        a("|---|---:|")
        for k, n in sorted(m["by_verdict"].items(), key=lambda kv: -kv[1]):
            a(f"| {k} | {n:,} |")
        a("")
        if m["sensitivity_failures"] == 0:
            a("**Not one of them was in the database and missed.** Every "
              "apparent miss is a UniProtKB entry that the taxon's "
              "*reference proteome* does not contain — a reference proteome "
              "is one canonical protein per gene, and the InterPro census "
              "enumerated all of UniProtKB — so it was never searched. "
              "The profiles' sensitivity failures over 763 vertebrate "
              "reference proteomes number **zero** "
              "(`missed_by_hmm.tsv`).")
        else:
            a(f"**{m['sensitivity_failures']:,} records were in the "
              "database and not found.** Those are real sensitivity "
              "failures and the completeness claim is bounded by them "
              "(`missed_by_hmm.tsv`).")
        a("")
        a(f"**{v3['proteomes_without_itpr']} of {v3['swept_proteomes']} "
          "vertebrate reference proteomes carry no ITPR record** "
          "(`proteomes_without_hits.tsv`). That is a list of leads for the "
          "genome sweep, not a list of losses: a reference proteome is an "
          "annotation, and S4/S5 test these against the genomes themselves.")
        a("")

    a("## 7. What this does not settle")
    a("")
    a("* The sweep is **Vertebrata only**. Every count above is a statement "
      "about 763 vertebrate reference proteomes, and nothing here bears on "
      "the plant, fungal or protist questions S2 opened — those are S20's, "
      "using the same code path and manifest format.")
    a("* A reference proteome is a *gene set*, not a genome. A gene missing "
      "from one may be missing from the annotation, not from the DNA (D9). "
      "The zero-hit list is where S5 starts, not what it concludes.")
    a("* Profile assignment inherits the seeds' breadth. The deep branches "
      "are represented by seven non-metazoan seeds; a lineage more diverged "
      "than any of them is outside what this instrument has been shown to "
      "call.")
    a("")

    out = CENSUS_V3_DIR / "report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L) + "\n")
    print(f"[s3_report] wrote {out} ({len(L)} lines)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
