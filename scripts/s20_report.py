"""S20 — render results/s20_sweep/report.md purely from the committed tables.

D13: nothing here is hand-written. Every number is read out of a table or a
stats JSON this task wrote, so the report and the data cannot drift. A
missing table is said so rather than dropping the section silently.

This half carries the scope and the instrument; `s20_report_results.py`
carries the results and takes this module's loader and formatter, so the two
halves cannot read the same tables differently.

Run:  python3 scripts/s20_report.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_groups import (  # noqa: E402
    ALL_GROUPS, EUKARYOTE_GROUPS, EVALUE_PRIMARY, EVALUE_RELAXED, GROUPS,
    PROKARYOTE_GROUPS, SAMPLE_RULE,
)
from scripts.s20_lib import (  # noqa: E402
    CENSUS_V5_DIR, S20_DIR, group_paths, log, read_json, read_tsv,
)
from scripts.s20_report_results import (  # noqa: E402
    section_census, section_d14, section_fungi, section_jackhmmer,
    section_negatives, section_plants, section_range,
)
from scripts.s3_assign import (  # noqa: E402
    MIN_PROFILE_POSITIONS, MIN_SCORE, REL_MARGIN,
)
from scripts.s3_seed_spec import ITPR_SEEDS, RYR_SEEDS  # noqa: E402


def rows(path: Path) -> list[dict]:
    return read_tsv(path) if path.exists() else []


def table(headers: list[str], body: list[list]) -> str:
    out = ["| " + " | ".join(str(h) for h in headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in body]
    return "\n".join(out) + "\n"


def main() -> int:
    presence = rows(S20_DIR / "proteome_presence.tsv")
    taxa = {int(r["taxid"]): r for r in rows(S20_DIR / "taxonomy.tsv")}
    verdicts = rows(S20_DIR / "plant_fungal_verdicts.tsv")
    relaxed = rows(S20_DIR / "relaxed_hits.tsv")
    conv = rows(S20_DIR / "jackhmmer_convergence_s20.tsv")
    conv_json = read_json(S20_DIR / "jackhmmer_verdicts_s20.json")
    sweep_stats = read_json(S20_DIR / "sweep_stats.json") or {"groups": {}}
    panel = read_json(S20_DIR / "relaxed_panel.json") or {}
    v5 = read_json(CENSUS_V5_DIR / "census_v5_stats.json")
    changes = rows(CENSUS_V5_DIR / "call_changes.tsv")
    conflicts = rows(CENSUS_V5_DIR / "conflicts.tsv")
    by_group = rows(CENSUS_V5_DIR / "calls_by_group.tsv")

    L: list[str] = []
    A = L.append
    A("# S20 — the non-vertebrate sweep: the family's true range")
    A("")
    A(f"_Rendered from the committed tables on "
      f"{time.strftime('%Y-%m-%d %H:%M')} by `scripts/s20_report.py` (D13)._")
    A("")
    A("S3 swept the vertebrate reference proteomes and found where the three "
      "paralogs live. This task asks the opposite question — how far the "
      "family reaches — and the one it was set up to answer: the databases "
      "hold a few tens of plant and fungal records while *Arabidopsis* and "
      "*S. cerevisiae* hold none. Which of those is a fact about genomes and "
      "which about databases?")
    A("")

    # ------------------------------------------------------------- scope
    A("## The declared search space\n")
    total_p = total_s = total_r = 0
    body = []
    for g in ALL_GROUPS:
        st = sweep_stats["groups"].get(g, {})
        db = read_json(group_paths(g)["stats"]) or {}
        n_prot = db.get("n_swept_proteomes", 0)
        n_seq = db.get("n_seqs", 0)
        n_res = db.get("n_residues", 0)
        total_p += n_prot
        total_s += n_seq
        total_r += n_res
        sample = (db.get("sample") or {}).get("rule", "none")
        body.append([g, f"{n_prot:,}", f"{n_seq:,}",
                     f"{n_res / 1e9:.2f} G" if n_res else "—",
                     "every reference proteome" if sample == "none"
                     else sample,
                     "yes" if st else "not swept"])
    A(table(["group", "proteomes", "proteins", "residues", "sampling",
             "swept"], body))
    A(f"\n**{total_p:,} reference proteomes, {total_s:,} proteins, "
      f"{total_r / 1e9:.2f} G residues.** The four eukaryote groups "
      "partition Eukaryota with S3's `vertebrata`, so no proteome is swept "
      "twice and the two sweeps' denominators add.\n")
    A(f"\nBacteria are sampled and archaea are not: "
      f"{SAMPLE_RULE['largest_per_genus']}. Archaea's reference set is small "
      "enough to take whole, so no sampling caveat attaches to it.\n")

    unavailable = []
    for g in ALL_GROUPS:
        unavailable += rows(group_paths(g)["unavailable"])
    A(f"\n{len(unavailable)} proteome(s) UniProt lists are not published in "
      "the current release FTP tree and 404 permanently; they are excluded "
      "from the denominator above and recorded in "
      "`proteome_unavailable_<group>.tsv` with what they took out of it.\n")

    # -------------------------------------------------------- instrument
    A("\n## The instrument\n")
    A(f"The same two profiles S3 built and calibrated — `itpr.hmm` and "
      f"`ryr.hmm` — with the same margin assignment: a target is called only "
      f"when the winning profile beats the loser by more than "
      f"{REL_MARGIN:.0%} of its own score, clears {MIN_SCORE:.0f} bits, and "
      f"spans at least {MIN_PROFILE_POSITIONS} match states (D22). Reusing "
      "S3's instrument rather than building a new one is what makes the "
      "vertebrate and non-vertebrate numbers comparable at all.\n")
    A(f"\n**One search, two sensitivities.** Every `hmmsearch` ran at "
      f"`-E {EVALUE_RELAXED}` and the primary call was taken by filtering "
      f"the same output at E ≤ {EVALUE_PRIMARY}. `-E` is a reporting "
      "threshold and does not touch hmmsearch's acceleration filters, so the "
      "strict set is exactly what a strict run would have produced and the "
      "relaxed set is a superset of it from the same search — not a second "
      "experiment that might have differed some other way.\n")
    if panel.get("pfam"):
        A(f"\nThe relaxed panel adds the family's four Pfam domain models:\n")
        A(table(["Pfam", "name", "match states", "why it is in the panel"],
                [[p["pfam"], p.get("name", ""), p.get("leng", ""), p["why"]]
                 for p in panel["pfam"]]))

    body = []
    for g in ALL_GROUPS:
        st = sweep_stats["groups"].get(g)
        if not st:
            continue
        prof = st["profiles"]
        a = st.get("assignment", {}).get("calls", {})
        body.append([g,
                     f"{prof['itpr']['targets_primary']:,}",
                     f"{prof['ryr']['targets_primary']:,}",
                     f"{a.get('ITPR', 0):,}", f"{a.get('RYR', 0):,}",
                     f"{a.get('unassigned', 0):,}",
                     f"{(prof['itpr']['elapsed_s'] + prof['ryr']['elapsed_s']) / 60:.0f}"])
    if body:
        A(table(["group", "itpr.hmm targets", "ryr.hmm targets", "→ ITPR",
                 "→ RYR", "→ unassigned", "search min"], body))

    # ----------------------------------------------------------- results
    assignments = [r for g in ALL_GROUPS
                   for r in rows(S20_DIR / f"assignments_{g}.tsv")]
    section_range(A, table, presence, taxa, v5)
    seed_accs = {s[0] for s in ITPR_SEEDS} | {s[0] for s in RYR_SEEDS}
    section_d14(A, table, assignments, MIN_SCORE, REL_MARGIN,
                MIN_PROFILE_POSITIONS, seed_accs)
    section_plants(A, table, presence, taxa, verdicts)
    section_fungi(A, table, presence, taxa, verdicts)
    section_negatives(A, table, relaxed, taxa, assignments,
                      [g for g in ALL_GROUPS if GROUPS[g]["relaxed"]],
                      EVALUE_RELAXED, EVALUE_PRIMARY)
    section_jackhmmer(A, table, conv, conv_json,
                      rows(S20_DIR / "jackhmmer_model_composition_s20.tsv"),
                      rows(S20_DIR / "jackhmmer_iteration_only_s20.tsv"))
    section_census(A, table, v5, changes, conflicts, by_group)

    # ------------------------------------------------------------ figures
    A("\n## Figures\n")
    figs = sorted((S20_DIR / "figures").glob("*.png"))
    if figs:
        for f in figs:
            A(f"![{f.stem}](figures/{f.name})\n")
    else:
        A("_No figures rendered yet — run `scripts/s20_figures.py`._\n")

    A("\n---\n")
    A(f"Groups: {', '.join(EUKARYOTE_GROUPS)} (eukaryotes), "
      f"{', '.join(PROKARYOTE_GROUPS)} (prokaryotes).\n")

    out = S20_DIR / "report.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(L) + "\n")
    log("s20_report", f"{len(L)} blocks → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
