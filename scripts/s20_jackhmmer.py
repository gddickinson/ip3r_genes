"""S20 step 4 — jackhmmer to convergence from a seed native to each group.

S3 ran three jackhmmer seeds against the *vertebrate* database and found the
family core is seed-independent. S20 asks the complementary question: search
each non-vertebrate group iteratively, from a seed that lives in it, and see
whether the range the single-pass profile reports is the range an iterated
model reaches. An iterated search is the only way to make a completeness
claim about a lineage whose members are too diverged for one profile.

**The seed is derived, not chosen.** For each group the seed is the ITPR
call in that group's own sweep with the highest bit score against
`itpr.hmm`, among records whose profile match spans at least half the model
— so it is a full-length family member of that group rather than the longest
protein or a name that looked right. The pick, its score and what it beat
are written to `jackhmmer_seeds_s20.tsv` before anything runs, so a rerun on
a changed sweep is visible as a changed seed rather than a silently
different experiment. A group with no qualifying record is reported as
having none; that is a result about the group, not a failure.

**D10 is evaluated exactly as in S3** — `s3_kill.evaluate` on the per-round
included lists parsed out of the archived log, with the sister-family call
taken from that group's own assignment table. K1's failure mode (an
ITPR-seeded model drifting across the shared architecture into the RyRs) is
if anything more likely out here, where the two families' members are more
diverged from the profiles than vertebrate ones are.

Run:  python3 scripts/s20_jackhmmer.py --group fungi --cpu 8
      python3 scripts/s20_jackhmmer.py --all --parse-only
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from scripts.s20_evidence import gather_sequences  # noqa: E402
from scripts.s20_groups import EUKARYOTE_GROUPS, EVALUE_PRIMARY, GROUPS  # noqa: E402
from scripts.s20_lib import (  # noqa: E402
    S20_DIR, group_paths, live, log, read_tsv, s20_dirs, write_json, write_tsv,
)
from scripts.s3_hmm_lib import acc_key, parse_jackhmmer_log, write_fasta  # noqa: E402
from scripts.s3_kill import accepted_targets, evaluate  # noqa: E402

MAX_ITER = 10                 # K3's ceiling, S3's value
MIN_SEED_COVERAGE = 0.50      # the seed must span half the model

SEED_FIELDS = ["group", "accession", "species", "taxon_id", "gene", "length",
               "itpr_score", "itpr_coverage", "rel_margin", "n_qualifying",
               "runner_up", "runner_up_score", "rule"]
CONV_FIELDS = ["group", "seed_tag", "accession", "round", "new_targets",
               "n_included", "n_own", "n_sister", "n_offfamily",
               "sister_frac", "sister_rise", "growth", "converged",
               "verdict", "kill_rule"]
COMP_FIELDS = ["group", "clade", "profile_call", "targets",
               "found_by_single_pass", "iteration_only"]


def model_composition(group: str, targets: set[str], assign: dict,
                      taxa: dict) -> list[dict]:
    """What the accepted model is actually built from, by lineage.

    This is the completeness statement the iterated search exists to make:
    a model seeded in one lineage and iterated to convergence either reaches
    a neighbouring lineage or it does not. Targets the single profile pass
    never reported are counted separately, because those are exactly the
    records iteration was run to find — and if they all turn out to sit in
    the same lineage as the seed, the absence next door is not a
    sensitivity artefact.

    A target absent from the sweep's assignment table has no taxid there, so
    it is resolved from the group FASTA's own header — the accession is in
    the database that was searched, by construction.
    """
    from scripts.s20_evidence import extract
    from scripts.s3_hmm_lib import parse_uniprot_header

    def acc_of(name: str) -> str:
        parts = name.split("|")
        return acc_key(parts[1] if len(parts) >= 3 else name)

    by_acc = {acc_of(n): n for n in targets}
    unresolved = {a for a in by_acc if a not in assign}
    extra: dict[str, dict] = {}
    if unresolved:
        for a, (header, _) in extract(group_paths(group)["db"],
                                      unresolved).items():
            extra[a] = parse_uniprot_header(header)
        log("s20_jack", f"{group}: {len(extra)}/{len(unresolved)} "
                        "iteration-only targets resolved from the group FASTA")

    agg: dict[tuple[str, str], list[int]] = {}
    for acc in by_acc:
        row = assign.get(acc)
        taxid = int((row or extra.get(acc, {})).get("taxon_id") or 0)
        clade = taxa.get(taxid, {}).get("clade", "") or "unresolved"
        call = row["assignment"] if row else "not called by the single pass"
        cell = agg.setdefault((clade, call), [0, 0, 0])
        cell[0] += 1
        cell[1] += int(row is not None)
        cell[2] += int(row is None)
    return [{"group": group, "clade": clade, "profile_call": call,
             "targets": v[0], "found_by_single_pass": v[1],
             "iteration_only": v[2]}
            for (clade, call), v in sorted(agg.items(),
                                           key=lambda kv: -kv[1][0])]


def pick_seed(group: str) -> dict | None:
    """The stated rule: highest `itpr.hmm` score among ITPR calls covering
    at least half the model."""
    path = S20_DIR / f"assignments_{group}.tsv"
    if not path.exists():
        return None
    ok = [r for r in read_tsv(path)
          if r["assignment"] == "ITPR"
          and float(r["itpr_coverage"] or 0) >= MIN_SEED_COVERAGE]
    if not ok:
        return None
    ok.sort(key=lambda r: -float(r["itpr_score"] or 0))
    best, runner = ok[0], (ok[1] if len(ok) > 1 else None)
    return {
        "group": group, "accession": best["accession"],
        "species": best["species"], "taxon_id": best["taxon_id"],
        "gene": best["gene"], "length": best["length"],
        "itpr_score": best["itpr_score"],
        "itpr_coverage": best["itpr_coverage"],
        "rel_margin": best["rel_margin"], "n_qualifying": len(ok),
        "runner_up": runner["accession"] if runner else "",
        "runner_up_score": runner["itpr_score"] if runner else "",
        "rule": f"highest itpr.hmm bit score with coverage >= "
                f"{MIN_SEED_COVERAGE}",
    }


def run_cmd(cmd: list[str]) -> float:
    log("s20_jack", "$ " + " ".join(cmd))
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log("s20_jack", f"  rc={proc.returncode} in {dt:.0f}s")
    if proc.returncode != 0:
        sys.stderr.write((proc.stderr or proc.stdout)[-3000:])
        raise SystemExit("jackhmmer failed")
    return dt


def run_group(group: str, seed: dict, cpu: int, raw_dir: Path) -> float:
    """One jackhmmer run: this group's seed against this group's own DB."""
    from scripts.s20_verdicts import fasta_sources
    tag = f"s20_{group}"
    seed_fa = raw_dir / f"jackhmmer_{tag}_seed.faa"
    if not seed_fa.exists():
        seqs = gather_sequences({seed["accession"]}, fasta_sources())
        if seed["accession"] not in seqs:
            raise SystemExit(f"{group}: seed {seed['accession']} has no "
                             "sequence in any archived FASTA")
        header, sequence = seqs[seed["accession"]]
        write_fasta(seed_fa, [(header, sequence)])

    logf = raw_dir / f"jackhmmer_{tag}.log"
    if logf.exists() and logf.stat().st_size:
        log("s20_jack", f"{group}: reusing archived log {logf.name}")
        return 0.0
    tmp_log = logf.with_suffix(".log.part")
    dt = run_cmd(["jackhmmer", "--cpu", str(cpu), "-N", str(MAX_ITER),
                  "-E", EVALUE_PRIMARY, "--incE", EVALUE_PRIMARY, "--noali",
                  "--domtblout", str(raw_dir / f"jackhmmer_{tag}.domtblout"),
                  "-o", str(tmp_log),
                  str(seed_fa), str(group_paths(group)["db"])])
    tmp_log.rename(logf)     # atomic, for the same reason S5b made miniprot
    return dt                # atomic: a truncated log parses as a short run


def parse_group(group: str, seed: dict, raw_dir: Path,
                elapsed_s: float | None) -> tuple[list[dict], dict]:
    """Archived log → per-round rows + the D10 verdict for one group."""
    tag = f"s20_{group}"
    logf = raw_dir / f"jackhmmer_{tag}.log"
    if not logf.exists():
        raise SystemExit(f"no archived log for {group}: {logf}")
    conv = parse_jackhmmer_log(logf)
    rounds = conv["rounds"]

    calls = {acc_key(r["accession"]): r["assignment"]
             for r in read_tsv(S20_DIR / f"assignments_{group}.tsv")}

    def call_of(name: str) -> str:
        parts = name.split("|")
        return calls.get(acc_key(parts[1] if len(parts) >= 3 else name),
                         "unassigned")

    by_name = {t: call_of(t)
               for rd in rounds for t in rd.get("included", [])}
    verdict = evaluate(rounds, by_name, own_family="ITPR",
                       converged=bool(conv["converged"]))
    verdict.update({"group": group, "tag": tag,
                    "accession": seed["accession"],
                    "converged": conv["converged"], "n_rounds": len(rounds),
                    "accepted_targets": len(
                        accepted_targets(rounds, verdict["accepted_rounds"]))})
    if elapsed_s:
        verdict["elapsed_s"] = round(elapsed_s, 1)

    (raw_dir / f"jackhmmer_{tag}_rounds.json").write_text(json.dumps(
        {"group": group, "tag": tag, "accession": seed["accession"],
         "converged": conv["converged"], "rounds": rounds}, indent=1))

    rows = []
    for pr in verdict["per_round"]:
        rd = next((r for r in rounds if r["round"] == pr["round"]), {})
        rows.append({
            "group": group, "seed_tag": tag, "accession": seed["accession"],
            "round": pr["round"], "new_targets": rd.get("new_targets", ""),
            "n_included": pr.get("n_included", ""),
            "n_own": pr.get("n_own", ""), "n_sister": pr.get("n_sister", ""),
            "n_offfamily": pr.get("n_offfamily", ""),
            "sister_frac": pr.get("sister_frac", ""),
            "sister_rise": pr.get("sister_rise", ""),
            "growth": pr.get("growth", ""),
            "converged": conv["converged"],
            "verdict": verdict["verdict"], "kill_rule": verdict.get("rule", ""),
        })
    log("s20_jack", f"{group}: {len(rounds)} rounds, "
                    f"converged={conv['converged']}, "
                    f"verdict={verdict['verdict']} "
                    f"({verdict.get('rule') or 'no rule fired'})")
    return rows, verdict


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--group", action="append", choices=list(GROUPS))
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--cpu", type=int, default=8)
    ap.add_argument("--parse-only", action="store_true")
    args = ap.parse_args()

    groups = list(EUKARYOTE_GROUPS) if args.all else (args.group or [])
    if not groups:
        ap.error("give --group (repeatable) or --all")
    _, raw_dir = s20_dirs()

    seeds, no_seed = [], []
    for g in groups:
        s = pick_seed(g)
        if s is None:
            no_seed.append(g)
            log("s20_jack", f"{g}: no ITPR call covering "
                            f">= {MIN_SEED_COVERAGE} of the model — no seed")
        else:
            seeds.append(s)
            log("s20_jack", f"{g}: seed {s['accession']} "
                            f"({s['species']}, {s['itpr_score']} bits, "
                            f"1 of {s['n_qualifying']} qualifying)")
    write_tsv(S20_DIR / "jackhmmer_seeds_s20.tsv", SEED_FIELDS, seeds)

    from scripts.s20_taxa import load_table
    taxa = load_table()
    rows, verdicts, comp = [], {}, []
    for i, seed in enumerate(seeds):
        g = seed["group"]
        live([("fetch proteomes", True), ("sweep", True),
              ("relaxed panel", True), ("plant/fungal verdicts", True)]
             + [(f"jackhmmer {s['group']}", s["group"] in
                 [x["group"] for x in seeds[:i]]) for s in seeds]
             + [("census v5", False)])
        dt = None if args.parse_only else run_group(g, seed, args.cpu, raw_dir)
        r, v = parse_group(g, seed, raw_dir, dt)
        rows += r
        verdicts[g] = v
        data = json.loads(
            (raw_dir / f"jackhmmer_s20_{g}_rounds.json").read_text())
        assign = {acc_key(a["accession"]): a
                  for a in read_tsv(S20_DIR / f"assignments_{g}.tsv")}
        comp += model_composition(
            g, accepted_targets(data["rounds"], v["accepted_rounds"]),
            assign, taxa)

    write_tsv(S20_DIR / "jackhmmer_convergence_s20.tsv", CONV_FIELDS, rows)
    write_tsv(S20_DIR / "jackhmmer_model_composition_s20.tsv",
              COMP_FIELDS, comp)
    write_json(S20_DIR / "jackhmmer_verdicts_s20.json",
               {"groups_without_seed": no_seed,
                "min_seed_coverage": MIN_SEED_COVERAGE,
                "max_iter": MAX_ITER, "evalue": EVALUE_PRIMARY,
                "verdicts": {g: {k: v for k, v in val.items()
                                 if k != "per_round"}
                             for g, val in verdicts.items()}})
    return 0


if __name__ == "__main__":
    sys.exit(main())
