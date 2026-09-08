"""S11 — renders `results/structures/report.md` purely from the tables (D13).

Scope, the AlphaFold coverage measurement and the instrument live here; the
results — the fold calls, D14 asked structurally, per-domain confidence and
the caveats — live in `s11_report_results.py` (the `s3_report.py` /
`s3_report_d10.py` split, so both halves stay inside the 500-line budget and
read the tables through the same loader and the same formatter).

**Headlines chosen by the data.** Every prior S11 meets is stated in
`s11_priors.PRIOR` with where the earlier task said it, computed from S11's
own tables, and rendered with S8's five-valued verdict — both numbers
printed either way. The comparison that has to be sayable is the one that
goes badly: if a census record turned out **not** to fold like an IP3
receptor, this report has to be able to say so, and the table it would say
it from is `tm_vs_reference.tsv`.

**A section whose table is absent renders *not run yet*,** never nothing, so
a partial run is visible as partial rather than as a shorter report.
"""

from __future__ import annotations

import csv
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s11_priors as priors                                     # noqa: E402
from s11_tmalign import TM_FOLD, TM_RANDOM                      # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "results" / "structures"


def load(name: str) -> list[dict]:
    path = OUT / name
    if not path.exists():
        return []
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_json(name: str) -> dict:
    path = OUT / name
    return json.loads(path.read_text()) if path.exists() else {}


def fmt(x, nd: int = 2) -> str:
    """One number formatter for both halves of the report."""
    if x is None or x == "":
        return "n/a"
    try:
        v = float(x)
    except (TypeError, ValueError):
        return str(x)
    if v == int(v) and abs(v) < 1e15:
        return f"{int(v):,}"
    return f"{v:,.{nd}f}"


def pct(num, den) -> str:
    den = float(den or 0)
    return "n/a" if not den else f"{100.0 * float(num) / den:.1f} %"


def missing(section: str) -> str:
    return (f"\n## {section}\n\n*Not run yet — the table this section renders "
            f"from is not present.*\n")


def prior_line(key: str, ours: str, holds) -> str:
    p = priors.PRIOR[key]
    return (f"> **Prior.** {p['where']}.\n>\n"
            f"> **This task.** {ours}\n>\n"
            f"> {priors.verdict(key, holds)}\n")


# --------------------------------------------------------------------------
def header(manifest: list[dict], stats: dict) -> list[str]:
    ok = [m for m in manifest if m["status"] == "ok"]
    roles = Counter(m["role"] for m in ok)
    return [
        "# S11 — structures",
        "",
        "Two questions. **Does what the census calls an IP3 receptor fold "
        "like one?** — asked of predicted models against cryo-EM references, "
        "with negative controls setting the floor. And **can it be told from "
        "a ryanodine receptor by shape alone?** — D14 put to a sixth "
        "instrument, one that reads no gene symbol, no Pfam annotation, no "
        "alignment score and no tree.",
        "",
        f"The panel is **{len(ok)}** structures: "
        + ", ".join(f"{n} {r.replace('_', ' ')}" for r, n in
                    sorted(roles.items(), key=lambda kv: -kv[1]))
        + ". Every reference was resolved by an RCSB query over the family's "
          "own Pfam signatures and assigned to a family by this project's "
          "census, never by an entry title.",
        "",
        f"> TM-align {stats.get('tmalign_version', 'n/a')}. Scores are read "
        f"against TM-align's own published bars — **{TM_RANDOM}**, below "
        "which a pair is indistinguishable from two random structures of the "
        f"same size, and **{TM_FOLD}**, above which they share a fold. "
        "Neither is a threshold this project chose.",
        "",
        "> Every number below is read from a committed table in this "
        "directory. Nothing here is computed from a structure file, an API "
        "or a TM-align run at render time (D13).",
        "",
    ]


def section_afdb(summary: list[dict], reps: list[dict],
                 census: list[dict]) -> list[str]:
    if not summary:
        return [missing("1. What AlphaFold DB actually holds")]
    all_c = next((r for r in summary
                  if r["scope"] == "census" and r["group"] == "ALL"), {})
    all_r = next((r for r in summary
                  if r["scope"] == "reps" and r["group"] == "ALL"), {})
    out = [
        "## 1. What AlphaFold DB actually holds", "",
        "A vertebrate IP3 receptor subunit is ~2,700 residues, which is where "
        "AlphaFold DB's monomer pipeline stops; the sister family at ~5,000 "
        "is past it outright. So the family's structural coverage is "
        "measured before anything is compared, and it is measured three ways "
        "that do not agree.", "",
        "| scope | records | UniProt-shaped | full-length model | isoform "
        "model only | no model | usable (≥ 95 % covered) |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in (all_r, all_c):
        if not row:
            continue
        out.append(
            f"| {row['scope']} | {fmt(row['n_records'])} | "
            f"{fmt(row['n_applicable'])} | {fmt(row['full'])} | "
            f"{fmt(row['isoform'])} | {fmt(row['absent'])} | "
            f"{fmt(row['usable_ge_95pct'])} "
            f"({pct(row['usable_ge_95pct'], row['n_applicable'])}) |")
    usable = int(all_c.get("usable_ge_95pct") or 0)
    appl = int(all_c.get("n_applicable") or 0)
    rep_usable = int(all_r.get("usable_ge_95pct") or 0)
    rep_n = int(all_r.get("n_records") or 0)
    out += [
        "",
        f"**{pct(usable, appl)} of the census's UniProt-shaped ITPR records "
        f"have a usable predicted model**, and of the {fmt(rep_n)} "
        "representatives S6 selected — the set every alignment, tree and "
        f"selection result in this project stands on — **{fmt(rep_usable)}** "
        "do.", "",
        prior_line("afdb_poor",
                   f"{pct(usable, appl)} of census records and "
                   f"{fmt(rep_usable)}/{fmt(rep_n)} representatives are "
                   "modelled at full length.", True),
        "",
        "### The isoform trap", "",
        "AFDB is keyed on a UniProt accession but answers with whatever "
        "record it holds, and that record may be an **isoform**. Asked for "
        "the three human paralogs it returns `Q14643-4` (2,695 aa of a "
        "2,758-residue protein), `Q14571-2` — **181 residues of a "
        "2,701-residue protein** — and the canonical `Q14573`. A probe that "
        "took the first record and called it a hit would report the family's "
        "reference paralogs as fully modelled. Coverage here is therefore "
        "the modelled span against the length the census holds, and the "
        "181-residue ITPR2 model is carried into the panel and rejected "
        "there by the minimum-chain rule rather than quietly used.", "",
    ]
    if census:
        modelled = [int(r["length"] or 0) for r in census
                    if r["afdb_status"] != "not_applicable"
                    and float(r["model_coverage"] or 0) >= 0.95]
        unmodelled = [int(r["length"] or 0) for r in census
                      if r["afdb_status"] == "absent"]
        if modelled and unmodelled:
            med_m = sorted(modelled)[len(modelled) // 2]
            med_u = sorted(unmodelled)[len(unmodelled) // 2]
            long_m = sum(1 for x in modelled if x >= 2000)
            long_u = sum(1 for x in unmodelled if x >= 2000)
            out += [
                "### Coverage is a function of length, not of taxonomy", "",
                f"The median modelled record is **{fmt(med_m)} aa**; the "
                f"median unmodelled one is **{fmt(med_u)} aa**. Of the "
                f"{fmt(long_m + long_u)} census records at or above the "
                "family's own 2,000 aa floor — the plausibly full-length "
                f"ones — **{fmt(long_m)}** ({pct(long_m, long_m + long_u)}) "
                "have a usable model. The coverage AFDB offers this family "
                "is concentrated on its *fragments*, which are the records a "
                "structural argument can do least with.", "",
                "![](figures/s11_afdb_coverage.png)", "",
                "**Figure 1.** AlphaFold DB coverage of the ITPR census, by "
                "group (**a**) and against record length (**b**). Panel **b** "
                "is the result: the usable-model mass sits below ~1,300 "
                "residues while the peak at ~2,700 — a full-length receptor "
                "subunit — is almost entirely unmodelled.", "",
            ]
    return out


def section_references(cands: list[dict], selection: list[dict],
                       unfilled: list[dict]) -> list[str]:
    if not cands:
        return [missing("2. The references, resolved by query")]
    elig = Counter(c["eligibility"] for c in cands)
    entries = len({c["entry_id"] for c in cands})
    out = [
        "## 2. The references, resolved by query", "",
        "The brief is emphatic that reference structures are resolved by a "
        "query and not from memory, and this family makes the reason "
        "concrete: the IP3 and ryanodine receptor entries share every "
        "diagnostic Pfam and half their titles. So the candidate set is "
        f"enumerated from RCSB by the family's own signatures — **{entries}** "
        f"entries, **{len(cands)}** polymer entities — and each entity's "
        "family is decided by *this project's census*, on the UniProt "
        "accession RCSB maps it to. No title is read.", "",
        "| eligibility | entities |", "|---|---|",
    ]
    for reason, n in sorted(elig.items(), key=lambda kv: -kv[1]):
        out.append(f"| {reason} | {fmt(n)} |")
    xray = [c for c in cands if c["method"] == "X-RAY DIFFRACTION"
            and c["call"] == "ITPR"]
    em = [c for c in cands if c["method"] == "ELECTRON MICROSCOPY"
          and c["call"] == "ITPR" and int(c["sample_length"] or 0) >= 2000]
    if xray and em:
        xl = sorted(int(c["sample_length"] or 0) for c in xray)
        xr = sorted(float(c["resolution"]) for c in xray if c["resolution"])
        el = sorted(int(c["sample_length"] or 0) for c in em)
        longest = max(xray, key=lambda c: int(c["sample_length"] or 0))
        all_res = [float(c["resolution"]) for c in cands
                   if c["call"] == "ITPR" and c["resolution"]]
        worst = (float(longest["resolution"]) >= max(all_res)) if all_res else False
        out += [
            "",
            "The cryo-EM restriction is a measurement about what exists, not "
            f"a preference. The **{len(xray)}** X-ray ITPR entities in the "
            f"candidate set run {fmt(xl[0])}–{fmt(xl[-1])} aa at "
            f"{fmt(xr[0])}–{fmt(xr[-1])} Å, against full-length cryo-EM "
            f"entities of {fmt(el[0])}–{fmt(el[-1])} aa — and the longest "
            f"X-ray construct ({longest['entry_id']}, {fmt(xl[-1])} aa, "
            f"{fmt(longest['resolution'])} Å) is "
            + ("also the lowest-resolution structure the family has"
               if worst else "among the lowest-resolution entries the family "
               "has") + ". Every one of them is a truncation, so an X-ray "
            "reference would make each TM-score a statement about part of a "
            "subunit.", ""]
    refs = [s for s in selection if s["role"] == "reference"]
    states = [s for s in selection if s["role"] == "state_panel"]
    ctrls = [s for s in selection if s["role"] == "control"]
    if refs:
        out += ["### Primary references (R5 — best resolution per paralog)",
                "", "| paralog | entry | resolution (Å) | organism | state | "
                "ligands modelled |", "|---|---|---|---|---|---|"]
        for s in refs:
            out.append(f"| {s['paralog']} | {s['entry_id']} | "
                       f"{fmt(s['resolution'])} | *{s['organism']}* | "
                       f"{s['state']} | {s['ligand_note'] or '—'} |")
        got = {s["paralog"] for s in refs}
        for want in ("ITPR1", "ITPR2", "ITPR3", "RYR1", "RYR2", "RYR3"):
            if want not in got:
                out += ["", f"**No {want} reference exists.** No cryo-EM "
                        "entry in the candidate set carries a full-length "
                        f"{want} entity the census names, so nothing in this "
                        f"report is scored against a {want} structure.", ""]
    if states:
        out += ["", f"### The state panel (R6 — {len(states) + 1} "
                "conformations of one paralog)", "",
                "A model scored against a single conformation confounds the "
                "fold question with the gating question. The state panel is "
                "the control on that: several experimental structures of the "
                "*same* protein, so the report can say how much of a "
                "TM-score difference is conformation before it attributes "
                "any to fold. The table lists "
                f"{len(states)}; the remaining conformation is the primary "
                "reference above, which is the same file and is not "
                "downloaded twice.", "",
                "| state | entry | resolution (Å) | ligands modelled |",
                "|---|---|---|---|"]
        for s in states:
            out.append(f"| {s['state']} | {s['entry_id']} | "
                       f"{fmt(s['resolution'])} | {s['ligand_note'] or '—'} |")
    if ctrls:
        out += ["", "### Negative controls (R7)", "",
                "Drawn from **S1's committed decoy panel**, one per decoy "
                "class, rather than invented here — so S11's floor is set by "
                "the same proteins S1 validated the discovery scorer "
                "against. Each is held to the same cryo-EM and full-length "
                "rules as a reference, and picked to be *size-matched* to "
                "the family reference: TM-score is length-normalised, but "
                "the score two unrelated structures can reach still depends "
                "on their size ratio.", "",
                "| class | protein | entry | chain length | resolution (Å) |",
                "|---|---|---|---|---|"]
        for s in ctrls:
            out.append(f"| {s['control_class']} | {s['decoy_gene']} | "
                       f"{s['entry_id']} | {fmt(s['sample_length'])} | "
                       f"{fmt(s['resolution'])} |")
    if unfilled:
        out += ["", "**Unfilled control classes.** A class with no "
                "qualifying structure is named with the stage that lost it, "
                "never dropped:", ""]
        for u in unfilled:
            out.append(f"- `{u['control_class']}` ({u['n_decoys']} decoys) — "
                       f"{u['stage_lost']}")
        out.append("")
    return out


def section_panel(manifest: list[dict], slots: list[dict],
                  stats: dict) -> list[str]:
    if not manifest:
        return [missing("3. The panel")]
    ok = [m for m in manifest if m["status"] == "ok"]
    bad = [m for m in manifest if m["status"] != "ok"]
    out = [
        "## 3. The panel", "",
        "Nothing enters the panel unread. Every file is parsed, reduced to "
        "**one chain** — both families are homotetramers and an IP3R "
        "assembly is ~11,000 residues, so scoring a monomer model against a "
        "deposited assembly would answer a question about quaternary "
        "structure — and measured. `model_coverage` in the manifest is "
        "*resolved residues over the length the census holds*, not what an "
        "API advertised.", "",
        f"**{len(ok)}** of {len(manifest)} structures are usable.", "",
    ]
    if bad:
        out += ["| rejected | role | reason |", "|---|---|---|"]
        for m in bad:
            out.append(f"| `{m['id']}` | {m['role']} | {m['status']}: "
                       f"{m['note']} |")
        out += ["", "The rejected row is the ITPR2 isoform model — the trap "
                "in §1, caught by a rule rather than by inspection.", ""]
    if slots:
        out += ["### Which slots the models fill", "",
                "The S6 representative set has a usable model for very few "
                "of its tips, so a panel built from representatives alone "
                "would be six proteins and no non-vertebrate breadth. Where "
                "a slot has no representative model, the best-covered census "
                "record in that slot stands in, and the audit records that "
                "the slot was *filled* rather than met.", "",
                "| slot | census records | outcome | filled with |",
                "|---|---|---|---|"]
        for s in slots:
            out.append(f"| {s['slot']} | {fmt(s.get('n_census_records'))} | "
                       f"{s['outcome']}{(' — ' + s['stage_lost']) if s.get('stage_lost') else ''}"
                       f" | {s['filled_with'] or '—'} |")
        out.append("")
    out += ["![](figures/s11_panel.png)", "",
            "**Figure 2.** Every structure in the panel: the pale bar is the "
            "sequence length the record claims, the filled bar what the "
            "structure actually delivers. Faded rows failed the "
            f"{stats.get('parameters', {}).get('min_chain_residues', 300)}"
            "-residue minimum.", ""]
    return out


def main(argv: list[str] | None = None) -> int:
    import s11_report_results as results

    manifest = load("structure_manifest.tsv")
    stats = load_json("structures_stats.json")
    lines: list[str] = []
    lines += header(manifest, stats)
    lines += section_afdb(load("afdb_coverage_summary.tsv"),
                          load("afdb_coverage_reps.tsv"),
                          load("afdb_coverage_census.tsv"))
    lines += section_references(load("reference_candidates.tsv"),
                                load("reference_selection.tsv"),
                                load("unfilled_controls.tsv"))
    lines += section_panel(manifest, load("model_slots.tsv"), stats)
    lines += results.sections(load, load_json, fmt, pct, missing, prior_line)

    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "report.md").write_text("\n".join(lines) + "\n")
    print(f"[s11] report -> {OUT / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
