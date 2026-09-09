"""S22 stage 8 — the lineage test with the sequences the question needs.

Stage 7 asks the brief's comparative question on the alignment this project
already had, and answers it with two tips.  That is a statement about
msa_v2's sampling, not about the hypothesis: msa_v2 was built to span four
kingdoms with 134 representatives, and the taxa that carry an IP3 receptor
without a PI-PLC are not the taxa a diversity rule picks.

The sequences exist.  Every one of the 662 non-vertebrate reference
proteomes that carries an ITPR was swept in S20 and its best family record
is named in `proteome_presence.tsv`, so the panel is the whole population
the pathway question is asked of rather than a sample of it.

Three things make the measurement usable at 20-30 % identity.

**The statistic is the paired difference and nothing else.**  A
*Phytophthora* receptor is divergent from the human reference everywhere;
what the hypothesis predicts is that the gap between the ligand core and the
pore widens.  Both numbers come from one pairwise alignment of one sequence,
so a bad alignment damages both alike.

**Both modules must be covered.**  A record resolving the pore and not the
core would enter as an extreme ligand-core divergence, which is what a
fragmentary gene model looks like.  The bar is the same one stage 7 uses and
the dropped records name the module that lost them.

**The comparison is made inside a clade.**  PI-PLC absence is not scattered
across the eukaryotes at random — it is concentrated in the oomycetes and
the early-diverging fungi — so a pooled test of absent against present is
partly a test of oomycete against everything else.  The stratified test asks
the question *inside* each phylum where both cells exist, which is the only
form of it that is about the pathway rather than about the clade.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(PROJECT_ROOT / "scripts"))

import s22_lib as L  # noqa: E402
import s22_modules as M  # noqa: E402
import s17_lib as S17  # noqa: E402
from src.utils.data_root import require_data_root  # noqa: E402

MIN_MODULE_COVERAGE = 0.50
MIN_CELL = 3            # a phylum needs this many records in *both* cells


def proteome_file(group: str, upid: str, taxid: str) -> Path | None:
    d = require_data_root() / "proteomes" / group
    p = d / f"{upid}_{taxid}.fasta.gz"
    if p.exists():
        return p
    hits = sorted(d.glob(f"{upid}_*.fasta.gz"))
    return hits[0] if hits else None


def fetch_sequence(path: Path, accession: str) -> str | None:
    """One record out of one proteome file, by accession."""
    want = f"|{accession}|"
    seq: list[str] = []
    keep = False
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if keep:
                    break
                keep = want in line or line[1:].split()[0] == accession
            elif keep:
                seq.append(line.strip())
    return "".join(seq) or None


def taxonomy() -> dict[str, dict]:
    return {r["taxid"]: r for r in L.read_tsv(L.S20_DIR / "taxonomy.tsv")}


def module_residues() -> tuple[set[int], set[int]]:
    mods = L.read_tsv(L.OUT_DIR / "module_map.tsv")
    return (M.residues(mods, L.STRUCTURE_REF, M.PRIMARY["ligand_core"]),
            M.residues(mods, L.STRUCTURE_REF, M.PRIMARY["pore_module"]))


def measure(ref: str, seq: str, core: set[int], pore: set[int]) -> dict:
    """Core and pore identity of one record against the human reference."""
    a, b = S17.mafft_pair(ref, seq)
    out: dict = {}
    resi = 0
    hit = {"core": [0, 0, 0], "pore": [0, 0, 0]}   # n_ref, covered, identical
    for ca, cb in zip(a, b):
        if ca == "-":
            continue
        resi += 1
        which = "core" if resi in core else ("pore" if resi in pore else None)
        if which is None:
            continue
        h = hit[which]
        h[0] += 1
        if cb != "-":
            h[1] += 1
            if cb == ca:
                h[2] += 1
    for which in ("core", "pore"):
        n_ref, cov, ident = hit[which]
        out[f"{which}_n_ref"] = n_ref
        out[f"{which}_n_covered"] = cov
        out[f"{which}_coverage"] = cov / n_ref if n_ref else 0.0
        out[f"{which}_identity"] = ident / cov if cov else None
    return out


def build_panel() -> list[dict]:
    presence = {r["upid"]: r for r in
                L.read_tsv(L.S20_DIR / "proteome_presence.tsv")}
    plc = {r["upid"]: r for r in L.read_tsv(L.OUT_DIR / "plc_repertoire.tsv")}
    tax = taxonomy()
    core, pore = module_residues()
    ref = S17.uniprot_fasta(L.REFERENCES[L.STRUCTURE_REF][1])
    targets = [r for r in presence.values() if r["itpr_status"] == "present"]
    L.log(f"deep lineage panel: {len(targets)} ITPR-carrying proteomes")
    rows: list[dict] = []
    for i, r in enumerate(targets, 1):
        if i % 50 == 0:
            L.log(f"  {i}/{len(targets)}")
            L.live("deep_lineage", [(f"aligned {i}/{len(targets)}", False)])
        p = plc.get(r["upid"])
        t = tax.get(r["taxid"], {})
        rec = {"upid": r["upid"], "organism": r["organism"],
               "taxid": r["taxid"], "group": r["group"],
               "kingdom": t.get("kingdom", ""), "phylum": t.get("phylum", ""),
               "class": t.get("class", ""), "order": t.get("order", ""),
               "n_itpr": r["n_itpr"], "accession": r["best_itpr_acc"],
               "protein_count": r["protein_count"],
               "n_pi_plc": p["n_pi_plc"] if p else "",
               "n_plc_x": p["n_plc_x_only_domain"] if p else "",
               "n_plc_y": p["n_plc_y_only_domain"] if p else "",
               "plc_status": p["plc_status"] if p else "not_searched"}
        path = proteome_file(r["group"], r["upid"], r["taxid"])
        seq = fetch_sequence(path, r["best_itpr_acc"]) if path else None
        if not seq:
            rec.update({"seq_length": 0, "in_test": False,
                        "drop_reason": "sequence not recovered"})
            rows.append(rec)
            continue
        rec["seq_length"] = len(seq)
        rec.update({k: L.fmt(v) if isinstance(v, float) else v
                    for k, v in measure(ref, seq, core, pore).items()})
        ok = (rec["core_identity"] not in ("", None)
              and rec["pore_identity"] not in ("", None)
              and float(rec["core_coverage"]) >= MIN_MODULE_COVERAGE
              and float(rec["pore_coverage"]) >= MIN_MODULE_COVERAGE)
        rec["in_test"] = ok
        rec["delta_core_minus_pore"] = (
            L.fmt(float(rec["core_identity"]) - float(rec["pore_identity"]))
            if ok else "")
        if not ok:
            lost = ("core" if rec["core_identity"] in ("", None)
                    or float(rec["core_coverage"]) < MIN_MODULE_COVERAGE
                    else "pore")
            rec["drop_reason"] = f"{lost} coverage below {MIN_MODULE_COVERAGE}"
        else:
            rec["drop_reason"] = ""
        rows.append(rec)
    return rows


def stratified(rows: list[dict]) -> list[dict]:
    """The test inside each phylum where both cells are populated."""
    live = [r for r in rows if r["in_test"] is True or r["in_test"] == "True"]
    by_phylum: dict[str, dict[str, list[float]]] = {}
    for r in live:
        cell = by_phylum.setdefault(r["phylum"] or "unclassified",
                                    {"absent": [], "present": []})
        if r["plc_status"] in cell:
            cell[r["plc_status"]].append(float(r["delta_core_minus_pore"]))
    out: list[dict] = []
    for phylum, cell in sorted(by_phylum.items()):
        a, b = cell["absent"], cell["present"]
        testable = len(a) >= MIN_CELL and len(b) >= MIN_CELL
        mw = L.mann_whitney(a, b) if testable else {"p": None, "cles": None}
        out.append({
            "stratum": phylum, "n_plc_absent": len(a), "n_plc_present": len(b),
            "testable": testable,
            "median_absent": L.fmt(L.median(a)),
            "median_present": L.fmt(L.median(b)),
            "shift": L.fmt((L.median(a) - L.median(b)) if a and b else None),
            "cles": L.fmt(mw["cles"]), "p_mannwhitney": L.fmt(mw["p"], 6),
            "note": "" if testable else
                    f"needs {MIN_CELL} records in both cells",
        })
    qs = L.benjamini_hochberg([float(r["p_mannwhitney"])
                               if r["p_mannwhitney"] else None for r in out])
    for r, q in zip(out, qs):
        r["q_mannwhitney"] = L.fmt(q, 6)
    # the pooled test, reported with its confound named
    a = [float(r["delta_core_minus_pore"]) for r in live
         if r["plc_status"] == "absent"]
    b = [float(r["delta_core_minus_pore"]) for r in live
         if r["plc_status"] == "present"]
    mw = L.mann_whitney(a, b)
    out.append({
        "stratum": "POOLED (all phyla)", "n_plc_absent": len(a),
        "n_plc_present": len(b), "testable": bool(a and b),
        "median_absent": L.fmt(L.median(a)),
        "median_present": L.fmt(L.median(b)),
        "shift": L.fmt((L.median(a) - L.median(b)) if a and b else None),
        "cles": L.fmt(mw["cles"]), "p_mannwhitney": L.fmt(mw["p"], 6),
        "q_mannwhitney": "", "note": "confounded by clade — PI-PLC absence "
                                     "is concentrated in two phyla",
    })
    return out


def covariates(rows: list[dict]) -> list[dict]:
    """What else moves the paired difference, measured rather than dismissed.

    The paired statistic controls for how divergent a tip is *only if* the
    two modules degrade at the same rate with distance.  They do not have
    to: the pore is transmembrane and structurally constrained, so at long
    evolutionary distance it can stay alignable while the core drifts, and
    the paired difference would then be a measure of divergence wearing the
    hypothesis's clothes.  This is S15a's confounder table (D47) applied to
    the same kind of claim.
    """
    from s15_lib import spearman
    live = [r for r in rows if r["in_test"] is True or r["in_test"] == "True"]
    d = [float(r["delta_core_minus_pore"]) for r in live]
    out = []
    for name, get in (("pore_identity", lambda r: float(r["pore_identity"])),
                      ("core_identity", lambda r: float(r["core_identity"])),
                      ("protein_count", lambda r: float(r["protein_count"] or 0)),
                      ("n_itpr", lambda r: float(r["n_itpr"] or 0)),
                      ("seq_length", lambda r: float(r["seq_length"] or 0))):
        xs = [get(r) for r in live]
        rho, pv, n = spearman(xs, d)
        ab = [get(r) for r in live if r["plc_status"] == "absent"]
        pr = [get(r) for r in live if r["plc_status"] == "present"]
        out.append({"covariate": name, "n": n, "rho_vs_delta": L.fmt(rho),
                    "p": L.fmt(pv, 6),
                    "median_plc_absent": L.fmt(L.median(ab)),
                    "median_plc_present": L.fmt(L.median(pr))})
    return out


#: Two records are divergence-matched when their pore identity — a proxy for
#: how far the whole protein has moved from the reference — differs by less
#: than this.  The bar is deliberately tight: a loose one would re-import the
#: confound it exists to remove.
MATCH_TOLERANCE = 0.03
MAX_MATCHES = 3


def divergence_matched(rows: list[dict]) -> tuple[list[dict], dict]:
    """The test again, with each PLC-absent record matched on divergence.

    PLC-absent proteomes sit at a median pore identity far below the
    PLC-present ones, so the unmatched comparison is partly a comparison of
    distant against near.  Each absent record is matched to the nearest
    PLC-present records within `MATCH_TOLERANCE` of its own pore identity and
    scored against their median: a within-pair difference, sign-tested, with
    the records that find no match reported rather than dropped silently.
    """
    live = [r for r in rows if r["in_test"] is True or r["in_test"] == "True"]
    pres = [r for r in live if r["plc_status"] == "present"]
    pairs: list[dict] = []
    for r in (x for x in live if x["plc_status"] == "absent"):
        pi = float(r["pore_identity"])
        near = sorted(
            ((abs(float(c["pore_identity"]) - pi), c["upid"], c) for c in pres
             if abs(float(c["pore_identity"]) - pi) <= MATCH_TOLERANCE),
            key=lambda t: (t[0], t[1]))
        chosen = [c for _, _, c in near[:MAX_MATCHES]]
        rec = {"upid": r["upid"], "organism": r["organism"],
               "phylum": r["phylum"], "pore_identity": r["pore_identity"],
               "delta": r["delta_core_minus_pore"],
               "n_matches": len(chosen),
               "matched_upids": ";".join(c["upid"] for c in chosen),
               "matched_median_pore_identity":
                   L.fmt(L.median([float(c["pore_identity"]) for c in chosen])),
               "matched_median_delta":
                   L.fmt(L.median([float(c["delta_core_minus_pore"])
                                   for c in chosen])),
               "difference": "", "note": ""}
        if chosen:
            rec["difference"] = L.fmt(
                float(r["delta_core_minus_pore"])
                - L.median([float(c["delta_core_minus_pore"]) for c in chosen]))
        else:
            rec["note"] = ("no PLC-present record within "
                           f"{MATCH_TOLERANCE} pore identity")
        pairs.append(rec)
    diffs = [float(p["difference"]) for p in pairs if p["difference"] != ""]
    st = L.sign_test(diffs)
    w = L.wilcoxon(diffs)
    lo, hi = L.bootstrap_ci(diffs)
    summary = {
        "n_plc_absent_in_test": len(pairs),
        "n_matched": len(diffs),
        "n_unmatched": len(pairs) - len(diffs),
        "median_difference": L.fmt(L.median(diffs)),
        "mean_difference": L.fmt(L.mean(diffs)),
        "ci95_lo": L.fmt(lo), "ci95_hi": L.fmt(hi),
        "n_core_more_relaxed": st["n_neg"], "n_core_less_relaxed": st["n_pos"],
        "n_ties": st["n_ties"],
        "p_sign": L.fmt(st["p"], 6), "p_wilcoxon": L.fmt(w["p"], 6),
        "match_tolerance": MATCH_TOLERANCE, "max_matches": MAX_MATCHES,
    }
    return pairs, summary


def ryr_control(rows: list[dict]) -> list[dict]:
    """The positive control, measured through *this* instrument.

    Stage 7 scores the ryanodine receptors on msa_v2's four-kingdom
    alignment; the matched test is a pairwise alignment to the human
    reference, and a control measured on a different instrument cannot
    calibrate this one.  So the same six RyR sequences are put through the
    same pairwise measurement — same reference, same module residues, same
    coverage bar — and their paired difference is what the null in
    `deep_lineage_matched_summary.tsv` has to be read against.
    """
    core, pore = module_residues()
    ref = S17.uniprot_fasta(L.REFERENCES[L.STRUCTURE_REF][1])
    _aln, reps = L.msa_v2()
    out: list[dict] = []
    for r in reps:
        if r["group"] != "RYR":
            continue
        seq = S17.uniprot_fasta(r["accession"])
        rec = {"tip": r["label"], "accession": r["accession"],
               "species": r["species"], "paralog": r["paralog"],
               "seq_length": len(seq)}
        m = measure(ref, seq, core, pore)
        rec.update({k: L.fmt(v) if isinstance(v, float) else v
                    for k, v in m.items()})
        ok = (m["core_identity"] is not None and m["pore_identity"] is not None
              and m["core_coverage"] >= MIN_MODULE_COVERAGE
              and m["pore_coverage"] >= MIN_MODULE_COVERAGE)
        rec["in_test"] = ok
        rec["delta_core_minus_pore"] = (
            L.fmt(m["core_identity"] - m["pore_identity"]) if ok else "")
        out.append(rec)
    itpr = [float(r["delta_core_minus_pore"]) for r in rows
            if (r["in_test"] is True or r["in_test"] == "True")]
    ryr = [float(r["delta_core_minus_pore"]) for r in out if r["in_test"]]
    mw = L.mann_whitney(ryr, itpr)
    for r in out:            # the summary row carries two extra columns, and
        r.setdefault("shift_vs_itpr", "")   # write_tsv takes its header from
        r.setdefault("p_mannwhitney", "")   # the first row (S22 T-series' point)
    out.append({"tip": "SUMMARY", "accession": "", "species": "",
                "paralog": "", "seq_length": "",
                "core_n_ref": "", "core_n_covered": "", "core_coverage": "",
                "core_identity": L.fmt(L.median(
                    [float(r["core_identity"]) for r in out if r["in_test"]])),
                "pore_n_ref": "", "pore_n_covered": "", "pore_coverage": "",
                "pore_identity": L.fmt(L.median(
                    [float(r["pore_identity"]) for r in out if r["in_test"]])),
                "in_test": len(ryr),
                "delta_core_minus_pore": L.fmt(L.median(ryr)),
                "shift_vs_itpr": L.fmt((L.median(ryr) or 0) - (L.median(itpr) or 0)),
                "p_mannwhitney": L.fmt(mw["p"], 6)})
    return out


def matched_power(pairs: list[dict], iters: int = 4000,
                  seed: int = 20220422) -> list[dict]:
    """What shift the *matched* test could have detected, by simulation.

    The differences are resampled from their own observed distribution with
    a shift added and the same two-sided sign test run, so the answer is in
    the units the result is reported in and carries the observed spread
    rather than an assumed one.
    """
    import random
    diffs = [float(p["difference"]) for p in pairs if p["difference"] != ""]
    if len(diffs) < 5:
        return [{"n_pairs": len(diffs), "shift": "", "power": "",
                 "note": "too few matched pairs to simulate"}]
    rng = random.Random(seed)
    out = []
    for shift in (0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.10):
        rej = 0
        for _ in range(iters):
            sample = [rng.choice(diffs) + shift for _ in range(len(diffs))]
            p = L.sign_test(sample)["p"]
            if p is not None and p < 0.05:
                rej += 1
        out.append({"n_pairs": len(diffs), "shift": L.fmt(shift),
                    "power": L.fmt(rej / iters), "note": ""})
    return out


def power(rows: list[dict], iters: int = 3000, seed: int = 20220422) -> list[dict]:
    import random
    live = [r for r in rows if r["in_test"] is True or r["in_test"] == "True"]
    ref = [float(r["delta_core_minus_pore"]) for r in live
           if r["plc_status"] == "present"]
    n_abs = sum(1 for r in live if r["plc_status"] == "absent")
    if n_abs == 0 or len(ref) < 5:
        return [{"n_test_group": n_abs, "n_reference": len(ref), "shift": "",
                 "power": "", "note": "the test could not be run"}]
    rng = random.Random(seed)
    out = []
    for shift in (0.0, 0.01, 0.02, 0.03, 0.05, 0.10, 0.15, 0.20):
        rej = 0
        for _ in range(iters):
            a = [rng.choice(ref) + shift for _ in range(n_abs)]
            b = [rng.choice(ref) for _ in range(len(ref))]
            p = L.mann_whitney(a, b)["p"]
            if p is not None and p < 0.05:
                rej += 1
        out.append({"n_test_group": n_abs, "n_reference": len(ref),
                    "shift": L.fmt(shift), "power": L.fmt(rej / iters),
                    "note": ""})
    return out


def main() -> int:
    L.OUT_DIR.mkdir(parents=True, exist_ok=True)
    L.live("deep_lineage", [("align panel", False), ("stratified test", False)])
    rows = build_panel()
    cols = ["upid", "organism", "taxid", "group", "kingdom", "phylum", "class",
            "order", "n_itpr", "accession", "protein_count", "seq_length",
            "n_pi_plc", "n_plc_x", "n_plc_y", "plc_status",
            "core_n_ref", "core_n_covered", "core_coverage", "core_identity",
            "pore_n_ref", "pore_n_covered", "pore_coverage", "pore_identity",
            "delta_core_minus_pore", "in_test", "drop_reason"]
    L.write_tsv(L.OUT_DIR / "deep_lineage_panel.tsv", rows, cols)
    n_live = sum(1 for r in rows if r["in_test"] is True)
    L.log(f"panel: {len(rows)} proteomes, {n_live} enter the test")

    st = stratified(rows)
    L.write_tsv(L.OUT_DIR / "deep_lineage_test.tsv", st, list(st[0]))
    cv = covariates(rows)
    L.write_tsv(L.OUT_DIR / "deep_lineage_covariates.tsv", cv, list(cv[0]))
    pairs, summary = divergence_matched(rows)
    L.write_tsv(L.OUT_DIR / "deep_lineage_matched.tsv", pairs, list(pairs[0]))
    L.write_tsv(L.OUT_DIR / "deep_lineage_matched_summary.tsv", [summary],
                list(summary))
    L.log(f"divergence-matched: {summary['n_matched']} of "
          f"{summary['n_plc_absent_in_test']} matched, median difference "
          f"{summary['median_difference']}, p_sign={summary['p_sign']}")
    rc = ryr_control(rows)
    L.write_tsv(L.OUT_DIR / "deep_lineage_ryr_control.tsv", rc, list(rc[0]))
    mp = matched_power(pairs)
    L.write_tsv(L.OUT_DIR / "deep_lineage_matched_power.tsv", mp, list(mp[0]))
    pw = power(rows)
    L.write_tsv(L.OUT_DIR / "deep_lineage_power.tsv", pw, list(pw[0]))
    for r in st:
        if r["testable"]:
            L.log(f"{r['stratum']}: {r['n_plc_absent']} absent vs "
                  f"{r['n_plc_present']} present, shift={r['shift']}, "
                  f"p={r['p_mannwhitney']}")
    L.live("deep_lineage", [("align panel", True), ("stratified test", True)])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
