"""Every committed S18 table, plus the by-source breakdown and the stats file.

Two aggregations live here rather than in the audit modules because they are
*comparisons* and the audit is a measurement:

**The by-source breakdown (D9).** RefSeq and submitter-deposited GenBank gene
sets are not comparable evidence, so the state counts are split by which
archive serves the annotation and each state is tested with Fisher's exact,
BH-corrected across the family of states. The comparison carries a confounder
in the same table rather than in a caveat: `GCA_` assemblies are on average
less contiguous, and a locus on a contig too short to hold the gene cannot be
annotated completely by anybody. `by_source()` therefore reports the raw
contrast *and* the contrast restricted to loci above D4's bar, which is the
control that decides whether a source difference is a source difference.

**The RyR control.** Every state count is reported for the sister family
beside the ITPR cells, from the same assemblies and the same instrument. A
family whose loci are as badly recorded as its sister's has an annotation
problem; one whose loci are worse has a *family-specific* one, and only the
control can tell those apart.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s18_lib as L                                               # noqa: E402
import s18_locus_rules as R                                       # noqa: E402


def state_counts(rows: list[dict], keys: list[str]) -> list[dict]:
    """State counts grouped by `keys`, every state present as an explicit 0.

    A missing row reads as *not measured*; a zero is a result (S15b's rule).
    """
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        groups.setdefault(tuple(r.get(k, "") for k in keys), []).append(r)
    out = []
    for gk, rs in sorted(groups.items(), key=lambda kv: [str(x) for x in kv[0]]):
        d = dict(zip(keys, gk))
        d["n_loci"] = len(rs)
        for st in R.STATES:
            n = sum(1 for r in rs if r["state"] == st)
            d[f"n_{st}"] = n
            d[f"frac_{st}"] = L.frac(n, len(rs))
        out.append(d)
    return out


def by_source(rows: list[dict]) -> list[dict]:
    """D9's contrast, twice: over every scorable locus and above D4's bar."""
    out = []
    for scope, sel in (
            ("all scorable loci", lambda r: True),
            ("loci on a contig that spans the gene (D4)",
             lambda r: str(r.get("contig_spans_gene")) in ("1", "True", "true"))):
        pool = [r for r in rows if r["gene_set"] == "present" and sel(r)]
        ref = [r for r in pool if r["source"] == "RefSeq"]
        gb = [r for r in pool if r["source"] == "GenBank"]
        rows_here = []
        for st in R.STATES:
            if st in ("no_gene_set", "gff_unavailable"):
                continue
            a = sum(1 for r in ref if r["state"] == st)
            c = sum(1 for r in gb if r["state"] == st)
            rows_here.append({
                "scope": scope, "state": st,
                "n_refseq": a, "n_refseq_total": len(ref),
                "frac_refseq": L.frac(a, len(ref)),
                "n_genbank": c, "n_genbank_total": len(gb),
                "frac_genbank": L.frac(c, len(gb)),
                "p_fisher": round(L.fisher_2x2(a, len(ref) - a,
                                               c, len(gb) - c), 6),
            })
        ps = [r["p_fisher"] for r in rows_here]
        for r, q in zip(rows_here, L.benjamini_hochberg(ps)):
            r["q_bh"] = q
        out.extend(rows_here)
    return out


def family_vs_control(rows: list[dict]) -> list[dict]:
    """The family against its sister, from the same assemblies and instrument.

    The audit's headline is a failure rate, and a failure rate means nothing
    without knowing whether it is this family's or the genome annotation's.
    The RyR control answers it: the same 274 gene sets, the same alignment
    pipeline, the same five states. A family recorded as badly as its sister
    has an annotation problem; one recorded worse has a *family-specific* one,
    and only this comparison separates them.

    Run twice — over every scorable locus and over the loci on a contig that
    can carry the gene — because RyR genes are 1.8x longer and therefore fail
    D4's bar more often, which would make the sister family look worse for a
    reason that has nothing to do with how it is annotated.
    """
    out = []
    fails = ("unannotated", "noncoding", "fragmentary", "split")
    for scope, sel in (
            ("all scorable loci", lambda r: True),
            ("loci on a contig that spans the gene (D4)",
             lambda r: str(r.get("contig_spans_gene")) in ("1", "True", "true"))):
        pool = [r for r in rows if r["gene_set"] == "present" and sel(r)]
        itpr = [r for r in pool if not r["is_control"]]
        ryr = [r for r in pool if r["is_control"]]
        a = sum(1 for r in itpr if r["state"] in fails)
        c = sum(1 for r in ryr if r["state"] in fails)
        out.append({
            "scope": scope, "measure": "any annotation failure",
            "n_itpr": len(itpr), "n_itpr_failing": a,
            "frac_itpr_failing": L.frac(a, len(itpr)),
            "n_ryr": len(ryr), "n_ryr_failing": c,
            "frac_ryr_failing": L.frac(c, len(ryr)),
            "p_fisher": round(L.fisher_2x2(a, len(itpr) - a,
                                           c, len(ryr) - c), 6)})
        for st in R.STATES:
            if st in ("no_gene_set", "gff_unavailable"):
                continue
            a = sum(1 for r in itpr if r["state"] == st)
            c = sum(1 for r in ryr if r["state"] == st)
            out.append({
                "scope": scope, "measure": st,
                "n_itpr": len(itpr), "n_itpr_failing": a,
                "frac_itpr_failing": L.frac(a, len(itpr)),
                "n_ryr": len(ryr), "n_ryr_failing": c,
                "frac_ryr_failing": L.frac(c, len(ryr)),
                "p_fisher": round(L.fisher_2x2(a, len(itpr) - a,
                                               c, len(ryr) - c), 6)})
    ps = [r["p_fisher"] for r in out]
    for r, q in zip(out, L.benjamini_hochberg(ps)):
        r["q_bh"] = q
    return out


def verdict_counts(rows: list[dict], field: str, keys: list[str]
                   ) -> list[dict]:
    """Name / symbol verdict counts, every verdict present as an explicit 0."""
    groups: dict[tuple, list[dict]] = {}
    for r in rows:
        groups.setdefault(tuple(r.get(k, "") for k in keys), []).append(r)
    out = []
    for gk, rs in sorted(groups.items(), key=lambda kv: [str(x) for x in kv[0]]):
        d = dict(zip(keys, gk))
        d["field"] = field
        d["n"] = len(rs)
        seen = set(R.NAME_VERDICTS) | {"paralog_not_callable"} | \
            {r.get(field, "") for r in rs}
        for v in sorted(seen):
            if not v:
                continue
            n = sum(1 for r in rs if r.get(field) == v)
            d[f"n_{v}"] = n
            d[f"frac_{v}"] = L.frac(n, len(rs))
        out.append(d)
    return out


def headline(locus_rows: list[dict], protein_rows: list[dict],
             corrections: list[dict], zero_rows: list[dict],
             calibration: dict) -> dict:
    """The numbers the report is not allowed to compute for itself (D13)."""
    scorable = [r for r in locus_rows
                if r["gene_set"] == "present" and not r["is_control"]]
    ctrl = [r for r in locus_rows
            if r["gene_set"] == "present" and r["is_control"]]
    spans = [r for r in scorable
             if str(r.get("contig_spans_gene")) in ("1", "True", "true")]

    def share(pool, sts):
        return L.frac(sum(1 for r in pool if r["state"] in sts), len(pool))

    fails = ("unannotated", "noncoding", "fragmentary", "split")
    vert = [r for r in protein_rows if r["paralog_askable"]]
    return {
        "n_genomes": len({r["accession"] for r in locus_rows}),
        "n_loci": len(locus_rows),
        "n_itpr_loci": len([r for r in locus_rows if not r["is_control"]]),
        "n_scorable_itpr_loci": len(scorable),
        "n_control_loci": len(ctrl),
        "n_loci_no_gene_set": len([r for r in locus_rows
                                   if r["state"] == "no_gene_set"]),
        "n_loci_cds_unavailable": len([r for r in locus_rows
                                       if r["state"] == "cds_unavailable"]),
        "frac_itpr_complete": share(scorable, ("complete",)),
        "frac_itpr_any_failure": share(scorable, fails),
        "frac_control_complete": share(ctrl, ("complete",)),
        "frac_control_any_failure": share(ctrl, fails),
        "frac_itpr_failure_above_d4": share(spans, fails),
        "n_itpr_unannotated": sum(1 for r in scorable
                                  if r["state"] == "unannotated"),
        "n_itpr_noncoding": sum(1 for r in scorable
                                if r["state"] == "noncoding"),
        "calibration_bar": calibration.get("bar"),
        "calibration_median": calibration.get("median"),
        "calibration_frac_below_bar": calibration.get("frac_below_bar"),
        "calibration_n": calibration.get("n"),
        "n_protein_records": len(protein_rows),
        "n_protein_vertebrate": len(vert),
        "n_protein_seq_family_disagrees": sum(
            1 for r in protein_rows
            if r["seq_family"] in ("ITPR", "RYR")
            and r["seq_family"] != r["census_call"]),
        "n_name_wrong_family": sum(1 for r in protein_rows
                                   if r["name_verdict"] == "wrong_family"),
        "n_name_family_ambiguous": sum(1 for r in protein_rows
                                       if r["name_verdict"] == "family_ambiguous"),
        "n_symbol_wrong_paralog": sum(
            1 for r in protein_rows
            if r["symbol_verdict"] == "correct_family_wrong_paralog"),
        "n_symbol_placeholder": sum(1 for r in protein_rows
                                    if r["symbol_verdict"] == "placeholder"),
        "n_symbol_absent": sum(1 for r in protein_rows
                               if r["symbol_verdict"] == "absent"),
        "n_paralog_not_callable": sum(
            1 for r in protein_rows
            if r["symbol_verdict"] == "paralog_not_callable"),
        "n_corrections": len(corrections),
        "n_corrections_high": sum(1 for r in corrections
                                  if r["priority"] == "high"),
        "n_corrections_vetoed": sum(1 for r in corrections if r["vetoed"]),
        "n_zero_hit_proteomes": len(zero_rows),
        "n_zero_hit_gene_caller": sum(1 for r in zero_rows
                                      if r["verdict"] == "gene_caller_missed_it"),
    }


def write_stats(paths: dict[str, Path], params: dict, self_test: str,
                head: dict) -> Path:
    out = L.OUT / "annotation_audit_stats.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({
        "task": "S18", "params": params, "self_test": self_test,
        "headline": head,
        "tables": {k: {"path": str(v.relative_to(L.PROJECT)),
                       "sha256": L.sha256(v)}
                   for k, v in sorted(paths.items()) if v and v.exists()},
    }, indent=1))
    return out
