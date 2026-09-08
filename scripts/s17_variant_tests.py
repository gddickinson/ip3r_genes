"""S17 stage 4b — the tests that turn the variant harvest into a result.

Kept apart from `s17_variants.py` (which only fetches and parses) so the
analysis can be read without wading through E-utilities plumbing.

* `constraint_test` — is the per-site score a *classifier*? Scored against
  ClinVar's own labelling for each of the four conservation layers. Two
  contrasts, not one, because ClinVar's benign set is itself ascertainment-
  biased: **P/LP vs B/LB** is the question a clinician asks, and **P/LP vs the
  whole protein** is the control for it — if the benign set were simply drawn
  from unconstrained regions the first AUC would be high and the second flat.
  This is also the arbitration stage 2 promised: if the ~250-orthologue `deep`
  layer does not beat the 11-18 sequence `shallow` one, the depth was not worth
  having.
* `paralog_audit` — do the labelled positions survive in the other two
  paralogs? Each labelled human position is carried through msa_v2 (the
  alignment S6 built and S7 validated, and the only one containing all three
  references) and asked whether the amino acid is the same, against two
  controls: the benign positions, and every aligned column of the pair.
* `vus_stratification` — the resource's actual output. 88 % of this family's
  ClinVar missense record is `uncertain significance`, and a per-site score is
  what a constraint resource offers such a variant. Each VUS is placed against
  the *labelled* distributions of the same gene, and the table reports how many
  sit above the P/LP median and below the B/LB median. It is deliberately a
  stratification and not a call: it says where a variant sits on an axis the
  labelled variants separate on, which is a different claim from pathogenicity.
* `variants_by_element` — where in the architecture each class of variant sits,
  with the element's share of the protein beside it, because "31 % of
  pathogenic variants are in the channel" means nothing until you know the
  channel is 10 % of the protein.

Every p-value is written through `s17_lib.fmt`, at `%.6g`: a p of 2e-17
rendered at four decimal places reads back as `0.0000`, which hides how strong
a claim is rather than how weak (S15a's rule).
"""

from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402

LAYERS = ("deep", "shallow", "vert", "family")

#: A layer is only scored on residues where it is reliable. For `deep` that is
#: S17's own occupancy flag; the msa_v2 layers carry their own occupancy.
MIN_OCCUPANCY = 0.50


def _auc(pos: list[float], neg: list[float]) -> tuple[float, float]:
    """ROC AUC of `pos` over `neg`, with the Mann-Whitney p it comes from."""
    from scipy.stats import mannwhitneyu
    if len(pos) < 3 or len(neg) < 3:
        return float("nan"), float("nan")
    r = mannwhitneyu(pos, neg, alternative="greater")
    return r.statistic / (len(pos) * len(neg)), float(r.pvalue)


def _score(site: dict, layer: str) -> float | None:
    v = site.get(f"{layer}_jsd", "")
    if v in ("", None):
        return None
    occ = site.get(f"{layer}_occupancy", "")
    if occ not in ("", None) and float(occ) < MIN_OCCUPANCY:
        return None
    return float(v)


def _scorable(sites: dict[str, dict[int, dict]], gene: str, resi: int,
              layers=LAYERS) -> bool:
    s = sites[gene].get(resi)
    return bool(s) and all(_score(s, ly) is not None for ly in layers)


def constraint_test(variants: list[dict], sites: dict[str, dict[int, dict]]
                    ) -> list[dict]:
    """Does constraint separate pathogenic from benign, layer by layer?

    Three contrasts per layer, and the third is what makes the layer comparison
    an answer rather than four separate ones:

    * **P/LP vs B/LB** — the question a clinician asks.
    * **P/LP vs whole protein** — the control for ClinVar's own ascertainment.
      If the benign set were simply drawn from unconstrained regions the first
      AUC would be high and this one flat.
    * **P/LP vs B/LB, all layers scorable** — the same test restricted to the
      positions where **every** layer has a reliable score. The layers do not
      cover the same residues (the msa_v2 layers lose gappy columns the deep
      alignment fills, and vice versa), so ranking four AUCs measured on four
      slightly different variant sets would compare the sets as much as the
      layers. This contrast holds the set fixed.
    """
    out = []
    for gene in list(L.PARALOGS) + ["POOLED"]:
        rows = [v for v in variants if v["source"] == "clinvar"
                and (gene == "POOLED" or v["gene"] == gene)]
        genes = L.PARALOGS if gene == "POOLED" else (gene,)
        for layer in LAYERS:
            def vals(b, matched=False):
                got, seen = [], set()
                for v in rows:
                    key = (v["gene"], int(v["resi"]))
                    if v["class_bucket"] != b or key in seen:
                        continue
                    resi = int(v["resi"])
                    if matched and not _scorable(sites, v["gene"], resi):
                        continue
                    s = sites[v["gene"]].get(resi)
                    x = _score(s, layer) if s else None
                    if x is not None:
                        seen.add(key)
                        got.append(x)
                return got

            protein = [x for g in genes
                       for x in (_score(s, layer) for s in sites[g].values())
                       if x is not None]
            path, ben = vals("P/LP"), vals("B/LB")
            mpath, mben = vals("P/LP", True), vals("B/LB", True)
            for contrast, pos, neg in (
                    ("P/LP vs B/LB", path, ben),
                    ("P/LP vs whole protein", path, protein),
                    ("P/LP vs B/LB, all layers scorable", mpath, mben)):
                auc, p = _auc(pos, neg)
                out.append({
                    "gene": gene, "layer": layer, "contrast": contrast,
                    "n_positive": len(pos), "n_negative": len(neg),
                    "mean_positive": round(sum(pos) / len(pos), 4) if pos else "",
                    "mean_negative": round(sum(neg) / len(neg), 4) if neg else "",
                    "auc": round(auc, 4) if auc == auc else "",
                    "p_mannwhitney": L.fmt(p) if p == p else "",
                })
    return out


def paralog_audit(variants: list[dict], out_dir: Path
                  ) -> tuple[list[dict], list[dict]]:
    """Do the positions human genetics has labelled survive in the sisters?

    Carried through msa_v2, because that is the alignment S6 built and S7
    validated and the only one containing all three references; the per-paralog
    deep alignments contain one paralog each and cannot answer a cross-paralog
    question.
    """
    from scipy.stats import fisher_exact

    msa = L.read_fasta(L.MSA_DIR / "aln.fasta")
    per_variant, summary = [], []
    pooled: dict[tuple[str, str], list[int]] = {}

    for gene in L.PARALOGS:
        label, acc, _d = L.REFERENCES[gene]
        ref_seq = L.uniprot_fasta(acc)
        msa_seq = msa[label].replace("-", "")
        fwd = L.transfer_positions(msa_seq, ref_seq)
        back = {v: k for k, v in fwd.items()}
        r2c = L.residue_to_col(msa[label])

        for other in L.PARALOGS:
            if other == gene:
                continue
            olabel, _oacc, _od = L.REFERENCES[other]
            orow = msa[olabel]
            c2ro = L.col_to_residue(orow)

            bg_same = bg_tot = 0
            for ca, cb in zip(msa[label], orow):
                if ca in "-." or cb in "-.":
                    continue
                bg_tot += 1
                bg_same += 1 if ca == cb else 0

            counts = {"P/LP": [0, 0], "B/LB": [0, 0], "VUS": [0, 0]}
            seen: set[tuple[int, str]] = set()
            for v in variants:
                b = v["class_bucket"]
                if v["gene"] != gene or b not in counts:
                    continue
                key = (int(v["resi"]), b)
                if key in seen:
                    continue                      # positions, not alleles
                seen.add(key)
                resi = int(v["resi"])
                col = r2c.get(back.get(resi, -1))
                oaa = orow[col] if col is not None else None
                aligned = oaa is not None and oaa not in "-."
                identical = aligned and oaa == v["ref_aa"]
                counts[b][0] += 1 if identical else 0
                counts[b][1] += 1 if aligned else 0
                per_variant.append({
                    "gene": gene, "other": other, "resi": resi,
                    "ref_aa": v["ref_aa"], "class_bucket": b,
                    "msa_col": col if col is not None else "",
                    "other_resi": c2ro.get(col, "") if col is not None else "",
                    "other_aa": oaa or "-", "aligned": aligned,
                    "identical": identical,
                })
            for b, (same, tot) in counts.items():
                if tot == 0:
                    continue
                odds, p = fisher_exact([[same, tot - same],
                                        [bg_same, bg_tot - bg_same]])
                summary.append({
                    "gene": gene, "other": other,
                    "contrast": f"{b} vs all aligned columns",
                    "class_bucket": b, "n_positions": tot, "n_identical": same,
                    "frac_identical": round(same / tot, 4),
                    "comparator_frac": round(bg_same / bg_tot, 4),
                    "comparator_n": bg_tot, "odds_ratio": round(float(odds), 3),
                    "p_fisher": L.fmt(float(p)),
                })
            (ps, pt), (bs, bt) = counts["P/LP"], counts["B/LB"]
            if pt and bt:
                odds, p = fisher_exact([[ps, pt - ps], [bs, bt - bs]])
                summary.append({
                    "gene": gene, "other": other, "contrast": "P/LP vs B/LB",
                    "class_bucket": "P/LP", "n_positions": pt, "n_identical": ps,
                    "frac_identical": round(ps / pt, 4),
                    "comparator_frac": round(bs / bt, 4), "comparator_n": bt,
                    "odds_ratio": round(float(odds), 3),
                    "p_fisher": L.fmt(float(p)),
                })
                acc4 = pooled.setdefault(("POOLED", other), [0, 0, 0, 0])
                acc4[0] += ps
                acc4[1] += pt
                acc4[2] += bs
                acc4[3] += bt

    for (gene, other), (ps, pt, bs, bt) in sorted(pooled.items()):
        if not (pt and bt):
            continue
        odds, p = fisher_exact([[ps, pt - ps], [bs, bt - bs]])
        summary.append({
            "gene": gene, "other": other, "contrast": "P/LP vs B/LB",
            "class_bucket": "P/LP", "n_positions": pt, "n_identical": ps,
            "frac_identical": round(ps / pt, 4),
            "comparator_frac": round(bs / bt, 4), "comparator_n": bt,
            "odds_ratio": round(float(odds), 3), "p_fisher": L.fmt(float(p)),
        })
    L.write_tsv(out_dir / "paralog_variant_positions.tsv", per_variant)
    return per_variant, summary


def vus_stratification(variants: list[dict],
                       sites: dict[str, dict[int, dict]]) -> list[dict]:
    """Where the uncertain variants sit on the axis the labelled ones separate on.

    Reported per gene *and* per layer, because the layers do not agree about
    which is the best classifier and a resource that quoted one number would be
    hiding that. The thresholds are the labelled distributions' own medians, so
    nothing here is a cut chosen to make a count.
    """
    import statistics as st

    out = []
    for gene in L.PARALOGS:
        for layer in LAYERS:
            def vals(bucket):
                got, seen = [], set()
                for v in variants:
                    key = int(v["resi"])
                    if (v["gene"] != gene or v["class_bucket"] != bucket
                            or key in seen):
                        continue
                    s = sites[gene].get(key)
                    x = _score(s, layer) if s else None
                    if x is not None:
                        seen.add(key)
                        got.append(x)
                return got

            path, ben, vus = vals("P/LP"), vals("B/LB"), vals("VUS")
            if not vus or not path or not ben:
                continue
            p_med, b_med = st.median(path), st.median(ben)
            out.append({
                "gene": gene, "layer": layer, "n_vus": len(vus),
                "n_pathogenic": len(path), "n_benign": len(ben),
                "median_pathogenic": round(p_med, 4),
                "median_benign": round(b_med, 4),
                "median_vus": round(st.median(vus), 4),
                "n_vus_above_pathogenic_median": sum(1 for x in vus if x >= p_med),
                "frac_vus_above_pathogenic_median": round(
                    sum(1 for x in vus if x >= p_med) / len(vus), 4),
                "n_vus_below_benign_median": sum(1 for x in vus if x <= b_med),
                "frac_vus_below_benign_median": round(
                    sum(1 for x in vus if x <= b_med) / len(vus), 4),
            })
    return out


def variants_by_element(variants: list[dict]) -> list[dict]:
    """Where each class of variant sits, against the element's own share.

    An element's share of the *protein* is the right denominator: without it a
    count in the largest element always looks like an enrichment. The residue
    denominators come from the same `element` column the variants were joined
    on, so the two cannot disagree.
    """
    from scipy.stats import fisher_exact

    out = []
    for gene in L.PARALOGS:
        acc = L.REFERENCES[gene][1]
        sites = L.read_tsv(L.OUT_DIR / f"constraint_{gene}_{acc}.tsv")
        n_res = Counter(r["element"] for r in sites)
        total_res = sum(n_res.values())
        for b in ("P/LP", "B/LB", "VUS"):
            pos = {(int(v["resi"])) for v in variants
                   if v["gene"] == gene and v["class_bucket"] == b}
            if not pos:
                continue
            by_el = Counter(r["element"] for r in sites
                            if int(r["resi"]) in pos)
            tot = sum(by_el.values())
            for el, n in sorted(by_el.items(), key=lambda x: -x[1]):
                odds, p = fisher_exact(
                    [[n, tot - n], [n_res[el] - n, total_res - n_res[el] - (tot - n)]])
                out.append({
                    "gene": gene, "class_bucket": b, "element": el,
                    "n_positions": n, "frac_of_class": round(n / tot, 4),
                    "n_residues": n_res[el],
                    "frac_of_protein": round(n_res[el] / total_res, 4),
                    "enrichment": round((n / tot) / (n_res[el] / total_res), 3)
                    if n_res[el] else "",
                    "odds_ratio": round(float(odds), 3), "p_fisher": L.fmt(float(p)),
                })
    return out
