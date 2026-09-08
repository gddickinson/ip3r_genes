"""S17 stage 5 — per-site dN/dS, mapped onto the same architecture.

S9 answered the whole-gene question: every paralog is far below neutrality
(one-ratio omega 0.024-0.043) and no site model found a positively selected
site. What it did not answer is *where* the constraint sits. This stage runs
HyPhy **FEL** on S9's three codon alignments — a per-site maximum-likelihood
estimate of the synonymous (alpha) and non-synonymous (beta) rate — and carries
every site onto the human reference protein, so selection and conservation are
reported on one set of coordinates and one set of elements.

Two limits, both stated in the report rather than worked around:

* **Power is per-element, not per-site.** These alignments have 19/13/19
  sequences and S9 measured median dS of 4.6-13.5 — deeply saturated, *within*
  a paralog and not only between them. A single codon's omega is noisy; the
  distribution over a 288-residue channel domain is not. Per-site calls are
  reported with Benjamini-Hochberg FDR and treated as a screen.
* **The coordinate carry is validated, not assumed.** The codon alignment is
  trimAl-trimmed, so its reference row is a *subsequence* of the reference
  protein. It is realigned to the full protein and every mapped site must agree
  in amino acid; a site that does not is written with no residue rather than
  with a plausible wrong one, and the agreement rate goes into the output.

    python scripts/s17_fel.py            # run + parse (cached)
    python scripts/s17_fel.py --force
"""

from __future__ import annotations

import argparse
import json
import statistics as st
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_domains as D                                        # noqa: E402
import s17_lib as L                                            # noqa: E402

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

#: HyPhy's own upper bound on a per-site rate. A site that reaches it is not a
#: site with a huge synonymous rate; it is a site the data cannot identify.
RATE_BOUND = 10000.0
#: Ceiling above which alpha is treated as unidentifiable for forming a ratio.
ALPHA_USABLE_MAX = 100.0


def ref_tip(paralog: str) -> str:
    """The human tip's short PAML code, read from S9's committed mapping.

    Never typed: `tip_codes.tsv` exists precisely so no later task has to guess
    how a 60-character label was abbreviated.
    """
    want = L.REFERENCES[paralog][0]
    for r in L.read_tsv(L.SEL_DIR / "tip_codes.tsv"):
        if r["label"] == want:
            return r["code"]
    raise SystemExit(f"no tip code for {want} in tip_codes.tsv")


def hyphy_bin() -> str:
    """HyPhy is env-resident here (D18); PATH first, then the recorded env."""
    import shutil
    return shutil.which("hyphy") or "/opt/anaconda3/envs/piezo1/bin/hyphy"


def run_fel(paralog: str, out_dir: Path, force: bool = False) -> Path:
    dest = out_dir / "fel" / f"{paralog}.FEL.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not force:
        return dest
    aln = L.SEL_DIR / f"codon_{paralog}.fasta"
    tree = L.SEL_DIR / f"tree_{paralog}.nwk"
    log = dest.with_suffix(".log")
    t0 = time.time()
    with log.open("w") as fl:
        subprocess.run([hyphy_bin(), "fel", "--alignment", str(aln),
                        "--tree", str(tree), "--branches", "All",
                        "--ci", "Yes", "--output", str(dest)],
                       stdout=fl, stderr=subprocess.STDOUT, check=True)
    print(f"[s17] FEL {paralog}: {time.time() - t0:.0f} s")
    return dest


def _bh(pvals: list[float]) -> list[float]:
    """Benjamini-Hochberg adjusted p-values, order preserved."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    adj, prev = [1.0] * n, 1.0
    for rank, i in enumerate(reversed(order), start=1):
        k = n - rank + 1
        prev = min(prev, pvals[i] * n / k)
        adj[i] = min(1.0, prev)
    return adj


def _translate_row(row: str) -> str:
    aa = []
    for i in range(0, len(row) - 2, 3):
        cod = row[i:i + 3].upper()
        if cod == "---":
            continue
        aa.append(CODON_TABLE.get(cod, "X"))
    return "".join(aa)


def parse(paralog: str, out_dir: Path) -> tuple[list[dict], dict]:
    d = json.loads((out_dir / "fel" / f"{paralog}.FEL.json").read_text())
    headers = [h[0] for h in d["MLE"]["headers"]]
    content = d["MLE"]["content"]["0"]
    idx = {name: i for i, name in enumerate(headers)}

    aln = L.read_fasta(L.SEL_DIR / f"codon_{paralog}.fasta")
    row = aln[ref_tip(paralog)]
    trimmed_prot = _translate_row(row)
    site_to_trim, t = {}, 0
    for s in range(len(row) // 3):
        if row[s * 3:s * 3 + 3] != "---":
            t += 1
            site_to_trim[s + 1] = t

    _label, acc, _desc = L.REFERENCES[paralog]
    ref_seq = L.uniprot_fasta(acc)
    xfer = L.transfer_positions(trimmed_prot, ref_seq)

    elem = D.residue_element(out_dir).get(paralog, {})
    sites = D.site_index(out_dir).get(paralog, {})

    pvals = [float(r[idx["p-value"]]) for r in content]
    adj = _bh(pvals)

    rows, matched, mapped = [], 0, 0
    for i, r in enumerate(content):
        alpha, beta = float(r[idx["alpha"]]), float(r[idx["beta"]])
        tpos = site_to_trim.get(i + 1)
        resi = xfer.get(tpos) if tpos else None
        if resi:
            mapped += 1
            if ref_seq[resi - 1] == trimmed_prot[tpos - 1]:
                matched += 1
            else:
                resi = None                  # refuse an unvalidated coordinate
        rows.append({
            "paralog": paralog, "site": i + 1, "resi": resi or "",
            "ref_acc": acc, "aa": ref_seq[resi - 1] if resi else "",
            "element": elem.get(resi, "") if resi else "",
            "ip3_contact": (resi in sites.get("ip3_contact", set())) if resi else "",
            "filter_lining": (resi in sites.get("filter_lining", set())) if resi else "",
            "gate_lining": (resi in sites.get("gate_lining", set())) if resi else "",
            "alpha": round(alpha, 4), "beta": round(beta, 4),
            "omega": round(beta / alpha, 4) if alpha > 0 else "",
            "lrt": round(float(r[idx["LRT"]]), 4),
            "p_value": L.fmt(pvals[i]), "q_value": L.fmt(adj[i]),
            "alpha_at_bound": alpha >= RATE_BOUND,
            "verdict": ("purifying" if adj[i] < 0.05 and beta < alpha else
                        "diversifying" if adj[i] < 0.05 and beta > alpha
                        else "neutral_or_undetermined"),
        })
    stats = {"paralog": paralog, "n_sites": len(rows), "n_mapped": mapped,
             "n_aa_agree": matched,
             "aa_agreement": round(matched / mapped, 4) if mapped else 0.0,
             "n_alpha_at_bound": sum(1 for x in rows if x["alpha_at_bound"]),
             "n_purifying_q05": sum(1 for x in rows if x["verdict"] == "purifying"),
             "n_diversifying_q05": sum(1 for x in rows
                                       if x["verdict"] == "diversifying")}
    return rows, stats


def by_element(rows: list[dict]) -> list[dict]:
    """Element-level selection, using statistics that survive dS saturation.

    The obvious aggregate — sum of beta over sum of alpha — is **not** usable
    here. S9 measured median dS of 4.6-13.5 on these alignments and FEL duly
    pins alpha at its upper bound at a share of sites; the sum is then dominated
    by sites whose synonymous rate is unidentifiable rather than large, and
    every element's omega collapses towards zero. Three statistics instead:

    * `median_beta` — the non-synonymous rate, which needs no ratio at all and
      is comparable **between elements of the same alignment** (FEL normalises
      tree scale per alignment, so it is not comparable across paralogs; S9's
      codeml omega is what compares those);
    * `frac_purifying_q05` — the share of sites FEL calls significantly
      constrained, a count that one site cannot inflate;
    * `median_omega_usable` — per-site beta/alpha over the sites where alpha is
      identifiable, with `n_usable_alpha` beside it so a reader can see how
      much of the element it rests on.
    """
    agg: dict[tuple[str, str], list] = {}
    for r in rows:
        if not r["element"]:
            continue
        keys = [(r["paralog"], r["element"]), (r["paralog"], "WHOLE_PROTEIN")]
        for flag in ("ip3_contact", "filter_lining", "gate_lining"):
            if r[flag] is True:
                keys.append((r["paralog"], flag))
        for key in keys:
            a = agg.setdefault(key, [[], [], 0, 0])
            a[0].append(r["beta"])
            if 0 < r["alpha"] < ALPHA_USABLE_MAX:
                a[1].append(r["beta"] / r["alpha"])
            a[2] += 1 if r["verdict"] == "purifying" else 0
            a[3] += 1 if r["alpha_at_bound"] else 0
    out = []
    for (paralog, element), (betas, om, npur, nbound) in sorted(agg.items()):
        n = len(betas)
        out.append({
            "paralog": paralog, "element": element, "n_sites": n,
            "median_beta": round(st.median(betas), 4),
            "mean_beta": round(st.mean(betas), 4),
            "frac_purifying_q05": round(npur / n, 4) if n else "",
            "n_usable_alpha": len(om),
            "median_omega_usable": round(st.median(om), 4) if om else "",
            "n_alpha_at_bound": nbound,
        })
    return out


def run(out_dir: Path = L.OUT_DIR, force: bool = False) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    all_rows, all_stats = [], []
    for paralog in L.PARALOGS:
        run_fel(paralog, out_dir, force=force)
        rows, stats = parse(paralog, out_dir)
        all_rows.extend(rows)
        all_stats.append(stats)
        print(f"[s17] {paralog}: {stats['n_sites']} sites, {stats['n_mapped']} "
              f"mapped (aa agreement {stats['aa_agreement']:.1%}), "
              f"{stats['n_purifying_q05']} purifying / "
              f"{stats['n_diversifying_q05']} diversifying at q<0.05")
    L.write_tsv(out_dir / "fel_sites.tsv", all_rows)
    L.write_tsv(out_dir / "fel_status.tsv", all_stats)
    L.write_tsv(out_dir / "selection_by_element.tsv", by_element(all_rows))
    return {"n_sites": len(all_rows)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    run(a.out, force=a.force)


if __name__ == "__main__":
    main()
