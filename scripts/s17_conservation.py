"""S17 stage 3 — per-site constraint in human ITPR1, ITPR2 and ITPR3.

Four conservation layers are computed on the same residues, because they answer
different questions and a single number would blur them:

* **`deep`** — within one paralog, the 240-260 vertebrate orthologues stage 2
  built. "Has *this* gene tolerated change at this position?" This is the layer
  a variant-interpretation resource needs, and the only one deep enough to
  score a column from more than a dozen observations.
* **`vert`** — msa_v2's vertebrate ITPR tips, all three paralogs together.
  "Was this position already fixed before the duplications?"
* **`family`** — every msa_v2 ITPR tip, protists and plants included: the
  eukaryote-wide floor, the positions the fold cannot do without. **The six RyR
  tips are excluded.** They are in msa_v2 as S6's outgroup, and a layer
  containing them would score what the *superfamily* conserves — a different
  question, and one D14 keeps separate everywhere else in this project.
* **`shallow`** — msa_v2 restricted to one paralog (15/11/18 tips). Not a
  claim, a **control**: it is what S17 would have had without stage 2, and the
  ClinVar test in `s17_variant_tests.py` uses the difference between `shallow`
  and `deep` to show empirically whether the extra depth buys discrimination or
  noise.

Every layer is sequence-weighted (Henikoff & Henikoff) and reported with the
column occupancy it is conditional on, so a score computed from four sequences
in a gappy region cannot be mistaken for a well-supported one.

Coordinates are the UniProt canonical numbering of each human reference. All
three msa_v2 reference tips *are* the canonical UniProt proteins here, so the
transfer to msa_v2 coordinates is the identity — but it is computed rather than
assumed, and the module says which it used.

    python scripts/s17_conservation.py
"""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_domains as D                                        # noqa: E402
import s17_lib as L                                            # noqa: E402

#: Occupancy below which a column's score is reported but flagged unreliable.
MIN_OCCUPANCY = 0.50

LAYERS = ("deep", "shallow", "vert", "family")

#: The within-protein control every element is tested against: the sequence
#: between the named domains, plus the two termini. A protein-wide mean would
#: be the wrong control — it contains the elements being tested.
#:
#: The membership test is a **prefix for linkers and an exact match for the
#: termini**, because `nterm_trefoil` is a domain whose name begins with
#: `nterm`: a `startswith` test on both quietly swept PF08709 — 225 residues of
#: the element the ligand question is about — into the control, and every other
#: element was then being compared against a set containing one of them
#: (`s17_test_constraint.py` T8).
CONTROL_LINKER_PREFIX = "linker_"
CONTROL_TERMINI = ("nterm", "cterm")


def _is_control(element: str) -> bool:
    return (element.startswith(CONTROL_LINKER_PREFIX)
            or element in CONTROL_TERMINI)


def _msa_subsets() -> dict[str, list[str]]:
    """Named subsets of msa_v2 labels, tree-corrected (S7) and RyR-free."""
    groups = L.msa_groups()
    out: dict[str, list[str]] = {p: [] for p in L.PARALOGS}
    vert, family = [], []
    for label, g in sorted(groups.items()):
        if g == "RYR":
            continue
        family.append(label)
        if g in L.PARALOGS:
            out[g].append(label)
            vert.append(label)
        elif g == "vertebrate_basal":
            vert.append(label)
    out["_vert"], out["_family"] = vert, family
    return out


def _layer_scores(alignment: dict[str, str], members: list[str],
                  ref_row: str) -> tuple[list[dict], dict[int, int]]:
    sub = {m: alignment[m] for m in members if m in alignment}
    if ref_row not in sub:                       # the reference must be scored
        sub[ref_row] = alignment[ref_row]
    return L.profile(sub), L.residue_to_col(sub[ref_row])


def build(out_dir: Path = L.OUT_DIR) -> dict:
    msa = L.read_fasta(L.MSA_DIR / "aln.fasta")
    subsets = _msa_subsets()
    elem = D.residue_element(out_dir)
    sites = D.site_index(out_dir)

    summary: dict[str, dict] = {}
    all_rows: list[dict] = []

    for paralog in L.PARALOGS:
        msa_label, acc, desc = L.REFERENCES[paralog]
        ref_seq = L.uniprot_fasta(acc)

        deep = L.read_fasta(out_dir / f"aln_{paralog}.fasta")
        deep_ref = f"REF|{paralog}|{acc}"
        deep_stats = L.profile(deep)
        deep_col = L.residue_to_col(deep[deep_ref])

        msa_ref_seq = msa[msa_label].replace("-", "")
        identical = msa_ref_seq == ref_seq
        xfer = L.transfer_positions(msa_ref_seq, ref_seq)
        back = {v: k for k, v in xfer.items()}
        msa_layers = {name: _layer_scores(msa, members, msa_label)
                      for name, members in (("shallow", subsets[paralog]),
                                            ("vert", subsets["_vert"]),
                                            ("family", subsets["_family"]))}

        rows = []
        psites = sites.get(paralog, {})
        for resi in range(1, len(ref_seq) + 1):
            row = {
                "paralog": paralog, "ref_acc": acc, "reference": desc,
                "resi": resi, "aa": ref_seq[resi - 1],
                "element": elem.get(paralog, {}).get(resi, "unassigned"),
                "ip3_contact": resi in psites.get("ip3_contact", set()),
                "filter_lining": resi in psites.get("filter_lining", set()),
                "gate_lining": resi in psites.get("gate_lining", set()),
            }
            c = deep_col.get(resi)
            s = deep_stats[c] if c is not None else None
            row.update({
                "deep_col": c if c is not None else "",
                "deep_jsd": s["jsd"] if s else "",
                "deep_entropy": s["entropy_norm"] if s else "",
                "deep_modal": s["modal"] if s else "",
                "deep_frac_modal": s["frac_modal"] if s else "",
                "deep_occupancy": s["occupancy"] if s else "",
                "deep_n": s["n_seq"] if s else 0,
                "deep_reliable": bool(s and s["occupancy"] >= MIN_OCCUPANCY),
            })
            mres = back.get(resi)
            for name in ("shallow", "vert", "family"):
                stats, col_of = msa_layers[name]
                c2 = col_of.get(mres) if mres else None
                s2 = stats[c2] if c2 is not None else None
                row[f"{name}_jsd"] = s2["jsd"] if s2 else ""
                row[f"{name}_occupancy"] = s2["occupancy"] if s2 else ""
                row[f"{name}_n"] = s2["n_seq"] if s2 else 0
            row["msa_col"] = (msa_layers["vert"][1].get(mres, "") if mres else "")
            rows.append(row)

        L.write_tsv(out_dir / f"constraint_{paralog}_{acc}.tsv", rows)
        all_rows.extend(rows)
        ok = [r for r in rows if r["deep_reliable"]]
        summary[paralog] = {
            "n_residues": len(rows), "n_reliable": len(ok),
            "median_deep_jsd": round(statistics.median(
                [r["deep_jsd"] for r in ok]), 4) if ok else 0.0,
            "n_seq": len(deep), "msa_ref_identical_to_uniprot": identical,
        }
        print(f"[s17] {paralog}: {len(rows)} residues, {len(ok)} with occupancy "
              f">= {MIN_OCCUPANCY:.0%}, median deep JSD "
              f"{summary[paralog]['median_deep_jsd']}, {len(deep)} orthologues"
              + ("" if identical else "  [msa tip != UniProt canonical]"))

    _by_element(out_dir, all_rows)
    _metric_controls(out_dir, all_rows)
    _functional_sites(out_dir, all_rows)
    _paralog_identity(out_dir, msa, elem, sites)
    _layer_sizes(out_dir, subsets, summary)
    return summary


def _layer_sizes(out_dir: Path, subsets: dict[str, list[str]],
                 summary: dict) -> None:
    rows = []
    for paralog in L.PARALOGS:
        rows.append({"layer": "deep", "paralog": paralog,
                     "n_sequences": summary[paralog]["n_seq"],
                     "note": "S5 gene models + curated, one locus per genome"})
        rows.append({"layer": "shallow", "paralog": paralog,
                     "n_sequences": len(subsets[paralog]),
                     "note": "msa_v2, this paralog only — the depth control"})
    rows.append({"layer": "vert", "paralog": "all three",
                 "n_sequences": len(subsets["_vert"]),
                 "note": "msa_v2 vertebrate ITPR tips incl. vertebrate_basal"})
    rows.append({"layer": "family", "paralog": "all",
                 "n_sequences": len(subsets["_family"]),
                 "note": "every msa_v2 ITPR tip; the 6 RyR outgroup tips excluded"})
    L.write_tsv(out_dir / "layer_sizes.tsv", rows)


def _mw(a: list[float], b: list[float]) -> str:
    from scipy.stats import mannwhitneyu
    if len(a) < 5 or len(b) < 5:
        return ""
    return L.fmt(float(mannwhitneyu(a, b, alternative="greater").pvalue))


def _by_element(out_dir: Path, rows: list[dict]) -> None:
    """Per-element constraint, tested against the protein's own linkers.

    The question is not "is the channel conserved" — everything in a 2,700
    residue receptor is somewhat conserved — but whether an element is *more*
    constrained than the same protein's inter-domain sequence, which controls
    for that paralog's own overall level. The elements are also ranked against
    each other, which is what answers the brief's ligand-core-versus-pore
    question.
    """
    out = []
    for paralog in L.PARALOGS:
        pr = [r for r in rows if r["paralog"] == paralog and r["deep_reliable"]]
        ctrl = [r["deep_jsd"] for r in pr if _is_control(r["element"])]
        whole = [r["deep_jsd"] for r in pr]
        for name in sorted({r["element"] for r in pr}):
            vals = [r["deep_jsd"] for r in pr if r["element"] == name]
            if not vals:
                continue
            out.append({
                "paralog": paralog, "element": name, "n_sites": len(vals),
                "is_control": _is_control(name),
                "mean_jsd": round(statistics.mean(vals), 4),
                "median_jsd": round(statistics.median(vals), 4),
                "vs_whole_protein": round(
                    statistics.mean(vals) - statistics.mean(whole), 4),
                "control_mean_jsd": round(statistics.mean(ctrl), 4) if ctrl else "",
                "p_greater_than_linkers": (
                    "" if _is_control(name) else _mw(vals, ctrl)),
                # the composition-free metric, beside the JSD on every row:
                # JSD is measured against a background frequency table, so a
                # conserved Leu column scores lower than a conserved Trp one,
                # and a transmembrane element is made of common residues. A
                # claim that survives in both columns is not a metric artefact.
                "mean_frac_modal": round(statistics.mean(
                    [r["deep_frac_modal"] for r in pr if r["element"] == name]), 4),
                "mean_entropy": round(statistics.mean(
                    [r["deep_entropy"] for r in pr if r["element"] == name]), 4),
                "mean_occupancy": round(statistics.mean(
                    [r["deep_occupancy"] for r in pr if r["element"] == name]), 4),
            })
    L.write_tsv(out_dir / "constraint_by_element.tsv", out)


def _metric_controls(out_dir: Path, rows: list[dict]) -> None:
    """The three things that could make an element's score mean something else.

    Written for every element rather than only where a result is surprising,
    because a control produced after seeing the answer is not one.

    * **composition** — mean BLOSUM62 background frequency of each column's
      modal residue. JSD is a divergence *from that background*, so an element
      built of common amino acids scores lower at equal conservation. Reported
      beside the two composition-free metrics (`frac_modal`, `entropy_norm`).
    * **occupancy** — a gappy element is scored from fewer sequences.
    * **gene models** — the deep alignment is ~95 % miniprot translations, and
      the exon-dense transmembrane region is where a mis-placed boundary would
      do most damage. So every element is recomputed on the **curated subset
      alone** (the UniProt/Ensembl proteins and the reference), which no gene
      caller produced. An element whose score moves between the two columns is
      an alignment result; one that does not is a sequence result.
    """
    out = []
    for paralog in L.PARALOGS:
        acc = L.REFERENCES[paralog][1]
        deep = L.read_fasta(out_dir / f"aln_{paralog}.fasta")
        curated = {k: v for k, v in deep.items()
                   if k.startswith("REF|") or not ("GCF_" in k or "GCA_" in k)}
        cur_stats = L.profile(curated)
        cur_col = L.residue_to_col(curated[f"REF|{paralog}|{acc}"])
        pr = [r for r in rows if r["paralog"] == paralog and r["deep_reliable"]]
        for name in sorted({r["element"] for r in pr}):
            sub = [r for r in pr if r["element"] == name]
            cur = [cur_stats[cur_col[r["resi"]]] for r in sub
                   if r["resi"] in cur_col]
            cur = [c for c in cur if c["occupancy"] >= MIN_OCCUPANCY]
            out.append({
                "paralog": paralog, "element": name, "n_sites": len(sub),
                "is_control": _is_control(name),
                "mean_jsd": round(statistics.mean(
                    [r["deep_jsd"] for r in sub]), 4),
                "mean_frac_modal": round(statistics.mean(
                    [r["deep_frac_modal"] for r in sub]), 4),
                "mean_entropy": round(statistics.mean(
                    [r["deep_entropy"] for r in sub]), 4),
                "mean_occupancy": round(statistics.mean(
                    [r["deep_occupancy"] for r in sub]), 4),
                "mean_background_freq_of_modal": round(statistics.mean(
                    [L.BLOSUM62_BG.get(r["deep_modal"], 0.05) for r in sub]), 5),
                "n_curated": len(curated),
                "n_sites_curated": len(cur),
                "curated_mean_jsd": round(statistics.mean(
                    [c["jsd"] for c in cur]), 4) if cur else "",
                "curated_mean_frac_modal": round(statistics.mean(
                    [c["frac_modal"] for c in cur]), 4) if cur else "",
            })
    L.write_tsv(out_dir / "metric_controls.tsv", out)


def _functional_sites(out_dir: Path, rows: list[dict]) -> None:
    """The measured residues, each against two controls.

    Against the whole protein, and — the sharper one — against the rest of its
    *own* element. A gate residue being more constrained than the average
    residue of a 2,700-aa receptor is nearly guaranteed; being more constrained
    than the rest of the channel domain is not.
    """
    out = []
    for paralog in L.PARALOGS:
        pr = [r for r in rows if r["paralog"] == paralog and r["deep_reliable"]]
        whole = [r["deep_jsd"] for r in pr]
        for flag in ("ip3_contact", "filter_lining", "gate_lining"):
            hit = [r for r in pr if r[flag]]
            if not hit:
                continue
            vals = [r["deep_jsd"] for r in hit]
            elems = sorted({r["element"] for r in hit})
            rest = [r["deep_jsd"] for r in pr
                    if r["element"] in elems and not r[flag]]
            out.append({
                "paralog": paralog, "site_class": flag, "n_sites": len(vals),
                "elements": ";".join(elems),
                "mean_jsd": round(statistics.mean(vals), 4),
                "median_jsd": round(statistics.median(vals), 4),
                "whole_protein_mean": round(statistics.mean(whole), 4),
                "own_element_mean": round(statistics.mean(rest), 4) if rest else "",
                "p_greater_than_protein": _mw(vals, whole),
                "p_greater_than_own_element": _mw(vals, rest),
                "n_invariant": sum(1 for r in hit if r["deep_frac_modal"] == 1.0),
            })
    L.write_tsv(out_dir / "functional_site_constraint.tsv", out)


def _paralog_identity(out_dir: Path, msa: dict[str, str],
                      elem: dict[str, dict[int, str]],
                      sites: dict[str, dict[str, set[int]]]) -> None:
    """Per-element sequence identity *between* the three paralogs.

    A second, independent way to ask which parts the family has kept alike over
    the ~500 Myr since 2R. Identity is measured over mutually covered columns
    only (S6's fragment-aware rule), and the element frame is named on every row
    because the three proteins do not have the same length.

    **On `aln.fasta`, not on `trimmed.fasta`, and that is why S6's number for
    the same pair is higher.** S6 measured its identity matrices on the
    trimAl product, which kept 1,797 of 11,777 columns; S17 cannot, because
    trimAl's job is to remove exactly the divergent and gappy columns, and the
    element this task's headline rests on is the most divergent and gappiest in
    the protein. So each row carries **both**: the identity over all mutually
    covered columns, and the same pair measured S6's way, with the number of
    element residues trimAl retained. A per-element identity read off the
    trimmed alignment would be a measurement of what survived trimming.
    """
    trimmed = L.read_fasta(L.MSA_DIR / "trimmed.fasta")
    out = []
    pairs = [("ITPR1", "ITPR2"), ("ITPR1", "ITPR3"), ("ITPR2", "ITPR3")]
    for a, b in pairs:
        la, aacc, _ = L.REFERENCES[a]
        lb, _bacc, _ = L.REFERENCES[b]
        row_a, row_b = msa[la], msa[lb]
        ref_a = L.uniprot_fasta(aacc)
        fwd = L.transfer_positions(row_a.replace("-", ""), ref_a)
        c2r = L.col_to_residue(row_a)
        buckets: dict[str, list[int]] = {}
        for c, (ca, cb) in enumerate(zip(row_a, row_b)):
            if ca in "-." or cb in "-.":
                continue
            resi = fwd.get(c2r.get(c, -1))
            if resi is None:
                continue
            same = 1 if ca == cb else 0
            buckets.setdefault(elem.get(a, {}).get(resi, "unassigned"),
                               []).append(same)
            for flag, positions in sites.get(a, {}).items():
                if resi in positions:
                    buckets.setdefault(flag, []).append(same)
            buckets.setdefault("_whole", []).append(same)
        whole = buckets.pop("_whole", [])
        wid = sum(whole) / len(whole) if whole else None
        # the same pair, measured the way S6 measured it
        ta, tb = trimmed.get(la, ""), trimmed.get(lb, "")
        tcov = [(x, y) for x, y in zip(ta, tb)
                if x not in "-." and y not in "-."]
        tid = (sum(1 for x, y in tcov if x == y) / len(tcov)) if tcov else None
        # which residues of the frame protein trimAl kept at all
        kept = _trimmed_residues(row_a, ta)
        for name, vals in sorted(buckets.items()):
            if len(vals) < 5:
                continue
            ident = sum(vals) / len(vals)
            in_elem = [r for r, e in elem.get(a, {}).items() if e == name] \
                if name in elem.get(a, {}).values() else \
                sorted(sites.get(a, {}).get(name, set()))
            out.append({
                "pair": f"{a}_vs_{b}", "frame": a, "element": name,
                "n_columns": len(vals), "identity": round(ident, 4),
                "whole_protein_identity": round(wid, 4) if wid else "",
                "delta_vs_whole": round(ident - wid, 4) if wid else "",
                "whole_protein_identity_trimmed": round(tid, 4) if tid else "",
                "n_residues_kept_by_trimal": sum(1 for r in in_elem if r in kept),
                "n_residues_in_element": len(in_elem),
            })
    L.write_tsv(out_dir / "paralog_identity_by_element.tsv", out)


def _trimmed_residues(aln_row: str, trimmed_row: str) -> set[int]:
    """Residues of one sequence that trimAl kept.

    trimAl deletes columns, so the surviving residues are a prefix-preserving
    subsequence: walking both rows in step and matching residue characters
    recovers which of the sequence's own positions survived.
    """
    kept: set[int] = set()
    seq = [c for c in aln_row if c not in "-."]
    tseq = [c for c in trimmed_row if c not in "-."]
    i = 0
    for ch in tseq:
        while i < len(seq) and seq[i] != ch:
            i += 1
        if i < len(seq):
            kept.add(i + 1)
            i += 1
    return kept


def run(out_dir: Path = L.OUT_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    return build(out_dir)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    run(ap.parse_args().out)


if __name__ == "__main__":
    main()
