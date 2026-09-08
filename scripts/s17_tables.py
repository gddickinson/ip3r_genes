"""S17's stats file and the numbers the report is not allowed to compute itself.

`constraint_stats.json` records the rules, the parameters, the self-test status
and the SHA-256 of every committed table (S6/S8/S9/S10's discipline), and
`headline()` returns the load-bearing numbers read off those tables, so the
report renders them rather than deriving them a second way (D13).
"""

from __future__ import annotations

import hashlib
import json
import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_conservation as C                                   # noqa: E402
import s17_domains as D                                        # noqa: E402
import s17_lib as L                                            # noqa: E402
import s17_orthologs as O                                      # noqa: E402
import s17_variant_tests as VT                                 # noqa: E402

TABLES = [
    "domain_map.tsv", "functional_sites.tsv", "orthologs_manifest.tsv",
    "ortholog_shape.tsv", "layer_sizes.tsv", "constraint_by_element.tsv",
    "metric_controls.tsv", "functional_site_constraint.tsv",
    "paralog_identity_by_element.tsv", "variants.tsv", "variant_parsing.tsv",
    "clinvar_transcripts.tsv", "curated_variant_control.tsv",
    "variant_constraint_test.tsv", "paralog_variant_audit.tsv",
    "paralog_variant_positions.tsv", "variant_by_element.tsv",
    "vus_stratification.tsv", "uniprot_topology.tsv", "fel_sites.tsv",
    "fel_status.tsv", "selection_by_element.tsv", "painted_structures.tsv",
]

#: Per-paralog constraint tables, named from the reference accessions so the
#: manifest cannot drift from `REFERENCES`.
TABLES += [f"constraint_{p}_{a}.tsv" for p, (_l, a, _d) in L.REFERENCES.items()]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _get(rows: list[dict], **eq) -> dict | None:
    for r in rows:
        if all(str(r.get(k, "")) == str(v) for k, v in eq.items()):
            return r
    return None


def headline(out_dir: Path = L.OUT_DIR) -> dict:
    """The numbers the report quotes, each read off one committed table."""
    h: dict = {}
    layers = L.read_tsv(out_dir / "layer_sizes.tsv")
    h["layer_sizes"] = {(r["layer"], r["paralog"]): int(r["n_sequences"])
                        for r in layers}
    h["deep_n"] = {r["paralog"]: int(r["n_sequences"]) for r in layers
                   if r["layer"] == "deep"}

    elem = L.read_tsv(out_dir / "constraint_by_element.tsv")
    h["elements"] = elem
    h["element_rank"] = {
        p: sorted([(float(r["mean_jsd"]), r["element"]) for r in elem
                   if r["paralog"] == p and r["is_control"] == "False"],
                  reverse=True)
        for p in L.PARALOGS}

    sites = L.read_tsv(out_dir / "functional_site_constraint.tsv")
    h["functional_sites"] = sites

    ident = L.read_tsv(out_dir / "paralog_identity_by_element.tsv")
    h["gate_identity"] = sorted({float(r["identity"]) for r in ident
                                 if r["element"] == "gate"})
    h["loop_identity"] = sorted(float(r["identity"]) for r in ident
                                if r["element"] == "luminal_loop")
    h["whole_identity"] = sorted({float(r["whole_protein_identity"])
                                  for r in ident
                                  if r["whole_protein_identity"] != ""})
    h["whole_identity_trimmed"] = sorted(
        {float(r["whole_protein_identity_trimmed"]) for r in ident
         if r["whole_protein_identity_trimmed"] != ""})
    # S6's own committed matrix, read rather than quoted, so the check in
    # §6.2 compares against the file and not against a number typed here.
    import csv as _csv
    with (L.MSA_DIR / "identity_covered.tsv").open() as fh:
        m = list(_csv.reader(fh, delimiter="\t"))
    col = {n: i for i, n in enumerate(m[0])}
    labels = [L.REFERENCES[q][0] for q in L.PARALOGS]
    s6 = set()
    for row in m[1:]:
        if row[0] in labels:
            for other in labels:
                if other != row[0] and other in col:
                    s6.add(round(float(row[col[other]]), 4))
    h["s6_identity_covered"] = sorted(s6)
    h["trimal_retention"] = {
        r["element"]: (int(r["n_residues_kept_by_trimal"]),
                       int(r["n_residues_in_element"]))
        for r in ident if r["frame"] == "ITPR1"}

    tests = L.read_tsv(out_dir / "variant_constraint_test.tsv")
    h["auc_matched"] = {
        r["layer"]: (r["auc"], r["p_mannwhitney"], int(r["n_positive"]),
                     int(r["n_negative"]))
        for r in tests if r["gene"] == "POOLED"
        and r["contrast"] == "P/LP vs B/LB, all layers scorable"}
    h["auc_protein"] = {r["layer"]: r["auc"] for r in tests
                        if r["gene"] == "POOLED"
                        and r["contrast"] == "P/LP vs whole protein"}
    h["best_layer"] = max(h["auc_matched"], key=lambda k: float(
        h["auc_matched"][k][0] or 0)) if h["auc_matched"] else ""

    variants = L.read_tsv(out_dir / "variants.tsv")
    buckets: dict[tuple[str, str], int] = {}
    for v in variants:
        if v["source"] != "clinvar":
            continue
        buckets[(v["gene"], v["class_bucket"])] = \
            buckets.get((v["gene"], v["class_bucket"]), 0) + 1
    h["buckets"] = buckets
    h["n_clinvar"] = sum(1 for v in variants if v["source"] == "clinvar")
    h["n_vus"] = sum(n for (_g, b), n in buckets.items() if b == "VUS")
    h["frac_vus"] = round(h["n_vus"] / h["n_clinvar"], 4) if h["n_clinvar"] else 0

    parsing = L.read_tsv(out_dir / "variant_parsing.tsv")
    h["n_dropped"] = sum(int(r[k]) for r in parsing for k in
                         ("no_protein_change", "no_transcript",
                          "unmappable_position", "ref_aa_mismatch"))
    h["transcripts"] = L.read_tsv(out_dir / "clinvar_transcripts.tsv")
    h["curated_control"] = L.read_tsv(out_dir / "curated_variant_control.tsv")
    h["paralog_audit"] = L.read_tsv(out_dir / "paralog_variant_audit.tsv")
    h["vus"] = L.read_tsv(out_dir / "vus_stratification.tsv")

    # The two residue-level variants S0 could cite, with their own scores —
    # so the report quotes a measured number rather than pointing at a file.
    cited = {}
    for r in h["curated_control"]:
        if r["verdict"] != "recovered":
            continue
        acc = L.REFERENCES[r["gene"]][1]
        row = _get(L.read_tsv(out_dir / f"constraint_{r['gene']}_{acc}.tsv"),
                   resi=r["resi"])
        if row:
            cited[(r["gene"], int(r["resi"]))] = row
    h["cited_variants"] = cited

    shape = L.read_tsv(out_dir / "ortholog_shape.tsv")
    h["n_shape_dropped"] = sum(1 for r in shape if r["kept"] == "False")
    h["shape_bar"] = shape[0]["bar"] if shape else ""
    h["shape_dropped"] = [r for r in shape if r["kept"] == "False"]

    algn = out_dir / "align_stats.json"
    h["alignments"] = (json.loads(algn.read_text()).get("alignments", [])
                       if algn.exists() else [])

    dm = L.read_tsv(out_dir / "domain_map.tsv")
    h["domain_map"] = dm
    h["transfer_checks"] = {(r["paralog"], r["element"]): r["check"]
                            for r in dm if r["kind"] == "structural"}

    for name, key in (("fel_status.tsv", "fel_status"),
                      ("selection_by_element.tsv", "fel_elements"),
                      ("painted_structures.tsv", "painted")):
        path = out_dir / name
        h[key] = L.read_tsv(path) if path.exists() else []
    return h


def build(out_dir: Path = L.OUT_DIR) -> dict:
    import s17_test_constraint as T
    stats = {
        "task": "S17",
        "references": {p: a for p, (_l, a, _d) in L.REFERENCES.items()},
        "structure_reference": L.STRUCTURE_REF,
        "rules": {
            "orthologs": {
                "min_coverage": O.MIN_COVERAGE,
                "min_identity": O.min_identity(),
                "min_aligned_aa": O.MIN_ALIGNED_AA,
                "one_locus_per_genome_x_paralog": True,
                "s15a_lesion_exclusion": True,
            },
            "conservation": {
                "metric": "Jensen-Shannon divergence vs BLOSUM62 background "
                          "(Capra & Singh 2007), sequence-weighted "
                          "(Henikoff & Henikoff 1994)",
                "min_occupancy": C.MIN_OCCUPANCY,
                "layers": list(C.LAYERS),
                "within_protein_control": "linkers and termini only",
            },
            "domains": {
                "pfam_source": "InterPro coordinates measured per accession "
                               "(S0 domain_coords.tsv) — no transfer",
                "structural_source": "PDB 6DQN measured in S0 — transferred "
                                     "with an anchor test",
                "primary_order": list(D.PRIMARY_ORDER),
            },
            "variants": {
                "sources": ["ClinVar missense SNV", "UniProt natural variant"],
                "numbering": "every cited transcript's own translated CDS "
                             "aligned to the canonical; the reference amino "
                             "acid must match or the variant is dropped",
                "layers_tested": list(VT.LAYERS),
            },
            "fel": {"tool": "HyPhy FEL", "fdr": "Benjamini-Hochberg",
                    "alpha_usable_max": 100.0},
        },
        "self_test": {"n_tests": len(T.TESTS),
                      "passed": T.self_test(verbose=False)},
        "tables": {},
    }
    for name in TABLES:
        path = out_dir / name
        if path.exists():
            stats["tables"][name] = {"sha256": _sha(path),
                                     "bytes": path.stat().st_size}
    (out_dir / "constraint_stats.json").write_text(json.dumps(stats, indent=2))
    return stats


def run(out_dir: Path = L.OUT_DIR) -> dict:
    s = build(out_dir)
    print(f"[s17] stats: {len(s['tables'])} tables hashed, self-test "
          f"{'passed' if s['self_test']['passed'] else 'FAILED'}")
    return s


def main() -> None:
    run()


if __name__ == "__main__":
    main()
