"""S17 stage 4 — the human variant harvest, and the numbering check it rests on.

Two things happen here, and the second is what turns a table of scores into a
resource.

**Is the per-site constraint score any good?** A per-site score is only a
variant-interpretation resource if it separates variants that matter from
variants that do not. ClinVar supplies that labelling for free —
pathogenic/likely-pathogenic against benign/likely-benign across ITPR1, ITPR2
and ITPR3. The score is then scored as a *classifier* (ROC AUC, Mann-Whitney)
for each of the four conservation layers. That is also the arbitration stage 2
promised: if the ~250-orthologue `deep` layer does not beat the 11-18 sequence
`shallow` one, the extra depth was not worth having and the report has to say
so.

**Do the positions human genetics has labelled survive in the other paralogs?**
This family's three paralogs are 74-79 % identical and their clinical records
are wildly unequal — ITPR1 carries a dominant ataxia burden, ITPR2 one reported
recessive phenotype. Carrying every labelled position onto the other two
paralogs through msa_v2 and asking whether the amino acid survived is what
separates "ITPR2 tolerates these positions" from "nobody has looked at ITPR2".

Numbering is checked, never assumed. ClinVar does not file variants in UniProt
numbering and does not file them all on one transcript, so **every cited
transcript's own translated CDS is fetched, aligned to the canonical, positions
transferred through that alignment, and the reference amino acid then required
to match**. A variant that cannot be placed is dropped and counted, never
renumbered: a silent transcript-numbering mismatch produces a resource that
looks fine and is wrong everywhere.

The harvest also carries a **positive control it can fail**: S0's curated
`disease_sites.tsv` names two variants at a specific residue — ITPR3
p.Thr1424Met and p.Arg2524Cys — and both must come back from ClinVar at that
residue with that reference amino acid, or the harvest is missing known
variants and says so.

    python scripts/s17_variants.py
    python scripts/s17_variants.py --relabel   # re-join scores, do not re-fetch
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402
import s17_variant_tests as T                                  # noqa: E402

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
EMAIL = "george.dickinson@gmail.com"
TOOL = "ip3r_genes"

THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*", "Sec": "U",
}

P_RE = re.compile(r"\(p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})\)")
NMV_RE = re.compile(r"(NM_\d+\.\d+)")

#: The join key S0's curated table uses for a variant it localises to a residue.
CURATED_RE = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")


# ------------------------------------------------------------------ fetching

def _get(url: str, tries: int = 5) -> bytes:
    for i in range(tries):
        try:
            return urllib.request.urlopen(url, timeout=120).read()
        except urllib.error.HTTPError as e:
            if e.code in (429, 500, 502, 503) and i < tries - 1:
                time.sleep(2 ** i)
                continue
            raise
        except Exception:
            if i < tries - 1:
                time.sleep(2 ** i)
                continue
            raise
    raise RuntimeError(url)


def clinvar_ids(gene: str) -> list[str]:
    cache = L.cache_dir() / f"clinvar_ids_{gene}.json"
    if cache.exists():
        return json.loads(cache.read_text())
    term = (f'{gene}[gene] AND "missense variant"[molecular consequence] '
            f'AND "single nucleotide variant"[Type of variation]')
    q = urllib.parse.urlencode({"db": "clinvar", "term": term, "retmode": "json",
                                "retmax": 100000, "email": EMAIL, "tool": TOOL})
    time.sleep(0.4)
    ids = json.loads(_get(EUTILS + "esearch.fcgi?" + q))["esearchresult"]["idlist"]
    cache.write_text(json.dumps(ids))
    return ids


def clinvar_records(gene: str) -> list[dict]:
    """ClinVar missense records, cached on disk in batches of 300."""
    ids = clinvar_ids(gene)
    out: list[dict] = []
    for i in range(0, len(ids), 300):
        cache = L.cache_dir() / f"clinvar_sum_{gene}_{i}.json"
        if not cache.exists():
            q = urllib.parse.urlencode(
                {"db": "clinvar", "retmode": "json",
                 "id": ",".join(ids[i:i + 300]), "email": EMAIL, "tool": TOOL})
            time.sleep(0.4)
            cache.write_bytes(_get(EUTILS + "esummary.fcgi?" + q))
        res = json.loads(cache.read_text()).get("result", {})
        for uid in res.get("uids", []):
            out.append(res[uid])
    return out


def transcript_protein(nm: str) -> str:
    """Translated CDS of a RefSeq transcript (cached)."""
    path = L.cache_dir() / f"{nm}.faa"
    if not path.exists():
        q = urllib.parse.urlencode({"db": "nuccore", "id": nm,
                                    "rettype": "fasta_cds_aa", "retmode": "text",
                                    "email": EMAIL, "tool": TOOL})
        time.sleep(0.4)
        path.write_bytes(_get(EUTILS + "efetch.fcgi?" + q))
    seq, started = [], False
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if started:
                break
            started = True
            continue
        if started:
            seq.append(line.strip())
    return "".join(seq).rstrip("*")


def transcript_map(nm: str, ref_seq: str) -> tuple[dict[int, int], dict]:
    """Transcript-protein position -> UniProt canonical position."""
    prot = transcript_protein(nm)
    if not prot:
        return {}, {"transcript": nm, "protein_len": 0, "n_mapped": 0,
                    "identical_to_canonical": False}
    m = ({i: i for i in range(1, len(prot) + 1)} if prot == ref_seq
         else L.transfer_positions(prot, ref_seq))
    return m, {"transcript": nm, "protein_len": len(prot), "n_mapped": len(m),
               "identical_to_canonical": prot == ref_seq}


# -------------------------------------------------------------------- parsing

def bucket(desc: str) -> str:
    d = (desc or "").lower()
    if "conflict" in d:
        return "conflicting"
    if "pathogenic" in d:
        return "P/LP"
    if "benign" in d:
        return "B/LB"
    if "uncertain" in d:
        return "VUS"
    return "other"


def parse_clinvar(gene: str, ref_seq: str) -> tuple[list[dict], dict, list[dict]]:
    rows = []
    drops = {"no_protein_change": 0, "no_transcript": 0,
             "unmappable_position": 0, "ref_aa_mismatch": 0}
    maps: dict[str, dict[int, int]] = {}
    tstats: dict[str, dict] = {}
    for rec in clinvar_records(gene):
        title = rec.get("title", "")
        nm, m = NMV_RE.search(title), P_RE.search(title)
        if not m:
            drops["no_protein_change"] += 1
            continue
        if not nm:
            drops["no_transcript"] += 1
            continue
        nmv = nm.group(1)
        if nmv not in maps:
            maps[nmv], tstats[nmv] = transcript_map(nmv, ref_seq)
            tstats[nmv].update({"gene": gene, "n_variants": 0, "n_kept": 0})
        tstats[nmv]["n_variants"] += 1
        ref1 = THREE_TO_ONE.get(m.group(1))
        alt1 = THREE_TO_ONE.get(m.group(3))
        tpos = int(m.group(2))
        if not ref1 or not alt1:
            drops["no_protein_change"] += 1
            continue
        pos = maps[nmv].get(tpos)
        if pos is None:
            drops["unmappable_position"] += 1
            continue
        if ref_seq[pos - 1] != ref1:
            drops["ref_aa_mismatch"] += 1
            continue
        tstats[nmv]["n_kept"] += 1
        g = rec.get("germline_classification", {}) or {}
        desc = g.get("description", "") or ""
        traits = "; ".join(t.get("trait_name", "")
                           for s in (g.get("trait_set") or []) for t in [s])[:200]
        rows.append({
            "gene": gene, "source": "clinvar",
            "variation_id": rec.get("accession", ""),
            "transcript": nmv, "transcript_resi": tpos,
            "resi": pos, "ref_aa": ref1, "alt_aa": alt1,
            "classification": desc, "review_status": g.get("review_status", ""),
            "condition": traits, "class_bucket": bucket(desc),
        })
    return rows, drops, list(tstats.values())


def uniprot_bucket(desc: str) -> str:
    """UniProt writes a disease variant as `in <ACRONYM>; ...`.

    Everything else in the block — "found in a patient with...", "uncertain
    significance", polymorphism notes — is not a disease assertion and is not
    counted as one.
    """
    d = (desc or "").strip()
    if "uncertain" in d.lower() or "unknown pathological" in d.lower():
        return "other"
    return "P/LP" if re.match(r"in [A-Z0-9]{2,};?", d) else "other"


def uniprot_variants(gene: str, acc: str, ref_seq: str) -> list[dict]:
    rows = []
    for f in L.uniprot_json(acc, "ft_variant").get("features", []):
        if f["type"] != "Natural variant":
            continue
        loc = f.get("location", {})
        if loc.get("start", {}).get("value") != loc.get("end", {}).get("value"):
            continue
        pos = loc["start"]["value"]
        alt = f.get("alternativeSequence", {}) or {}
        orig = alt.get("originalSequence", "")
        new = (alt.get("alternativeSequences") or [""])[0]
        if len(orig) != 1 or len(new) != 1:
            continue
        if pos < 1 or pos > len(ref_seq) or ref_seq[pos - 1] != orig:
            continue
        desc = f.get("description", "")
        rows.append({
            "gene": gene, "source": "uniprot", "transcript": acc,
            "transcript_resi": pos, "variation_id": f.get("featureId", ""),
            "resi": pos, "ref_aa": orig, "alt_aa": new,
            "classification": desc[:160], "review_status": "uniprot_curated",
            "condition": desc[:160], "class_bucket": uniprot_bucket(desc),
        })
    return rows


def uniprot_topology(gene: str, acc: str) -> list[dict]:
    rows = []
    for f in L.uniprot_json(acc, "ft_transmem,ft_topo_dom").get("features", []):
        if f["type"] not in ("Transmembrane", "Topological domain"):
            continue
        rows.append({"gene": gene, "acc": acc, "feature": f["type"],
                     "start": f["location"]["start"]["value"],
                     "end": f["location"]["end"]["value"],
                     "description": f.get("description", "")})
    return rows


def curated_control(variants: list[dict]) -> list[dict]:
    """S0's literature-anchored variants must be in the harvest.

    A positive control the harvest can fail: `disease_sites.tsv` localises two
    entries to a named residue, and if a query that returns ~1,750 records does
    not contain the two variants the review cites, it is not enumerating what
    it claims to.
    """
    out = []
    for r in L.read_tsv(L.S0_FIG_DIR / "disease_sites.tsv"):
        if r["kind"] != "point":
            continue
        m = CURATED_RE.search(r["label"])
        if not m:
            out.append({"gene": r["gene"], "label": r["label"], "resi": r["where"],
                        "expected_aa": "", "found": "", "n_clinvar_at_resi": "",
                        "verdict": "no residue-level HGVS in the curated label"})
            continue
        aa1, pos, _alt = THREE_TO_ONE.get(m.group(1), ""), int(m.group(2)), m.group(3)
        hits = [v for v in variants if v["gene"] == r["gene"]
                and int(v["resi"]) == pos]
        exact = [v for v in hits if v["ref_aa"] == aa1]
        out.append({
            "gene": r["gene"], "label": r["label"], "resi": pos,
            "expected_aa": aa1, "found": bool(exact),
            "n_clinvar_at_resi": len(hits),
            "sources": ";".join(sorted({v["source"] for v in exact})),
            "verdict": ("recovered" if exact else
                        "position present, reference amino acid differs" if hits
                        else "NOT RECOVERED"),
        })
    return out


# ---------------------------------------------------------------------- driver

def _site_index(out_dir: Path) -> dict[str, dict[int, dict]]:
    return {p: {int(r["resi"]): r for r in L.read_tsv(
        out_dir / f"constraint_{p}_{acc}.tsv")}
        for p, (_l, acc, _d) in L.REFERENCES.items()}


JOIN_COLS = ("element", "ip3_contact", "filter_lining", "gate_lining",
             "deep_jsd", "deep_n", "deep_occupancy", "deep_reliable",
             "shallow_jsd", "vert_jsd", "family_jsd")


def _join(variants: list[dict], sites: dict[str, dict[int, dict]]) -> None:
    for v in variants:
        s = sites[v["gene"]].get(int(v["resi"]), {})
        for k in JOIN_COLS:
            v[k] = s.get(k, "")


def run(out_dir: Path = L.OUT_DIR) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    sites = _site_index(out_dir)

    variants, drop_rows, topo, tstats = [], [], [], []
    for gene, (_label, acc, _desc) in L.REFERENCES.items():
        ref_seq = L.uniprot_fasta(acc)
        cv, drops, ts = parse_clinvar(gene, ref_seq)
        up = uniprot_variants(gene, acc, ref_seq)
        variants.extend(cv + up)
        topo.extend(uniprot_topology(gene, acc))
        tstats.extend(ts)
        drop_rows.append({"gene": gene, "n_kept": len(cv), **drops})
        print(f"[s17] {gene}: {len(cv)} ClinVar missense kept, "
              f"{len(up)} UniProt variants; dropped {drops}")
    L.write_tsv(out_dir / "clinvar_transcripts.tsv", tstats)

    ctrl = curated_control(variants)
    L.write_tsv(out_dir / "curated_variant_control.tsv", ctrl)
    for r in ctrl:
        print(f"[s17] curated control {r['gene']} {r['label'][:40]}: {r['verdict']}")

    _join(variants, sites)
    L.write_tsv(out_dir / "variants.tsv", variants)
    L.write_tsv(out_dir / "variant_parsing.tsv", drop_rows)
    L.write_tsv(out_dir / "uniprot_topology.tsv", topo)

    return _tests(out_dir, variants, sites)


def _tests(out_dir: Path, variants: list[dict],
           sites: dict[str, dict[int, dict]]) -> dict:
    test = T.constraint_test(variants, sites)
    L.write_tsv(out_dir / "variant_constraint_test.tsv", test)
    for r in test:
        if r["contrast"] == "P/LP vs B/LB":
            print(f"[s17] {r['gene']:<7} {r['layer']:<8} AUC {r['auc']} "
                  f"(P/LP {r['n_positive']} vs B/LB {r['n_negative']}, "
                  f"p {r['p_mannwhitney']})")

    _pv, audit = T.paralog_audit(variants, out_dir)
    L.write_tsv(out_dir / "paralog_variant_audit.tsv", audit)
    for r in audit:
        if r["contrast"] == "P/LP vs B/LB":
            print(f"[s17] {r['gene']}->{r['other']} {r['contrast']:<16} "
                  f"{r['n_identical']}/{r['n_positions']} vs "
                  f"{r['comparator_frac']:.0%}, p {r['p_fisher']}")

    L.write_tsv(out_dir / "variant_by_element.tsv",
                T.variants_by_element(variants))
    strat = T.vus_stratification(variants, sites)
    L.write_tsv(out_dir / "vus_stratification.tsv", strat)
    for r in strat:
        if r["layer"] == "family":
            print(f"[s17] {r['gene']} VUS: {r['n_vus']} scored, "
                  f"{r['frac_vus_above_pathogenic_median']:.0%} at or above the "
                  f"P/LP median, {r['frac_vus_below_benign_median']:.0%} at or "
                  f"below the B/LB median (family layer)")
    return {"n_variants": len(variants), "n_tests": len(test)}


def relabel(out_dir: Path = L.OUT_DIR) -> dict:
    """Re-join the constraint columns onto the harvest already committed.

    A change to the domain map must not be allowed to re-harvest ClinVar: the
    labels would then move for two reasons at once — the new boundaries and
    whatever ClinVar has accessioned since — and the two would be
    indistinguishable. Any number that moves under `--relabel` moved because of
    the boundaries.
    """
    sites = _site_index(out_dir)
    variants = L.read_tsv(out_dir / "variants.tsv")
    _join(variants, sites)
    L.write_tsv(out_dir / "variants.tsv", variants)
    print(f"[s17] re-labelled {len(variants)} variants from the current map")
    return _tests(out_dir, variants, sites)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=L.OUT_DIR)
    ap.add_argument("--relabel", action="store_true")
    args = ap.parse_args()
    relabel(args.out) if args.relabel else run(args.out)


if __name__ == "__main__":
    main()
