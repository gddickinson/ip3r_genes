"""S9 step 2 — codon alignment (PAL2NAL) + pruned trees for codeml/HyPhy.

Takes S6's protein alignment and the validated CDS set, and emits
everything the model-based tests need:

    results/selection/
      prot_sub.fasta        protein MSA restricted to CDS-backed tips
      cds_sub.fasta         matching CDS (same order, same names)
      codon_aln.fasta       PAL2NAL codon alignment (full columns)
      codon_trimmed.fasta   trimAl -automated1 columns, applied codon-aware
      codon_trimmed.phy     PAML sequential phylip (short tip codes)
      tip_codes.tsv         short code <-> S6 label <-> paralog set
      tree_all.nwk          the S7 topology pruned to these tips
      codon_<set>.fasta/.phy + tree_<set>.nwk   per-paralog subsets
      codon_aln_stats.md

Two guards, both of which have to be here rather than downstream:

  * PAL2NAL is the declared tool, and its output is cross-checked against
    an independent in-house protein->codon mapping. A disagreement aborts.
  * A residue whose codon was masked in step 1 is written as `X` in the
    protein rows too. Otherwise the protein says one thing and the codon
    beneath it says `NNN`, and pal2nal is being asked to reconcile a pair
    that does not agree — which it does by dropping the sequence, quietly.

Run:  python scripts/s9_codon_aln.py
"""

from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_lib import Node, parse_newick, to_newick  # noqa: E402
from s9_cds_lib import (ALN_FASTA, OUT_DIR, PARALOGS, ROOTED_NWK,  # noqa: E402
                        load_reps, read_fasta, tool_bin, write_fasta)
from s9_sets import paralog_of  # noqa: E402

#: Fewest tips a set may have and still be estimated on its own.
MIN_SET_TIPS = 4
#: Fewest curated (non-genome-model) tips before the sensitivity subset is
#: worth writing — below this the comparison is between two small samples
#: and says more about sampling than about the gene models.
MIN_CURATED_TIPS = 5

TAG = {"ITPR1": "T1", "ITPR2": "T2", "ITPR3": "T3", "background": "BG"}


def short_code(label: str, group: str, used: set[str]) -> str:
    """A <= 20-char PAML/HyPhy-safe tip code, unique within the set."""
    parts = label.split("_")
    tag = TAG.get(group, "XX")
    genus = parts[1][:6] if len(parts) > 1 else "sp"
    species = parts[2][:5] if len(parts) > 2 else ""
    base = "".join(c for c in f"{tag}_{genus}_{species}".rstrip("_")
                   if c.isalnum() or c == "_")
    code, n = base, 1
    while code in used:
        n += 1
        code = f"{base}{n}"
    used.add(code)
    return code


def prune(node: Node, keep: set[str]) -> Node | None:
    """Copy of `node` containing only `keep` leaves; unary nodes collapsed."""
    if node.is_leaf:
        return node if node.name in keep else None
    kids = [k for k in (prune(c, keep) for c in node.children) if k is not None]
    if not kids:
        return None
    if len(kids) == 1:
        kids[0].length = (kids[0].length or 0.0) + (node.length or 0.0)
        return kids[0]
    out = Node(name=node.name, length=node.length)
    out.children = kids
    return out


def relabel(node: Node, mapping: dict[str, str]) -> None:
    for n in node.walk():
        if n.is_leaf and n.name in mapping:
            n.name = mapping[n.name]


def strip_gap_columns(rows: dict[str, str]) -> dict[str, str]:
    labels = list(rows)
    keep = [i for i in range(len(rows[labels[0]]))
            if any(rows[l][i] != "-" for l in labels)]
    return {l: "".join(rows[l][i] for i in keep) for l in labels}


def strip_gap_codons(rows: dict[str, str]) -> dict[str, str]:
    """Drop codon columns that are all-gap across the subset — the whole
    triplet, never a single nucleotide, or the reading frame shifts."""
    labels = list(rows)
    ncod = len(rows[labels[0]]) // 3
    keep = [c for c in range(ncod)
            if any(rows[l][c * 3:c * 3 + 3] != "---" for l in labels)]
    return {l: "".join(rows[l][c * 3:c * 3 + 3] for c in keep) for l in labels}


def mask_protein(aln_row: str, cds: str) -> str:
    """Aligned protein row with `X` wherever its codon was masked to NNN."""
    out, k = [], 0
    for aa in aln_row:
        if aa == "-":
            out.append("-")
            continue
        codon = cds[k * 3:k * 3 + 3]
        out.append("X" if "N" in codon.upper() else aa)
        k += 1
    return "".join(out)


def inhouse_codon_map(prot_aln: str, cds: str) -> str:
    """Independent protein-alignment -> codon-alignment mapping."""
    out, k = [], 0
    for aa in prot_aln:
        if aa == "-":
            out.append("---")
        else:
            out.append(cds[k * 3:k * 3 + 3].ljust(3, "N"))
            k += 1
    return "".join(out)


def trimal_columns(prot_fasta: Path, log: Path) -> list[int]:
    proc = subprocess.run(
        [tool_bin("trimal"), "-in", str(prot_fasta), "-automated1",
         "-colnumbering"],
        capture_output=True, text=True)
    log.write_text(proc.stdout + proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"trimal failed: {proc.stderr[:300]}")
    text = proc.stdout
    if "#ColumnsMap" in text:
        text = text.split("#ColumnsMap", 1)[1]
    nums = [t.strip() for t in text.replace("\n", " ").split(",")]
    return [int(t) for t in nums if t.strip().lstrip("-").isdigit()]


def write_phylip(path: Path, rows: dict[str, str]) -> None:
    labels = list(rows)
    n, ln = len(labels), len(rows[labels[0]])
    with open(path, "w") as fh:
        fh.write(f" {n} {ln}\n")
        for lab in labels:
            fh.write(f"{lab}  {rows[lab]}\n")


def subset_rows(labels: list[str], trimmed: dict[str, str],
                code_of: dict[str, str]) -> dict[str, str]:
    """The rows of one selection subset, **in `labels` order**.

    Never through a set. Python hashes strings with a per-process seed, so
    a dict built by walking `{code_of[l] for l in labels}` writes its rows
    in a different order on every run: the alignment is the same alignment
    and the likelihood is the same likelihood, but the file's SHA-256
    changes and the artefact stops being reproducible. S8's T11 caught
    exactly this in a ranked table; here it rewrote a phylip file
    underneath a running codeml job.
    """
    return strip_gap_codons({code_of[l]: trimmed[code_of[l]] for l in labels})


def write_subset(name: str, labels: list[str], trimmed: dict[str, str],
                 code_of: dict[str, str]) -> int:
    """A per-set codon alignment + its pruned tree. Returns codon count."""
    rows = subset_rows(labels, trimmed, code_of)
    write_fasta(OUT_DIR / f"codon_{name}.fasta", list(rows.items()))
    write_phylip(OUT_DIR / f"codon_{name}.phy", rows)
    sub = prune(parse_newick(ROOTED_NWK.read_text()), set(labels))
    relabel(sub, code_of)
    (OUT_DIR / f"tree_{name}.nwk").write_text(
        to_newick(sub, with_support=False) + "\n")
    return len(next(iter(rows.values()))) // 3


def main() -> None:
    # The self-test runs before anything is written (the pattern
    # `s5_bait_screen.self_test()` set): a rule that is not exercised on
    # every build is a rule nobody knows still holds.
    from s9_test_codon import self_test  # noqa: PLC0415  (cyclic at import)
    failed = [n for n, ok, _ in self_test() if not ok]
    if failed:
        raise SystemExit("S9 self-test FAILED: " + ", ".join(failed))
    print()

    cds_all = read_fasta(OUT_DIR / "cds.fasta")
    aln = read_fasta(ALN_FASTA)
    reps = {r["label"]: r for r in load_reps()}
    groups = {k: v["group"] for k, v in reps.items()}
    routes = {}
    with open(OUT_DIR / "cds_status.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            routes[row["label"]] = row["route"]

    tips = [l for l in aln if l in cds_all]
    assign, why = paralog_of(tips)
    # A tip the tree placed in no paralog clade is background: it keeps the
    # whole-tree branch lengths honest and is in no foreground.
    setof = {l: assign.get(l, "background") for l in tips}
    tips.sort(key=lambda l: (setof[l], l))
    print(f"tips with a validated CDS: {len(tips)}")
    for l in tips:
        if setof[l] != "background" and groups.get(l) != setof[l]:
            print(f"  tree assignment: {l} ({groups.get(l)}) -> {setof[l]}")
    n_bg = sum(1 for l in tips if setof[l] == "background")
    print(f"  background (in no paralog clade of the S7 tree): {n_bg}")

    used: set[str] = set()
    code_of = {l: short_code(l, setof[l], used) for l in tips}

    prot_rows = strip_gap_columns(
        {code_of[l]: mask_protein(aln[l], cds_all[l]) for l in tips})
    cds_rows = {code_of[l]: cds_all[l] for l in tips}
    write_fasta(OUT_DIR / "prot_sub.fasta", list(prot_rows.items()))
    write_fasta(OUT_DIR / "cds_sub.fasta", list(cds_rows.items()))

    # ---- PAL2NAL --------------------------------------------------------
    proc = subprocess.run(
        [tool_bin("pal2nal.pl"), str(OUT_DIR / "prot_sub.fasta"),
         str(OUT_DIR / "cds_sub.fasta"), "-output", "fasta", "-codontable", "1"],
        capture_output=True, text=True)
    (OUT_DIR / "pal2nal.log").write_text(proc.stderr[:20000])
    if proc.returncode != 0 or not proc.stdout.strip():
        raise SystemExit(f"pal2nal failed (rc={proc.returncode}): "
                         f"{proc.stderr[:500]}")
    (OUT_DIR / "codon_aln.fasta").write_text(proc.stdout)
    codon = read_fasta(OUT_DIR / "codon_aln.fasta")
    if len(codon) != len(prot_rows):
        raise SystemExit(f"pal2nal returned {len(codon)} of {len(prot_rows)} "
                         "sequences — it dropped some, which is silent")

    # ---- independent cross-check ---------------------------------------
    bad = []
    for code in prot_rows:
        mine = inhouse_codon_map(prot_rows[code], cds_rows[code])
        theirs = codon.get(code, "")
        if len(mine) != len(theirs):
            bad.append(f"{code}: length {len(theirs)} vs {len(mine)}")
            continue
        diff = sum(1 for a, b in zip(mine.upper(), theirs.upper()) if a != b)
        if diff:
            bad.append(f"{code}: {diff} nt differ from the in-house mapping")
    if bad:
        raise SystemExit("pal2nal cross-check FAILED:\n  " + "\n  ".join(bad))
    print(f"pal2nal cross-check passed for {len(prot_rows)} sequences")

    # ---- trimAl columns, applied codon-aware ----------------------------
    cols = trimal_columns(OUT_DIR / "prot_sub.fasta", OUT_DIR / "trimal.log")
    trimmed = {c: "".join(codon[c][i * 3:i * 3 + 3] for i in cols)
               for c in codon}
    write_fasta(OUT_DIR / "codon_trimmed.fasta", list(trimmed.items()))
    write_phylip(OUT_DIR / "codon_trimmed.phy", trimmed)

    # ---- trees -----------------------------------------------------------
    pruned = prune(parse_newick(ROOTED_NWK.read_text()), set(tips))
    relabel(pruned, code_of)
    (OUT_DIR / "tree_all.nwk").write_text(
        to_newick(pruned, with_support=False) + "\n")

    ncod = len(next(iter(codon.values()))) // 3
    stats = ["# S9 codon alignment", "",
             f"- tips with a validated CDS: **{len(tips)}**",
             f"- codon alignment: {ncod} codons",
             f"- trimAl -automated1 kept: **{len(cols)} codons** "
             f"({100 * len(cols) / ncod:.1f} %)",
             f"- tips the S7 tree places in no paralog clade "
             f"(background only): {n_bg}",
             "", "| set | tips | codons (trimmed subset) | note |",
             "|---|---|---|---|"]
    for para in PARALOGS:
        sub = [l for l in tips if setof[l] == para]
        if len(sub) < MIN_SET_TIPS:
            stats.append(f"| {para} | {len(sub)} | — | too few, skipped |")
            continue
        n = write_subset(para, sub, trimmed, code_of)
        n_model = sum(1 for l in sub if routes[l] == "miniprot")
        stats.append(f"| {para} | {len(sub)} | {n} | "
                     f"{n_model} genome gene models |")
    # Sensitivity subsets: the same paralog without the genome gene models.
    # Those carry masked frameshift/stop codons, and a paralog's ω resting
    # partly on them has to be shown not to depend on them.
    for para in PARALOGS:
        sub = [l for l in tips if setof[l] == para]
        cur = [l for l in sub if routes[l] != "miniprot"]
        if len(cur) == len(sub) or len(cur) < MIN_CURATED_TIPS:
            continue
        n = write_subset(f"{para}_curated", cur, trimmed, code_of)
        stats.append(f"| {para} (curated CDS only) | {len(cur)} | {n} | "
                     "sensitivity subset |")

    with open(OUT_DIR / "tip_codes.tsv", "w") as fh:
        fh.write("code\tlabel\tspecies\tclass\tcensus_group\tset\troute"
                 "\twhy\n")
        for l in tips:
            row = reps.get(l, {})
            fh.write(f"{code_of[l]}\t{l}\t{row.get('species', '')}"
                     f"\t{row.get('class', '')}\t{groups.get(l, '?')}"
                     f"\t{setof[l]}\t{routes.get(l, '?')}\t"
                     f"{why.get(l, 'no paralog clade contains it')}\n")
    (OUT_DIR / "codon_aln_stats.md").write_text("\n".join(stats) + "\n")
    print("\n".join(stats))


if __name__ == "__main__":
    main()
