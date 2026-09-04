"""s4_notes.py — render results/genome_manifest_notes.md from the manifest.

D13 applied to the scope document: the denominator's prose is generated from
the rows it describes, so the manifest and the sentence declaring what it
covers cannot drift. Split from `s4_build_manifest.py` to keep both files
short; it is pure formatting over data the driver already has.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

# Assemblies above this get listed individually: at ~1 GB of zip per 3 Gbp
# they dominate the storage plan and S5's per-genome runtime.
BIG_GENOME_GBP = 10.0

REASON_BLURB = {
    "zero_hit_proteome": "a swept reference proteome with no ITPR record at "
                         "all",
    "missing_paralog": "fewer ITPR records than the vertebrate paralog count",
    "fragment_only": "no ITPR record reaching the family's length floor",
    "anchor": "a reference species the project is calibrated on",
}


def _reason_class(reason: str) -> str:
    return reason.split(":")[0]


def write_notes(path: Path, totals: dict, rows: list[dict], order_reps: dict,
                no_order: list[dict], unresolved: list[dict],
                margin_rows: list[dict], margins: list[dict],
                free_gb: float, params: dict) -> None:
    reason_counts: Counter = Counter()
    for r in rows:
        for reason in r["reasons"]:
            reason_counts[_reason_class(reason)] += 1
    classes: Counter = Counter(
        asm.get("vclass") or "?" for asm in order_reps.values())
    margin_classes: Counter = Counter(
        m.get("vclass") or "?" for m in margin_rows)
    by_source = Counter(
        "RefSeq" if r["accession"].startswith("GCF_") else "GenBank"
        for r in rows)
    unannotated = [r for r in rows if not r["annotated"]]
    big = sorted((r for r in rows if r["total_length"] / 1e9 > BIG_GENOME_GBP),
                 key=lambda r: -r["total_length"])

    a = [
        "# Genome manifest notes (S4) — the declared denominator", "",
        "Every absence claim this project makes is a claim about the "
        "assemblies listed in `genome_manifest.tsv`, and about nothing else.",
        "", "## The scope rule", "",
        "**One best reference assembly per vertebrate order, UNION every "
        "species the protein-level census leaves undecided.**", "",
        "Order representatives are ranked *annotated > RefSeq > assembly "
        "level > scaffold N50* (D9: a RefSeq gene set and a submitter "
        "GenBank gene set are not comparable evidence, so the annotation "
        "source is ranked on and then carried into the manifest rather than "
        "discarded).", "",
        "The margin species are **derived from the committed census tables**, "
        "not hand-listed — the rule and its threshold are in "
        "`scripts/s4_manifest_lib.py:derive_margins`, so the denominator can "
        "be rebuilt from the census and cannot drift from it:", "",
        "| Rule | Threshold | Species | What the genome sweep decides |",
        "|---|---|---:|---|",
    ]
    rule_rows = [
        ("`zero_hit_proteome`", "0 ITPR records in the reference proteome",
         "whether the gene is absent or merely unannotated"),
        ("`missing_paralog`",
         f"< {params['expected_paralogs']} ITPR records",
         "whether a paralog is lost or the proteome is shallow"),
        ("`fragment_only`",
         f"longest record < {params['min_length_aa']:,} aa "
         "(`family.MIN_LENGTH_AA`)",
         "whether the gene model is broken or the protein is truly short"),
        ("`anchor`", "the calibration species", "nothing — it is the control"),
    ]
    counts = Counter(_reason_class(reason) for m in margins
                     for reason in m["reasons"])
    for label, thresh, decides in rule_rows:
        key = label.strip("`")
        a.append(f"| {label} | {thresh} | {counts.get(key, 0)} | {decides} |")
    a += [
        "",
        f"A species can fire several rules; the union is **{len(margins)} "
        f"species**, of which {len(margin_rows)} resolved to an assembly.", "",
        "`missing_paralog` and `fragment_only` are *questions, not findings*. "
        "S3 measured copy number tracking annotation depth — bird reference "
        "proteomes hold a median 13,894 proteins against Mammalia's 34,127 — "
        "so a low count is exactly what a genome search exists to resolve. "
        "That is why these species are in the denominator rather than in a "
        "results table.", "",
        "## What the manifest contains", "",
        f"- **{totals['n']} genomes**, from {len(order_reps)} vertebrate "
        f"orders and {len(margin_rows)} margin species (a species that is "
        "both is one row with both reasons).",
        "- Rows per reason: " + ", ".join(
            f"`{k}` {v}" for k, v in sorted(reason_counts.items(),
                                            key=lambda kv: -kv[1])),
        "- Orders per class: " + ", ".join(
            f"{k} {v}" for k, v in sorted(classes.items())),
        "- Margin species per class: " + ", ".join(
            f"{k} {v}" for k, v in sorted(margin_classes.items(),
                                          key=lambda kv: -kv[1])),
        "- Annotation source: " + ", ".join(
            f"{k} {v}" for k, v in sorted(by_source.items())) +
        f"; **{len(unannotated)} carry no gene set at all** (D9 — an "
        "annotation claim may not be quoted for these).",
        f"- Species in the reference dump with no NCBI order rank: "
        f"{len(no_order)} (not order-eligible; the margin rules can still "
        "pull them in).", "",
        "## Storage", "",
        f"- Total sequence: **{totals['total_bp'] / 1e9:,.1f} Gbp**.",
        f"- Uncompressed FASTA ≈ **{totals['total_fasta_gb']:.0f} GB**; "
        f"download as zip ≈ **{totals['total_zip_gb']:.0f} GB** "
        f"(at {params['zip_factor']:.0%} of sequence length).",
        f"- Free on the data root at build time: **{free_gb:,.0f} GB** — "
        + ("sufficient for the whole manifest resident at once."
           if free_gb > totals["total_fasta_gb"] + totals["total_zip_gb"]
           else "**not** sufficient for the whole manifest at once; run "
                "`fetch_genomes.py --delete-after-search`."),
        "",
    ]

    if big:
        a += [f"### Assemblies over {BIG_GENOME_GBP:.0f} Gbp", "",
              "These dominate both the download and S5's per-genome runtime.",
              "", "| Accession | Species | Gbp | Reasons |", "|---|---|---:|---|"]
        for r in big:
            a.append(f"| {r['accession']} | {r['organism']} | "
                     f"{r['total_length'] / 1e9:.1f} | "
                     f"{';'.join(r['reasons'])} |")
        a.append("")

    if unresolved:
        a += ["### Margin species with no assembly at NCBI (excluded)", "",
              "Excluded from the denominator, and therefore outside every "
              "claim it supports.", ""]
        for m in unresolved:
            a.append(f"- **{m['organism']}** (taxid {m['taxid']}, "
                     f"{';'.join(m['reasons'])}): {' | '.join(m['notes'])}")
        a.append("")

    if no_order:
        by_class: dict[str, list[str]] = defaultdict(list)
        for asm in no_order:
            by_class[asm.get("vclass") or "?"].append(asm["organism"])
        a += ["### Reference species without an NCBI order rank", "",
              "NCBI taxonomy places these in *incertae sedis* buckets with no "
              "ranked order, so they cannot represent one. The margin rules "
              "still reach them.", ""]
        for cls, names in sorted(by_class.items()):
            a.append(f"- **{cls}** ({len(names)}): "
                     + ", ".join(sorted(names)[:12])
                     + (" …" if len(names) > 12 else ""))
        a.append("")

    a += ["## Reproducing this", "", "```",
          "python3 scripts/s4_build_manifest.py            # offline from cache",
          "python3 scripts/s4_build_manifest.py --refresh-assemblies",
          "```", "",
          "Raw `datasets` dumps are archived under "
          "`<data_root>/raw_api/ncbi_datasets/`, so the manifest re-derives "
          "offline from the same bytes it was first built from.", ""]
    path.write_text("\n".join(a))
