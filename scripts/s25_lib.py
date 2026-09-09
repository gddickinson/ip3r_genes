"""S25 shared declarations — the thesis: chapters, assignment, paths.

The thesis is the *long* form of this project: the manuscript compresses
109,243 words of committed task reports into 17,807, and almost none of the
reasoning behind the numbers survives that compression. This module declares
what the long form is made of.

Two things are declared here rather than written into prose, because both are
guards:

* **`CHAPTERS`** — the chapter order. A chapter file missing from `thesis/`
  is a build error, as a missing section is in `s14_assemble.py`.
* **`ASSIGNMENT`** — one row per committed results directory, naming the
  chapter it is *primary* in, the rule that placed it there, and the chapter
  it might plausibly have gone to instead with the reason it did not. A
  results directory that exists and is not assigned is a build error: that is
  rule **T6** (a chapter that only exists to hold leftovers is an appendix)
  enforced rather than asserted.

The rules themselves are written out in `thesis/chapter_rules.md`, which was
committed before any chapter prose was written (brief step 1). S26 carves the
same body of results into papers and starts from this table.

Stdlib only, and it imports `s14_lib` for page geometry, TSV/JSON/SHA-256
helpers: the thesis and the paper must measure a figure the same way.
"""

from __future__ import annotations

from pathlib import Path

import s14_lib as lib

ROOT = lib.ROOT
RESULTS = lib.RESULTS
TH = ROOT / "thesis"
TH_FIGS = TH / "figures"
DOCS_FIGS = ROOT / "docs" / "figures"

#: Page geometry is the manuscript's, so a figure drawn for the paper is
#: placed at the width it was drawn at here too (D19).
W_FULL = lib.W_FULL
W_HALF = lib.W_HALF

TITLE = ("A genome-scale census of the inositol 1,4,5-trisphosphate "
         "receptor family")
SUBTITLE = ("Range, origin, retention, constraint and record quality across "
            "503 genomes and 7,691 reference proteomes")
#: The document's author. This thesis was written by Claude; the project it
#: reports was directed by the correspondent named in the front matter, who
#: is also the point of contact, since the author has no address of its own.
AUTHOR = "Claude (Opus 5, Anthropic)"
CORRESPONDENCE = "george.dickinson@gmail.com"

# ------------------------------------------------------------------ chapters

#: (n, filename, short title). `n` is None for unnumbered front/back matter.
#: The order here is the order they are stitched.
CHAPTERS: list[tuple[int | None, str, str]] = [
    (None, "00_frontmatter.md", "Front matter"),
    (1, "01_introduction.md", "The receptor, and the question"),
    (2, "02_baseline.md", "The baseline, and what it could be trusted to say"),
    (3, "03_separation.md", "Two families, one architecture"),
    (4, "04_vertebrate_sweep.md",
     "The vertebrate genomic sweep, and what the search was worth"),
    (4, "04b_search_worth.md",
     "The vertebrate genomic sweep, and what the search was worth "
     "(continued)"),
    (5, "05_eukaryote_range.md", "How far the family reaches"),
    (6, "06_phylogeny.md", "The alignment and the tree"),
    (7, "07_origin.md", "Where ITPR1, ITPR2 and ITPR3 came from"),
    (7, "07b_origin_duplication.md",
     "Where ITPR1, ITPR2 and ITPR3 came from (continued)"),
    (8, "08_gene_architecture.md",
     "The gene: exons, introns, and what a fragment is"),
    (9, "09_retention.md", "Counting a loss that never happened"),
    (10, "10_selection.md", "Selection across the vertebrate family"),
    (11, "11_channel.md", "The channel: shape, constraint and the clinic"),
    (11, "11b_channel_variants.md",
     "The channel: shape, constraint and the clinic (continued)"),
    (12, "12_ligand_site.md", "The part the ryanodine receptors do not share"),
    (13, "13_annotation.md", "An archive that cannot find the gene"),
    (13, "13b_annotation_records.md",
     "An archive that cannot find the gene (continued)"),
    (14, "14_methods.md", "Methods, and the reasoning behind them"),
    (14, "14b_methods_controls.md",
     "Methods, and the reasoning behind them (continued)"),
    (15, "15_discussion.md", "General discussion"),
    (None, "16_appendix_corrections.md", "Appendix A — the correction list"),
    (None, "17_appendix_controls.md", "Appendix B — the control inventory"),
    (None, "18_appendix_sensitivity.md", "Appendix C — sensitivity tables"),
    (None, "19_appendix_references.md", "Appendix D — the reference audit"),
    (None, "20_appendix_assignment.md",
     "Appendix E — chapters, papers and the results tree"),
    (None, "21_references.md", "References"),
]

CHAPTER_FILES = [c[1] for c in CHAPTERS]
#: filename -> chapter number (None for front/back matter).
CHAPTER_NUM = {c[1]: c[0] for c in CHAPTERS}

#: chapter number -> its files, in order. A long chapter is written as two or
#: three files so that no source file exceeds the project's 500-line budget;
#: only the first carries the chapter's `#` heading, and the build treats them
#: as one chapter throughout — one figure sequence, one set of claims.
FILES_OF_CHAPTER: dict[int, list[str]] = {}
for _n, _fn, _t in CHAPTERS:
    if _n is not None:
        FILES_OF_CHAPTER.setdefault(_n, []).append(_fn)

#: chapter number -> its first file, for messages that name one.
FILE_OF_CHAPTER = {n: f[0] for n, f in FILES_OF_CHAPTER.items()}

SECTION_REFERENCES = "21_references.md"

# ---------------------------------------------------------------- assignment

#: The rule vocabulary. Written out in full in `thesis/chapter_rules.md`;
#: named here so every assignment row carries the rule that placed it and a
#: rule cannot be invented at assignment time.
RULES = {
    "T1": "one thread of argument, from question to answer",
    "T2": "the instrument is explained where it is first used",
    "T3": "one sitting: four to twelve figures, readable end to end",
    "T4": "a results directory is primary in exactly one chapter",
    "T5": "a measured dead end belongs with the instrument it was measured on",
    "T6": "a grouping that only holds leftovers is an appendix, not a chapter",
    "T7": "every figure is copied from a committed results directory, "
          "never re-plotted",
}

#: results-tree entry -> (chapter, rule, what it is, where else it could have
#: gone and why it did not). Every directory under `results/` must appear
#: here or in `ASSIGNMENT_EXCLUDED`; `s25_assign.py` fails otherwise.
ASSIGNMENT: dict[str, tuple[int, str, str, str]] = {
    "2026-08-18_091322_itpr1_itpr2_itpr3": (
        2, "T2", "S0 app smoke-test bundle, human ITPR1/2/3",
        "Chapter 14 holds the toolchain, but these bundles are the evidence "
        "for a claim made in Chapter 2 — that the three public databases "
        "answer at all — so they sit with the claim."),
    "2026-08-18_092259_itpr1a_itpr1b_itpr2": (
        2, "T2", "S0 app smoke-test bundle, zebrafish",
        "As above; the zebrafish panel is the second half of the same test."),
    "2026-08-18_123301_itpr1_itpr2_itpr3": (
        2, "T2", "S0 app smoke-test bundle, re-run after the Ensembl fix",
        "As above; it is the after half of a before/after measurement."),
    "2026-08-18_123618_s0_smoke_human": (
        2, "T2", "S0 app smoke-test bundle, human single-gene",
        "As above."),
    "s0_baseline": (
        2, "T1", "the literature audit, the database snapshot and the "
                 "review's figure data",
        "Its figure data also feeds Chapter 1, but the audit — 19 claims "
        "scored, two struck — is a result of its own and Chapter 2 is about "
        "that result."),
    "benchmark_controls": (
        3, "T2", "the S1 positive/negative control benchmark",
        "Chapter 14 collects the negative controls as a body of work; this "
        "one is not a constructed unit test but the measurement that "
        "established the sister-family margin, so it belongs where the "
        "margin is derived."),
    "census_v2": (
        3, "T2", "the InterPro enumeration and the architecture call",
        "Chapter 5 reports the range the census implies; the call itself is "
        "the instrument and is built here."),
    "census_v3": (
        3, "T1", "the profile sweep merged with the architecture call",
        "Chapter 4 uses census v3 to derive its baits; v3 is where the two "
        "instruments are reconciled, which is this chapter's thread."),
    "hmm_sweep": (
        3, "T2", "the two profiles, their seeds, calibration and sweep",
        "Chapter 5 re-uses the same profiles outside the vertebrates; they "
        "are built once, here."),
    "genome_manifest.tsv": (
        4, "T2", "the declared vertebrate assembly denominator (S4)",
        "Chapter 9 rests on this denominator for its retention claim, but a "
        "denominator is declared where it is built."),
    "genome_manifest_notes.md": (
        4, "T2", "the rendered scope document for that manifest",
        "As above."),
    "s5_baits": (
        4, "T2", "the 38-bait panel, its seven rules, screen and intron "
                 "calibration",
        "Chapter 5 builds a second panel by different rules; the two are "
        "separate instruments and each is explained where it is used."),
    "genome_ledger": (
        4, "T1", "the 309-genome sweep ledger and its margins",
        "Chapter 9 counts losses off the same sweep, but through S15's "
        "restated states rather than the ledger's, so the ledger is primary "
        "here."),
    "census_v4": (
        4, "T2", "census v3 plus the sweep's 1,058 vertebrate gene models",
        "Chapter 13 asks how databases hold those loci; the merge itself is "
        "part of the sweep."),
    "methods": (
        4, "T2", "S19 — the search measured rather than estimated",
        "S26's proposal splits S19 across two papers (the retention paper "
        "takes the false-negative rate, the archive paper the contribution "
        "half). T4 forbids that here: the false-negative rate is a property "
        "of the search, so the whole of S19 sits with the search and "
        "Chapters 9 and 13 cite it. This is the largest departure between "
        "the two groupings and is recorded in Appendix E."),
    "s20_sweep": (
        5, "T1", "the eukaryote-wide proteome sweep and the plant/fungal "
                 "chase",
        "Chapter 3 owns the profiles; this chapter owns what they found "
        "outside the vertebrates."),
    "s23_baits": (
        5, "T2", "the non-vertebrate bait panel and its measured length band",
        "Chapter 4 holds the vertebrate panel; B2 replaces the ryanodine "
        "control with the MIR-domain sharer, which makes this a different "
        "instrument and not a reuse."),
    "s23_scope": (
        5, "T1", "the 194-genome non-vertebrate sweep and the absence claims",
        "Chapter 9's absence claims are about vertebrate paralogues; these "
        "are about whole clades, and the two use different bars."),
    "census_v5": (
        5, "T2", "census v4 plus the eukaryote-wide proteome sweep",
        "It is the edition the range claim is read off."),
    "census_v6": (
        5, "T2", "census v5 plus the non-vertebrate genomic sweep",
        "The last edition; every later chapter reads it, none builds it."),
    "msa_v2": (
        6, "T2", "the representative set, the alignment and its trim",
        "Chapters 10, 11 and 12 all stand on `trimmed.fasta`; it is built "
        "once, here."),
    "phylogeny": (
        6, "T1", "the ML tree, the model search and the AU sister test",
        "Chapter 7 dates the duplications the tree implies; the tree itself "
        "is this chapter."),
    "synteny": (
        7, "T2", "flanking-gene orthology and the consensus paralog caller",
        "Chapter 9 uses the caller on trace regions and finds it cannot "
        "reach them; the caller is calibrated here."),
    "reconciliation": (
        7, "T1", "gene-tree/species-tree reconciliation and the dated "
                 "placement",
        "Chapter 6 supplies the gene tree; the species tree is an input "
        "declared here."),
    "duplication": (
        7, "T1", "the 2R and 3R tests, and the copy-number landscape",
        "Chapter 8 measures tandem duplicates as a by-product of gene "
        "structure; the duplication history itself is this chapter."),
    "gene_architecture": (
        8, "T1", "exon and intron structure across the sweep",
        "Chapter 13 asks whether a fragmentary annotation is real; this "
        "chapter is where the exon evidence that answers it is measured."),
    "loss_dynamics": (
        9, "T2", "the character matrix and the loss instrument",
        "Chapter 4's ledger records the same cells; S15 restates every one "
        "of them under a rule chain built so that exactly one state can be "
        "counted as a loss, and that restatement is this chapter."),
    "loss_counts": (
        9, "T1", "the Dollo counts and the sensitivity matrix",
        "It is the second half of the same measurement."),
    "selection": (
        10, "T1", "the codon alignment and every model-based selection test",
        "Chapter 11 reads FEL rates onto the structure; the codon alignment "
        "and the codeml suite are built here."),
    "structures": (
        11, "T2", "the AlphaFold coverage probe and the TM-align panel",
        "Chapter 12 measures the ligand pocket on six depositions; the "
        "structural family call and the panel are built here."),
    "constraint": (
        11, "T1", "per-residue constraint, the variant classifier and the "
                  "painted structures",
        "Chapter 12 re-uses the same per-residue tables for two modules; "
        "they are computed here."),
    "ligand_site": (
        12, "T1", "the ligand core against the pore, and the upstream "
                  "pathway",
        "Chapter 11 owns constraint as a quantity; this chapter owns the one "
        "comparison that is about what the ryanodine receptors do not "
        "share."),
    "annotation_bugs": (
        13, "T1", "two annotation failures proved at the exon",
        "Chapter 8 measures exon structure family-wide; these two cases are "
        "about what the annotation does with it."),
    "expression": (
        13, "T2", "junction-spanning RNA-seq over the recovered loci",
        "It exists to answer a question Chapter 13 raises and no other."),
    "annotation_audit": (
        13, "T1", "how the family is recorded, and by whom",
        "Chapter 4 measures whether the search found the gene; this chapter "
        "measures whether the archive did."),
    "supplementary": (
        14, "T5", "the six supplementary figures and the figure audit",
        "Its figures illustrate Chapters 6, 10, 11 and 12 and are placed "
        "there; its *result* is the audit — 26 findings, 16 legends "
        "corrected — which is a statement about how the figures were made, "
        "so the directory is primary in the methods chapter."),
    "toolchain_manifest.txt": (
        14, "T2", "the recorded version of every external tool",
        "Every chapter shells out to something in it; it is stated once."),
}

#: Entries under `results/` deliberately not assigned, with the reason.
ASSIGNMENT_EXCLUDED: dict[str, str] = {
    "session_live.json": "the dashboard's live progress panel — transient "
                         "session state, rewritten by whichever driver is "
                         "running, and not a result",
}


def results_entries() -> list[str]:
    """Every entry directly under `results/`, as the assignment sees them."""
    return sorted(p.name for p in RESULTS.iterdir()
                  if not p.name.startswith("."))


def chapter_title(n: int) -> str:
    for num, _fn, title in CHAPTERS:
        if num == n:
            return title
    raise KeyError(n)


#: The appendices are the methods chapter's back matter: each one is a table
#: the methods chapter refers to, written out. A number stated in Appendix B
#: or E is stated by that material, so the claims ledger's appearance check
#: reads them as part of chapter 14 rather than requiring the chapter body to
#: repeat every number its own appendices carry. They are listed explicitly so
#: that this is a decision and not a side effect of the file naming.
APPENDIX_FILES_OF_CHAPTER: dict[int, list[str]] = {
    14: ["16_appendix_corrections.md", "17_appendix_controls.md",
         "18_appendix_sensitivity.md", "19_appendix_references.md",
         "20_appendix_assignment.md"],
}


def chapter_text(n: int) -> str:
    """Every file of a chapter, concatenated. '' if none is written yet."""
    names = (FILES_OF_CHAPTER.get(n, [])
             + APPENDIX_FILES_OF_CHAPTER.get(n, []))
    return "\n".join((TH / f).read_text(encoding="utf-8")
                      for f in names if (TH / f).exists())


def read_tsv(path: Path) -> list[dict[str, str]]:
    return lib.read_tsv(path)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    lib.write_tsv(path, rows, fields)
