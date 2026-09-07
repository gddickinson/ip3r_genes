"""S9 — the codeml job list, and what each job is asking.

Every model-based test S9 makes is one `CodemlJob` here, so the suite can
be read in one place and a job cannot exist without a stated question.

    m0_<p>              one-ratio ω within paralog p — "is this a gene
                        under purifying selection", the headline per cell
    m0_<p>_curated      the same without the genome gene models, so the
                        answer can be shown not to rest on masked codons
    pair_<p>            runmode -2 pairwise ML dN/dS; also the dS
                        saturation diagnostic that qualifies every
                        cross-paralog comparison
    m0_all              one-ratio over the whole vertebrate family tree —
                        the null for every branch test
    two_ratio_<p>       paralog p's whole clade as foreground ($1)
    bs_null_<p>         branch-site model A null on p's stem (#1), ω fixed
    bs_alt_<p>_w<x>     the alternative, restarted from initial ω = x
    m1a/m2a, m7/m8_<p>  site models within p

**Branch-site model A is restarted by construction, not by patch.** A
nested alternative cannot have a lower optimum than its own null, yet
codeml reaches one routinely on alignments this size — the PIEZO project
hit it and had to add restarts afterwards. Running the alternative from
four initial ω from the start makes "the best of several optima" the
reported number rather than a repair, and the spread across restarts is
itself committed, so a run that is still stuck is visible.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from scripts.s7_lib import parse_newick  # noqa: E402
from s9_cds_lib import OUT_DIR, PARALOGS  # noqa: E402
from s9_codeml_lib import (CodemlJob, labelled_newick,  # noqa: E402
                           newick_no_lengths, unroot)

CODEML_DIR = OUT_DIR / "codeml"

#: Initial ω for the branch-site alternative. codeml's own documented
#: remedy for a local optimum is several starting points; 1.5 is the value
#: the manual suggests and the other three bracket it on both sides.
BS_INIT_OMEGAS = (0.5, 1.5, 2.5, 4.0)

#: Site models. M7/M8 use ten beta categories and are the slow pair.
SITE_MODELS_FAST = (("m1a", 1, 3), ("m2a", 2, 3))
SITE_MODELS_SLOW = (("m7", 7, 10), ("m8", 8, 10))


def tip_sets() -> dict[str, str]:
    """short code -> selection set, from the committed tip table."""
    out: dict[str, str] = {}
    with open(OUT_DIR / "tip_codes.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            out[row["code"]] = row["set"]
    return out


def write_tree(path: Path, text: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text + "\n")
    return path


def _paralog_tree(para: str, suffix: str = ""):
    seq = OUT_DIR / f"codon_{para}{suffix}.phy"
    nwk = OUT_DIR / f"tree_{para}{suffix}.nwk"
    if not seq.exists() or not nwk.exists():
        return None, None
    return seq, unroot(parse_newick(nwk.read_text()))


def build_jobs(skip_slow: bool = False) -> list[CodemlJob]:
    jobs: list[CodemlJob] = []
    sets = tip_sets()

    # ---- per-paralog one-ratio, pairwise, sensitivity, site models ------
    for para in PARALOGS:
        seq, tree = _paralog_tree(para)
        if seq is None:
            continue
        for name, settings, desc in (
            (f"m0_{para}", {"model": 0, "NSsites": 0, "omega": 0.4},
             f"one-ratio ω within {para}"),
        ):
            wd = CODEML_DIR / name
            jobs.append(CodemlJob(
                name=name, seqfile=seq,
                treefile=write_tree(wd / "tree.nwk", newick_no_lengths(tree)),
                workdir=wd, settings=settings, description=desc))
        wd = CODEML_DIR / f"pair_{para}"
        jobs.append(CodemlJob(
            name=f"pair_{para}", seqfile=seq,
            treefile=write_tree(wd / "tree.nwk", newick_no_lengths(tree)),
            workdir=wd,
            settings={"runmode": -2, "model": 0, "NSsites": 0, "omega": 0.4},
            pairwise=True,
            description=f"pairwise ML dN/dS within {para}"))

        cseq, ctree = _paralog_tree(para, "_curated")
        if cseq is not None:
            wd = CODEML_DIR / f"m0_{para}_curated"
            jobs.append(CodemlJob(
                name=f"m0_{para}_curated", seqfile=cseq,
                treefile=write_tree(wd / "tree.nwk", newick_no_lengths(ctree)),
                workdir=wd, settings={"model": 0, "NSsites": 0, "omega": 0.4},
                description=f"one-ratio ω within {para}, curated CDS only "
                            "(no genome gene models)"))

        models = list(SITE_MODELS_FAST) + ([] if skip_slow
                                           else list(SITE_MODELS_SLOW))
        for name, ns, ncat in models:
            wd = CODEML_DIR / f"{name}_{para}"
            jobs.append(CodemlJob(
                name=f"{name}_{para}", seqfile=seq,
                treefile=write_tree(wd / "tree.nwk", newick_no_lengths(tree)),
                workdir=wd,
                settings={"model": 0, "NSsites": ns, "omega": 0.4,
                          "ncatG": ncat},
                description=f"site model {name.upper()} within {para}"))

    # ---- whole-tree branch models ---------------------------------------
    seq_all = OUT_DIR / "codon_trimmed.phy"
    tree_all = unroot(parse_newick((OUT_DIR / "tree_all.nwk").read_text()))
    wd = CODEML_DIR / "m0_all"
    jobs.append(CodemlJob(
        name="m0_all", seqfile=seq_all,
        treefile=write_tree(wd / "tree.nwk", newick_no_lengths(tree_all)),
        workdir=wd, settings={"model": 0, "NSsites": 0, "omega": 0.4},
        description="one-ratio over the whole vertebrate family tree "
                    "(the branch-model null)"))

    for para in PARALOGS:
        fg = {c for c, p in sets.items() if p == para}
        if not fg:
            continue
        wd = CODEML_DIR / f"two_ratio_{para}"
        jobs.append(CodemlJob(
            name=f"two_ratio_{para}", seqfile=seq_all,
            treefile=write_tree(wd / "tree.nwk",
                                labelled_newick(tree_all, fg, "clade")),
            workdir=wd, settings={"model": 2, "NSsites": 0, "omega": 0.4},
            description=f"two-ratio: the whole {para} clade as foreground"))

        wd = CODEML_DIR / f"bs_null_{para}"
        jobs.append(CodemlJob(
            name=f"bs_null_{para}", seqfile=seq_all,
            treefile=write_tree(wd / "tree.nwk",
                                labelled_newick(tree_all, fg, "stem")),
            workdir=wd,
            settings={"model": 2, "NSsites": 2, "fix_omega": 1, "omega": 1},
            description=f"branch-site model A null on the {para} stem"))
        for w0 in BS_INIT_OMEGAS:
            name = f"bs_alt_{para}_w{w0:g}"
            wd = CODEML_DIR / name
            jobs.append(CodemlJob(
                name=name, seqfile=seq_all,
                treefile=write_tree(wd / "tree.nwk",
                                    labelled_newick(tree_all, fg, "stem")),
                workdir=wd,
                settings={"model": 2, "NSsites": 2, "fix_omega": 0,
                          "omega": w0},
                description=f"branch-site model A on the {para} stem, "
                            f"initial ω = {w0:g}"))
    return jobs
