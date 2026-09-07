"""S7 — the report's membership half: what the tree says tip by tip.

Split out of `s7_report_results.py` to keep every module under the
project's 500-line budget, and split *here* because these three sections
answer one kind of question — not "what is the family's shape" but "is
this particular record what its name says, and do the sets S6 handed
S7 as stress tests survive contact with a tree":

  § the tips whose census label the tree overturns, and the reciprocal
    best hits that settle each one with an instrument that has never
    seen this alignment
  § the cyclostome loci — S6's open question
  § the 3R duplicate pairs — S6's stress test on its own naming

Like `s3_report_d10.py`, it takes the caller's loader, table renderer
and formatter rather than importing its own, so the two halves of the
report cannot read the same table differently.
"""

from __future__ import annotations

from s7_priors import PRIOR


def _sup(r: dict) -> str:
    a, u = r.get("alrt", ""), r.get("ufboot", "")
    return f"{a}/{u}" if (a and u) else (a or u or "—")


def _strong(r: dict, min_alrt=80.0, min_uf=95.0) -> bool:
    try:
        return float(r.get("alrt") or 0) >= min_alrt and \
            float(r.get("ufboot") or 0) >= min_uf
    except (TypeError, ValueError):
        return False

def _hit(title: str) -> str:
    """A BLAST subject title, safe in a markdown cell.

    UniProt titles are `db|ACC|NAME description`, and an unescaped pipe
    silently splits the row into extra columns — the table renders, just
    with every later cell shifted one place left.
    """
    return title.replace("|", "\\|")


# --------------------------------------------------------- the relabels

def _section_relabels(load, table, num, reps) -> list[str]:
    audit = load("membership_audit.tsv")
    conf = load("naming_conflicts.tsv")
    rbh = load("naming_rbh.tsv")
    if not audit:
        return []
    n_re = sum(1 for a in audit if a["rule"] == "reassigned")
    n_free = sum(1 for a in audit if a["rule"] == "unconstrained")
    n_core = sum(1 for a in audit if a["rule"] == "core")
    L = ["### § The tips whose census label the tree overturns", "",
         "Three rules, each a positive test, each recording the number it "
         "fired on (`s7_constraints.py`, `membership_audit.tsv`): a tip "
         "inside its own group's largest pure clade keeps its label "
         "(`core`); a tip outside it whose smallest containing clade holds "
         "at least three members of exactly one *other* paralog, at a node "
         "clearing SH-aLRT ≥ 80 **and** UFBoot ≥ 95, is moved there "
         "(`reassigned`); anything else is left out of all three "
         "constraint sets (`unconstrained`) rather than assumed either "
         "way.", ""]
    L += table(["rule", "tips"],
               [["`core` — label and tree agree", n_core],
                ["`reassigned` — tree overrules the label", n_re],
                ["`unconstrained` — tree cannot place it", n_free]])
    if conf:
        L += ["", "Every tip in the last two rows, with what put it "
              "there:", ""]
        L += table(["tip", "census label", "rule", "tree places it with",
                    "clade size", "SH-aLRT/UFBoot", "why"],
                   [[f"`{c['label'][:52]}`", c["census_group"], c["rule"],
                     c["assigned_to"] or "—", c["clade_size"] or "—",
                     _sup(c), c["why"]] for c in conf])
    if rbh:
        L += ["", "**Each naming call the tree contradicts is settled "
              "independently by reciprocal best hits** — the candidate "
              "against the human reference proteome, then that human "
              "protein back against the candidate's own species. RBH "
              "knows nothing about this alignment, this model or this "
              "tree:", ""]
        # What `agrees` compares depends on what the tree offered — a
        # reassignment to test, or only a refusal to place the tip — so
        # the column is headed by what was actually compared rather than
        # by "agrees with tree", which would be false for every
        # `unconstrained` row (there is no tree call there to agree
        # with). `naming_rbh.tsv` carries the choice per row.
        L += table(["tip", "census label", "tree rule", "forward best hit",
                    "reciprocal", "RBH call", "RBH upholds"],
                   [[f"`{r['label'][:40]}`", r.get("census_group", ""),
                     r.get("tree_rule", ""), _hit(r.get("forward_hit", "")),
                     r.get("reciprocal_ok", ""), r.get("rbh_call", ""),
                     (r.get("agrees", "") + " ("
                      + {"tree_says": "the tree's reassignment",
                         "census_group": "the census name"}.get(
                             r.get("compared_against", ""), "—") + ")")]
                    for r in rbh])
        upheld = sum(1 for r in rbh if r.get("agrees") == "yes")
        free_up = sum(1 for r in rbh if r.get("agrees") == "yes"
                      and r.get("compared_against") == "census_group")
        L += ["", f"**{upheld} of {len(rbh)}** disputed tips are upheld by "
              "RBH. "
              + (f"For {free_up} of them the tree offered no alternative "
                 "placement, only a refusal to place the tip, so what RBH "
                 "upholds there is the **census name** — the record is "
                 "correctly named and the conflict is this tree's own "
                 "uncertainty about where to hang it, not an annotation "
                 "error." if free_up else ""), ""]
    elif conf:
        L += ["", "*RBH verification has not been run — "
              "`python3 scripts/s7_rbh.py`.*", ""]
    return L


# ---------------------------------------------------------- cyclostomes

def p_all(rows: list[dict]) -> str:
    """The clade sizes those sets sit in, as a phrase."""
    sizes = sorted({r.get("mrca_size", "") for r in rows if r.get("mrca_size")})
    return sizes[0] + "-tip" if len(sizes) == 1 else "-/".join(sizes) + "-tip"


def _set_name(s: str) -> str:
    """The row's own name for the table.

    A species carries a common name in brackets that adds nothing here;
    a `cyclostome-only clade (N species)` row carries its *count* there,
    and two such rows stripped of it read as one clade listed twice.
    """
    if s.startswith("cyclostome-only clade"):
        return s.replace(" (", ", ").rstrip(")")
    return s.split(" (")[0]


def _strong_all(rows: list[dict], min_alrt=80.0, min_uf=95.0) -> bool:
    """Every row clears both of IQ-TREE's thresholds.

    A conclusion drawn from a set of clades is only as good as its
    weakest member, so the branch that reads them as separate lineages
    has to be gated on all of them, not on their average.
    """
    return bool(rows) and all(_strong(r, min_alrt, min_uf) for r in rows)


def _section_cyclostome(load, table, num) -> list[str]:
    place = load("cyclostome_placement.tsv")
    mono = load("cyclostome_monophyly.tsv")
    if not place:
        return []
    p = PRIOR["cyclostome_nearest"]
    L = ["### § The cyclostome trio — expansion, or three 1:1 orthologs?",
         "",
         f"Prior: **all six loci nearest {p['value']}** — {p['where']}.", "",
         "The tree can separate the two readings that identity could not. "
         "Three 1:1 orthologs from 2R put each locus *inside* a different "
         "paralog clade; a lineage-specific expansion puts all of them "
         "together, outside all three. A third arrangement — several "
         "cyclostome-only clades, each holding both species — is neither, "
         "and is what the maximal pure-cyclostome pieces below test for.",
         ""]
    # The composition is printed, not just the home. `home` is empty
    # whenever the nearest neighbourhood is mixed, and a column of
    # dashes reads as "the tree says nothing" when what the tree
    # actually says is in the mixture.
    L += table(["locus", "species", "tree places it with", "clade size",
                "what else is in that clade", "SH-aLRT/UFBoot"],
               [[f"`{r['label'][-28:]}`", (r["species"] or "").split(" (")[0],
                 r["home"] or "—", r["clade_size"] or "—",
                 r.get("composition") or "—", _sup(r)] for r in place])
    homes = {r["home"] for r in place if r["home"]}
    all_set = next((r for r in mono if r["set"] == "all cyclostome loci"),
                   None)
    if all_set and all_set["monophyletic"] == "yes":
        L += [f"**All {all_set['n']} cyclostome loci form a single clade** "
              f"({_sup(all_set)}), outside the three vertebrate paralogs. "
              f"That is the expansion answer: the copies are cyclostome "
              f"duplicates, not 2R ohnologs, and the ITPR1 they all leaned "
              f"towards in S6 is the artefact of ITPR1 being the nearest "
              f"*outside* thing, not evidence of orthology. The prior is "
              f"**confirmed as an artefact of the metric** — S6 said "
              f"identity could not tell these apart and it could not.", ""]
    elif len(homes) >= 3:
        L += [f"**The loci resolve into {len(homes)} different paralog "
              f"neighbourhoods ({', '.join(sorted(homes))})** — the 1:1 "
              f"orthology answer, and a direct contradiction of the "
              f"identity ranking, which put every one of them nearest "
              f"{p['value']}. The prior is **contradicted**, which is the "
              f"outcome S6 predicted was possible and could not test.", ""]
    # Neither of S6's two readings, but not therefore nothing: the
    # maximal pure-cyclostome pieces are computed in `s7_analyse.py`
    # precisely because "not one clade" and "scattered across the
    # paralogs" are different results, and a piece holding both species
    # is a duplication older than the hagfish/lamprey split.
    pieces = [r for r in mono if r["set"].startswith("cyclostome-only")
              and r["monophyletic"] == "yes"]
    covered = sum(int(r["n"]) for r in pieces)
    shared = [r for r in pieces if not r["set"].startswith(
        "cyclostome-only clade (1 ")]
    n_all = int(all_set["n"]) if all_set else 0
    if pieces and covered == n_all and shared and _strong_all(pieces):
        L += [f"**Every one of the {n_all} loci sits in a cyclostome-only "
              f"clade** — {len(pieces)} of them, each well supported, and "
              f"{len(shared)} holding *both* hagfish and lamprey. So the "
              f"loci are neither one expansion nor three 1:1 ohnologs: "
              f"they are **{len(pieces)} anciently separate cyclostome "
              f"lineages**, and a clade spanning the two species is a "
              f"duplication older than the hagfish/lamprey split rather "
              f"than a copy either genome made on its own. The prior — "
              f"every locus nearest {p['value']} — is **contradicted as a "
              f"reading**: identity ranked all six against the same "
              f"paralog because it can only measure distance to the "
              f"labelled clades, and these loci are outside all of them. "
              f"Which side of the vertebrate duplication each lineage "
              f"attaches to is the next table, and S8's synteny is what "
              f"would settle it *(pending: S8)*.", ""]
    else:
        where = ", ".join(sorted(homes)) or \
            "nothing at the required support"
        L += [f"The loci neither form one clade nor spread across three "
              f"paralogs — they resolve with {where}. **Unresolved**: "
              f"this alignment does not settle the cyclostome question "
              f"either, and S8's synteny is the next instrument.", ""]
    if mono and len(mono) > 1:
        L += ["Per species, in case the expansion is not shared:", ""]
        L += table(["set", "loci", "form a clade", "SH-aLRT/UFBoot"],
                   [[_set_name(r["set"]), r["n"],
                     "**yes**" if r["monophyletic"] == "yes" else "no",
                     _sup(r)] for r in mono])
    return L


# ------------------------------------------------------- the 3R stress test

def _section_duplicates(load, table, num) -> list[str]:
    rows = load("duplicate_pairs.tsv")
    if not rows:
        return []
    p = PRIOR["duplicate_pairs"]
    ok = sum(1 for r in rows if r["sisters"] == "yes")
    L = ["### § S6's stress test: the same species, the same paralog, twice",
         "",
         f"Prior: **{p['value']}** — {p['where']}.", ""]
    L += table(["paralog", "species", "copies", "sisters on the tree",
                "SH-aLRT/UFBoot", "smallest clade holding them all",
                "what else is in it"],
               [[r["paralog"], (r["species"] or "").split(" (")[0],
                 r["n"], "**yes**" if r["sisters"] == "yes" else "no",
                 _sup(r), r.get("mrca_size", "—"),
                 r.get("mrca_intruders") or "nothing"] for r in rows])
    verdict = "confirmed" if ok == len(rows) else (
        "contradicted" if ok == 0 else "partly confirmed")

    # "Not sisters" is not one finding but two, and the difference is
    # what the copies' smallest containing clade holds. A duplication
    # *older than the species* — which is what teleost 3R is — puts each
    # copy with the other species' matching copy, so same-species copies
    # are expected **not** to be sisters and their MRCA clade is expected
    # to hold other tips of the same paralog and nothing else. That is a
    # different result from copies pulled apart by a paralog or a
    # non-family tip, which is the naming or alignment failure S6 had in
    # mind. Reporting them as one number would score the expected
    # outcome as a failure.
    broken = [r for r in rows if r["sisters"] != "yes"]
    same_par = [r for r in broken
                if r.get("mrca_intruders")
                and all(i.split(":")[0] == r["paralog"]
                        for i in r["mrca_intruders"].split(";") if i)]
    other = [r for r in broken if r not in same_par]
    L += [f"**{ok} of {len(rows)}** same-species same-paralog sets come "
          f"back as clades — against the prior, that reads "
          f"**{verdict}**.", ""]
    if same_par:
        L += [f"But {len(same_par)} of the {len(broken)} that do not are "
              f"broken **only by other tips of the same paralog**: the "
              f"copies sit in a small {p_all(same_par)} clade made "
              f"entirely of that paralog, with each copy nearer another "
              f"species' copy than its own genome's. That is what a "
              f"duplication *older than the species* looks like, and the "
              f"teleost 3R is exactly that — so for these the prior was "
              f"the wrong expectation rather than the tree the wrong "
              f"answer. S6 put them in as a stress test on the naming; "
              f"the naming passes, and what fails is the assumption that "
              f"co-orthologs of a shared duplication should be sisters.",
              ""]
    if other:
        L += [f"**{len(other)} set(s) are broken by something else** — "
              + "; ".join(f"{(r['species'] or '').split(' (')[0]} "
                          f"{r['paralog']} by {r['mrca_intruders'] or 'a non-family tip'}"
                          for r in other)
              + ". That is the failure S6's test was looking for: the "
                "alignment or the naming, and it is named here rather "
                "than counted.", ""]
    return L
