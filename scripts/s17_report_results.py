"""The results half of S17's report, split to keep both under 500 lines.

Takes the caller's `headline()` dict so the two halves cannot read the tables
differently (the `s3_report.py` / `s3_report_d10.py` pattern).

This half carries the architecture results — the element ranking and the ligand
site. The variant half, the per-site selection, the structures and the caveats
are in `s17_report_variants.py`, so both stay inside the 500-line budget.

**Headlines chosen by the data.** Every prior is stated in `s17_priors.PRIOR`
with where the earlier task said it, computed on S17's own tables, and rendered
with both numbers printed either way. The comparison that has to be sayable is
the one that goes badly — here, that the deep within-paralog layer S17 built
its instrument around is *not* the best classifier, and the variant half says
so with the number.
"""

from __future__ import annotations

import statistics as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import s17_lib as L                                            # noqa: E402
import s17_priors as P                                         # noqa: E402


def A(*parts: str) -> str:
    return "\n".join(parts)


def _f(x, nd=3) -> str:
    try:
        return f"{float(x):.{nd}f}"
    except (TypeError, ValueError):
        return str(x) if x not in (None, "") else "—"


def _el(h: dict, paralog: str, element: str, key: str):
    for r in h["elements"]:
        if r["paralog"] == paralog and r["element"] == element:
            return r.get(key, "")
    return ""


def s_elements(h: dict) -> str:
    """§5 — the element ranking, and the element that was hiding inside one."""
    rows = []
    order = ["gate", "selectivity_filter", "RIH_assoc", "RIH_N",
             "nterm_trefoil", "RIH_C", "channel", "MIR", "luminal_loop"]
    for e in order:
        cells = []
        for p in L.PARALOGS:
            j, f = _el(h, p, e, "mean_jsd"), _el(h, p, e, "mean_frac_modal")
            cells.append(f"{_f(j)} / {_f(f)}" if j != "" else "—")
        pv = [_el(h, p, e, "p_greater_than_linkers") for p in L.PARALOGS]
        n = _el(h, "ITPR1", e, "n_sites")
        rows.append(f"| `{e}` | {n} | " + " | ".join(cells) + " | "
                    + "; ".join(x for x in pv if x) + " |")
    ctrl = [float(r["mean_jsd"]) for r in h["elements"]
            if r["is_control"] == "True"]
    ctrl_f = [float(r["mean_frac_modal"]) for r in h["elements"]
              if r["is_control"] == "True"]
    loop = {p: _el(h, p, "luminal_loop", "mean_jsd") for p in L.PARALOGS}
    ch = {p: _el(h, p, "channel", "mean_jsd") for p in L.PARALOGS}

    return A(
        "## 5. Which parts of the receptor are intolerant",
        "",
        "Every element is tested against the **same protein's own linkers** — "
        "the inter-domain sequence and the two termini — rather than against "
        "its whole-protein mean, which would contain the element being tested. "
        f"Pooled linker means: JSD {_f(st.mean(ctrl))}, modal-residue fraction "
        f"{_f(st.mean(ctrl_f))}.",
        "",
        "| element | n (ITPR1) | ITPR1 | ITPR2 | ITPR3 | p vs linkers |",
        "|---|---|---|---|---|---|",
        "| | | *JSD / modal fraction* | | | |",
        *rows,
        "",
        "**The gate and the selectivity filter are the most constrained "
        "elements of the protein, on both metrics and in all three "
        "paralogs.** RIH-associated is the most constrained *domain*. Every "
        "named element except one sits above its own protein's linker mean.",
        "",
        "### 5.1 The exception is inside the channel, and finding it changed §5",
        "",
        f"`luminal_loop` scores JSD {_f(loop['ITPR1'])} / {_f(loop['ITPR2'])} "
        f"/ {_f(loop['ITPR3'])} against a linker mean of {_f(st.mean(ctrl))} — "
        "by a wide margin the least conserved element in the receptor, on the "
        "composition-free metric as well as the JSD.",
        "",
        "Before it was separated out, the Pfam channel domain read as the "
        "*least* constrained named element in ITPR1, below the linkers "
        "(JSD 0.697, p = 0.96 for the one-sided test) — which is a strange "
        "thing to report about the pore of an ion channel, and the sort of "
        "result a report should not print without asking what else could "
        "produce it. Three explanations were checked and all three were "
        "measured rather than argued:",
        "",
        "1. **A metric artefact.** JSD is a divergence from a background "
        "frequency table and a TM domain is built of common residues. Ruled "
        "out: the composition-free modal fraction gave the channel 0.824 "
        "against the linkers' 0.824–0.835, the same picture.",
        "2. **A gene-model artefact.** The deep alignment is ~95 % miniprot "
        "translations and the exon-dense TM region is where a mis-placed "
        "boundary would do most damage. Ruled out: recomputed on the curated "
        "subset alone, the channel still sat below the linkers.",
        "3. **An unresolved element.** Confirmed. The channel's mean was the "
        "average of the most conserved stretch in the protein and the least, "
        "and a bar of that mean is a number no residue has.",
        "",
        f"With the loop separated, `channel` is {_f(ch['ITPR1'])} / "
        f"{_f(ch['ITPR2'])} / {_f(ch['ITPR3'])} and clears the linker control "
        f"(p = {_el(h, 'ITPR1', 'channel', 'p_greater_than_linkers')} in "
        "ITPR1). **The pore module is more constrained than the linkers; the "
        "50 residues of luminal loop inside it are the most variable "
        "sequence in the receptor.** Those two live about fifty residues "
        "apart in the same Pfam domain.",
        "",
        "This is what the geometric definition in §2 was for. Had the loop "
        "boundary been drawn on the conservation profile, this section would "
        "be circular; drawn on the membrane's own axial span in 6DQN, it is "
        "an independent prediction that landed on the dip.",
    )


def s_ligand(h: dict) -> str:
    """§6 — the ligand site and the pore, which is the brief's own question."""
    rows = []
    for r in h["functional_sites"]:
        rows.append(
            f"| {r['paralog']} | `{r['site_class']}` | {r['n_sites']} | "
            f"{r['elements']} | {_f(r['mean_jsd'])} | "
            f"{_f(r['whole_protein_mean'])} | {_f(r['own_element_mean'])} | "
            f"{r['p_greater_than_protein'] or '—'} | "
            f"{r['p_greater_than_own_element'] or '—'} | {r['n_invariant']} |")
    ident = h["gate_identity"]
    loop_id = h["loop_identity"]
    whole_id = h["whole_identity"]
    trimmed_id = h["whole_identity_trimmed"]
    s6_id = h["s6_identity_covered"]
    ret = h["trimal_retention"]
    loop_keep = ret.get("luminal_loop", (0, 0))
    gate_keep = ret.get("gate", (0, 0))
    filt_keep = ret.get("selectivity_filter", (0, 0))
    return A(
        "## 6. The ligand site against the pore",
        "",
        "This is the brief's question, and §2.1 changed how it can be asked: "
        "the Pfam signature named for IP₃ binding contains none of the "
        "measured contacts, so the ligand site here is the **ten residues "
        "within 4.5 Å of IP₃ in 6DQN**, not a domain.",
        "",
        "Each site class is tested against two controls. Against the whole "
        "protein, and — the sharper one — against the rest of its **own** "
        "element. A gate residue beating the average residue of a "
        "2,700-residue receptor is nearly guaranteed; beating the rest of the "
        "channel domain is not.",
        "",
        "| paralog | sites | n | in | mean JSD | protein | own element | "
        "p vs protein | p vs own element | invariant |",
        "|---|---|---|---|---|---|---|---|---|---|",
        *rows,
        "",
        "**The measured IP₃ contacts are more constrained than the rest of "
        "the domains that carry them**, in all three paralogs — which is a "
        "statement about ten residues and is reported as such. The filter and "
        "gate sets are two residues each: their means are given because the "
        "element-level test in §5 is where those elements are actually "
        "powered, and a p-value on n = 2 is not offered.",
        "",
        "### 6.1 What the three copies have kept identical",
        "",
        "A second, independent instrument on the same question: sequence "
        "identity *between* the paralogs, measured over mutually covered "
        "columns only (S6's fragment-aware rule). It needs no alignment depth "
        "and no conservation metric at all.",
        "",
        f"- **the gate is {'identical' if ident == [1.0] else 'nearly identical'} "
        f"in all three pairs** ({', '.join(_f(x, 2) for x in ident)})",
        f"- the luminal loop is {', '.join(_f(x, 2) for x in loop_id)}",
        f"- the whole protein is {', '.join(_f(x, 2) for x in whole_id)} "
        f"(over all mutually covered columns — see §6.2)",
        "",
        "So the two extremes of paralog divergence in this receptor are "
        "**five residues that have not changed since 2R** and **fifty that "
        "retain 13–31 % identity**, and they sit inside the same domain. "
        "Whatever the three copies were free to differ in, it was not the "
        "gate.",
        "",
        "### 6.2 The number S6 published is a different measurement, and it checks out",
        "",
        "S6 reports between-paralog covered identity at 0.741–0.791 and this "
        f"section reports {', '.join(_f(x, 2) for x in whole_id)}. The two are "
        "not in conflict and neither is wrong: **S6 measured on "
        "`trimmed.fasta`**, which kept 1,797 of the alignment's 11,777 "
        "columns (and S6's headline is a group mean, where these are the "
        "human pairs). So each row of `paralog_identity_by_element.tsv` "
        "carries the same pair measured S6's way as well, and it reproduces "
        "S6's own committed matrix exactly: "
        f"{', '.join(_f(x, 3) for x in trimmed_id)} here against "
        f"{', '.join(_f(x, 3) for x in s6_id)} in "
        "`msa_v2/identity_covered.tsv`.",
        "",
        "S17 cannot use the trimmed alignment, and the retention column says "
        "why. trimAl's job is to delete divergent and gappy columns, and the "
        "element this section's headline rests on is the most divergent and "
        f"gappiest in the receptor: trimAl keeps "
        f"{loop_keep[0]} of the luminal loop's {loop_keep[1]} residues, "
        f"against {gate_keep[0]}/{gate_keep[1]} for the gate and "
        f"{filt_keep[0]}/{filt_keep[1]} for the filter. A per-element identity "
        "read off the trimmed alignment would be a measurement of what "
        "survived trimming.",
        "",
        P.line("within_paralog_identity", "orthogonal",
               f"measured S6's way S17 reproduces S6's number exactly; "
               f"measured over all mutually covered columns it is "
               f"{_f(min(whole_id), 3)}–{_f(max(whole_id), 3)}, about 0.09 "
               f"lower, and the whole difference is trimAl. These are the "
               f"same metric on two different column sets, so the comparison "
               f"is not a corroboration — what is new is the per-element "
               f"breakdown, a range from {_f(min(loop_id), 2)} to 1.00 inside "
               f"one protein, which the trimmed alignment cannot carry"),
        "",
        P.line("ryr_shares_the_domains", "confirmed",
               "the ranking puts RIH-associated and RIH among the most "
               "constrained domains, and both are shared with the ryanodine "
               "receptor — so the parts of an IP₃ receptor that cannot change "
               "are largely the parts that do not distinguish it from its "
               "sister family. The elements that *are* diagnostic of ligand "
               "gating are the ten contact residues, not a domain"),
    )


def build(h: dict) -> str:
    import s17_report_variants as V
    return A(s_elements(h), "", s_ligand(h), "", V.build(h))
