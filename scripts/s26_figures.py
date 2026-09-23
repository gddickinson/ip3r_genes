"""S26 stage `figures`: number, check and copy one paper's figures.

Figures are numbered **in order of first mention in the paper's body**,
separately for main, Extended Data and Supplementary figures. The source text
cites a figure by its slug (`{fig:tree}`, `{fig:tree}b`) and never by a
number, so a paper can be reordered without renumbering anything by hand,
and the order check S14c had to add to the manuscript is true by
construction here.

Each figure is copied with `s14_figures._copy_one`, the manuscript's own
copier, so a paper figure is byte-identical to the committed one and is never
re-plotted (D13, D19). The stage fails on:

* a missing source file;
* a figure whose source is not primary in this paper (**P5**: to show a
  sibling's result, cite the sibling);
* fewer than four or more than seven main figures (**P4**);
* a declared figure that no body sentence mentions (R3), a `{fig:}`
  reference to a slug the paper does not declare, and a duplicated slug;
* a declared figure with no legend, a legend for an undeclared figure, or a
  legend in the wrong section for its kind;
* a multi-file figure whose legend does not name each of its panels.
"""

from __future__ import annotations

import re
import sys

import s14_figures
import s26_assign
import s26_lib as L

PREFIX = {"main": "Fig", "ed": "ExtDataFig", "supp": "SuppFig"}
LEGEND_SECTION = {"main": "05_figure_legends.md", "ed": "06_extended_data.md",
                  "supp": "06_extended_data.md"}


def numbering(pid: str) -> dict[str, tuple[str, int]]:
    """slug -> (kind, number), by first mention in the body sections."""
    kinds = {slug: kind for slug, kind, _s, _c in L.figures(pid)}
    order: list[str] = []
    for slug in L.FIG_REF.findall(L.body_text(pid)):
        if slug in kinds and slug not in order:
            order.append(slug)
    order += [s for s in kinds if s not in order]      # reported, not hidden
    out, count = {}, {k: 0 for k in L.KINDS}
    for slug in order:
        count[kinds[slug]] += 1
        out[slug] = (kinds[slug], count[kinds[slug]])
    return out


def legends(pid: str) -> dict[str, tuple[str, str]]:
    """slug -> (section file, legend paragraph)."""
    out = {}
    for name in sorted(L.LEGEND_SECTIONS):
        for para in L.section_text(pid, name).split("\n\n"):
            m = L.LEGEND_HEAD.match(para.strip())
            if m:
                out[m.group(1)] = (name, para.strip())
    return out


def check(pid: str) -> list[str]:
    fails = []
    figs = L.figures(pid)
    slugs = [f[0] for f in figs]
    for s in {s for s in slugs if slugs.count(s) > 1}:
        fails.append(f"{pid}: figure slug {s!r} declared twice")
    n_main = sum(1 for f in figs if f[1] == "main")
    if not L.MAIN_MIN <= n_main <= L.MAIN_MAX:
        fails.append(f"{pid}: {n_main} main figures; P4 allows "
                     f"{L.MAIN_MIN} to {L.MAIN_MAX}")
    body_refs = set(L.FIG_REF.findall(L.body_text(pid)))
    all_text = "\n".join(L.section_text(pid, n) for n in L.SECTIONS)
    for ref in set(L.FIG_REF.findall(all_text)) - set(slugs):
        fails.append(f"{pid}: {{fig:{ref}}} refers to no declared figure")
    legs = legends(pid)
    for slug, kind, stems, _cap in figs:
        if kind not in L.KINDS:
            fails.append(f"{pid}: figure {slug} has unknown kind {kind!r}")
            continue
        if slug not in body_refs:
            fails.append(f"{pid}: figure {slug} is placed but no body "
                         f"sentence refers to it (R3)")
        for stem in stems:
            owner = s26_assign.primary_paper(f"results/{stem}.png")
            if owner != pid:
                fails.append(f"{pid}: figure {slug} draws results/{stem}, "
                             f"which is primary in {owner!r} (P5)")
        if slug not in legs:
            fails.append(f"{pid}: figure {slug} has no legend")
            continue
        section, para = legs[slug]
        if section != LEGEND_SECTION[kind]:
            fails.append(f"{pid}: legend of {kind} figure {slug} is in "
                         f"{section}, expected {LEGEND_SECTION[kind]}")
        if len(stems) > 1:
            for i in range(len(stems)):
                letter = chr(97 + i)
                if not re.search(rf"\({letter}\)|\*\*{letter}\*\*|"
                                 rf"\b{letter},|\b{letter}\)", para):
                    fails.append(f"{pid}: legend of {slug} does not name "
                                 f"panel ({letter}) of its {len(stems)} "
                                 f"files")
                    break
    for slug in set(legs) - set(slugs):
        fails.append(f"{pid}: legend for undeclared figure {slug}")
    return fails


def run(pid: str) -> int:
    fails = check(pid)
    nums = numbering(pid)
    out_dir = L.paper_dir(pid) / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows, missing = [], 0
    for slug, kind, stems, cap in L.figures(pid):
        _k, n = nums[slug]
        for i, stem in enumerate(stems):
            suffix = "" if len(stems) == 1 else chr(97 + i)
            got = s14_figures._copy_one(
                stem, f"{PREFIX[kind]}{n}{suffix}_{slug}", dest_dir=out_dir)
            if not got:
                missing += 1
            for r in got:
                r.update(kind=kind, number=n, slug=slug, caption_stub=cap)
                r["deposited_as"] = r["deposited_as"].split(
                    f"papers/{pid}/", 1)[-1]
            rows.extend(got)
    keep = {r["deposited_as"].split("/")[-1] for r in rows}
    for stale in sorted(out_dir.iterdir()):
        if stale.is_file() and stale.name not in keep:
            stale.unlink()
    L.write_tsv(L.paper_dir(pid) / "figure_manifest.tsv", rows,
                ["kind", "number", "slug", "figure", "format", "width_in",
                 "height_in", "source", "deposited_as", "bytes", "sha256",
                 "caption_stub"])
    for f in fails:
        print(f"  [FAIL] {f}", file=sys.stderr)
    counts = {k: len({r['slug'] for r in rows if r['kind'] == k})
              for k in L.KINDS}
    print(f"[s26 figures {pid}] {counts['main']} main, {counts['ed']} "
          f"Extended Data, {counts['supp']} Supplementary ({len(rows)} "
          f"files, {missing} missing): {len(fails)} failure(s)")
    return 1 if fails or missing else 0
