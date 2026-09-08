"""s15b_sensitivity.py — the deliverable: which settings manufacture a loss.

S15a's matrix holds no `absent` cell, so there is no loss to place and the
count is zero.  What is worth reporting is therefore the shape of that
zero: how far the rules have to be moved before a loss appears, which rule
has to move, and how many losses each move buys.

The grid is the brief's four axes.  Three of them change the character:
the **evidence** ladder (`s15b_coding.EVIDENCE_LEVELS`), D45's
`paralog_unassignable` rule and D4's contiguity rule.  The fourth,
**branch lengths**, cannot change a parsimony count at all — Dollo counts
edges — so it is carried as a column and its invariance is reported as a
result rather than left as an omission.  It moves the Mk fits, and only
those.

The **coding axis** is D46's: `family` asks whether the assembly carries
the family at all and is the primary, `paralog` asks per cell and is the
sensitivity axis, because the `co_trace` population measured how
unreliable per-fragment paralog attribution is in a shattered assembly.

`manufactured_losses()` reports the same result at cell resolution: every
genome x cell that ever reaches `absent`, the state and rule that held it
at S15a's operating point, and the loosest setting that releases it.  A
row there is not a candidate loss; it is a measurement of what one
particular decision is worth in losses.
"""

from __future__ import annotations

import s15b_coding as coding
import s15b_dollo as dollo
import s15b_lib as lib
import s15b_mk as mk

CODINGS = ("family", "paralog")

#: where Dollo's single gain goes, per coding, and why.
#:
#: The **family** character is pinned at the root: S20 swept 6,928
#: non-vertebrate eukaryotic reference proteomes and S23 194 non-vertebrate
#: genomes, and both found IP3 receptors well outside the vertebrates, so
#: the family was present at the root of any vertebrate tree and a
#: vertebrate clade carrying none has lost it.  The **paralog** characters
#: are not pinned: S13 places the ITPR2/ITPR3 duplication on the
#: gnathostome stem, so a cyclostome that has neither never had either, and
#: pinning would score the origin of the paralogs as two losses.
PIN_GAIN = {"family": True, "paralog": False}


def _characters(recoded: list[dict], coding_name: str):
    """`[(character_name, {tip: state})]` under one coding."""
    if coding_name == "family":
        return [("ITPR_family", coding.family_states(recoded))]
    return [(c, s) for c, s in sorted(coding.paralog_states(recoded).items())]


def _state_counts(states: dict[str, str]) -> dict:
    present = sum(1 for s in states.values()
                  if s in coding.PRESENT_STATES or s == "present")
    absent = sum(1 for s in states.values() if s == "absent")
    return dict(n_tips_scored=len(states), n_present=present, n_absent=absent,
                n_undecided=len(states) - present - absent)


def dollo_counts(rows: list[dict], root: lib.Node) -> list[dict]:
    """Dollo over every (coding, character, setting).  No lengths read."""
    out = []
    for setting in coding.settings_grid():
        recoded = coding.recode_all(rows, **setting)
        for coding_name in CODINGS:
            for name, states in _characters(recoded, coding_name):
                res = dollo.dollo(root, coding.to_character(states),
                                  pin_gain_at_root=PIN_GAIN[coding_name])
                row = dict(coding=coding_name, character=name,
                           evidence=setting["evidence"],
                           use_r5=int(setting["use_r5"]),
                           use_r6=int(setting["use_r6"]),
                           is_base=int(setting == coding.BASE_SETTING))
                row.update(_state_counts(states))
                row.update(dollo_losses_max=res["n_max"],
                           dollo_losses_min=res["n_min"],
                           gain_node=res["gain_node"],
                           gain_rule=res["gain_rule"],
                           mrca_of_present=res["mrca_of_present"],
                           n_unknown_tips=res["n_unknown"],
                           resolution_note=res["note"])
                out.append(row)
    return out


def loss_edges(rows: list[dict], root: lib.Node) -> list[dict]:
    """Every loss edge Dollo places, under every setting that places one."""
    out = []
    for setting in coding.settings_grid():
        recoded = coding.recode_all(rows, **setting)
        for coding_name in CODINGS:
            for name, states in _characters(recoded, coding_name):
                res = dollo.dollo(root, coding.to_character(states),
                                  pin_gain_at_root=PIN_GAIN[coding_name])
                for ls in res["losses"]:
                    out.append(dict(
                        coding=coding_name, character=name,
                        evidence=setting["evidence"],
                        use_r5=int(setting["use_r5"]),
                        use_r6=int(setting["use_r6"]),
                        node=ls["node"], parent=ls["parent"],
                        parent_degree=ls["parent_degree"],
                        parent_is_polytomy=ls["parent_is_polytomy"],
                        siblings_lost_under_parent=ls[
                            "siblings_lost_under_parent"],
                        n_tips_lost=ls["n_tips_lost"],
                        n_tips_unknown=ls["n_tips_unknown"],
                        tips=ls["tips"]))
    return out


def manufactured_losses(rows: list[dict]) -> list[dict]:
    """Per genome x cell: does any setting make it a loss, and which.

    The columns that carry the argument are `base_state` / `base_rule` —
    what held the cell at S15a's operating point — and `needs_*`, the
    knobs that have to move before it reads as absent.
    """
    settings = list(coding.settings_grid())
    base = {(r["accession"], r["cell"]): d
            for r, d in zip(rows, coding.recode_all(
                rows, **coding.BASE_SETTING))}
    hits: dict[tuple[str, str], list[dict]] = {}
    for setting in settings:
        for r, d in zip(rows, coding.recode_all(rows, **setting)):
            if d["state"] == "absent":
                hits.setdefault((r["accession"], r["cell"]), []).append(setting)
    order = {name: i for i, name in enumerate(coding.EVIDENCE_NAMES)}
    out = []
    for r in rows:
        key = (r["accession"], r["cell"])
        got = hits.get(key)
        if not got:
            continue
        loosest = min(got, key=lambda s: (order[s["evidence"]],
                                          -int(s["use_r5"]), -int(s["use_r6"])))
        out.append(dict(
            accession=r["accession"], organism=r["organism"],
            vclass=r["vclass"], cell=r["cell"],
            ledger_status=r["ledger_status"],
            base_state=base[key]["state"], base_rule=base[key]["rule"],
            base_reason=r["reason"],
            recon_coverage=r["recon_coverage"],
            recon_n_contigs=r["recon_n_contigs"],
            n_spare_itpr_loci=r["n_spare_itpr_loci"],
            contig_spans_gene=r["contig_spans_gene"],
            n_settings_absent=len(got), n_settings=len(settings),
            loosest_evidence=loosest["evidence"],
            needs_r5_off=int(not loosest["use_r5"]),
            needs_r6_off=int(not loosest["use_r6"]),
            released_by=_released_by(loosest)))
    out.sort(key=lambda d: (-d["n_settings_absent"], d["accession"],
                            d["cell"]))
    return out


def _released_by(setting: dict) -> str:
    bits = []
    if setting["evidence"] != coding.BASE_SETTING["evidence"]:
        bits.append(f"evidence={setting['evidence']}")
    if not setting["use_r5"]:
        bits.append("R5 off (D45)")
    if not setting["use_r6"]:
        bits.append("R6 off (D4)")
    return "; ".join(bits) or "S15a's own setting"


def matrix(rows: list[dict], root: lib.Node,
           cal: dict[str, float] | None = None) -> tuple[list[dict], list[dict]]:
    """The four-axis matrix and the Mk fits behind its rate columns.

    The Dollo columns repeat across the branch-length axis by
    construction; that repetition is the point, and the report says so.
    Mk results are memoised on the character vector, so an evidence level
    that leaves the character unchanged costs nothing.
    """
    cal = cal if cal is not None else lib.load_calibrations()
    cache: dict[tuple, dict] = {}
    grid, fits = [], []
    for scheme in lib.BL_SCHEMES:
        lib.apply_lengths(root, lib.branch_lengths(root, scheme, cal))
        for setting in coding.settings_grid():
            recoded = coding.recode_all(rows, **setting)
            for coding_name in CODINGS:
                for name, states in _characters(recoded, coding_name):
                    char = coding.to_character(states)
                    res = dollo.dollo(root, char,
                                      pin_gain_at_root=PIN_GAIN[coding_name])
                    sig = (scheme, tuple(sorted(
                        (k, -1 if v is None else v) for k, v in char.items())))
                    row = dict(coding=coding_name, character=name,
                               evidence=setting["evidence"],
                               use_r5=int(setting["use_r5"]),
                               use_r6=int(setting["use_r6"]),
                               branch_lengths=scheme,
                               is_base=int(setting == coding.BASE_SETTING
                                           and scheme == "unit"))
                    row.update(_state_counts(states))
                    row.update(dollo_losses_max=res["n_max"],
                               dollo_losses_min=res["n_min"],
                               gain_node=res["gain_node"],
                               gain_rule=res["gain_rule"])
                    if sig not in cache:
                        got = {}
                        for model in mk.MODELS:
                            f = mk.fit(root, char, model)
                            got[model] = f
                            fits.append(dict(
                                coding=coding_name, character=name,
                                evidence=setting["evidence"],
                                first_seen_note="one row per distinct "
                                                "character; the fits are "
                                                "memoised on the character "
                                                "vector",
                                use_r5=int(setting["use_r5"]),
                                use_r6=int(setting["use_r6"]),
                                branch_lengths=scheme,
                                n_absent=row["n_absent"], **f))
                        cache[sig] = got
                    got = cache[sig]
                    best = _best_model(got)
                    row.update(
                        mk_fitted=int(got["ER"]["fitted"]),
                        mk_refusal=got["ER"]["reason"],
                        mk_best_model=best,
                        mk_loss_rate=(got[best]["rate"] if best else ""),
                        mk_irrev_loglik=got["irreversible"]["loglik"],
                        mk_irrev_rate=got["irreversible"]["rate"])
                    grid.append(row)
    return grid, fits


def _best_model(got: dict[str, dict]) -> str:
    usable = [(v["aic"], k) for k, v in got.items()
              if v["fitted"] and isinstance(v["aic"], float)]
    return min(usable)[1] if usable else ""


def profile_rows(root: lib.Node, rows: list[dict],
                 cal: dict[str, float] | None = None) -> list[dict]:
    """The likelihood profile behind the Mk section's refusal.

    Computed on the *primary* character at S15a's own setting, under every
    branch-length scheme, so `monotone decreasing` is a number a reader can
    check rather than a claim.
    """
    cal = cal if cal is not None else lib.load_calibrations()
    base = coding.recode_all(rows, **coding.BASE_SETTING)
    char = coding.to_character(coding.family_states(base))
    out = []
    for scheme in lib.BL_SCHEMES:
        lib.apply_lengths(root, lib.branch_lengths(root, scheme, cal))
        for model in mk.MODELS:
            axes = ("loss", "gain") if model == "ARD" else ("loss",)
            for axis in axes:
                prof = mk.profile(root, char, model, axis=axis)
                for rate, ll in zip(prof["grid"], prof["loglik"]):
                    out.append(dict(branch_lengths=scheme, model=model,
                                    axis=axis,
                                    rate=f"{rate:.6g}", loglik=round(ll, 6),
                                    shape=prof["shape"],
                                    argmax=f"{prof['argmax']:.6g}",
                                    at_boundary=prof["at_boundary"]))
    return out
