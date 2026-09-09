"""S25 — the manuscript claims this thesis also makes, remapped to chapters.

The thesis states most of the paper's load-bearing numbers and many more. A
number quoted in both documents must be checked once, so those claims are
**carried** from `s14_claims` rather than restated: same identifier, same
source table, same operation, same expected value — only the chapter is new.
Copying them into a second ledger would create exactly the drift the ledger
exists to prevent.

`CARRIED` names, per chapter, the manuscript claim identifiers this thesis
also makes. A claim is carried only if the chapter's text actually states its
value, which `s25_claims._in_chapter` re-checks on every build, so this table
cannot quietly acquire a claim the thesis does not make.

The manuscript's ledger has 276 rows. The ones not carried here are of two
kinds: rows whose check is a phrase in a report rather than a value (those
cannot be held to the appearance rule, which is about values), and numbers the
thesis simply does not quote. Which rows those are is computable, because both
ledgers are committed and this table is a filter over one of them.

A carried identifier that no longer exists in the manuscript's ledger raises
at import, so the two documents cannot drift apart by deletion either.
"""

from __future__ import annotations

import s14_claims

#: chapter -> manuscript claim ids the thesis also makes.
CARRIED: dict[int, list[str]] = {
    3: [
        "C04", "C05", "C246b"
    ],
    4: [
        "C01", "C25", "C26", "C27", "C28", "C30", "C31", "C32",
        "C33", "C34", "C35", "C37", "C42", "C43", "C44", "C45",
        "C46", "C47", "C48"
    ],
    5: [
        "C02", "C03", "C07", "C08", "C09", "C11", "C12", "C13",
        "C14", "C15", "C16", "C17", "C19", "C20", "C21", "C24"
    ],
    6: [
        "C52", "C53", "C57", "C60"
    ],
    7: [
        "C51", "C63", "C64", "C66", "C67", "C68", "C69", "C70",
        "C76", "C77", "C78", "C79", "C80", "C81", "C82", "C83",
        "C84"
    ],
    8: [
        "C181", "C182", "C183", "C187", "C190", "C193", "C194",
        "C195", "C202", "C205", "C206", "C207", "C208", "C209",
        "C214", "C215", "C216", "C217", "C218", "C219", "C220",
        "C221", "C260", "C261", "C263", "C265", "C266", "C267",
        "C268"
    ],
    9: [
        "C85", "C86", "C88", "C89", "C90", "C91", "C92", "C93",
        "C94", "C95", "C96", "C97", "C98", "C99", "C110"
    ],
    11: [
        "C116", "C117", "C122", "C127", "C128", "C133", "C134",
        "C136", "C138", "C143", "C180", "C222", "C223", "C229",
        "C230", "C239", "C247", "C248"
    ],
    12: [
        "C251", "C252", "C253", "C254", "C255", "C246c", "C248c"
    ],
    13: [
        "C148", "C149", "C150", "C157", "C158", "C159", "C160",
        "C161", "C162", "C163", "C164", "C165", "C166", "C168",
        "C169", "C173", "C175", "C262", "C264"
    ],
    14: [
        "C176", "C177", "C178", "C179"
    ],
}

_BY_ID = {c["id"]: c for c in s14_claims.CLAIMS}

_missing = [i for ids in CARRIED.values() for i in ids if i not in _BY_ID]
if _missing:                                              # pragma: no cover
    raise ImportError("carried claim ids not in the manuscript ledger: "
                      + ", ".join(sorted(_missing)))

CLAIMS: list[dict] = [
    dict(_BY_ID[i], chapter=ch, section=f"ch{ch}")
    for ch in sorted(CARRIED) for i in CARRIED[ch]
]
