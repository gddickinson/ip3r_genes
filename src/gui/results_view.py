"""Results table + filter row.

A sortable `ttk.Treeview` showing one row per ProteinVariant. A text filter
at the top does substring matching across all visible columns; clicking a
column header sorts ascending / descending. The currently-selected variant is
exposed via `selected()` for the details pane.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from typing import Callable, Optional

from ..core.models import ProteinVariant


COLUMNS = [
    ("source",      "Source",      80),
    ("gene_symbol", "Gene",        100),
    ("accession",   "Accession",   170),
    ("species",     "Species",     180),
    ("length_aa",   "Length",      70),
    ("description", "Description", 480),
]


class ResultsView(ttk.Frame):
    def __init__(
        self,
        master: tk.Misc,
        on_select: Optional[Callable[[ProteinVariant], None]] = None,
    ) -> None:
        super().__init__(master, padding=(8, 0, 8, 8))
        self.on_select = on_select
        self._variants: dict[str, ProteinVariant] = {}  # iid -> variant
        self._all_iids: list[str] = []
        self._filter_var = tk.StringVar()
        self._sort_state: dict[str, bool] = {}  # col -> ascending?
        self._context_menu: Optional[tk.Menu] = None
        self._build()

    def _build(self) -> None:
        bar = ttk.Frame(self)
        bar.pack(fill="x", pady=(0, 4))
        ttk.Label(bar, text="Filter:").pack(side="left")
        e = ttk.Entry(bar, textvariable=self._filter_var)
        e.pack(side="left", fill="x", expand=True, padx=(4, 8))
        e.bind("<KeyRelease>", lambda _e: self._apply_filter())
        self.count_var = tk.StringVar(value="0 rows")
        ttk.Label(bar, textvariable=self.count_var).pack(side="right")

        cols = [c[0] for c in COLUMNS]
        self.tree = ttk.Treeview(self, columns=cols, show="headings", selectmode="browse")
        for col, label, width in COLUMNS:
            self.tree.heading(col, text=label, command=lambda c=col: self._sort(c))
            anchor = "e" if col == "length_aa" else "w"
            self.tree.column(col, width=width, anchor=anchor, stretch=(col == "description"))

        vsb = ttk.Scrollbar(self, orient="vertical", command=self.tree.yview)
        hsb = ttk.Scrollbar(self, orient="horizontal", command=self.tree.xview)
        self.tree.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        self.tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")
        hsb.pack(side="bottom", fill="x")

        self.tree.bind("<<TreeviewSelect>>", self._on_select)

    # ---- public API ------------------------------------------------------

    def add(self, variant: ProteinVariant) -> None:
        iid = variant.key()
        # de-dup across sources gracefully — same accession from same DB only once
        if iid in self._variants:
            return
        values = (
            variant.source,
            variant.gene_symbol,
            variant.accession,
            variant.species,
            variant.length_aa if variant.length_aa is not None else "",
            variant.description,
        )
        self.tree.insert("", "end", iid=iid, values=values)
        self._variants[iid] = variant
        self._all_iids.append(iid)
        self._update_count()

    def extend(self, variants: list[ProteinVariant]) -> None:
        for v in variants:
            self.add(v)

    def clear(self) -> None:
        for iid in self.tree.get_children():
            self.tree.delete(iid)
        self._variants.clear()
        self._all_iids.clear()
        self._update_count()

    def all_variants(self) -> list[ProteinVariant]:
        return list(self._variants.values())

    def selected(self) -> Optional[ProteinVariant]:
        sel = self.tree.selection()
        if not sel:
            return None
        return self._variants.get(sel[0])

    def add_context_action(self, label: str, callback: Callable[[ProteinVariant], None]) -> None:
        """Add a right-click menu action; `callback` receives the clicked row's variant."""
        if self._context_menu is None:
            self._context_menu = tk.Menu(self, tearoff=False)
            # macOS reports right-click as Button-2, elsewhere Button-3.
            self.tree.bind("<Button-2>", self._popup_context)
            self.tree.bind("<Button-3>", self._popup_context)
            self.tree.bind("<Control-Button-1>", self._popup_context)
        self._context_menu.add_command(
            label=label, command=lambda: self._run_context(callback))

    def _popup_context(self, event: tk.Event) -> None:
        iid = self.tree.identify_row(event.y)
        if not iid or self._context_menu is None:
            return
        self.tree.selection_set(iid)
        self._context_menu.tk_popup(event.x_root, event.y_root)

    def _run_context(self, callback: Callable[[ProteinVariant], None]) -> None:
        v = self.selected()
        if v:
            callback(v)

    # ---- internals -------------------------------------------------------

    def _on_select(self, _event: tk.Event) -> None:
        v = self.selected()
        if v and self.on_select:
            self.on_select(v)

    def _apply_filter(self) -> None:
        needle = self._filter_var.get().strip().lower()
        # detach everything then reattach the matching subset
        for iid in self._all_iids:
            self.tree.detach(iid)
        kept = 0
        for iid in self._all_iids:
            if not needle:
                self.tree.move(iid, "", "end")
                kept += 1
                continue
            v = self._variants[iid]
            hay = " ".join([
                v.source, v.gene_symbol, v.accession, v.species,
                str(v.length_aa or ""), v.description,
            ]).lower()
            if needle in hay:
                self.tree.move(iid, "", "end")
                kept += 1
        self._update_count(kept)

    def _update_count(self, visible: Optional[int] = None) -> None:
        total = len(self._all_iids)
        if visible is None:
            visible = total
        if visible == total:
            self.count_var.set(f"{total} rows")
        else:
            self.count_var.set(f"{visible} / {total} rows")

    def _sort(self, col: str) -> None:
        ascending = not self._sort_state.get(col, False)
        self._sort_state[col] = ascending

        def key(iid: str):
            v = self._variants[iid]
            value = getattr(v, col, "") if col != "length_aa" else (v.length_aa or 0)
            if col == "length_aa":
                return value
            return str(value).lower()

        ordered = sorted(self._all_iids, key=key, reverse=not ascending)
        self._all_iids = ordered
        self._apply_filter()  # re-attach in new order honoring filter
