"""Detail pane for the currently-selected protein variant.

Shows a labelled summary, the URL (clickable: copied to clipboard on click),
and a scrollable sequence preview if a sequence has been fetched.
"""

from __future__ import annotations

import tkinter as tk
import webbrowser
from tkinter import ttk
from typing import Optional

from ..core.models import ProteinVariant


class DetailsView(ttk.Frame):
    def __init__(self, master: tk.Misc) -> None:
        super().__init__(master, padding=8)
        self._variant: Optional[ProteinVariant] = None
        self._build()

    def _build(self) -> None:
        self.fields: dict[str, tk.StringVar] = {}
        labels = [
            ("Source",      "source"),
            ("Gene",        "gene_symbol"),
            ("Accession",   "accession"),
            ("Species",     "species"),
            ("Taxon ID",    "taxon_id"),
            ("Length (aa)", "length_aa"),
            ("Transcript",  "transcript_id"),
            ("Gene ID",     "gene_id"),
            ("Description", "description"),
        ]
        grid = ttk.Frame(self)
        grid.pack(fill="x")
        for i, (label, attr) in enumerate(labels):
            ttk.Label(grid, text=label + ":", anchor="e", width=14).grid(row=i, column=0, sticky="ne", padx=(0, 6), pady=1)
            var = tk.StringVar()
            self.fields[attr] = var
            lbl = ttk.Label(grid, textvariable=var, anchor="w", wraplength=520, justify="left")
            lbl.grid(row=i, column=1, sticky="we", pady=1)
        grid.columnconfigure(1, weight=1)

        link_row = ttk.Frame(self)
        link_row.pack(fill="x", pady=(6, 4))
        ttk.Label(link_row, text="URL:", width=14, anchor="e").pack(side="left", padx=(0, 6))
        self.url_label = ttk.Label(link_row, text="", foreground="#1a4dbf", cursor="hand2")
        self.url_label.pack(side="left", fill="x", expand=True)
        self.url_label.bind("<Button-1>", self._open_url)

        ttk.Label(self, text="Sequence preview:").pack(anchor="w", pady=(8, 2))
        self.seq_text = tk.Text(self, height=10, wrap="none", font=("Menlo", 10))
        self.seq_text.pack(fill="both", expand=True)
        self.seq_text.configure(state="disabled")

    def show(self, variant: Optional[ProteinVariant]) -> None:
        self._variant = variant
        if variant is None:
            for v in self.fields.values():
                v.set("")
            self.url_label.configure(text="")
            self._set_seq("")
            return
        self.fields["source"].set(variant.source)
        self.fields["gene_symbol"].set(variant.gene_symbol)
        self.fields["accession"].set(variant.accession)
        self.fields["species"].set(variant.species)
        self.fields["taxon_id"].set(str(variant.taxon_id) if variant.taxon_id else "")
        self.fields["length_aa"].set(str(variant.length_aa) if variant.length_aa else "")
        self.fields["transcript_id"].set(variant.transcript_id)
        self.fields["gene_id"].set(variant.gene_id)
        self.fields["description"].set(variant.description)
        self.url_label.configure(text=variant.url)
        if variant.sequence:
            preview = self._format_sequence(variant.sequence)
            self._set_seq(preview)
        else:
            self._set_seq("(sequence not fetched — enable 'Fetch sequences' and re-search)")

    def _format_sequence(self, seq: str, width: int = 60) -> str:
        chunks = [seq[i : i + width] for i in range(0, len(seq), width)]
        return "\n".join(f"{i*width+1:>6}  {c}" for i, c in enumerate(chunks))

    def _set_seq(self, text: str) -> None:
        self.seq_text.configure(state="normal")
        self.seq_text.delete("1.0", "end")
        self.seq_text.insert("1.0", text)
        self.seq_text.configure(state="disabled")

    def _open_url(self, _event: tk.Event) -> None:
        if self._variant and self._variant.url:
            webbrowser.open(self._variant.url)
