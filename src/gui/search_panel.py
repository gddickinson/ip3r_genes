"""Search input panel.

Lets the user enter one or more gene symbols (comma-separated), pick a species
from a combobox (or leave blank for cross-species), toggle which databases to
query, and load one of the bundled presets.

Emits its result via a callback `on_search(query)` rather than holding a
reference to the orchestrator — keeps this widget testable in isolation.
"""

from __future__ import annotations

import json
import tkinter as tk
from pathlib import Path
from tkinter import ttk
from typing import Callable, Optional

from ..core.models import SearchQuery
from ..utils.species import all_species, resolve_species


class SearchPanel(ttk.Frame):
    def __init__(
        self,
        master: tk.Misc,
        on_search: Callable[[SearchQuery], None],
        on_cancel: Optional[Callable[[], None]] = None,
        presets_dir: Optional[Path] = None,
        on_special_preset: Optional[Callable[[str, dict], None]] = None,
    ) -> None:
        super().__init__(master, padding=8)
        self.on_search = on_search
        self.on_cancel = on_cancel
        self.on_special_preset = on_special_preset
        self.presets_dir = presets_dir or (Path(__file__).resolve().parents[2] / "presets")
        self.loaded_preset: dict = {}   # raw JSON of the last-loaded preset

        self.gene_var = tk.StringVar()
        self.species_var = tk.StringVar(value="(any)")
        self.max_var = tk.IntVar(value=50)
        self.include_seq_var = tk.BooleanVar(value=False)
        self.source_vars = {
            "NCBI": tk.BooleanVar(value=True),
            "Ensembl": tk.BooleanVar(value=True),
            "UniProt": tk.BooleanVar(value=True),
            "Compara": tk.BooleanVar(value=False),
            "AlphaFold": tk.BooleanVar(value=False),
            "Foldseek": tk.BooleanVar(value=False),
            "BLAST": tk.BooleanVar(value=False),
        }
        self.preset_var = tk.StringVar(value="")

        self._build()

    def _build(self) -> None:
        # Row 0 — gene + species
        ttk.Label(self, text="Gene symbol(s):").grid(row=0, column=0, sticky="w")
        gene_entry = ttk.Entry(self, textvariable=self.gene_var, width=32)
        gene_entry.grid(row=0, column=1, sticky="we", padx=(4, 12))
        gene_entry.bind("<Return>", lambda _e: self._fire())

        ttk.Label(self, text="Species:").grid(row=0, column=2, sticky="w")
        species_values = ["(any)"] + [f"{s.common} — {s.scientific}" for s in all_species()]
        self.species_combo = ttk.Combobox(
            self, textvariable=self.species_var, values=species_values, width=32
        )
        self.species_combo.grid(row=0, column=3, sticky="we", padx=(4, 0))

        # Row 1 — DB toggles + sequence toggle + max
        srcs = ttk.Frame(self)
        srcs.grid(row=1, column=0, columnspan=4, sticky="we", pady=(8, 0))
        ttk.Label(srcs, text="Sequence DBs:").pack(side="left")
        for name in ("NCBI", "Ensembl", "UniProt"):
            ttk.Checkbutton(srcs, text=name, variable=self.source_vars[name]).pack(side="left", padx=4)
        ttk.Checkbutton(srcs, text="Compara (paralogs)", variable=self.source_vars["Compara"]).pack(side="left", padx=4)
        ttk.Separator(srcs, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Label(srcs, text="Structure:").pack(side="left")
        ttk.Checkbutton(srcs, text="AlphaFold DB", variable=self.source_vars["AlphaFold"]).pack(side="left", padx=4)
        ttk.Separator(srcs, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Label(srcs, text="Bait:").pack(side="left")
        ttk.Checkbutton(srcs, text="BLAST*", variable=self.source_vars["BLAST"]).pack(side="left", padx=4)
        ttk.Checkbutton(srcs, text="Foldseek*", variable=self.source_vars["Foldseek"]).pack(side="left", padx=4)
        ttk.Separator(srcs, orient="vertical").pack(side="left", fill="y", padx=8)
        ttk.Checkbutton(srcs, text="Fetch sequences", variable=self.include_seq_var).pack(side="left", padx=4)
        ttk.Label(srcs, text="Max/source:").pack(side="left", padx=(12, 2))
        ttk.Spinbox(srcs, from_=1, to=500, textvariable=self.max_var, width=5).pack(side="left")

        hint = ttk.Label(
            self,
            text="Compara mines Ensembl gene trees for paralogs regardless of naming (unnamed genes = candidates). "
                 "* BLAST (sequence) and Foldseek (structure) search with the longest fetched sequence per gene as "
                 "bait — require 'Fetch sequences'. Async: BLAST 2–10 min, Foldseek 1–3 min.",
            foreground="#666",
        )
        hint.grid(row=3, column=0, columnspan=4, sticky="w", pady=(4, 0))

        # Row 2 — presets + buttons
        action = ttk.Frame(self)
        action.grid(row=2, column=0, columnspan=4, sticky="we", pady=(8, 0))
        ttk.Label(action, text="Preset:").pack(side="left")
        self.preset_combo = ttk.Combobox(action, textvariable=self.preset_var, values=self._list_presets(), width=24)
        self.preset_combo.pack(side="left", padx=4)
        ttk.Button(action, text="Load", command=self._load_preset).pack(side="left", padx=2)
        ttk.Button(action, text="Search", command=self._fire).pack(side="right")
        if self.on_cancel:
            ttk.Button(action, text="Cancel", command=self.on_cancel).pack(side="right", padx=4)

        self.columnconfigure(1, weight=1)
        self.columnconfigure(3, weight=1)

    # ---- presets ---------------------------------------------------------

    def _list_presets(self) -> list[str]:
        if not self.presets_dir.exists():
            return []
        return sorted(p.stem for p in self.presets_dir.glob("*.json"))

    def _load_preset(self) -> None:
        name = self.preset_var.get().strip()
        if not name:
            return
        self.load_preset(name)

    def load_preset(self, name: str) -> bool:
        """Programmatic preset load — used by the GUI button and `run.py --preset`.
        Returns True if the preset was found and applied.
        """
        path = self.presets_dir / f"{name}.json"
        if not path.exists():
            return False
        try:
            data = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            return False
        if data.get("mode"):
            # Not a form-fill preset — special modes ("domain_scan",
            # "exhaustive") are routed to the Analysis tools.
            if self.on_special_preset:
                self.preset_var.set(name)
                self.on_special_preset(name, data)
                return True
            return False
        self.loaded_preset = data
        self.preset_var.set(name)
        self.gene_var.set(", ".join(data.get("genes", [])))
        sp = data.get("species", "")
        self.species_var.set(sp if sp else "(any)")
        for src_name in self.source_vars:
            self.source_vars[src_name].set(src_name in data.get("sources", list(self.source_vars)))
        if "max_per_source" in data:
            try:
                self.max_var.set(int(data["max_per_source"]))
            except (TypeError, ValueError):
                pass
        if "include_sequence" in data:
            self.include_seq_var.set(bool(data["include_sequence"]))
        return True

    def fire(self) -> None:
        """Public alias for programmatic 'Search' click."""
        self._fire()

    # ---- search submission ----------------------------------------------

    def _parse_species(self) -> tuple[str, Optional[int]]:
        raw = self.species_var.get().strip()
        if not raw or raw == "(any)":
            return "", None
        # Combobox values come in "Common — Scientific" form; resolve to scientific
        if "—" in raw:
            scientific = raw.split("—", 1)[1].strip()
        else:
            scientific = raw
        info = resolve_species(scientific)
        if info:
            return info.scientific, info.taxon_id
        return scientific, None

    def _fire(self) -> None:
        symbols = [s.strip() for s in self.gene_var.get().split(",") if s.strip()]
        if not symbols:
            return
        species, taxon = self._parse_species()
        sources = [k for k, v in self.source_vars.items() if v.get()]
        if not sources:
            return
        query = SearchQuery(
            gene_symbols=symbols,
            species=species,
            taxon_id=taxon,
            sources=sources,
            max_results_per_source=max(1, int(self.max_var.get() or 50)),
            include_sequence=bool(self.include_seq_var.get()),
        )
        self.on_search(query)
