"""Modal option dialogs for the Analysis menu.

Each dialog is a small modal Toplevel; construct it and read `.result`
afterwards — a dict of the chosen options, or None if the user cancelled.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import ttk
from typing import Optional


class _ModalDialog(tk.Toplevel):
    """OK/Cancel modal base. Subclasses build rows in `body()` and return
    their options dict from `collect()` (returning None keeps it open)."""

    def __init__(self, master: tk.Misc, title: str) -> None:
        super().__init__(master)
        self.title(title)
        self.resizable(False, False)
        self.result: Optional[dict] = None
        frame = ttk.Frame(self, padding=12)
        frame.pack(fill="both", expand=True)
        frame.columnconfigure(1, weight=1)
        self._next_row = 0
        self.body(frame)
        btns = ttk.Frame(frame)
        btns.grid(row=99, column=0, columnspan=2, sticky="e", pady=(12, 0))
        ttk.Button(btns, text="Cancel", command=self.destroy).pack(side="right", padx=(6, 0))
        ttk.Button(btns, text="OK", command=self._ok).pack(side="right")
        self.bind("<Return>", lambda _e: self._ok())
        self.bind("<Escape>", lambda _e: self.destroy())
        self.transient(master.winfo_toplevel())
        self.grab_set()
        self.wait_window(self)

    def body(self, frame: ttk.Frame) -> None:  # pragma: no cover - abstract
        raise NotImplementedError

    def collect(self) -> Optional[dict]:  # pragma: no cover - abstract
        raise NotImplementedError

    def _ok(self) -> None:
        result = self.collect()
        if result is None:
            return
        self.result = result
        self.destroy()

    # -- layout helpers ----------------------------------------------------

    def row(self, frame: ttk.Frame, label: str, widget: tk.Widget) -> None:
        ttk.Label(frame, text=label).grid(row=self._next_row, column=0, sticky="w", pady=2)
        widget.grid(row=self._next_row, column=1, sticky="we", padx=(8, 0), pady=2)
        self._next_row += 1

    def note(self, frame: ttk.Frame, text: str) -> None:
        lbl = ttk.Label(frame, text=text, foreground="#666", wraplength=420, justify="left")
        lbl.grid(row=self._next_row, column=0, columnspan=2, sticky="w", pady=(6, 0))
        self._next_row += 1


class AnalysisDialog(_ModalDialog):
    """Options for the MSA / tree / clusters / mutations pipeline."""

    def __init__(self, master: tk.Misc, n_with_seq: int) -> None:
        self.n_with_seq = n_with_seq
        self.identity_var = tk.DoubleVar(master, value=0.5)
        self.conservation_var = tk.DoubleVar(master, value=0.6)
        self.mafft_var = tk.BooleanVar(master, value=False)
        self.max_var = tk.IntVar(master, value=60)
        super().__init__(master, "Run sequence analysis")

    def body(self, frame: ttk.Frame) -> None:
        self.row(frame, "Cluster identity threshold:",
                 ttk.Spinbox(frame, from_=0.1, to=0.95, increment=0.05,
                             textvariable=self.identity_var, width=8))
        self.row(frame, "Conservation threshold:",
                 ttk.Spinbox(frame, from_=0.1, to=0.95, increment=0.05,
                             textvariable=self.conservation_var, width=8))
        self.row(frame, "Max variants in MSA:",
                 ttk.Spinbox(frame, from_=3, to=200, textvariable=self.max_var, width=8))
        self.row(frame, "", ttk.Checkbutton(frame, text="Use MAFFT if installed (better MSA)",
                                            variable=self.mafft_var))
        self.note(frame, f"{self.n_with_seq} variant(s) have sequences and will be analysed "
                         "(longest kept if over the max). Builds MSA → identity matrix → "
                         "NJ tree → clusters → mutation calls → conservation profile.")

    def collect(self) -> Optional[dict]:
        try:
            ident = float(self.identity_var.get())
            cons = float(self.conservation_var.get())
            max_n = int(self.max_var.get())
        except (tk.TclError, ValueError):
            return None
        if not (0.0 < ident < 1.0 and 0.0 < cons < 1.0 and max_n >= 3):
            return None
        return {
            "identity_threshold": ident,
            "conservation_threshold": cons,
            "use_mafft": bool(self.mafft_var.get()),
            "max_variants": max_n,
        }


class DiscoveryDialog(_ModalDialog):
    """Options for the novel-paralog discovery scorer."""

    def __init__(self, master: tk.Misc, default_known: list[str], has_analysis: bool) -> None:
        self.has_analysis = has_analysis
        self.known_var = tk.StringVar(master, value=", ".join(default_known))
        self.interpro_var = tk.BooleanVar(master, value=False)
        super().__init__(master, "Discover novel paralogs")

    def body(self, frame: ttk.Frame) -> None:
        self.row(frame, "Known paralog gene symbols:",
                 ttk.Entry(frame, textvariable=self.known_var, width=36))
        self.row(frame, "", ttk.Checkbutton(
            frame, text="Query InterPro for Pfam domain signatures (slower, network)",
            variable=self.interpro_var))
        if not self.has_analysis:
            self.note(frame, "Tip: run 'Run sequence analysis…' first — its identity matrix "
                             "powers the sequence-outlier and cluster-exclusion evidence "
                             "(worth up to 30 of the 100 points).")

    def collect(self) -> Optional[dict]:
        known = [s.strip() for s in self.known_var.get().split(",") if s.strip()]
        if not known:
            return None
        return {"known_paralogs": known, "interpro": bool(self.interpro_var.get())}


class InvestigateDialog(_ModalDialog):
    """Options for the single-accession deep-dive investigation."""

    def __init__(self, master: tk.Misc, default_accession: str = "") -> None:
        self.acc_var = tk.StringVar(master, value=default_accession)
        self.panel_var = tk.StringVar(master, value="Q92508, Q9H5I5")
        self.foldseek_var = tk.BooleanVar(master, value=False)
        super().__init__(master, "Investigate accession")

    def body(self, frame: ttk.Frame) -> None:
        self.row(frame, "UniProt accession:", ttk.Entry(frame, textvariable=self.acc_var, width=24))
        self.row(frame, "Comparison panel (accessions):",
                 ttk.Entry(frame, textvariable=self.panel_var, width=36))
        self.row(frame, "", ttk.Checkbutton(
            frame, text="Also run Foldseek structural search (adds 1–3 min)",
            variable=self.foldseek_var))
        self.note(frame, "Gathers seven lines of evidence (domain architecture, UniProt "
                         "cross-refs, AlphaFold structure, phylogenetic placement vs the "
                         "panel, synteny, literature) and renders a case file with a "
                         "verdict. Takes ~30–90 s. Default panel: human ITPR1/2/3.")

    def collect(self) -> Optional[dict]:
        acc = self.acc_var.get().strip()
        if not acc:
            return None
        return {
            "accession": acc,
            "panel_csv": self.panel_var.get().strip(),
            "foldseek": bool(self.foldseek_var.get()),
        }


class SelectionTestDialog(_ModalDialog):
    """Options for the dN/dS selection-pressure test."""

    def __init__(self, master: tk.Misc, default_accession: str = "") -> None:
        self.acc_var = tk.StringVar(master, value=default_accession)
        self.partner_var = tk.StringVar(master, value="")
        self.max_partners_var = tk.IntVar(master, value=5)
        super().__init__(master, "Selection test (dN/dS)")

    def body(self, frame: ttk.Frame) -> None:
        self.row(frame, "Ensembl protein/gene ID:",
                 ttk.Entry(frame, textvariable=self.acc_var, width=28))
        self.row(frame, "Partner ID (blank = auto):",
                 ttk.Entry(frame, textvariable=self.partner_var, width=28))
        self.row(frame, "Max auto-partners to try:",
                 ttk.Spinbox(frame, from_=1, to=10, textvariable=self.max_partners_var, width=6))
        self.note(frame, "Computes pairwise dN/dS (Nei–Gojobori) from Ensembl CDS. "
                         "dN/dS ≪ 1 = purifying selection = functional gene; ≈ 1 = "
                         "pseudogene-like drift. Auto mode tries the candidate's own "
                         "closest Compara orthologs and keeps the least-saturated "
                         "comparison. Works for Ensembl IDs only (ENS…P/T/G). ~30–90 s.")

    def collect(self) -> Optional[dict]:
        acc = self.acc_var.get().strip()
        if not acc:
            return None
        try:
            k = max(1, int(self.max_partners_var.get()))
        except (tk.TclError, ValueError):
            k = 5
        return {"accession": acc, "partner": self.partner_var.get().strip(),
                "max_partners": k}


class FoldCheckDialog(_ModalDialog):
    """Options for the ESMFold segment-fold check."""

    def __init__(self, master: tk.Misc, seq_len: int) -> None:
        self.seq_len = seq_len
        self.where_var = tk.StringVar(master, value="cterm")
        self.len_var = tk.IntVar(master, value=min(350, seq_len))
        super().__init__(master, "Fold check (ESMFold)")

    def body(self, frame: ttk.Frame) -> None:
        seg = ttk.Frame(frame)
        ttk.Radiobutton(seg, text="C-terminal (pore + tail — most conserved)",
                        variable=self.where_var, value="cterm").pack(anchor="w")
        ttk.Radiobutton(seg, text="N-terminal",
                        variable=self.where_var, value="nterm").pack(anchor="w")
        self.row(frame, "Segment:", seg)
        self.row(frame, "Segment length (aa, ≤400):",
                 ttk.Spinbox(frame, from_=50, to=400, textvariable=self.len_var, width=6))
        self.note(frame, f"Sequence is {self.seq_len} aa. ESMFold's public API folds "
                         "≤400 aa, no MSA needed — it works for proteins AlphaFold DB "
                         "doesn't cover (unnamed gene models, fresh BLAST hits). "
                         "Reports mean pLDDT; >70 = confident fold. ~30–120 s.")

    def collect(self) -> Optional[dict]:
        try:
            length = max(50, min(400, int(self.len_var.get())))
        except (tk.TclError, ValueError):
            return None
        return {"where": self.where_var.get(), "length": length}


class DomainScanDialog(_ModalDialog):
    """Options for the domain-bait scan (novel-member hunt)."""

    def __init__(self, master: tk.Misc, presets: list[str], default: str = "") -> None:
        self.presets = presets
        self.preset_var = tk.StringVar(master, value=default or (presets[0] if presets else ""))
        self.save_var = tk.BooleanVar(master, value=True)
        super().__init__(master, "Domain-bait scan")

    def body(self, frame: ttk.Frame) -> None:
        combo = ttk.Combobox(frame, textvariable=self.preset_var,
                             values=self.presets, state="readonly", width=30)
        self.row(frame, "Scan preset:", combo)
        self.row(frame, "", ttk.Checkbutton(
            frame, text="Save results bundle (analysis + discovery + report)",
            variable=self.save_var))
        self.note(frame, "Enumerates every protein carrying the family's Pfam signature "
                         "via InterPro, filters out known paralogs, fetches sequences for "
                         "the top candidates, then runs analysis + discovery scoring. The "
                         "hits replace the current results table. Takes several minutes.")

    def collect(self) -> Optional[dict]:
        preset = self.preset_var.get().strip()
        if not preset:
            return None
        return {"preset": preset, "save_results": bool(self.save_var.get())}
