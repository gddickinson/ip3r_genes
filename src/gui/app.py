"""Main application window.

Lays out the three sub-widgets (search panel on top, results table in the
middle, details pane on the right), owns the SearchOrchestrator, polls the
worker-result queue with `after()`, and handles export menu commands.
"""

from __future__ import annotations

import queue
import tkinter as tk
import webbrowser
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from typing import Optional

from ..core.cache import DiskCache
from ..core.models import ProteinVariant, SearchQuery, SearchResult, SearchStatus
from ..core.search import SearchOrchestrator
from ..utils.exporters import write_csv, write_fasta, write_json
from ..utils.results_writer import make_bundle_dir, write_bundle
from .details_view import DetailsView
from .results_view import ResultsView
from .search_panel import SearchPanel
from .tools import ToolsController


POLL_MS = 100


class MainWindow:
    def __init__(self, root: tk.Tk, project_root: Path, email: str = "") -> None:
        self.root = root
        self.project_root = project_root
        self.email = email
        self.root.title("Protein Variant Finder")
        self.root.geometry("1280x780")

        cache_dir = project_root / "cache"
        self.cache = DiskCache(cache_dir)
        self.orchestrator = SearchOrchestrator(cache=self.cache, email=email)
        self._pending_sources: set[str] = set()
        self.last_query: Optional[SearchQuery] = None
        self.last_results: list[SearchResult] = []
        self.tools = ToolsController(self)
        self._build_menu()
        self._build_layout()
        self._build_status()
        self.root.after(POLL_MS, self._drain_queue)

    # ---- layout ----------------------------------------------------------

    def _build_menu(self) -> None:
        menubar = tk.Menu(self.root)
        filem = tk.Menu(menubar, tearoff=False)
        filem.add_command(label="Save results bundle (+ report)…", command=self._save_bundle)
        filem.add_separator()
        filem.add_command(label="Export FASTA…", command=lambda: self._export("fasta"))
        filem.add_command(label="Export CSV…",   command=lambda: self._export("csv"))
        filem.add_command(label="Export JSON…",  command=lambda: self._export("json"))
        filem.add_separator()
        filem.add_command(label="Clear cache", command=self._clear_cache)
        filem.add_command(label="Quit", command=self.root.destroy)
        menubar.add_cascade(label="File", menu=filem)

        self.tools.attach_menu(menubar)

        helpm = tk.Menu(menubar, tearoff=False)
        helpm.add_command(label="About", command=self._about)
        menubar.add_cascade(label="Help", menu=helpm)
        self.root.config(menu=menubar)

    def _build_layout(self) -> None:
        self.search = SearchPanel(
            self.root,
            on_search=self._on_search,
            on_cancel=self._on_cancel,
            presets_dir=self.project_root / "presets",
            on_special_preset=lambda name, data: self.tools.run_special_preset(name, data),
        )
        self.search.pack(fill="x")

        paned = ttk.PanedWindow(self.root, orient="horizontal")
        paned.pack(fill="both", expand=True)

        self.results = ResultsView(paned, on_select=self._on_select)
        self.details = DetailsView(paned)
        paned.add(self.results, weight=3)
        paned.add(self.details, weight=2)
        self.results.add_context_action(
            "Investigate this accession…", lambda _v: self.tools.investigate())
        self.results.add_context_action(
            "Selection test (dN/dS)…", lambda _v: self.tools.run_selection_test())
        self.results.add_context_action(
            "Fold check (ESMFold)…", lambda _v: self.tools.run_fold_check())

    def _build_status(self) -> None:
        self.status_var = tk.StringVar(value="Ready.")
        bar = ttk.Frame(self.root, padding=(8, 2))
        bar.pack(fill="x")
        ttk.Label(bar, textvariable=self.status_var).pack(side="left")

    # ---- events ----------------------------------------------------------

    def _on_search(self, query: SearchQuery) -> None:
        self.results.clear()
        self.details.show(None)
        self._pending_sources = set(query.sources)
        self.last_query = query
        self.last_results = []
        self.tools.reset()
        n = self.orchestrator.start(query, on_done=lambda: None)
        if n == 0:
            self.status("No databases selected.")
            return
        self.status(f"Querying {n} source(s): {', '.join(query.sources)} for {', '.join(query.gene_symbols)}…")

    def _on_cancel(self) -> None:
        self.orchestrator.cancel()
        self.status("Cancellation requested — workers finish current call.")

    def _on_select(self, variant) -> None:
        self.details.show(variant)

    def load_variants(
        self,
        variants: list[ProteinVariant],
        query: Optional[SearchQuery] = None,
        results: Optional[list[SearchResult]] = None,
    ) -> None:
        """Replace the table contents programmatically (e.g. domain-bait scan)."""
        self.results.clear()
        self.details.show(None)
        self._pending_sources = set()
        self.results.extend(variants)
        if query is not None:
            self.last_query = query
        self.last_results = list(results or [])

    def _drain_queue(self) -> None:
        try:
            while True:
                result: SearchResult = self.orchestrator.queue.get_nowait()
                self._apply_result(result)
        except queue.Empty:
            pass
        self.root.after(POLL_MS, self._drain_queue)

    def _apply_result(self, result: SearchResult) -> None:
        self._pending_sources.discard(result.source)
        self.last_results.append(result)
        if result.status == SearchStatus.OK:
            self.results.extend(result.variants)
            note = f" {result.message}" if result.message else ""
            msg = f"{result.source}: {len(result.variants)} hit(s) in {result.elapsed_s:.1f}s{note}"
        elif result.status == SearchStatus.EMPTY:
            msg = f"{result.source}: no hits"
        else:
            msg = f"{result.source}: ERROR — {result.message}"
        if self._pending_sources:
            tail = f"  (waiting on {', '.join(sorted(self._pending_sources))})"
        else:
            tail = f"  · {len(self.results.all_variants())} total rows"
        self.status(msg + tail)

    def status(self, text: str) -> None:
        self.status_var.set(text)

    # ---- export / housekeeping ------------------------------------------

    def _export(self, fmt: str) -> None:
        variants = self.results.all_variants()
        if not variants:
            messagebox.showinfo("Export", "No results to export.")
            return
        ext = {"fasta": ".fasta", "csv": ".csv", "json": ".json"}[fmt]
        path = filedialog.asksaveasfilename(
            defaultextension=ext,
            filetypes=[(fmt.upper(), f"*{ext}"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            if fmt == "fasta":
                n = write_fasta(variants, path)
                msg = f"Wrote {n} FASTA record(s).\nVariants without sequences were skipped — re-search with 'Fetch sequences' checked."
            elif fmt == "csv":
                n = write_csv(variants, path)
                msg = f"Wrote {n} rows to CSV."
            else:
                n = write_json(variants, path)
                msg = f"Wrote {n} rows to JSON."
            messagebox.showinfo("Export", msg)
        except OSError as e:
            messagebox.showerror("Export failed", str(e))

    def _save_bundle(self) -> None:
        variants = self.results.all_variants()
        if not variants or not self.last_query:
            messagebox.showinfo("Save results bundle", "No results to save — run a search first.")
            return
        results_root = self.project_root / "results"
        results_root.mkdir(exist_ok=True)
        label = "_".join(self.last_query.gene_symbols[:3]) or "results"
        bundle_dir = make_bundle_dir(results_root, label)
        notes = "Saved from Protein Variant Finder GUI."
        summary = write_bundle(bundle_dir, variants, self.last_query, self.last_results, notes=notes)
        report_paths = self.tools.augment_bundle(bundle_dir, variants, self.last_query, notes=notes)
        extras = []
        if self.tools.analysis:
            extras.append("  • analysis/ (MSA, tree, clusters, mutations)")
        if self.tools.discovery:
            extras.append("  • discovery/ (candidates.tsv + report)")
        msg = (
            f"Wrote {summary['n_variants']} variants to:\n{summary['dir']}\n\n"
            "Bundle contains:\n"
            "  • results.csv / .json / .fasta\n"
            "  • metadata.json (query + per-source stats)\n"
            "  • summary.md (human-readable report)\n"
            + ("\n".join(extras) + "\n" if extras else "")
            + "  • report.md / report.html (+ figures)"
        )
        self.status(f"Saved bundle to {bundle_dir}")
        if messagebox.askyesno("Saved", msg + "\n\nOpen report.html in your browser?"):
            webbrowser.open(Path(report_paths["report_html"]).resolve().as_uri())

    def _clear_cache(self) -> None:
        n = self.cache.clear()
        self.status(f"Cleared {n} cache file(s).")

    def _about(self) -> None:
        messagebox.showinfo(
            "Protein Variant Finder",
            "A GUI tool for enumerating protein isoforms / orthologs across\n"
            "NCBI, Ensembl, UniProt, AlphaFold DB, and Foldseek.\n\n"
            "The Analysis menu adds MSA/tree/cluster analysis, novel-paralog\n"
            "discovery, domain-bait scans, and deep-dive investigations.\n\n"
            "First use case: the IP3 receptor (ITPR) calcium-release family.",
        )


def run(
    project_root: Optional[Path] = None,
    email: str = "",
    preset: str = "",
    auto_search: bool = False,
) -> None:
    root = tk.Tk()
    try:
        ttk.Style().theme_use("clam")
    except tk.TclError:
        pass
    project_root = project_root or Path(__file__).resolve().parents[2]
    win = MainWindow(root, project_root=project_root, email=email)
    if preset:
        loaded = win.search.load_preset(preset)
        if loaded and auto_search:
            # Defer until the window is mapped so status updates render.
            root.after(200, win.search.fire)
    root.mainloop()
