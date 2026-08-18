"""Analysis-menu controller.

Gives the GUI access to everything the headless CLI can do: the sequence
analysis pipeline, the novel-paralog discovery scorer (with optional
InterPro Pfam lookup), the domain-bait scan, single-accession deep-dive
investigations, and analysis/discovery/report artefacts inside saved
results bundles.

All heavy work runs through `TaskRunner` on a worker thread; the GUI never
blocks. Completed results are cached on the controller so "Save results
bundle" persists them alongside the raw hits, exactly like the CLI's
`--analyze --discover --save-results` path.
"""

from __future__ import annotations

import json
import tkinter as tk
import webbrowser
from pathlib import Path
from tkinter import filedialog, messagebox
from typing import Optional

from ..analysis import analyse
from ..analysis.esmfold import fold_segment
from ..analysis.pipeline import AnalysisResult, _label_for, write_analysis
from ..analysis.presence import build_presence_matrix
from ..analysis.selection import selection_test
from ..databases.interpro import batch_fetch_family_signatures
from ..discovery import (
    DiscoveryConfig,
    DiscoveryReport,
    build_signature_set,
    discover_novel_paralogs,
    run_domain_scan,
    run_exhaustive_hunt,
    write_discovery,
)
from ..investigation import InvestigationOptions, Investigator, write_investigation
from ..investigation.synthesis import assess_signals, render_case_file, verdict_from_signals
from ..utils.report import write_report
from ..utils.results_writer import make_bundle_dir
from .dialogs import (
    AnalysisDialog,
    DiscoveryDialog,
    DomainScanDialog,
    FoldCheckDialog,
    InvestigateDialog,
    SelectionTestDialog,
)
from .task_runner import TaskRunner
from .text_window import TextWindow


class ToolsController:
    """Owns the Analysis menu, its background tasks, and cached results."""

    def __init__(self, app) -> None:  # app: MainWindow (avoids circular import)
        self.app = app
        self.runner = TaskRunner(app.root)
        self.analysis: Optional[AnalysisResult] = None
        self.discovery: Optional[DiscoveryReport] = None
        self.presence = None                      # last PresenceReport
        self.novel_labels: set[str] = set()       # analysis labels to highlight

    def attach_menu(self, menubar: tk.Menu) -> None:
        m = tk.Menu(menubar, tearoff=False)
        m.add_command(label="Run sequence analysis…", command=self.run_analysis)
        m.add_command(label="Discover novel paralogs…", command=self.run_discovery)
        m.add_command(label="Domain-bait scan (hunt novel members)…", command=self.run_domain_scan)
        m.add_command(label="Exhaustive hunt (all sources + all analyses)…",
                      command=self.run_exhaustive)
        m.add_separator()
        m.add_command(label="Alignment & tree viewer", command=self.show_phylo_viewer)
        m.add_separator()
        m.add_command(label="Investigate accession (deep dive)…", command=self.investigate)
        m.add_command(label="Selection test (dN/dS)…", command=self.run_selection_test)
        m.add_command(label="Fold check (ESMFold)…", command=self.run_fold_check)
        m.add_command(label="Presence/absence matrix", command=self.show_presence_matrix)
        m.add_separator()
        m.add_command(label="View last analysis summary", command=self.show_analysis)
        m.add_command(label="View last discovery report", command=self.show_discovery)
        menubar.add_cascade(label="Analysis", menu=m)

    def reset(self) -> None:
        """A new search invalidates results computed for the old table."""
        self.analysis = None
        self.discovery = None
        self.presence = None
        self.novel_labels = set()

    # ---- shared helpers --------------------------------------------------

    def _guard_busy(self) -> bool:
        busy = self.runner.busy_with()
        if busy:
            messagebox.showinfo("Busy", f"A background task is already running: {busy}.")
            return True
        return False

    def _task_failed(self, tb: str) -> None:
        self.app.status("Background task failed — see details window.")
        TextWindow(self.app.root, "Task failed", tb, save_name="error.txt",
                   width=760, height=420)

    # ---- sequence analysis ----------------------------------------------

    def run_analysis(self) -> None:
        if self._guard_busy():
            return
        variants = self.app.results.all_variants()
        n_seq = sum(1 for v in variants if v.sequence)
        if not variants:
            messagebox.showinfo("Analysis", "No results to analyse — run a search first.")
            return
        if n_seq < 2:
            messagebox.showinfo(
                "Analysis",
                "Need at least 2 variants with sequences.\n"
                "Re-search with 'Fetch sequences' checked.")
            return
        dlg = AnalysisDialog(self.app.root, n_with_seq=n_seq)
        if not dlg.result:
            return
        opts = dlg.result

        def job(progress):
            progress(f"Aligning {min(n_seq, opts['max_variants'])} sequence(s) — this can take a while…")
            return analyse(
                variants,
                identity_threshold=opts["identity_threshold"],
                conservation_threshold=opts["conservation_threshold"],
                use_mafft=opts["use_mafft"],
                max_variants=opts["max_variants"],
            )

        if self.runner.run("sequence analysis", job, self._analysis_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Running sequence analysis on {n_seq} sequence(s)…")

    def _analysis_done(self, result: AnalysisResult) -> None:
        self.analysis = result
        self.app.status(
            f"Analysis done: {result.n_analyzed}/{result.n_input} variants analysed. "
            "Included in the next saved bundle.")
        self.show_analysis()

    def show_analysis(self) -> None:
        if not self.analysis:
            messagebox.showinfo("Analysis", "No analysis yet — run 'Run sequence analysis…' first.")
            return
        TextWindow(self.app.root, "Sequence analysis", self.analysis.text_summary(),
                   save_name="analysis_summary.txt",
                   extra_actions=[("Save artefacts…", self._save_analysis_artefacts)])

    def _save_analysis_artefacts(self) -> None:
        target = filedialog.askdirectory(title="Folder for the analysis/ artefacts")
        if not target or not self.analysis:
            return
        written = write_analysis(Path(target) / "analysis", self.analysis)
        messagebox.showinfo("Analysis", f"Wrote {len(written)} file(s) to {Path(target) / 'analysis'}")

    # ---- novel-paralog discovery ----------------------------------------

    def run_discovery(self) -> None:
        if self._guard_busy():
            return
        variants = self.app.results.all_variants()
        if not variants:
            messagebox.showinfo("Discovery", "No results to score — run a search first.")
            return
        # Honour the loaded preset's discovery tuning (same keys the
        # headless CLI reads), so a preset-driven method is repeatable
        # from the GUI.
        preset = getattr(self.app.search, "loaded_preset", {}) or {}
        default_known = list(preset.get("known_paralogs", [])) \
            or (list(self.app.last_query.gene_symbols) if self.app.last_query else [])
        dlg = DiscoveryDialog(self.app.root, default_known=default_known,
                              has_analysis=self.analysis is not None)
        if not dlg.result:
            return
        known = dlg.result["known_paralogs"]
        do_interpro = dlg.result["interpro"]
        config_kwargs = dict(preset.get("discovery", {}))
        analysis = self.analysis

        def job(progress):
            domain_hits: dict = {}
            if do_interpro:
                accs = sorted({v.accession for v in variants if v.source == "UniProt"})
                progress(f"InterPro: fetching Pfam signatures for {len(accs)} UniProt entrie(s)…")
                domain_hits = batch_fetch_family_signatures(accs)
            label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                         for v in variants}
            sset = None
            if analysis and analysis.msa:
                sset = build_signature_set(variants, analysis.msa, label_map, known)
                if sset:
                    progress(f"Derived {len(sset.signatures)} family-signature blocks "
                             "from the known-paralog MSA…")
            progress("Scoring candidates…")
            return discover_novel_paralogs(
                variants,
                DiscoveryConfig(known_paralogs=known, **config_kwargs),
                analysis_label_for=label_map,
                distances=(analysis.distances if analysis else None),
                domain_hits=domain_hits,
                foldseek_hits=[v for v in variants if v.source == "Foldseek"],
                signature_set=sset,
            )

        if self.runner.run("novel-paralog discovery", job, self._discovery_done,
                           self._task_failed, self.app.status):
            self.app.status("Running novel-paralog discovery…")

    def _discovery_done(self, report: DiscoveryReport) -> None:
        self.discovery = report
        self.novel_labels = self._novel_labels_from_discovery()
        bins = report.by_verdict()
        self.app.status(
            f"Discovery done: {len(bins['promising'])} promising, "
            f"{len(bins['worth manual review'])} worth review, {len(bins['weak'])} weak.")
        self.show_discovery()

    def show_discovery(self) -> None:
        if not self.discovery:
            messagebox.showinfo("Discovery", "No discovery report yet — run 'Discover novel paralogs…' first.")
            return
        TextWindow(self.app.root, "Novel-paralog discovery", self.discovery.text_summary(),
                   save_name="discovery_report.txt",
                   extra_actions=[("Save TSV + report…", self._save_discovery_artefacts)])

    def _save_discovery_artefacts(self) -> None:
        target = filedialog.askdirectory(title="Folder for the discovery/ artefacts")
        if not target or not self.discovery:
            return
        paths = write_discovery(Path(target) / "discovery", self.discovery)
        messagebox.showinfo("Discovery", "Wrote:\n" + "\n".join(paths.values()))

    # ---- selection test (dN/dS) ------------------------------------------

    def run_selection_test(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        default = ""
        if sel:
            # Compara rows carry an Ensembl protein id as accession.
            default = sel.accession if sel.accession.startswith("ENS") else \
                sel.raw.get("protein_id", "")
        dlg = SelectionTestDialog(self.app.root, default_accession=default)
        if not dlg.result:
            return
        acc, partner, k = (dlg.result["accession"], dlg.result["partner"],
                           dlg.result["max_partners"])

        def job(progress):
            return selection_test(acc, partner_id=partner, max_partners=k,
                                  log=progress)

        if self.runner.run(f"selection test on {acc}", job,
                           self._selection_done, self._task_failed,
                           lambda m: self.app.status(f"dN/dS {acc}: {m}")):
            self.app.status(f"Running selection test on {acc}…")

    def _selection_done(self, report) -> None:
        self.app.status(
            f"Selection test done: dN/dS "
            f"{'n/a' if report.result.ratio is None else f'{report.result.ratio:.2f}'}"
            f" — {report.result.verdict()}")
        TextWindow(self.app.root, f"Selection test — {report.query_id}",
                   report.summary(), save_name=f"dnds_{report.query_id}.txt",
                   width=780, height=360)

    # ---- fold check (ESMFold) --------------------------------------------

    def run_fold_check(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        if not sel:
            messagebox.showinfo("Fold check", "Select a row first.")
            return
        if not sel.sequence:
            messagebox.showinfo("Fold check",
                                "The selected row has no sequence — re-search with "
                                "'Fetch sequences' checked.")
            return
        dlg = FoldCheckDialog(self.app.root, seq_len=len(sel.sequence))
        if not dlg.result:
            return
        seq, where, length = sel.sequence, dlg.result["where"], dlg.result["length"]
        label = f"{sel.gene_symbol} ({sel.accession})"

        def job(progress):
            progress("Submitting segment to ESMFold (30–120 s)…")
            return fold_segment(seq, where=where, length=length)

        if self.runner.run(f"ESMFold check on {sel.accession}", job,
                           lambda r: self._fold_done(label, r),
                           self._task_failed, self.app.status):
            self.app.status(f"ESMFold: folding {length} aa of {label}…")

    def _fold_done(self, label: str, result) -> None:
        self.app.status(f"ESMFold done: mean pLDDT {result.mean_plddt:.1f} — "
                        f"{result.verdict()}")

        def save_pdb() -> None:
            path = filedialog.asksaveasfilename(
                defaultextension=".pdb", initialfile="esmfold_segment.pdb",
                filetypes=[("PDB", "*.pdb"), ("All files", "*.*")])
            if path:
                Path(path).write_text(result.pdb_text)

        TextWindow(self.app.root, f"ESMFold — {label}",
                   result.summary(), save_name="esmfold_check.txt",
                   extra_actions=[("Save PDB…", save_pdb)],
                   width=760, height=300)

    # ---- presence/absence matrix -----------------------------------------

    def show_presence_matrix(self) -> None:
        variants = self.app.results.all_variants()
        if not variants:
            messagebox.showinfo("Presence/absence", "No results — run a search first.")
            return
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                     for v in variants}
        report = build_presence_matrix(
            variants,
            distances=(self.analysis.distances if self.analysis else None),
            analysis_label_for=label_map,
        )
        if not self.analysis:
            report.notes.append(
                "Tip: run 'Run sequence analysis…' first so unnamed cross-species "
                "orthologs merge into homolog clusters instead of singletons.")
        TextWindow(self.app.root, "Family presence/absence matrix",
                   report.summary(), save_name="presence_matrix.txt")

    # ---- exhaustive hunt -------------------------------------------------

    def run_special_preset(self, name: str, data: dict) -> None:
        """Dispatch presets carrying a 'mode' key (loaded via the search
        panel's preset picker) to the matching tool."""
        mode = data.get("mode", "")
        if mode == "domain_scan":
            self.run_domain_scan(preset=name)
        elif mode == "exhaustive":
            self.run_exhaustive(preset=name)
        else:
            messagebox.showinfo("Preset", f"Unknown preset mode: {mode!r}")

    def run_exhaustive(self, preset: str = "") -> None:
        if self._guard_busy():
            return
        presets_dir = self.app.project_root / "presets"
        candidates = []
        for p in sorted(presets_dir.glob("*.json")):
            try:
                if json.loads(p.read_text()).get("mode") == "exhaustive":
                    candidates.append(p.stem)
            except (OSError, json.JSONDecodeError):
                continue
        if not candidates:
            messagebox.showinfo("Exhaustive hunt",
                                'No presets with "mode": "exhaustive" found in presets/.')
            return
        name = preset if preset in candidates else candidates[0]
        preset_data = json.loads((presets_dir / f"{name}.json").read_text())
        if not messagebox.askyesno(
                "Exhaustive hunt",
                f"Run the full hunt from preset '{name}'?\n\n"
                "Harvests every ITPR-like sequence (Compara mine + ortholog "
                "expansion + InterPro), saves the census, and runs every "
                "analysis. Takes 5–15 minutes against live APIs."):
            return
        project_root, email = self.app.project_root, self.app.email

        def job(progress):
            return run_exhaustive_hunt(
                project_root=project_root, preset_data=preset_data,
                preset_name=name, email=email, save_results=True,
                log=progress,
            )

        if self.runner.run("exhaustive hunt", job, self._exhaustive_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Exhaustive hunt '{name}' running — Compara mine first…")

    def _exhaustive_done(self, out: dict) -> None:
        variants = out.get("variants") or []
        if not variants:
            self.app.status("Exhaustive hunt finished with no sequences.")
            return
        self.analysis = out.get("analysis")
        self.discovery = out.get("discovery_report")
        self.presence = out.get("presence_report")
        self.novel_labels = out.get("novel_labels") or set()
        self.app.load_variants(variants, query=out.get("query"), results=out.get("results"))
        msg = (f"Exhaustive hunt done: {len(variants)} ITPR-like protein(s) in census. "
               f"Bundle: {out.get('dir', '?')}")
        self.app.status(msg)
        self.show_phylo_viewer()
        html = out.get("report_html")
        if html and messagebox.askyesno("Exhaustive hunt", msg + "\n\nOpen the HTML report?"):
            webbrowser.open(Path(html).resolve().as_uri())

    # ---- alignment & tree viewer -----------------------------------------

    def show_phylo_viewer(self) -> None:
        if not self.analysis or not self.analysis.msa:
            messagebox.showinfo(
                "Alignment & tree",
                "No alignment yet — run 'Run sequence analysis…' (or the "
                "exhaustive hunt) first.")
            return
        from .phylo_view import PhyloWindow
        PhyloWindow(self.app.root, analysis=self.analysis,
                    novel_labels=getattr(self, "novel_labels", set()) or
                    self._novel_labels_from_discovery())

    def _novel_labels_from_discovery(self) -> set[str]:
        if not self.discovery:
            return set()
        variants = self.app.results.all_variants()
        label_map = {f"{v.source}|{v.accession}|{v.gene_symbol}": _label_for(v)
                     for v in variants}
        return {label_map.get(c.label, "") for c in self.discovery.candidates
                if c.score >= 40} - {""}

    # ---- domain-bait scan ------------------------------------------------

    def run_domain_scan(self, preset: str = "") -> None:
        if self._guard_busy():
            return
        presets_dir = self.app.project_root / "presets"
        scan_presets = []
        for p in sorted(presets_dir.glob("*.json")):
            try:
                if json.loads(p.read_text()).get("mode") == "domain_scan":
                    scan_presets.append(p.stem)
            except (OSError, json.JSONDecodeError):
                continue
        if not scan_presets:
            messagebox.showinfo("Domain-bait scan",
                                'No presets with "mode": "domain_scan" found in presets/.')
            return
        dlg = DomainScanDialog(self.app.root, scan_presets, default=preset)
        if not dlg.result:
            return
        name = dlg.result["preset"]
        save = dlg.result["save_results"]
        preset_data = json.loads((presets_dir / f"{name}.json").read_text())
        project_root, email = self.app.project_root, self.app.email

        def job(progress):
            return run_domain_scan(
                project_root=project_root, preset_data=preset_data, preset_name=name,
                email=email, save_results=save, analyze=True, discover=True,
                log=progress,
            )

        if self.runner.run("domain-bait scan", job, self._domain_scan_done,
                           self._task_failed, self.app.status):
            self.app.status(f"Domain-bait scan '{name}' running — enumerating InterPro…")

    def _domain_scan_done(self, out: dict) -> None:
        variants = out.get("variants") or []
        if not variants:
            self.app.status("Domain scan finished with no candidates — see preset filters.")
            return
        self.analysis = out.get("analysis")
        self.discovery = out.get("discovery_report")
        self.app.load_variants(variants, query=out.get("query"), results=out.get("results"))
        msg = f"Domain scan done: {len(variants)} candidate(s) loaded into the table."
        if out.get("dir"):
            msg += f"  Bundle: {out['dir']}"
        self.app.status(msg)
        if self.discovery:
            self.show_discovery()
        html = out.get("report_html")
        if html and messagebox.askyesno("Domain-bait scan", msg + "\n\nOpen the HTML report?"):
            webbrowser.open(Path(html).resolve().as_uri())

    # ---- deep-dive investigation ----------------------------------------

    def investigate(self) -> None:
        if self._guard_busy():
            return
        sel = self.app.results.selected()
        dlg = InvestigateDialog(self.app.root,
                                default_accession=sel.accession if sel else "")
        if not dlg.result:
            return
        acc = dlg.result["accession"]
        panel = [(f"panel_{a}", a)
                 for a in (s.strip() for s in dlg.result["panel_csv"].split(",")) if a]
        options = InvestigationOptions(
            panel=panel or list(InvestigationOptions().panel),
            run_foldseek=dlg.result["foldseek"],
            email=self.app.email,
        )
        seq = sel.sequence if (sel and sel.accession == acc) else ""

        def job(progress):
            result = Investigator(acc, options, candidate_sequence=seq).run(on_progress=progress)
            signals = assess_signals(result)
            verdict, pct = verdict_from_signals(signals)
            return result, verdict, pct

        if self.runner.run(f"investigation of {acc}", job, self._investigation_done,
                           self._task_failed,
                           lambda msg: self.app.status(f"Investigating {acc}: {msg}")):
            self.app.status(f"Investigating {acc}…")

    def _investigation_done(self, payload) -> None:
        result, verdict, pct = payload
        self.app.status(f"Investigation of {result.accession} done — verdict: {verdict} ({pct}/100).")

        def save() -> None:
            results_root = self.app.project_root / "results"
            results_root.mkdir(exist_ok=True)
            bundle = make_bundle_dir(results_root, f"investigate_{result.accession}")
            paths = write_investigation(bundle / "investigations" / result.accession, result)
            messagebox.showinfo("Investigation", "Saved:\n" + "\n".join(paths.values()))

        TextWindow(self.app.root,
                   f"Case file — {result.accession}  ({verdict}, {pct}/100)",
                   render_case_file(result),
                   save_name=f"case_file_{result.accession}.md",
                   extra_actions=[("Save artefacts…", save)])

    # ---- bundle augmentation --------------------------------------------

    def augment_bundle(self, bundle_dir: Path, variants, query, notes: str = "") -> dict:
        """Write analysis/, discovery/, and report.md|.html into a bundle,
        mirroring the CLI's --analyze/--discover/--save-results path.
        Returns the report paths."""
        discovery_text = ""
        if self.analysis:
            write_analysis(bundle_dir / "analysis", self.analysis)
        if self.discovery:
            discovery_text = self.discovery.text_summary()
            write_discovery(bundle_dir / "discovery", self.discovery)
        return write_report(bundle_dir, variants, query,
                            analysis=self.analysis,
                            discovery_summary=discovery_text,
                            notes=notes)
