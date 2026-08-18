"""Alignment & phylogenetics viewer.

A Toplevel with two panes built from the cached `AnalysisResult`:

    * top — the NJ tree drawn as a phylogram on a Canvas (x = cumulative
      branch length, tips evenly spaced). Novel-candidate tips (the
      discovery scorer's ≥40-point labels) are red/bold.
    * bottom — an alignment overview: one row per tree tip, columns binned
      to fit the width, colored by per-column conservation (dark = highly
      conserved, light = variable, white = gap).

Buttons: save the tree as PNG (matplotlib, same highlighting), save the
alignment FASTA, save the Newick.
"""

from __future__ import annotations

import tkinter as tk
from pathlib import Path
from tkinter import filedialog, messagebox, ttk

from ..analysis.figures import HAVE_MPL, tree_figure
from ..analysis.pipeline import AnalysisResult
from ..analysis.tree import build_nj_tree

HIGHLIGHT = "#c0392b"
NORMAL = "#1c2733"
BRANCH = "#5a6b7a"


class PhyloWindow(tk.Toplevel):
    def __init__(self, master: tk.Misc, analysis: AnalysisResult,
                 novel_labels: set[str] | None = None) -> None:
        super().__init__(master)
        self.title("Alignment & tree")
        self.geometry("1100x760")
        self.analysis = analysis
        self.novel = novel_labels or set()
        self._build()

    # ---- layout ----------------------------------------------------------

    def _build(self) -> None:
        info = ttk.Frame(self, padding=(8, 6))
        info.pack(fill="x")
        n_novel = len(self.novel)
        ttk.Label(info, text=(
            f"{self.analysis.n_analyzed} sequences aligned · "
            f"{len(self.analysis.msa[0].aligned) if self.analysis.msa else 0} columns · "
            f"{n_novel} novel candidate(s) highlighted in red"
        )).pack(side="left")
        ttk.Button(info, text="Save Newick…", command=self._save_newick).pack(side="right", padx=2)
        ttk.Button(info, text="Save alignment FASTA…", command=self._save_fasta).pack(side="right", padx=2)
        ttk.Button(info, text="Save tree PNG…", command=self._save_png).pack(side="right", padx=2)

        paned = ttk.PanedWindow(self, orient="vertical")
        paned.pack(fill="both", expand=True)

        tree_frame = ttk.Frame(paned)
        self.tree_canvas = tk.Canvas(tree_frame, bg="white", highlightthickness=0)
        tvsb = ttk.Scrollbar(tree_frame, orient="vertical", command=self.tree_canvas.yview)
        self.tree_canvas.configure(yscrollcommand=tvsb.set)
        self.tree_canvas.pack(side="left", fill="both", expand=True)
        tvsb.pack(side="right", fill="y")
        paned.add(tree_frame, weight=3)

        aln_frame = ttk.Frame(paned)
        self.aln_canvas = tk.Canvas(aln_frame, bg="white", highlightthickness=0)
        avsb = ttk.Scrollbar(aln_frame, orient="vertical", command=self.aln_canvas.yview)
        self.aln_canvas.configure(yscrollcommand=avsb.set)
        self.aln_canvas.pack(side="left", fill="both", expand=True)
        avsb.pack(side="right", fill="y")
        paned.add(aln_frame, weight=2)

        self.tip_order = self._draw_tree()
        self._draw_alignment(self.tip_order)

    # ---- tree drawing ----------------------------------------------------

    def _draw_tree(self) -> list[str]:
        """Phylogram on the canvas; returns tip labels in display order."""
        c = self.tree_canvas
        if not self.analysis.distances:
            c.create_text(20, 20, anchor="nw", text="(no distance matrix)")
            return []
        tree = build_nj_tree(self.analysis.distances)
        if tree is None:
            c.create_text(20, 20, anchor="nw", text="(tree unavailable)")
            return []
        tips = tree.get_terminals()
        row_h, pad, label_w = 22, 30, 340
        width = 1050
        depths = tree.depths()
        if not any(depths.values()):                    # zero branch lengths
            depths = tree.depths(unit_branch_lengths=True)
        max_depth = max(depths.values()) or 1.0
        x_scale = (width - label_w - 2 * pad) / max_depth

        ys: dict = {}
        for i, tip in enumerate(tips):
            ys[tip] = pad + i * row_h

        def assign_y(clade) -> float:
            if clade in ys:
                return ys[clade]
            child_ys = [assign_y(ch) for ch in clade.clades]
            ys[clade] = sum(child_ys) / len(child_ys)
            return ys[clade]

        assign_y(tree.root)

        def draw(clade, x_parent: float) -> None:
            x = pad + depths[clade] * x_scale
            y = ys[clade]
            c.create_line(x_parent, y, x, y, fill=BRANCH, width=1.4)
            if clade.is_terminal():
                label = str(clade.name or "")
                novel = label in self.novel
                c.create_text(
                    x + 6, y, anchor="w", text=("★ " if novel else "") + label,
                    fill=HIGHLIGHT if novel else NORMAL,
                    font=("Menlo", 10, "bold" if novel else "normal"),
                )
            else:
                child_ys = [ys[ch] for ch in clade.clades]
                c.create_line(x, min(child_ys), x, max(child_ys), fill=BRANCH, width=1.4)
                for ch in clade.clades:
                    draw(ch, x)

        draw(tree.root, pad)
        c.configure(scrollregion=(0, 0, width, pad * 2 + len(tips) * row_h))
        return [str(t.name or "") for t in tips]

    # ---- alignment overview ----------------------------------------------

    def _draw_alignment(self, tip_order: list[str]) -> None:
        c = self.aln_canvas
        rows = {r.label: r.aligned for r in self.analysis.msa}
        ordered = [l for l in tip_order if l in rows] or list(rows)
        if not ordered:
            c.create_text(20, 20, anchor="nw", text="(no alignment)")
            return
        n_cols = len(next(iter(rows.values())))
        cons = self.analysis.conservation or [0.5] * n_cols
        n_bins = min(160, n_cols)
        bin_w, row_h, label_w, pad = 5, 14, 340, 8
        per_bin = max(1, n_cols // n_bins)

        # Per-bin conservation → color ramp (light → dark blue).
        def shade(value: float) -> str:
            v = max(0.0, min(1.0, value))
            r = int(232 - 160 * v)
            g = int(240 - 140 * v)
            b = int(248 - 90 * v)
            return f"#{r:02x}{g:02x}{b:02x}"

        c.create_text(label_w + pad, 4, anchor="nw", fill="#66727e",
                      font=("Menlo", 9),
                      text="alignment overview — dark = conserved, white = gap")
        for ri, label in enumerate(ordered):
            y = 22 + ri * row_h
            novel = label in self.novel
            c.create_text(label_w - 4, y + row_h / 2, anchor="e", text=label,
                          fill=HIGHLIGHT if novel else NORMAL,
                          font=("Menlo", 9, "bold" if novel else "normal"))
            aligned = rows[label]
            for b in range(n_bins):
                start = b * per_bin
                seg = aligned[start:start + per_bin]
                if not seg or all(ch == "-" for ch in seg):
                    continue
                mean_cons = sum(cons[start:start + per_bin]) / max(1, len(seg))
                occupancy = 1 - seg.count("-") / len(seg)
                c.create_rectangle(
                    label_w + pad + b * bin_w, y,
                    label_w + pad + (b + 1) * bin_w - 1, y + row_h - 2,
                    fill=shade(mean_cons * occupancy), outline="")
        c.configure(scrollregion=(0, 0, label_w + pad * 2 + n_bins * bin_w,
                                  30 + len(ordered) * row_h))

    # ---- save actions ----------------------------------------------------

    def _save_png(self) -> None:
        if not HAVE_MPL:
            messagebox.showinfo("Save tree", "matplotlib not installed.", parent=self)
            return
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".png", initialfile="tree_highlighted.png",
            filetypes=[("PNG", "*.png")])
        if not path:
            return
        tree = build_nj_tree(self.analysis.distances) if self.analysis.distances else None
        out = tree_figure(tree, Path(path), title="Family NJ tree",
                          highlight_labels=self.novel)
        if out:
            messagebox.showinfo("Save tree", f"Wrote {out}", parent=self)

    def _save_fasta(self) -> None:
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".fasta", initialfile="alignment.fasta",
            filetypes=[("FASTA", "*.fasta")])
        if not path:
            return
        with open(path, "w") as f:
            for r in self.analysis.msa:
                f.write(f">{r.label}\n{r.aligned}\n")
        messagebox.showinfo("Save alignment", f"Wrote {len(self.analysis.msa)} rows.", parent=self)

    def _save_newick(self) -> None:
        path = filedialog.asksaveasfilename(
            parent=self, defaultextension=".newick", initialfile="tree.newick",
            filetypes=[("Newick", "*.newick"), ("All files", "*.*")])
        if not path:
            return
        Path(path).write_text(self.analysis.newick + "\n")
        messagebox.showinfo("Save Newick", f"Wrote {path}", parent=self)
