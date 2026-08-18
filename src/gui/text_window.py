"""Reusable scrollable text window for tool output.

Shows monospace text (analysis summaries, discovery reports, case files) in
a Toplevel with a "Save as…" button and optional extra action buttons.
"""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Optional


class TextWindow(tk.Toplevel):
    def __init__(
        self,
        master: tk.Misc,
        title: str,
        text: str,
        save_name: str = "output.txt",
        extra_actions: Optional[list[tuple[str, Callable[[], None]]]] = None,
        width: int = 920,
        height: int = 620,
    ) -> None:
        super().__init__(master)
        self.title(title)
        self.geometry(f"{width}x{height}")
        self._text = text
        self._save_name = save_name

        body = ttk.Frame(self, padding=8)
        body.pack(fill="both", expand=True)

        txt = tk.Text(body, wrap="none", font=("Menlo", 11))
        vsb = ttk.Scrollbar(body, orient="vertical", command=txt.yview)
        hsb = ttk.Scrollbar(body, orient="horizontal", command=txt.xview)
        txt.configure(yscrollcommand=vsb.set, xscrollcommand=hsb.set)
        txt.grid(row=0, column=0, sticky="nsew")
        vsb.grid(row=0, column=1, sticky="ns")
        hsb.grid(row=1, column=0, sticky="we")
        body.rowconfigure(0, weight=1)
        body.columnconfigure(0, weight=1)
        txt.insert("1.0", text)
        txt.configure(state="disabled")

        btns = ttk.Frame(self, padding=(8, 0, 8, 8))
        btns.pack(fill="x")
        for label, cb in (extra_actions or []):
            ttk.Button(btns, text=label, command=cb).pack(side="left", padx=(0, 6))
        ttk.Button(btns, text="Close", command=self.destroy).pack(side="right")
        ttk.Button(btns, text="Save as…", command=self._save).pack(side="right", padx=6)

    def _save(self) -> None:
        path = filedialog.asksaveasfilename(
            parent=self,
            initialfile=self._save_name,
            defaultextension=".txt",
            filetypes=[("Text", "*.txt"), ("Markdown", "*.md"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            with open(path, "w") as f:
                f.write(self._text)
        except OSError as e:
            messagebox.showerror("Save failed", str(e), parent=self)
