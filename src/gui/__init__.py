"""tkinter/ttk GUI for the Protein Variant Finder.

Contents:
    app.py          — MainWindow (Tk root), menus, status bar, queue polling
    search_panel.py — gene + species + DB checkbox input panel
    results_view.py — sortable/filterable Treeview of ProteinVariant rows
    details_view.py — read-only detail pane for the selected row
    tools.py        — ToolsController: Analysis menu (analysis, discovery,
                      domain scan, investigation, bundle/report augmentation)
    dialogs.py      — modal option dialogs for the Analysis-menu tools
    task_runner.py  — one-at-a-time worker-thread bridge (queue + after())
    text_window.py  — reusable scrollable output window with Save-as
"""

from .app import MainWindow, run

__all__ = ["MainWindow", "run"]
