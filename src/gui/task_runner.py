"""Generic background-task bridge for the GUI.

Long-running work (network fetches, MSA building, report figures) must never
run on the Tk main thread. `TaskRunner` runs one job at a time on a worker
thread and marshals progress / completion back to the main thread via a
`queue.Queue` drained with `root.after()` — the same pattern the search
orchestrator uses.
"""

from __future__ import annotations

import queue
import threading
import tkinter as tk
import traceback
from typing import Any, Callable, Optional

POLL_MS = 100


class TaskRunner:
    """Runs one named background job at a time for a Tk application."""

    def __init__(self, root: tk.Misc) -> None:
        self.root = root
        self._queue: "queue.Queue[tuple[str, Any]]" = queue.Queue()
        self._task_name: Optional[str] = None
        self._on_done: Optional[Callable[[Any], None]] = None
        self._on_error: Optional[Callable[[str], None]] = None
        self._on_progress: Optional[Callable[[str], None]] = None
        self.root.after(POLL_MS, self._poll)

    def busy_with(self) -> Optional[str]:
        """Name of the running task, or None when idle."""
        return self._task_name

    def run(
        self,
        name: str,
        fn: Callable[[Callable[[str], None]], Any],
        on_done: Callable[[Any], None],
        on_error: Optional[Callable[[str], None]] = None,
        on_progress: Optional[Callable[[str], None]] = None,
    ) -> bool:
        """Start `fn(progress)` on a worker thread. Returns False if busy.

        `fn` receives a `progress(message)` callable that is safe to call
        from the worker; messages arrive on the main thread via `on_progress`.
        """
        if self._task_name is not None:
            return False
        self._task_name = name
        self._on_done = on_done
        self._on_error = on_error
        self._on_progress = on_progress

        def progress(message: str) -> None:
            self._queue.put(("progress", str(message)))

        def worker() -> None:
            try:
                result = fn(progress)
            except Exception:
                self._queue.put(("error", traceback.format_exc(limit=8)))
            else:
                self._queue.put(("done", result))

        threading.Thread(target=worker, daemon=True, name=f"task-{name}").start()
        return True

    def _poll(self) -> None:
        try:
            while True:
                kind, payload = self._queue.get_nowait()
                if kind == "progress":
                    if self._on_progress:
                        self._on_progress(payload)
                    continue
                on_done, on_error = self._on_done, self._on_error
                self._task_name = None
                self._on_done = self._on_error = self._on_progress = None
                if kind == "done" and on_done:
                    on_done(payload)
                elif kind == "error" and on_error:
                    on_error(payload)
        except queue.Empty:
            pass
        self.root.after(POLL_MS, self._poll)
