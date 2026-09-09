"""s19_lib — S19's shared surface.

Two modules behind one name, because the whole task imports `s19_lib` and the
500-line budget does not fit both halves: `s19_lib_base` carries the paths,
the control definition, the ledger join, the census and the small statistics,
and `s19_universe` carries the swept accession index, which is the one part
that reads bulk data.
"""

from __future__ import annotations

from s19_lib_base import *                        # noqa: F401,F403
from s19_lib_base import (ALL_GROUPS, BAIT_MANIFEST, CELLS,   # noqa: F401
                          CENSUS, CONTROL_CELL, FOUND_STATUSES, HMM_DIR,
                          JACK_RUNS, LEDGER, MANIFEST, MATRIX, METHOD_LABEL,
                          METHOD_ORDER, METHOD_UNIVERSE, NONVERT_GROUPS,
                          OUT_DIR, PARALOGS, PRESENT_STATES, PROJECT_ROOT,
                          RESULTS, S20_DIR, S23_DIR, UNDECIDABLE_STATES,
                          V1_DELTA, acc_key, bait_meta, cache_dir,
                          census_rows, contiguity_bar, control_cells,
                          convergence_rounds, data_root, fmt, fnum,
                          is_found, jackhmmer_log, jackhmmer_targets,
                          jackhmmer_verdicts, ledger_cells, log, median,
                          out_dir, read_json, read_tsv, sweep_dir, wilson,
                          write_json, write_tsv)
from s19_universe import (build_universe, header_accession,   # noqa: F401
                          universe_intersect, universe_size, universe_taxids)
