"""africas2s.select_modes vs the modes CPT chose (find_modes.py), per case.

  python check_select_modes.py <case>

Runs the search with the fit's parity settings and CPT's cyclic CV, with
tau-b and with CPT's tau; the over-fitting check is on in both (the default).
"""
import sys
import time
import warnings

warnings.filterwarnings("ignore")
from africas2s.cv import loyo_cyclic
from africas2s.methods.cca import select_modes
from africas2s_side import PARITY, inputs

models, obs, meta, case = inputs(sys.argv[1])
years = [int(y) for y in obs["year"].values]
for name, (hcst, _) in models.items():
    for tau in ("b", "cpt"):
        t0 = time.time()
        xe, ye, cc, g, *_ = select_modes(hcst, obs, years, window=5, method_kwargs=PARITY,
                                         cv=loyo_cyclic, tau=tau)
        print(f"RESULT {sys.argv[1]} {name} tau={tau}: ({xe}, {ye}, {cc})  goodness {g:+.4f}"
              f"  ({time.time() - t0:.0f} s)", flush=True)
