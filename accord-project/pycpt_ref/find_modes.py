"""Infer the EOF truncations CPT selected from pyCPT's output.

  python find_modes.py <case>

pyCPT's output shows the number of CCA modes but writes out every EOF it
computed. For each model, fit africas2s (parity settings) with every
(x_eof, y_eof) pair at CPT's CCA count and report the pair whose cross-validated
hindcast is closest to pyCPT's. A near-zero difference identifies CPT's choice
and shows the two implementations agree given the same modes.
"""

import sys
import warnings

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore")
import africas2s
from africas2s_side import PARITY, inputs

models, obs, meta, case = inputs(sys.argv[1])
clim = obs.mean("year")
wet = clim >= 0.1 * meta["season_days"]

for name, pair in models.items():
    px = xr.open_dataset(case / "out" / f"px_{name}.nc", decode_times=False)
    scores = px.x_cca_scores.where(px.x_cca_scores > -900)
    cca = int(scores.notnull().any("T").sum())
    ref = xr.open_dataset(case / "out" / f"hcst_{name}.nc", decode_times=False)["deterministic"]
    ref = ref.where(ref > -900).rename(Y="lat", X="lon", T="year").sortby("lat")
    ref = ref.assign_coords(year=obs["year"].values).where(wet)
    scale = float(np.sqrt((ref ** 2).mean()))
    results = []
    for xe in range(cca, 11):
        for ye in range(cca, 11):
            args = dict(x_eof_modes=xe, y_eof_modes=ye, cca_modes=cca,
                        crossvalidation_window=5, **PARITY)
            try:
                res = africas2s.seasonal_mme({"precip": {name: pair}}, obs, method="cca",
                                             cpt_args=args, cv="loyo_cyclic",
                                             probability_aggregation="cpt_per_model",
                                             forecast_year=meta["forecast_year"], verbose=False)
            except Exception:
                continue
            cv = next(iter(res.per_model_cv_hindcasts.values())) - clim
            rel = float(np.sqrt(((cv - ref) ** 2).mean())) / scale
            results.append((rel, xe, ye))
    results.sort()
    best = ", ".join(f"x={x} y={y}: {r:.4f}" for r, x, y in results[:4])
    print(f"{name}: CCA modes {cca}. Closest (RMS diff / RMS pyCPT): {best}", flush=True)
    print(f"MODES {name}={results[0][1]},{results[0][2]},{cca}", flush=True)
