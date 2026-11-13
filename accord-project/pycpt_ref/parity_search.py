"""CPT's CCA mode search, step by step, to see which steps africas2s needs to match it.

  python parity_search.py <case>

Ranks every (x_eof, y_eof, cca) combination by cross-validated goodness, using
africas2s's CCA with the CPT-parity settings and CPT's cyclic CV, under three rules:

  A  mean Kendall tau-b over grid cells (what africas2s select_modes does)
  B  A, plus CPT's over-fitting check: an (x_eof, y_eof) pair is dropped if any
     fold has a canonical correlation of 1 (cca.F95: mu > 1 - tol, icco == 0)
  C  B, with CPT's own tau: 2 * AFC - 1, where tied pairs score one half
     (scores.F95: k_tau via two_afc_cc)

CPT's actual choice comes from find_modes.py.
"""

import sys
import time
import warnings

import numpy as np
from scipy.stats import kendalltau

warnings.filterwarnings("ignore")
from africas2s.cv import loyo_cyclic
from africas2s.methods.cca import CCAMethod
from africas2s_side import PARITY, inputs

TOL = np.sqrt(np.finfo(float).eps)  # CPT data_numbers.F95
models, obs, meta, case = inputs(sys.argv[1])
years = list(obs["year"].values)
obs_flat = obs.values.reshape(len(years), -1)
cells = np.flatnonzero(~np.isnan(obs_flat).any(axis=0))


def tau_b(pred):
    taus = [kendalltau(pred[:, j], obs_flat[:, j]).statistic for j in cells]
    return float(np.nanmean(taus))


def tau_cpt(pred):
    """Mean over cells of CPT's k_tau: pairs tied in obs or forecast score 0.5."""
    total = 0.0
    for j in cells:
        do = np.sign(obs_flat[:, None, j] - obs_flat[None, :, j])
        df = np.sign(pred[:, None, j] - pred[None, :, j])
        iu = np.triu_indices(len(years), 1)
        do, df = do[iu], df[iu]
        afc = np.where((do == 0) | (df == 0), 0.5, np.abs(do + df) / 2.0).mean()
        total += 2.0 * afc - 1.0
    return total / len(cells)


for name, (hcst, _) in models.items():
    t0 = time.time()
    combos = [(xe, ye, cc) for xe in range(1, 11) for ye in range(1, 11)
              for cc in range(1, min(xe, ye) + 1)]
    preds = np.full((len(combos), len(years), obs_flat.shape[1]), np.nan)
    overfit = set()
    for train, test in loyo_cyclic(years, window=5):
        i = years.index(test)
        x_train, y_train = hcst.sel(year=train), obs.sel(year=train)
        x_test = hcst.sel(year=[test]).isel(year=0, drop=True)
        for ci, (xe, ye, cc) in enumerate(combos):
            try:
                m = CCAMethod(x_eof_modes=xe, y_eof_modes=ye, cca_modes=cc, **PARITY)
                m.fit(x_train, y_train)
                preds[ci, i] = m.predict(x_test).mean("member").values.reshape(-1)
                if cc == min(xe, ye) and np.any(m.mu_ > 1 - TOL):
                    overfit.add((xe, ye))
            except Exception:
                pass
    ok = [ci for ci in range(len(combos)) if np.isfinite(preds[ci][:, cells]).all()]

    def best(score, drop_overfit):
        rows = sorted(((score(preds[ci]), combos[ci]) for ci in ok
                       if not (drop_overfit and combos[ci][:2] in overfit)), reverse=True)
        return ", ".join(f"{c} {g:+.4f}" for g, c in rows[:3])

    print(f"{name}: {len(overfit)} of 100 (x, y) pairs over-fitted   ({time.time() - t0:.0f} s)")
    print(f"   A tau-b, no check    : {best(tau_b, False)}")
    print(f"   B tau-b + check      : {best(tau_b, True)}")
    print(f"   C CPT tau + check    : {best(tau_cpt, True)}", flush=True)
