"""Compare the pyCPT reference with africas2s on one case.

  python compare.py <case>

Reads cases/<case>/out (pyCPT) and cases/<case>/africas2s_out (every tag found).
Cells with an observed season under 0.1 mm/day are left out.
"""

import json
import sys
import warnings

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore")
import africas2s
from africas2s_side import inputs

models, obs, meta, case = inputs(sys.argv[1])
YEARS = obs["year"].values
clim = obs.mean("year")
wet = clim >= 0.1 * meta["season_days"]
obs_anom = obs - clim


def cpt(name):
    """A pyCPT output file with africas2s dimension names and NaN for CPT's -999."""
    ds = xr.open_dataset(case / "out" / f"{name}.nc", decode_times=False)
    ds = ds.where(ds > -900).rename(Y="lat", X="lon").sortby("lat")
    if "T" in ds.dims:
        ds = ds.assign_coords(
            T=YEARS[: ds.sizes["T"]] if ds.sizes["T"] > 1 else [meta["forecast_year"]]).rename(T="year")
    if "C" in ds.dims:
        ds = ds.assign_coords(C=[0, 1, 2]).rename(C="tercile")
    return ds


def n_modes(da):
    da = da.where(da > -900)
    return int(da.notnull().any([d for d in da.dims if d != "Mode"]).sum())


def cpt_modes(name):
    px = xr.open_dataset(case / "out" / f"px_{name}.nc", decode_times=False)
    py = xr.open_dataset(case / "out" / f"py_{name}.nc", decode_times=False)
    return [n_modes(px.x_eof_scores), n_modes(py.y_eof_scores), n_modes(px.x_cca_scores)]


def cell_corr(a, b):
    a, b = xr.align(a, b, join="inner")
    return float(xr.corr(a, b, dim="year").where(wet).mean())


def prob_scores(probs):
    probs = probs.transpose("year", "tercile", "lat", "lon").where(wet)
    rep = africas2s.skill(probs, obs, metrics=["rpss", "generalized_roc"], spatial=True)
    return (float(rep.spatial["rpss"].where(wet).mean()),
            float(rep.spatial["generalized_roc"].where(wet).mean()))


tags = sorted(p.name[: -len("_modes.json")] for p in (case / "africas2s_out").glob("*_modes.json"))
print(f"case {sys.argv[1]}: {meta['season']} {meta['first_year']}-{meta['last_year']}, "
      f"forecast {meta['forecast_year']}, {int(wet.sum())} wet cells of {wet.size}")

for name in models:
    ref_h, ref_f = cpt(f"hcst_{name}"), cpt(f"fcst_{name}")
    ref_det = ref_h["deterministic"]
    ref_fc = ref_f["deterministic"].squeeze("year", drop=True).where(wet)
    ref_pf = (ref_f["probabilistic"].squeeze("year", drop=True) / 100.0).where(wet).transpose("tercile", "lat", "lon")
    rp, rg = prob_scores(ref_h["probabilistic"] / 100.0)
    print(f"\n=== {name}   pyCPT CCA modes = {cpt_modes(name)[2]} (EOF counts are not reported; see find_modes.py)")
    print(f"  pyCPT: CV hindcast Pearson with obs {cell_corr(ref_det, obs_anom):+.3f}"
          f" | in-sample probabilistic RPSS {rp:+.3f} GROC {rg:.3f}")
    for tag in tags:
        info = json.loads((case / "africas2s_out" / f"{tag}_modes.json").read_text())
        get = lambda kind: xr.open_dataarray(case / "africas2s_out" / f"{tag}_{kind}_{name}.nc")
        cv_anom = get("hcst") - clim           # pyCPT returns anomalies (tailoring="Anomaly")
        fc_anom = (get("fcst") - clim).where(wet)
        pr_fc = get("pr_fcst").transpose("tercile", "lat", "lon").where(wet)
        ap, ag = prob_scores(get("pr_hcst"))
        dp = np.abs(pr_fc - ref_pf)
        same = float((pr_fc.fillna(-1).argmax("tercile") == ref_pf.fillna(-1).argmax("tercile")).where(wet).mean())
        diff = (cv_anom - ref_det).where(wet)
        rel = float(np.sqrt((diff ** 2).mean()) / np.sqrt((ref_det.where(wet) ** 2).mean()))
        print(f"  africas2s [{tag}] modes {info['modes'][name]} parity={info['parity']}:")
        print(f"      CV hindcast: corr with pyCPT {cell_corr(cv_anom, ref_det):.4f}, RMS diff / RMS pyCPT {rel:.3f},"
              f" Pearson with obs {cell_corr(cv_anom, obs_anom):+.3f}, CV RPSS {ap:+.3f} GROC {ag:.3f}")
        print(f"      {meta['forecast_year']} forecast: deterministic pattern corr {float(xr.corr(fc_anom, ref_fc)):.4f},"
              f" probabilities mean abs diff {float(dp.mean()):.4f} max {float(dp.max()):.3f},"
              f" same leading category {same:.0%}")
