"""africas2s CCA on the same input files the pyCPT reference uses.

  python africas2s_side.py <case> <tag> [--no-parity] [MODEL=x_eof,y_eof,cca ...]

MODEL=... fixes the modes for that model; without it africas2s's CPT-style
search chooses them. The CPT-parity settings are on unless --no-parity.
Writes cases/<case>/africas2s_out/<tag>_*.nc.
"""

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

warnings.filterwarnings("ignore")
import africas2s

HERE = Path(__file__).resolve().parent
PARITY = dict(standardize=True, lat_weights="cpt", sign_convention="cpt", y_reconstruction="cpt")


def load(path, first_year):
    """A pyCPT-format file as (year, lat, lon), years counted from ``first_year``."""
    ds = xr.open_dataset(path, decode_times=False)
    da = list(ds.data_vars.values())[0]
    da = da.rename(T="year", Y="lat", X="lon").assign_coords(
        year=first_year + np.arange(da.sizes["T"]))
    da = da.drop_vars([c for c in da.coords if c not in da.dims]).sortby("lat")
    return da.where(da > -900)


def inputs(case):
    case = HERE / "cases" / case
    meta = json.loads((case / "case.json").read_text())
    y0, y1, fy = meta["first_year"], meta["last_year"], meta["forecast_year"]
    obs = load(case / "inputs" / "obs.nc", y0).sel(year=slice(y0, y1))
    models = {}
    for name in meta["models"]:
        hcst = load(case / "inputs" / f"{name}_hcst.nc", y0).expand_dims(member=[0])
        fcst = load(case / "inputs" / f"{name}_fcst.nc", fy).expand_dims(member=[0])
        models[name] = (hcst.transpose("year", "member", "lat", "lon"),
                        fcst.transpose("year", "member", "lat", "lon"))
    return models, obs, meta, case


def main():
    case_name, tag = sys.argv[1], sys.argv[2]
    parity = "--no-parity" not in sys.argv
    fixed = dict(arg.split("=") for arg in sys.argv[3:] if "=" in arg)
    models, obs, meta, case = inputs(case_name)
    out = case / "africas2s_out"
    out.mkdir(exist_ok=True)
    chosen = {}
    for name, pair in models.items():
        if name in fixed:
            xe, ye, cc = (int(v) for v in fixed[name].split(","))
            args = dict(x_eof_modes=xe, y_eof_modes=ye, cca_modes=cc, crossvalidation_window=5)
        else:
            args = dict(mode_selection="cpt", x_eof_range=(1, 10), y_eof_range=(1, 10),
                        cca_range=(1, 10), crossvalidation_window=5)
        if parity:
            args.update(PARITY)
        res = africas2s.seasonal_mme({"precip": {name: pair}}, obs, method="cca", cpt_args=args,
                                     cv="loyo_cyclic" if parity else "loyo",
                                     probability_aggregation="cpt_per_model",
                                     forecast_year=meta["forecast_year"], verbose=False)
        next(iter(res.per_model_cv_hindcasts.values())).to_netcdf(out / f"{tag}_hcst_{name}.nc")
        next(iter(res.per_model_forecasts.values())).to_netcdf(out / f"{tag}_fcst_{name}.nc")
        res.tercile_forecast.to_netcdf(out / f"{tag}_pr_fcst_{name}.nc")
        res.tercile_cv.to_netcdf(out / f"{tag}_pr_hcst_{name}.nc")
        method = next(iter(res.per_model_methods.values()), None)
        chosen[name] = [int(getattr(method, a, -1) or -1)
                        for a in ("x_eof_modes", "y_eof_modes", "cca_modes")]
        print(f"{tag} {name}: modes {chosen[name]} parity={parity}", flush=True)
    (out / f"{tag}_modes.json").write_text(json.dumps(dict(modes=chosen, parity=parity)))


if __name__ == "__main__":
    main()
