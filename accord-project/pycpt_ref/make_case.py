"""Write a comparison case: model and observation fields in pyCPT's input format.

  python make_case.py horn ../configs/horn_c3s.yml     # from an osop_acmad config
  python make_case.py sahara --osop-test               # OSOP's bundled converter outputs

Creates cases/<name>/inputs/*.nc and cases/<name>/case.json. Model fields are
ensemble means, as OSOP hands them to pyCPT; precipitation is a season total in mm.
"""

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))


def to_pycpt(da, cfg, days, name):
    """(year, lat, lon) season mean in mm/day -> pyCPT layout (T, Y, X) season total in mm."""
    first, last = cfg.target_months[0], cfg.target_months[-1]
    starts, ends, inits = [], [], []
    for year in da["year"].values:
        y0 = int(year) + cfg.season_year_offset
        start = pd.Timestamp(y0, first, 1)
        end = pd.Timestamp(y0 + (1 if last < first else 0), last, 1) + pd.offsets.MonthEnd(0)
        starts.append(start), ends.append(end), inits.append(pd.Timestamp(int(year), cfg.month, 1))
    starts, ends = pd.DatetimeIndex(starts), pd.DatetimeIndex(ends)
    out = (da * days).rename(year="T", lat="Y", lon="X").sortby("Y", ascending=False)
    out = out.assign_coords(T=starts + (ends - starts) / 2)
    out = out.assign_coords(Ti=("T", starts), Tf=("T", ends))
    if name != "prcp":  # models carry their initialisation date; observations must not
        out = out.assign_coords(S=("T", pd.DatetimeIndex(inits)))
    out["Y"].attrs = {"units": "degrees_north", "standard_name": "latitude"}
    out["X"].attrs = {"units": "degrees_east", "standard_name": "longitude"}
    out.attrs = {"units": "mm", "missing": -999.0}
    return out.transpose("T", "Y", "X").astype("float64").rename(name)


def write(da, path):
    enc = {da.name: {"_FillValue": np.float64(-999.0), "dtype": "float64"}}
    da.to_dataset().to_netcdf(path, encoding=enc)


def from_config(case, config_path):
    from osop_acmad.config import load_config
    from osop_acmad.data import _season_days, load_inputs

    cfg = load_config(config_path)
    inp = load_inputs(cfg)
    days = _season_days(cfg)
    inputs = case / "inputs"
    write(to_pycpt(inp.obs, cfg, days, "prcp"), inputs / "obs.nc")
    models = []
    for centre, (hcst, fcst) in inp.predictors.items():
        name = centre.upper()
        write(to_pycpt(hcst.mean("member"), cfg, days, "aprod"), inputs / f"{name}_hcst.nc")
        write(to_pycpt(fcst.mean("member"), cfg, days, "aprod"), inputs / f"{name}_fcst.nc")
        models.append(name)
    n, w, s, e = cfg.predictor_area
    y0, y1 = cfg.hindcast_years
    return dict(models=models, predictor_extent=dict(west=w, east=e, south=s, north=n),
                first_year=y0, last_year=y1, forecast_year=cfg.forecast_year,
                season=cfg.season, season_days=days, source=str(config_path), sources=inp.sources)


def from_osop_test(case):
    data = HERE.parents[1] / "lib/osop/test/test_data/test_pycptf"
    tail = "monthly_mean_5_234_40.10.10.40_total_precipitation.pycpt.nc"
    inputs = case / "inputs"
    shutil.copy(data / "era5_total_precipitation_1993-2016_monthly_5_123_40.10.10.40.pycpt.nc",
                inputs / "obs.nc")
    for name, stem in {"ECMWF": "pycpt_ecmwf_51", "UKMO": "pycpt_ukmo_604"}.items():
        shutil.copy(data / f"{stem}_1993-2016_{tail}", inputs / f"{name}_hcst.nc")
        shutil.copy(data / f"{stem}_2025-2025_{tail}", inputs / f"{name}_fcst.nc")
    return dict(models=["ECMWF", "UKMO"],
                predictor_extent=dict(west=10.0, east=40.0, south=10.0, north=40.0),
                first_year=1993, last_year=2016, forecast_year=2025, season="JJA",
                season_days=92, source="OSOP bundled test data (pycpt_convert.py output)")


def main():
    name, arg = sys.argv[1], sys.argv[2]
    case = HERE / "cases" / name
    (case / "inputs").mkdir(parents=True, exist_ok=True)
    meta = from_osop_test(case) if arg == "--osop-test" else from_config(case, arg)
    (case / "case.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote {case}: {meta['models']} {meta['season']}")


if __name__ == "__main__":
    main()
