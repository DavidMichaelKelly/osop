"""Load model hindcasts, forecasts and observations as season means.

Two backends return the same shapes:

* ``acmaddl``: fetch from source through acmadDL (C3S models need a CDS key).
* ``grib``: read the GRIB files an OSOP run leaves in its downloads directory,
  so an existing OSOP download can be reused with no network at all.

Shapes: hindcast ``(year, member, lat, lon)``, forecast the same with a single
year, observations ``(year, lat, lon)``. ``year`` is always the initialisation year,
also for seasons that run into the next calendar year. Latitude ascends.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from .config import Config, ConfigError


@dataclass
class Inputs:
    predictors: dict  # centre -> (hindcast, forecast or None) over the predictor box
    targets: dict     # centre -> (hindcast, forecast or None) over the target box
    obs: xr.DataArray
    sources: dict     # centre / "obs" -> where the data came from


def load_inputs(cfg: Config) -> Inputs:
    return _load_grib(cfg) if cfg.backend == "grib" else _load_acmaddl(cfg)


def crop(da: xr.DataArray, region) -> xr.DataArray:
    """Cut ``da`` to ``[south, north, west, east]``, switching longitude convention if needed."""
    south, north, west, east = region
    lon = da["lon"]
    if west < 0 and float(lon.max()) > 180:
        da = da.assign_coords(lon=((lon + 180) % 360) - 180).sortby("lon")
    elif east > 180 and float(lon.min()) < 0:
        da = da.assign_coords(lon=lon % 360).sortby("lon")
    da = da.sortby("lat")
    return da.sel(lat=slice(south, north), lon=slice(west, east))


def _contains(outer, inner) -> bool:
    os_, on, ow, oe = outer
    is_, in_, iw, ie = inner
    lon_ok = (ow <= iw and ie <= oe) or (ow == 0 and oe >= 359)
    return os_ <= is_ and in_ <= on and lon_ok


# --------------------------------------------------------------------------- acmadDL

def product_for(cfg: Config, centre: str) -> str:
    """acmadDL product for an OSOP centre + CDS system number.

    OSOP names the system in its config; acmadDL pins one system per catalog
    entry. So look the pair up, and fail loudly when the catalog has no entry
    for that system instead of quietly fetching a different one.
    """
    if centre in cfg.products:
        return cfg.products[centre]
    from acmaddl import catalog

    seen = {}
    for name in catalog.list_products(include_deprecated=False):
        try:
            entry = catalog.get(name)
        except KeyError:
            continue
        if (entry.get("adapter") != "cds"
                or entry.get("cds_dataset") != "seasonal-monthly-single-levels"
                or entry.get("cds_model") != centre):
            continue
        seen[name] = str(entry.get("system"))
        if str(entry.get("system")) == cfg.systems[centre]:
            return name
    raise ConfigError(
        f"acmadDL has no monthly C3S product for centre {centre!r} system "
        f"{cfg.systems[centre]}. Catalog entries for this centre: {seen or 'none'}. "
        f"Change the system number under services, or set rebuild.products.{centre}."
    )


def _season_days(cfg: Config) -> int:
    return sum(pd.Period(f"2001-{m:02d}").days_in_month for m in cfg.target_months)


def _canonical(da: xr.DataArray, cfg: Config) -> xr.DataArray:
    if da.attrs.get("units") == "mm":
        # acmadDL returns model precipitation as a season total; observations come
        # as mm/day. Put both in mm/day (non-leap day count).
        da = (da / _season_days(cfg)).assign_attrs(units="mm/day")
    if "member" not in da.dims:
        da = da.expand_dims(member=[0])
    stray = [c for c in ("number", "member", "spatial_ref") if c in da.coords and c not in da.dims]
    return da.drop_vars(stray, errors="ignore").transpose("year", "member", "lat", "lon").sortby("lat")


def _load_acmaddl(cfg: Config) -> Inputs:
    import acmaddl

    var, season = cfg.variable, cfg.season
    y0, y1 = cfg.hindcast_years
    init_year = cfg.forecast_year if cfg.do_forecast else y1
    init = f"{init_year}-{cfg.month:02d}"
    wide = cfg.calibrate and cfg.predictor_region != cfg.region
    reuse_wide = wide and _contains(cfg.predictor_region, cfg.region)

    def fetch_pair(product, region):
        common = dict(init=init, target=season, region=list(region), seasonal="mean",
                      boundary="center", year_index=True, verbose=False,
                      grid_res=cfg.model_grid_res)
        hcst = _canonical(acmaddl.fetch(product, var, hindcast=(y0, y1), **common)[var], cfg)
        fcst = None
        if cfg.do_forecast:
            fcst = _canonical(acmaddl.fetch(product, var, **common)[var], cfg)
            fcst = fcst.sel(year=[init_year])
        return hcst, fcst

    predictors, targets, sources = {}, {}, {}
    for centre in cfg.centres:
        product = product_for(cfg, centre)
        sources[centre] = f"acmadDL {product}"
        print(f"[data] {centre}: fetching {product} {var} {season} init {init}")
        if not wide:
            pair = fetch_pair(product, cfg.region)
            predictors[centre] = targets[centre] = pair
            continue
        predictors[centre] = fetch_pair(product, cfg.predictor_region)
        if reuse_wide:
            targets[centre] = tuple(None if x is None else crop(x, cfg.region)
                                    for x in predictors[centre])
        else:
            targets[centre] = fetch_pair(product, cfg.region)

    off = cfg.season_year_offset
    print(f"[data] obs: fetching {cfg.obs} {var} {season}")
    # `months=` prunes the download to the season's months; `target=` alone fetches
    # all twelve (288 CHIRPS files instead of 72). Seasons that cross the calendar
    # year keep `target=`, which is what labels them by their starting year.
    months = cfg.target_months
    wraps = months != sorted(months)
    when = dict(target=season) if wraps else dict(months=months)
    if cfg.obs_request_interval:
        # Extra pacing for rate-limited servers (CHIRPS bans an IP for 4 hours above
        # ~2 requests/s). acmadDL enforces its own per-product floor regardless.
        when["request_interval"] = cfg.obs_request_interval
    # A season that crosses New Year needs the months after the last labelled year.
    obs = acmaddl.fetch(cfg.obs, var, region=cfg.region,
                        hindcast=(y0 + off, y1 + off + (1 if wraps else 0)), seasonal="mean",
                        grid_res=cfg.obs_grid_res, boundary="center", verbose=False,
                        **when)[var]
    if obs.attrs.get("units") == "mm/month":
        # e.g. CHIRPS: the season mean of monthly totals. Total over the season's days.
        obs = (obs * len(months) / _season_days(cfg)).assign_attrs(units="mm/day")
    obs = obs.assign_coords(year=obs["year"].values - off).transpose("year", "lat", "lon").sortby("lat")
    obs = obs.sel(year=slice(y0, y1))
    obs = obs.drop_vars([c for c in obs.coords if c not in obs.dims])  # e.g. a scalar `member`
    if obs.sizes["year"] < y1 - y0 + 1:
        print(f"[data] obs: only {obs.sizes['year']} of {y1 - y0 + 1} hindcast years available")
    sources["obs"] = f"acmadDL {cfg.obs}"
    return Inputs(predictors, targets, obs, sources)


# --------------------------------------------------------------------------- OSOP GRIB files

def _start_dim(path) -> str:
    """"indexing_time" for lagged ensembles (UKMO), "time" otherwise. As osop.util.get_tindex."""
    import eccodes

    with open(path, "rb") as handle:
        gid = eccodes.codes_grib_new_from_file(handle)
        try:
            eccodes.codes_get(gid, "indexingDate")
            return "indexing_time"
        except eccodes.KeyValueNotFoundError:
            return "time"
        finally:
            eccodes.codes_release(gid)


def _to_units(da: xr.DataArray, cfg: Config, *, model: bool) -> xr.DataArray:
    if cfg.variable == "precip":
        # model: mean rate in m/s. ERA5 monthly means: m per day.
        da = da * (86400.0 * 1000.0 if model else 1000.0)
        return da.assign_attrs(units="mm/day")
    return (da - 273.15).assign_attrs(units="C")


def read_model_grib(path, cfg: Config) -> xr.DataArray:
    """Season mean of an OSOP hindcast or forecast GRIB (unweighted mean of the lead months)."""
    start = _start_dim(path)
    ds = xr.open_dataset(path, engine="cfgrib", decode_timedelta=False,
                         backend_kwargs=dict(time_dims=("forecastMonth", start), indexpath=""))
    da = ds[list(ds.data_vars)[0]].load()
    da = da.sel(forecastMonth=cfg.leads).mean("forecastMonth")
    da = da.rename(latitude="lat", longitude="lon", number="member")
    if start in da.dims:
        years = pd.to_datetime(da[start].values).year
        da = da.assign_coords({start: years}).rename({start: "year"})
        da = da.transpose("year", "member", "lat", "lon")
    else:
        year = int(pd.to_datetime(da[start].values).year)
        da = da.expand_dims(year=[year]).transpose("year", "member", "lat", "lon")
    keep = set(da.dims)
    da = da.drop_vars([c for c in da.coords if c not in keep])
    return _to_units(da.sortby("lat"), cfg, model=True)


def read_obs_grib(path, cfg: Config) -> xr.DataArray:
    """Season mean of an OSOP ERA5 GRIB, labelled by initialisation year."""
    # decode_timedelta=False: recent xarray/pandas pairs trip an assertion decoding `step`.
    ds = xr.open_dataset(path, engine="cfgrib", decode_timedelta=False,
                         backend_kwargs=dict(indexpath=""))
    name = list(ds.data_vars)[0]
    da = ds[name].drop_vars("step", errors="ignore").load()
    valid = pd.to_datetime(da["valid_time"].values)
    da = da.assign_coords(time=valid.year * 12 + (valid.month - 1))
    da = da.rename(latitude="lat", longitude="lon")
    y0, y1 = cfg.hindcast_years
    seasons, years = [], []
    for year in range(y0, y1 + 1):
        wanted = [year * 12 + (cfg.month - 1) + (lead - 1) for lead in cfg.leads]
        if not all(w in da["time"].values for w in wanted):
            continue
        seasons.append(da.sel(time=wanted).mean("time"))
        years.append(year)
    if not seasons:
        raise ValueError(f"{path} holds none of the months needed for {cfg.season}")
    out = xr.concat(seasons, dim=pd.Index(years, name="year"))
    keep = set(out.dims)
    out = out.drop_vars([c for c in out.coords if c not in keep])
    return _to_units(out.transpose("year", "lat", "lon").sortby("lat"), cfg, model=False)


def _first_existing(*candidates) -> Path | None:
    for path in candidates:
        if path.is_file():
            return path
    return None


def _load_grib(cfg: Config) -> Inputs:
    if not cfg.grib_dir:
        raise ConfigError("rebuild.grib_dir is required when rebuild.backend is 'grib'")
    d = Path(cfg.grib_dir)
    y0, y1 = cfg.hindcast_years
    fy = cfg.forecast_year

    def name(prefix, area_str, centre, years):
        return d / (f"{prefix}{centre}_{cfg.systems[centre]}_{years}_monthly_mean_"
                    f"{cfg.month}_{cfg.leads_str}_{area_str}_{cfg.osop_variable}.grib")

    def pair(centre, primary, fallback):
        """(hindcast, forecast) from the first of two file families that exists."""
        out = []
        for years in (f"{y0}-{y1}", f"{fy}-{fy}"):
            is_forecast = years != f"{y0}-{y1}"
            if is_forecast and not cfg.do_forecast:
                out.append((None, None))
                continue
            path = _first_existing(name(*primary, centre, years), name(*fallback, centre, years))
            if path is None:
                raise FileNotFoundError(
                    f"No GRIB for {centre} {years} in {d}: looked for "
                    f"{name(*primary, centre, years).name} and {name(*fallback, centre, years).name}")
            out.append((read_model_grib(path, cfg), path))
        (hcst, hpath), (fcst, _) = out
        return hcst, fcst, hpath

    target_family = ("", cfg.area_str)
    predictor_family = ("pycpt_", cfg.predictor_area_str)
    predictors, targets, sources = {}, {}, {}
    for centre in cfg.centres:
        th, tf, tpath = pair(centre, target_family, predictor_family)
        targets[centre] = tuple(None if x is None else crop(x, cfg.region) for x in (th, tf))
        ph, pf, ppath = pair(centre, predictor_family, target_family)
        predictors[centre] = (ph, pf)
        sources[centre] = f"GRIB {ppath.name}" if ppath == tpath else f"GRIB {tpath.name} + {ppath.name}"
        print(f"[data] {centre}: {sources[centre]}")

    obs_var = cfg.osop_variable
    obs_path = d / f"era5_{obs_var}_{y0}-{y1}_monthly_{cfg.month}_{cfg.obs_leads_str}_{cfg.area_str}.grib"
    if not obs_path.is_file():
        raise FileNotFoundError(f"No observation GRIB: {obs_path}")
    obs = read_obs_grib(obs_path, cfg)
    sources["obs"] = f"GRIB {obs_path.name}"
    print(f"[data] obs: {sources['obs']}")
    return Inputs(predictors, targets, obs, sources)


def to_grid(da: xr.DataArray, like: xr.DataArray) -> xr.DataArray:
    """Bilinear interpolation of ``da`` onto the lat/lon of ``like`` (no-op if already there)."""
    if (da.sizes["lat"] == like.sizes["lat"] and da.sizes["lon"] == like.sizes["lon"]
            and np.allclose(da["lat"], like["lat"]) and np.allclose(da["lon"], like["lon"])):
        return da
    return da.interp(lat=like["lat"], lon=like["lon"], method="linear")
