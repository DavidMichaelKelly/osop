"""Run the OSOP workflow on acmadDL + africas2s.

Products, per centre and for the weighted multi-model ensemble (MME):

* ``member_count``: OSOP's uncalibrated product (members per hindcast tercile).
* ``cca``: africas2s CCA with OSOP's pyCPT settings. Replaces the pyCPT step.
* ``ereg``: africas2s ensemble regression, for comparison.

Every product has a forecast map and a cross-validated hindcast scored against
the observations, so the methods can be compared like for like.
"""

from __future__ import annotations

import json
import time
import warnings
from pathlib import Path

import numpy as np
import xarray as xr

from . import products
from .config import Config
from .data import Inputs, load_inputs, to_grid

# OSOP's fixed pyCPT settings (lib/osop/run_pycpt.py), in africas2s spelling.
OSOP_CCA_ARGS = dict(
    mode_selection="cpt",          # cross-validated search, as CPT does
    mode_selection_tau="cpt",      # rank by CPT's own Kendall tau (ties as in CPT)
    x_eof_range=(1, 10),
    y_eof_range=(1, 10),
    cca_range=(1, 10),
    crossvalidation_window=5,
    # africas2s's opt-in CPT-parity behaviour. Against a real pyCPT run on OSOP's
    # test data these raise the agreement of the cross-validated hindcasts from
    # about 0.3 to 0.73 (see pycpt_ref/).
    standardize=True,
    lat_weights="cpt",
    sign_convention="cpt",
    y_reconstruction="cpt",
    # OSOP also sets tailoring="Anomaly". Left out on purpose: in africas2s 0.3.0
    # it turns the deterministic forecast into an anomaly that is then compared
    # with absolute tercile boundaries, which puts nearly every cell in "below".
    # In CPT the setting only changes how the deterministic map is presented.
)
SKILL_METRICS = ["rpss", "generalized_roc"]


def _ordered(probs: xr.DataArray) -> xr.DataArray:
    lead = ("year", "tercile") if "year" in probs.dims else ("tercile",)
    return probs.transpose(*lead, "lat", "lon")


def combine(per_centre: dict, weights: dict) -> xr.DataArray:
    """Weighted mean of tercile probabilities; a centre with no value at a cell abstains."""
    num, den = None, None
    for centre, probs in per_centre.items():
        present = probs.notnull()
        term, wgt = probs.fillna(0.0) * weights[centre], present * weights[centre]
        num = term if num is None else num + term
        den = wgt if den is None else den + wgt
    out = (num / den.where(den > 0))
    return _ordered(out / out.sum("tercile", skipna=False))


def _member_count(cfg: Config, inp: Inputs, grid: xr.DataArray):
    forecast, cv, native = {}, {}, {}
    for centre, (hcst, fcst) in inp.targets.items():
        thr = products.tercile_thresholds(hcst)
        hind = products.hindcast_probabilities(hcst, thr)
        native[centre] = {"thresholds": thr, "hindcast": hind}
        cv[centre] = _ordered(to_grid(hind, grid))
        if fcst is not None:
            fc = products.forecast_probabilities(fcst, thr)
            native[centre]["forecast"] = fc
            forecast[centre] = _ordered(to_grid(fc, grid))
    return forecast, cv, native


def _cca(cfg: Config, inp: Inputs):
    import africas2s

    args = {**OSOP_CCA_ARGS, **cfg.cca_args}
    forecast, cv, info = {}, {}, {}
    for centre, (hcst, fcst) in inp.predictors.items():
        t0 = time.time()
        res = africas2s.seasonal_mme(
            {cfg.variable: {centre: (hcst, fcst)}}, inp.obs, method="cca",
            cpt_args=args, probability_aggregation=cfg.probability_aggregation,
            cv="loyo_cyclic",  # CPT's wraparound leave-out window
            forecast_year=cfg.forecast_year if fcst is not None else None,
            verbose=False)  # no skill_metrics: an explicit list raises in africas2s 0.3.0
        cv[centre] = _ordered(res.tercile_cv)
        if fcst is not None:
            forecast[centre] = _ordered(res.tercile_forecast)
        modes = {k: v for k, v in res.metadata.items() if "mode" in k}
        info[centre] = {"seconds": round(time.time() - t0, 1), "modes": modes}
        print(f"[cca] {centre}: {info[centre]['seconds']} s {modes}")
    return forecast, cv, info


def _ereg(cfg: Config, inp: Inputs):
    import africas2s

    obs = inp.obs
    on_obs = {c: (to_grid(h, obs), None if f is None else to_grid(f, obs))
              for c, (h, f) in inp.targets.items()}
    kw = dict(method="ereg", return_components=True,
              clip_negative=cfg.variable == "precip")
    cv_res = africas2s.calibrate({c: (h, None) for c, (h, _) in on_obs.items()}, obs,
                                 cv="loyo", cv_window=5, **kw)
    cv = {c: _ordered(p) for c, p in cv_res.per_model.items()}
    forecast = {}
    if cfg.do_forecast:
        fc_res = africas2s.calibrate(on_obs, obs, forecast_year=cfg.forecast_year, **kw)
        forecast = {c: _ordered(p) for c, p in fc_res.per_model.items()}
    return forecast, cv


def _score(cv: xr.DataArray, obs: xr.DataArray):
    import africas2s

    cv = to_grid(cv, obs)
    years = np.intersect1d(cv["year"].values, obs["year"].values)
    report = africas2s.skill(cv.sel(year=years), obs.sel(year=years),
                             metrics=SKILL_METRICS, spatial=True)
    return report


def run(cfg: Config, inp: Inputs | None = None) -> dict:
    """Run every configured product. Returns the in-memory results; also writes files."""
    warnings.filterwarnings("ignore", category=FutureWarning)
    inp = inp or load_inputs(cfg)
    obs = inp.obs
    first = next(iter(inp.targets.values()))[0]

    forecasts: dict[str, dict] = {}   # method -> {centre or "mme": (tercile, lat, lon)}
    hindcasts: dict[str, dict] = {}   # method -> {centre or "mme": (year, tercile, lat, lon)}
    extra: dict = {"sources": inp.sources}

    if "member_count" in cfg.methods:
        forecasts["member_count"], hindcasts["member_count"], native = _member_count(cfg, inp, first)
        extra["member_count_native"] = native
    if cfg.calibrate and "cca" in cfg.methods:
        forecasts["cca"], hindcasts["cca"], extra["cca"] = _cca(cfg, inp)
    if cfg.calibrate and "ereg" in cfg.methods:
        forecasts["ereg"], hindcasts["ereg"] = _ereg(cfg, inp)

    if cfg.mme and len(cfg.centres) > 1:
        for store in (forecasts, hindcasts):
            for method, per_centre in store.items():
                if per_centre:
                    ref = next(iter(per_centre.values()))
                    same = {c: to_grid(p, ref) for c, p in per_centre.items()}
                    per_centre["mme"] = combine(same, cfg.weights)

    # --- dry mask: where the observed season is close to rainless, terciles are
    # degenerate (both boundaries near zero) and regression-based methods report
    # "above normal" every year. Blank those cells in every product and score.
    if cfg.variable == "precip" and cfg.dry_mask is not None:
        clim = obs.mean("year")
        extra["dry_cells"] = round(float((clim < cfg.dry_mask).mean()), 3)
        for store in (forecasts, hindcasts):
            for per_centre in store.values():
                for centre, probs in per_centre.items():
                    wet = to_grid(clim, probs) >= cfg.dry_mask
                    per_centre[centre] = probs.where(wet)

    rows, skill_maps = [], {}
    if cfg.do_hindcast:
        for method, per_centre in hindcasts.items():
            for centre, cv in per_centre.items():
                report = _score(cv, obs)
                rows.append({"method": method, "centre": centre,
                             **{m: round(float(report.scores[m]), 4) for m in SKILL_METRICS}})
                skill_maps[(method, centre)] = report.spatial
    _write(cfg, forecasts, hindcasts, rows, skill_maps, extra)
    return {"forecasts": forecasts, "hindcasts": hindcasts, "skill": rows,
            "skill_maps": skill_maps, "inputs": inp, "extra": extra}


def _stem(cfg: Config, method: str, centre: str) -> str:
    system = cfg.systems.get(centre, "mme")
    return (f"{centre}_{system}_{cfg.month}_{cfg.leads_str}_{cfg.area_str}_"
            f"{cfg.osop_variable}.{method}")


def _write(cfg, forecasts, hindcasts, rows, skill_maps, extra):
    for method, per_centre in hindcasts.items():
        for centre, cv in per_centre.items():
            cv.rename("probability").to_netcdf(
                cfg.out("hindcast", "products") / f"{_stem(cfg, method, centre)}.hindcast_probs.nc")
    for method, per_centre in forecasts.items():
        for centre, fc in per_centre.items():
            fc.rename("probability").to_netcdf(
                cfg.out("forecast", "products")
                / f"{_stem(cfg, method, centre)}.{cfg.forecast_year}.forecast_probs.nc")
    for centre, native in extra.get("member_count_native", {}).items():
        native["thresholds"].rename("threshold").to_netcdf(
            cfg.out("hindcast", "products")
            / f"{_stem(cfg, 'member_count', centre)}.tercile_thresholds.nc")

    scores_dir = cfg.out("hindcast", "scores")
    if rows:
        lines = ["method,centre," + ",".join(SKILL_METRICS)]
        lines += [f"{r['method']},{r['centre']}," + ",".join(str(r[m]) for m in SKILL_METRICS)
                  for r in rows]
        (scores_dir / "skill_table.csv").write_text("\n".join(lines) + "\n")
        for (method, centre), maps in skill_maps.items():
            xr.Dataset({m: maps[m] for m in SKILL_METRICS if m in maps}).to_netcdf(
                scores_dir / f"{_stem(cfg, method, centre)}.skill.nc")

    summary = {
        "season": cfg.season, "init_month": cfg.month, "forecast_year": cfg.forecast_year,
        "hindcast_years": list(cfg.hindcast_years), "area_nwse": cfg.area,
        "centres": cfg.centres, "weights": cfg.weights, "sources": extra["sources"],
        "cca": extra.get("cca"), "dry_mask_mm_day": cfg.dry_mask,
        "dry_cell_share": extra.get("dry_cells"), "skill": rows,
    }
    Path(cfg.paths["base"]).mkdir(parents=True, exist_ok=True)
    (Path(cfg.paths["base"]) / "run_summary.json").write_text(json.dumps(summary, indent=2, default=str))

    try:
        _plots(cfg, forecasts, skill_maps)
    except Exception as exc:  # a plotting failure must not lose the numbers above
        print(f"[plots] skipped: {type(exc).__name__}: {exc}")


def _plots(cfg, forecasts, skill_maps):
    import africas2s
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n, w, s, e = cfg.area
    style_name = "icpac" if cfg.variable == "precip" else "icpac-temperature"
    style = africas2s.TercileStyle.named(style_name, extent=(w, e, s, n))
    order = [*cfg.centres, "mme"]
    labels = {"member_count": "member count", "cca": "CCA", "ereg": "eReg"}

    panels = {f"{c.upper()} {labels[m]}": forecasts[m][c]
              for m in forecasts for c in order if c in forecasts[m]}
    if panels:
        ncols = len([c for c in order if c in next(iter(forecasts.values()))])
        fig = africas2s.plot_matrix(
            panels, style=style, ncols=ncols, legend_detailed=True,
            variable_kind=cfg.variable,
            suptitle=f"{cfg.season} {cfg.forecast_year} tercile outlook, init month {cfg.month}")
        fig.savefig(cfg.out("forecast", "plots") / "forecast_by_method.png", dpi=130,
                    bbox_inches="tight")
        plt.close(fig)

    for metric, (vmin, vmax, center, cmap) in {"rpss": (-0.3, 0.3, 0.0, "RdBu"),
                                               "generalized_roc": (0.2, 0.8, 0.5, "RdBu")}.items():
        maps = {f"{c.upper()} {labels[m]}": sm[metric]
                for (m, c), sm in skill_maps.items() if metric in sm}
        if not maps:
            continue
        ncols = len({c for (_, c) in skill_maps})
        fig = africas2s.plot_matrix(
            maps, style=style, ncols=ncols, cmap=cmap, vmin=vmin, vmax=vmax,
            cbar_label=metric, suptitle=f"{cfg.season} hindcast skill: {metric}")
        fig.savefig(cfg.out("hindcast", "plots") / f"skill_{metric}.png", dpi=130,
                    bbox_inches="tight")
        plt.close(fig)
