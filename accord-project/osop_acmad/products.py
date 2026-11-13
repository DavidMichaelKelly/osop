"""OSOP's uncalibrated "member count" tercile product, ported to plain xarray.

Mirrors ``osop.compute_products_func.prob_terc`` (hindcast) and
``osop.compare_terciles.nth_month`` (forecast), including their inequality
conventions, so the two can be compared number for number.

Inputs are season means shaped ``(year, member, lat, lon)``; a forecast holds
a single year. Outputs carry ``tercile = [0, 1, 2]`` for
below / normal / above, as fractions.
"""

from __future__ import annotations

import xarray as xr

TERCILE = [0, 1, 2]


def tercile_thresholds(hindcast: xr.DataArray) -> xr.DataArray:
    """Lower and upper tercile boundaries over all hindcast years and members.

    In-sample, as in OSOP: the year being scored is part of its own climatology.
    """
    q = hindcast.quantile([1 / 3.0, 2 / 3.0], dim=["year", "member"])
    return q.rename(quantile="bound").assign_coords(bound=["lower", "upper"])


def _stack(below, normal, above, members: int) -> xr.DataArray:
    probs = xr.concat([below, normal, above], dim="tercile") / float(members)
    return probs.assign_coords(tercile=TERCILE)


def hindcast_probabilities(hindcast: xr.DataArray, thresholds: xr.DataArray) -> xr.DataArray:
    """Fraction of members in each tercile, per hindcast year: ``(tercile, year, lat, lon)``."""
    lo = thresholds.sel(bound="lower", drop=True)
    hi = thresholds.sel(bound="upper", drop=True)
    below = (hindcast <= lo).sum("member")
    normal = ((hindcast > lo) & (hindcast <= hi)).sum("member")
    above = (hindcast > hi).sum("member")
    return _stack(below, normal, above, hindcast.sizes["member"])


def forecast_probabilities(forecast: xr.DataArray, thresholds: xr.DataArray) -> xr.DataArray:
    """Fraction of forecast members in each hindcast tercile: ``(tercile, lat, lon)``.

    OSOP uses strict inequalities here, so a member exactly on a boundary is
    counted in no category. That cannot happen with real float data.
    """
    if "year" in forecast.dims:
        forecast = forecast.squeeze("year", drop=True)
    lo = thresholds.sel(bound="lower", drop=True)
    hi = thresholds.sel(bound="upper", drop=True)
    below = (forecast < lo).sum("member")
    normal = ((forecast > lo) & (forecast < hi)).sum("member")
    above = (forecast > hi).sum("member")
    return _stack(below, normal, above, forecast.sizes["member"])


def weighted_mme(per_centre: dict[str, xr.DataArray], weights: dict[str, float]) -> xr.DataArray:
    """OSOP's multi-model product: weighted mean of the per-centre fields.

    ``weights`` must sum to 1 over the centres given.
    """
    total = sum(weights[c] for c in per_centre)
    out = None
    for centre, field in per_centre.items():
        term = field * (weights[centre] / total)
        out = term if out is None else out + term
    return out
