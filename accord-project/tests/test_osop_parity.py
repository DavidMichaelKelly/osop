"""The ported member-count product must equal what OSOP's own code produces.

Runs OSOP's library functions (not copies of them) on the GRIB files OSOP ships
as test data, then compares with ``osop_acmad.products`` on the same files.
"""

import shutil

import numpy as np
import pytest
import xarray as xr

from conftest import TEST_DATA
from osop_acmad import products

AREA = "40.10.10.40"
VAR = "total_precipitation"
CENTRES = {"ecmwf": 51, "ukmo": 604}


def _bname(centre, system, years):
    return f"{centre}_{system}_{years}_monthly_mean_5_234_{AREA}_{VAR}"


@pytest.fixture(scope="module")
def osop_run(tmp_path_factory):
    """OSOP's hindcast products, forecast products and MME, written to a temp dir."""
    from osop.compare_terciles import compute_forecast, mme_products
    from osop.compute_products_func import calc_products, calc_products_mme

    root = tmp_path_factory.mktemp("osop_ref")
    dl, hc, fc = root / "downloads", root / "hc_products", root / "fc_products"
    for d in (dl, hc, fc):
        d.mkdir()
    # OSOP skips a download when the file exists, so naming the shipped files the
    # way its scripts expect is enough to run it offline.
    for centre, system in CENTRES.items():
        for years in ("1993-2016", "2025-2025"):
            shutil.copy(TEST_DATA / f"pycpt_{_bname(centre, system, years)}.grib",
                        dl / f"{_bname(centre, system, years)}.grib")

    base = dict(start_month=5, leads=[2, 3, 4], leads_str="234", area_str=AREA,
                var=VAR, hc_var=VAR, hcstarty=1993, hcendy=2016, fcstarty=2025, fcendy=2025)
    for centre, system in CENTRES.items():
        config = dict(base, origin=centre, system=system, systemhc=system, systemfc=system)
        calc_products(config, str(dl), str(hc))
        compute_forecast(config, str(dl), str(hc), str(fc))

    def services():  # fresh each call: OSOP mutates the dict
        return {"ecmwf": [51, 0.5], "ukmo": [604, 0.5], "mme": [1, 0]}

    mme_cfg = dict(base, origin="mme", systemfc=1, systemhc=1)
    calc_products_mme(services(), mme_cfg, str(hc))
    mme_products(services(), mme_cfg, str(fc))
    return hc, fc


def _osop_hindcast(hc, centre, system):
    ds = xr.open_dataset(hc / f"{_bname(centre, system, '1993-2016')}.nm.tercile_probs.nc")
    da = ds[list(ds.data_vars)[0]].squeeze("forecastMonth", drop=True)
    return da.transpose("category", "start_date", "lat", "lon").sortby("lat").values


def _osop_forecast(fc, centre, system):
    ds = xr.open_dataset(fc / f"{_bname(centre, system, '2025-2025')}.nm.forecast_percentages.nc")
    ds = ds.squeeze(drop=True).sortby("lat")
    return np.stack([ds[k].transpose("lat", "lon").values for k in ("lower", "middle", "higher")]) / 100.0


@pytest.mark.parametrize("centre", list(CENTRES))
def test_hindcast_probabilities_match_osop(osop_run, inputs, centre):
    hc, _ = osop_run
    hcst, _ = inputs.targets[centre]
    mine = products.hindcast_probabilities(hcst, products.tercile_thresholds(hcst))
    mine = mine.transpose("tercile", "year", "lat", "lon").values
    np.testing.assert_allclose(mine, _osop_hindcast(hc, centre, CENTRES[centre]), atol=1e-12)


def _assert_forecast_parity(mine, ref, members):
    """Equal, except that a member within float32 rounding of a boundary may flip.

    Seen only in desert cells with ~1e-5 mm/day, where GRIB-packed member values
    sit on the hindcast tercile boundary. OSOP thresholds in m/s, this package in
    mm/day, and the unit conversion decides which side such a member lands on.
    """
    diff = np.abs(mine - ref)
    assert diff.max() <= 1.0 / members + 1e-12
    assert (diff > 1e-12).mean() < 0.005


@pytest.mark.parametrize("centre", list(CENTRES))
def test_forecast_probabilities_match_osop(osop_run, inputs, centre):
    _, fc = osop_run
    hcst, fcst = inputs.targets[centre]
    mine = products.forecast_probabilities(fcst, products.tercile_thresholds(hcst))
    _assert_forecast_parity(mine.transpose("tercile", "lat", "lon").values,
                            _osop_forecast(fc, centre, CENTRES[centre]), fcst.sizes["member"])


def test_mme_matches_osop(osop_run, inputs, cfg):
    hc, fc = osop_run
    hind, fore = {}, {}
    for centre, (hcst, fcst) in inputs.targets.items():
        thr = products.tercile_thresholds(hcst)
        hind[centre] = products.hindcast_probabilities(hcst, thr)
        fore[centre] = products.forecast_probabilities(fcst, thr)
    mine_h = products.weighted_mme(hind, cfg.weights).transpose("tercile", "year", "lat", "lon")
    mine_f = products.weighted_mme(fore, cfg.weights).transpose("tercile", "lat", "lon")
    np.testing.assert_allclose(mine_h.values, _osop_hindcast(hc, "mme", 1), atol=1e-12)
    fewest = min(f.sizes["member"] for _, f in inputs.targets.values())
    _assert_forecast_parity(mine_f.values, _osop_forecast(fc, "mme", 1), fewest)
