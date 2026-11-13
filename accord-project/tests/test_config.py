"""Config loading: OSOP's own YAML files must load unchanged."""

import pytest

from conftest import OSOP_ROOT, REBUILD
from osop_acmad.config import ConfigError, load_config
from osop_acmad.data import product_for


def test_osop_stock_config_loads():
    cfg = load_config(OSOP_ROOT / "osop_config.yml")
    assert cfg.season == "DJF"            # November start, leads 2-4
    assert cfg.season_year_offset == 0    # December is in the initialisation year
    assert cfg.region == [21.0, 34.0, 22.0, 38.0]
    assert cfg.centres == ["dwd", "ukmo", "ecmwf"] and cfg.mme
    assert sum(cfg.weights.values()) == pytest.approx(1.0)
    assert cfg.backend == "acmaddl"


def test_season_in_next_calendar_year(tmp_path):
    text = (OSOP_ROOT / "osop_config.yml").read_text().replace("month: 11", "month: 12")
    path = tmp_path / "dec.yml"
    path.write_text(text)
    cfg = load_config(path)
    assert cfg.season == "JFM" and cfg.season_year_offset == 1


def test_c3s_centres_map_to_acmaddl_products():
    cfg = load_config(OSOP_ROOT / "osop_config.yml")
    # ukmo is left out: OSOP's config names system 604 and acmadDL carries 610.
    assert {c: product_for(cfg, c) for c in ("dwd", "ecmwf")} == {
        "dwd": "c3s/dwd", "ecmwf": "c3s/ecmwf-monthly"}


def test_unknown_system_fails_loudly():
    cfg = load_config(OSOP_ROOT / "osop_config.yml")
    cfg.systems["ukmo"] = "9999"
    with pytest.raises(ConfigError, match="system 9999"):
        product_for(cfg, "ukmo")


def test_offline_config_uses_grib_backend():
    cfg = load_config(REBUILD / "configs" / "offline_test.yml")
    assert cfg.backend == "grib" and cfg.season == "JJA"
