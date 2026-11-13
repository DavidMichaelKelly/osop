"""Shared fixtures: the offline config and OSOP's own library on the path."""

import sys
from pathlib import Path

import pytest

REBUILD = Path(__file__).resolve().parents[1]
OSOP_ROOT = REBUILD.parent
TEST_DATA = OSOP_ROOT / "lib" / "osop" / "test" / "test_data" / "test_pycptf"

sys.path.insert(0, str(REBUILD))
sys.path.insert(0, str(OSOP_ROOT / "lib"))


@pytest.fixture(scope="session")
def cfg():
    from osop_acmad.config import load_config

    return load_config(REBUILD / "configs" / "offline_test.yml")


@pytest.fixture(scope="session")
def inputs(cfg):
    from osop_acmad.data import load_inputs

    return load_inputs(cfg)
