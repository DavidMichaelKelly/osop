"""Command line: ``python -m osop_acmad --config path/to/osop_config.yml``."""

from __future__ import annotations

import argparse
import sys

from .config import ConfigError, load_config
from .pipeline import SKILL_METRICS, run


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="osop_acmad", description="Run the OSOP workflow on acmadDL + africas2s")
    parser.add_argument("--config", required=True, help="OSOP YAML configuration file")
    args = parser.parse_args(argv)
    try:
        cfg = load_config(args.config)
        result = run(cfg)
    except ConfigError as exc:
        print(f"Configuration error: {exc}", file=sys.stderr)
        return 2

    if result["skill"]:
        print(f"\n{cfg.season} hindcast skill vs observations, {cfg.hindcast_years[0]}-{cfg.hindcast_years[1]}")
        print(f"{'method':<14}{'centre':<10}" + "".join(f"{m:>18}" for m in SKILL_METRICS))
        for row in result["skill"]:
            print(f"{row['method']:<14}{row['centre']:<10}"
                  + "".join(f"{row[m]:>18.3f}" for m in SKILL_METRICS))
    print(f"\nOutputs under {cfg.paths['base']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
