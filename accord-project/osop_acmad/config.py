"""Load an OSOP YAML config, plus an optional ``rebuild:`` section.

The required sections (workflow, parameters, paths, centres, services) are
OSOP's own, so one file drives both OSOP's ``run_all.py`` and this package.
OSOP ignores the extra ``rebuild:`` section.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

MONTH_INITIALS = "JFMAMJJASOND"
VARIABLES = {"total_precipitation": "precip", "2m_temperature": "temp"}
HINDCAST_YEARS = (1993, 2016)  # OSOP's hindcast period: the C3S common period.


class ConfigError(Exception):
    """The configuration file is missing something or holds a bad value."""


@dataclass
class Config:
    month: int                     # forecast initialisation month, 1-12
    leads: list[int]               # OSOP leads: lead 1 is the initialisation month
    area: list[float]              # target box [north, west, south, east]
    predictor_area: list[float]    # predictor box for CCA, same order
    area_str: str                  # OSOP file-name spelling of the boxes
    predictor_area_str: str
    osop_variable: str             # "total_precipitation" or "2m_temperature"
    location: str
    forecast_year: int
    do_hindcast: bool
    do_forecast: bool
    calibrate: bool                # OSOP's ``pycpt`` switch
    centres: list[str]             # without "mme"
    mme: bool
    systems: dict[str, str]        # centre -> CDS system number
    weights: dict[str, float]      # centre -> MME weight, normalised to sum 1
    paths: dict
    hindcast_years: tuple[int, int] = HINDCAST_YEARS
    # --- rebuild-only options ---
    backend: str = "acmaddl"       # "acmaddl" (live fetch) or "grib" (OSOP downloads on disk)
    grib_dir: str | None = None
    obs: str = "obs/era5"          # acmadDL observation product
    obs_grid_res: float | None = 1.0
    obs_request_interval: float | None = None  # seconds between obs file opens; None = acmadDL's floor
    model_grid_res: float | None = None
    products: dict[str, str] = field(default_factory=dict)  # centre -> acmadDL product override
    methods: list[str] = field(default_factory=lambda: ["member_count", "cca", "ereg"])
    cca_args: dict = field(default_factory=dict)
    probability_aggregation: str = "cpt_per_model"
    dry_mask: float | None = 0.1   # mm/day; cells with a drier observed season are blanked
    config_path: Path | None = None

    @property
    def variable(self) -> str:
        """acmadDL / africas2s variable name."""
        return VARIABLES[self.osop_variable]

    @property
    def target_months(self) -> list[int]:
        return [(self.month - 1 + lead - 1) % 12 + 1 for lead in self.leads]

    @property
    def season(self) -> str:
        """Season code of the target months, e.g. "JJA"."""
        return "".join(MONTH_INITIALS[m - 1] for m in self.target_months)

    @property
    def season_year_offset(self) -> int:
        """Calendar years between the initialisation and the season's first month."""
        return (self.month - 1 + self.leads[0] - 1) // 12

    @property
    def region(self) -> list[float]:
        """Target box in acmadDL order [south, north, west, east]."""
        n, w, s, e = self.area
        return [s, n, w, e]

    @property
    def predictor_region(self) -> list[float]:
        n, w, s, e = self.predictor_area
        return [s, n, w, e]

    @property
    def leads_str(self) -> str:
        return "".join(str(lead) for lead in self.leads)

    @property
    def obs_leads_str(self) -> str:
        return "".join(str(lead - 1) for lead in self.leads)

    def out(self, phase: str, kind: str) -> Path:
        path = Path(self.paths[phase][kind])
        path.mkdir(parents=True, exist_ok=True)
        return path


def _bbox(text, name) -> list[float]:
    parts = [p.strip() for p in str(text).split(",") if p.strip()]
    if len(parts) != 4:
        raise ConfigError(f"{name} must hold 4 comma-separated numbers, got: {text}")
    return [float(p) for p in parts]


def load_config(path) -> Config:
    path = Path(path).expanduser().resolve()
    if not path.exists():
        raise ConfigError(f"Config file not found: {path}")
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    for section in ("workflow", "parameters", "paths", "centres", "services"):
        if section not in raw:
            raise ConfigError(f"Missing required section: {section}")

    params, workflow, extra = raw["parameters"], raw["workflow"], raw.get("rebuild") or {}
    variable = str(params["variable"])
    if variable not in VARIABLES:
        raise ConfigError(f"Unsupported variable: {variable}. Allowed: {sorted(VARIABLES)}")
    month = int(params["month"])
    if not 1 <= month <= 12:
        raise ConfigError(f"Month must be 1-12, got: {month}")
    leads = [int(x) for x in str(params["leads"]).split(",") if x.strip()]
    if not leads or min(leads) <= 0:
        raise ConfigError(f"Leads must be positive integers, got: {params['leads']}")
    if leads != list(range(leads[0], leads[0] + len(leads))):
        raise ConfigError(f"Leads must be consecutive months, got: {params['leads']}")

    area = _bbox(params["area"], "parameters.area")
    predictor_area = (_bbox(params["predictor_area"], "parameters.predictor_area")
                      if params.get("predictor_area") else area)

    selected = [str(c) for c in raw["centres"]]
    centres = [c for c in selected if c != "mme"]
    if not centres:
        raise ConfigError("centres must name at least one modelling centre")
    services = raw["services"]
    systems, weights = {}, {}
    for centre in centres:
        if centre not in services:
            raise ConfigError(f"centre {centre!r} has no entry under services")
        entry = services[centre]
        entry = list(entry) if isinstance(entry, (list, tuple)) else [entry]
        systems[centre] = str(entry[0])
        weights[centre] = float(entry[1]) if len(entry) > 1 else 1.0
    total = sum(weights.values())
    if total <= 0:
        raise ConfigError("MME weights must sum to a positive number")
    weights = {c: w / total for c, w in weights.items()}

    base = str(raw["paths"]["base"])
    if not Path(base).is_absolute():
        base = str(Path.cwd() / base)
    paths = {"base": base}
    for key, value in raw["paths"].items():
        if key == "base":
            continue
        paths[key] = ({k: str(v).format(base=base) for k, v in value.items()}
                      if isinstance(value, dict) else str(value).format(base=base))

    backend = str(extra.get("backend", "acmaddl"))
    if backend not in ("acmaddl", "grib"):
        raise ConfigError(f"rebuild.backend must be 'acmaddl' or 'grib', got: {backend}")
    grib_dir = extra.get("grib_dir")
    if grib_dir and not Path(grib_dir).is_absolute():
        grib_dir = str((path.parent / grib_dir).resolve())
    hindcast_years = tuple(extra.get("hindcast_years", HINDCAST_YEARS))

    do_forecast = bool(workflow.get("forecast", False))
    if do_forecast and "forecast_year" not in params:
        raise ConfigError("parameters.forecast_year is required when forecast is enabled")

    return Config(
        month=month, leads=leads, area=area, predictor_area=predictor_area,
        area_str=str(params["area"]).replace(",", "."),
        predictor_area_str=str(params.get("predictor_area") or params["area"]).replace(",", "."),
        osop_variable=variable, location=str(params.get("location", "None")),
        forecast_year=int(params.get("forecast_year", 0)),
        do_hindcast=bool(workflow.get("hindcast", True)), do_forecast=do_forecast,
        calibrate=str(params.get("pycpt", True)).strip().lower() in ("1", "true", "yes", "y", "on"),
        centres=centres, mme="mme" in selected, systems=systems, weights=weights, paths=paths,
        hindcast_years=(int(hindcast_years[0]), int(hindcast_years[1])),
        backend=backend, grib_dir=grib_dir,
        obs=str(extra.get("obs", "obs/era5")),
        obs_grid_res=extra.get("obs_grid_res", 1.0),
        obs_request_interval=extra.get("obs_request_interval"),
        model_grid_res=extra.get("model_grid_res"),
        products=dict(extra.get("products") or {}),
        methods=list(extra.get("methods", ["member_count", "cca", "ereg"])),
        cca_args=dict(extra.get("cca_args") or {}),
        probability_aggregation=str(extra.get("probability_aggregation", "cpt_per_model")),
        dry_mask=extra.get("dry_mask", 0.1),
        config_path=path,
    )
