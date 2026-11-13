"""pyCPT reference run on one comparison case. Runs inside the pycpt container.

Makes the same pyCPT calls, with the same settings, as OSOP's
``lib/osop/run_pycpt.py::calibrate`` (settings copied verbatim from there), and
saves what pyCPT returns as NetCDF so it can be compared with africas2s.

The case directory (see make_case.py) is mounted at /case; outputs go to /case/out.
"""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import pycpt
import xarray as xr

CASE = Path("/case")
OUT = CASE / "out"
meta = json.loads((CASE / "case.json").read_text())
names = meta["models"]

MOS = "CCA"
cpt_args = {  # lib/osop/run_pycpt.py, lines 232-243
    "transform_predictand": None,
    "tailoring": "Anomaly",
    "cca_modes": (1, 10),
    "x_eof_modes": (1, 10),
    "y_eof_modes": (1, 10),
    "validation": "crossvalidation",
    "drymask_threshold": None,
    "skillmask_threshold": None,
    "crossvalidation_window": 5,
    "synchronous_predictors": True,
}


def first_var(path):
    ds = xr.open_dataset(path)
    return list(ds.data_vars.values())[0]


def save(obj, name):
    if isinstance(obj, (list, tuple)):
        for model, item in zip(names, obj):
            save(item, f"{name}_{model}")
        return
    try:
        obj.to_netcdf(OUT / f"{name}.nc")
        print(f"saved {name}.nc", flush=True)
    except Exception as exc:  # keep going: partial output is still useful
        print(f"could not save {name}: {type(exc).__name__}: {exc}", flush=True)


Y = first_var(CASE / "inputs" / "obs.nc")
hind = [first_var(CASE / "inputs" / f"{m}_hcst.nc") for m in names]
fore = [first_var(CASE / "inputs" / f"{m}_fcst.nc") for m in names]
print("pycpt", pycpt.__version__, "| Y", dict(Y.sizes), "| X", [dict(h.sizes) for h in hind], flush=True)

OUT.mkdir(parents=True, exist_ok=True)
domain_dir = pycpt.setup(OUT / "case", meta["predictor_extent"])
hcsts, fcsts, skill, pxs, pys = pycpt.evaluate_models(
    hind, MOS, Y, fore, cpt_args, domain_dir, names, False)
for label, obj in (("hcst", hcsts), ("fcst", fcsts), ("skill", skill), ("px", pxs), ("py", pys)):
    save(obj, label)

if len(names) > 1:
    det_fcst, pr_fcst, pev_fcst, nextgen_skill = pycpt.construct_mme(
        fcsts, hcsts, Y, names, names, cpt_args, domain_dir)
    for label, obj in (("mme_det_fcst", det_fcst), ("mme_pr_fcst", pr_fcst),
                       ("mme_pev_fcst", pev_fcst), ("mme_skill", nextgen_skill)):
        save(obj, label)
print("done", flush=True)
