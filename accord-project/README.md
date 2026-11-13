# osop_acmad: the OSOP workflow on acmadDL + africas2s

An experiment: keep OSOP's YAML-driven workflow, replace its engine.
acmadDL supplies the data, africas2s the statistics and maps. No pyCPT,
cptcore, CPT Fortran binary or IRI conda channel.

Status: working prototype. Runs OSOP's stock settings against live C3S data
(with UKMO moved from system 604 to 610: see
[Data-source caution](#data-source-caution)) and
matches OSOP's own code for the uncalibrated product. The CCA reproduces real
pyCPT to about 1% when given the same mode counts, but africas2s's automatic
mode search does not pick the same counts as CPT: see
[Comparison with pyCPT](#comparison-with-pycpt). See [What is and is not verified](#what-is-and-is-not-verified).

## Run it

```bash
cd accord-project
uv venv --python 3.13 .venv
uv pip install --python .venv/bin/python \
    "acmadDL>=0.3.2" "africas2s>=0.4.1" \
    pyyaml xskillscore xarray-regrid pytest

# Offline, on the GRIB files OSOP ships as test data (about 2 minutes):
.venv/bin/python -m osop_acmad --config configs/offline_test.yml

# Live through acmadDL, no credentials (NMME + CMAP):
.venv/bin/python -m osop_acmad --config configs/nmme_live_test.yml

# OSOP's stock config with UKMO on system 610 and the forecast switched on
# (needs ~/.cdsapirc; about 20 minutes the first time, half of it CCA):
.venv/bin/python -m osop_acmad --config configs/egypt_c3s.yml

# OSOP's own ../osop_config.yml needs the same UKMO change (604 to 610) to run.

.venv/bin/python -m pytest tests -q
```

Outputs land in OSOP's directory layout under `paths.base`:
`forecast/products`, `forecast/plots`, `hindcast/products`, `hindcast/scores`
(including `skill_table.csv`), `hindcast/plots`, plus `run_summary.json`.

## What it does

The config is OSOP's (`workflow`, `parameters`, `paths`, `centres`,
`services`). An optional `rebuild:` section, which OSOP ignores, adds the
options below.

For each centre and for the weighted multi-model ensemble it produces a
forecast tercile map and a cross-validated hindcast for three methods:

| Product | What it is | Replaces |
|---|---|---|
| `member_count` | Members per hindcast tercile, uncalibrated | `compute_products_func.py`, `compare_terciles.py` |
| `cca` | africas2s CCA with OSOP's pyCPT settings | `run_pycpt.py` + pyCPT + CPT |
| `ereg` | africas2s ensemble regression | new, for comparison |

Every hindcast is scored against the observations with RPSS and generalized
ROC, so the three methods can be compared on the same inputs.

| `rebuild:` key | Default | Meaning |
|---|---|---|
| `backend` | `acmaddl` | `acmaddl` fetches from source; `grib` reads an OSOP downloads directory |
| `grib_dir` | none | Directory of OSOP GRIB files, for `backend: grib` |
| `obs` | `obs/era5` | acmadDL observation product |
| `obs_grid_res`, `model_grid_res` | `1.0`, none | Regrid on fetch, in degrees |
| `products` | none | `centre: acmadDL product`, bypassing the C3S lookup |
| `methods` | all three | Which products to run |
| `cca_args` | OSOP's | Overrides for the CCA settings, e.g. `mode_selection: north` |
| `hindcast_years` | `[1993, 2016]` | OSOP hard-codes this period |
| `dry_mask` | `0.1` | mm/day. Cells with a drier observed season are blanked in every product and score; `null` disables |

## OSOP step to module

| OSOP step | OSOP code | Here |
|---|---|---|
| Download hindcast, forecast, ERA5 | `get_any_hindcast.py`, `get_era5.py` | `data.py`: `acmaddl.fetch` |
| GRIB to pyCPT NetCDF | `pycpt_convert.py` (1,044 lines) | not needed |
| Member-count terciles and MME | `compute_products_func.py`, `compare_terciles.py` | `products.py` (65 lines) |
| pyCPT CCA and MME | `run_pycpt.py`, pycpt, cptcore, CPT | `pipeline.py`: `africas2s.seasonal_mme` |
| Scores | `compute_scores_func.py`, xskillscore | `pipeline.py`: `africas2s.skill` |
| Plots | `plot_verify.py`, `ens_plotting.py` | `pipeline.py`: `africas2s.plot_matrix` |
| Driver | `run_all.py`, subprocess per step | `pipeline.py`, `__main__.py`, one process |

The package is about 700 lines against OSOP's 6,600.

OSOP names a CDS system number per centre; acmadDL pins one system per catalog
entry. `data.product_for` looks the pair up in the acmadDL catalog and stops
with a clear message when there is no entry for that system, so a config never
silently runs a different model version.

## What is and is not verified

Verified:

- **Member-count parity with OSOP.** `tests/test_osop_parity.py` runs OSOP's
  own library functions on OSOP's test GRIBs and compares. Hindcast
  probabilities are identical for ECMWF, UKMO and the MME. Forecast
  probabilities are identical except in under 0.5% of cells, where one member
  differs: desert cells with about 0.00001 mm/day where a member sits within
  float32 rounding of a tercile boundary.
- **Season aggregation.** The season means read from GRIB correlate at 0.99999
  with the output of OSOP's `pycpt_convert.py` for both models and for ERA5.
- **OSOP's stock config** loads unchanged. Its `ecmwf 51` and `dwd 22` map to
  `c3s/ecmwf-monthly` and `c3s/dwd`; its `ukmo 604` stops with a clear message,
  because acmadDL carries system 610 only.
- **The whole pipeline runs offline** on OSOP's test data: three methods, two
  centres plus MME, forecast maps, skill table and skill maps.
- **OSOP's stock settings run end to end on live C3S data**
  (`configs/egypt_c3s.yml`): DWD 22, UKMO 610 and ECMWF 51 hindcasts and
  November 2026 forecasts over the 40N-40S predictor band, ERA5 observations,
  all through acmadDL.
- **The same path runs with credential-free products**
  (`configs/nmme_live_test.yml`: three NMME models, CMAP observations).
- **2 m temperature runs end to end** (`configs/horn_c3s_t2m.yml`): values
  arrive in degrees C, all three
  methods, forecast and skill. Not compared with OSOP or pyCPT output.

Not verified:

- **Automatic CCA mode selection** matches CPT's in the four models compared
  (see below), with africas2s #36 (version 0.4.1). Four models in two regions
  is a small sample.
- The pyCPT comparison was made with UKMO 604 and a 2025 forecast, before
  acmadDL moved to 610. It has not been repeated on 610.
- OSOP's 1-month products, Brier score, reliability diagrams and correlation
  maps are not ported. Only the season (`nm`) products are.

## Findings worth passing on

1. `tailoring="Anomaly"` breaks `africas2s.seasonal_mme` (0.3.0). OSOP sets it
   for pyCPT, where it only changes how the deterministic map is presented. In
   africas2s the forecast becomes an anomaly that is then compared with
   absolute tercile boundaries: 74% "below" everywhere on the test data. This
   package leaves the setting out.
2. `seasonal_mme(skill_metrics=[...])` raises for any explicit list, because
   the same list is applied to the deterministic and the probabilistic
   forecast and no metric accepts both. The example in the africas2s skill
   docs passes such a list.
3. `acmaddl.fetch(obs, target="OND", seasonal="mean")` downloads all twelve
   months of CHIRPS (288 files for 24 years); `months=[10, 11, 12]` downloads
   72. This package passes `months=` whenever the season stays in one year.
4. Calibrated methods report "above normal" over desert. Where most seasons
   are rainless both tercile boundaries sit near zero, and CCA and eReg then
   give "above" 53-58% of the time in every hindcast year (southern Egypt in
   winter: 68% of cell-seasons under 0.01 mm/day). OSOP runs pyCPT with
   `drymask_threshold: None`, so its pyCPT maps probably show the same thing;
   that is untested. This package masks cells under 0.1 mm/day by default.
5. With OSOP's pyCPT settings the CCA is the weakest of the three methods in
   both C3S cases below, and below 0.5 generalized ROC over Egypt. The
   cross-validated mode search picks 4 to 10 modes from 24 years at a goodness
   score near zero. `mode_selection: north` gave much the same skill over Egypt
   in 1.5 s per model instead of 3.5 minutes.
6. africas2s's CPT-parity settings are off by default and are required to
   match pyCPT (next section). This package switches them on. In africas2s
   0.3.0 the mode search ignored them and so picked different modes from CPT;
   fixed in africas2s #34.
7. pyCPT's probabilistic hindcast scores are in-sample. Its own source says so
   (`pycpt/notebook.py`: "this produces unrealistically optimistic values").
   The RPSS and ROC maps OSOP's pyCPT step plots are therefore not comparable
   with a cross-validated score.
8. OSOP's forecast member count uses strict inequalities on both sides, so a
   member exactly on a boundary is counted in no category. Harmless with real
   data, but the three categories then do not sum to 100%.

## Comparison with pyCPT

`pycpt_ref/` compares africas2s CCA with real pyCPT 2.10.4 (CPT binary in the
Docker image from `accord-research/pycpt`) on identical input files. The pyCPT
side makes the same calls with the same settings as OSOP's `run_pycpt.py`.

Two cases, both ECMWF 51 and UKMO 604 against ERA5, hindcast 1993-2016,
forecast 2025:

- `horn`: Horn of Africa, September start, October to December. Real skill.
  `configs/horn_c3s.yml` has since moved to UKMO 610 and a 2026 forecast, so
  `make_case.py horn` now builds a different case from the one tabulated here.
- `sahara`: OSOP's bundled test data, 10-40N 10-40E, June to August. No skill.

```bash
cd pycpt_ref
python make_case.py horn ../configs/horn_c3s.yml   # or: make_case.py sahara --osop-test
./run.sh horn                                      # pyCPT, in Docker
python find_modes.py horn                          # which EOF counts CPT selected
python africas2s_side.py horn cpt_modes ECMWF=1,1,1 UKMO=9,2,2
python compare.py horn
```

### Given the same modes, africas2s reproduces pyCPT

With africas2s's CPT-parity settings on and CPT's mode counts:

| Case, model | CPT modes (X, Y, CCA) | CV hindcast: RMS difference / RMS | CV hindcast correlation | 2025 probabilities: mean abs. difference | Same leading category |
|---|---|---|---|---|---|
| Horn, ECMWF | 1, 1, 1 | 1.5% | 0.9999 | 0.004 | 99% |
| Horn, UKMO | 9, 2, 2 | 0.2% | 1.0000 | 0.003 | 100% |
| Sahara, ECMWF | 10, 4, 2 | 1.1% | 0.9997 | 0.003 | 99% |
| Sahara, UKMO | 7, 9, 2 | 1.0% | 0.9996 | 0.003 | 99% |

Correlation of the cross-validated hindcast with observations is the same to
three decimals in all four (for example +0.583 for Horn UKMO in both).

pyCPT does not report the EOF counts CPT selected: its output files hold every
EOF computed, ten here, whatever was selected. `find_modes.py` infers them by
fitting every pair; in each model exactly one pair matches, at 0.2-1.5%
against 7-29% for the next best.

### The parity settings are required

Same modes, africas2s defaults instead of the parity settings:

| Case, model | CV hindcast: RMS difference / RMS | CV hindcast correlation | Same leading category |
|---|---|---|---|
| Horn, ECMWF | 32% | 0.956 | 94% |
| Horn, UKMO | 40% | 0.867 | 83% |
| Sahara, ECMWF | 135% | 0.424 | 68% |
| Sahara, UKMO | 101% | 0.467 | 80% |

The parity settings are `standardize=True`, `lat_weights="cpt"`,
`sign_convention="cpt"`, `y_reconstruction="cpt"` and `cv="loyo_cyclic"`. This
package passes them by default.

### Mode selection is where the two differ

africas2s's CPT-style search (`mode_selection="cpt"`) chose CPT's modes in one
model of four:

| Case, model | CPT | africas2s 0.3.0 search | Search redone with parity settings |
|---|---|---|---|
| Horn, ECMWF | 1, 1, 1 | 1, 1, 1 | 10, 10, 3 (CPT's choice third, 0.012 behind) |
| Horn, UKMO | 9, 2, 2 | 2, 9, 2 | 9, 2, 2 |
| Sahara, ECMWF | 10, 4, 2 | 8, 10, 7 | 10, 4, 2 |
| Sahara, UKMO | 7, 9, 2 | 10, 10, 6 | 6, 9, 2 (CPT's choice second, 0.001 behind) |

The stock search (`africas2s.methods.cca.select_modes`) ranks combinations
using a CCA with default settings and plain leave-out cross-validation, so it
ranks a different CCA from the one CPT ranks. `parity_search.py` repeats the
search with the parity settings and CPT's cyclic cross-validation: it then
picks CPT's modes in two of four and places CPT's choice within 0.012 of the
top in the other two.

This is fixed in africas2s `main` (#34, after 0.3.0): `select_modes` takes the
fit's settings and CV scheme, and `seasonal_mme` passes them, so
`mode_selection="cpt"` now gives the third column above from a normal run.

The remaining two differences have two causes, both found in CPT's source
(17.8.3) and confirmed by test:

1. **CPT excludes over-fitted combinations.** If any cross-validation fold of
   an (X modes, Y modes) pair has a canonical correlation of 1, CPT drops that
   pair (`cca.F95`, "Check for over-fitting"). With 19 training years, 10 + 10
   and 10 + 9 modes fit exactly. This is why the search picked 10, 10, 3 for
   Horn ECMWF.
2. **CPT's Kendall tau treats ties differently.** CPT computes tau as
   2 x AFC - 1 (`scores.F95`, `k_tau`, `two_afc_cc`); scipy's tau-b corrects
   for ties. They differ only where a cell has tied values, such as several
   rainless years. This decided Sahara UKMO (7, 9, 2 against 6, 9, 2, a margin
   of 0.0002 to 0.001).

Both are in africas2s from #36 (version 0.4.1): the over-fitting exclusion is on by default, and CPT's tau is chosen
with `cpt_args={"mode_selection_tau": "cpt"}`, which this package sets.
`pycpt_ref/check_select_modes.py` runs the search:

| Case, model | CPT | Search, tau-b | Search, CPT's tau |
|---|---|---|---|
| Horn, ECMWF | 1, 1, 1 | 1, 1, 1 | 1, 1, 1 |
| Horn, UKMO | 9, 2, 2 | 9, 2, 2 | 9, 2, 2 |
| Sahara, ECMWF | 10, 4, 2 | 10, 4, 2 | 10, 4, 2 |
| Sahara, UKMO | 7, 9, 2 | 6, 9, 2 | 7, 9, 2 |

Caveats: four models; the Docker image runs a later CPT than the source read;
in africas2s 0.4.0 and earlier the `mode_selection_tau` setting is silently
ignored and there is no over-fitting exclusion, so the third column of the
earlier table applies.

### pyCPT's probabilistic skill is in-sample

pyCPT's probabilistic hindcasts scored RPSS +0.19 / +0.26 and generalized ROC
0.69 / 0.78 in the Horn case. africas2s with the same modes, cross-validated,
scored +0.13 / +0.16 and 0.66 / 0.70. The difference is not skill: pyCPT fits
those hindcasts in-sample, and says so in its source (`pycpt/notebook.py`:
"this produces unrealistically optimistic values"). The probabilistic skill
maps OSOP's pyCPT step plots overstate skill by that margin.

## Results

RPSS and generalized ROC of the cross-validated hindcast, 1993-2016, dry cells
masked. The member-count rows are scored in-sample, as OSOP produces them
(each year's tercile boundaries include that year), which flatters them
slightly. CCA and eReg are cross-validated with a five-year window. Each table
is one season in one region: enough to show the pipeline works, not enough to
rank methods.

### OSOP's stock config: Egypt, December to February, C3S, ERA5

`configs/egypt_c3s.yml`, November start, 21-34N 22-38E, predictor band 40N-40S,
DWD 22, UKMO 610 and ECMWF 51, forecast 2026. 53% of cells are masked as dry.
CCA modes chosen by the search: DWD 10, 3, 3; UKMO 10, 3, 3; ECMWF 8, 2, 1.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| member count | DWD | -0.014 | 0.514 |
| member count | UKMO | 0.024 | 0.565 |
| member count | ECMWF | 0.012 | 0.566 |
| member count | MME | 0.034 | 0.574 |
| CCA | DWD | -0.071 | 0.459 |
| CCA | UKMO | -0.059 | 0.475 |
| CCA | ECMWF | -0.061 | 0.460 |
| CCA | MME | -0.058 | 0.457 |
| eReg | DWD | -0.041 | 0.553 |
| eReg | UKMO | -0.014 | 0.589 |
| eReg | ECMWF | -0.023 | 0.578 |
| eReg | MME | 0.002 | 0.594 |

### OSOP test data, offline: 10-40N 10-40E, June to August

May start, ECMWF and UKMO from OSOP's bundled GRIBs, ERA5, forecast 2025. The
predictor box equals the target box. 59% of cells are masked as dry.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| member count | ECMWF | -0.056 | 0.541 |
| member count | UKMO | 0.026 | 0.587 |
| member count | MME | 0.007 | 0.574 |
| CCA | ECMWF | -0.058 | 0.505 |
| CCA | UKMO | -0.045 | 0.524 |
| CCA | MME | -0.039 | 0.529 |
| eReg | ECMWF | -0.071 | 0.548 |
| eReg | UKMO | -0.020 | 0.604 |
| eReg | MME | -0.029 | 0.589 |

### NMME, no credentials: Horn of Africa, October to December

September start, 5S-15N 33-48E, CMAP as truth, forecast 2024,
`mode_selection: north`. No cells masked.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| member count | CFSv2 | 0.008 | 0.576 |
| member count | GEOS-S2S | -0.045 | 0.593 |
| member count | SPEAR | 0.095 | 0.633 |
| member count | MME | 0.134 | 0.646 |
| CCA | CFSv2 | 0.025 | 0.568 |
| CCA | GEOS-S2S | 0.042 | 0.568 |
| CCA | SPEAR | 0.081 | 0.606 |
| CCA | MME | 0.079 | 0.606 |
| eReg | CFSv2 | 0.028 | 0.610 |
| eReg | GEOS-S2S | 0.043 | 0.627 |
| eReg | SPEAR | 0.122 | 0.682 |
| eReg | MME | 0.106 | 0.690 |

The GEOS-S2S member-count forecast is "above normal" at 70-90% almost everywhere:
its 2024 forecast mean is 50% wetter than its hindcast mean, a model-version
drift that counting members against hindcast terciles cannot correct.

### Precipitation against ERA5: Horn of Africa, October to December, C3S

`configs/horn_c3s.yml`, September start, ECMWF 51 and UKMO 610, forecast 2026.
2% of cells are masked as dry. CCA modes chosen by the search: ECMWF 1, 1, 1;
UKMO 8, 10, 3.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| Member count | ECMWF | 0.179 | 0.669 |
| Member count | UKMO | 0.203 | 0.681 |
| Member count | MME | 0.211 | 0.682 |
| CCA | ECMWF | 0.132 | 0.657 |
| CCA | UKMO | 0.121 | 0.662 |
| CCA | MME | 0.142 | 0.669 |
| eReg | ECMWF | 0.168 | 0.696 |
| eReg | UKMO | 0.192 | 0.710 |
| eReg | MME | 0.200 | 0.711 |

### 2 m temperature: Horn of Africa, October to December, C3S, ERA5

`configs/horn_c3s_t2m.yml`, September start, ECMWF 51 and UKMO 610, forecast
2026. CCA modes chosen by the search: ECMWF 5, 5, 3; UKMO 9, 7, 4.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| Member count | ECMWF | 0.190 | 0.693 |
| Member count | UKMO | 0.048 | 0.656 |
| Member count | MME | 0.167 | 0.684 |
| CCA | ECMWF | 0.091 | 0.648 |
| CCA | UKMO | 0.126 | 0.669 |
| CCA | MME | 0.120 | 0.666 |
| eReg | ECMWF | 0.174 | 0.705 |
| eReg | UKMO | 0.131 | 0.668 |
| eReg | MME | 0.181 | 0.708 |

The 2026 forecast is "above normal" nearly everywhere in every method (MME
probability 78-96% averaged over the box). That is probably mostly the warming
trend against the 1993-2016 reference period, not a seasonal signal. Part of
the hindcast skill is trend too.

### Precipitation against CHIRPS: Horn of Africa, October to December, C3S

September start, ECMWF 51 only, forecast 2025, scored against CHIRPS v2 on land
at 1 degree. CCA modes chosen by the search: 1, 1, 1. The config is not
included, because it has not been rerun since 12 November 2026 and a first run
is a slow CHIRPS download (see the caution below). It is `configs/horn_c3s.yml`
with `obs: obs/chirps-v2-monthly`, `obs_request_interval: 5.0`, ECMWF as the
only centre and `forecast_year: 2025`.

| Method | Centre | RPSS | Generalized ROC |
|---|---|---|---|
| Member count | ECMWF | 0.111 | 0.623 |
| CCA | ECMWF | 0.035 | 0.597 |
| eReg | ECMWF | 0.085 | 0.649 |

## Data-source caution

**CHIRPS (`data.chc.ucsb.edu`).** CHC publishes its limits in the README of its
own downloader (https://github.com/chc-ucsb/chc-data-download, "Downloading
Politely"): an address that keeps sending more than about 2 requests a second
is blocked for several hours, more than 20 open connections get HTTP 429, and
the advice is about 1 request a second from one connection.

acmadDL up to 0.3.1 did not observe this. It paced file opens, but GDAL's
`/vsicurl/` sent about 60 requests per CHIRPS file within a fraction of a second (mostly
probes for sidecar files that do not exist). A 72-file pull at 6 seconds per
file got this address blocked. acmadDL #27 (version 0.3.2) replaces that with its own reader: 3 or 4 paced requests per file.

With that reader the 72-file Horn pull completed, in 2 h 46 min. The server
answered 43 requests with a plain HTTP 429, at spacings from 5 to 30 seconds
between requests and about as often at 30 as at 5. Each cleared on a retry and
no block followed. That is far below the published limit and the cause is not
known. The reader waits each 429 out and stops the pull if one persists.

Do not retry a CHIRPS pull straight after a block: the retry lands inside it.
The result is cached, so later runs of the same config make no requests.

**UKMO system number.** acmadDL's `c3s/ukmo` is system 610, which has
forecasts from April 2026 only (604 covered March 2025 to February 2026). The
configs here name 610 and a 2026 forecast. A config that names 604, including
OSOP's own `osop_config.yml`, stops with "no monthly C3S product for centre
'ukmo' system 604". acmadDL's cache key does not include the system number: a
machine that cached UKMO data under 604 is served that data as 610 until those
cache entries are removed.

Long runs: start them under `caffeinate -i`. A laptop that sleeps mid-download
leaves the CDS client waiting on a dead connection indefinitely.
