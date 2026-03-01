# Running and Configuration

## Build and run the standard example

From the repo root:

```bash
python -m pip install -e .
python examples/Run_PySpeedy_1week.py
```

That example uses:

- `pyspeedy/data/model_config.yml` for dates, output cadence, and default output directory
- `model.set_bc()` with the packaged example BC and SST-anomaly files
- `DiagnosticCheck` and `XarrayExporter` callbacks

With the current defaults the run is:

- start date: `1980-01-01`
- end date: `1980-01-08`
- one full week because `end_date` is exclusive
- one output file per day

## Output files

`XarrayExporter` now prefixes filenames with both the physics-case tag and the spectral tag.

Examples:

- `SPEEDY_T30_1980-01-02_0000.nc`
- `HS_T30_1980-01-02_0000.nc`

## Configuration layers

There are three different kinds of user choices in the current repo.

### 1. Compile-time model choices

These live in the `model` section of `pyspeedy/data/model_config.yml`.
If you change them, you must rebuild the extension.

Main controls:

- `trunc`: spectral truncation, for example `30` for `T30`
- `ix`, `iy`: horizontal Gaussian-grid dimensions
- `kx`: number of sigma levels
- `nsteps`: timesteps per day
- `nstrad`: shortwave-radiation cadence in timesteps
- `sppt_on`: compile-time SPPT switch
- `rob`, `wil`, `alph`: time-integration controls

### 2. Runtime default choices

These live in the `run` section of `pyspeedy/data/model_config.yml`.
They affect the example scripts and callback defaults, not the compiled Fortran geometry.

Main controls:

- `start_date`
- `end_date`
- `output_dir`
- `history_interval`
- `diag_interval`
- `verbose_output`
- `output_vars`

The `held_suarez` section is also runtime-configurable and does not require a rebuild.
It controls the benchmark constants such as:

- equilibrium temperature parameters
- boundary-layer sigma threshold
- Newtonian-cooling timescales
- Rayleigh-drag timescale

### 3. Per-run Python choices

These are passed directly when you build a `Speedy` object or callbacks.

Examples:

- choose dates at construction
- choose `physics_mode="speedy"` or `physics_mode="held_suarez"`
- choose a BC file and SST anomaly file in `set_bc(...)`
- choose callback intervals and output variables
- choose ensemble size with `SpeedyEns`
- modify state variables directly through `model["name"]`

## Minimal custom run

```python
from datetime import datetime

from pyspeedy import Speedy
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter

model = Speedy(
    start_date=datetime(1980, 1, 1),
    end_date=datetime(1980, 1, 8),
)

model.set_bc()

callbacks = [
    DiagnosticCheck(interval=36),
    XarrayExporter(interval=36, output_dir="./data", verbose=True),
]

model.run(callbacks=callbacks)
```

## Held-Suarez run

The repo now supports a dry Held-Suarez mode directly from Python:

```python
from datetime import datetime

from pyspeedy import Speedy

model = Speedy(
    start_date=datetime(1980, 1, 1),
    end_date=datetime(1980, 1, 8),
    physics_mode="held_suarez",
)

model.set_bc()
```

In that mode:

- the model applies Held-Suarez Newtonian cooling and Rayleigh drag
- moist physics and radiative/surface-flux tendencies are bypassed
- humidity is kept inert
- `set_bc()` generates flat default lower-boundary data if you do not provide your own files

## Boundary-condition requirements

The current `set_bc()` path expects these variables:

- invariant: `orog`, `lsm`, `vegl`, `vegh`, `alb`
- monthly climatologies: `stl`, `snowd`, `swl1`, `swl2`, `swl3`, `icec`, `sst`
- monthly SST anomalies from a separate dataset: `ssta`

In the packaged files:

- BC climatologies are stored with a `time=12` monthly dimension
- SST anomalies are stored monthly over a long time axis and interpolated in time by the sea model

## Native versus remapped BCs

The shipped BC and SST-anomaly files are `T30`.

Current behavior in `pyspeedy/speedy.py`:

- if the data already match the configured grid, pySPEEDY prints that no interpolation is applied
- if the configured grid differs, pySPEEDY remaps the data to the configured lon/lat grid before initialization

The current regridding is mask-aware for land-only and sea-only fields.

## Direct state edits from Python

The wrapper supports direct reads and writes of state variables.

Examples:

```python
temp = model["t_grid"]
model["t_grid"] = temp + 1.0
model.grid2spectral()
```

Useful mutable scalar flags include:

- `increase_co2`
- `held_suarez_mode`
- `land_coupling_flag`
- `sst_anomaly_coupling_flag`

Held-Suarez tuning scalars are also exposed through the model state if you want to override the
YAML defaults programmatically:

- `hs_trefc`
- `hs_delta_ty`
- `hs_delta_theta_z`
- `hs_tmin`
- `hs_sigma_b`
- `hs_tau_a_days`
- `hs_tau_s_days`
- `hs_tau_f_days`

## Ensemble runs

For multiple members, use `SpeedyEns`.
`XarrayExporter` writes each member to a separate subdirectory so outputs do not collide.

## Practical advice

If you want the least friction:

- keep the compiled grid at `T30`
- use the packaged BCs
- only change run dates and output intervals at first

If you want to change resolution:

- edit the `model` section of `model_config.yml`
- rebuild
- verify the BC-interpolation message on startup
- run a short smoke test before a long integration
