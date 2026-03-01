# Model Overview

## What pySPEEDY is

pySPEEDY is a Python wrapper around a refactored Fortran implementation of the SPEEDY
intermediate-complexity atmospheric general circulation model. The Python layer handles
configuration, I/O, callbacks, and model-state access. The Fortran layer does the actual
time integration and physical parameterizations.

At a high level:

- Python allocates a model state, loads boundary conditions, and drives the run loop.
- Fortran owns the prognostic state and performs the dynamics and physics.
- The model can expose state variables back to Python through generated getters and setters.

## Repository map

The main pieces are:

- `pyspeedy/speedy.py`: Python `Speedy` and `SpeedyEns` APIs.
- `pyspeedy/callbacks.py`: output/export and diagnostic callbacks.
- `pyspeedy/config.py`: shared runtime/config helpers.
- `pyspeedy/data/model_config.yml`: editable defaults for model and run settings.
- `speedy.f90/`: Fortran core.
- `examples/Run_PySpeedy_1week.py`: current single-run example.

The main Fortran modules are:

- `speedy.f90/initialization.f90`: full model initialization.
- `speedy.f90/tendencies.f90`: dynamics plus physics tendency assembly.
- `speedy.f90/time_stepping.f90`: leapfrog and filtering.
- `speedy.f90/physics.f90`: moist physics, radiation, surface fluxes, PBL tendencies.
- `speedy.f90/forcing.f90`: daily forcing updates.
- `speedy.f90/land_model.f90` and `speedy.f90/sea_model.f90`: slab land/sea updates.

## State variables

The core prognostic atmospheric variables are:

- vorticity
- divergence
- temperature
- log of normalized surface pressure
- tracers, currently specific humidity

In Python these live primarily in spectral space, but the wrapper exposes grid-space variables
such as:

- `u_grid`
- `v_grid`
- `t_grid`
- `q_grid`
- `phi_grid`
- `ps_grid`

The full generated registry is in `pyspeedy/data/model_state.json`.

## Default grid and vertical coordinate

With the current repo defaults:

- horizontal truncation: `T30`
- Gaussian grid: `96` longitudes by `48` latitudes
- vertical coordinate: hydrostatic sigma coordinate
- vertical levels: `8`

The sigma half levels are defined in the geometry module and the full-level sigma coordinates
are exposed back to Python as `lev`.

## What happens during a run

The normal single-run flow is:

1. Create `Speedy(start_date=..., end_date=...)`.
2. Call `model.set_bc(...)` to load monthly climatological BCs and SST anomalies.
3. Fortran initialization builds geometry, spectral transforms, implicit operators, boundary masks,
   prognostic initial conditions, and coupler state.
4. `model.run(callbacks=...)` advances one timestep at a time until `current_date < end_date`
   is no longer true. In other words, `end_date` is exclusive.
5. Callbacks can export output, accumulate checkpoints, or run diagnostic checks.

## Initialization details

The current initialization path starts from a reference atmosphere at rest rather than from a
restart file. The reference state is described in `speedy.f90/prognostics.f90`:

- zero wind and zero tracer perturbations
- troposphere with a fixed lapse-rate structure
- stratosphere capped at `216 K`
- surface pressure consistent with that reference temperature profile
- humidity initialized from a reference relative humidity profile

## Boundary-condition data flow

The current boundary-condition interface expects:

- invariant fields on `(lon, lat)`:
  - `orog`, `lsm`, `vegl`, `vegh`, `alb`
- monthly climatologies on `(lon, lat, time)` with `time=12`:
  - `stl`, `snowd`, `swl1`, `swl2`, `swl3`, `icec`, `sst`
- monthly SST anomalies on `(lon, lat, time)` for the simulation period:
  - `ssta`

Important current behavior:

- The shipped example BC and SST anomaly datasets are `T30`.
- If the configured model grid differs, `pyspeedy/speedy.py` remaps them to the configured grid.
- The code now prints whether interpolation is applied.

## Daily versus sub-daily updates

The atmosphere advances every timestep, but some coupled fields update on slower schedules:

- shortwave radiation is evaluated every `nstrad` timesteps
- forcing fields and land/sea exchanges are updated once per model day

That split is important when you change timestep settings or design idealized experiments.
