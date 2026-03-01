# Physics and Coupling

## Physics sequence in the current model

`speedy.f90/physics.f90` applies the physical parameterizations in this order:

1. convert the current spectral state to grid-point fields
2. diagnose thermodynamic variables such as relative humidity and saturation humidity
3. deep convection
4. large-scale condensation
5. shortwave radiation
6. longwave radiation
7. surface fluxes
8. vertical diffusion and shallow convection
9. optional SPPT perturbation of the physical tendency contribution

## Deep convection

Deep convection is implemented in `speedy.f90/convection.f90` as a simplified Tiedtke-style
mass-flux scheme.

The scheme diagnoses unstable columns and then computes:

- cloud-base mass flux
- vertical moisture and dry-static-energy fluxes
- convective precipitation
- cloud-top level used by later cloud and radiation logic

## Large-scale condensation

`speedy.f90/large_scale_condensation.f90` relaxes humidity toward a sigma-dependent relative
humidity threshold and converts the resulting moisture sink into latent heating and precipitation.

That module is one of the most explicit in the repo and is a good starting point if you want to
change the moist physics.

## Radiation

Radiation is split into:

- `speedy.f90/shortwave_radiation.f90`
- `speedy.f90/longwave_radiation.f90`

Current characteristics:

- shortwave is not necessarily computed every timestep; it is gated by `nstrad`
- cloud properties from the moist physics feed back into radiative transmissivity
- surface and top-of-atmosphere radiative fluxes are stored in the model state

The state exposes fields such as:

- `tsr`
- `ssrd`
- `ssr`
- `slrd`
- `slr`
- `olr`
- `tt_rsw`

## Surface fluxes

`speedy.f90/surface_fluxes.f90` computes:

- momentum fluxes
- sensible heat flux
- evaporation
- net surface heat flux
- surface longwave emission

It works over land, sea, and area-weighted mixed surfaces.

## Boundary layer and shallow convection

`speedy.f90/vertical_diffusion.f90` combines:

- shallow convection between the lowest layers
- moisture diffusion in stable conditions
- damping of near-superadiabatic lapse rates

This is where lower-tropospheric mixing is handled in the current physics package.

## Land model

The land model in `speedy.f90/land_model.f90` is a slab land-surface model with:

- climatological land temperature, snow depth, and soil-water inputs
- prognostic land-surface temperature anomaly
- optional land coupling through `land_coupling_flag`

Important current behavior:

- snow depth and soil-water availability remain climatologically prescribed
- only land temperature anomaly is prognostic

## Sea and ice model

The sea module in `speedy.f90/sea_model.f90` handles:

- prescribed SST climatology
- observed monthly SST anomalies
- slab mixed-layer and sea-ice bookkeeping

Current limitations matter:

- `sea_coupling_flag` is compiled as `0`
- more strongly coupled sea options are described in comments but marked as not supported

So the practical default today is:

- prescribed SST climatology
- optionally plus observed monthly SST anomaly forcing

## Coupling schedule

The coupler (`speedy.f90/coupler.f90`) initializes land and sea state at startup and then updates
them once per model day. The atmospheric timestep is much shorter, so land/sea boundary exchange is
not happening every sub-daily step.

## Runtime state switches a user can change from Python

Several useful experiment switches already exist in the mutable model state:

- `increase_co2`: turn on the built-in CO2 optical-thickness trend used in `forcing.f90`
- `land_coupling_flag`: use climatological land temperature or evolve the slab-land anomaly
- `sst_anomaly_coupling_flag`: use SST anomalies or climatological SST only

These are not all exposed through `model_config.yml`; they can be changed directly on a model
instance before or after initialization depending on the experiment.

Example:

```python
model["increase_co2"] = True
model["land_coupling_flag"] = False
model["sst_anomaly_coupling_flag"] = False
```

## What this means for experiment design

This codebase is best thought of as:

- a moist intermediate-complexity AGCM with prescribed or weakly coupled lower-boundary forcing
- plus a dedicated Held-Suarez dry benchmark mode

Aquaplanet still works by supplying custom BC files. Held-Suarez is now implemented as a separate
runtime forcing path rather than as a BC-only variant.
