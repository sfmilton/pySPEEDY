# Equations and Numerics

## Governing system

At the level of the current code, pySPEEDY solves the hydrostatic primitive-equation system on
the sphere in sigma coordinates using a spectral-transform dynamical core.

The prognostic variables are:

- spectral vorticity `zeta`
- spectral divergence `D`
- spectral temperature `T`
- spectral `log(ps / p0)`
- tracer fields `q_n`, with humidity as tracer `1`

The model is described in the code and existing docs as:

- spectral-transform
- vorticity-divergence form
- hydrostatic sigma-coordinate
- semi-implicit for gravity-wave terms

## Grid and spectral representations

The model switches between:

- spectral space for the actual prognostic update
- Gaussian lon/lat grid space for most physical parameterizations

This conversion is done through the spectral module:

- `vort2vel`
- `grid_vel2vort`
- `spec2grid`
- `grid2spec`

## Pressure, geopotential, and humidity

Surface pressure is advanced as the logarithm of normalized pressure. Grid-space pressure is then
recovered as:

`ps_grid = p0 * exp(ps)`

Geopotential is diagnosed hydrostatically from temperature and surface geopotential in
`speedy.f90/geopotential.f90`.

Humidity is carried as a tracer in spectral space and transformed to grid space when physical
tendencies are computed.

## Tendency assembly

`speedy.f90/tendencies.f90` builds the full tendencies in two stages:

1. grid-point tendencies from dynamics plus physical parameterizations
2. spectral tendencies for divergence, temperature, and pressure including the semi-implicit terms

The code computes:

- advection terms in grid space
- pressure-gradient and Coriolis contributions
- sigma-dot vertical-motion terms
- tracer advection
- then adds physical tendencies from `physics.f90`

## Time stepping

The model uses a leapfrog family timestep with Robert/Williams filtering.

The scheme in the code comment is:

`F_new = F(1) + dt * [T_dyn(F(J2)) + T_phy(F(J1))]`

with filtered updates to the two time levels. The relevant controls are:

- `rob`: Robert filter strength
- `wil`: Williams-filter parameter
- `alph`: semi-implicit weighting

The first step is special:

- half step for initialization
- full step to seed the leapfrog state
- then repeated `2 * delt` leapfrog steps

## Semi-implicit and diffusion terms

The implicit module adds the reference-state terms used to stabilize gravity waves.
Horizontal diffusion is applied spectrally in `speedy.f90/horizontal_diffusion.f90`.

Current numerical damping features include:

- scale-selective horizontal diffusion
- stronger stratospheric diffusion
- stratospheric zonal-wind damping
- optional SPPT multiplicative stochastic perturbations on physical tendencies

## Diagnostic equations implemented explicitly in code

The clearest explicit formulas in the repo are in the moist-physics modules.

Large-scale condensation is implemented as relaxation toward a sigma-dependent humidity threshold:

`dq/dt = -(q - RH(sigma) q_sat) / tau_lsc`

with latent-heating tendency:

`dT/dt = -(L / cp) * dq/dt`

and precipitation diagnosed from the vertically integrated moisture sink.

## What is not solved here

The current codebase does not include:

- nonhydrostatic dynamics
- interactive full ocean dynamics
- an out-of-the-box aquaplanet experiment mode

Aquaplanet still needs custom boundary conditions. Held-Suarez is now available as a dedicated
runtime dry-physics mode. See
[Idealized Cases](./idealized-cases.md).
