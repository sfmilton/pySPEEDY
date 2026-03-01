# pySPEEDY Markdown Notes

This directory contains repo-oriented Markdown notes that complement the existing Sphinx docs.
They focus on how the current codebase is wired today rather than on the historical SPEEDY
documentation alone.

Files in this directory:

- [Model Overview](./model-overview.md): package layout, run flow, state variables, and data flow.
- [Equations and Numerics](./equations-and-numerics.md): what system is solved and how it is advanced.
- [Physics and Coupling](./physics-and-coupling.md): moist physics, radiation, surface fluxes, land, sea, and coupling.
- [Running and Configuration](./running-and-configuration.md): build/run workflow and the main user choices.
- [Idealized Cases](./idealized-cases.md): current status of aquaplanet and Held-Suarez in this repo.

Current repo facts worth knowing up front:

- The active default configuration is `T30`, `96 x 48`, `8` sigma levels, `36` timesteps per day.
- Output files written by `XarrayExporter` are prefixed with the physics-case and spectral tags, for example
  `SPEEDY_T30_1980-01-02_0000.nc` or `HS_T30_1980-01-02_0000.nc`.
- The bundled boundary-condition and SST-anomaly files are native `T30` files.
  If the configured grid changes, pySPEEDY remaps those inputs to the configured lon/lat grid at load time.
- A custom aquaplanet experiment can be run today by supplying custom boundary-condition files.
- Held-Suarez is now available as a runtime physics mode with generated flat lower-boundary inputs.

For the original SPEEDY scientific description, see the link already referenced in the main docs:
<http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf>
