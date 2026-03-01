# Idealized Cases

## Current status

The current codebase does not support all idealized benchmarks equally.

| Case | Current status | Notes |
| --- | --- | --- |
| Standard moist SPEEDY with prescribed BCs | Supported | This is the default path. |
| Aquaplanet | Supported through custom BC files | No dedicated switch, but the current API is enough. |
| Held-Suarez | Supported | Use `physics_mode="held_suarez"` and the generated flat defaults or custom flat BCs. |

## Aquaplanet

### What aquaplanet means here

For this repo, an aquaplanet experiment means:

- no land
- no orography
- prescribed all-ocean lower boundary
- typically zonally symmetric SST
- usually no imposed SST anomaly field

There is no one-line `aquaplanet=True` mode today, but the current BC interface is flexible enough
to build one from a modified copy of the packaged BC files.

### Practical recipe

1. Start from `pyspeedy/data/example_bc.nc` and `pyspeedy/data/sst_anomaly.nc`.
2. Replace the land/sea mask with all sea.
3. Zero the topography.
4. Define a zonally symmetric SST climatology.
5. Zero out SST anomalies.
6. Run with `sst_anomaly_coupling_flag = False`.

### Example script

```python
from pathlib import Path
from datetime import datetime

import numpy as np
import xarray as xr

from pyspeedy import Speedy, example_bc_file, example_sst_anomaly_file
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter

work = Path("./aquaplanet_case")
work.mkdir(exist_ok=True)

bc = xr.load_dataset(example_bc_file(), engine="netcdf4")
ssta = xr.load_dataset(example_sst_anomaly_file(), engine="netcdf4")

# All ocean, no orography, no sea ice.
bc["lsm"] = xr.zeros_like(bc["lsm"])
bc["orog"] = xr.zeros_like(bc["orog"])
bc["icec"] = xr.zeros_like(bc["icec"])

# Keep all fields finite even if they are unused over ocean.
bc["vegh"] = xr.zeros_like(bc["vegh"])
bc["vegl"] = xr.zeros_like(bc["vegl"])
bc["alb"] = xr.zeros_like(bc["alb"]) + 0.06
bc["stl"] = xr.zeros_like(bc["stl"]) + 300.0
bc["snowd"] = xr.zeros_like(bc["snowd"])
bc["swl1"] = xr.zeros_like(bc["swl1"])
bc["swl2"] = xr.zeros_like(bc["swl2"])
bc["swl3"] = xr.zeros_like(bc["swl3"])

# Example zonally symmetric SST profile:
# warm tropics, cold poles, no zonal structure.
lat_rad = xr.DataArray(np.deg2rad(bc["lat"].values), dims=("lat",), coords={"lat": bc["lat"]})
sst_profile = 300.0 - 27.0 * np.sin(lat_rad) ** 2
bc["sst"] = xr.zeros_like(bc["sst"]) + sst_profile

# No anomaly forcing.
ssta["ssta"] = xr.zeros_like(ssta["ssta"])

bc_path = work / "aquaplanet_bc.nc"
ssta_path = work / "aquaplanet_ssta.nc"
bc.to_netcdf(bc_path)
ssta.to_netcdf(ssta_path)

model = Speedy(
    start_date=datetime(1980, 1, 1),
    end_date=datetime(1980, 1, 8),
)

model["sst_anomaly_coupling_flag"] = False
model.set_bc(bc_file=str(bc_path), sst_anomaly=str(ssta_path))

callbacks = [
    DiagnosticCheck(interval=36),
    XarrayExporter(interval=36, output_dir=str(work / "data"), verbose=True),
]

model.run(callbacks=callbacks)
```

### Notes on aquaplanet realism

This gives you an aquaplanet-style moist experiment with the current physics package, but it is not
the same thing as a dry benchmark or a standardized CMIP-style aquaplanet protocol.

The following are still active unless you change source code:

- moist physics
- cloud-radiation interaction
- slab sea/ice bookkeeping
- the default seasonal-cycle setting from the compiled config

## Held-Suarez

### Current status

Held-Suarez is now available as a real runtime mode in this repo.
The model-state registry exposes `held_suarez_mode`, and the Python API exposes it as
`physics_mode="held_suarez"`.

### What the current implementation does

When Held-Suarez mode is active:

- `speedy.f90/physics.f90` bypasses moist physics, radiation, surface fluxes, and PBL moist tendencies
- `speedy.f90/held_suarez.f90` applies:
  - Newtonian relaxation of temperature
  - Rayleigh drag on winds in the lower atmosphere
- humidity tendency is set to zero and the tracer is initialized dry
- `set_bc()` generates flat lower-boundary fields and zero SST anomalies if you do not supply your own inputs

### How to run it

The simplest path is:

```python
from datetime import datetime

from pyspeedy import Speedy
from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter

model = Speedy(
    start_date=datetime(1980, 1, 1),
    end_date=datetime(1980, 1, 8),
    physics_mode="held_suarez",
)

model.set_bc()

model.run(callbacks=[
    DiagnosticCheck(interval=36),
    XarrayExporter(interval=36, output_dir="./data", verbose=True),
])
```

There is also a dedicated example script in `examples/Run_PySpeedy_HeldSuarez.py`.

### Implementation notes

The current Held-Suarez mode is closest to the benchmark when you keep the lower boundary flat.
The generated defaults do that automatically:

- zero orography
- all ocean mask
- constant finite surface fields
- zero SST anomalies

You can still provide your own BC files, but then you are leaving the strict benchmark setup.

### What not to call Held-Suarez

These are not Held-Suarez:

- turning off SST anomalies only
- flattening the topography only
- making the planet all ocean while keeping moist physics on

Those can still be useful experiments, but they are not the standard Held-Suarez dry benchmark.
