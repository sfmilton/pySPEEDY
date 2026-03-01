===================
pySPEEDY user guide
===================

Description of the SPEEDY model
================================

SPEEDY is an atmospheric global circulation model that uses the spectral dynamical core developed
at the Geophysical Fluid Dynamics Laboratory.
Its main characteristics are:

- A spectral-transform model in the vorticity-divergence form with a semi-implicit treatment of gravity waves.
- Hydrostatic σ-coordinate in the vertical coordinate.
- The principal model prognostic variables are vorticity, divergence, temperature, and the logarithm of surface pressure.
- Humidity is advected by the dynamical core, and with its sources and sinks are determined physical parametrizations.
- The default horizontal resolution in this repository is a triangular spectral truncation at total wavenumber 30 (
  T30, approximately 3.75 x 3.75 degree resolution).
  This corresponds to a [Gaussian grid](https://en.wikipedia.org/wiki/Gaussian_grid) of 96 (longitude) by 48 (latitude) points.
- For the vertical coordinate, eight levels are used with boundaries at σ values of 0, 0.05, 0.14, 0.26, 0.42, 0.60, 0.77, 0.90 and 1.
  The prognostic variables (except log(ps)) are specified at the σ levels in between the latter boundaries, namely at
  σ 0.025, 0.095, 0.20, 0.34, 0.51, 0.685, 0.835 and 0.9
- Time step: 40min. That is, 36 time steps per day.

For additional details, check the SPEEDY documentation
`here <http://users.ictp.it/~kucharsk/speedy_description/km_ver41_appendixA.pdf>`_.

The boundary conditions
-----------------------

The boundary conditions needed to run SPEEDY are:

- Time invariant fields (lon, lat):
  - orog: Orographic height [m]
  - lsm: Land sea mask fraction. Values between 0 and 1.
  - vegl: Low vegetation cover (fraction). Values between 0 and 1.
  - vegh: High vegetation cover (fraction). Values between 0 and 1.
  - alb: Annual-mean albedo (fraction). Values between 0 and 1.
- Monthly-average climatological fiel for each month of the year (lon, lat, month):
  - stl: Land surface temp (top-layer) [degK].
  - snowd: Snow depth [kg/m2]
  - swl1: Soil wetness (layer 1) [vol. fraction, 0-1]
  - swl2: Soil wetness (layer 2) [vol. fraction, 0-1]
  - swl3: Soil wetness (layer 3) [vol. fraction, 0-1]
  - icec: Sea-ice concentration (fraction). Values between 0 and 1.
  - sst: Sea surface temperature [degK].
- Anomaly fields (lon, lat, day):
  - ssta: Sea surface temperature anomaly [degK].

The exact shapes for the invariant fields are `(lon, lat, month) = (96, 48, 12)`.
In contrast, the anomaly field (SST anomaly) needs to be provided for each day of the simulation period.
For example, for a one-month forecast (30 days), the shape of the anomaly field is (96, 48, 30).

By default, the pySPEEDY package includes the example boundary conditions initially included in the SPEEDY.f90 package.
These fields were derived from the ERA-interim re-analysis using the 1979-2008 period.
In addition, the default boundary conditions include the monthly SST anomalies for the 1979-01-01 to 2013-12-01 period.
If the configured model grid differs from the bundled data grid, pySPEEDY will regrid the bundled inputs to the
configured lon/lat coordinates before initializing the model.

The initial conditions
----------------------

The current version of pySPEEDY initializes all spectral variables at a reference atmosphere at rest with the following properties:
Troposphere: T = 288°K at the surface. Constant temperature lapse rate.
Stratosphere: T = 216 °K, lapse rate = 0.
Pressure field consistent with the temperature field.
Humidity: qv = RHref * q_saturation(288K, 1013hPa) with RHref=0.7.
For this initialization, the following spectral variables are initialized:
- "vor [mx, nx, levels, 2]": Vorticity in the spectral space. Complex.
- "div [mx, nx, levels, 2]": Divergence in the spectral space. Complex.
- "t [mx, nx, levels, 2]": Temperature in the spectral space. Complex.
- "ps [mx, nx, 2]": Log of (normalised) surface pressure.
- "phi [mx, nx, levels]": Atmospheric geopotential. Real.
- "phis [mx, nx]": Surface geopotential. Real.
- "tr [mx, nx, levels, 2, ntr]": Tracers. Currently, it only contains humidity. Real.

where mx=31, nx=32, levels=8, ntr=1 (number of tracers).
The "2" in the last dimension indicates that the variable is complex (0-real/1-imaginary).


Model State variables
=====================

A complete list of the model state variables in the pySPEEDY model is presented `in this
interactive table <./model_state_def.html>`_.
The values of these variables for a particular SPEEDY instance can be accessed/modified the bracket getter/setters
in python. For example::

    model = Speedy() # Create an speedy instance
    model.set_bc()

    # Read the temperature
    temp = model["t_grid"] # This returns a copy of the t_grid fortran array.

    # Add a small increment to the temperature everywhere.
    model["t_grid"]=temp + 0.001



Running the model
=================

To run the model, we need to create an SPEEDY instance, set the boundary conditions, and finally run the model.

The next example shows how to run a 1 week simulation using the default boundary and initial conditions::

  from pyspeedy import Speedy
  from pyspeedy.callbacks import DiagnosticCheck, XarrayExporter
  from pyspeedy.config import load_config

  config = load_config()
  run_config = config.run

  # Create an instance of the speedy model using the defaults in pyspeedy/data/model_config.yml.
  model = Speedy(
      start_date=run_config.start_date,
      end_date=run_config.end_date,
  )
  # At this point, the model state is "empty".

  # To initialize the model, we need to define its boundary conditions first.
  # This function will set the default boundary conditions derived from the ERA reanalysis.
  model.set_bc()

  callbacks = [
      DiagnosticCheck(interval=run_config.diag_interval),
      XarrayExporter(
          interval=run_config.history_interval,
          output_dir=run_config.output_dir,
          variables=run_config.output_vars,
          verbose=run_config.verbose_output,
      ),
  ]

  # Print the names of output variables that will be saved.
  # Note that the variables shown next are in the grid space (not the spectral space)
  print(callbacks[1].variables)

  # Run the model
  model.run(callbacks=callbacks)
  # After the model is run, the model state will keep the last values of the last integration step.


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
