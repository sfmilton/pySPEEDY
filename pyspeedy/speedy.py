"""
The Speedy model
================

.. autosummary::
    :toctree: ./generated/

    Speedy
    SpeedyEns
"""

import os

import xarray as xr
import numpy as np
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import json

from pyspeedy import (
    _speedy,  # noqa
    example_bc_file,
    example_sst_anomaly_file,
    PACKAGE_DATA_DIR,
    DEFAULT_OUTPUT_VARS,
)
from pyspeedy.config import grid_coordinates, load_config
from pyspeedy.error_codes import ERROR_CODES

_DEFAULT_PARAMS = dict(
    start_date=datetime(1982, 1, 1),
    end_date=datetime(1982, 1, 2),
)

# Make this dict immutable.
_DEFAULT_PARAMS = tuple(_DEFAULT_PARAMS.items())
with open(PACKAGE_DATA_DIR / "model_state.json") as fp:
    MODEL_STATE_DEF = json.load(fp)

DEFAULT_CONFIG = load_config()
MODEL_CONFIG = DEFAULT_CONFIG.model
HELD_SUAREZ_CONFIG = DEFAULT_CONFIG.held_suarez
VALID_PHYSICS_MODES = frozenset(("speedy", "held_suarez"))


def _grid_alignment(dataset):
    if "lon" not in dataset.coords or "lat" not in dataset.coords:
        return False, False
    target_lon, target_lat = grid_coordinates(MODEL_CONFIG)
    nominal_resolution_matches = (
        dataset.sizes["lon"] == len(target_lon) and dataset.sizes["lat"] == len(target_lat)
    )
    coordinate_matches = nominal_resolution_matches and np.allclose(
        dataset["lon"].values,
        target_lon,
        atol=1.0e-6,
    ) and np.allclose(
        dataset["lat"].values,
        target_lat,
        atol=1.0e-1,
    )
    return nominal_resolution_matches, coordinate_matches


def _grid_matches_config(dataset):
    _, coordinate_matches = _grid_alignment(dataset)
    return coordinate_matches


def _sort_lat_lon(dataset):
    if "lon" in dataset.coords:
        dataset = dataset.sortby("lon")
    if "lat" in dataset.coords:
        dataset = dataset.sortby("lat")
    return dataset


def _extend_longitudes(dataset):
    if "lon" not in dataset.coords:
        return dataset
    if dataset.sizes["lon"] == 0:
        return dataset

    lon_values = dataset["lon"].values
    if np.isclose(lon_values[-1], 360.0):
        return dataset

    edge = dataset.isel(lon=[0]).copy(deep=True)
    edge = edge.assign_coords(lon=np.array([360.0]))
    return xr.concat([dataset, edge], dim="lon")


def _extend_polar_caps(dataset):
    if "lat" not in dataset.coords:
        return dataset
    if dataset.sizes["lat"] == 0:
        return dataset

    pieces = [dataset]
    lat_values = dataset["lat"].values
    if lat_values[0] > -90.0:
        south = dataset.isel(lat=[0]).copy(deep=True)
        south = south.assign_coords(lat=np.array([-90.0]))
        pieces.insert(0, south)
    if lat_values[-1] < 90.0:
        north = dataset.isel(lat=[-1]).copy(deep=True)
        north = north.assign_coords(lat=np.array([90.0]))
        pieces.append(north)
    return xr.concat(pieces, dim="lat")


def _regrid_dataset(dataset):
    if "lon" not in dataset.coords or "lat" not in dataset.coords:
        return dataset

    target_lon, target_lat = grid_coordinates(MODEL_CONFIG)
    if _grid_matches_config(dataset):
        return dataset

    source_dataset = _sort_lat_lon(dataset)
    source_dataset = _extend_longitudes(source_dataset)
    source_dataset = _extend_polar_caps(source_dataset)

    target_lon_da = xr.DataArray(target_lon, dims="lon")
    target_lat_da = xr.DataArray(target_lat, dims="lat")

    linear = source_dataset.interp(lon=target_lon_da, lat=target_lat_da)
    nearest = source_dataset.interp(
        lon=target_lon_da,
        lat=target_lat_da,
        method="nearest",
    )
    return linear.combine_first(nearest)


def _interp_field(field, target_lon, target_lat, source_mask=None, fill_value=None):
    source_field = _sort_lat_lon(field)
    if source_mask is not None:
        source_mask = _sort_lat_lon(source_mask)
        source_field = source_field.where(source_mask)

    source_field = _extend_longitudes(source_field)
    source_field = _extend_polar_caps(source_field)

    target_lon_da = xr.DataArray(target_lon, dims="lon")
    target_lat_da = xr.DataArray(target_lat, dims="lat")

    linear = source_field.interp(lon=target_lon_da, lat=target_lat_da)
    nearest = source_field.interp(
        lon=target_lon_da,
        lat=target_lat_da,
        method="nearest",
    )
    result = linear.combine_first(nearest)
    if fill_value is not None:
        result = result.fillna(fill_value)
    return result


def _regrid_bc_dataset(dataset):
    if "lon" not in dataset.coords or "lat" not in dataset.coords:
        return dataset

    target_lon, target_lat = grid_coordinates(MODEL_CONFIG)
    target_ds = xr.Dataset(coords={"lon": target_lon, "lat": target_lat})

    source_lsm = _sort_lat_lon(dataset["lsm"])
    land_mask = source_lsm > 0.0
    sea_mask = source_lsm < 1.0

    target_ds["lsm"] = _interp_field(source_lsm, target_lon, target_lat, fill_value=0.0).clip(0.0, 1.0)
    for field_name in ("orog", "alb", "vegh", "vegl"):
        target_ds[field_name] = _interp_field(dataset[field_name], target_lon, target_lat)

    target_ds["stl"] = _interp_field(dataset["stl"], target_lon, target_lat, source_mask=land_mask, fill_value=273.0)
    target_ds["snowd"] = _interp_field(dataset["snowd"], target_lon, target_lat, source_mask=land_mask, fill_value=0.0)
    target_ds["swl1"] = _interp_field(dataset["swl1"], target_lon, target_lat, source_mask=land_mask, fill_value=0.0)
    target_ds["swl2"] = _interp_field(dataset["swl2"], target_lon, target_lat, source_mask=land_mask, fill_value=0.0)
    target_ds["swl3"] = _interp_field(dataset["swl3"], target_lon, target_lat, source_mask=land_mask, fill_value=0.0)
    target_ds["sst"] = _interp_field(dataset["sst"], target_lon, target_lat, source_mask=sea_mask, fill_value=273.0)
    target_ds["icec"] = _interp_field(dataset["icec"], target_lon, target_lat, source_mask=sea_mask, fill_value=0.0).clip(0.0, 1.0)

    if "time" in dataset.coords:
        target_ds = target_ds.assign_coords(time=dataset["time"])
    return target_ds


def _regrid_ssta_dataset(dataset, land_sea_mask=None):
    if "lon" not in dataset.coords or "lat" not in dataset.coords:
        return dataset

    target_lon, target_lat = grid_coordinates(MODEL_CONFIG)
    if land_sea_mask is None:
        sea_mask = None
    else:
        sea_mask = _sort_lat_lon(land_sea_mask) < 1.0

    target_ds = xr.Dataset(coords={"lon": target_lon, "lat": target_lat, "time": dataset["time"]})
    target_ds["ssta"] = _interp_field(
        dataset["ssta"],
        target_lon,
        target_lat,
        source_mask=sea_mask,
        fill_value=0.0,
    )
    return target_ds


def _assert_finite_fields(dataset, field_names, dataset_name):
    for field_name in field_names:
        values = np.asarray(dataset[field_name].values)
        if np.isfinite(values).all():
            continue
        bad_points = int(np.size(values) - np.count_nonzero(np.isfinite(values)))
        raise RuntimeError(
            f"{dataset_name} contains {bad_points} non-finite values after regridding in field "
            f"'{field_name}'."
        )


def _iter_month_starts(start_date, end_date):
    current = start_date.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    last = end_date.replace(day=1, hour=0, minute=0, second=0, microsecond=0)
    while current <= last:
        yield current
        current += relativedelta(months=1)


def _build_held_suarez_bc_dataset():
    lon, lat = grid_coordinates(MODEL_CONFIG)
    months = np.arange(12, dtype=np.int32)
    coords = {"lon": lon, "lat": lat, "time": months}

    surface = xr.DataArray(np.zeros((MODEL_CONFIG.ix, MODEL_CONFIG.il), dtype=np.float64), dims=("lon", "lat"), coords={"lon": lon, "lat": lat})
    monthly = xr.DataArray(
        np.zeros((MODEL_CONFIG.ix, MODEL_CONFIG.il, 12), dtype=np.float64),
        dims=("lon", "lat", "time"),
        coords=coords,
    )

    dataset = xr.Dataset(coords=coords)
    dataset["orog"] = surface
    dataset["lsm"] = surface
    dataset["vegl"] = surface
    dataset["vegh"] = surface
    dataset["alb"] = surface + 0.06
    dataset["stl"] = monthly + 300.0
    dataset["snowd"] = monthly
    dataset["swl1"] = monthly
    dataset["swl2"] = monthly
    dataset["swl3"] = monthly
    dataset["icec"] = monthly
    dataset["sst"] = monthly + 300.0
    return dataset


def _build_held_suarez_ssta_dataset(start_date, end_date):
    lon, lat = grid_coordinates(MODEL_CONFIG)
    ssta_start = start_date.replace(day=1, hour=0, minute=0, second=0, microsecond=0) - relativedelta(months=1)
    ssta_end = end_date.replace(day=1, hour=0, minute=0, second=0, microsecond=0) + relativedelta(months=1)
    month_starts = list(_iter_month_starts(ssta_start, ssta_end))

    zeros = np.zeros((MODEL_CONFIG.ix, MODEL_CONFIG.il, len(month_starts)), dtype=np.float64)
    surface = np.zeros((MODEL_CONFIG.ix, MODEL_CONFIG.il), dtype=np.float64)
    return xr.Dataset(
        data_vars={
            "ssta": (("lon", "lat", "time"), zeros),
            "lsm": (("lon", "lat"), surface),
        },
        coords={"lon": lon, "lat": lat, "time": month_starts},
    )


if _speedy is None:  # pragma: no cover - exercised only before build/install
    raise ModuleNotFoundError(
        "The compiled pySPEEDY extension is not available. Run `python -m pip install -e .` first."
    )


class Speedy:
    """
    Speedy model.
    """

    def __init__(
        self,
        start_date=datetime(1982, 1, 1),
        end_date=datetime(1982, 1, 2),
        physics_mode="speedy",
        member=None,
    ):
        """
        Constructor. Initializes the model.

        Parameters
        ----------
        start_date: datetime
            Model start date.
        end_date: datetime
            Model end date.
        physics_mode: str
            Physics configuration. Supported values are ``"speedy"`` and ``"held_suarez"``.
        member: int
            Member number. If None, the run is not considered part of an ensemble run.
            This value is used by the callback functions to determine how to treat the simulation
            (a single isolated run, or an ensemble member).
        """
        self._start_date = None
        self._end_date = None
        self.member_id = member  # member number
        self.is_ensemble_member = self.member_id is not None

        # Allocate the model state
        self._state_cnt = _speedy.modelstate_init()
        self.set_params(start_date=start_date, end_date=end_date)
        self._physics_mode = None
        self.set_physics_mode(physics_mode)

        self._initialized_bc = False
        self._initialized_ssta = False
        self.current_date = self.start_date

    def __del__(self):
        """Clean up."""
        _speedy.modelstate_close(self._state_cnt)
        _speedy.controlparams_close(self._control_cnt)

        self._dealloc_date(self._start_date)
        self._dealloc_date(self._end_date)

    def set_params(
        self,
        start_date=datetime(1982, 1, 1),
        end_date=datetime(1982, 1, 2),
    ):
        """
        Set the model's control parameters.

        Parameters
        ----------
        start_date: datetime
            Model start date.
        end_date: datetime
            Model end date.
        """
        self.start_date = start_date
        self.end_date = end_date

        if self.start_date > self.end_date:
            raise ValueError("The start date should be lower than the en date.")

        self._control_cnt = _speedy.controlparams_init(self._start_date, self._end_date)

        self._model_date = None

        self.current_date = start_date

        self.n_months = (
            (self.end_date.year - self.start_date.year) * 12
            + (self.end_date.month - self.start_date.month)
            + 1
        )

    @staticmethod
    def _dealloc_date(container):
        """Deallocate a datetime fortran object."""
        if container is not None:
            _speedy.close_datetime(container)

    def __getitem__(self, var_name):
        """Getter for state variables."""
        _getter = getattr(_speedy, f"get_{var_name}", None)
        if _getter is None:
            raise AttributeError(f"The state variable '{var_name}' does not exist.")

        time_dim = MODEL_STATE_DEF[var_name]["time_dim"]
        if time_dim:
            return _getter(self._state_cnt, getattr(self, time_dim))

        return _getter(self._state_cnt)

    def get_shape(self, var_name):
        """Get the shape of an state variable."""
        _getter = getattr(_speedy, f"get_{var_name}_shape", None)
        if _getter is None:
            raise AttributeError(
                f"The 'get-shape' method for the state variable "
                f"{var_name}' does not exist."
            )
        return tuple(_getter(self._state_cnt))

    def __setitem__(self, var_name, value):
        """Setter for state variables."""
        _setter = getattr(_speedy, f"set_{var_name}")
        if _setter is None:
            raise AttributeError(
                f"The setter for the state variable '{var_name}' does not exist."
            )

        is_array_func = getattr(_speedy, f"is_array_{var_name}")
        if is_array_func():
            if self.get_shape(var_name) != value.shape:
                raise ValueError("Array shape missmatch")
            value = np.asfortranarray(value)

            time_dim = MODEL_STATE_DEF[var_name]["time_dim"]
            if time_dim:
                return _setter(self._state_cnt, value, getattr(self, time_dim))

            return _setter(self._state_cnt, value)

        return _setter(self._state_cnt, value)

    @staticmethod
    def _get_fortran_date(container):
        """Get a datetime object from a fortran datetime"""
        return datetime(*_speedy.get_datetime(container))

    @staticmethod
    def _set_fortran_date(container, date_value):
        """Create a python datetime object from a fortran datetime"""
        Speedy._dealloc_date(container)
        if isinstance(date_value, datetime):
            return _speedy.create_datetime(
                date_value.year,
                date_value.month,
                date_value.day,
                date_value.hour,
                date_value.minute,
            )
        else:
            raise TypeError("The input value is not a datetime object.")

    def get_current_step(self):
        """Return the current step in the simulation."""
        return self["current_step"]

    @property
    def physics_mode(self):
        return self._physics_mode

    @property
    def case_tag(self):
        if self.physics_mode == "held_suarez":
            return "HS"
        return "SPEEDY"

    @property
    def output_tag(self):
        return f"{self.case_tag}_{MODEL_CONFIG.spectral_tag}"

    def _apply_held_suarez_defaults(self):
        self["hs_trefc"] = HELD_SUAREZ_CONFIG.equilibrium_surface_temperature
        self["hs_delta_ty"] = HELD_SUAREZ_CONFIG.equator_pole_temperature_contrast
        self["hs_delta_theta_z"] = HELD_SUAREZ_CONFIG.vertical_temperature_contrast
        self["hs_tmin"] = HELD_SUAREZ_CONFIG.minimum_equilibrium_temperature
        self["hs_sigma_b"] = HELD_SUAREZ_CONFIG.boundary_layer_sigma
        self["hs_tau_a_days"] = HELD_SUAREZ_CONFIG.upper_air_relaxation_days
        self["hs_tau_s_days"] = HELD_SUAREZ_CONFIG.boundary_layer_relaxation_days
        self["hs_tau_f_days"] = HELD_SUAREZ_CONFIG.rayleigh_drag_days
        self["hs_min_pressure_ratio"] = HELD_SUAREZ_CONFIG.minimum_pressure_ratio

    def set_physics_mode(self, physics_mode):
        """Set the physics mode for the current model state before initialization."""
        normalized = physics_mode.lower()
        if normalized not in VALID_PHYSICS_MODES:
            raise ValueError(
                f"Unsupported physics_mode '{physics_mode}'. Expected one of: "
                f"{', '.join(sorted(VALID_PHYSICS_MODES))}."
            )

        self["held_suarez_mode"] = normalized == "held_suarez"
        if normalized == "held_suarez":
            self["land_coupling_flag"] = False
            self["sst_anomaly_coupling_flag"] = False
            self["increase_co2"] = False
            self._apply_held_suarez_defaults()
        self._physics_mode = normalized

    @property
    def start_date(self):
        return self._get_fortran_date(self._start_date)

    @start_date.setter
    def start_date(self, value):
        self._start_date = self._set_fortran_date(self._start_date, value)

    @property
    def current_date(self):
        return self._get_fortran_date(self._model_date)

    @current_date.setter
    def current_date(self, value):
        self._model_date = self._set_fortran_date(self._model_date, value)

    @property
    def end_date(self):
        return self._get_fortran_date(self._end_date)

    @end_date.setter
    def end_date(self, value):
        self._end_date = self._set_fortran_date(self._end_date, value)

    def set_bc(self, bc_file=None, sst_anomaly=None):
        """
        Set the model boundary conditions from a Netcdf file.
        If no file is provided, the default boundary conditions from the original SPEEDY model are used.

        The boundary conditions file (`bc_file`) should contain the following fields:

        - Time invariant (lon, lat):
            - orog: Orographic height [m]
            - lsm: Land sea mask fraction. Values between 0 and 1.
            - vegl: Low vegetation cover (fraction). Values between 0 and 1.
            - vegh: High vegetation cover (fraction). Values between 0 and 1.
            - alb: Annual-mean albedo (fraction). Values between 0 and 1.
        - Climatological values for each month of the year (lon, lat, month):
            - stl: Land surface temp (top-layer) [degK].
            - snowd: Snow depth [kg/m2]
            - swl1: Soil wetness (layer 1) [vol. fraction, 0-1]
            - swl2: Soil wetness (layer 2) [vol. fraction, 0-1]
            - swl3: Soil wetness (layer 3) [vol. fraction, 0-1]
            - icec: Sea-ice concentration (fraction). Values between 0 and 1.
            - sst: Sea surface temperature [degK].

        In addition to the climatological fields, the SPEEDY model requires the
        sea surface temperature anomalies with respect to the the climatologies (lon, lat, day):
        These anomalies are loaded from the file specified in the `sst_anomaly` keyword, and it should contain the
        following field:

        - Anomalies fields (lon, lat, day):
            - ssta: Sea surface temperature anomaly [degK].

        By default, the anomalies available in the original SPEEDY model are used.
        Note that the anomaly fields should cover the simulation period.

        Notes
        -----
        The exact shapes for the invariant fields are:
        (lon, lat, month) = (ix, 2 * iy, 12)

        The SPEEDY boundary conditions are included in the pySPEEDY package.
        This climatology was derived from the ERA interim re-analysis using the 1979-2008 period.
        Also, the example data provided in the pySPEEDY package include the monthly SST anomalies
        from 1979-01-01 to 2013-12-01 (Y/m/d). If the input grid does not match the configured
        model grid, the dataset is regridded before the boundary conditions are applied.
        """

        if self._initialized_bc:
            raise RuntimeError(
                "The model was already initialized. Create a new instance if you need different boundary conditions."
            )

        if self.physics_mode == "held_suarez":
            if bc_file is None:
                bc_file = _build_held_suarez_bc_dataset()
            if sst_anomaly is None:
                sst_anomaly = _build_held_suarez_ssta_dataset(self.start_date, self.end_date)

        self._set_sst_anomalies(sst_anomaly=sst_anomaly)

        # In the model state, the variables follow the lon/lat dimension ordering.
        if bc_file is None:
            bc_file = example_bc_file()

        if isinstance(bc_file, str):
            if not os.path.isfile(bc_file):
                raise RuntimeError(
                    "The boundary conditions file does not exist.\n" f"File: {bc_file}"
                )
            ds = xr.load_dataset(bc_file, engine="netcdf4")
        elif isinstance(bc_file, xr.Dataset):
            ds = bc_file
        else:
            raise TypeError(f"Unsupported bc_file input: {type(bc_file)}")
        nominal_resolution_matches, coordinate_matches = _grid_alignment(ds)
        if coordinate_matches:
            print(
                f"Boundary conditions already match configured {MODEL_CONFIG.spectral_tag} grid; "
                "no interpolation applied."
            )
        elif nominal_resolution_matches:
            print(
                "Boundary conditions use the same nominal resolution as the configured "
                f"{MODEL_CONFIG.spectral_tag} grid, but lon/lat coordinates differ; remapping "
                "to the model grid."
            )
        else:
            print(
                "Interpolating boundary conditions from "
                f"{ds.sizes['lon']}x{ds.sizes['lat']} to configured "
                f"{MODEL_CONFIG.spectral_tag} grid ({MODEL_CONFIG.ix}x{MODEL_CONFIG.il})."
            )
        ds = _regrid_bc_dataset(ds)
        _assert_finite_fields(
            ds,
            (
                "orog",
                "lsm",
                "alb",
                "vegh",
                "vegl",
                "stl",
                "snowd",
                "swl1",
                "swl2",
                "swl3",
                "sst",
                "icec",
            ),
            "Boundary conditions dataset",
        )

        self["orog"] = ds["orog"].values
        self["fmask_orig"] = ds["lsm"].values
        self["alb0"] = ds["alb"].values

        self["veg_high"] = ds["vegh"].values
        self["veg_low"] = ds["vegl"].values

        self["stl12"] = ds["stl"].values
        self["snowd12"] = ds["snowd"].values

        self["soil_wc_l1"] = ds["swl1"].values

        self["soil_wc_l2"] = ds["swl2"].values
        self["soil_wc_l3"] = ds["swl3"].values

        self["sst12"] = ds["sst"].values

        self["sea_ice_frac12"] = ds["icec"].values

        error_code = _speedy.init(self._state_cnt, self._control_cnt)
        if error_code < 0:
            raise RuntimeError(ERROR_CODES[error_code])

        self.spectral2grid()  # Convert some spectral variables to the grid space.
        for field_name in ("u_grid", "v_grid", "t_grid", "q_grid", "phi_grid", "ps_grid"):
            values = np.asarray(self[field_name])
            if np.isfinite(values).all():
                continue
            raise RuntimeError(
                f"Model state contains non-finite values immediately after initialization in '{field_name}'."
            )
        self._initialized_bc = True

    def _set_sst_anomalies(self, sst_anomaly=None):
        """
        Load SST anomalies from netcdf file.

        **Important**: The SST anomalies need to be set before setting the ICs!

        Only the times between the simulation's start and end date are loaded.

        See the :meth:`set_bc` documentation for additional details on the expected fields.
        """

        if self._initialized_ssta:
            raise RuntimeError(
                "The SST anomaly was already initialized."
                " Create a new instance if you need different boundary conditions."
            )
        if sst_anomaly is None:
            if self.physics_mode == "held_suarez":
                sst_anomaly = _build_held_suarez_ssta_dataset(self.start_date, self.end_date)
            else:
                sst_anomaly = example_sst_anomaly_file()

        if isinstance(sst_anomaly, str):
            if not os.path.isfile(sst_anomaly):
                raise RuntimeError(
                    "The SST anomaly file does not exist.\n" f"File: {sst_anomaly}"
                )
            ds = xr.load_dataset(sst_anomaly)
        elif isinstance(sst_anomaly, xr.Dataset):
            ds = sst_anomaly
        else:
            raise TypeError(f"Unsupported sst_anomaly input: {type(sst_anomaly)}")

        land_sea_mask = None
        if "lsm" in ds:
            land_sea_mask = ds["lsm"]
        elif os.path.isfile(example_bc_file()):
            land_sea_mask = xr.load_dataset(example_bc_file(), engine="netcdf4")["lsm"]

        nominal_resolution_matches, coordinate_matches = _grid_alignment(ds)
        if coordinate_matches:
            print(
                f"SST anomalies already match configured {MODEL_CONFIG.spectral_tag} grid; "
                "no interpolation applied."
            )
        elif nominal_resolution_matches:
            print(
                "SST anomalies use the same nominal resolution as the configured "
                f"{MODEL_CONFIG.spectral_tag} grid, but lon/lat coordinates differ; remapping "
                "to the model grid."
            )
        else:
            print(
                "Interpolating SST anomalies from "
                f"{ds.sizes['lon']}x{ds.sizes['lat']} to configured "
                f"{MODEL_CONFIG.spectral_tag} grid ({MODEL_CONFIG.ix}x{MODEL_CONFIG.il})."
            )
        ds = _regrid_ssta_dataset(ds, land_sea_mask=land_sea_mask)

        #########################################################################
        # Select only the times in the dataset between the start and the end date

        # At each timestep, Speedy uses a 3-month window for the computations.
        # Correct the start/end dates to account for those dates.
        start_date = self.start_date.replace(
            day=1, hour=0, minute=0, microsecond=0
        ) - relativedelta(months=1)

        end_date = self.end_date.replace(
            day=1, hour=0, minute=0, microsecond=0
        ) + relativedelta(months=1, days=1)

        ds = ds.loc[
            dict(
                lon=slice(None),
                lat=slice(None),
                time=slice(start_date, end_date),
            )
        ]

        expected_months = (
            (end_date.year - start_date.year) * 12
            + (end_date.month - start_date.month)
            + 1
        )

        missing_months = expected_months - len(ds["time"])

        if missing_months > 0:
            raise RuntimeError(
                f"{missing_months} months are missing in the SST anomalies file for the period: "
                + start_date.strftime("%Y/%m/%d")
                + " , "
                + end_date.strftime("%Y/%m/%d")
                + ".\n "
            )

        _assert_finite_fields(ds, ("ssta",), "SST anomaly dataset")

        _speedy.modelstate_init_sst_anom(self._state_cnt, expected_months - 2)
        self["sst_anom"] = ds["ssta"]
        self._initialized_ssta = True

    def run(self, callbacks=None):
        """
        Run the model between the start and the end date (defined in the instance `start_date` and
        `end_date` attributes).

        Before running the model, the boundary conditions need to be set (see the :meth:`set_bc` method).

        Parameters
        ----------
        callbacks: iterable of callback functions
            Sequence of callback functions to be call every `interval` time steps. The interval is specified in each
            callback instance.
        """
        if callbacks is None:
            callbacks = list()
        if not self._initialized_bc:
            raise RuntimeError(
                "The SPEEDY model was not initialized. Call the `set_bc` method to initialize the model."
            )

        self.current_date = self.start_date
        dt_step = timedelta(seconds=MODEL_CONFIG.timestep_seconds)

        while self.current_date < self.end_date:
            error_code = _speedy.step(self._state_cnt, self._control_cnt)
            if error_code < 0:
                raise RuntimeError(ERROR_CODES[error_code])
            self.current_date += dt_step

            for callback in callbacks:
                callback(self)

    def grid2spectral(self):
        """Transform the grid u, v, t, qv, ps, and phi fields to the spectral domain."""
        _speedy.transform_grid2spectral(self._state_cnt)  # noqa

    def spectral2grid(self):
        """Transform the spectral u, v, t, qv, ps, and phi fields to the grid domain."""
        _speedy.transform_spectral2grid(self._state_cnt)  # noqa

    def to_dataframe(self, variables=None):
        """
        Return an xarray DataFrame with the current model state.
        """
        if variables is None:
            variables = DEFAULT_OUTPUT_VARS

        self.spectral2grid()
        data_vars = dict()

        for var in variables:
            # Append time dimension
            dims = MODEL_STATE_DEF[var]["nc_dims"] + ["time"]
            var_data = self[var][..., None].astype("float32")
            if self.is_ensemble_member:
                # Append ensemble dimension if needed.
                dims = dims + ["ens"]
                var_data = var_data[..., None]

            data_vars[MODEL_STATE_DEF[var]["alt_name"]] = (dims, var_data)

        coords = dict(
            lon=self["lon"],
            lat=self["lat"],
            lev=self["lev"],
            time=[self.current_date],
        )
        if self.is_ensemble_member:
            coords["ens"] = [self.member_id]

        output_ds = xr.Dataset(
            data_vars=data_vars,
            coords=coords,
        )

        encoding = dict()
        for var in list(variables) + list(coords.keys()):
            if var in ("time", "ens"):
                encoding[var] = {"dtype": "int32"}
                continue
            alt_name = MODEL_STATE_DEF[var]["alt_name"]
            if MODEL_STATE_DEF[var]["units"] is not None:
                output_ds[alt_name].attrs["units"] = MODEL_STATE_DEF[var]["units"]
            output_ds[alt_name].attrs["long_name"] = MODEL_STATE_DEF[var]["desc"]
            output_ds[alt_name].attrs["standard_name"] = MODEL_STATE_DEF[var][
                "std_name"
            ]
            encoding[alt_name] = {"dtype": "float32", "zlib": True}

        output_ds["lat"].attrs["axis"] = "Y"
        output_ds["lon"].attrs["axis"] = "X"
        output_ds["time"].attrs["axis"] = "T"
        output_ds["time"].attrs["standard_name"] = "time"

        # Reorder the dimensions to follow NWP output standards.
        # - dimensions: ("time", ["ensemble"], "lev", "lat", "lon").
        # - vertical levels increasing with height.
        if self.is_ensemble_member:
            sorted_dims = ("time", "ens", "lev", "lat", "lon")
        else:
            sorted_dims = ("time", "lev", "lat", "lon")
        output_ds = output_ds.reindex(lev=output_ds.lev[::-1]).transpose(*sorted_dims)
        return output_ds

    def check(self):
        """Run a diagnostic check for the model state."""
        error_code = _speedy.check(self._state_cnt)
        if error_code < 0:
            raise RuntimeError(ERROR_CODES[error_code])


class SpeedyEns:
    """
    Ensemble of Speedy model instances.
    """

    def __init__(
        self,
        num_of_members,
        start_date=datetime(1982, 1, 1),
        end_date=datetime(1982, 1, 2),
        physics_mode="speedy",
    ):
        """
        Constructor. Initializes the ensemble of Speedy instances.

        Parameters
        ----------
        num_of_members: int
            Ensemble size.
        start_date: datetime
            Model start date.
        end_date: datetime
            Model end date.
        physics_mode: str
            Physics configuration shared by all ensemble members.
        """
        self.n_members = num_of_members
        self.members = [
            Speedy(
                start_date=start_date,
                end_date=end_date,
                physics_mode=physics_mode,
                member=mem_num,
            )
            for mem_num in range(num_of_members)
        ]
        self.current_date = self.members[
            0
        ].current_date  # Current date for the ensemble run.

    def __iter__(self):
        """Iterate over the ensemble members"""
        return iter(self.members)

    def __len__(self):
        return self.n_members

    @property
    def case_tag(self):
        return self.members[0].case_tag

    @property
    def output_tag(self):
        return self.members[0].output_tag

    def set_params(
        self, start_date=datetime(1982, 1, 1), end_date=datetime(1982, 1, 2)
    ):
        """
        Set the model's control parameters.

        See :meth:`~Speedy.set_params` for additional details.
        """
        for member in self:
            member.set_params(start_date=start_date, end_date=end_date)

        self.current_date = start_date

    def to_dataframe(self, variables=None):
        """
        Return an xarray DataFrame with the current model state.
        """
        member_dfs = []
        for member in self:
            member_dfs.append(member.to_dataframe(variables=variables))
        return xr.merge(member_dfs, join="outer", compat="no_conflicts")

    def run(self, callbacks=None):
        """
        Run each ensemble member between the start and the end dates

        Before running the model, the boundary conditions need to be set for each members
        (see the :meth:`set_bc` method).

        To initialize all the member with the default parameters, run:

        >>> speedy_ensemble = SpeedyEns()
        >>> for member in speedy_ensemble:
        >>>     member.set_bc()

        Parameters
        ----------
        callbacks: iterable of callback functions
            Sequence of callback functions to be call every `interval` time steps. The interval is specified in each
            callback instance.
        """
        if callbacks is None:
            callbacks = []

        end_date = self.members[0].end_date
        dt_step = timedelta(seconds=MODEL_CONFIG.timestep_seconds)

        state_cnts = np.zeros(self.n_members, dtype=np.int64)
        control_cnts = np.zeros(self.n_members, dtype=np.int64)
        for m, member in enumerate(self):
            state_cnts[m] = member._state_cnt  # noqa
            control_cnts[m] = member._control_cnt  # noqa

        while self.current_date < end_date:
            error_codes = _speedy.parallel_step(state_cnts, control_cnts)
            self.current_date += dt_step

            if (error_codes < 0).any():
                msg = ""
                for n, code in enumerate(error_codes):
                    msg += f"Member{n}: {ERROR_CODES[code]}\n"
                raise RuntimeError(msg)

            # Update current date in all members
            for member in self:
                member.current_date = self.current_date

            for callback in callbacks:
                callback(self)

    def get_current_step(self):
        """Return the current step in the simulation."""
        return self.members[0]["current_step"]
