"""
Callbacks module
================

A callback is an object can perform a particular task at the end of the model time step.
The following callbacks are available in pySPEEDY:

.. autosummary::
    :toctree: ../generated/

    BaseCallback
    DiagnosticCheck
    RuntimeSummary
    ModelCheckpoint
    DailyDiagnosticsExporter
    XarrayExporter
"""

import copy
import glob
import os

from datetime import datetime

import numpy as np
import xarray as xr
from pyspeedy import (
    _speedy,  # noqa
    DEFAULT_OUTPUT_VARS,
    MODEL_STATE_DEF,
)
from pyspeedy.config import load_config
from pyspeedy.pressure_levels import sigma_level_data_vars, sigma_to_pressure_dataset
from pyspeedy.speedy import SpeedyEns, Speedy

DEFAULT_RUN_CONFIG = load_config().run
DEFAULT_MODEL_CONFIG = load_config().model
SEASON_MONTHS = {
    "DJF": ((-1, 12), (0, 1), (0, 2)),
    "MAM": ((0, 3), (0, 4), (0, 5)),
    "JJA": ((0, 6), (0, 7), (0, 8)),
    "SON": ((0, 9), (0, 10), (0, 11)),
}
DEFAULT_DAILY_TENDENCY_VARS = (
    "ttend_dyn_mean",
    "ttend_phy_mean",
    "ttend_cnv_mean",
    "ttend_lsc_mean",
    "ttend_sw_mean",
    "ttend_lw_mean",
    "ttend_pbl_mean",
    "ttend_hs_mean",
    "qtend_dyn_mean",
    "qtend_phy_mean",
    "qtend_cnv_mean",
    "qtend_lsc_mean",
    "qtend_pbl_mean",
    "utend_dyn_mean",
    "vtend_dyn_mean",
    "utend_phy_mean",
    "vtend_phy_mean",
    "utend_pbl_mean",
    "vtend_pbl_mean",
    "utend_gwd_mean",
    "vtend_gwd_mean",
    "utend_hs_mean",
    "vtend_hs_mean",
    "ustr_sfc_mean",
    "vstr_sfc_mean",
    "cloud_cover_mean",
    "stratiform_cloud_cover_mean",
    "total_cloud_top_pressure_mean",
    "conv_cloud_top_pressure_mean",
    "column_water_vapor_mean",
    "precip_mean",
    "evap_mean",
    "toa_sw_down_mean",
    "toa_sw_up_mean",
    "toa_sw_net_mean",
    "olr_mean",
    "surface_lh_flux_mean",
    "surface_sh_flux_mean",
    "surface_sw_down_mean",
    "surface_sw_up_mean",
    "surface_sw_net_mean",
    "surface_lw_down_mean",
    "surface_lw_up_mean",
    "surface_lw_net_mean",
    "soil_avail_water",
)


def _save_dataset(dataset, output_file_path, print_msg):
    vars_and_coords = list(dataset.variables.keys()) + list(dataset.coords.keys())
    encoding = dict()
    for var in vars_and_coords:
        if var in ("time", "ens"):
            continue

    print_msg(f"Saving model output at: {output_file_path}.")
    dataset.to_netcdf(output_file_path, encoding=encoding)


def _pressure_output_path(output_file_path):
    root, ext = os.path.splitext(output_file_path)
    return f"{root}_p{ext}"


def _append_time_slice(existing_dataset, current_dataset):
    if existing_dataset is None:
        return current_dataset.sortby("time")

    return current_dataset.sortby("time").combine_first(existing_dataset.sortby("time")).sortby("time")


def _load_dataset(dataset_path):
    with xr.open_dataset(dataset_path) as dataset:
        return dataset.load()


def _month_key(model_date):
    return model_date.year, model_date.month


def _season_key(model_date):
    if model_date.month == 12:
        return model_date.year + 1, "DJF"
    if model_date.month in (1, 2):
        return model_date.year, "DJF"
    if model_date.month in (3, 4, 5):
        return model_date.year, "MAM"
    if model_date.month in (6, 7, 8):
        return model_date.year, "JJA"
    return model_date.year, "SON"


def _season_month_keys(season_key):
    season_year, season = season_key
    return tuple((season_year + year_offset, month) for year_offset, month in SEASON_MONTHS[season])


def _time_value_to_datetime(value):
    return datetime.strptime(
        np.datetime_as_string(np.datetime64(value, "s"), unit="s"),
        "%Y-%m-%dT%H:%M:%S",
    )


def _dataset_time_bounds(dataset):
    time_values = np.sort(np.asarray(dataset["time"].values))
    return _time_value_to_datetime(time_values[0]), _time_value_to_datetime(time_values[-1])


def _monthly_filename(dataset, suffix, output_tag=None):
    year = int(dataset["time"].dt.year.isel(time=0).item())
    month = int(dataset["time"].dt.month.isel(time=0).item())
    start_day = int(dataset["time"].dt.day.min().item())
    end_day = int(dataset["time"].dt.day.max().item())
    file_name = f"{year:04d}-{month:02d}_d{start_day:02d}-d{end_day:02d}{suffix}"
    if output_tag is not None:
        file_name = f"{output_tag}_{file_name}"
    return file_name


def _seasonal_filename(dataset, suffix, output_tag=None):
    start_time, end_time = _dataset_time_bounds(dataset)
    season_year, season = _season_key(end_time)
    file_name = (
        f"{season_year:04d}-{season}_m{start_time.month:02d}d{start_time.day:02d}"
        f"-m{end_time.month:02d}d{end_time.day:02d}{suffix}"
    )
    if output_tag is not None:
        file_name = f"{output_tag}_{file_name}"
    return file_name


def _existing_monthly_output(output_dir, month_key, suffix, output_tag=None):
    prefix = f"{output_tag}_" if output_tag is not None else ""
    year, month = month_key
    pattern = os.path.join(
        output_dir,
        f"{prefix}{year:04d}-{month:02d}_d??-d??{suffix}",
    )
    matches = glob.glob(pattern)
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)


def _matching_monthly_outputs(output_dir, month_key, suffix, output_tag=None):
    prefix = f"{output_tag}_" if output_tag is not None else ""
    year, month = month_key
    pattern = os.path.join(
        output_dir,
        f"{prefix}{year:04d}-{month:02d}_d??-d??{suffix}",
    )
    return glob.glob(pattern)


def _existing_seasonal_output(output_dir, season_key, suffix, output_tag=None):
    prefix = f"{output_tag}_" if output_tag is not None else ""
    season_year, season = season_key
    pattern = os.path.join(
        output_dir,
        f"{prefix}{season_year:04d}-{season}_m??d??-m??d??{suffix}",
    )
    matches = glob.glob(pattern)
    if not matches:
        return None
    return max(matches, key=os.path.getmtime)


def _matching_seasonal_outputs(output_dir, season_key, suffix, output_tag=None):
    prefix = f"{output_tag}_" if output_tag is not None else ""
    season_year, season = season_key
    pattern = os.path.join(
        output_dir,
        f"{prefix}{season_year:04d}-{season}_m??d??-m??d??{suffix}",
    )
    return glob.glob(pattern)


def _pressure_source_variables(variables):
    needs_pressure_levels = any("lev" in (MODEL_STATE_DEF[var]["nc_dims"] or []) for var in variables)
    if not needs_pressure_levels:
        return None
    if "ps_grid" in variables:
        return tuple(variables)
    return tuple(variables) + ("ps_grid",)


def _pressure_companion_dataset(model_instance, model_df, variables, pressure_levels):
    pressure_source_variables = _pressure_source_variables(variables)
    if pressure_source_variables is None:
        return None

    if pressure_source_variables == tuple(variables):
        pressure_source_df = model_df
    else:
        pressure_source_df = model_instance.to_dataframe(variables=pressure_source_variables)

    pressure_df = sigma_to_pressure_dataset(
        pressure_source_df,
        pressure_levels_hpa=pressure_levels,
    )
    if "ps_grid" not in variables and "ps" in pressure_df.data_vars:
        pressure_df = pressure_df.drop_vars("ps")
    if not sigma_level_data_vars(pressure_df):
        return None
    return pressure_df


def _write_monthly_output(
    callback,
    model_instance,
    model_df,
    suffix,
    pressure_df=None,
):
    os.makedirs(callback.output_dir, exist_ok=True)

    month_key = _month_key(model_instance.current_date)
    output_tag = model_instance.output_tag if callback.include_output_tag else None

    if getattr(callback, "_active_month_key", None) != month_key:
        callback._active_month_key = month_key
        callback._monthly_dataset = None
        callback._monthly_pressure_dataset = None
        callback._monthly_output_path = None
        callback._monthly_pressure_output_path = None

    existing_output_path = _existing_monthly_output(
        callback.output_dir,
        month_key,
        suffix,
        output_tag=output_tag,
    )
    if existing_output_path is not None:
        # Reload the current on-disk monthly file before each append so
        # multiple callbacks can safely merge variables into the same output.
        callback._monthly_dataset = _load_dataset(existing_output_path)
        callback._monthly_output_path = existing_output_path
        pressure_output_path = _pressure_output_path(existing_output_path)
        if pressure_df is not None and os.path.exists(pressure_output_path):
            callback._monthly_pressure_dataset = _load_dataset(pressure_output_path)
            callback._monthly_pressure_output_path = pressure_output_path

    callback._monthly_dataset = _append_time_slice(callback._monthly_dataset, model_df)
    output_file_name = _monthly_filename(
        callback._monthly_dataset,
        suffix,
        output_tag=output_tag,
    )
    output_file_path = os.path.join(callback.output_dir, output_file_name)
    previous_output_path = callback._monthly_output_path
    _save_dataset(callback._monthly_dataset, output_file_path, callback.print_msg)
    if (
        previous_output_path is not None
        and previous_output_path != output_file_path
        and os.path.exists(previous_output_path)
    ):
        os.remove(previous_output_path)
    callback._monthly_output_path = output_file_path

    if pressure_df is None:
        return

    callback._monthly_pressure_dataset = _append_time_slice(
        callback._monthly_pressure_dataset,
        pressure_df,
    )
    pressure_output_path = _pressure_output_path(output_file_path)
    previous_pressure_output_path = callback._monthly_pressure_output_path
    _save_dataset(
        callback._monthly_pressure_dataset,
        pressure_output_path,
        callback.print_msg,
    )
    if (
        previous_pressure_output_path is not None
        and previous_pressure_output_path != pressure_output_path
        and os.path.exists(previous_pressure_output_path)
    ):
        os.remove(previous_pressure_output_path)
    callback._monthly_pressure_output_path = pressure_output_path


def _mean_dataset(dataset, aggregation):
    start_time, end_time = _dataset_time_bounds(dataset)
    mean_dataset = dataset.mean(dim="time", skipna=True, keep_attrs=True).expand_dims(
        time=[np.datetime64(end_time)]
    )
    mean_dataset.attrs = dict(dataset.attrs)
    mean_dataset.attrs["aggregation"] = aggregation
    mean_dataset.attrs["period_start"] = start_time.strftime("%Y-%m-%d")
    mean_dataset.attrs["period_end"] = end_time.strftime("%Y-%m-%d")
    mean_dataset.attrs["sample_count"] = int(dataset.sizes["time"])
    if "time" in dataset.coords:
        mean_dataset["time"].attrs = dict(dataset["time"].attrs)
    for var_name in mean_dataset.data_vars:
        mean_dataset[var_name].attrs = dict(dataset[var_name].attrs)
    return mean_dataset


def _write_derived_mean_output(
    dataset,
    output_dir,
    output_file_name,
    print_msg,
    previous_output_path=None,
    stale_output_paths=None,
):
    if dataset is None:
        return None

    output_file_path = os.path.join(output_dir, output_file_name)
    _save_dataset(dataset, output_file_path, print_msg)
    if (
        previous_output_path is not None
        and previous_output_path != output_file_path
        and os.path.exists(previous_output_path)
    ):
        os.remove(previous_output_path)
    if stale_output_paths is not None:
        for stale_output_path in stale_output_paths:
            if stale_output_path != output_file_path and os.path.exists(stale_output_path):
                os.remove(stale_output_path)
    return output_file_path


def _load_seasonal_daily_dataset(callback, season_key, output_tag, pressure=False):
    seasonal_dataset = None
    active_month_key = getattr(callback, "_active_month_key", None)
    active_month_dataset = getattr(
        callback,
        "_monthly_pressure_dataset" if pressure else "_monthly_dataset",
        None,
    )

    for month_key in _season_month_keys(season_key):
        if month_key == active_month_key and active_month_dataset is not None:
            month_dataset = active_month_dataset
        else:
            month_output_path = _existing_monthly_output(
                callback.output_dir,
                month_key,
                ".nc",
                output_tag=output_tag,
            )
            if month_output_path is None:
                continue
            if pressure:
                month_output_path = _pressure_output_path(month_output_path)
                if not os.path.exists(month_output_path):
                    continue
            month_dataset = _load_dataset(month_output_path)
        seasonal_dataset = _append_time_slice(seasonal_dataset, month_dataset)

    return seasonal_dataset


def _write_period_mean_outputs(callback, output_tag):
    if callback._monthly_dataset is not None and callback.write_monthly_means:
        existing_monthly_mean_paths = _matching_monthly_outputs(
            callback.output_dir,
            callback._active_month_key,
            "_monthly_mean.nc",
            output_tag=output_tag,
        )
        previous_monthly_mean_path = _existing_monthly_output(
            callback.output_dir,
            callback._active_month_key,
            "_monthly_mean.nc",
            output_tag=output_tag,
        )
        previous_monthly_pressure_mean_path = _existing_monthly_output(
            callback.output_dir,
            callback._active_month_key,
            "_monthly_mean_p.nc",
            output_tag=output_tag,
        )
        monthly_mean = _mean_dataset(callback._monthly_dataset, "monthly_mean")
        monthly_mean_path = _write_derived_mean_output(
            monthly_mean,
            callback.output_dir,
            _monthly_filename(callback._monthly_dataset, "_monthly_mean.nc", output_tag=output_tag),
            callback.print_msg,
            previous_output_path=previous_monthly_mean_path,
            stale_output_paths=existing_monthly_mean_paths,
        )
        if callback._monthly_pressure_dataset is not None:
            existing_monthly_pressure_mean_paths = _matching_monthly_outputs(
                callback.output_dir,
                callback._active_month_key,
                "_monthly_mean_p.nc",
                output_tag=output_tag,
            )
            monthly_pressure_mean = _mean_dataset(callback._monthly_pressure_dataset, "monthly_mean")
            _write_derived_mean_output(
                monthly_pressure_mean,
                callback.output_dir,
                os.path.basename(_pressure_output_path(monthly_mean_path)),
                callback.print_msg,
                previous_output_path=previous_monthly_pressure_mean_path,
                stale_output_paths=existing_monthly_pressure_mean_paths,
            )

    if not callback.write_seasonal_means:
        return

    season_key = _season_key(callback.current_output_date)
    seasonal_daily_dataset = _load_seasonal_daily_dataset(callback, season_key, output_tag, pressure=False)
    if seasonal_daily_dataset is None:
        return

    existing_seasonal_mean_paths = _matching_seasonal_outputs(
        callback.output_dir,
        season_key,
        "_seasonal_mean.nc",
        output_tag=output_tag,
    )
    previous_seasonal_mean_path = _existing_seasonal_output(
        callback.output_dir,
        season_key,
        "_seasonal_mean.nc",
        output_tag=output_tag,
    )
    previous_seasonal_pressure_mean_path = _existing_seasonal_output(
        callback.output_dir,
        season_key,
        "_seasonal_mean_p.nc",
        output_tag=output_tag,
    )
    seasonal_mean = _mean_dataset(seasonal_daily_dataset, "seasonal_mean")
    seasonal_mean_path = _write_derived_mean_output(
        seasonal_mean,
        callback.output_dir,
        _seasonal_filename(seasonal_daily_dataset, "_seasonal_mean.nc", output_tag=output_tag),
        callback.print_msg,
        previous_output_path=previous_seasonal_mean_path,
        stale_output_paths=existing_seasonal_mean_paths,
    )

    seasonal_pressure_daily_dataset = _load_seasonal_daily_dataset(
        callback,
        season_key,
        output_tag,
        pressure=True,
    )
    if seasonal_pressure_daily_dataset is None:
        return

    seasonal_pressure_mean = _mean_dataset(seasonal_pressure_daily_dataset, "seasonal_mean")
    existing_seasonal_pressure_mean_paths = _matching_seasonal_outputs(
        callback.output_dir,
        season_key,
        "_seasonal_mean_p.nc",
        output_tag=output_tag,
    )
    _write_derived_mean_output(
        seasonal_pressure_mean,
        callback.output_dir,
        os.path.basename(_pressure_output_path(seasonal_mean_path)),
        callback.print_msg,
        previous_output_path=previous_seasonal_pressure_mean_path,
        stale_output_paths=existing_seasonal_pressure_mean_paths,
    )


def _area_weighted_global_mean(field_2d, lat):
    weights = np.cos(np.deg2rad(lat))[None, :]
    valid = np.isfinite(field_2d)
    weighted_sum = np.sum(np.where(valid, field_2d, 0.0) * weights)
    total_weight = np.sum(weights * valid)
    if total_weight == 0.0:
        return np.nan
    return float(weighted_sum / total_weight)


class BaseCallback:
    """
    Base callback class.
    """

    def __init__(self, *args, **kwargs):
        """
        Constructor.

        Parameters
        ----------
        interval: int
            Interval, in time steps, for which the callback should be applied.
        verbose: bool
            If true, print debug and progress messages.
        spinup_date: datetime or None:
            End date of the spinup period. During the spinup the callback is ignored.
        """
        self.verbose = kwargs.pop("verbose", False)
        self.interval = kwargs.pop("interval", 1)
        self.spinup_date = kwargs.pop("spinup_date", None)

    def skip_flag(self, model_instance):
        """
        Return True when the callback execution is skipped for this time step if
        - model_date < spinup_date
        - current_step % interval !=0
        """
        if self.spinup_date is not None:
            if model_instance.current_date < self.spinup_date:
                # Do not save during spinup time
                return True

        return model_instance.get_current_step() % self.interval != 0

    def print_msg(self, msg):
        """Print debug message if `verbose` was set to True."""
        if self.verbose:
            print(msg)

    def copy(self):
        """Create a copy of the instance."""
        return copy.deepcopy(self)

    def __call__(self, model_instance):
        """Object call."""
        pass


class DiagnosticCheck(BaseCallback):
    """
    Callback used to check on that the prognostic variables are inside reasonable ranges.
    """

    def __init__(self, interval=None):
        """
        Constructor.

        Parameters
        ----------
        interval: int
            Interval, in time steps, for which the diagnostic checkes are run.
        """
        if interval is None:
            interval = DEFAULT_MODEL_CONFIG.nsteps
        super().__init__(interval=interval)

    def __call__(self, model_instance):
        """
        Object call.

        Run diagnostic check if needed.
        """
        if self.skip_flag(model_instance):
            # Only run tests every `interval` steps.
            return

        if isinstance(model_instance, Speedy):
            # If it is a single run, we convert it to a list to make iterable.
            model_instance = [model_instance]

        for _member in model_instance:
            # If the test fails, it will raise a runtime exception.
            _member.check()


class RuntimeSummary(BaseCallback):
    """
    Print a compact runtime health summary for long integrations.

    The summary includes:
    - current model date
    - global-mean lowest-model-level air temperature
    - global-mean daily precipitation
    - global-mean daily net TOA radiation
    - max absolute zonal wind
    - minimum specific humidity
    """

    def __init__(self, interval=None, verbose=True, spinup_date=None):
        if interval is None:
            interval = DEFAULT_RUN_CONFIG.diag_interval
        if interval % DEFAULT_MODEL_CONFIG.nsteps != 0:
            raise ValueError(
                "RuntimeSummary requires an interval that is an integer multiple "
                f"of one model day ({DEFAULT_MODEL_CONFIG.nsteps} timesteps)."
            )
        super().__init__(verbose=verbose, interval=interval, spinup_date=spinup_date)

    def __call__(self, model_instance):
        if self.skip_flag(model_instance):
            return

        model_instance.spectral2grid()

        lat = np.asarray(model_instance["lat"], dtype=float)
        lev = np.asarray(model_instance["lev"], dtype=float)
        lowest_level_idx = int(np.argmax(lev))

        t_grid = np.asarray(model_instance["t_grid"], dtype=float)
        q_grid = np.asarray(model_instance["q_grid"], dtype=float)
        u_grid = np.asarray(model_instance["u_grid"], dtype=float)
        precip_day = np.asarray(model_instance["precip_mean"], dtype=float)
        toa_net_day = np.asarray(model_instance["toa_sw_net_mean"], dtype=float) - np.asarray(
            model_instance["olr_mean"],
            dtype=float,
        )

        global_mean_t_low = _area_weighted_global_mean(t_grid[:, :, lowest_level_idx], lat)
        global_mean_precip = _area_weighted_global_mean(precip_day, lat)
        global_mean_toa_net = _area_weighted_global_mean(toa_net_day, lat)
        max_abs_u = float(np.nanmax(np.abs(u_grid)))
        min_q = float(np.nanmin(q_grid) * 1000.0)

        summary = (
            f"{model_instance.current_date:%Y-%m-%d}  "
            f"Tlow={global_mean_t_low:.1f} K  "
            f"P={global_mean_precip:.2f} mm/day  "
            f"TOA_net={global_mean_toa_net:.1f} W/m^2  "
            f"max|u|={max_abs_u:.1f} m/s  "
            f"min(q)={min_q:.2f} g/kg"
        )
        self.print_msg(summary)


class ModelCheckpoint(BaseCallback):
    """
    Callback used save selected variables of the model state at the current time step.
    Note that this is not the full model state, only a small subset
    of it. We will refer to this saved state as "Checkpoint".

    Each instance of this callback will keep internally a time series of the checkpoints saved during the the model run
    in an xarray DataFrame following CF conventions.

    Notes
    -----
    The variables are saved in the lat/lon grid space (not the spectral domain).
    Spectral variables are not supported by this callback.
    """

    def __init__(
        self,
        interval=None,
        verbose=False,
        spinup_date=None,
        variables=None,
        output_dir="./",
    ):
        """
        Parameters
        ----------
        interval: int
            Interval, in time steps, for which the model variables are saved.
        verbose: bool
            If true, print debug and progress messages.
        spinup_date: datetime or None:
            End date of the spinup period. During the spinup the output files are not saved.
        variables: list, tuple
            List of variables to save
        output_dir: str
            Path to folder where the output files are stored.
        interval: int
            History interval in timesteps every which the output files are saved.
        spinup_date: datetime or None
            Model spinup date. From `start_date` to `spinup_date` the callbacks functions are not called.
        """
        if variables is None:
            variables = DEFAULT_OUTPUT_VARS
        if interval is None:
            interval = DEFAULT_RUN_CONFIG.history_interval
        self.variables = variables
        self.output_dir = output_dir
        self.history_interval = interval
        super().__init__(verbose=verbose, interval=interval, spinup_date=spinup_date)
        self.dataframe = None

    def __call__(self, model_instance):
        """
        Object call.

        Export the model data to file using xarray.
        """
        if self.skip_flag(model_instance):
            # Only save files at every history_interval steps.
            return

        model_df = model_instance.to_dataframe(variables=self.variables)
        if self.dataframe is None:
            self.dataframe = model_df
        else:
            self.dataframe = xr.merge(
                (self.dataframe, model_df), join="outer", compat="no_conflicts"
            )


class XarrayExporter(BaseCallback):
    """
    Callback used to create an xarray dataset with selected variables at the current model time step.

    The variables are saved in the lat/lon grid space (not the spectral domain).
    Spectral variables are not supported by this callback.

    For ensemble runs, each output file is saved in a different subdirectory (named "member###").
    This ensure that concurrent calls to the callback from different ensemble members do not interfere with each other.
    """

    def __init__(
        self,
        interval=None,
        verbose=False,
        spinup_date=None,
        variables=None,
        output_dir="./",
        filename_fmt=None,
        export_pressure_levels=True,
        pressure_levels=None,
    ):
        """
        Parameters
        ----------
        interval: int
            Interval, in time steps, for which the model variables are saved.
        verbose: bool
            If true, print debug and progress messages.
        spinup_date: datetime or None:
            End date of the spinup period. During the spinup the model checkpoints are not saved.
        variables: list, tuple
            List of variables to save
        output_dir: str
            Path to folder where the output files are stored.
        filename_fmt: str
            Optional format string used to generate one file per callback invocation. When omitted, outputs are
            accumulated into one NetCDF per calendar month named with the covered day range, for example
            ``..._1982-01_d02-d31.nc``.
        export_pressure_levels: bool
            If true, write a companion NetCDF file with sigma-level variables interpolated to pressure levels.
        pressure_levels: array-like or None
            Pressure levels in hPa used for the companion export. If None, use the package defaults.
        interval: int
            History interval in timesteps every which the output files are saved.
        spinup_date: datetime or None
            Model spinup date. From `start_date` to `spinup_date` the callbacks functions are not called.
        """
        if variables is None:
            variables = DEFAULT_OUTPUT_VARS
        if interval is None:
            interval = DEFAULT_RUN_CONFIG.history_interval
        self.include_output_tag = filename_fmt is None
        self.use_monthly_files = filename_fmt is None
        if filename_fmt is None:
            filename_fmt = "%Y-%m-%d_%H%M.nc"
        self.variables = variables
        self.output_dir = output_dir
        self.filename_fmt = filename_fmt
        self.export_pressure_levels = export_pressure_levels
        self.pressure_levels = pressure_levels
        self.history_interval = interval
        self._active_month_key = None
        self._monthly_dataset = None
        self._monthly_pressure_dataset = None
        self._monthly_output_path = None
        self._monthly_pressure_output_path = None
        super().__init__(verbose=verbose, interval=interval, spinup_date=spinup_date)

    def __call__(self, model_instance):
        """
        Object call.

        Export the model data to file using xarray.
        """
        if self.skip_flag(model_instance):
            # Only save files at every history_interval steps.
            return

        model_df = model_instance.to_dataframe(variables=self.variables)

        pressure_df = None
        if self.export_pressure_levels:
            pressure_df = _pressure_companion_dataset(
                model_instance,
                model_df,
                self.variables,
                self.pressure_levels,
            )

        if self.use_monthly_files:
            _write_monthly_output(
                self,
                model_instance,
                model_df,
                ".nc",
                pressure_df=pressure_df,
            )
            return

        file_name = model_instance.current_date.strftime(self.filename_fmt)
        if self.include_output_tag:
            file_name = f"{model_instance.output_tag}_{file_name}"
        os.makedirs(self.output_dir, exist_ok=True)
        output_file_path = os.path.join(self.output_dir, file_name)
        _save_dataset(model_df, output_file_path, self.print_msg)

        if pressure_df is None:
            return
        _save_dataset(pressure_df, _pressure_output_path(output_file_path), self.print_msg)


class DailyDiagnosticsExporter(BaseCallback):
    """
    Callback used to export daily-mean tendencies together with selected
    hydro, radiation, and surface-flux diagnostics.

    By default, the daily diagnostics are merged into the same monthly sigma-level
    and pressure-level files written by ``XarrayExporter``. The callback can also
    write derived monthly-mean and seasonal-mean files from those unified daily outputs.
    """

    def __init__(
        self,
        interval=None,
        verbose=False,
        spinup_date=None,
        variables=None,
        output_dir="./",
        filename_fmt=None,
        export_pressure_levels=True,
        pressure_levels=None,
        write_monthly_means=True,
        write_seasonal_means=True,
    ):
        if variables is None:
            variables = DEFAULT_DAILY_TENDENCY_VARS
        if interval is None:
            interval = DEFAULT_MODEL_CONFIG.nsteps
        elif interval != DEFAULT_MODEL_CONFIG.nsteps:
            raise ValueError(
                "DailyDiagnosticsExporter requires interval to equal one model day "
                f"({DEFAULT_MODEL_CONFIG.nsteps} timesteps)."
            )
        self.include_output_tag = filename_fmt is None
        self.use_monthly_files = filename_fmt is None
        if filename_fmt is None:
            filename_fmt = "%Y-%m-%d_%H%M.nc"
        self.variables = tuple(variables)
        self.output_dir = output_dir
        self.filename_fmt = filename_fmt
        self.export_pressure_levels = export_pressure_levels
        self.pressure_levels = pressure_levels
        self.write_monthly_means = write_monthly_means and self.use_monthly_files
        self.write_seasonal_means = write_seasonal_means and self.use_monthly_files
        self._active_month_key = None
        self._monthly_dataset = None
        self._monthly_pressure_dataset = None
        self._monthly_output_path = None
        self._monthly_pressure_output_path = None
        self.current_output_date = None
        super().__init__(verbose=verbose, interval=interval, spinup_date=spinup_date)

    def __call__(self, model_instance):
        if self.skip_flag(model_instance):
            return

        self.current_output_date = model_instance.current_date
        model_df = model_instance.to_dataframe(variables=self.variables)
        pressure_df = None
        if self.export_pressure_levels:
            pressure_df = _pressure_companion_dataset(
                model_instance,
                model_df,
                self.variables,
                self.pressure_levels,
            )

        if self.use_monthly_files:
            _write_monthly_output(
                self,
                model_instance,
                model_df,
                ".nc",
                pressure_df=pressure_df,
            )
            if self.write_monthly_means or self.write_seasonal_means:
                output_tag = model_instance.output_tag if self.include_output_tag else None
                _write_period_mean_outputs(self, output_tag)
            return

        file_name = model_instance.current_date.strftime(self.filename_fmt)
        if self.include_output_tag:
            file_name = f"{model_instance.output_tag}_{file_name}"
        os.makedirs(self.output_dir, exist_ok=True)
        output_file_path = os.path.join(self.output_dir, file_name)
        _save_dataset(model_df, output_file_path, self.print_msg)

        if pressure_df is None:
            return
        _save_dataset(pressure_df, _pressure_output_path(output_file_path), self.print_msg)


# Backwards-compatible alias for older scripts.
DailyTemperatureTendencyExporter = DailyDiagnosticsExporter
