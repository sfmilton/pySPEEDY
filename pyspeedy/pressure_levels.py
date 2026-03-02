"""Utilities for converting sigma-level output to pressure levels."""

from __future__ import annotations

import numpy as np
import xarray as xr


DEFAULT_PRESSURE_LEVELS_HPA = np.array(
    [25.0, 100.0, 200.0, 300.0, 500.0, 700.0, 850.0, 950.0],
    dtype=np.float32,
)


def sigma_level_data_vars(dataset: xr.Dataset) -> tuple[str, ...]:
    """Return exported variables that depend on the sigma-level coordinate."""
    return tuple(name for name, data_array in dataset.data_vars.items() if "lev" in data_array.dims)


def _interp_profile_to_pressure(values, sigma_levels, surface_pressure_pa, pressure_levels_hpa):
    source_pressures_hpa = np.asarray(surface_pressure_pa, dtype=np.float64) * np.asarray(
        sigma_levels,
        dtype=np.float64,
    ) / 100.0
    source_values = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(source_pressures_hpa) & np.isfinite(source_values)

    if np.count_nonzero(valid) < 2:
        return np.full(pressure_levels_hpa.shape, np.nan, dtype=np.float32)

    source_pressures_hpa = source_pressures_hpa[valid]
    source_values = source_values[valid]
    order = np.argsort(source_pressures_hpa)
    interpolated = np.interp(
        pressure_levels_hpa,
        source_pressures_hpa[order],
        source_values[order],
        left=np.nan,
        right=np.nan,
    )
    return interpolated.astype(np.float32, copy=False)


def sigma_to_pressure_dataset(
    dataset: xr.Dataset,
    pressure_levels_hpa=None,
    surface_pressure_var: str = "ps",
) -> xr.Dataset:
    """Convert all sigma-level data variables in a dataset to pressure levels."""
    if pressure_levels_hpa is None:
        pressure_levels_hpa = DEFAULT_PRESSURE_LEVELS_HPA
    pressure_levels_hpa = np.asarray(pressure_levels_hpa, dtype=np.float32)

    sigma_vars = sigma_level_data_vars(dataset)
    if not sigma_vars:
        return dataset.copy(deep=True)

    if "lev" not in dataset.coords:
        raise ValueError("Dataset does not contain a 'lev' coordinate.")
    if surface_pressure_var not in dataset.data_vars:
        raise ValueError(
            f"Dataset does not contain the surface pressure variable '{surface_pressure_var}'."
        )

    pressure_dataset = dataset.drop_vars(list(sigma_vars)).copy(deep=True)
    pressure_dataset = pressure_dataset.assign_coords(lev=("lev", pressure_levels_hpa))
    pressure_dataset.attrs = dataset.attrs.copy()
    sigma_levels = xr.DataArray(
        np.asarray(dataset["lev"].values, dtype=np.float64),
        dims=("lev",),
        coords={"lev": dataset["lev"]},
    )
    surface_pressure = dataset[surface_pressure_var]

    for var_name in sigma_vars:
        data_array = dataset[var_name]
        converted = xr.apply_ufunc(
            _interp_profile_to_pressure,
            data_array,
            sigma_levels,
            surface_pressure,
            input_core_dims=[["lev"], ["lev"], []],
            output_core_dims=[["lev"]],
            vectorize=True,
            dask="parallelized",
            output_dtypes=[np.float32],
            dask_gufunc_kwargs={"output_sizes": {"lev": pressure_levels_hpa.size}},
            kwargs={"pressure_levels_hpa": pressure_levels_hpa},
        ).transpose(*data_array.dims)
        converted = converted.assign_coords(lev=pressure_levels_hpa)
        converted.attrs = data_array.attrs.copy()
        pressure_dataset[var_name] = converted

    lev_attrs = dict(dataset["lev"].attrs)
    lev_attrs.update(
        {
            "long_name": "air_pressure",
            "standard_name": "air_pressure",
            "units": "hPa",
            "positive": "down",
            "axis": "Z",
        }
    )
    pressure_dataset["lev"].attrs = lev_attrs
    return pressure_dataset
