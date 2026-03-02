from argparse import ArgumentParser

import numpy as np
import xarray as xr


def _format_stat(value):
    if value is None or not np.isfinite(value):
        return "nan"
    return f"{value:.4g}"


parser = ArgumentParser(description="Summarize tendency NetCDF contents and magnitudes.")
parser.add_argument("filename", type=str, help="Path to tendency NetCDF file")
parser.add_argument(
    "--variables",
    nargs="*",
    default=None,
    help="Optional subset of variables to summarize",
)
args = parser.parse_args()

with xr.open_dataset(args.filename) as dataset:
    print(f"File: {args.filename}")
    print(f"Dimensions: {dict(dataset.sizes)}")
    if "time" in dataset.coords and dataset.sizes.get("time", 0) > 0:
        start_time = str(dataset["time"].values[0])
        end_time = str(dataset["time"].values[-1])
        print(f"Time range: {start_time} -> {end_time}")

    variable_names = list(dataset.data_vars)
    if args.variables is not None:
        requested = list(args.variables)
        missing = [name for name in requested if name not in dataset.data_vars]
        if missing:
            print(f"Missing variables: {', '.join(missing)}")
        variable_names = [name for name in requested if name in dataset.data_vars]

    for variable_name in sorted(variable_names):
        data_array = dataset[variable_name]
        values = np.asarray(data_array.values, dtype=np.float64)
        finite = np.isfinite(values)
        finite_count = int(np.count_nonzero(finite))
        total_count = int(values.size)
        units = data_array.attrs.get("units", "-")

        if finite_count == 0:
            min_value = max_value = max_abs = None
        else:
            finite_values = values[finite]
            min_value = float(np.min(finite_values))
            max_value = float(np.max(finite_values))
            max_abs = float(np.max(np.abs(finite_values)))

        print(
            f"{variable_name}: dims={data_array.dims}, units={units}, "
            f"finite={finite_count}/{total_count}, "
            f"min={_format_stat(min_value)}, max={_format_stat(max_value)}, "
            f"maxabs={_format_stat(max_abs)}"
        )
