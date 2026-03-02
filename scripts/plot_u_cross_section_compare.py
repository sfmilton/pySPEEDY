#!/usr/bin/env python
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import xarray as xr


def _parse_args():
    parser = argparse.ArgumentParser(
        description="Plot a zonal-mean zonal-wind cross-section comparison between two monthly-mean files."
    )
    parser.add_argument("gwd_file", help="Monthly-mean pressure-level NetCDF with GWD enabled.")
    parser.add_argument("nogwd_file", help="Monthly-mean pressure-level NetCDF with GWD disabled.")
    parser.add_argument(
        "--output",
        default="data/u_cross_section_gwd_vs_nogwd_january_1981.png",
        help="Path to the output PNG.",
    )
    return parser.parse_args()


def _open_zonal_mean_u(path):
    ds = xr.open_dataset(path)
    if "u" not in ds:
        raise KeyError(f"{path} does not contain variable 'u'.")
    u = ds["u"]
    if "lon" not in u.dims or "lat" not in u.dims or "lev" not in u.dims:
        raise ValueError(f"{path} variable 'u' must have lon/lat/lev dimensions.")
    if "time" in u.dims:
        u = u.isel(time=-1)
    return u.mean(dim="lon")


def _rounded_levels(max_abs, step):
    upper = step * max(1, int(np.ceil(max_abs / step)))
    return np.arange(-upper, upper + step, step)


def _difference_levels(max_abs):
    if max_abs <= 0.1:
        step = 0.02
    elif max_abs <= 0.25:
        step = 0.05
    elif max_abs <= 0.5:
        step = 0.1
    elif max_abs <= 1.0:
        step = 0.2
    else:
        step = 0.5
    return _rounded_levels(max_abs, step)


def main():
    args = _parse_args()

    gwd = _open_zonal_mean_u(args.gwd_file)
    nogwd = _open_zonal_mean_u(args.nogwd_file)
    diff = gwd - nogwd

    lat = gwd["lat"].values
    lev = gwd["lev"].values

    main_levels = _rounded_levels(
        max(float(np.nanmax(np.abs(gwd.values))), float(np.nanmax(np.abs(nogwd.values)))),
        5.0,
    )
    diff_levels = _difference_levels(float(np.nanmax(np.abs(diff.values))))
    if diff_levels.size < 3:
        diff_levels = np.array([-1.0, 0.0, 1.0])

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True, constrained_layout=True)

    panels = (
        (gwd, "GWD On", "RdBu_r", main_levels),
        (nogwd, "GWD Off", "RdBu_r", main_levels),
        (diff, "Difference (On - Off)", "BrBG", diff_levels),
    )

    for ax, (field, title, cmap, levels) in zip(axes, panels):
        contour = ax.contourf(lat, lev, field.values, levels=levels, cmap=cmap, extend="both")
        ax.contour(lat, lev, field.values, levels=levels[::2], colors="k", linewidths=0.4, alpha=0.5)
        ax.set_title(title)
        ax.set_xlabel("Latitude")
        ax.set_xticks(np.arange(-90, 91, 30))
        ax.invert_yaxis()
        fig.colorbar(contour, ax=ax, pad=0.02, shrink=0.9, label="m/s")

    axes[0].set_ylabel("Pressure (hPa)")
    fig.suptitle("January 1981 Zonal-Mean Zonal Wind", fontsize=14)

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    print(output_path)


if __name__ == "__main__":
    main()
