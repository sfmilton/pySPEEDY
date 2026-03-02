from argparse import ArgumentParser
from os.path import splitext
import xarray as xr

from pyspeedy.pressure_levels import DEFAULT_PRESSURE_LEVELS_HPA, sigma_level_data_vars, sigma_to_pressure_dataset

p_levels = DEFAULT_PRESSURE_LEVELS_HPA


# Parse command line arguments
parser = ArgumentParser(description="Converts the given file from sigma level to pressure level")
parser.add_argument("filename", type=str, help="File name of input")
args = parser.parse_args()

with xr.open_dataset(args.filename) as dataset:
    if "ps" not in dataset.data_vars:
        print(f"Couldn't find surface pressure variable 'ps' in file {args.filename}")
        raise SystemExit
    sigma_vars = sigma_level_data_vars(dataset)
    if not sigma_vars:
        print(f"No sigma-level variables found in {args.filename}")
        raise SystemExit

    print("Processing the following sigma-level variables...")
    print(", ".join(sigma_vars))
    print(
        f"Processing {args.filename} from sigma levels\n"
        f"{', '.join(str(level) for level in dataset['lev'].values)}\n"
        f"to pressure levels\n"
        f"{', '.join(str(level) + ' hPa' for level in p_levels)}"
    )
    pressure_dataset = sigma_to_pressure_dataset(dataset)
    new_filename = f"{splitext(args.filename)[0]}_p.nc"
    print(f"Saving to {new_filename}")
    pressure_dataset.to_netcdf(new_filename)
