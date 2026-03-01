import os
import tempfile
from datetime import datetime, timedelta

import pytest
import xarray as xr

from pyspeedy.callbacks import XarrayExporter
from pyspeedy.config import load_config
from pyspeedy.speedy import Speedy, SpeedyEns

start_dates = (
    # Run twice the same date to check if the globals variables in the library are modified.
    (datetime(1982, 1, 1), datetime(1982, 1, 2)),
    (datetime(1982, 1, 1), datetime(1982, 1, 2)),
    # Run for several days
    (datetime(1982, 1, 1), datetime(1982, 1, 4)),
)

export_variables = (
    ["u_grid", "v_grid"],
    ["t_grid", "q_grid"],
    ["phi_grid", "ps_grid"],
    ["precnv", "precls"],
)


@pytest.mark.parametrize("start_date, end_date", start_dates)
def test_speedy_run(start_date, end_date):
    """Run speedy and verify the exported dataset matches the configured grid."""

    config = load_config()

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = f"{model.output_tag}_{end_date.strftime('%Y-%m-%d_%H%M.nc')}"
        model.set_bc()
        model.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir)])

        model_file = os.path.join(tmp_work_dir, file_name)
        model_ds = xr.open_dataset(model_file)
        assert model_ds.sizes["lon"] == config.model.ix
        assert model_ds.sizes["lat"] == config.model.il
        assert model_ds.sizes["lev"] == config.model.kx
        assert model_ds.sizes["time"] == 1
        for data_var in model_ds.data_vars:
            assert model_ds[data_var].notnull().all().item()


def test_speedy_concurrent():
    """Run serially 2 speedy models for 4 days, alternating each model run every one day of simulation."""
    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 4)
    ndays = 3

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        tmp_work_dir1 = os.path.join(tmp_work_dir, "run1")
        tmp_work_dir2 = os.path.join(tmp_work_dir, "run2")

        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = f"{model.output_tag}_{end_date.strftime('%Y-%m-%d_%H%M.nc')}"
        model.set_bc()

        # Create another speedy instance
        model2 = Speedy(start_date=start_date, end_date=end_date)
        model2.set_bc()

        for day in range(ndays):
            model.start_date = start_date + timedelta(days=day)
            model.end_date = start_date + timedelta(days=day + 1)
            model.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir1)])

            model2.start_date = start_date + timedelta(days=day)
            model2.end_date = start_date + timedelta(days=day + 1)
            model2.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir2)])

        model_file = os.path.join(tmp_work_dir1, file_name)
        model_ds = xr.open_dataset(model_file)

        model_file = os.path.join(tmp_work_dir2, file_name)
        model2_ds = xr.open_dataset(model_file)
        xr.testing.assert_allclose(model_ds, model2_ds, rtol=1e-06, atol=0)


def test_ens_speedy():
    """Test the SpeedyEns class."""
    num_of_members = 3
    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)
    model_ens = SpeedyEns(num_of_members, start_date=start_date, end_date=end_date)
    file_name = f"{model_ens.output_tag}_{end_date.strftime('%Y-%m-%d_%H%M.nc')}"
    for member in model_ens:
        member.set_bc()
    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model_ens.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir)])

        model_ens_ds = xr.open_dataset(os.path.join(tmp_work_dir, file_name))
        for m, member in enumerate(model_ens):
            xr.testing.assert_allclose(
                member.to_dataframe().squeeze(dim="ens", drop=True),
                model_ens_ds.sel(ens=m).drop_vars("ens"),
                rtol=1e-06,
                atol=0,
            )


def test_exceptions():
    """Test that certain exceptions are raised."""
    model = Speedy(start_date=datetime(1982, 1, 1), end_date=datetime(1982, 1, 2))
    model.set_bc()
    model.run()

    # Force a failure in the diagnostic check
    t = model["t"]
    t[:] = 0
    model["t"] = t
    with pytest.raises(RuntimeError):
        model.check()

    with pytest.raises(ValueError):
        Speedy(
            start_date=datetime(1982, 1, 1),
            end_date=datetime(1982, 1, 2),
            physics_mode="not-a-mode",
        )


def test_held_suarez_run():
    """Run the Held-Suarez mode and verify the exported dataset is finite."""

    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(
            start_date=start_date,
            end_date=end_date,
            physics_mode="held_suarez",
        )
        file_name = f"{model.output_tag}_{end_date.strftime('%Y-%m-%d_%H%M.nc')}"
        model.set_bc()
        model.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir)])

        model_file = os.path.join(tmp_work_dir, file_name)
        model_ds = xr.open_dataset(model_file)
        assert model_ds["q"].notnull().all().item()
        assert float(abs(model_ds["q"]).max().item()) == pytest.approx(0.0)


@pytest.mark.parametrize("variables", export_variables)
def test_speedy_variable_export(variables):
    """Run speedy and compare the output with a reference."""

    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = f"{model.output_tag}_{end_date.strftime('%Y-%m-%d_%H%M.nc')}"
        model.set_bc()

        exporter = XarrayExporter(output_dir=tmp_work_dir, variables=variables)
        model.run(callbacks=[exporter])

        model_file = os.path.join(tmp_work_dir, file_name)
        model_ds = xr.open_dataset(model_file)

        assert set((v.replace("_grid", "") for v in variables)) == set(model_ds.keys())
