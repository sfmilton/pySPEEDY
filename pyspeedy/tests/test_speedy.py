import os
import subprocess
import sys
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from pyspeedy.callbacks import (
    DailyDiagnosticsExporter,
    DailyTemperatureTendencyExporter,
    RuntimeSummary,
    XarrayExporter,
)
from pyspeedy.config import load_config
from pyspeedy.pressure_levels import DEFAULT_PRESSURE_LEVELS_HPA
from pyspeedy.speedy import Speedy, SpeedyEns

REPO_ROOT = Path(__file__).resolve().parents[2]

DAILY_TEMPERATURE_TENDENCY_EXPORT_VARS = {
    "t_tend_dyn",
    "t_tend_phy",
    "t_tend_cnv",
    "t_tend_lsc",
    "t_tend_sw",
    "t_tend_lw",
    "t_tend_pbl",
    "t_tend_hs",
}

DAILY_MOISTURE_TENDENCY_EXPORT_VARS = {
    "q_tend_dyn",
    "q_tend_phy",
    "q_tend_cnv",
    "q_tend_lsc",
    "q_tend_pbl",
}

DAILY_MOMENTUM_TENDENCY_EXPORT_VARS = {
    "u_tend_dyn",
    "v_tend_dyn",
    "u_tend_phy",
    "v_tend_phy",
    "u_tend_pbl",
    "v_tend_pbl",
    "u_tend_gwd",
    "v_tend_gwd",
    "u_tend_hs",
    "v_tend_hs",
}

DAILY_STRESS_EXPORT_VARS = {
    "u_stress",
    "v_stress",
}

DAILY_CLOUD_EXPORT_VARS = {
    "cloud_cover",
    "stratiform_cloud_cover",
    "total_cloud_top_pressure",
    "conv_cloud_top_pressure",
    "column_water_vapor",
}

DAILY_HYDROLOGY_EXPORT_VARS = {
    "precip",
    "evap",
}

DAILY_LAND_EXPORT_VARS = {
    "soil_avail_water",
}

DAILY_TOA_RADIATION_EXPORT_VARS = {
    "toa_sw_down",
    "toa_sw_up",
    "toa_sw_net",
    "olr",
}

DAILY_SURFACE_FLUX_EXPORT_VARS = {
    "surface_lh_flux",
    "surface_sh_flux",
    "surface_sw_down",
    "surface_sw_up",
    "surface_sw_net",
    "surface_lw_down",
    "surface_lw_up",
    "surface_lw_net",
}

DAILY_TENDENCY_EXPORT_VARS = (
    DAILY_TEMPERATURE_TENDENCY_EXPORT_VARS
    | DAILY_MOISTURE_TENDENCY_EXPORT_VARS
    | DAILY_MOMENTUM_TENDENCY_EXPORT_VARS
    | DAILY_STRESS_EXPORT_VARS
    | DAILY_CLOUD_EXPORT_VARS
    | DAILY_HYDROLOGY_EXPORT_VARS
    | DAILY_LAND_EXPORT_VARS
    | DAILY_TOA_RADIATION_EXPORT_VARS
    | DAILY_SURFACE_FLUX_EXPORT_VARS
)

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


def _monthly_file_name(model, year, month, start_day, end_day, suffix=".nc"):
    return (
        f"{model.output_tag}_{year:04d}-{month:02d}_d{start_day:02d}-d{end_day:02d}{suffix}"
    )


def _seasonal_file_name(
    model,
    season_year,
    season,
    start_month,
    start_day,
    end_month,
    end_day,
    suffix="_seasonal_mean.nc",
):
    return (
        f"{model.output_tag}_{season_year:04d}-{season}_m{start_month:02d}d{start_day:02d}"
        f"-m{end_month:02d}d{end_day:02d}{suffix}"
    )


@pytest.mark.parametrize("start_date, end_date", start_dates)
def test_speedy_run(start_date, end_date):
    """Run speedy and verify the exported dataset matches the configured grid."""

    config = load_config()

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        output_start_day = start_date.day + 1
        output_end_day = end_date.day
        file_name = _monthly_file_name(
            model,
            end_date.year,
            end_date.month,
            output_start_day,
            output_end_day,
        )
        pressure_file_name = _monthly_file_name(
            model,
            end_date.year,
            end_date.month,
            output_start_day,
            output_end_day,
            suffix="_p.nc",
        )
        model.set_bc()
        model.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir)])

        model_file = os.path.join(tmp_work_dir, file_name)
        pressure_file = os.path.join(tmp_work_dir, pressure_file_name)
        model_ds = xr.open_dataset(model_file)
        pressure_ds = xr.open_dataset(pressure_file)
        assert model_ds.sizes["lon"] == config.model.ix
        assert model_ds.sizes["lat"] == config.model.il
        assert model_ds.sizes["lev"] == config.model.kx
        assert model_ds.sizes["time"] == (end_date - start_date).days
        for data_var in model_ds.data_vars:
            assert model_ds[data_var].notnull().all().item()
        assert pressure_ds.sizes["lon"] == config.model.ix
        assert pressure_ds.sizes["lat"] == config.model.il
        assert pressure_ds.sizes["lev"] == len(DEFAULT_PRESSURE_LEVELS_HPA)
        assert pressure_ds.sizes["time"] == (end_date - start_date).days
        assert np.allclose(pressure_ds["lev"].values, DEFAULT_PRESSURE_LEVELS_HPA)
        assert pressure_ds["lev"].attrs["standard_name"] == "air_pressure"
        assert np.isfinite(pressure_ds["t"]).any().item()


def test_speedy_concurrent():
    """Run serially 2 speedy models for 4 days, alternating each model run every one day of simulation."""
    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 4)
    ndays = 3

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        tmp_work_dir1 = os.path.join(tmp_work_dir, "run1")
        tmp_work_dir2 = os.path.join(tmp_work_dir, "run2")

        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = _monthly_file_name(model, 1982, 1, 2, 4)
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
        assert model_ds.sizes["time"] == ndays

        model_file = os.path.join(tmp_work_dir2, file_name)
        model2_ds = xr.open_dataset(model_file)
        xr.testing.assert_allclose(model_ds, model2_ds, rtol=1e-06, atol=0)


def test_ens_speedy():
    """Test the SpeedyEns class."""
    num_of_members = 3
    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)
    model_ens = SpeedyEns(num_of_members, start_date=start_date, end_date=end_date)
    file_name = _monthly_file_name(model_ens, 1982, 1, 2, 2)
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

    with pytest.raises(ValueError):
        DailyDiagnosticsExporter(interval=1)

    with pytest.raises(ValueError):
        RuntimeSummary(interval=1)


def test_runtime_summary_output(capsys):
    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)

    model = Speedy(start_date=start_date, end_date=end_date)
    model.set_bc()
    model.run(callbacks=[RuntimeSummary(interval=load_config().model.nsteps, verbose=True)])

    captured = capsys.readouterr()
    assert "Tlow=" in captured.out
    assert "P=" in captured.out
    assert "TOA_net=" in captured.out
    assert "max|u|=" in captured.out
    assert "min(q)=" in captured.out


def test_gwd_defaults_loaded_from_config():
    """A new model instance should pick up the YAML-backed GWD defaults."""

    config = load_config()
    model = Speedy(start_date=datetime(1982, 1, 1), end_date=datetime(1982, 1, 2))

    assert bool(model["orographic_gwd_enabled"]) is config.gravity_wave_drag.enabled
    assert float(model["gwd_time_scale_days"]) == pytest.approx(config.gravity_wave_drag.time_scale_days)
    assert float(model["gwd_oro_threshold_m"]) == pytest.approx(config.gravity_wave_drag.orography_threshold_m)
    assert float(model["gwd_oro_scale_m"]) == pytest.approx(config.gravity_wave_drag.orography_scale_m)
    assert float(model["gwd_launch_sigma"]) == pytest.approx(config.gravity_wave_drag.launch_sigma)


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
        file_name = _monthly_file_name(model, 1982, 1, 2, 2)
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
        file_name = _monthly_file_name(model, 1982, 1, 2, 2)
        model.set_bc()

        exporter = XarrayExporter(output_dir=tmp_work_dir, variables=variables)
        model.run(callbacks=[exporter])

        model_file = os.path.join(tmp_work_dir, file_name)
        model_ds = xr.open_dataset(model_file)

        assert set((v.replace("_grid", "") for v in variables)) == set(model_ds.keys())


def test_pressure_export_skipped_when_no_sigma_variables():
    """Do not write a pressure-level companion file when only 2D fields are exported."""

    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = _monthly_file_name(model, 1982, 1, 2, 2)
        pressure_file_name = _monthly_file_name(model, 1982, 1, 2, 2, suffix="_p.nc")
        model.set_bc()

        exporter = XarrayExporter(output_dir=tmp_work_dir, variables=["precnv", "precls"])
        model.run(callbacks=[exporter])

        assert os.path.exists(os.path.join(tmp_work_dir, file_name))
        assert not os.path.exists(os.path.join(tmp_work_dir, pressure_file_name))


def test_daily_temperature_tendency_export():
    """Merge daily diagnostics into the main monthly output file."""

    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)
    config = load_config()

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = _monthly_file_name(model, 1982, 1, 2, 2)
        pressure_file_name = _monthly_file_name(model, 1982, 1, 2, 2, suffix="_p.nc")
        model.set_bc()

        model.run(
            callbacks=[
                XarrayExporter(output_dir=tmp_work_dir),
                DailyDiagnosticsExporter(output_dir=tmp_work_dir),
            ]
        )

        tendency_ds = xr.open_dataset(os.path.join(tmp_work_dir, file_name))
        pressure_ds = xr.open_dataset(os.path.join(tmp_work_dir, pressure_file_name))

        assert DAILY_TENDENCY_EXPORT_VARS.issubset(tendency_ds.data_vars)
        assert {"u", "v", "t", "q", "phi", "ps"}.issubset(tendency_ds.data_vars)
        assert not os.path.exists(os.path.join(tmp_work_dir, _monthly_file_name(model, 1982, 1, 2, 2, suffix="_tend.nc")))
        assert not os.path.exists(
            os.path.join(tmp_work_dir, _monthly_file_name(model, 1982, 1, 2, 2, suffix="_tend_p.nc"))
        )
        for var_name in ("u", "v", "t", "q", "phi", "ps"):
            assert tendency_ds[var_name].notnull().all().item()
        for var_name in DAILY_TEMPERATURE_TENDENCY_EXPORT_VARS:
            assert tendency_ds[var_name].attrs["units"] == "K/day"
            assert tendency_ds[var_name].notnull().all().item()
        for var_name in DAILY_MOISTURE_TENDENCY_EXPORT_VARS:
            assert tendency_ds[var_name].attrs["units"] == "g/kg/day"
            assert tendency_ds[var_name].notnull().all().item()
        for var_name in DAILY_MOMENTUM_TENDENCY_EXPORT_VARS:
            assert tendency_ds[var_name].attrs["units"] == "m/s/day"
            assert tendency_ds[var_name].notnull().all().item()
        for var_name in DAILY_STRESS_EXPORT_VARS:
            assert tendency_ds[var_name].notnull().all().item()
        assert tendency_ds["cloud_cover"].attrs["units"] == "1"
        assert tendency_ds["stratiform_cloud_cover"].attrs["units"] == "1"
        assert tendency_ds["total_cloud_top_pressure"].attrs["units"] == "hPa"
        assert tendency_ds["conv_cloud_top_pressure"].attrs["units"] == "hPa"
        assert tendency_ds["column_water_vapor"].attrs["units"] == "mm"
        assert tendency_ds["cloud_cover"].notnull().all().item()
        assert tendency_ds["stratiform_cloud_cover"].notnull().all().item()
        assert tendency_ds["column_water_vapor"].notnull().all().item()
        assert np.isfinite(tendency_ds["total_cloud_top_pressure"]).any().item()
        assert np.isfinite(tendency_ds["conv_cloud_top_pressure"]).any().item()
        for var_name in DAILY_HYDROLOGY_EXPORT_VARS:
            assert tendency_ds[var_name].attrs["units"] == "mm/day"
            assert tendency_ds[var_name].notnull().all().item()
        assert tendency_ds["soil_avail_water"].attrs["units"] == "1"
        assert tendency_ds["soil_avail_water"].notnull().all().item()
        for var_name in DAILY_TOA_RADIATION_EXPORT_VARS | DAILY_SURFACE_FLUX_EXPORT_VARS:
            assert tendency_ds[var_name].attrs["units"] == "W/m^2"
            assert tendency_ds[var_name].notnull().all().item()
        assert pressure_ds.sizes["lev"] == len(DEFAULT_PRESSURE_LEVELS_HPA)
        for var_name in DAILY_TENDENCY_EXPORT_VARS:
            assert np.isfinite(pressure_ds[var_name]).any().item()
        assert float(tendency_ds["cloud_cover"].min().item()) >= 0.0
        assert float(tendency_ds["cloud_cover"].max().item()) <= 1.0
        assert float(tendency_ds["stratiform_cloud_cover"].min().item()) >= 0.0
        assert float(tendency_ds["stratiform_cloud_cover"].max().item()) <= 1.0
        assert float(tendency_ds["column_water_vapor"].min().item()) >= 0.0
        assert float(tendency_ds["soil_avail_water"].min().item()) >= 0.0
        assert float(tendency_ds["soil_avail_water"].max().item()) <= 1.0
        assert float(np.nanmin(tendency_ds["total_cloud_top_pressure"].values)) >= 0.0
        assert float(np.nanmin(tendency_ds["conv_cloud_top_pressure"].values)) >= 0.0
        assert float(tendency_ds["precip"].min().item()) >= 0.0
        for dataset in (tendency_ds, pressure_ds):
            toa_sw_budget = dataset["toa_sw_down"] - dataset["toa_sw_up"] - dataset["toa_sw_net"]
            surface_sw_budget = dataset["surface_sw_down"] - dataset["surface_sw_up"] - dataset["surface_sw_net"]
            surface_lw_budget = dataset["surface_lw_down"] - dataset["surface_lw_up"] - dataset["surface_lw_net"]
            assert float(np.abs(toa_sw_budget).max().item()) < 1.0e-4
            assert float(np.abs(surface_sw_budget).max().item()) < 1.0e-4
            assert float(np.abs(surface_lw_budget).max().item()) < 1.0e-4
        assert np.allclose(
            tendency_ds["q_tend_cnv"] + tendency_ds["q_tend_lsc"] + tendency_ds["q_tend_pbl"],
            tendency_ds["q_tend_phy"],
            atol=1.0e-4,
        )
        if config.gravity_wave_drag.enabled:
            assert float(np.abs(tendency_ds["u_tend_gwd"]).max().item()) > 0.0
            assert float(np.abs(tendency_ds["v_tend_gwd"]).max().item()) > 0.0
        else:
            assert float(np.abs(tendency_ds["u_tend_gwd"]).max().item()) == pytest.approx(0.0)
            assert float(np.abs(tendency_ds["v_tend_gwd"]).max().item()) == pytest.approx(0.0)

        summary = subprocess.run(
            [
                sys.executable,
                str(REPO_ROOT / "scripts" / "inspect_tendency.py"),
                os.path.join(tmp_work_dir, file_name),
            ],
            check=True,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
        )
        assert "u_stress" in summary.stdout
        assert "u_tend_phy" in summary.stdout
        assert "q_tend_phy" in summary.stdout
        assert "cloud_cover" in summary.stdout
        assert "soil_avail_water" in summary.stdout
        assert "surface_sw_net" in summary.stdout


def test_daily_temperature_tendency_export_with_gwd():
    """Enable GWD and verify the exported momentum split includes it."""

    start_date = datetime(1982, 1, 1)
    end_date = datetime(1982, 1, 2)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        file_name = _monthly_file_name(model, 1982, 1, 2, 2)
        model["orographic_gwd_enabled"] = True
        model["gwd_oro_threshold_m"] = 0.0
        model["gwd_time_scale_days"] = 2.0
        model["gwd_launch_sigma"] = 0.9
        model.set_bc()

        exporter = DailyDiagnosticsExporter(output_dir=tmp_work_dir)
        model.run(callbacks=[exporter])

        tendency_ds = xr.open_dataset(os.path.join(tmp_work_dir, file_name))

        assert np.isfinite(tendency_ds["u_tend_gwd"]).any().item()
        assert np.isfinite(tendency_ds["v_tend_gwd"]).any().item()
        assert float(np.abs(tendency_ds["u_tend_gwd"]).max().item()) > 0.0
        assert float(np.abs(tendency_ds["v_tend_gwd"]).max().item()) > 0.0
        xr.testing.assert_allclose(
            tendency_ds["u_tend_phy"],
            tendency_ds["u_tend_pbl"] + tendency_ds["u_tend_gwd"],
            rtol=1e-5,
            atol=1e-5,
        )
        xr.testing.assert_allclose(
            tendency_ds["v_tend_phy"],
            tendency_ds["v_tend_pbl"] + tendency_ds["v_tend_gwd"],
            rtol=1e-5,
            atol=1e-5,
        )


def test_daily_temperature_tendency_exporter_alias():
    assert DailyTemperatureTendencyExporter is DailyDiagnosticsExporter


def test_daily_diagnostics_write_monthly_and_seasonal_means():
    start_date = datetime(1982, 1, 30)
    end_date = datetime(1982, 2, 3)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        january_mean_file = _monthly_file_name(model, 1982, 1, 31, 31, suffix="_monthly_mean.nc")
        february_mean_file = _monthly_file_name(model, 1982, 2, 1, 3, suffix="_monthly_mean.nc")
        january_pressure_mean_file = _monthly_file_name(model, 1982, 1, 31, 31, suffix="_monthly_mean_p.nc")
        february_pressure_mean_file = _monthly_file_name(model, 1982, 2, 1, 3, suffix="_monthly_mean_p.nc")
        seasonal_mean_file = _seasonal_file_name(model, 1982, "DJF", 1, 31, 2, 3)
        seasonal_pressure_mean_file = _seasonal_file_name(
            model,
            1982,
            "DJF",
            1,
            31,
            2,
            3,
            suffix="_seasonal_mean_p.nc",
        )
        model.set_bc()

        model.run(
            callbacks=[
                XarrayExporter(output_dir=tmp_work_dir),
                DailyDiagnosticsExporter(output_dir=tmp_work_dir),
            ]
        )

        january_mean_ds = xr.open_dataset(os.path.join(tmp_work_dir, january_mean_file))
        february_mean_ds = xr.open_dataset(os.path.join(tmp_work_dir, february_mean_file))
        seasonal_mean_ds = xr.open_dataset(os.path.join(tmp_work_dir, seasonal_mean_file))

        assert os.path.exists(os.path.join(tmp_work_dir, january_pressure_mean_file))
        assert os.path.exists(os.path.join(tmp_work_dir, february_pressure_mean_file))
        assert os.path.exists(os.path.join(tmp_work_dir, seasonal_pressure_mean_file))

        assert january_mean_ds.attrs["aggregation"] == "monthly_mean"
        assert february_mean_ds.attrs["aggregation"] == "monthly_mean"
        assert seasonal_mean_ds.attrs["aggregation"] == "seasonal_mean"
        assert int(january_mean_ds.attrs["sample_count"]) == 1
        assert int(february_mean_ds.attrs["sample_count"]) == 3
        assert int(seasonal_mean_ds.attrs["sample_count"]) == 4
        assert january_mean_ds.sizes["time"] == 1
        assert seasonal_mean_ds.sizes["time"] == 1
        assert "u" in seasonal_mean_ds.data_vars
        assert "t_tend_dyn" in seasonal_mean_ds.data_vars


def test_monthly_export_rolls_over_at_calendar_month():
    start_date = datetime(1982, 1, 30)
    end_date = datetime(1982, 2, 3)

    with tempfile.TemporaryDirectory() as tmp_work_dir:
        model = Speedy(start_date=start_date, end_date=end_date)
        january_file = _monthly_file_name(model, 1982, 1, 31, 31)
        february_file = _monthly_file_name(model, 1982, 2, 1, 3)
        january_pressure_file = _monthly_file_name(model, 1982, 1, 31, 31, suffix="_p.nc")
        february_pressure_file = _monthly_file_name(model, 1982, 2, 1, 3, suffix="_p.nc")
        model.set_bc()

        model.run(callbacks=[XarrayExporter(output_dir=tmp_work_dir)])

        january_ds = xr.open_dataset(os.path.join(tmp_work_dir, january_file))
        february_ds = xr.open_dataset(os.path.join(tmp_work_dir, february_file))

        assert january_ds.sizes["time"] == 1
        assert february_ds.sizes["time"] == 3
        assert os.path.exists(os.path.join(tmp_work_dir, january_pressure_file))
        assert os.path.exists(os.path.join(tmp_work_dir, february_pressure_file))
