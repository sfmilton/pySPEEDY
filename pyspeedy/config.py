"""Shared model and run configuration helpers."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime
from functools import lru_cache
import json
from pathlib import Path

import numpy as np

from .generated_model_config import MODEL_CONFIG as GENERATED_MODEL_CONFIG


PACKAGE_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG_PATH = PACKAGE_ROOT / "data" / "model_config.yml"


@dataclass(frozen=True)
class ModelConfig:
    trunc: int
    ix: int
    iy: int
    kx: int
    ntr: int
    t_levs: int
    aux_dim: int
    nsteps: int
    iseasc: int
    nstrad: int
    sppt_on: bool
    rob: float
    wil: float
    alph: float

    @property
    def il(self) -> int:
        return 2 * self.iy

    @property
    def mx(self) -> int:
        return self.trunc + 1

    @property
    def nx(self) -> int:
        return self.trunc + 2

    @property
    def lon_step(self) -> float:
        return 360.0 / self.ix

    @property
    def timestep_seconds(self) -> float:
        return 86400.0 / self.nsteps

    @property
    def spectral_tag(self) -> str:
        return f"T{self.trunc}"


@dataclass(frozen=True)
class RunConfig:
    start_date: datetime
    end_date: datetime
    output_dir: str
    history_interval: int
    diag_interval: int
    verbose_output: bool
    output_vars: tuple[str, ...] | None


@dataclass(frozen=True)
class HeldSuarezConfig:
    equilibrium_surface_temperature: float
    equator_pole_temperature_contrast: float
    vertical_temperature_contrast: float
    minimum_equilibrium_temperature: float
    boundary_layer_sigma: float
    upper_air_relaxation_days: float
    boundary_layer_relaxation_days: float
    rayleigh_drag_days: float
    minimum_pressure_ratio: float


@dataclass(frozen=True)
class PySpeedyConfig:
    model: ModelConfig
    run: RunConfig
    held_suarez: HeldSuarezConfig


def _coerce_datetime(value: date | datetime | str) -> datetime:
    if isinstance(value, datetime):
        return value
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day)
    if isinstance(value, str):
        return datetime.fromisoformat(value)
    raise TypeError(f"Unsupported datetime value: {value!r}")


def _load_yaml(path: Path) -> dict:
    with path.open() as stream:
        filtered = "\n".join(
            line for line in stream.read().splitlines() if not line.lstrip().startswith("#")
        )
        data = json.loads(filtered)
    if not isinstance(data, dict):
        raise ValueError(f"Invalid config file: {path}")
    return data


@lru_cache(maxsize=None)
def load_config(path: str | Path | None = None) -> PySpeedyConfig:
    config_path = Path(path) if path is not None else DEFAULT_CONFIG_PATH
    raw = _load_yaml(config_path)
    model = ModelConfig(**GENERATED_MODEL_CONFIG)
    run = RunConfig(
        start_date=_coerce_datetime(raw["run"]["start_date"]),
        end_date=_coerce_datetime(raw["run"]["end_date"]),
        output_dir=raw["run"]["output_dir"],
        history_interval=raw["run"]["history_interval"],
        diag_interval=raw["run"]["diag_interval"],
        verbose_output=raw["run"].get("verbose_output", False),
        output_vars=tuple(raw["run"]["output_vars"]) if raw["run"].get("output_vars") else None,
    )
    held_suarez_raw = raw.get("held_suarez", {})
    held_suarez = HeldSuarezConfig(
        equilibrium_surface_temperature=held_suarez_raw.get("equilibrium_surface_temperature", 315.0),
        equator_pole_temperature_contrast=held_suarez_raw.get("equator_pole_temperature_contrast", 60.0),
        vertical_temperature_contrast=held_suarez_raw.get("vertical_temperature_contrast", 10.0),
        minimum_equilibrium_temperature=held_suarez_raw.get("minimum_equilibrium_temperature", 200.0),
        boundary_layer_sigma=held_suarez_raw.get("boundary_layer_sigma", 0.7),
        upper_air_relaxation_days=held_suarez_raw.get("upper_air_relaxation_days", 40.0),
        boundary_layer_relaxation_days=held_suarez_raw.get("boundary_layer_relaxation_days", 4.0),
        rayleigh_drag_days=held_suarez_raw.get("rayleigh_drag_days", 1.0),
        minimum_pressure_ratio=held_suarez_raw.get("minimum_pressure_ratio", 1.0e-4),
    )
    return PySpeedyConfig(model=model, run=run, held_suarez=held_suarez)


def gaussian_latitudes(iy: int) -> np.ndarray:
    """Return the model latitudes in degrees, matching the Fortran geometry module."""
    il = 2 * iy
    latitudes = np.empty(il, dtype=np.float64)
    for j in range(1, iy + 1):
        jj = il + 1 - j
        sia_half = np.cos(np.pi * (j - 0.25) / (il + 0.5))
        angle = np.degrees(np.arcsin(sia_half))
        latitudes[j - 1] = -angle
        latitudes[jj - 1] = angle
    return latitudes


def longitudes(ix: int) -> np.ndarray:
    """Return the model longitudes in degrees."""
    return np.arange(ix, dtype=np.float64) * (360.0 / ix)


def grid_coordinates(model_config: ModelConfig | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return target longitudes and latitudes for the configured model grid."""
    cfg = model_config or load_config().model
    return longitudes(cfg.ix), gaussian_latitudes(cfg.iy)
