"""
pySpeedy main module
====================

.. autosummary::
    :toctree: ../generated/

    example_bc_file
    example_sst_anomaly_file
    load_config
    Speedy
    SpeedyEns
"""
from pathlib import Path

PACKAGE_DATA_DIR = Path(__file__).parent / "data"
_IMPORT_ERROR = None

from .config import load_config  # noqa

DEFAULT_OUTPUT_VARS = (
    "u_grid",
    "v_grid",
    "t_grid",
    "q_grid",
    "phi_grid",
    "ps_grid",
)


def example_bc_file():
    """Returns the Path to the example bc file."""
    return str(PACKAGE_DATA_DIR / "example_bc.nc")


def example_sst_anomaly_file():
    """Returns the Path to the example SST anomaly file."""
    return str(PACKAGE_DATA_DIR / "sst_anomaly.nc")


try:
    from .speedy_driver import speedy_driver as _speedy  # noqa
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only before build/install
    _speedy = None
    _IMPORT_ERROR = exc
else:
    from .speedy import Speedy, SpeedyEns, MODEL_STATE_DEF  # noqa


def __getattr__(name):
    if name in {"Speedy", "SpeedyEns", "MODEL_STATE_DEF"} and _IMPORT_ERROR is not None:
        raise ModuleNotFoundError(
            "The compiled pySPEEDY extension is not available. Run `python -m pip install -e .` first."
        ) from _IMPORT_ERROR
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
