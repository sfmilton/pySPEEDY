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
import importlib
import importlib.util
import json
from pathlib import Path
import sys

PACKAGE_DATA_DIR = Path(__file__).parent / "data"
PACKAGE_ROOT = PACKAGE_DATA_DIR.parent
REPO_ROOT = PACKAGE_ROOT.parent
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


def _expected_speedy_driver_symbols():
    state_def = json.loads((PACKAGE_DATA_DIR / "model_state.json").read_text())
    symbols = set()
    for var_name, var_meta in state_def.items():
        if "class(" in (var_meta["dtype"] or "").lower():
            continue
        symbols.add(f"get_{var_name}")
        symbols.add(f"set_{var_name}")
        symbols.add(f"is_array_{var_name}")
        if var_meta["dims"] is not None:
            symbols.add(f"get_{var_name}_shape")
    return symbols


def _speedy_driver_missing_symbols(module):
    driver = module.speedy_driver
    return sorted(
        symbol for symbol in _expected_speedy_driver_symbols() if not hasattr(driver, symbol)
    )


def _load_speedy_driver_module(module_path=None):
    module_name = f"{__name__}.speedy_driver"
    if module_path is None:
        return importlib.import_module(".speedy_driver", __name__)

    sys.modules.pop(module_name, None)
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise ModuleNotFoundError(f"Could not load extension module from {module_path}.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _candidate_speedy_driver_paths():
    build_dir = REPO_ROOT / "build"
    if not build_dir.exists():
        return []
    return sorted(
        build_dir.glob("*/speedy_driver*.so"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )


def _import_speedy_driver():
    module = _load_speedy_driver_module()
    missing = _speedy_driver_missing_symbols(module)
    if not missing:
        return module.speedy_driver

    local_path = getattr(module, "__file__", "<unknown>")
    for candidate in _candidate_speedy_driver_paths():
        try:
            candidate_module = _load_speedy_driver_module(candidate)
        except Exception:
            continue
        candidate_missing = _speedy_driver_missing_symbols(candidate_module)
        if not candidate_missing:
            return candidate_module.speedy_driver

    preview = ", ".join(missing[:5])
    raise ModuleNotFoundError(
        "The compiled pySPEEDY extension is out of sync with the Python sources. "
        f"Loaded stale extension: {local_path}. Missing wrapper symbols such as: {preview}. "
        "Rebuild the extension and ensure the rebuilt `speedy_driver*.so` is the one being imported."
    )


try:
    _speedy = _import_speedy_driver()  # noqa
except ModuleNotFoundError as exc:  # pragma: no cover - exercised only before build/install
    _speedy = None
    _IMPORT_ERROR = exc
else:
    from .speedy import Speedy, SpeedyEns, MODEL_STATE_DEF  # noqa


def __getattr__(name):
    if name in {"Speedy", "SpeedyEns", "MODEL_STATE_DEF"} and _IMPORT_ERROR is not None:
        raise ModuleNotFoundError(str(_IMPORT_ERROR)) from _IMPORT_ERROR
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
