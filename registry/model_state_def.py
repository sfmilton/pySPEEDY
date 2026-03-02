from collections import defaultdict

from jinja2 import FileSystemLoader, Environment
from pathlib import Path
import json
import re

try:
    # Used to export the model state vars into an excel file
    import xlsxwriter

    XLSXWRITER_IMPORTED = True
except ImportError:
    XLSXWRITER_IMPORTED = False

try:
    # Used to export the model state vars into an excel file
    import pandas as pd

    PANDAS_IMPORTED = True
except ImportError:
    PANDAS_IMPORTED = False

THIS_FILE_DIR = Path(__file__).parent
SOURCES_DIR = (THIS_FILE_DIR / "../speedy.f90").resolve()
DOCS_DIR = (THIS_FILE_DIR / "../docs").resolve()
PYSPEEDY_DATA_DIR = (THIS_FILE_DIR / "../pyspeedy/data").resolve()

NC_DIMS_LUT = {"ix": "lon", "il": "lat", "kx": "lev"}


class VarDef:
    def __init__(
        self,
        name,
        dtype,
        dims,
        desc,
        std_name=None,
        units=None,
        time_dim=None,
        alt_name=None,
        value=None,
    ):
        """
        If time_dim is not None, the variable is not allocated during the
        initialization since it depends on the duration of the simulation.

        Parameters
        ----------
        name: str
            SPEEDY.f90 variable name.
        dtype: str
            Fortran data type used in the SPEEDY.f90 model.
            E.g.: "complex(8)"
        dims: str or None
            Variable dimensions in the SPEEDY.f90 model. E.g.: "(mx, nx, kx, t_levs)"
        desc: str
            Variable description.
        std_name: str
            CF standard name associated with the variable. If std_name is not provided, the `name` is used.
        units: str
            Variable units
        time_dim: str
            Time dimension. Only required if the variable has a time dimension.
        alt_name: str
            Alternative name used for the variable. For example, name used in the NETCDF files.
        value: str
            Default value used for the variable initialization.
        """
        self.name = name
        self.dtype = dtype
        self.dims = dims
        self.nc_dims = None

        if dims is not None:
            for dim, nc_dim in NC_DIMS_LUT.items():
                dims = re.sub(r"\b%s\b" % dim, nc_dim, dims)
            self.nc_dims = dims.replace(" ", "").replace("(", "").replace(")", "")
            self.nc_dims = [s for s in self.nc_dims.split(",")]

        self.desc = desc
        self.time_dim = time_dim
        self.units = units
        self.alt_name = name
        self.std_name = name

        if alt_name is not None:
            self.alt_name = alt_name
        if std_name is None:
            self.std_name = name

        self.is_module_instance = "class" in dtype.lower()
        self.value = value

    @property
    def dimension(self):
        if self.dims:
            dimension = ", ".join(":" * self.ndim)
            return f"dimension({dimension})"

    @property
    def dimension_args(self):
        if self.dims:
            return self.dims.replace("(", "").replace(")", "")

    @property
    def dimension_args_declaration(self):
        if self.dims:
            return f"integer, intent(in) :: {self.dimension_args}"

    @property
    def ndim(self):
        if self.dims:
            return len(self.dims.split(","))

    def __repr__(self):
        return str({k: v for k, v in self.__dict__.items() if not k.startswith("_")})


model_state = [
    ########################################
    # Model integration control variables
    ########################################
    VarDef("current_step", "integer", None, "Current model step."),
    VarDef(
        "ttend_daily_count",
        "integer",
        None,
        "Number of time steps accumulated in the current daily tendency mean.",
    ),
    VarDef(
        "cloud_daily_count",
        "integer",
        None,
        "Number of radiation-step samples accumulated in the current daily cloud diagnostics mean.",
    ),
    ########################################
    # Prognostic variables (spectral domain)
    ########################################
    VarDef("vor", "complex(8)", "(mx, nx, kx, t_levs)", "Vorticity"),
    VarDef("div", "complex(8)", "(mx, nx, kx, t_levs)", "Divergence"),
    VarDef("t", "complex(8)", "(mx, nx, kx, t_levs)", "Temperature", "K"),
    VarDef(
        "ps",
        "complex(8)",
        "(mx, nx, t_levs)",
        "Log of (normalised) surface pressure",
        "(p_s/p0)",
    ),
    VarDef(
        "tr",
        "complex(8)",
        "(mx, nx, kx, t_levs,ntr)",
        "Tracers (tr(1): specific humidity in g/kg)",
    ),
    VarDef(
        "phi",
        "complex(8)",
        "(mx, nx, kx)",
        "Atmospheric geopotential",
        std_name="geopotential_height",
        units="m",
    ),
    VarDef("phis", "complex(8)", "(mx, nx)", "Surface geopotential"),
    ####################################
    # Prognostic variables (grid domain)
    ####################################
    VarDef(
        "u_grid",
        "real(8)",
        "(ix, il, kx)",
        "eastward_wind",
        units="m/s",
        alt_name="u",
        std_name="eastward_wind",
    ),
    VarDef(
        "v_grid",
        "real(8)",
        "(ix, il, kx)",
        "northward_wind",
        units="m/s",
        alt_name="v",
        std_name="northward_wind",
    ),
    VarDef(
        "t_grid",
        "real(8)",
        "(ix, il, kx)",
        "air_temperature",
        units="K",
        alt_name="t",
        std_name="air_temperature",
    ),
    VarDef(
        "q_grid", "real(8)", "(ix, il, kx)", "specific_humidity", "Kg/Kg", alt_name="q"
    ),
    VarDef(
        "phi_grid",
        "real(8)",
        "(ix, il, kx)",
        "geopotential_height",
        "m",
        alt_name="phi",
    ),
    VarDef(
        "ps_grid", "real(8)", "(ix, il)", "surface_air_pressure", "Pa", alt_name="ps"
    ),
    VarDef(
        "ttend_dyn_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily dynamical temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_phy_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily physical temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_cnv_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily convective temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_lsc_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily large-scale condensation temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_sw_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily shortwave-radiative temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_lw_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily longwave-radiative temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_pbl_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily boundary-layer and surface-flux temperature tendency",
        units="K/s",
    ),
    VarDef(
        "ttend_hs_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily Held-Suarez temperature tendency",
        units="K/s",
    ),
    VarDef(
        "qtend_dyn_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily dynamical specific humidity tendency",
        units="g/kg/s",
    ),
    VarDef(
        "qtend_phy_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily physical specific humidity tendency",
        units="g/kg/s",
    ),
    VarDef(
        "qtend_cnv_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily convective specific humidity tendency",
        units="g/kg/s",
    ),
    VarDef(
        "qtend_lsc_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily large-scale-condensation specific humidity tendency",
        units="g/kg/s",
    ),
    VarDef(
        "qtend_pbl_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily boundary-layer and surface-flux specific humidity tendency",
        units="g/kg/s",
    ),
    VarDef(
        "utend_dyn_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily dynamical eastward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "vtend_dyn_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily dynamical northward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "utend_phy_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily physical eastward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "vtend_phy_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily physical northward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "utend_pbl_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily boundary-layer and surface-flux eastward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "vtend_pbl_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily boundary-layer and surface-flux northward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "utend_gwd_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily orographic gravity-wave-drag eastward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "vtend_gwd_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily orographic gravity-wave-drag northward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "utend_hs_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily Held-Suarez eastward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "vtend_hs_accum",
        "real(8)",
        "(ix, il, kx)",
        "Accumulated daily Held-Suarez northward-wind tendency",
        units="m/s^2",
    ),
    VarDef(
        "ustr_sfc_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily area-weighted surface eastward stress diagnostic",
    ),
    VarDef(
        "vstr_sfc_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily area-weighted surface northward stress diagnostic",
    ),
    VarDef(
        "cloud_cover_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily total cloud-cover diagnostic",
    ),
    VarDef(
        "stratiform_cloud_cover_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily stratiform cloud-cover diagnostic",
    ),
    VarDef(
        "total_cloud_top_pressure_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily total cloud-top pressure diagnostic",
        units="hPa",
    ),
    VarDef(
        "total_cloud_top_count",
        "real(8)",
        "(ix, il)",
        "Number of valid samples accumulated in the daily total cloud-top pressure diagnostic",
    ),
    VarDef(
        "conv_cloud_top_pressure_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily convective cloud-top pressure diagnostic",
        units="hPa",
    ),
    VarDef(
        "conv_cloud_top_count",
        "real(8)",
        "(ix, il)",
        "Number of valid samples accumulated in the daily convective cloud-top pressure diagnostic",
    ),
    VarDef(
        "column_water_vapor_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily column water vapor diagnostic",
        units="mm",
    ),
    VarDef(
        "precip_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily total precipitation diagnostic",
        units="g/(m^2 s)",
    ),
    VarDef(
        "evap_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily surface evaporation diagnostic",
        units="g/(m^2 s)",
    ),
    VarDef(
        "toa_sw_down_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily incoming top-of-atmosphere shortwave radiation",
        units="W/m^2",
    ),
    VarDef(
        "toa_sw_up_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily outgoing top-of-atmosphere shortwave radiation",
        units="W/m^2",
    ),
    VarDef(
        "toa_sw_net_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily net top-of-atmosphere shortwave radiation (downward positive)",
        units="W/m^2",
    ),
    VarDef(
        "olr_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily outgoing longwave radiation at the top of atmosphere",
        units="W/m^2",
    ),
    VarDef(
        "surface_lh_flux_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily upward latent heat flux at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_sh_flux_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily upward sensible heat flux at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_sw_down_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily downward shortwave radiation at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_sw_up_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily upward shortwave radiation at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_sw_net_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily net shortwave radiation at the surface (downward positive)",
        units="W/m^2",
    ),
    VarDef(
        "surface_lw_down_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily downward longwave radiation at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_lw_up_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily upward longwave radiation at the surface",
        units="W/m^2",
    ),
    VarDef(
        "surface_lw_net_accum",
        "real(8)",
        "(ix, il)",
        "Accumulated daily net longwave radiation at the surface (downward positive)",
        units="W/m^2",
    ),
    VarDef(
        "ttend_dyn_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean dynamical temperature tendency",
        units="K/day",
        alt_name="t_tend_dyn",
    ),
    VarDef(
        "ttend_phy_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean physical temperature tendency",
        units="K/day",
        alt_name="t_tend_phy",
    ),
    VarDef(
        "ttend_cnv_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean convective temperature tendency",
        units="K/day",
        alt_name="t_tend_cnv",
    ),
    VarDef(
        "ttend_lsc_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean large-scale condensation temperature tendency",
        units="K/day",
        alt_name="t_tend_lsc",
    ),
    VarDef(
        "ttend_sw_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean shortwave-radiative temperature tendency",
        units="K/day",
        alt_name="t_tend_sw",
    ),
    VarDef(
        "ttend_lw_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean longwave-radiative temperature tendency",
        units="K/day",
        alt_name="t_tend_lw",
    ),
    VarDef(
        "ttend_pbl_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean boundary-layer and surface-flux temperature tendency",
        units="K/day",
        alt_name="t_tend_pbl",
    ),
    VarDef(
        "ttend_hs_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean Held-Suarez temperature tendency",
        units="K/day",
        alt_name="t_tend_hs",
    ),
    VarDef(
        "qtend_dyn_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean dynamical specific humidity tendency",
        units="g/kg/day",
        alt_name="q_tend_dyn",
    ),
    VarDef(
        "qtend_phy_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean physical specific humidity tendency",
        units="g/kg/day",
        alt_name="q_tend_phy",
    ),
    VarDef(
        "qtend_cnv_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean convective specific humidity tendency",
        units="g/kg/day",
        alt_name="q_tend_cnv",
    ),
    VarDef(
        "qtend_lsc_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean large-scale-condensation specific humidity tendency",
        units="g/kg/day",
        alt_name="q_tend_lsc",
    ),
    VarDef(
        "qtend_pbl_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean boundary-layer and surface-flux specific humidity tendency",
        units="g/kg/day",
        alt_name="q_tend_pbl",
    ),
    VarDef(
        "utend_dyn_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean dynamical eastward-wind tendency",
        units="m/s/day",
        alt_name="u_tend_dyn",
    ),
    VarDef(
        "vtend_dyn_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean dynamical northward-wind tendency",
        units="m/s/day",
        alt_name="v_tend_dyn",
    ),
    VarDef(
        "utend_phy_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean physical eastward-wind tendency",
        units="m/s/day",
        alt_name="u_tend_phy",
    ),
    VarDef(
        "vtend_phy_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean physical northward-wind tendency",
        units="m/s/day",
        alt_name="v_tend_phy",
    ),
    VarDef(
        "utend_pbl_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean boundary-layer and surface-flux eastward-wind tendency",
        units="m/s/day",
        alt_name="u_tend_pbl",
    ),
    VarDef(
        "vtend_pbl_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean boundary-layer and surface-flux northward-wind tendency",
        units="m/s/day",
        alt_name="v_tend_pbl",
    ),
    VarDef(
        "utend_gwd_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean orographic gravity-wave-drag eastward-wind tendency",
        units="m/s/day",
        alt_name="u_tend_gwd",
    ),
    VarDef(
        "vtend_gwd_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean orographic gravity-wave-drag northward-wind tendency",
        units="m/s/day",
        alt_name="v_tend_gwd",
    ),
    VarDef(
        "utend_hs_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean Held-Suarez eastward-wind tendency",
        units="m/s/day",
        alt_name="u_tend_hs",
    ),
    VarDef(
        "vtend_hs_mean",
        "real(8)",
        "(ix, il, kx)",
        "Daily-mean Held-Suarez northward-wind tendency",
        units="m/s/day",
        alt_name="v_tend_hs",
    ),
    VarDef(
        "ustr_sfc_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean area-weighted surface eastward stress diagnostic",
        alt_name="u_stress",
    ),
    VarDef(
        "vstr_sfc_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean area-weighted surface northward stress diagnostic",
        alt_name="v_stress",
    ),
    VarDef(
        "cloud_cover_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean total cloud cover",
        units="1",
        alt_name="cloud_cover",
    ),
    VarDef(
        "stratiform_cloud_cover_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean stratiform cloud cover",
        units="1",
        alt_name="stratiform_cloud_cover",
    ),
    VarDef(
        "total_cloud_top_pressure_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean total cloud-top pressure",
        units="hPa",
        alt_name="total_cloud_top_pressure",
    ),
    VarDef(
        "conv_cloud_top_pressure_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean convective cloud-top pressure",
        units="hPa",
        alt_name="conv_cloud_top_pressure",
    ),
    VarDef(
        "column_water_vapor_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean column water vapor",
        units="mm",
        alt_name="column_water_vapor",
    ),
    VarDef(
        "precip_mean",
        "real(8)",
        "(ix, il)",
        "Daily total precipitation",
        units="mm/day",
        alt_name="precip",
    ),
    VarDef(
        "evap_mean",
        "real(8)",
        "(ix, il)",
        "Daily total evaporation",
        units="mm/day",
        alt_name="evap",
    ),
    VarDef(
        "toa_sw_down_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean incoming top-of-atmosphere shortwave radiation",
        units="W/m^2",
        alt_name="toa_sw_down",
    ),
    VarDef(
        "toa_sw_up_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean outgoing top-of-atmosphere shortwave radiation",
        units="W/m^2",
        alt_name="toa_sw_up",
    ),
    VarDef(
        "toa_sw_net_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean net top-of-atmosphere shortwave radiation (downward positive)",
        units="W/m^2",
        alt_name="toa_sw_net",
    ),
    VarDef(
        "olr_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean outgoing longwave radiation at the top of atmosphere",
        units="W/m^2",
        alt_name="olr",
    ),
    VarDef(
        "surface_lh_flux_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean upward latent heat flux at the surface",
        units="W/m^2",
        alt_name="surface_lh_flux",
    ),
    VarDef(
        "surface_sh_flux_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean upward sensible heat flux at the surface",
        units="W/m^2",
        alt_name="surface_sh_flux",
    ),
    VarDef(
        "surface_sw_down_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean downward shortwave radiation at the surface",
        units="W/m^2",
        alt_name="surface_sw_down",
    ),
    VarDef(
        "surface_sw_up_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean upward shortwave radiation at the surface",
        units="W/m^2",
        alt_name="surface_sw_up",
    ),
    VarDef(
        "surface_sw_net_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean net shortwave radiation at the surface (downward positive)",
        units="W/m^2",
        alt_name="surface_sw_net",
    ),
    VarDef(
        "surface_lw_down_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean downward longwave radiation at the surface",
        units="W/m^2",
        alt_name="surface_lw_down",
    ),
    VarDef(
        "surface_lw_up_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean upward longwave radiation at the surface",
        units="W/m^2",
        alt_name="surface_lw_up",
    ),
    VarDef(
        "surface_lw_net_mean",
        "real(8)",
        "(ix, il)",
        "Daily-mean net longwave radiation at the surface (downward positive)",
        units="W/m^2",
        alt_name="surface_lw_net",
    ),
    ################################################
    # Auxiliary variables used by the physic schemes
    ################################################
    VarDef(
        "precnv", "real(8)", "(ix, il)", "Convective precipitation, total", "g/(m^2 s)"
    ),
    VarDef(
        "precls", "real(8)", "(ix, il)", "Large-scale precipitation, total", "g/(m^2 s)"
    ),
    VarDef(
        "snowcv",
        "real(8)",
        "(ix, il)",
        "Convective precipitation, snow only",
        "g/(m^2 s)",
    ),
    VarDef(
        "snowls",
        "real(8)",
        "(ix, il)",
        "Large-scale precipitation, snow only",
        "g/(m^2 s)",
    ),
    VarDef("cbmf", "real(8)", "(ix, il)", "Cloud-base mass flux"),
    VarDef(
        "tsr", "real(8)", "(ix, il)", "Top-of-atmosphere shortwave radiation (downward)"
    ),
    VarDef(
        "ssrd", "real(8)", "(ix, il)", "Surface shortwave radiation (downward-only)"
    ),
    VarDef("ssr", "real(8)", "(ix, il)", "Surface shortwave radiation (net downward)"),
    VarDef("slrd", "real(8)", "(ix, il)", "Surface longwave radiation (downward-only)"),
    VarDef("slr", "real(8)", "(ix, il)", "Surface longwave radiation (net upward)"),
    VarDef("olr", "real(8)", "(ix, il)", "Outgoing longwave radiation (upward)"),
    # Third dimension -> 1:land, 2:sea, 3: weighted average
    VarDef("slru", "real(8)", "(ix, il,aux_dim)", "Surface longwave emission (upward)"),
    VarDef("ustr", "real(8)", "(ix, il,aux_dim)", "U-stress"),
    VarDef("vstr", "real(8)", "(ix, il,aux_dim)", "Vstress"),
    VarDef("shf", "real(8)", "(ix, il,aux_dim)", "Sensible heat flux"),
    VarDef("evap", "real(8)", "(ix, il,aux_dim)", "Evaporation", "g/(m^2 s)"),
    VarDef("hfluxn", "real(8)", "(ix, il,aux_dim)", "Net heat flux into surface"),
    #
    # Saved computations
    VarDef(
        "tt_rsw",
        "real(8)",
        "(ix, il,kx)",
        "Flux of short-wave radiation absorbed in each atmospheric layer",
    ),
    ###########################
    # Boundary module variables
    ###########################
    VarDef("phi0", "real(8)", "(ix, il)", "Unfiltered surface geopotential"),
    VarDef("orog", "real(8)", "(ix, il)", "Orography", "m"),
    VarDef("phis0", "real(8)", "(ix, il)", "Spectrally-filtered surface geopotential"),
    VarDef("alb0", "real(8)", "(ix, il)", "Bare-land annual-mean albedo"),
    VarDef("forog", "real(8)", "(ix, il)", "Orographic factor for land surface drag"),
    #################################
    # Surface fluxes module variables
    #################################
    VarDef("fmask_orig", "real(8)", "(ix, il)", "Original (fractional) land-sea mask"),
    ###############################
    # Geopotential module variables
    ###############################
    VarDef("xgeop1", "real(8)", "(kx)", "Constant 1 for hydrostatic equation"),
    VarDef("xgeop2", "real(8)", "(kx)", "Constant 2 for hydrostatic equation"),
    #############################
    # Land model module variables
    #############################
    VarDef(
        "stl12",
        "real(8)",
        "(ix, il, 12)",
        "Land surface temperature monthly-mean climatology",
    ),
    VarDef(
        "snowd12",
        "real(8)",
        "(ix, il, 12)",
        "Snow depth (water equivalent) monthly-mean climatology",
    ),
    VarDef(
        "soilw12",
        "real(8)",
        "(ix, il, 12)",
        "Soil water availability monthly-mean climatology",
    ),
    VarDef("veg_low", "real(8)", "(ix, il)", "Low vegetation fraction"),
    VarDef("veg_high", "real(8)", "(ix, il)", "High vegetation fraction"),
    VarDef("soil_wc_l1", "real(8)", "(ix, il, 12)", "Soil water content: Layer 1"),
    VarDef("soil_wc_l2", "real(8)", "(ix, il, 12)", "Soil water content: Layer 2"),
    VarDef("soil_wc_l3", "real(8)", "(ix, il, 12)", "Soil water content: Layer 3"),
    ############################
    # Sea model module variables
    ############################
    VarDef("sst12", "real(8)", "(ix, il, 12)", "Sea/ice surface temperature", "K"),
    VarDef("sea_ice_frac12", "real(8)", "(ix, il, 12)", "Sea ice fraction"),
    VarDef(
        "sst_anom",
        "real(8)",
        "(ix, il, 0:n_months+1)",
        "Observed SST anomaly (input).",
        time_dim="n_months",
    ),
    ############################
    # Shortware radiation module
    ############################
    VarDef(
        "increase_co2",
        "logical",
        None,
        " Flag for CO2 optical thickness increase",
        value=".false.",
    ),
    VarDef(
        "compute_shortwave",
        "logical",
        None,
        "Flag for shortwave radiation routine (turned on and off in main loop depending on the value of nstrad)",
        value=".true.",
    ),
    VarDef(
        "held_suarez_mode",
        "logical",
        None,
        "Flag for Held-Suarez dry-physics forcing",
        value=".false.",
    ),
    VarDef(
        "hs_trefc",
        "real(8)",
        None,
        "Held-Suarez equilibrium surface temperature",
        value="315.0",
    ),
    VarDef(
        "hs_delta_ty",
        "real(8)",
        None,
        "Held-Suarez equator-to-pole temperature contrast",
        value="60.0",
    ),
    VarDef(
        "hs_delta_theta_z",
        "real(8)",
        None,
        "Held-Suarez vertical temperature contrast",
        value="10.0",
    ),
    VarDef(
        "hs_tmin",
        "real(8)",
        None,
        "Held-Suarez minimum equilibrium temperature",
        value="200.0",
    ),
    VarDef(
        "hs_sigma_b",
        "real(8)",
        None,
        "Held-Suarez boundary-layer sigma threshold",
        value="0.7",
    ),
    VarDef(
        "hs_tau_a_days",
        "real(8)",
        None,
        "Held-Suarez upper-air relaxation time in days",
        value="40.0",
    ),
    VarDef(
        "hs_tau_s_days",
        "real(8)",
        None,
        "Held-Suarez boundary-layer relaxation time in days",
        value="4.0",
    ),
    VarDef(
        "hs_tau_f_days",
        "real(8)",
        None,
        "Held-Suarez Rayleigh-drag time in days",
        value="1.0",
    ),
    VarDef(
        "hs_min_pressure_ratio",
        "real(8)",
        None,
        "Held-Suarez minimum pressure ratio used in log(p/p0)",
        value="1.0e-4",
    ),
    VarDef(
        "orographic_gwd_enabled",
        "logical",
        None,
        "Flag for the optional orographic gravity-wave-drag scheme",
        value=".false.",
    ),
    VarDef(
        "gwd_time_scale_days",
        "real(8)",
        None,
        "Reference drag time scale for the optional orographic gravity-wave-drag scheme",
        value="10.0",
    ),
    VarDef(
        "gwd_oro_threshold_m",
        "real(8)",
        None,
        "Orography threshold for the optional orographic gravity-wave-drag scheme",
        value="500.0",
    ),
    VarDef(
        "gwd_oro_scale_m",
        "real(8)",
        None,
        "Orography scale height used to saturate the optional orographic gravity-wave-drag scheme",
        value="1500.0",
    ),
    VarDef(
        "gwd_launch_sigma",
        "real(8)",
        None,
        "Sigma level above which the optional orographic gravity-wave-drag profile is applied",
        value="0.7",
    ),
    VarDef(
        "air_absortivity_co2",
        "real(8)",
        None,
        "Absorptivity of air in CO2 band",
        value="6.0",
    ),
    VarDef("flux_solar_in", "real(8)", "(ix, il)", "Flux of incoming solar radiation"),
    VarDef(
        "flux_ozone_lower",
        "real(8)",
        "(ix, il)",
        "Flux absorbed by ozone (lower stratosphere)",
    ),
    VarDef(
        "flux_ozone_upper",
        "real(8)",
        "(ix, il)",
        "Flux absorbed by ozone (upper stratosphere)",
    ),
    VarDef(
        "zenit_correction",
        "real(8)",
        "(ix, il)",
        "Zenith angle correction to (downward) absorptivity",
    ),
    VarDef(
        "stratospheric_correction",
        "real(8)",
        "(ix, il)",
        "Stratospheric correction for polar night",
    ),
    VarDef(
        "qcloud_equiv", "real(8)", "(ix, il)", " Equivalent specific humidity of clouds"
    ),
    ###################
    # Land model module
    ###################
    VarDef("rhcapl", "real(8)", "(ix, il)", "1/heat capacity (land)"),
    VarDef("cdland", "real(8)", "(ix, il)", " 1/dissipation time (land)"),
    VarDef(
        "stlcl_obs", "real(8)", "(ix, il)", "Climatological land surface temperature"
    ),
    VarDef(
        "snowdcl_obs",
        "real(8)",
        "(ix, il)",
        "Climatological snow depth (water equivalent)",
    ),
    VarDef(
        "soilwcl_obs", "real(8)", "(ix, il)", "Climatological soil water availability"
    ),
    VarDef("land_temp", "real(8)", "(ix, il)", "Land surface temperature"),
    VarDef("snow_depth", "real(8)", "(ix, il)", "Snow depth (water equivalent)"),
    VarDef("soil_avail_water", "real(8)", "(ix, il)", "Soil water availability", units="1"),
    VarDef("stl_lm", "real(8)", "(ix, il)", "Land-model surface temperature"),
    VarDef("fmask_land", "real(8)", "(ix, il)", "Fraction of land"),
    VarDef("bmask_land", "real(8)", "(ix, il)", " Binary land mask"),
    VarDef(
        "land_coupling_flag",
        "logical",
        None,
        "Flag for land-coupling (0: off, 1: on)",
        value=".true.",
    ),
    ###
    # Sea model
    ####
    VarDef("rhcaps", "real(p)", "(ix, il)", "1./heat_capacity (sea)"),
    VarDef("rhcapi", "real(p)", "(ix, il)", "1./heat_capacity (ice)"),
    VarDef("cdsea", "real(p)", "(ix, il)", "1./dissip_time (sea)"),
    VarDef("cdice", "real(p)", "(ix, il)", "1./dissip_time (ice)"),
    VarDef("fmask_sea", "real(p)", "(ix, il)", "Fraction of sea"),
    VarDef("bmask_sea", "real(p)", "(ix, il)", "Binary sea mask"),
    VarDef("deglat_s", "real(p)", "(il)", "Grid latitudes"),
    VarDef("hfseacl", "real(p)", "(ix, il)", "Annual-mean heat flux into sea sfc."),
    VarDef("sstom12", "real(p)", "(ix, il, 12)", "Ocean model SST climatology"),
    VarDef("sstcl_ob", "real(p)", "(ix, il)", "Observed clim. SST"),
    VarDef("sicecl_ob", "real(p)", "(ix, il)", "Clim. sea ice fraction"),
    VarDef("ticecl_ob", "real(p)", "(ix, il)", "Clim. sea ice temperature"),
    VarDef("sstan_ob", "real(p)", "(ix, il)", "Daily observed SST anomaly"),
    VarDef("sstcl_om", "real(p)", "(ix, il)", "Ocean model clim. SST"),
    VarDef("sst_am", "real(p)", "(ix, il)", "SST (full-field)"),
    VarDef("sstan_am", "real(p)", "(ix, il)", "SST anomaly"),
    VarDef("sice_am", "real(p)", "(ix, il)", "Sea ice fraction"),
    VarDef("tice_am", "real(p)", "(ix, il)", "Sea ice temperature"),
    VarDef("sst_om", "real(p)", "(ix, il)", "Ocean model SST"),
    VarDef("sice_om", "real(p)", "(ix, il)", "Model sea ice fraction"),
    VarDef("tice_om", "real(p)", "(ix, il)", "Model sea ice temperature"),
    VarDef("ssti_om", "real(p)", "(ix, il)", "Model SST + sea ice temp."),
    VarDef(
        "wsst_ob", "real(p)", "(ix, il)", "Weight for obs. SST anomaly in coupled runs"
    ),
    VarDef(
        "sst_anomaly_coupling_flag",
        "logical",
        None,
        "Weight for obs. SST anomaly in coupled runs",
        value=".true.",
    ),
    ###################
    # mod_radcon module
    ###################
    VarDef(
        "ablco2_ref", "real(8)", None, "Initial absorptivity of air in CO2 band (t=t0)"
    ),
    VarDef(
        "fband",
        "real(8)",
        "(100:400,4)",
        "Energy fraction emitted in each LW band = f(T)",
    ),
    VarDef(
        "alb_land",
        "real(8)",
        "(ix,il)",
        "Daily-mean albedo over land (bare-land + snow)",
    ),
    VarDef(
        "alb_sea",
        "real(8)",
        "(ix,il)",
        "Daily-mean albedo over sea  (open sea + sea ice)",
    ),
    VarDef("alb_surface", "real(8)", "(ix,il)", "Combined surface albedo (land + sea)"),
    VarDef("snowc", "real(8)", "(ix,il)", "Effective snow cover (fraction)"),
    VarDef(
        "rad_flux", "real(8)", "(ix,il,4)", "Radiative flux in different spectral bands"
    ),
    VarDef(
        "rad_tau2", "real(8)", "(ix,il,kx,4)", "Transmissivity of atmospheric layers"
    ),
    VarDef(
        "rad_st4a",
        "real(8)",
        "(ix,il,kx,2)",
        "Blackbody emission from full and half atmospheric levels",
    ),
    VarDef("rad_strat_corr", "real(8)", "(ix,il,2)", "Stratospheric correction term"),
    #############
    # Coordinates
    #############
    VarDef(
        "lon", "real", "(ix)", "longitude", units="degrees_east", std_name="longitude"
    ),
    VarDef(
        "lat", "real", "(il)", "latitude", units="degrees_north", std_name="latitude"
    ),
    VarDef(
        "lev",
        "real",
        "(kx)",
        "Vertical sigma coordinate",
        std_name="atmosphere_sigma_coordinate",
    ),
    #
    # Module instances
    #
    VarDef(
        "mod_geometry",
        "class(ModGeometry_t)",
        None,
        "Geometry module instance",
    ),
    VarDef(
        "mod_spectral",
        "class(ModSpectral_t)",
        None,
        "Spectral module instance",
    ),
    VarDef(
        "mod_implicit",
        "class(ModImplicit_t)",
        None,
        "Implicit module instance",
    ),
]


def build_fortran_sources():
    """Create the sources for the Python interface and the model state."""
    state_arrays = [
        var for var in model_state if var.dims and not var.is_module_instance
    ]
    state_scalars = [
        var for var in model_state if var.dims is None and not var.is_module_instance
    ]
    state_modules = [var for var in model_state if var.is_module_instance]

    file_loader = FileSystemLoader(THIS_FILE_DIR / "templates")
    env = Environment(loader=file_loader, trim_blocks=True, lstrip_blocks=True)
    template = env.get_template("model_state.f90.j2")
    output_file = str(SOURCES_DIR / "model_state.f90")
    template.stream(
        state_arrays=state_arrays,
        state_scalars=state_scalars,
        state_modules=state_modules,
    ).dump(output_file)
    print(f"Saved source: {output_file} file.")

    template = env.get_template("speedy_driver.f90.j2")
    output_file = str(SOURCES_DIR / "speedy_driver.f90")
    template.stream(state_arrays=state_arrays, state_scalars=state_scalars).dump(
        str(SOURCES_DIR / "speedy_driver.f90")
    )
    print(f"Saved source: {output_file} file.")


def export_model_state_html(file_path=None):
    """
    Export state variables description to HTML for docs.
    """
    file_loader = FileSystemLoader(THIS_FILE_DIR / "templates")
    env = Environment(loader=file_loader, trim_blocks=True, lstrip_blocks=True)
    template = env.get_template("model_state_def.html")
    if file_path is None:
        file_path = str(DOCS_DIR / "model_state_def.html")
    template.stream(model_state=model_state).dump(file_path)
    return file_path


def export_model_state_json():
    """Export state variables description in JSON format."""
    data2json = {
        var.name: dict(
            dtype=var.dtype,
            dims=var.dims,
            desc=var.desc,
            time_dim=var.time_dim,
            units=var.units,
            nc_dims=var.nc_dims,
            alt_name=var.alt_name,
            std_name=var.std_name,
        )
        for var in model_state
    }

    file_path = PYSPEEDY_DATA_DIR / "model_state.json"
    with open(file_path, "w") as outfile:
        json.dump(data2json, outfile, indent=4)

    print(f"Saved state definition: {file_path} file.")
    return file_path


def export_model_state_excel():
    """Export state variables description in Excel format."""
    if XLSXWRITER_IMPORTED and PANDAS_IMPORTED:
        _data = defaultdict(list)
        for var in model_state:
            _data["name"].append(var.name)
            _data["dtype"].append(var.dtype)
            _data["dims"].append(var.dims)
            _data["nc_dims"].append(var.dims)
            _data["desc"].append(var.desc)
            _data["units"].append(var.units)
            _data["time_dim"].append(var.time_dim)
            _data["alt_name"].append(var.alt_name)
            _data["std_name"].append(var.std_name)

        my_dataframe = pd.DataFrame(data=_data)

        file_path = PYSPEEDY_DATA_DIR / "model_state.xlsx"

        writer = pd.ExcelWriter(file_path, engine="xlsxwriter")
        sheetname = "state_variables"
        my_dataframe.to_excel(writer, sheet_name=sheetname, index=False)
        # Adjust the columns size
        worksheet = writer.sheets[sheetname]  # pull worksheet object
        for idx, col in enumerate(my_dataframe):  # loop through all columns
            series = my_dataframe[col]
            max_len = (
                max(
                    (
                        series.astype(str).map(len).max(),  # len of largest item
                        len(str(series.name)),  # len of column name/header
                    )
                )
                + 1
            )  # adding a little extra space
            worksheet.set_column(idx, idx, max_len)  # set column width
        writer.save()
        return file_path


if __name__ == "__main__":
    build_fortran_sources()
    export_model_state_json()
    export_model_state_excel()
