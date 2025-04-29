#!/usr/bin/env python3
"""
The MarsVars executable is for performing variable manipulations in
existing files. Most often, it is used to derive and add variables to
existing files, but it also differentiates variables with respect to
(w.r.t) the Z axis, column-integrates variables, converts aerosol
opacities from opacity per Pascal to opacity per meter, removes and
extracts variables from files, and enables scaling variables or editing
variable names, units, etc.

The executable requires:

    * ``[input_file]``           The file to be transformed

and optionally accepts:

    * ``[-add --add_variable]``          Derive and add variable to file
    * ``[-zdiff --differentiate_wrt_z]`` Differentiate variable w.r.t. Z axis
    * ``[-col --column_integrate]``      Column-integrate variable
    * ``[-zd --zonal_detrend]``          Subtract zonal mean from variable
    * ``[-to_dz --dp_to_dz]``            Convert aerosol opacity op/Pa -> op/m
    * ``[-to_dp --dz_to_dp]``            Convert aerosol opacity op/m -> op/Pa
    * ``[-rm --remove_variable]``        Remove variable from file
    * ``[-extract --extract_copy]``      Copy variable to new file
    * ``[-edit --edit_variable]``        Edit variable attributes or scale it

Third-party Requirements:

    * ``numpy``
    * ``netCDF4``
    * ``argparse``
    * ``os``
    * ``subprocess``
    * ``matplotlib``

"""

# Make print statements appear in color
from amescap.Script_utils import (
    Yellow, Cyan, Red, Nclr, Green, Blue
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import subprocess   # Run command-line commands
import warnings     # Suppress errors triggered by NaNs
import re           # Regular expressions
import matplotlib
import numpy as np
from netCDF4 import Dataset

# Force matplotlib NOT to load Xwindows backend
matplotlib.use("Agg")

# Load amesCAP modules
from amescap.FV3_utils import (
    fms_press_calc, fms_Z_calc, dvar_dh, cart_to_azimut_TR, mass_stream,
    zonal_detrend, spherical_div, spherical_curl, frontogenesis
)
from amescap.Script_utils import (
    check_file_tape, FV3_file_type, filter_vars,
    get_longname_unit, ak_bk_loader, except_message
)
from amescap.Ncdf_wrapper import Ncdf

# ======================================================
#                  DEFINITIONS
# ======================================================

import functools
import traceback

def debug_wrapper(func):
    """
    A decorator that wraps a function with error handling based on the 
    --debug flag.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global debug
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if debug:
                # In debug mode, show the full traceback
                print(f"{Red}ERROR in {func.__name__}: {str(e)}{Nclr}")
                traceback.print_exc()
            else:
                # In normal mode, show a clean error message
                print(f"{Red}ERROR in {func.__name__}: {str(e)}\nUse "
                      f"--debug for more information.{Nclr}")
            return 1  # Error exit code
    return wrapper

# List of supported variables for [-add --add_variable]
cap_str = " (derived w/CAP)"

master_list = {
    'curl': [
    "Relative vorticity",
        'Hz',
        ['ucomp', 'vcomp'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'div': [
        "Wind divergence",
        'Hz',
        ['ucomp', 'vcomp'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'DP': [
        "Layer thickness (P)",
        'Pa',
        ['ps', 'temp'],
        ['pfull']
    ],
    'dst_mass_mom': [
        "Dust MMR",
        'kg/kg',
        ['dzTau', 'temp'],
        ['pfull']
    ],
    'DZ': [
        "Layer thickness (Z)",
        'm',
        ['ps', 'temp'],
        ['pfull']
    ],
    'dzTau': [
        "Dust extinction rate",
        'km-1',
        ['dst_mass_mom', 'temp'],
        ['pfull']
    ],
    'fn': [
        "Frontogenesis",
        'K/m/s',
        ['ucomp', 'vcomp', 'theta'],
        ['pstd', 'zstd', 'zagl']
    ],
    'ice_mass_mom': [
        "Ice MMR",
        'kg/kg',
        ['izTau', 'temp'],
        ['pfull']
    ],
    'izTau': [
        "Ice extinction rate",
        'km-1',
        ['ice_mass_mom', 'temp'],
        ['pfull']
    ],
    'N': [
        "Brunt Vaisala freq.",
        'rad/s',
        ['ps', 'temp'],
        ['pfull']
    ],
    'pfull3D': [
        "Mid-layer pressure",
        'Pa',
        ['ps', 'temp'],
        ['pfull']
    ],
    'rho': [
        "Density",
        'kg/m^3',
        ['ps', 'temp'],
        ['pfull']
    ],
    'Ri': [
        "Richardson number",
        'none',
        ['ps', 'temp', 'ucomp', 'vcomp'],
        ['pfull']
    ],
    'scorer_wl': [
        "Scorer horiz. λ = 2π/√(l^2)",
        'm',
        ['ps', 'temp', 'ucomp'],
        ['pfull']
    ],
    'Tco2': [
        "CO2 condensation temperature",
        'K',
        ['ps', 'temp'],
        ['pfull', 'pstd']
    ],
    'theta': [
        "Potential temperature",
        'K',
        ['ps', 'temp'],
        ['pfull']
    ],
    'Vg_sed': [
        "Sedimentation rate",
        'm/s',
        ['dst_mass_mom', 'dst_num_mom', 'temp'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'w': [
        "vertical wind",
        'm/s',
        ['ps', 'temp', 'omega'],
        ['pfull']
    ],
    'w_net': [
        "Net vertical wind [w-Vg_sed]",
        'm/s',
        ['Vg_sed', 'w'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'wdir': [
        "Wind direction",
        'degree',
        ['ucomp', 'vcomp'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'wspeed': [
        "Wind speed",
        'm/s',
        ['ucomp', 'vcomp'],
        ['pfull', 'pstd', 'zstd', 'zagl']
    ],
    'zfull': [
        "Mid-layer altitude AGL",
        'm',
        ['ps', 'temp'],
        ['pfull']
    ],
    'ax': [
        "Zonal wave-mean flow forcing",
        'm/s^2',
        ['ucomp', 'w', 'rho'],
        ['pstd', 'zstd', 'zagl']
    ],
    'ay': [
        "Merid. wave-mean flow forcing",
        'm/s^2',
        ['vcomp', 'w', 'rho'],
        ['pstd', 'zstd', 'zagl']
    ],
    'ek': [
        "Wave kinetic energy",
        'J/kg',
        ['ucomp', 'vcomp'],
        ['pstd', 'zstd', 'zagl']
    ],
    'ep': [
        "Wave potential energy",
        'J/kg',
        ['temp'],
        ['pstd', 'zstd', 'zagl']
    ],
    'msf': [
        "Mass stream function",
        '1.e8 kg/s',
        ['vcomp'],
        ['pstd', 'zstd', 'zagl']
    ],
    'mx': [
        "Zonal momentum flux, vertical",
        'J/kg',
        ['ucomp', 'w'],
        ['pstd', 'zstd', 'zagl']
    ],
    'my': [
        "Merid. momentum flux, vertical",
        'J/kg',
        ['vcomp', 'w'],
        ['pstd', 'zstd', 'zagl']
    ],
    'tp_t': [
        "Normalized temperature perturbation",
        'None',
        ['temp'],
        ['pstd', 'zstd', 'zagl']
    ],
}

def add_help(var_list):
    help_text = (f"{'VARIABLE':9s} {'DESCRIPTION':33s} {'UNIT':11s} "
                 f"{'REQUIRED VARIABLES':24s} {'SUPPORTED FILETYPES'}"
                 f"\n{Cyan}")
    for var in var_list.keys():
        longname, unit, reqd_var, compat_files = var_list[var]
        reqd_var_fmt = ", ".join([f"{rv}" for rv in reqd_var])
        compat_file_fmt = ", ".join([f"{cf}" for cf in compat_files])
        help_text += (
            f"{var:9s} {longname:33s} {unit:11s} {reqd_var_fmt:24s} "
            f"{compat_file_fmt}\n"
            )
    return(help_text)

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    prog=('MarsVars'),
    description=(
        f"{Yellow}Enables derivations, manipulations, or removal of "
        f"variables from netCDF files.\n"
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('input_file', nargs='+',
    type=argparse.FileType('rb'),
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument('-add', '--add_variable', nargs='+', default=[],
    help=(
        f"Add a new variable to file. "
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"Variables that can be added are listed below.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -add rho\n{Yellow}\n"
        f"{add_help(master_list)}\n"
        f"{Yellow}NOTE: MarsVars offers some support on interpolated\n"
        f"files, particularly if ``pfull3D`` and ``zfull`` are added \n"
        f"to the file before interpolation.\n\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-zdiff', '--differentiate_wrt_z', nargs='+',
    default=[],
    help=(
        f"Differentiate a variable w.r.t. the Z axis.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"*Requires a variable with a vertical dimension*\n"
        f"A new variable ``d_dz_var`` in [Unit/m] will be added to the "
        f"file.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -zdiff dst_mass_mom\n"
        f"  {Blue}d_dz_dst_mass_mom is derived and added to the file"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-col', '--column_integrate', nargs='+', default=[],
    help=(
        f"Integrate a variable through the column.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"*Requires a variable with a vertical dimension*\n"
        f"A new variable (``var_col``) in [kg/m2] will be added to the "
        f"file.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -col dst_mass_mom\n"
        f"{Blue}(derive and add dst_mass_mom_col to the file)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-zd', '--zonal_detrend', nargs='+', default=[],
    help=(
        f"Detrend a variable by substracting its zonal mean value.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"A new a variable (``var_p``) (for prime) will be added to the"
        f" file.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -zd temp\n"
        f"{Blue}(temp_p is added to the file)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-to_dz', '--dp_to_dz', nargs='+', default=[],
    help=(
        f"Convert aerosol opacity [op/Pa] to [op/m]. "
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"Requires ``DP`` & ``DZ`` to be present in the file already.\n"
        f"A new variable (``[variable]_dp_to_dz``) is added to the "
        f"file.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -to_dz temp\n"
        f"{Blue}(temp_dp_to_dz is added to the file)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-to_dp', '--dz_to_dp', nargs='+', default=[],
    help=(
        f"Convert aerosol opacity [op/m] to [op/Pa]. "
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"Requires ``DP`` & ``DZ`` to be present in the file already.\n"
        f"A new variable (``[variable]_dz_to_dp``) is added to the "
        f"file.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -to_dp temp\n"
        f"{Blue}(temp_dz_to_dp is added to the file)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-rm', '--remove_variable', nargs='+', default=[],
    help=(
        f"Remove a variable from a file.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -rm ps"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-extract', '--extract_copy', nargs='+', default=[],
    help=(
        f"Copy one or more variables from a file into a new file of "
        f"the same name with the appended extension: '_extract'.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -extract ps temp\n"
        f"{Blue}(Creates 01336.atmos_average_extract.nc containing ps "
        f"and temp\nplus their dimensional variables [lat, "
        f"lon, lev, etc.])"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-edit', '--edit_variable', default=None,
    help=(
        f"Edit a variable's attributes or scale its values.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"Requires the use of one or more of the following flags:\n"
        f"``-rename``\n``-longname``\n``-unit``\n``-multiply``\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -edit ps -rename ps_mbar "
        f"-multiply 0.01 -longname 'Pressure in mb' -unit 'mbar'"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

# To be used jointly with --edit
parser.add_argument('-rename', '--rename', type=str, default=None,
    help=(
        f"Rename a variable. Requires ``-edit``.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -edit ps -rename ps_mbar\n"
        f"{Nclr}\n\n"
    )
)

# To be used jointly with --edit
parser.add_argument('-longname', '--longname', type=str, default=None,
    help=(
        f"Change a variable's 'longname' attribute. Requires ``-edit``.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -edit ps -longname "
        f"'Pressure scaled to mb'"
        f"{Nclr}\n\n"
    )
)

# To be used jointly with --edit
parser.add_argument('-unit', '--unit', type=str, default=None,
    help=(
        f"Change a variable's unit text. Requires ``-edit``.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -edit ps -unit 'mbar'"
        f"{Nclr}\n\n"
    )
)

# To be used jointly with --edit
parser.add_argument('-multiply', '--multiply', type=float, default=None,
    help=(
        f"Scale a variable's values. Requires ``-edit``.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -edit ps -multiply 0.01"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug',  action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsVars 01336.atmos_average.nc -add rho --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()
debug = args.debug

if args.input_file:
    for file in args.input_file:
        if not re.search(".nc", file.name):
            parser.error(f"{Red}{file.name} is not a netCDF file{Nclr}")
            exit()

if args.rename and args.edit_variable is None:
    parser.error(f"{Red}The -rename argument requires -edit to be used "
                 f"with it (e.g., MarsVars 01336.atmos_average.nc "
                 f"-edit ps -rename ps_mbar)"
                 f"{Nclr}")
    exit()

if args.longname and args.edit_variable is None:
    parser.error(f"{Red}The -longname argument requires -edit to be "
                 f"used with it (e.g., MarsVars 01336.atmos_average.nc "
                 f"-edit ps -longname 'Pressure scaled to mb')"
                 f"{Nclr}")
    exit()

if args.unit and args.edit_variable is None:
    parser.error(f"{Red}The -unit argument requires -edit to be used "
                 f"with it (e.g., MarsVars 01336.atmos_average.nc "
                 f"-edit ps -unit 'mbar')"
                 f"{Nclr}")
    exit()

if args.multiply and args.edit_variable is None:
    parser.error(f"{Red}The -multiply argument requires -edit to be "
                 f"used with it (e.g., MarsVars 01336.atmos_average.nc "
                 f"-edit ps -multiply 0.01)"
                 f"{Nclr}")
    exit()

# ======================================================================
# TODO : If only one timestep, reshape from
#       (lev, lat, lon) to (t, lev, lat, lon)

# Fill values for NaN. np.NaN, raises errors when running runpinterp
fill_value = 0.

# Define constants
global rgas, psrf, Tpole, g, R, Rd, rho_air, rho_dst, rho_ice
global Qext_dst, Qext_ice, n0, S0, T0, Cp, Na, amu, amu_co2, mass_co2
global sigma, M_co2, N, C_dst, C_ice

rgas = 189.  # Gas const. CO2 [J/kg/K or m^2/s^2/K]
psrf = 610.  # Mars surface pressure [Pa or kg/m/s^2]
Tpole = 150.  # Polar temperature [K]
g = 3.72  # Gravitational constant for Mars [m/s^2]
R = 8.314  # Universal gas constant [J/mol/K]
Rd = 192.0  # R for dry air on Mars [J/kg/K]
rho_air = psrf/(rgas*Tpole)  # Air density (ρ) [kg/m^3]
rho_dst = 2500.  # Dust particle ρ [kg/m^3]
# rho_dst = 3000  # Dust particle ρ [kg/m^3] (Kleinbohl, 2009)
rho_ice = 900  # Ice particle ρ [kg/m^3] (Heavens, 2010)
Qext_dst = 0.35  # Dust extinction efficiency (MCS) (Kleinbohl, 2009)
Qext_ice = 0.773  # Ice extinction efficiency (MCS) (Heavens, 2010)
Reff_dst = 1.06  # Effective dust particle radius [µm] (Kleinbohl, 2009)
Reff_ice = 1.41  # Effective ice particle radius [µm] (Heavens, 2010)
n0 = 1.37*1.e-5  # Sutherland's law [N-s/m^2]
S0 = 222  # Sutherland's law [K]
T0 = 273.15  # Sutherland's law [K]
Cp = 735.0  # [J/K]
Na = 6.022*1.e23  # Avogadro's number [per mol]
Kb = R/Na  # Boltzmann constant [m^2*kg/s^2/K]
amu = 1.66054*1.e-27  # Atomic mass Unit [kg/amu]
amu_co2 = 44.0  # Molecular mass of CO2 [amu]
mass_co2 = amu_co2*amu  # Mass of 1 CO2 particle [kg]
sigma = 0.63676  # Gives effective variance = 0.5 (Dust)
M_co2 = 0.044  # Molar mass of CO2 [kg/mol]
N = 0.01  # For wave potential energy calc. [rad/s]

# For mmr <-> extinction rate calculations:
C_dst = (4/3) * (rho_dst/Qext_dst) * Reff_dst # = 12114.286 [m-2]
C_ice = (4/3) * (rho_ice/Qext_ice) * Reff_ice # = 2188.874 [m-2]

# ===========================

def compute_p_3D(ps, ak, bk, shape_out):
    """
    Compute the 3D pressure at layer midpoints.

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param ak: Vertical coordinate pressure value (Pa)
    :type ak: array [phalf]

    :param bk: Vertical coordinate sigma value (None)
    :type bk: array [phalf]

    :param shape_out: Determines how to handle the dimensions of p_3D.
        If ``len(time) = 1`` (one timestep), ``p_3D`` is returned as
        [1, lev, lat, lon] as opposed to [lev, lat, lon]
    :type shape_out: float

    :raises:

    :return: ``p_3D`` The full 3D pressure array (Pa)
    :rtype: array [time, lev, lat, lon]

    """
    p_3D = fms_press_calc(ps, ak, bk, lev_type="full")
    # Swap dimensions 0 and 1 (time and lev)
    p_3D = p_3D.transpose(lev_T)
    return p_3D.reshape(shape_out)

# =====================================================================
def compute_rho(p_3D, temp):
    """
    Compute density.

    :param p_3D: Pressure (Pa)
    :type p_3D: array [time, lev, lat, lon]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :raises:

    :return: Density (kg/m^3)
    :rtype: array [time, lev, lat, lon]

    """
    return p_3D / (rgas*temp)

# =====================================================================
def compute_xzTau(q, temp, lev, const, f_type):
    """
    Compute the dust or ice extinction rate.
    Adapted from Heavens et al. (2011) observations from MCS (JGR).

    :param q: Dust or ice mass mixing ratio (ppm)
    :type q: array [time, lev, lat, lon]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
    :type lev: array [lev]

    :param const: Dust or ice constant
    :type const: array

    :param f_type: The FV3 file type: diurn, daily, or average
    :type f_stype: str

    :raises:

    :return: ``xzTau`` Dust or ice extinction rate (km-1)
    :rtype: array [time, lev, lat, lon]

    """
    if f_type == "diurn":
        # Handle diurn files 
        PT = np.repeat(
            lev,
            (q.shape[0] * q.shape[1] * q.shape[3] * q.shape[4])
        )
        PT = np.reshape(
            PT,
            (q.shape[2], q.shape[0], q.shape[1], q.shape[3], q.shape[4])
        )
        # (lev, tim, tod, lat, lon) -> (tim, tod, lev, lat, lon)
        P = PT.transpose((1, 2, 0, 3, 4))
    else:
        # For average and daily files, ensure proper broadcasting across all times
        # Create a properly sized pressure field with correct time dimension
        P = np.zeros_like(q)
        
        # Fill P with the appropriate pressure level for each vertical index
        for z in range(len(lev)):
            if len(q.shape) == 4:  # Standard [time, lev, lat, lon] format
                P[:, z, :, :] = lev[z]
            else:
                # Handle other shapes appropriately
                P[..., z, :, :] = lev[z]
                
        # PT = np.repeat(lev, (q.shape[0] * q.shape[2] * q.shape[3]))
        # PT = np.reshape(
        #     PT,
        #     (q.shape[1], q.shape[0], q.shape[2], q.shape[3])
        # )
        # # Swap dimensions 0 and 1 (time and lev)
        # P = PT.transpose(lev_T)

    rho_z = P / (Rd*temp)
    # Converts mass mixing ratio (q) from kg/kg -> ppm (mg/kg)
    # Converts extinction (xzTau) from m-1 -> km-1
    xzTau = (rho_z * (q*1.e6)/const) * 1000
    return xzTau

# =====================================================================
def compute_mmr(xTau, temp, lev, const, f_type):
    """
    Compute the dust or ice mixing ratio.
    Adapted from Heavens et al. (2011) observations from MCS (JGR).

    :param xTau: Dust or ice extinction rate (km-1)
    :type xTau: array [time, lev, lat, lon]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :param lev: Vertical coordinate (e.g., pstd) (e.g., Pa)
    :type lev: array [lev]

    :param const: Dust or ice constant
    :type const: array

    :param f_type: The FV3 file type: diurn, daily, or average
    :type f_stype: str

    :raises:

    :return: ``q``, Dust or ice mass mixing ratio (ppm)
    :rtype: array [time, lev, lat, lon]

    """
    if f_type == "diurn":
        # Handle diurnal files
        PT = np.repeat(lev,(xTau.shape[0] * xTau.shape[1]
                            * xTau.shape[3] * xTau.shape[4]))
        PT = np.reshape(PT, (xTau.shape[2], xTau.shape[0],
                             xTau.shape[1], xTau.shape[3],
                             xTau.shape[4]))
        # (lev, tim, tod, lat, lon) -> (tim, tod, lev, lat, lon)
        P = PT.transpose((1, 2, 0, 3, 4))
    else:
        # For average and daily files, create properly broadcast pressure array
        P = np.zeros_like(xTau)
        
        # Fill P with the appropriate pressure level for each vertical index
        for z in range(len(lev)):
            if len(xTau.shape) == 4:  # Standard [time, lev, lat, lon] format
                P[:, z, :, :] = lev[z]
            else:
                # Handle other shapes appropriately
                P[..., z, :, :] = lev[z]
        # PT = np.repeat(lev, (xTau.shape[0] * xTau.shape[2]
        #                      * xTau.shape[3]))
        # PT = np.reshape(PT,(xTau.shape[1], xTau.shape[0],
        #                     xTau.shape[2], xTau.shape[3]))
        # # Swap dimensions 0 and 1 (time and lev)
        # P = PT.transpose(lev_T)

    rho_z = P / (Rd*temp)
    # Converts extinction (xzTau) from km-1 -> m-1
    # Converts mass mixing ratio (q) from ppm (kg/kg) -> mg/kg
    q = (const * (xTau/1000) / rho_z) / 1.e6
    return q

# =====================================================================
def compute_Vg_sed(xTau, nTau, temp):
    """
    Calculate the sedimentation rate of the dust.

    :param xTau: Dust or ice MASS mixing ratio (ppm)
    :type xTau: array [time, lev, lat, lon]

    :param nTau: Dust or ice NUMBER mixing ratio (None)
    :type nTau: array [time, lev, lat, lon]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :raises:

    :return: ``Vg`` Dust sedimentation rate (m/s)
    :rtype: array [time, lev, lat, lon]

    """
    r0 = (((3.*xTau) / (4.*np.pi*rho_dst*nTau)) ** (1/3)
          * np.exp(-3 * (sigma**2) / 2))
    Rp = r0 * np.exp(3.5*sigma**2)
    c = (2/9) * rho_dst * (Rp)**2 * g
    eta = n0 * ((temp/T0)**(3/2)) * ((T0+S0)/(temp+S0))
    v = np.sqrt((3*Kb*temp) / mass_co2)
    mfp = (2*eta) / (rho_air*v)
    Kn = mfp / Rp
    alpha = 1.246 + 0.42*np.exp(-0.87/Kn)
    Vg = c*(1 + alpha*Kn)/eta
    return Vg

# =====================================================================
def compute_w_net(Vg, wvar):
    """
    Computes the net vertical wind, which is the vertical wind (w)
    minus the sedimentation rate (``Vg_sed``)::

        w_net = w - Vg_sed

    :param Vg: Dust sedimentation rate (m/s)
    :type Vg: array [time, lev, lat, lon]

    :param wvar: Vertical wind (m/s)
    :type wvar: array [time, lev, lat, lon]

    :raises:

    :return: `w_net` Net vertical wind speed (m/s)
    :rtype: array [time, lev, lat, lon]

    """
    w_net = np.subtract(wvar, Vg)
    return w_net

# =====================================================================
def compute_theta(p_3D, ps, temp, f_type):
    """
    Compute the potential temperature.

    :param p_3D: The full 3D pressure array (Pa)
    :type p_3D: array [time, lev, lat, lon]

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :param f_type: The FV3 file type: diurn, daily, or average
    :type f_type: str

    :raises:

    :return: Potential temperature (K)
    :rtype: array [time, lev, lat, lon]

    """
    theta_exp = R / (M_co2*Cp)
    # Broadcast dimensions
    if f_type == "diurn":
        # (time, tod, lat, lon) -> (time, tod, 1, lat, lon)
        ps_shape = [ps.shape[0], ps.shape[1], 1, ps.shape[2],
                    ps.shape[3]]
    else:
        # (time, lat, lon) -> (time, 1, lat, lon)
        ps_shape = [ps.shape[0], 1, ps.shape[1], ps.shape[2]]

    return temp * (np.reshape(ps, ps_shape)/p_3D) ** theta_exp

# =====================================================================
def compute_w(rho, omega):
    """
    Compute the vertical wind using the omega equation.

    Under hydrostatic balance, omega is proportional to the vertical
    wind velocity (``w``)::

        omega = dp/dt = (dp/dz)(dz/dt) = (dp/dz) * w

    Under hydrostatic equilibrium::

        dp/dz = -rho * g

    So ``omega`` can be calculated as::

        omega = -rho * g * w

    :param rho: Atmospheric density (kg/m^3)
    :type rho: array [time, lev, lat, lon]

    :param omega: Rate of change in pressure at layer midpoint (Pa/s)
    :type omega: array [time, lev, lat, lon]

    :raises:

    :return: vertical wind (m/s)
    :rtype: array [time, lev, lat, lon]

    """
    return -omega / (rho*g)

# =====================================================================
def compute_zfull(ps, ak, bk, temp):
    """
    Calculate the altitude of the layer midpoints above ground level.

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param ak: Vertical coordinate pressure value (Pa)
    :type ak: array [phalf]

    :param bk: Vertical coordinate sigma value (None)
    :type bk: array [phalf]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :raises:

    :return: ``zfull`` (m)
    :rtype: array [time, lev, lat, lon]

    """
    zfull = fms_Z_calc(ps, ak, bk, temp.transpose(lev_T),
                       topo=0., lev_type="full")

    # .. note:: lev_T swaps dims 0 & 1, ensuring level is the first
    # dimension for the calculation

    zfull = zfull.transpose(lev_T_out)
    return zfull

# =====================================================================
def compute_zhalf(ps, ak, bk, temp):
    """
    Calculate the altitude of the layer interfaces above ground level.

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param ak: Vertical coordinate pressure value (Pa)
    :type ak: array [phalf]

    :param bk: Vertical coordinate sigma value (None)
    :type bk: array [phalf]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :raises:

    :return: ``zhalf`` (m)
    :rtype: array [time, lev, lat, lon]

    """
    zhalf = fms_Z_calc(ps, ak, bk, temp.transpose(lev_T),
                       topo=0., lev_type="half")

    # .. note:: lev_T swaps dims 0 & 1, ensuring level is the first
    # dimension for the calculation

    zhalf = zhalf.transpose(lev_T_out)
    return zhalf

# =====================================================================
def compute_DZ_full_pstd(pstd, temp, ftype="average"):
    """
    Calculate the thickness of a layer from the midpoint of the
    standard pressure levels (``pstd``).

    In this context, ``pfull=pstd`` with the layer interfaces
    defined somewhere in between successive layers::

        --- Nk --- TOP       ========  phalf
        --- Nk-1 ---
                             --------  pfull = pstd    ^
                                                       | DZ_full_pstd
                             ========  phalf           |
        --- 1 ---            --------  pfull = pstd    v
        --- 0 --- SFC        ========  phalf
                              / / / /

    :param pstd: Vertical coordinate (pstd; Pa)
    :type pstd: array [lev]

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :param f_type: The FV3 file type: diurn, daily, or average
    :type f_stype: str

    :raises:

    :return: DZ_full_pstd, Layer thicknesses (Pa)
    :rtype: array [time, lev, lat, lon]

    """
    # Determine whether the lev dimension is located at i = 1 or i = 2
    if ftype == "diurn":
        axis = 2
    else:
        axis = 1

    # Make lev the first dimension, swapping it with time
    temp = np.swapaxes(temp, 0, axis)

    # Create a new shape = [1, 1, 1, 1]
    new_shape = [1 for i in range(0, len(temp.shape))]

    # Make the first dimesion = the length of the lev dimension (pstd)
    new_shape[0] = len(pstd)

    # Reshape pstd according to new_shape
    pstd_reshaped = pstd.reshape(new_shape)

    # Ensure pstd is broadcast to match the shape of temp 
    # (along non-level dimensions)
    broadcast_shape = list(temp.shape)
    broadcast_shape[0] = len(pstd)  # Keep level dimension the same
    pstd_broadcast = np.broadcast_to(pstd_reshaped, broadcast_shape)
    
    # Compute thicknesses using avg. temperature of both layers
    DZ_full_pstd = np.zeros_like(temp)
    DZ_full_pstd[0:-1, ...] = (
        -rgas * 0.5 * (temp[1:, ...]+temp[0:-1, ...]) / g
        * np.log(pstd_broadcast[1:, ...]/pstd_broadcast[0:-1, ...])
    )

    # There is nothing to differentiate the last layer with, so copy
    # the second-to-last layer.
    DZ_full_pstd[-1, ...] = DZ_full_pstd[-2, ...]

    # .. note:: that unless you fine-tune the standard pressure levels to
    # match the model top, data is usually missing in the last few
    # layers.

    return np.swapaxes(DZ_full_pstd, 0, axis)

# =====================================================================
def compute_N(theta, zfull):
    """
    Calculate the Brunt Vaisala freqency.

    :param theta: Potential temperature (K)
    :type theta: array [time, lev, lat, lon]

    :param zfull: Altitude above ground level at the layer midpoint (m)
    :type zfull: array [time, lev, lat, lon]

    :raises:

    :return: ``N``, Brunt Vaisala freqency [rad/s]
    :rtype: array [time, lev, lat, lon]

    """
    # Differentiate theta w.r.t. zfull to obdain d(theta)/dz
    dtheta_dz = dvar_dh(theta.transpose(lev_T),
                        zfull.transpose(lev_T)).transpose(lev_T)

    # .. note:: lev_T swaps dims 0 & 1, ensuring level is the first
    # dimension for the differentiation

    # Calculate the Brunt Vaisala frequency
    N = np.sqrt(g/theta * dtheta_dz)

    return N

# =====================================================================
def compute_Tco2(P_3D):
    """
    Calculate the frost point of CO2.
    Adapted from Fannale (1982) - Mars: The regolith-atmosphere cap
    system and climate change. Icarus.

    :param P_3D: The full 3D pressure array (Pa)
    :type p_3D: array [time, lev, lat, lon]

    :raises:

    :return: CO2 frost point [K]
    :rtype: array [time, lev, lat, lon]

    """
    # Set some constants
    B = -3167.8 # K
    CO2_triple_pt_P = 518000 # Pa

    # Determine where the pressure < the CO2 triple point pressure
    condition = (P_3D < CO2_triple_pt_P)

    # If P < triple point, calculate temperature
    # modified vapor pressure curve equation
    temp_where_true = B/(np.log(0.01*P_3D) - 23.23)

    # If P > triple point, calculate temperature
    temp_where_false = 684.2 - 92.3*np.log(P_3D) + 4.32*np.log(P_3D)**2

    return np.where(condition, temp_where_true, temp_where_false)

# =====================================================================
def compute_scorer(N, ucomp, zfull):
    """
    Calculate the Scorer wavelength.

    :param N: Brunt Vaisala freqency (rad/s)
    :type N: float [time, lev, lat, lon]

    :param ucomp: Zonal wind (m/s)
    :type ucomp: array [time, lev, lat, lon]

    :param zfull: Altitude above ground level at the layer midpoint (m)
    :type zfull: array [time, lev, lat, lon]

    :raises:

    :return: ``scorer_wl`` Scorer horizontal wavelength (m)
    :rtype: array [time, lev, lat, lon]

    """
    # Differentiate U w.r.t. zfull TWICE to obdain d^2U/dz^2
    dUdz = dvar_dh(ucomp.transpose(lev_T),
                   zfull.transpose(lev_T)).transpose(lev_T)
    dUdz2 = dvar_dh(dUdz.transpose(lev_T),
                    zfull.transpose(lev_T)).transpose(lev_T)

    # .. note:: lev_T swaps dims 0 & 1, ensuring level is the first
    # dimension for the differentiation

    # Compute the scorer parameter I^2(z) (m-1)
    scorer_param = N**2/ucomp**2 - dUdz2/ucomp

    # Compute the wavelength
    # I = sqrt(I^2) = wavenumber (k)
    # wavelength (lambda) = 2pi/k
    scorer_wl = 2*np.pi/np.sqrt(scorer_param)

    return scorer_wl

# =====================================================================
def compute_DP_3D(ps, ak, bk, shape_out):
    """
    Calculate the thickness of a layer in pressure units.

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param ak: Vertical coordinate pressure value (Pa)
    :type ak: array [phalf]

    :param bk: Vertical coordinate sigma value (None)
    :type bk: array [phalf]

    :param shape_out: Determines how to handle the dimensions of DP_3D.
        If len(time) = 1 (one timestep), DP_3D is returned as
        [1, lev, lat, lon] as opposed to [lev, lat, lon]
    :type shape_out: float

    :raises:

    :return: ``DP`` Layer thickness in pressure units (Pa)
    :rtype: array [time, lev, lat, lon]

    """
    # Get the 3D pressure field from fms_press_calc
    p_half3D = fms_press_calc(ps, ak, bk, lev_type="half")
    # fms_press_calc will swap dimensions 0 and 1 so p_half3D has
    # dimensions = [lev, t, lat, lon]
    # Calculate the differences in pressure between each layer midpoint
    DP_3D = p_half3D[1:, ..., ] - p_half3D[0:-1, ...]

    # Swap dimensions 0 and 1, back to [t, lev, lat, lon]
    DP_3D = DP_3D.transpose(lev_T)

    DP = DP_3D.reshape(shape_out)
    return DP

# =====================================================================
def compute_DZ_3D(ps, ak, bk, temp, shape_out):
    """
    Calculate the thickness of a layer in altitude units.

    :param ps: Surface pressure (Pa)
    :type ps: array [time, lat, lon]

    :param ak: Vertical coordinate pressure value (Pa)
    :type ak: array [phalf]

    :param bk: Vertical coordinate sigma value (None)
    :type bk: array [phalf]

    :param shape_out: Determines how to handle the dimensions of DZ_3D.
        If len(time) = 1 (one timestep), DZ_3D is returned as
        [1, lev, lat, lon] as opposed to [lev, lat, lon]
    :type shape_out: float

    :raises:

    :return: ``DZ`` Layer thickness in altitude units (m)
    :rtype: array [time, lev, lat, lon]

    """

    # Get the 3D altitude field from fms_Z_calc
    z_half3D = fms_Z_calc(ps, ak, bk, temp.transpose(lev_T), topo=0.,
                          lev_type="half")
    # fms_press_calc will swap dimensions 0 and 1 so p_half3D has
    # dimensions = [lev, t, lat, lon]

    # Calculate the differences in pressure between each layer midpoint
    DZ_3D = z_half3D[0:-1, ...]-z_half3D[1:, ..., ]
    # .. note:: the reversed order: Z decreases with increasing levels

    # Swap dimensions 0 and 1, back to [t, lev, lat, lon]
    DZ_3D = DZ_3D.transpose(lev_T)

    DZ = DZ_3D.reshape(shape_out)

    return DZ

# =====================================================================
def compute_Ep(temp):
    """
    Calculate wave potential energy::

        Ep = 1/2 (g/N)^2 (temp'/temp)^2

    :param temp: Temperature (K)
    :type temp: array [time, lev, lat, lon]

    :raises:

    :return: ``Ep`` Wave potential energy (J/kg)
    :rtype: array [time, lev, lat, lon]

    """
    return 0.5*g**2*(zonal_detrend(temp)/(temp*N))**2

# =====================================================================
def compute_Ek(ucomp, vcomp):
    """
    Calculate wave kinetic energ::

        Ek = 1/2 (u'**2+v'**2)

    :param ucomp: Zonal wind (m/s)
    :type ucomp: array [time, lev, lat, lon]

    :param vcomp: Meridional wind (m/s)
    :type vcomp: array [time, lev, lat, lon]

    :raises:

    :return: ``Ek`` Wave kinetic energy (J/kg)
    :rtype: array [time, lev, lat, lon]

    """
    return 0.5*(zonal_detrend(ucomp)**2+zonal_detrend(vcomp)**2)

# =====================================================================
def compute_MF(UVcomp, w):
    """
    Calculate zonal or meridional momentum fluxes.

    :param UVcomp: Zonal or meridional wind (ucomp or vcomp)(m/s)
    :type UVcomp: array

    :param w: Vertical wind (m/s)
    :type w: array [time, lev, lat, lon]

    :raises:

    :return: ``u'w'`` or ``v'w'``, Zonal/meridional momentum flux (J/kg)
    :rtype: array [time, lev, lat, lon]

    """
    return zonal_detrend(UVcomp)*zonal_detrend(w)

# =====================================================================
def compute_WMFF(MF, rho, lev, interp_type):
    """
    Calculate the zonal or meridional wave-mean flow forcing::

        ax = -1/rho d(rho u'w')/dz
        ay = -1/rho d(rho v'w')/dz

    If interp_type == ``pstd``, then::

        [du/dz = (du/dp).(dp/dz)] > [du/dz = -rho*g * (du/dp)]

    where::

        dp/dz = -rho*g
        [du/dz = (du/dp).(-rho*g)] > [du/dz = -rho*g * (du/dp)]

    :param MF: Zonal/meridional momentum flux (J/kg)
    :type MF: array [time, lev, lat, lon]

    :param rho: Atmospheric density (kg/m^3)
    :type rho: array [time, lev, lat, lon]

    :param lev: Array for the vertical grid (zagl, zstd, pstd, or pfull)
    :type lev: array [lev]

    :param interp_type: The vertical grid type (``zagl``, ``zstd``,
        ``pstd``, or ``pfull``)
    :type interp_type: str

    :raises:

    :return: The zonal or meridional wave-mean flow forcing (m/s2)
    :rtype: array [time, lev, lat, lon]

    """
    # Differentiate the momentum flux (MF)
    darr_dz = dvar_dh((rho*MF).transpose(lev_T), lev).transpose(lev_T)
    # Manually swap dimensions 0 and 1 so lev_T has lev for first
    # dimension [lev, t, lat, lon] for the differentiation

    if interp_type == "pstd":
        # Computed du/dp, need to multiply by (-rho g) to obtain du/dz
        return g * darr_dz
    else:
        # zagl and zstd grids have levels in meters, so du/dz
        # is not multiplied by g.
        return -1/rho*darr_dz


# =====================================================================
def check_dependencies(f, var, master_list, add_missing=True, dependency_chain=None):
    """
    Check if all dependencies for a variable are present in the file,
    and optionally try to add missing dependencies.
    
    :param f: NetCDF file object
    :param var: Variable to check dependencies for
    :param master_list: Dictionary of supported variables and their dependencies
    :param add_missing: Whether to try adding missing dependencies (default: True)
    :param dependency_chain: List of variables in the current dependency chain (for detecting cycles)
    :return: True if all dependencies are present or successfully added, False otherwise
    """
    # Initialize dependency chain if None
    if dependency_chain is None:
        dependency_chain = []
    
    # Check if we're in a circular dependency
    if var in dependency_chain:
        print(f"{Red}Circular dependency detected: {' -> '.join(dependency_chain + [var])}{Nclr}")
        return False
    
    # Add current variable to dependency chain
    dependency_chain = dependency_chain + [var]
    
    if var not in master_list:
        print(f"{Red}Variable `{var}` is not in the master list of supported variables.{Nclr}")
        return False
    
    # Get the list of required variables for this variable
    required_vars = master_list[var][2]
    
    # Check each required variable
    missing_vars = []
    for req_var in required_vars:
        if req_var not in f.variables:
            missing_vars.append(req_var)
    
    if not missing_vars:
        # All dependencies are present
        return True
    
    if not add_missing:
        # Dependencies are missing but we're not adding them
        dependency_list = ", ".join(missing_vars)
        print(f"{Red}Missing dependencies for {var}: {dependency_list}{Nclr}")
        return False
    
    # Try to add missing dependencies
    successfully_added = []
    failed_to_add = []
    
    for missing_var in missing_vars:
        # Check if we can add this dependency (must be in master_list)
        if missing_var in master_list:
            # Recursively check dependencies for this variable, passing the current dependency chain
            if check_dependencies(f, missing_var, master_list, add_missing=True, dependency_chain=dependency_chain):
                # If dependencies are satisfied, try to add the variable
                try:
                    print(f"{Yellow}Dependency {missing_var} for {var} can be added{Nclr}")
                    # Get the file type and interpolation type
                    f_type, interp_type = FV3_file_type(f)
                    
                    # Check if the interpolation type is compatible with this variable
                    compat_file_fmt = ", ".join([cf for cf in master_list[missing_var][3]])
                    if interp_type not in master_list[missing_var][3]:
                        print(f"{Red}Cannot add {missing_var}: incompatible file type {interp_type}{Nclr}")
                        failed_to_add.append(missing_var)
                        continue
                    
                    # Mark it as successfully checked
                    successfully_added.append(missing_var)
                    
                except Exception as e:
                    print(f"{Red}Error checking dependency {missing_var}: {str(e)}{Nclr}")
                    failed_to_add.append(missing_var)
            else:
                # If dependencies for this dependency are not satisfied
                failed_to_add.append(missing_var)
        else:
            # This dependency is not in the master list, cannot be added
            print(f"{Red}Dependency {missing_var} for {var} is not in the master list and cannot be added automatically{Nclr}")
            failed_to_add.append(missing_var)
    
    # Check if all dependencies were added
    if not failed_to_add:
        return True
    else:
        # Some dependencies could not be added
        dependency_list = ", ".join(failed_to_add)
        print(f"{Red}Cannot add {var}: missing dependencies {dependency_list}{Nclr}")
        return False


# =====================================================================
def check_variable_exists(var_name, file_vars):
    """Check if a variable exists in the file, considering alternative naming conventions"""
    if var_name in file_vars:
        return True
        
    # Handle _micro/_mom naming variations
    if var_name.endswith('_micro'):
        alt_name = var_name.replace('_micro', '_mom')
        if alt_name in file_vars:
            return True
    elif var_name.endswith('_mom'):
        alt_name = var_name.replace('_mom', '_micro')
        if alt_name in file_vars:
            return True
    elif var_name == 'dst_num_mom':
        # Special case for dst_num_mom/dst_num_micro
        if 'dst_num_micro' in file_vars:
            return True
            
    return False


# =====================================================================
def get_existing_var_name(var_name, file_vars):
    """Get the actual variable name that exists in the file"""
    if var_name in file_vars:
        return var_name
        
    # Check alternative names
    if var_name.endswith('_micro'):
        alt_name = var_name.replace('_micro', '_mom')
        if alt_name in file_vars:
            return alt_name
    elif var_name.endswith('_mom'):
        alt_name = var_name.replace('_mom', '_micro')
        if alt_name in file_vars:
            return alt_name
    elif var_name == 'dst_num_mom':
        # Special case for dst_num_mom/dst_num_micro
        if 'dst_num_micro' in file_vars:
            return 'dst_num_micro'
            
    return var_name  # Return original if no match found


# =====================================================================
def process_add_variables(ifile, add_list, master_list, debug=False):
    """
    Process the list of variables to add, handling dependencies appropriately.
    
    :param ifile: Input file path
    :param add_list: List of variables to add
    :param master_list: Dictionary of supported variables and their dependencies
    :param debug: Whether to show debug information
    """
        # Create a topologically sorted list of variables to add
    variables_to_add = []
    already_in_file = []
    
    # First check if all requested variables already exist in the file
    with Dataset(ifile, "r", format="NETCDF4_CLASSIC") as f:
        file_vars = set(f.variables.keys())
        
        # Check if all requested variables are already in the file
        for var in add_list:
            if check_variable_exists(var, file_vars):
                existing_name = get_existing_var_name(var, file_vars)
                already_in_file.append((var, existing_name))
        
        # If all requested variables exist, report and exit
        if len(already_in_file) == len(add_list):
            if len(add_list) == 1:
                var, actual_var = already_in_file[0]
                if var == actual_var:
                    print(f"{Yellow}Variable '{var}' is already in the file.{Nclr}")
                else:
                    print(f"{Yellow}Variable '{var}' is already in the file (as '{actual_var}').{Nclr}")
            else:
                print(f"{Yellow}All requested variables are already in the file:{Nclr}")
                for var, actual_var in already_in_file:
                    if var == actual_var:
                        print(f"{Yellow}  - {var}{Nclr}")
                    else:
                        print(f"{Yellow}  - {var} (as '{actual_var}'){Nclr}")
            return
    
    def add_with_dependencies(var):
        """Helper function to add a variable and its dependencies to the list"""
        # Skip if already processed
        if var in variables_to_add:
            return True
        
        # Open the file to check dependencies
        with Dataset(ifile, "a", format="NETCDF4_CLASSIC") as f:
            file_vars = set(f.variables.keys())
            
            # Skip if already in file
            if check_variable_exists(var, file_vars):
                return True

            f_type, interp_type = FV3_file_type(f)
                
            # Check file compatibility
            if interp_type not in master_list[var][3]:
                compat_file_fmt = ", ".join(master_list[var][3])
                print(f"{Red}ERROR: Variable '{Yellow}{var}{Red}' can only be added to file type: {Yellow}{compat_file_fmt}{Nclr}")
                return False
                
            # Check each dependency
            all_deps_ok = True
            for dep in master_list[var][2]:
                # Skip if already in file (including alternative names)
                if check_variable_exists(dep, file_vars):
                    continue
                    
                # If dependency can be added, try to add it
                if dep in master_list:
                    if not add_with_dependencies(dep):
                        all_deps_ok = False
                        print(f"{Red}Cannot add {var}: Required dependency {dep} cannot be added{Nclr}")
                else:
                    # Cannot add this dependency
                    all_deps_ok = False
                    print(f"{Red}Cannot add {var}: Required dependency {dep} is not in the list of supported variables{Nclr}")
            
            if all_deps_ok:
                variables_to_add.append(var)
                return True
            else:
                return False
    
        # Check all requested variables
    for var in add_list:
        if var not in master_list:
            print(f"{Red}Variable '{var}' is not supported and cannot be added to the file.{Nclr}")
            continue
            
        # Skip if already in file
        with Dataset(ifile, "r", format="NETCDF4_CLASSIC") as f:
            if check_variable_exists(var, f.variables.keys()):
                existing_name = get_existing_var_name(var, f.variables.keys())
                if var != existing_name:
                    print(f"{Yellow}Variable '{var}' is already in the file (as '{existing_name}').{Nclr}")
                else:
                    print(f"{Yellow}Variable '{var}' is already in the file.{Nclr}")
                continue
        
        # Try to add the variable and its dependencies
        add_with_dependencies(var)
    
    # Now add the variables in the correct order
    for var in variables_to_add:
        try:
            f = Dataset(ifile, "a", format="NETCDF4_CLASSIC")
            
            # Skip if already in file (double-check)
            if check_variable_exists(var, f.variables.keys()):
                f.close()
                continue
                
            print(f"Processing: {var}...")
            
            # Define lev_T and lev_T_out for this file
            f_type, interp_type = FV3_file_type(f)
            
            if f_type == "diurn":
                lev_T = [2, 1, 0, 3, 4]
                lev_T_out = [1, 2, 0, 3, 4]
                lev_axis = 2
            else:
                lev_T = [1, 0, 2, 3]
                lev_T_out = lev_T
                lev_axis = 1
                
            # Make lev_T and lev_T_out available to compute functions
            globals()['lev_T'] = lev_T
            globals()['lev_T_out'] = lev_T_out
            
            # temp and ps are always required. Get dimension
            dim_out = f.variables["temp"].dimensions
            temp = f.variables["temp"][:]
            shape_out = temp.shape
            
            if interp_type == "pfull":
                # Load ak and bk for pressure calculation.
                # Usually required.
                ak, bk = ak_bk_loader(f)
                
                # level, ps, and p_3d are often required.
                lev = f.variables["pfull"][:]
                ps = f.variables["ps"][:]
                p_3D = compute_p_3D(ps, ak, bk, shape_out)

            elif interp_type == "pstd":
                # If file interpolated to pstd, calculate the 3D
                # pressure field.
                lev = f.variables["pstd"][:]

                # Create the right shape that includes all time steps
                rshp_shape = [1 for i in range(0, len(shape_out))]
                rshp_shape[0] = shape_out[0]  # Set the correct number of time steps
                rshp_shape[lev_axis] = len(lev)

                # Reshape and broadcast properly
                p_levels = lev.reshape([1, len(lev), 1, 1])
                p_3D = np.broadcast_to(p_levels, shape_out)

            else:
                try:
                    # If requested interp_type is zstd, or zagl,
                    # pfull3D is required before interpolation.
                    # Some computations (e.g. wind speed) do not
                    # require pfull3D and will work without it,
                    # so we use a try statement here.
                    p_3D = f.variables["pfull3D"][:]
                except:
                    pass

            if var == "dzTau":
                if "dst_mass_micro" in f.variables.keys():
                    q = f.variables["dst_mass_micro"][:]
                elif "dst_mass_mom" in f.variables.keys():
                    q = f.variables["dst_mass_mom"][:]
                OUT = compute_xzTau(q, temp, lev, C_dst, f_type)

            if var == "izTau":
                if "ice_mass_micro" in f.variables.keys():
                    q = f.variables["ice_mass_micro"][:]
                elif "ice_mass_mom" in f.variables.keys():
                    q = f.variables["ice_mass_mom"][:]
                OUT = compute_xzTau(q, temp, lev, C_ice, f_type)

            if var == "dst_mass_micro" or var == "dst_mass_mom":
                xTau = f.variables["dzTau"][:]
                OUT = compute_mmr(xTau, temp, lev, C_dst, f_type)

            if var == "ice_mass_micro" or var == "ice_mass_mom":
                xTau = f.variables["izTau"][:]
                OUT = compute_mmr(xTau, temp, lev, C_ice, f_type)

            if var == "Vg_sed":
                if "dst_mass_micro" in f.variables.keys():
                    xTau = f.variables["dst_mass_micro"][:]
                elif "dst_mass_mom" in f.variables.keys():
                    xTau = f.variables["dst_mass_mom"][:]
                if "dst_num_micro" in f.variables.keys():
                    nTau = f.variables["dst_num_micro"][:]
                elif "dst_num_mom" in f.variables.keys():
                    nTau = f.variables["dst_num_mom"][:]
                OUT = compute_Vg_sed(xTau, nTau, temp)

            if var == "w_net":
                Vg = f.variables["Vg_sed"][:]
                wvar = f.variables["w"][:]
                OUT = compute_w_net(Vg, wvar)

            if var == "pfull3D":
                OUT = p_3D

            if var == "DP":
                OUT = compute_DP_3D(ps, ak, bk, shape_out)

            if var == "rho":
                OUT = compute_rho(p_3D, temp)

            if var == "theta":
                OUT = compute_theta(p_3D, ps, temp, f_type)

            if var == "w":
                omega = f.variables["omega"][:]
                rho = compute_rho(p_3D, temp)
                OUT = compute_w(rho, omega)

            if var == "zfull":
                OUT = compute_zfull(ps, ak, bk, temp)

            if var == "DZ":
                OUT = compute_DZ_3D(ps, ak, bk, temp, shape_out)

            if var == "wspeed" or var == "wdir":
                ucomp = f.variables["ucomp"][:]
                vcomp = f.variables["vcomp"][:]
                theta, mag = cart_to_azimut_TR(ucomp, vcomp,
                                                mode="from")
                if var == "wdir":
                    OUT = theta
                if var == "wspeed":
                    OUT = mag

            if var == "N":
                theta = compute_theta(p_3D, ps, temp, f_type)
                zfull = compute_zfull(ps, ak, bk, temp)
                OUT = compute_N(theta, zfull)

            if var == "Ri":
                theta = compute_theta(p_3D, ps, temp, f_type)
                zfull = compute_zfull(ps, ak, bk, temp)
                N = compute_N(theta, zfull)

                ucomp = f.variables["ucomp"][:]
                vcomp = f.variables["vcomp"][:]
                du_dz = dvar_dh(
                    ucomp.transpose(lev_T),
                    zfull.transpose(lev_T)
                    ).transpose(lev_T)
                dv_dz = dvar_dh(
                    vcomp.transpose(lev_T),
                    zfull.transpose(lev_T)
                    ).transpose(lev_T)
                OUT = N**2 / (du_dz**2 + dv_dz**2)

            # .. note:: lev_T swaps dims 0 & 1, ensuring level is
            # the first dimension for the differentiation

            if var == "Tco2":
                OUT = compute_Tco2(p_3D)

            if var == "scorer_wl":
                ucomp = f.variables["ucomp"][:]
                theta = compute_theta(p_3D, ps, temp, f_type)
                zfull = compute_zfull(ps, ak, bk, temp)
                N = compute_N(theta, zfull)
                OUT = compute_scorer(N, ucomp, zfull)

            if var in ["div", "curl", "fn"]:
                lat = f.variables["lat"][:]
                lon = f.variables["lon"][:]
                ucomp = f.variables["ucomp"][:]
                vcomp = f.variables["vcomp"][:]

            if var == "div":
                OUT = spherical_div(ucomp, vcomp, lon, lat,
                                    R=3400*1000.,
                                    spacing="regular")

            if var == "curl":
                OUT = spherical_curl(ucomp, vcomp, lon, lat,
                                        R=3400*1000.,
                                        spacing="regular")

            if var == "fn":
                theta = f.variables["theta"][:]
                OUT = frontogenesis(ucomp, vcomp, theta, lon, lat,
                                    R=3400*1000.,
                                    spacing="regular")

            # ==================================================
            #               Interpolated Files
            # ==================================================
            # All interpolated files have the following
            if interp_type != "pfull":
                lev = f.variables[interp_type][:]

            # The next several variables can ONLY be added to
            # pressure interpolated files.
            if var == "msf":
                vcomp = f.variables["vcomp"][:]
                lat = f.variables["lat"][:]
                if f_type == "diurn":
                    # [lev, lat, t, tod, lon]
                    # -> [t, tod, lev, lat, lon]
                    # [0 1 2 3 4] -> [2 3 0 1 4]
                    OUT = mass_stream(
                        vcomp.transpose([2, 3, 0, 1, 4]), lat, lev,
                        type=interp_type
                    ).transpose([2, 3, 0, 1, 4])
                else:
                    OUT = mass_stream(
                        vcomp.transpose([1, 2, 3, 0]), lat, lev,
                        type=interp_type
                    ).transpose([3, 0, 1, 2])
                    # [t, lev, lat, lon]
                    # -> [lev, lat, lon, t]
                    # ->  [t, lev, lat, lon]
                    # [0 1 2 3] -> [1 2 3 0] -> [3 0 1 2]

            if var == "ep":
                OUT = compute_Ep(temp)

            if var == "ek":
                ucomp = f.variables["ucomp"][:]
                vcomp = f.variables["vcomp"][:]
                OUT = compute_Ek(ucomp, vcomp)

            if var == "mx":
                OUT = compute_MF(f.variables["ucomp"][:],
                                    f.variables["w"][:])

            if var == "my":
                OUT = compute_MF(f.variables["vcomp"][:],
                                    f.variables["w"][:])

            if var == "ax":
                mx = compute_MF(f.variables["ucomp"][:],
                                f.variables["w"][:])
                rho = f.variables["rho"][:]
                OUT = compute_WMFF(mx, rho, lev, interp_type)

            if var == "ay":
                my = compute_MF(f.variables["vcomp"][:],
                                f.variables["w"][:])
                rho = f.variables["rho"][:]
                OUT = compute_WMFF(my, rho, lev, interp_type)

            if var == "tp_t":
                OUT = zonal_detrend(temp)/temp

            if interp_type == "pfull":
                # Filter out NANs in the native files
                OUT[np.isnan(OUT)] = fill_value

            else:
                # Add NANs to the interpolated files
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore",
                                            category = RuntimeWarning)
                    OUT[OUT > 1.e30] = np.nan
                    OUT[OUT < -1.e30] = np.nan

            # Log the variable
            var_Ncdf = f.createVariable(var, "f4", dim_out)
            var_Ncdf.long_name = (master_list[var][0]
                                    + cap_str)
            var_Ncdf.units = master_list[var][1]
            var_Ncdf[:] = OUT

            # After successfully adding
            print(f"{Green}*** Variable '{var}' added successfully ***{Nclr}")
            f.close()
            
        except Exception as e:
            except_message(debug, e, var, ifile)
            
# ======================================================
#                  MAIN PROGRAM
# ======================================================

filepath = os.getcwd()

@debug_wrapper
def main():
    # Load all the .nc files
    file_list = [f.name for f in args.input_file]
    add_list = args.add_variable
    zdiff_list = args.differentiate_wrt_z
    zdetrend_list = args.zonal_detrend
    dp_to_dz_list = args.dp_to_dz
    dz_to_dp_list = args.dz_to_dp
    col_list = args.column_integrate
    remove_list = args.remove_variable
    extract_list = args.extract_copy
    edit_var = args.edit_variable
    debug = args.debug

    # An array to swap vertical axis forward and backward:
    # [1, 0, 2, 3] for [t, lev, lat, lon] and
    # [2, 1, 0, 3, 4] for [t, tod, lev, lat, lon]
    global lev_T
    # Reshape ``lev_T_out`` in zfull and zhalf calculation
    global lev_T_out

    # For all the files
    for ifile in file_list:
        # First check if file is on the disk (Lou only)
        # Create a wrapper object to pass to check_file_tape
        class FileWrapper:
            def __init__(self, name):
                self.name = name

        file_wrapper = FileWrapper(ifile)
        check_file_tape(file_wrapper)

        # ==============================================================
        # Remove Function
        # ==============================================================
        if remove_list:
            try:
                # If ncks is available, use it
                cmd_txt = "ncks --version"
                subprocess.check_call(cmd_txt, shell = True,
                                      stdout = open(os.devnull, "w"),
                                      stderr = open(os.devnull, "w"))
                print("Using ncks.")
                for ivar in remove_list:
                    print(f"Creating new file {ifile} without {ivar}:")
                    cmd_txt = f"ncks -C -O -x -v {ivar} {ifile} {ifile}"
                    try:
                        subprocess.check_call(
                            cmd_txt, shell = True,
                            stdout = open(os.devnull, "w"),
                            stderr = open(os.devnull, "w")
                        )
                    except Exception as e:
                        print(f"{e.__class__.__name__ }: "
                              f"{e.message}")

            except subprocess.CalledProcessError:
                # If ncks is not available, use internal method
                print("Using internal method.")
                f_IN = Dataset(ifile, "r", format = "NETCDF4_CLASSIC")
                ifile_tmp = f"{ifile[:-3]}_tmp.nc"
                Log = Ncdf(ifile_tmp, "Edited postprocess")
                Log.copy_all_dims_from_Ncfile(f_IN)
                Log.copy_all_vars_from_Ncfile(f_IN, remove_list)
                f_IN.close()
                Log.close()
                cmd_txt = f"mv {ifile_tmp} {ifile}"
                p = subprocess.run(cmd_txt, universal_newlines = True,
                                   shell = True)
                print(f"{Cyan}{ifile} was updated{Nclr}")

        # ==============================================================
        # Extract Function
        # ==============================================================
        if extract_list:
            f_IN = Dataset(ifile, "r", format = "NETCDF4_CLASSIC")

            # The variable to exclude
            exclude_list = filter_vars(f_IN,
                                       args.extract_copy,
                                       giveExclude = True)
            ifile_tmp = f"{ifile[:-3]}_extract.nc"
            Log = Ncdf(ifile_tmp, "Edited in postprocessing")
            Log.copy_all_dims_from_Ncfile(f_IN)
            Log.copy_all_vars_from_Ncfile(f_IN, exclude_list)
            f_IN.close()
            Log.close()
            print(f"{Cyan}complete{Nclr}\n")

        # ==============================================================
        #  Add Function
        # ==============================================================
        # If the list is not empty, load ak and bk for the pressure
        # calculation. ak and bk are always necessary.
        process_add_variables(ifile, add_list, master_list, debug)

        # First pass: Build the list of variables to add, including dependencies
        # for ivar in add_list:
        #     if ivar not in master_list.keys():
        #         # If the variable to be added is NOT supported
        #         print(f"{Red}Variable ``{ivar}`` is not supported and "
        #               f"cannot be added to the file.{Nclr}")
        #         failed_variables.append(ivar)
        #         continue
            
        #     try:
        #         with Dataset(ifile, "a", format="NETCDF4_CLASSIC") as f:
        #             # Check if variable already exists
        #             if ivar in f.variables:
        #                 print(f"{Yellow}Variable {ivar} already exists in the file.{Nclr}")
        #                 continue
                        
        #             # Check file compatibility
        #             f_type, interp_type = FV3_file_type(f)
        #             compat_file_fmt = ", ".join([cf for cf in master_list[ivar][3]])
                    
        #             if interp_type not in master_list[ivar][3]:
        #                 if "pfull" in compat_file_fmt or "phalf" in compat_file_fmt:
        #                     print(f"{Red}ERROR: Variable '{Yellow}{ivar}{Red}' can only be added to non-interpolated file(s){Nclr}")
        #                 else:
        #                     print(f"{Red}ERROR: Variable '{Yellow}{ivar}{Red}' can only be added to file(s) ending in: {Yellow}{compat_file_fmt}{Nclr}")
        #                 failed_variables.append(ivar)
        #                 continue
                    
        #             # Check dependencies
        #             if check_dependencies(f, ivar, master_list, add_missing=True):
        #                 # If all dependencies are present or can be added, add this variable to our list
        #                 if ivar not in variables_to_add:
        #                     variables_to_add.append(ivar)
                        
        #                 # Check if we need to add any missing dependencies
        #                 for dep_var in master_list[ivar][2]:
        #                     if dep_var not in f.variables and dep_var not in variables_to_add:
        #                         variables_to_add.append(dep_var)
        #             else:
        #                 # If dependencies cannot be satisfied
        #                 failed_variables.append(ivar)
            
        #     except Exception as e:
        #         print(f"{Red}Error checking dependencies for {ivar}: {str(e)}{Nclr}")
        #         failed_variables.append(ivar)
        
        # # Second pass: Add variables in order (dependencies first)
        # for ivar in variables_to_add:
        #     try:
        #         f = Dataset(ifile, "a", format = "NETCDF4_CLASSIC")
        #         f_type, interp_type = FV3_file_type(f)

        #         print(f"Processing: {ivar}...")

                # if interp_type == "pfull":
                #     # Load ak and bk for pressure calculation.
                #     # Usually required.
                #     ak, bk = ak_bk_loader(f)

                # # temp and ps are always required. Get dimension
                # dim_out = f.variables["temp"].dimensions
                # temp = f.variables["temp"][:]
                # shape_out = temp.shape

                # if f_type == "diurn":
                #     # [t, tod, lev, lat, lon]
                #     # -> [lev, tod, t, lat, lon]
                #     # -> [t, tod, lev, lat, lon]
                #     lev_T = [2, 1, 0, 3, 4]
                #     # [0 1 2 3 4] -> [2 1 0 3 4] -> [2 1 0 3 4]
                #     lev_T_out = [1, 2, 0, 3, 4]
                #     # In diurn file, level is the 3rd axis:
                #     # [t, tod, lev, lat, lon]
                #     lev_axis = 2
                # else:
                #     # [t, lev, lat, lon]
                #     # -> [lev, t, lat, lon]
                #     # -> [t, lev, lat, lon]
                #     lev_T = [1, 0, 2, 3]
                #     # [0 1 2 3] -> [1 0 2 3] -> [1 0 2 3]
                #     lev_T_out = lev_T
                #     # In average and daily files, level is the 2nd
                #     # axis = [t, lev, lat, lon]
                #     lev_axis = 1

                # # ==================================================
                # #               Non-Interpolated Files
                # # ==================================================

                # if interp_type == "pfull":
                #     # level, ps, and p_3d are often required.
                #     lev = f.variables["pfull"][:]
                #     ps = f.variables["ps"][:]
                #     p_3D = compute_p_3D(ps, ak, bk, shape_out)

                # elif interp_type == "pstd":
                #     # If file interpolated to pstd, calculate the 3D
                #     # pressure field.
                #     lev = f.variables["pstd"][:]

                #     # Create the right shape that includes all time steps
                #     rshp_shape = [1 for i in range(0, len(shape_out))]
                #     rshp_shape[0] = shape_out[0]  # Set the correct number of time steps
                #     rshp_shape[lev_axis] = len(lev)

                #     # p_3D = lev.reshape(rshp_shape)
                #     # Reshape and broadcast properly
                #     p_levels = lev.reshape([1, len(lev), 1, 1])
                #     p_3D = np.broadcast_to(p_levels, shape_out)

                # else:
                #     try:
                #         # If requested interp_type is zstd, or zagl,
                #         # pfull3D is required before interpolation.
                #         # Some computations (e.g. wind speed) do not
                #         # require pfull3D and will work without it,
                #         # so we use a try statement here.
                #         p_3D = f.variables["pfull3D"][:]
                #     except:
                #         pass

                # if ivar == "dzTau":
                #     if "dst_mass_micro" in f.variables.keys():
                #         q = f.variables["dst_mass_micro"][:]
                #     elif "dst_mass_mom" in f.variables.keys():
                #         q = f.variables["dst_mass_mom"][:]
                #     OUT = compute_xzTau(q, temp, lev, C_dst, f_type)

                # if ivar == "izTau":
                #     if "ice_mass_micro" in f.variables.keys():
                #         q = f.variables["ice_mass_micro"][:]
                #     elif "ice_mass_mom" in f.variables.keys():
                #         q = f.variables["ice_mass_mom"][:]
                #     OUT = compute_xzTau(q, temp, lev, C_ice, f_type)

                # if ivar == "dst_mass_micro" or ivar == "dst_mass_mom":
                #     xTau = f.variables["dzTau"][:]
                #     OUT = compute_mmr(xTau, temp, lev, C_dst, f_type)

                # if ivar == "ice_mass_micro" or ivar == "ice_mass_mom":
                #     xTau = f.variables["izTau"][:]
                #     OUT = compute_mmr(xTau, temp, lev, C_ice, f_type)

                # if ivar == "Vg_sed":
                #     if "dst_mass_micro" in f.variables.keys():
                #         xTau = f.variables["dst_mass_micro"][:]
                #     elif "dst_mass_mom" in f.variables.keys():
                #         xTau = f.variables["dst_mass_mom"][:]
                #     if "dst_num_micro" in f.variables.keys():
                #         nTau = f.variables["dst_num_micro"][:]
                #     elif "dst_num_mom" in f.variables.keys():
                #         nTau = f.variables["dst_num_mom"][:]
                #     OUT = compute_Vg_sed(xTau, nTau, temp)

                # if ivar == "w_net":
                #     Vg = f.variables["Vg_sed"][:]
                #     wvar = f.variables["w"][:]
                #     OUT = compute_w_net(Vg, wvar)

                # if ivar == "pfull3D":
                #     OUT = p_3D

                # if ivar == "DP":
                #     OUT = compute_DP_3D(ps, ak, bk, shape_out)

                # if ivar == "rho":
                #     OUT = compute_rho(p_3D, temp)

                # if ivar == "theta":
                #     OUT = compute_theta(p_3D, ps, temp, f_type)

                # if ivar == "w":
                #     omega = f.variables["omega"][:]
                #     rho = compute_rho(p_3D, temp)
                #     OUT = compute_w(rho, omega)

                # if ivar == "zfull":
                #     OUT = compute_zfull(ps, ak, bk, temp)

                # if ivar == "DZ":
                #     OUT = compute_DZ_3D(ps, ak, bk, temp, shape_out)

                # if ivar == "wspeed" or ivar == "wdir":
                #     ucomp = f.variables["ucomp"][:]
                #     vcomp = f.variables["vcomp"][:]
                #     theta, mag = cart_to_azimut_TR(ucomp, vcomp,
                #                                    mode="from")
                #     if ivar == "wdir":
                #         OUT = theta
                #     if ivar == "wspeed":
                #         OUT = mag

                # if ivar == "N":
                #     theta = compute_theta(p_3D, ps, temp, f_type)
                #     zfull = compute_zfull(ps, ak, bk, temp)
                #     OUT = compute_N(theta, zfull)

                # if ivar == "Ri":
                #     theta = compute_theta(p_3D, ps, temp, f_type)
                #     zfull = compute_zfull(ps, ak, bk, temp)
                #     N = compute_N(theta, zfull)

                #     ucomp = f.variables["ucomp"][:]
                #     vcomp = f.variables["vcomp"][:]
                #     du_dz = dvar_dh(
                #         ucomp.transpose(lev_T),
                #         zfull.transpose(lev_T)
                #         ).transpose(lev_T)
                #     dv_dz = dvar_dh(
                #         vcomp.transpose(lev_T),
                #         zfull.transpose(lev_T)
                #         ).transpose(lev_T)
                #     OUT = N**2 / (du_dz**2 + dv_dz**2)

                # # .. note:: lev_T swaps dims 0 & 1, ensuring level is
                # # the first dimension for the differentiation

                # if ivar == "Tco2":
                #     OUT = compute_Tco2(p_3D)

                # if ivar == "scorer_wl":
                #     ucomp = f.variables["ucomp"][:]
                #     theta = compute_theta(p_3D, ps, temp, f_type)
                #     zfull = compute_zfull(ps, ak, bk, temp)
                #     N = compute_N(theta, zfull)
                #     OUT = compute_scorer(N, ucomp, zfull)

                # if ivar in ["div", "curl", "fn"]:
                #     lat = f.variables["lat"][:]
                #     lon = f.variables["lon"][:]
                #     ucomp = f.variables["ucomp"][:]
                #     vcomp = f.variables["vcomp"][:]

                # if ivar == "div":
                #     OUT = spherical_div(ucomp, vcomp, lon, lat,
                #                         R=3400*1000.,
                #                         spacing="regular")

                # if ivar == "curl":
                #     OUT = spherical_curl(ucomp, vcomp, lon, lat,
                #                          R=3400*1000.,
                #                          spacing="regular")

                # if ivar == "fn":
                #     theta = f.variables["theta"][:]
                #     OUT = frontogenesis(ucomp, vcomp, theta, lon, lat,
                #                         R=3400*1000.,
                #                         spacing="regular")

                # # ==================================================
                # #               Interpolated Files
                # # ==================================================
                # # All interpolated files have the following
                # if interp_type != "pfull":
                #     lev = f.variables[interp_type][:]

                # # The next several variables can ONLY be added to
                # # pressure interpolated files.
                # if ivar == "msf":
                #     vcomp = f.variables["vcomp"][:]
                #     lat = f.variables["lat"][:]
                #     if f_type == "diurn":
                #         # [lev, lat, t, tod, lon]
                #         # -> [t, tod, lev, lat, lon]
                #         # [0 1 2 3 4] -> [2 3 0 1 4]
                #         OUT = mass_stream(
                #             vcomp.transpose([2, 3, 0, 1, 4]), lat, lev,
                #             type=interp_type
                #         ).transpose([2, 3, 0, 1, 4])
                #     else:
                #         OUT = mass_stream(
                #             vcomp.transpose([1, 2, 3, 0]), lat, lev,
                #             type=interp_type
                #         ).transpose([3, 0, 1, 2])
                #         # [t, lev, lat, lon]
                #         # -> [lev, lat, lon, t]
                #         # ->  [t, lev, lat, lon]
                #         # [0 1 2 3] -> [1 2 3 0] -> [3 0 1 2]

                # if ivar == "ep":
                #     OUT = compute_Ep(temp)

                # if ivar == "ek":
                #     ucomp = f.variables["ucomp"][:]
                #     vcomp = f.variables["vcomp"][:]
                #     OUT = compute_Ek(ucomp, vcomp)

                # if ivar == "mx":
                #     OUT = compute_MF(f.variables["ucomp"][:],
                #                      f.variables["w"][:])

                # if ivar == "my":
                #     OUT = compute_MF(f.variables["vcomp"][:],
                #                      f.variables["w"][:])

                # if ivar == "ax":
                #     mx = compute_MF(f.variables["ucomp"][:],
                #                     f.variables["w"][:])
                #     rho = f.variables["rho"][:]
                #     OUT = compute_WMFF(mx, rho, lev, interp_type)

                # if ivar == "ay":
                #     my = compute_MF(f.variables["vcomp"][:],
                #                     f.variables["w"][:])
                #     rho = f.variables["rho"][:]
                #     OUT = compute_WMFF(my, rho, lev, interp_type)

                # if ivar == "tp_t":
                #     OUT = zonal_detrend(temp)/temp

                # if interp_type == "pfull":
                #     # Filter out NANs in the native files
                #     OUT[np.isnan(OUT)] = fill_value

                # else:
                #     # Add NANs to the interpolated files
                #     with warnings.catch_warnings():
                #         warnings.simplefilter("ignore",
                #                               category = RuntimeWarning)
                #         OUT[OUT > 1.e30] = np.nan
                #         OUT[OUT < -1.e30] = np.nan

                # # Log the variable
                # var_Ncdf = f.createVariable(ivar, "f4", dim_out)
                # var_Ncdf.long_name = (master_list[ivar][0]
                #                       + cap_str)
                # var_Ncdf.units = master_list[ivar][1]
                # var_Ncdf[:] = OUT
                
                # print(f"{Green}*** Variable '{ivar}' added successfully ***{Nclr}")
                # f.close()

            # except Exception as e:
            #     except_message(debug, e, ivar, ifile)
                
        # ==============================================================
        #                   Vertical Differentiation
        # ==============================================================
        for idiff in zdiff_list:
            f = Dataset(ifile, "a", format = "NETCDF4_CLASSIC")
            f_type, interp_type = FV3_file_type(f)

            if interp_type == "pfull":
                ak, bk = ak_bk_loader(f)

            if idiff not in f.variables.keys():
                print(f"{Red}zdiff error: variable {idiff} is not "
                      f"present in {ifile}{Nclr}")
                f.close()
            else:
                print(f"Differentiating: {idiff}...")
                if f_type == "diurn":
                    lev_T = [2, 1, 0, 3, 4]
                else:
                    # If [t, lat, lon] -> [lev, t, lat, lon]
                    lev_T = [1, 0, 2, 3]
                try:
                    var = f.variables[idiff][:]
                    lname_text, unit_text = get_longname_unit(f, idiff)
                    # Remove the last ] to update the units (e.g [kg]
                    # to [kg/m])
                    new_unit = f"{unit_text[:-2]}/m]"
                    new_lname = f"vertical gradient of {lname_text}"
                    # temp and ps are always required. Get dimension
                    dim_out = f.variables["temp"].dimensions
                    if interp_type == "pfull":
                        if "zfull" in f.variables.keys():
                            zfull = f.variables["zfull"][:]
                        else:
                            temp = f.variables["temp"][:]
                            ps = f.variables["ps"][:]
                            # Z is the first axis
                            zfull = fms_Z_calc(
                                ps, ak, bk, temp.transpose(lev_T),
                                topo=0., lev_type="full"
                                ).transpose(lev_T)

                        # Average file: zfull = [lev, t, lat, lon]
                        # Diurn file: zfull = [lev, tod, t, lat, lon]
                        # Differentiate the variable w.r.t. Z:
                        darr_dz = dvar_dh(
                            var.transpose(lev_T),
                            zfull.transpose(lev_T)
                            ).transpose(lev_T)

                        # .. note:: lev_T swaps dims 0 & 1, ensuring level
                        # is the first dimension for the differentiation

                    elif interp_type == "pstd":
                        # If pstd, requires zfull
                        if "zfull" in f.variables.keys():
                            zfull = f.variables["zfull"][:]
                            darr_dz = dvar_dh(
                                var.transpose(lev_T),
                                zfull.transpose(lev_T)
                                ).transpose(lev_T)
                        else:
                            lev = f.variables[interp_type][:]
                            temp = f.variables["temp"][:]
                            dzfull_pstd = compute_DZ_full_pstd(
                                lev, temp
                                )
                            darr_dz = (dvar_dh(
                                var.transpose(lev_T)
                                ).transpose(lev_T) / dzfull_pstd)

                    elif interp_type in ["zagl", "zstd"]:
                        lev = f.variables[interp_type][:]
                        darr_dz = dvar_dh(
                            var.transpose(lev_T), lev
                            ).transpose(lev_T)
                    # .. note:: lev_T swaps dims 0 & 1, ensuring level is
                    # the first dimension for the differentiation

                    # Log the variable
                    var_Ncdf = f.createVariable(
                        f"d_dz_{idiff}", "f4", dim_out
                        )
                    var_Ncdf.long_name = (new_lname + cap_str)
                    var_Ncdf.units = new_unit
                    var_Ncdf[:] = darr_dz
                    f.close()

                    print(f"d_dz_{idiff}: {Green}Done{Nclr}")
                except Exception as e:
                    except_message(debug,e,idiff,ifile,pre="d_dz_")

        # ==============================================================
        #                       Zonal Detrending
        # ==============================================================
        for izdetrend in zdetrend_list:
            f = Dataset(ifile, "a", format = "NETCDF4_CLASSIC")
            f_type, interp_type = FV3_file_type(f)
            if izdetrend not in f.variables.keys():
                print(f"{Red}zdiff error: variable {izdetrend} is not "
                      f"in {ifile}{Nclr}")
                f.close()
            else:
                print(f"Detrending: {izdetrend}...")
                try:
                    var = f.variables[izdetrend][:]
                    lname_text, unit_text = get_longname_unit(
                        f, izdetrend
                        )
                    new_lname = f"zonal perturbation of {lname_text}"

                    # Get dimension
                    dim_out = f.variables[izdetrend].dimensions

                    # Log the variable
                    var_Ncdf = f.createVariable(
                        izdetrend+"_p", "f4", dim_out
                        )
                    var_Ncdf.long_name = (new_lname + cap_str)
                    var_Ncdf.units = unit_text
                    var_Ncdf[:] = zonal_detrend(var)
                    f.close()

                    print(f"{izdetrend}_p: {Green}Done{Nclr}")
                except Exception as e:
                    except_message(debug,e,izdetrend,ifile,ext="_p")

        # ==============================================================
        #           Opacity Conversion (dp_to_dz and dz_to_dp)
        # ==============================================================
        for idp_to_dz in dp_to_dz_list:
            # ========= Case 1: dp_to_dz
            f = Dataset(ifile, "a", format = "NETCDF4_CLASSIC")
            f_type, interp_type = FV3_file_type(f)
            if idp_to_dz not in f.variables.keys():
                print(f"{Red}dp_to_dz error: variable {idp_to_dz} is "
                      f"not in {ifile}{Nclr}")
                f.close()
            else:
                print("Converting: {idp_to_dz}...")

                try:
                    var = f.variables[idp_to_dz][:]
                    new_unit = (getattr(f.variables[idp_to_dz],
                                        "units", "")
                                + "/m")
                    new_lname = (getattr(f.variables[idp_to_dz],
                                         "long_name", "")
                                 + " rescaled to meter-1")
                    # Get dimension
                    dim_out = f.variables[idp_to_dz].dimensions

                    # Log the variable
                    var_Ncdf = f.createVariable(
                        f"{idp_to_dz}_dp_to_dz", "f4", dim_out
                        )
                    var_Ncdf.long_name = (new_lname + cap_str)
                    var_Ncdf.units = new_unit
                    var_Ncdf[:] = (var
                                   * f.variables["DP"][:]
                                   / f.variables["DZ"][:])
                    f.close()

                    print(f"{idp_to_dz}_dp_to_dz: {Green}Done{Nclr}")

                except Exception as e:
                    except_message(debug,e,idp_to_dz,ifile,ext="_dp_to_dz")

        for idz_to_dp in dz_to_dp_list:
            # ========= Case 2: dz_to_dp
            f = Dataset(ifile, "a", format = "NETCDF4_CLASSIC")
            f_type, interp_type = FV3_file_type(f)
            if idz_to_dp not in f.variables.keys():
                print(f"{Red}dz_to_dp error: variable {idz_to_dp} is "
                      f"not in {ifile}{Nclr}")
                f.close()
            else:
                print(f"Converting: {idz_to_dp}...")

                try:
                    var = f.variables[idz_to_dp][:]
                    new_unit = (getattr(f.variables[idz_to_dp],
                                        "units", "")
                                + "/m")
                    new_lname = (getattr(f.variables[idz_to_dp],
                                         "long_name", "")
                                 + " rescaled to Pa-1")
                    # Get dimension
                    dim_out = f.variables[idz_to_dp].dimensions

                    # Log the variable
                    var_Ncdf = f.createVariable(
                        f"{idz_to_dp}_dz_to_dp", "f4", dim_out
                        )
                    var_Ncdf.long_name = (new_lname + cap_str)
                    var_Ncdf.units = new_unit
                    var_Ncdf[:] = (var * f.variables["DP"][:] 
                                   / f.variables["DZ"][:])
                    f.close()

                    print(f"{idz_to_dp}_dz_to_dp: {Green}Done{Nclr}")

                except Exception as e:
                    except_message(debug,e,idz_to_dp,ifile,ext="_dz_to_dp")

        # ==============================================================
        #                    Column Integration
        # ==============================================================
        """
        Column-integrate the variable::

                            z_top
                            ⌠
            We have col=    ⌡ var (rho dz)
                            0

            with [(dp/dz) = (-rho g)] => [(rho dz) = (-dp/g)]

                        ___ p_sfc
                >  col = \
                        /__ var (dp/g)
                            p_top

        """
        for icol in col_list:
            f = Dataset(ifile, "a")
            f_type, interp_type = FV3_file_type(f)

            if interp_type == "pfull":
                ak, bk = ak_bk_loader(f)
            if icol not in f.variables.keys():
                print(f"{Red}column integration error: variable {icol} "
                      f"is not in {ifile}{Nclr}")
                f.close()
            else:
                print(f"Performing column integration: {icol}...")
                try:
                    var = f.variables[icol][:]
                    lname_text, unit_text = get_longname_unit(f, icol)
                    # turn "kg/kg" -> "kg/m2"
                    new_unit = f"{unit_text[:-3]}/m2"
                    new_lname = f"column integration of {lname_text}"
                    # temp and ps always required
                    # Get dimension
                    dim_in = f.variables["temp"].dimensions
                    shape_in = f.variables["temp"].shape
                    # TODO edge cases where time = 1
                    if f_type == "diurn":
                        # if [t, tod, lat, lon]
                        lev_T = [2, 1, 0, 3, 4]
                        # -> [lev, tod, t, lat, lon]
                        dim_out = tuple(
                            [dim_in[0], dim_in[1], dim_in[3], dim_in[4]]
                            )
                        # In diurn, lev is the 3rd axis (index 2):
                        # [t, tod, lev, lat, lon]
                        lev_axis = 2
                    else:
                        # if [t, lat, lon]
                        lev_T = [1, 0, 2, 3]
                        # -> [lev, t, lat, lon]
                        dim_out = tuple(
                            [dim_in[0], dim_in[2], dim_in[3]]
                            )
                        lev_axis = 1

                    ps = f.variables["ps"][:]
                    DP = compute_DP_3D(ps, ak, bk, shape_in)
                    out = np.sum(var * DP / g, axis=lev_axis)

                    # Log the variable
                    var_Ncdf = f.createVariable(
                        f"{icol}_col", "f4", dim_out
                        )
                    var_Ncdf.long_name = (new_lname + cap_str)
                    var_Ncdf.units = new_unit
                    var_Ncdf[:] = out
                    f.close()

                    print(f"{icol}_col: {Green}Done{Nclr}")

                except Exception as e:
                    except_message(debug,e,icol,ifile,ext="_col")
        if edit_var:
            f_IN = Dataset(ifile, "r", format = "NETCDF4_CLASSIC")
            ifile_tmp = f"{ifile[:-3]}_tmp.nc"
            Log = Ncdf(ifile_tmp, "Edited in postprocessing")
            Log.copy_all_dims_from_Ncfile(f_IN)

            # Copy all variables but this one
            Log.copy_all_vars_from_Ncfile(f_IN, exclude_var = edit_var)

            # Read value, longname, units, name, and log the new var
            var_Ncdf = f_IN.variables[edit_var]

            name_text = edit_var
            vals = var_Ncdf[:]
            dim_out = var_Ncdf.dimensions
            lname_text = getattr(var_Ncdf, "long_name", "")
            unit_text = getattr(var_Ncdf, "units", "")
            cart_text = getattr(var_Ncdf, "cartesian_axis", "")

            if args.rename:
                name_text = args.rename
            if args.longname:
                lname_text = args.longname
            if args.unit:
                unit_text = args.unit
            if args.multiply:
                vals *= args.multiply

            if cart_text == "":
                Log.log_variable(
                    name_text, vals, dim_out, lname_text, unit_text
                    )
            else:
                Log.log_axis1D(name_text, vals, dim_out, lname_text,
                               unit_text, cart_text)
            f_IN.close()
            Log.close()

            # Rename the new file
            cmd_txt = f"mv {ifile_tmp} {ifile}"
            subprocess.call(cmd_txt, shell = True)

            print(f"{Cyan}{ifile} was updated{Nclr}")

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
