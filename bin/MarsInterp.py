#!/usr/bin/env python3
"""
The MarsInterp executable is for interpolating files to pressure or
altitude coordinates. Options include interpolation to standard
pressure (``pstd``), standard altitude (``zstd``), altitude above
ground level (``zagl``), or a custom vertical grid.

The executable requires:

    * ``[input_file]``          The file to be transformed

and optionally accepts:

    * ``[-t --interp_type]``    Type of interpolation to perform (altitude, pressure, etc.)
    * ``[-v --vertical_grid]``  Specific vertical grid to interpolate to
    * ``[-incl --include]``     Variables to include in the new interpolated file
    * ``[-ext --extension]``    Custom extension for the new file
    * ``[-print --print_grid]`` Print the vertical grid to the screen


Third-party Requirements:

    * ``numpy``
    * ``netCDF4``
    * ``argparse``
    * ``os``
    * ``time``
    * ``matplotlib``
    * ``re``
    * ``functools``
    * ``traceback``
    * ``sys``
    * ``amescap``
"""

# Make print statements appear in color
from amescap.Script_utils import (
    Cyan, Red, Blue, Yellow, Nclr, Green, Cyan
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import time         # Monitor interpolation time
import re           # Regular expressions
import matplotlib
import numpy as np
from netCDF4 import Dataset
import functools    # For function decorators
import traceback    # For printing stack traces

# Force matplotlib NOT to load Xwindows backend
matplotlib.use("Agg")

# Load amesCAP modules
from amescap.FV3_utils import (
    fms_press_calc, fms_Z_calc, vinterp,find_n
)
from amescap.Script_utils import (
    check_file_tape, section_content_amescap_profile, find_tod_in_diurn,
    filter_vars, find_fixedfile, ak_bk_loader,
    read_variable_dict_amescap_profile
)
from amescap.Ncdf_wrapper import Ncdf


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


# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    prog=('MarsInterp'),
    description=(
        f"{Yellow}Performs a pressure interpolation on the vertical "
        f"coordinate of the netCDF file.{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('input_file', nargs='+',
    type=argparse.FileType('rb'),
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument('-t', '--interp_type', type=str, default='pstd',
    help=(
        f"Interpolation to standard pressure (pstd), standard altitude "
        f"(zstd), or altitude above ground level (zagl).\nWorks on "
        f"'daily', 'average', and 'diurn' files.\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_average.nc\n"
        f"> MarsInterp 01336.atmos_average.nc -t pstd\n"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-v', '--vertical_grid', type=str, default=None,
    help=(
        f"For use with ``-t``. Specify a custom vertical grid to "
        f"interpolate to.\n"
        f"Custom grids defined in ``amescap_profile``.\nFor first "
        f"time use, copy ``amescap_profile`` to your home directory:\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Cyan}cp path/to/amesCAP/mars_templates/amescap_profile "
        f"~/.amescap_profile\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_average.nc -t zstd -v phalf_mb"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-incl', '--include', nargs='+',
    help=(
        f"Only include the listed variables in the action. Dimensions "
        f"and 1D variables are always included.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_daily.nc -incl temp ps ts"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-print', '--print_grid', action='store_true',
    help=(
        f"Print the vertical grid to the screen.\n{Yellow}This does not "
        f"run the interpolation, it only prints grid information.\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_average.nc -t pstd -v pstd_default -print"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-ext', '--extension', type=str, default=None,
    help=(
        f"Must be paired with an argument listed above.\nInstead of "
        f"overwriting a file to perform a function, ``-ext``\ntells "
        f"CAP to create a new file with the extension name specified "
        f"here.\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_average.nc -t pstd -ext _my_pstd\n"
        f"{Blue}(Produces 01336.atmos_average_my_pstd.nc and "
        f"preserves all other files)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsInterp 01336.atmos_average.nc -t pstd --debug"
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


# ======================================================================
#                           DEFINITIONS
# ======================================================================

# TODO: If only one time step, reshape from (lev,lat,lon) to
# (time, lev, lat, lon).

# Fill values for NaN. Do not use np.NaN because it is deprecated and
# will raise issues when using runpinterp
fill_value = 0.

# Define constants
rgas = 189.     # J/(kg-K) -> m2/(s2 K)
g = 3.72        # m/s2
R = 8.314       # J/ mol. K
Cp = 735.0      # J/K
M_co2 = 0.044   # kg/mol

filepath = os.getcwd()

# ======================================================================
#                               MAIN PROGRAM
# ======================================================================


@debug_wrapper
def main():
    start_time   = time.time()
    debug        = args.debug
    # Load all of the netcdf files
    file_list    = file_list = [f.name for f in args.input_file]
    interp_type  = args.interp_type  # e.g. pstd
    custom_level = args.vertical_grid # e.g. p44
    grid_out     = args.print_grid

    # Create a namespace with numpy available
    namespace = {'np': np}

    # PRELIMINARY DEFINITIONS
    # =========================== pstd ===========================
    if interp_type == "pstd":
        longname_txt = "standard pressure"
        units_txt = "Pa"
        need_to_reverse = False
        interp_technic = "log"

        content_txt = section_content_amescap_profile("Pressure definitions for pstd")

        # Execute in controlled namespace
        exec(content_txt, namespace)

        if custom_level:
            lev_in = eval(f"np.array({custom_level})", namespace)
        else:
            lev_in = np.array(namespace['pstd_default'])

    # =========================== zstd ===========================
    elif interp_type == "zstd":
        longname_txt = "standard altitude"
        units_txt = "m"
        need_to_reverse = True
        interp_technic = "lin"

        content_txt = section_content_amescap_profile("Altitude definitions "
                                                      "for zstd")
        # Load all variables in that section
        exec(content_txt, namespace)

        if custom_level:
            lev_in = eval(f"np.array({custom_level})", namespace)
        else:
            lev_in = eval("np.array(zstd_default)", namespace)

        # The fixed file is necessary if pk, bk are not in the
        # requested file, orto load the topography if zstd output is
        # requested.
        name_fixed = find_fixedfile(file_list[0])
        try:
            f_fixed = Dataset(name_fixed, 'r')
            model=read_variable_dict_amescap_profile(f_fixed)
            zsurf = f_fixed.variables["zsurf"][:]
            f_fixed.close()
        except FileNotFoundError:
            print(f"{Red}***Error*** Topography (zsurf) is required for "
                  f"interpolation to zstd, but the file {name_fixed} "
                  f"cannot be not found{Nclr}")
            exit()

    # =========================== zagl ===========================
    elif interp_type == "zagl":
        longname_txt = "altitude above ground level"
        units_txt = "m"
        need_to_reverse = True
        interp_technic = "lin"

        content_txt = section_content_amescap_profile("Altitude definitions "
                                                      "for zagl")
        # Load all variables in that section
        exec(content_txt, namespace)

        if custom_level:
            lev_in = eval(f"np.array({custom_level})", namespace)
        else:
            lev_in = eval("np.array(zagl_default)", namespace)
    else:
        print(f"{Red}Interpolation interp_ {interp_type} is not supported, use "
              f"``pstd``, ``zstd`` or ``zagl``{Nclr}")
        exit()

    if grid_out:
        # Only print grid content and exit the code
        print(*lev_in)
        exit()

    for ifile in file_list:
        # First check if file is present on the disk (Lou only)
        check_file_tape(ifile)

        # Append extension, if any
        if args.extension:
            newname = (f"{filepath}/{ifile[:-3]}_{interp_type}_"
                       f"{args.extension}.nc")
        else:
            newname = (f"{filepath}/{ifile[:-3]}_{interp_type}.nc")


        # ==============================================================
        #                       Interpolation
        # ==============================================================

        fNcdf = Dataset(ifile, "r", format = "NETCDF4_CLASSIC")
        # Load pk, bk, and ps for 3D pressure field calculation.
        # Read the pk and bk for each file in case the vertical resolution has changed.
        model=read_variable_dict_amescap_profile(fNcdf)
        ak, bk = ak_bk_loader(fNcdf)
        ps = np.array(fNcdf.variables["ps"])

        ps = np.array(fNcdf.variables["ps"])

        if len(ps.shape) == 3:
            do_diurn = False
            tod_name = "not_used"
            # Put vertical axis first for 4D variable,
            # e.g., [time, lev, lat, lon] -> [lev, time, lat, lon]
            permut = [1, 0, 2, 3]
            # [0 1 2 3] -> [1 0 2 3]
        elif len(ps.shape) == 4:
            do_diurn = True
            # Find time_of_day variable name
            tod_name = find_tod_in_diurn(fNcdf)
            # Same for diurn files,
            # e.g., [time, time_of_day_XX, lev, lat, lon]
            # -> [lev, time_of_day_XX, time, lat, lon]
            permut = [2, 1, 0, 3, 4]
            # [0 1 2 3 4] -> [2 1 0 3 4]

        # Compute levels in the file, these are permutted arrays
        # Suppress "divide by zero" error
        with np.errstate(divide = "ignore", invalid = "ignore"):
            if interp_type == "pstd":
                # Permute by default dimension, e.g., lev is first
                L_3D_P = fms_press_calc(ps, ak, bk, lev_type = "full")

            elif interp_type == 'zagl':
                temp = fNcdf.variables["temp"][:]
                L_3D_P = fms_Z_calc(ps, ak, bk, temp.transpose(
                    permut), topo=0., lev_type='full')

            elif interp_type == 'zstd':
                temp = fNcdf.variables["temp"][:]
                # Expand the 'zsurf' array to the 'time' dimension
                zflat = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis=0)
                if do_diurn:
                    zflat = np.repeat(zflat[:, np.newaxis, :, :], ps.shape[1],
                                      axis = 1)

                L_3D_P = fms_Z_calc(ps, ak, bk, temp.transpose(permut),
                                    topo = zflat, lev_type = "full")

        fnew = Ncdf(newname, "Pressure interpolation using MarsInterp")

        # Copy existing DIMENSIONS other than pfull
        # Get all variables in the file
        # var_list=fNcdf.variables.keys()
        # Get the variables
        var_list = filter_vars(fNcdf, args.include)

        fnew.copy_all_dims_from_Ncfile(fNcdf, exclude_dim=["pfull"])
        # Add new vertical dimension
        fnew.add_dim_with_content(interp_type, lev_in, longname_txt, units_txt)

        #TODO :this is fine but FV3-specific, is there a more flexible approach?
        if "tile" in ifile:
            fnew.copy_Ncaxis_with_content(fNcdf.variables["grid_xt"])
            fnew.copy_Ncaxis_with_content(fNcdf.variables["grid_yt"])
        else:
            fnew.copy_Ncaxis_with_content(fNcdf.variables["lon"])
            fnew.copy_Ncaxis_with_content(fNcdf.variables["lat"])

        fnew.copy_Ncaxis_with_content(fNcdf.variables["time"])

        if do_diurn:
            fnew.copy_Ncaxis_with_content(fNcdf.variables[tod_name])

        # Re-use the indices for each file, speeds up the calculation
        compute_indices = True
        for ivar in var_list:
            if (fNcdf.variables[ivar].dimensions == ("time", "pfull", "lat",
                                                     "lon") or
                fNcdf.variables[ivar].dimensions == ("time", tod_name, "pfull",
                                                     "lat", "lon") or
                fNcdf.variables[ivar].dimensions == ("time", "pfull",
                                                     "grid_yt", "grid_xt")):
                if compute_indices:
                    print(f"{Cyan}Computing indices ...{Nclr}")
                    index = find_n(L_3D_P, lev_in,
                                   reverse_input = need_to_reverse)
                    compute_indices = False

                print(f"{Cyan}Interpolating: {ivar} ...{Nclr}")
                varIN = fNcdf.variables[ivar][:]
                # This with the loop suppresses "divide by zero" errors
                with np.errstate(divide = "ignore", invalid = "ignore"):
                    varOUT = vinterp(varIN.transpose(permut), L_3D_P, lev_in,
                                     type_int = interp_technic,
                                     reverse_input = need_to_reverse,
                                     masktop = True,
                                     index = index).transpose(permut)

                long_name_txt = getattr(fNcdf.variables[ivar], "long_name", "")
                units_txt = getattr(fNcdf.variables[ivar], "units", "")
                # long_name_txt=fNcdf.variables[ivar].long_name
                # units_txt=fNcdf.variables[ivar].units)

                if not do_diurn:
                    if "tile" in ifile:
                        fnew.log_variable(ivar, varOUT, ("time", interp_type,
                                                         "grid_yt", "grid_xt"),
                                          long_name_txt, units_txt)
                    else:
                        fnew.log_variable(ivar, varOUT, ("time", interp_type,
                                                         "lat", "lon"),
                                          long_name_txt, units_txt)
                else:
                    if "tile" in ifile:
                        fnew.log_variable(ivar, varOUT, ("time", tod_name,
                                                         interp_type,
                                                         "grid_yt", "grid_xt"),
                                          long_name_txt, units_txt)
                    else:
                        fnew.log_variable(ivar, varOUT, ("time", tod_name,
                                                         interp_type, "lat",
                                                         "lon"),
                                          long_name_txt, units_txt)
            else:

                #TODO logic could be improved over here
                if ivar not in ["time", "pfull", "lat",
                                "lon", 'phalf', 'ak', 'pk', 'bk',
                                "pstd", "zstd", "zagl",
                                tod_name, 'grid_xt', 'grid_yt']:

                    dim_list=fNcdf.dimensions.keys()

                    if 'pfull' not in fNcdf.variables[ivar].dimensions:
                        print(f"{Cyan}Copying over: {ivar}...")
                        if ivar in dim_list:
                            fnew.copy_Ncaxis_with_content(fNcdf.variables[ivar])
                        else:
                            fnew.copy_Ncvar(fNcdf.variables[ivar])

        print("\r ", end="")
        fNcdf.close()
        fnew.close()
        print(f"Completed in {(time.time() - start_time):3f} sec")


# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
