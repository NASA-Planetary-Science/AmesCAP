#!/usr/bin/env python3
"""
The MarsInterp executable is for interpolating files to pressure or \
altitude coordinates. Options include interpolation to standard \
pressure (``pstd``), standard altitude (``zstd``), altitude above \
ground level (``zagl``), or a custom vertical grid.

The executable requires:
    * ``[input_file]``          the file to be transformed

and optionally accepts:
    * ``[-t --type]``           type of interpolation to perform \
        (altitude, pressure, etc.)
    * ``[-l --level]``          specific vertical grid to interpolate to
    * ``[-include --include]``  variables to include in the new \
        interpolated file
    * ``[-e --ext]``            custom extension for the new file
    * ``[-g --grid]``           print the vertical grid chosen by \
        [-l --level] to the screen


Third-party Requirements:
    * ``numpy``
    * ``netCDF4``
    * ``argparse``
    * ``os``
    * ``time``
    * ``matplotlib``
"""

# make print statements appear in color
from amescap.Script_utils import (
    prCyan, prRed, Blue, Yellow,  NoColor, Green, Cyan
)

# load generic Python modules
import argparse     # Parse arguments
import os           # Access operating system functions
import time         # Monitor interpolation time
import matplotlib
import numpy as np
from netCDF4 import Dataset

# Force matplotlib NOT load Xwindows backend
matplotlib.use("Agg") 

# load amesCAP modules
from amescap.FV3_utils import (
    fms_press_calc, fms_Z_calc, vinterp,find_n
)
from amescap.Script_utils import (
    check_file_tape, section_content_amescap_profile, find_tod_in_diurn, 
    filter_vars, find_fixedfile, ak_bk_loader
)
from amescap.Ncdf_wrapper import Ncdf

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}MarsInterp, pressure interpolation on fixed "
        f"layers.{NoColor}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("input_file", nargs="+",
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument("-t", "--type", type=str, default="pstd",
    help=(
        f"Interpolation type. Accepts ``pstd``, ``zstd``, or "
        f"``zagl``.\n{Green}Usage:\n"
        f"> MarsInterp.py ****.atmos.average.nc\n"
        f"> MarsInterp.py ****.atmos.average.nc -t zstd\n"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-l", "--level", type=str, default=None,
    help=(
        f"Layer IDs as defined in ``~/.amescap_profile``. For first "
        f"time use, copy ``~/.amescap_profile`` to ``~/amesCAP``:\n"
        f"{Cyan}cp ~/amesCAP/mars_templates/amescap_profile "
        f"~/.amescap_profile\n"
        f"{Green}Usage:\n"
        f"> MarsInterp.py ****.atmos.average.nc -t pstd -l p44\n"
        f"> MarsInterp.py ****.atmos.average.nc -t zstd -l phalf_mb\n"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-include", "--include", nargs="+",
    help=(
        f"Only include the listed variables. Dimensions and 1D "
        f"variables are always included.\n"
        f"{Green}Usage:\n"
        f"> MarsInterp.py *.atmos_daily.nc --include ps ts temp\n"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-e", "--ext", type=str, default=None,
    help=(
        f"Append an extension (``_ext.nc``) to the output file instead"
        f" of replacing the existing file.\n"
        f"{Green}Usage:\n"
        f"> MarsInterp.py ****.atmos.average.nc -ext B\n"
        f"  {Blue}Produces ****.atmos.average_pstd_B.nc.\n"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-g", "--grid", action="store_true",
    help=(
        f"Output current grid information to standard output. This "
        f"will not run the interpolation.\n"
        f"{Green}Usage:\n"
        f"> MarsInterp.py ****.atmos.average.nc -t pstd -l p44 -g\n"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("--debug", action="store_true",
    help=(f"Debug flag: do not bypass errors.\n\n"))

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

def main():
    start_time   = time.time()
    debug        = parser.parse_args().debug
    # Load all of the netcdf files
    file_list    = parser.parse_args().input_file
    interp_type  = parser.parse_args().type  # e.g. pstd
    custom_level = parser.parse_args().level # e.g. p44
    grid_out     = parser.parse_args().grid

    # PRELIMINARY DEFINITIONS
    # =========================== pstd ===========================
    if interp_type == "pstd":
        longname_txt = "standard pressure"
        units_txt = "Pa"
        need_to_reverse = False
        interp_technic = "log"

        content_txt = section_content_amescap_profile(
            "Pressure definitions for pstd"
        )
        # Load all variables in that section
        exec(content_txt)  

        if custom_level:
            lev_in = eval(f"np.array({custom_level})")
        else:
            lev_in = eval("np.array(pstd_default)")

    # =========================== zstd ===========================
    elif interp_type == "zstd":
        longname_txt = "standard altitude"
        units_txt = "m"
        need_to_reverse = True
        interp_technic = "lin"

        content_txt = section_content_amescap_profile(
            "Altitude definitions for zstd"
        )
        # Load all variables in that section
        exec(content_txt)  

        if custom_level:
            lev_in = eval(f"np.array({custom_level})")
        else:
            lev_in = eval("np.array(zstd_default)")
            # Default levels, this is size 45

        # The fixed file is necessary if pk, bk are not in the 
        # requested file, orto load the topography if zstd output is 
        # requested.
        name_fixed = find_fixedfile(file_list[0])
        try:
            f_fixed = Dataset(name_fixed, "r")
            zsurf = f_fixed.variables["zsurf"][:]
            f_fixed.close()
        except FileNotFoundError:
            prRed("***Error*** Topography (zsurf) is required for \
                interpolation to zstd, but the file {name_fixed} \
                    cannot be not found")
            exit()

    # =========================== zagl ===========================
    elif interp_type == "zagl":
        longname_txt = "altitude above ground level"
        units_txt = "m"
        need_to_reverse = True
        interp_technic = "lin"

        content_txt = section_content_amescap_profile(
            "Altitude definitions for zagl"
        )
        # Load all variables in that section
        exec(content_txt)  

        if custom_level:
            lev_in = eval(f"np.array({custom_level})")
        else:
            lev_in = eval("np.array(zagl_default)")
    else:
        prRed(f"Interpolation type {interp_type} is not supported, \
            use ``pstd``, ``zstd`` or ``zagl``" )
        exit()
    
    if grid_out:
        # Only print grid content and exit the code
        print(*lev_in)
        exit()
    
    for ifile in file_list:
        # First check if file is present on the disk (Lou only)
        check_file_tape(ifile)

        # Append extension, if any
        if parser.parse_args().ext:
            newname = (f"{filepath}/{ifile[:-3]}_{interp_type}_\
                {parser.parse_args().ext}.nc")
        else:
            newname = (f"{filepath}/{ifile[:-3]}_{interp_type}.nc")

        # ==============================================================
        # ======================== Interpolation =======================
        # ==============================================================

        fNcdf = Dataset(ifile, "r", format="NETCDF4_CLASSIC")
        # Load pk, bk, and ps for 3D pressure field calculation.
        # Read the pk and bk for each file in case the vertical 
        # resolution has changed.
        ak, bk = ak_bk_loader(fNcdf)

        ps = np.array(fNcdf.variables["ps"])

        # For pstd only: Uncommenting the following line will use 
        # pfull as default layers:
        # if interp_type == "pstd":lev_in=fNcdf.variables["pfull"][::-1]
        
        if len(ps.shape) == 3:
            do_diurn = False
            tod_name = "not_used"
            # Put vertical axis first for 4D variable, 
            # e.g (time, lev, lat, lon) >>> (lev, time, lat, lon)
            permut = [1, 0, 2, 3]
            # ( 0 1 2 3 ) >>> ( 1 0 2 3 )
        elif len(ps.shape) == 4:
            do_diurn = True
            # Find time_of_day variable name
            tod_name = find_tod_in_diurn(fNcdf)
            # Same for diurn files, 
            # e.g (time, time_of_day_XX, lev, lat, lon) 
            # >>> (lev, time_of_day_XX, time, lat, lon)
            permut = [2, 1, 0, 3, 4]
            # ( 0 1 2 3 4) >>> ( 2 1 0 3 4 )

        # Compute levels in the file, these are permutted arrays
        # Suppress "divide by zero" error
        with np.errstate(divide="ignore", invalid="ignore"):
            if interp_type == "pstd":
                # Permute by default dimension, e.g lev is first
                L_3D_P = fms_press_calc(ps, ak, bk, lev_type="full")

            elif interp_type == "zagl":
                temp = fNcdf.variables["temp"][:]
                L_3D_P = fms_Z_calc(ps, ak, bk, 
                                    temp.transpose(permut), 
                                    topo=0.,
                                    lev_type="full")

            elif interp_type == "zstd":
                temp = fNcdf.variables["temp"][:]
                # Expand the zsurf array to the time dimension
                zflat = np.repeat(zsurf[np.newaxis, :], ps.shape[0], axis=0)
                if do_diurn:
                    zflat = np.repeat(zflat[:, np.newaxis, :, :],
                                      ps.shape[1],
                                      axis=1)

                L_3D_P = fms_Z_calc(ps, ak, bk, 
                                    temp.transpose(permut), 
                                    topo=zflat, 
                                    lev_type="full")

        fnew = Ncdf(newname, "Pressure interpolation using MarsInterp.py")

        # Copy existing DIMENSIONS other than pfull
        # Get all variables in the file
        # var_list=fNcdf.variables.keys()
        # Get the variables
        var_list = filter_vars(fNcdf, parser.parse_args().include)  

        fnew.copy_all_dims_from_Ncfile(fNcdf, exclude_dim=["pfull"])
        # Add new vertical dimension
        fnew.add_dim_with_content(interp_type, lev_in, longname_txt, units_txt)

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
            if (fNcdf.variables[ivar].dimensions == ("time", "pfull", 
                                                     "lat", "lon") or
                fNcdf.variables[ivar].dimensions == ("time", tod_name, "pfull", 
                                                     "lat", "lon") or
                fNcdf.variables[ivar].dimensions == ("time", "pfull", 
                                                     "grid_yt", "grid_xt")):
                if compute_indices:
                    prCyan("Computing indices ...")
                    index = find_n(L_3D_P, lev_in, 
                                   reverse_input=need_to_reverse)
                    compute_indices = False

                prCyan(f"Interpolating: {ivar} ...")
                varIN = fNcdf.variables[ivar][:]
                # This with the loop suppresses "divide by zero" errors
                with np.errstate(divide="ignore", invalid="ignore"):
                    varOUT = vinterp(
                        varIN.transpose(permut), L_3D_P, lev_in, 
                        type_int=interp_technic, reverse_input=need_to_reverse, 
                        masktop=True, index=index).transpose(permut)

                long_name_txt = getattr(fNcdf.variables[ivar], "long_name", "")
                units_txt = getattr(fNcdf.variables[ivar], "units", "")
                # long_name_txt=fNcdf.variables[ivar].long_name
                # units_txt=fNcdf.variables[ivar].units)

                if not do_diurn:
                    if "tile" in ifile:
                        fnew.log_variable(
                            ivar, varOUT, ("time", interp_type, 
                                           "grid_yt", "grid_xt"),
                            long_name_txt, units_txt)
                    else:
                        fnew.log_variable(
                            ivar, varOUT, ("time", interp_type, "lat", "lon"),
                            long_name_txt, units_txt)
                else:
                    if "tile" in ifile:
                        fnew.log_variable(
                            ivar, varOUT, ("time", tod_name, interp_type, 
                                           "grid_yt", "grid_xt"),
                            long_name_txt, units_txt)
                    else:
                        fnew.log_variable(
                            ivar, varOUT, ('time', tod_name, interp_type,
                                           'lat', 'lon'),
                            long_name_txt, units_txt)
            else:
                if ivar not in ['time', 'pfull', 'lat', 'lon', 'phalf',
                                'ak', 'pk', 'bk', 'pstd', 'zstd', 
                                'zagl', tod_name, 'grid_xt', 'grid_yt']:
                    #print("\r Copying over: %s..."%(ivar), end="")
                    prCyan(f"Copying over: {ivar}...")
                    fnew.copy_Ncvar(fNcdf.variables[ivar])

        print("\r ", end="")
        fNcdf.close()
        fnew.close()
        print(f"Completed in {(time.time() - start_time):3f} sec")

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    main()
