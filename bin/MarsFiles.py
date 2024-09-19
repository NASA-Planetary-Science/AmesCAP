#!/usr/bin/env python3
"""
The MarsFiles executable has functions for manipulating entire files.
The capabilities include time-shifting, binning, and regridding data,
as well as band pass filtering, tide analysis, zonal averaging, and
extracting variables from files.

The executable requires:
    * ``[input_file]``                  the file for manipulation

and optionally accepts:
    * ``[-fv3, --fv3]``                 produce MGCM ``fixed``,
        ``diurn``, ``average`` and ``daily`` files from Legacy output
    * ``[-c, --combine]``               Combine sequential files of
        the same type into one file
    * ``[-t, --tshift]``                apply a time-shift to
        ``diurn`` files
    * ``[-ba, --bin_average]``          bin MGCM ``daily`` files like
        ``average`` files
    * ``[-bd, --bin_diurn]``            bin MGCM ``daily`` files like
        ``diurn`` files
    * ``[-hp, --high_pass_filter]``     temporal filtering: high-pass
    * ``[-lp, --low_pass_filter]``      temporal filtering: low-pass
    * ``[-bp, --band_pass_filter]``     temporal filtering: band-pass
    * ``[-no_trend, --no_trend]``       filter and compute amplitudes
        only (use with filtering)
    * ``[-hpk, --high_pass_zonal]``     spatial filtering: high-pass
    * ``[-lpk, --low_pass_zonal]``      spatial filtering: low-pass
    * ``[-bpk, --band_pass_zonal]``     spatial filtering: band-pass
    * ``[-tidal, --tidal]``             extracts diurnal tide and its
        harmonics
    * ``[-reconstruct, --reconstruct]`` reconstructs the first N
        harmonics
    * ``[-norm, --normalize]``          provides ``-tidal`` result in
        percent amplitude
    * ``[-rs, --regrid_source]``        regrid a target file to match
        a source file
    * ``[-za, --zonal_avg]``            zonally average all variables
        in a file
    * ``[-include, --include]``         only include specific
        variables from the target file
    * ``[-e, --ext]``                   create a new file with a
        unique extension instead of overwriting current file

Third-party Requirements:
    * ``numpy``
    * ``netCDF4``
    * ``sys``
    * ``argparse``
    * ``os``
    * ``subprocess``
    * ``warnings``
"""

# Make print statements appear in color
from amescap.Script_utils import (Yellow, Cyan, Red, Blue, Yellow, Nclr, Green)

# Load generic Python Modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import subprocess   # Run command-line commands
import warnings     # Suppress errors triggered by NaNs
import numpy as np
from netCDF4 import Dataset

# Load amesCAP modules
from amescap.Ncdf_wrapper import (Ncdf, Fort)
from amescap.FV3_utils import (tshift, daily_to_average, daily_to_diurn)
from amescap.Script_utils import (find_tod_in_diurn, FV3_file_type, 
                                  filter_vars, regrid_Ncfile,
                                  get_longname_units, extract_path_basename)

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}MarsFiles is a file manager. Use it to modify a "
        f"netCDF file format.{Nclr}\n\n"
    ),
    formatter_class = argparse.RawTextHelpFormatter
)

parser.add_argument("input_file", nargs="+",
    help=(
        f"A netCDF file or list of netCDF files.\n\n"
    )
)

parser.add_argument("-fv3", "--fv3", nargs="+",
    help=(
        f"Produce MGCM ``fixed``, ``diurn``, ``average`` and "
        f"``daily`` files from Legacy output.\n"
        f"Available options are:\n"
        f"  - ``fixed``  : static fields (e.g., topography)\n"
        f"  - ``average``: 5-sol averages\n"
        f"  - ``daily``  : 5-sol continuous\n"
        f"  - ``diurn``  : 5-sol averages for each time of day\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py filename.nc -fv3 fixed\n"
        f"> MarsFiles.py filename.nc -fv3 fixed diurn"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-c", "--combine", action="store_true",
    help=(
        f"Combine sequential files of the same type into one file.\n"
        f"Works with all file types (``fixed``, ``average``, "
        f"``daily`` and ``diurn``).\n"
        f"{Yellow}Overwrites the first file in the series. "
        f"To override, use --ext.{Nclr}\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_average.nc --combine"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-split", "--split", nargs="+",
    help=(
        f"Extract values between min and max solar longitudes 0-360 [°]\n"
        f"This assumes all values in the file are from only one Mars Year.\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py 00668.atmos_average.nc --split 0 90 \n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-dim", "--dim", type=str, default = 'areo',
    help=(
        f"Flag to specify dimension to split on. Acceptable values are \n"
        f"areo, lat, lon, lev. For use with --split.\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py 00668.atmos_average.nc --split 0 90 --dim areo"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-t", "--tshift", nargs="?", const=999, type=str,
    help=(
    f"Apply a time-shift to {Yellow}``diurn``{Nclr}  files.\n"
        "Vertically interpolated ``diurn`` files OK.\n"
        f"{Yellow}Generates a new file ending in ``_T.nc``{Nclr}\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_diurn.nc --tshift\n"
        f"  {Blue}(outputs data for all 24 local times){Green}\n"
        f"> MarsFiles.py *.atmos_diurn.nc --tshift ``3 15``"
        f"\n"
        f"  {Blue}(outputs data for target local times only)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-ba", "--bin_average", nargs="?", const=5,type=int,
    help=(
        f"Bin MGCM ``daily`` files like ``average`` files.\n"
        f"{Yellow}Generates a new file ending in ``_to_average.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -ba\n"
        f"  {Blue}(NC, bin 5 days){Green}\n"
        f"> MarsFiles.py *.atmos_daily_pstd.nc -ba 10\n"
        f"  {Blue}(bin 10 days)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-bd", "--bin_diurn", action="store_true",
    help=(
        f"Bin MGCM ``daily`` files like ``diurn`` files.\n"
        f"May be used jointly with --bin_average.\n"
        f"{Yellow}Generates a new file ending in ``_to_diurn.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -bd\n"
        f"  {Blue}(default 5-day bin){Green}\n"
        f"> MarsFiles.py *.atmos_daily_pstd.nc -bd -ba 10\n"
        f"  {Blue}(10-day bin){Green}\n"
        f"> MarsFiles.py *.atmos_daily_pstd.nc -bd -ba 1\n"
        f"  {Blue}(No binning. Mimics raw Legacy output)"
        f"{Nclr}\n\n"
    )
)


parser.add_argument("-hpf", "--high_pass_filter", nargs="+", type=float,
    help=(
        f"Temporal filtering utilities: low-, high-, and "
        f"band-pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``_hpf.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -hpf 10.\n"
        f"  {Blue}(-hpf) --high_pass_filter sol_min "
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-lpf", "--low_pass_filter", nargs="+", type=float,
    help=(
        f"Temporal filtering utilities: low-, high-, and "
        f"band-pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``_lpf.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -lpf 0.5\n"
        f"  {Blue}(-lpf) --low_pass_filter sol_max "
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-bpf", "--band_pass_filter", nargs="+",
    help=(
        f"Temporal filtering utilities: low-, high-, and "
        f"band-pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``bpf.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -hpf 0.5 10.\n"
        f"  {Blue}(-bpf) --band_pass_filter sol_min sol max "
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-no_trend", "--no_trend", action="store_true",
    help=(
        f"Filter and compute amplitudes only.\n"
        f"For use with temporal filtering utilities (``-lpf``, "
        f"``-hpf``, ``-bpf``).\n"
        f"{Yellow}Generates a new file ending in ``_no_trend.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -hpf 10. --no_trend\n"
        f"> MarsFiles.py *.atmos_daily.nc -lpf 0.5 --no_trend\n"
        f"> MarsFiles.py *.atmos_daily.nc -hpf 0.5 10. --no_trend"
        f"{Nclr}\n\n"
    )
)

# Decomposition in zonal harmonics, disabled for initial CAP release:
parser.add_argument("-hpk", "--high_pass_zonal", nargs="+", type=int,
    help=(
        f"Spatial filtering utilities: low-, high-, and "
        f"band pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``_hpk.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -hpk 10 --no_trend\n"
        f"      {Blue}(-hpk) --high_pass_zonal kmin "
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-lpk", "--low_pass_zonal", nargs="+", type=int,
    help=(
        f"Spatial filtering utilities: low-, high-, and "
        f"band pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``_lpk.nc``\n"
        f"{Green}Usage:\n"
        f"    > MarsFiles.py *.atmos_daily.nc -lpk 20 --no_trend\n"
        f"      {Blue}(-lpk) --low_pass_zonal kmax "
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-bpk", "--band_pass_zonal", nargs="+",
    help=(
        f"Spatial filtering utilities: low-, high-, and "
        f"band pass filters.\n"
        f"Use ``--no_trend`` to compute amplitudes only.\n"
        f"Data detrended before filtering.\n"
        f"{Yellow}Generates a new file ending in ``_bpk.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -bpk 10 20 --no_trend\n"
        f"      {Blue}(-bpk) --band_pass_zonal kmin kmax"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-tidal", "--tidal", nargs="+", type=int,
    help=(
        f"Performs a tidal analyis on ``diurn`` files.\n"
        f"Extracts diurnal tide and its harmonics.\n"
        f"N = 1 diurnal, N = 2 semi-diurnal etc.\n"
        f"{Yellow}Generates a new file ending in ``_tidal.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_diurn.nc -tidal 4\n"
        f"  {Blue}(extracts 4 harmonics)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-reconstruct", "--reconstruct", action="store_true",
    help=(
        f"Reconstructs the first N harmonics.\n"
        f"{Yellow}Generates a new file ending in ``_reconstruct.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_diurn.nc -tidal 6 "
        f"--include ps temp --reconstruct"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-norm", "--normalize", action="store_true",
    help=(
        f"Provides result in percent amplitude.\n"
        f"{Yellow}Generates a new file ending in ``_norm.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_diurn.nc -tidal 6 "
        f"--include ps --normalize"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-rs", "--regrid_source", nargs="+",
    help=(
        f"Regrid a target file to match a source file.\n"
        f"Both source and target files should be vertically\n"
        f"interpolated to the same standard grid\n"
        f"(e.g. zstd, zagl, pstd, etc.).\n"
        f"{Yellow}Generates a new file ending in ``_regrid.nc``\n"
        f"{Green}Usage:\n"
        f"> MarsInterp.py *.atmos.average_pstd.nc -rs "
        f"simu2/00668.atmos_average_pstd.nc"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-za", "--zonal_avg", action="store_true",
    help=(
        f"Zonally average all variables in a file.\n"
        f"{Yellow}Generates a new file ending in ``_zonal_avg.nc``\n"
        f"{Green}Usage:\n"
        "> MarsFiles.py *.atmos_diurn.nc -za"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-include", "--include", nargs="+",
    help=(
        f"Flag to include only the variables listed after \n"
        f"-include in the target file.\n"
        f"All dimensional and 1D variables are always included.\n"
        f"{Yellow}Overwrites existing target file. To override, "
        f"use --ext.{Nclr}\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos_daily.nc -ba --include ps ts ucomp"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-e", "--ext", type=str, default = None,
    help=(
        f"Do not overwrite file. Append the extension provided \n"
        f"after --ext to the new file.\n"
        f"{Green}Usage:\n"
        f"> MarsFiles.py *.atmos.average.nc --combine --ext _combined\n"
        f"  {Blue}(produces *.atmos.average_combined.nc)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("--debug", action="store_true",
    help=(f"Debug flag: release the exceptions.\n\n"))

# ======================================================================
#                               DEFINITIONS
# ======================================================================

def combine_files(file_list, full_file_list):
    """
    Concatenates sequential output files in chronological order.

    :param file_list: list of file names
    :type file_list: list
    :param full_file_list: list of file names and full paths
    :type full_file_list: list
    """
    print(f"{Yellow}Using internal method for concatenation{Nclr}")

    # For fixed files, deleting all but the first file has the same
    # effect as combining files
    num_files = len(full_file_list)
    if (file_list[0][5:] == ".fixed.nc" and num_files >= 2):
        rm_cmd = "rm -f "
        for i in range(1, num_files):
            # 1-N files ensures file number 0 is preserved
            rm_cmd += f" {full_file_list[i]}"
        p = subprocess.run(rm_cmd, universal_newlines = True, shell = True)
        print(f"{Cyan}Cleaned all but {file_list[0]}{Nclr}")
        exit()

    print(f"{Cyan}Merging {num_files} files starting with {file_list[0]}..."
          f"{Nclr}")

    if parser.parse_args().include:
        # Exclude variables NOT listed after --include
        f = Dataset(file_list[0], "r")
        exclude_list = filter_vars(f, parser.parse_args().include, 
                                   giveExclude = True)
        f.close()
    else:
        exclude_list = []

    # Create a temporary file ending in _tmp.nc to work in
    tmp_file = f"{full_file_list[0][:-3]}_tmp.nc"
    Log = Ncdf(tmp_file, "Merged file")
    Log.merge_files_from_list(full_file_list, exclude_var=exclude_list)
    Log.close()

    # ----- Delete the files that were used for combine -----

    # First, rename temporary file for the final merged file
    #   For Legacy netCDF files, rename using initial and end Ls
    #   For MGCM netCDF files, rename to the first file in the list
    if file_list[0][:12] == "LegacyGCM_Ls":
        ls_ini = file_list[0][12:15]
        ls_end = file_list[-1][18:21]
        merged_file = f"LegacyGCM_Ls{ls_ini}_Ls{ls_end}.nc"
    else:
        merged_file = full_file_list[0]

    # Second, delete the files that were combined.
    # Apply the new name created above
    rm_cmd = "rm -f "
    for file in full_file_list:
        rm_cmd += f" {file}"

    cmd_txt = f"mv {tmp_file} {merged_file}"
    p = subprocess.run(rm_cmd, universal_newlines = True, shell = True)
    p = subprocess.run(cmd_txt, universal_newlines = True, shell = True)
    print(f"{Cyan}{merged_file} was created from a merge{Nclr}")
    
    return

def split_files(file_list, split_dim):
    if split_dim not in ['time', 'areo', 'lat', 'lon']:
        print(f"{Red}Split dimension must be one of the following:"
              f"    time, areo, lat, lon{Nclr}")
        exit()
        
    if split_dim == 'areo':
        split_dim = 'time'
    bounds = np.asarray(parser.parse_args().split).astype(float)
    
    if len(np.atleast_1d(bounds)) != 2:
        print(f"{Red}Requires two values: [lower_bound] [upper_bound]{Nclr}")
        exit()
        
    # Add path unless full path is provided
    if not ("/" in file_list[0]):
        input_file_name = f"{data_dir}/{file_list[0]}"
    else:
        input_file_name = file_list[0]
    original_date = (input_file_name.split('.')[0]).split('/')[-1]
    
    fNcdf = Dataset(input_file_name, 'r', format = 'NETCDF4_CLASSIC')
    var_list = filter_vars(fNcdf, parser.parse_args().include)

    dim_var = fNcdf.variables[split_dim][:]
    
    # Get file type (diurn, average, daily, etc.)
    f_type, _ = FV3_file_type(fNcdf)

    if f_type == 'diurn': 
        if split_dim == 'time':
            # size = areo (133, 24, 1)
            bounds_limit = np.squeeze(fNcdf.variables['areo'][:, 0, :]) % 360
        else:
            bounds_limit = np.squeeze(fNcdf.variables[split_dim][:, 0])
    else:
        if split_dim == 'time':
            # size = areo (133, 1)
            bounds_limit = np.squeeze(fNcdf.variables['areo'][:]) % 360
        else:
            bounds_limit = np.squeeze(fNcdf.variables[split_dim][:])

    lower_bound = np.argmin(np.abs(bounds[0] - bounds_limit))
    upper_bound = np.argmin(np.abs(bounds[1] - bounds_limit))

    if lower_bound == upper_bound:
        print(f"{Red}Warning, requested {split_dim} min, max ({bounds[0]}, "
              f"{bounds[1]}) are out of file range ({split_dim} = "
              f"{bounds_limit[0]:.1f}-{bounds_limit[-1]:.1f})")
        exit()

    dim_out = dim_var[lower_bound:upper_bound]
    print(f"{Yellow}dim_var = {dim_var}")
    print(f"{Yellow}dim_out = {dim_out}")
    print(f"{Yellow}dim_out.ndim = {dim_out.ndim}{Nclr}")

    fpath, fname = extract_path_basename(input_file_name)
    if split_dim == 'time':
        output_file_name = (f"{fpath}/{int(dim_out[0]):05d}{fname[5:-3]}_"
                            f"Ls{int(bounds[0]):03d}_{int(bounds[1]):03d}.nc")
    elif split_dim == 'lat':
        new_bounds = [str(abs(int(b)))+"S" if b < 0 else str(int(b))+"N" for b in bounds]
        print(f"{Yellow}bounds = {bounds[0]} {bounds[1]}")
        print(f"{Yellow}new_bounds = {new_bounds[0]} {new_bounds[1]}")
        output_file_name = (f"{fpath}/{original_date}{fname[5:-3]}_{split_dim}"
                            f"{new_bounds[0]}_{new_bounds[1]}.nc")
    else:
        output_file_name = (f"{fpath}/{original_date}{fname[5:-3]}_{split_dim}"
                            f"{int(bounds[0]):03d}_{int(bounds[1]):03d}.nc")
    
    print(f"{Cyan}new filename = {output_file_name}")
    Log = Ncdf(output_file_name)
    
    Log.copy_all_dims_from_Ncfile(fNcdf, exclude_dim = [split_dim])
    
    if split_dim == 'time':
        Log.add_dimension(split_dim, None)
    else:
        Log.add_dimension(split_dim, len(dim_out))
        print(f"{Yellow}len(split_dim) = {len(dim_out)}{Nclr}")
    
    if split_dim == 'time':
        Log.log_axis1D(variable_name = 'time', 
                       DATAin = dim_out, 
                       dim_name = 'time', 
                       longname_txt = 'sol number',
                       units_txt = 'days since 0000-00-00 00:00:00', 
                       cart_txt = 'T')
    elif split_dim == 'lat':
        Log.log_axis1D(variable_name = 'lat', 
                       DATAin = dim_out, 
                       dim_name = 'lat', 
                       longname_txt = 'latitude',
                       units_txt = 'degrees_N', 
                       cart_txt = 'T')
    elif split_dim == 'lon':
        Log.log_axis1D(variable_name = 'lon', 
                       DATAin = dim_out, 
                       dim_name = 'lon', 
                       longname_txt = 'longitude',
                       units_txt = 'degrees_E', 
                       cart_txt = 'T')
    
    # Loop over all variables in the file
    for ivar in var_list:
        varNcf = fNcdf.variables[ivar]
        if split_dim in varNcf.dimensions and ivar != split_dim:  
            # ivar is a dim of ivar but ivar is not ivar
            print(f'{Cyan}Processing: {ivar}...{Nclr}')
            if split_dim == 'time':
                var_out = varNcf[lower_bound:upper_bound, ...]
            elif split_dim == 'lat' and varNcf.ndim == 5:
                var_out = varNcf[:, :, :, lower_bound:upper_bound, :]
            elif split_dim == 'lat' and varNcf.ndim == 4:
                var_out = varNcf[:, :, lower_bound:upper_bound, :]
            elif split_dim == 'lat' and varNcf.ndim == 3:
                var_out = varNcf[:, lower_bound:upper_bound, :]
            elif split_dim == 'lat' and varNcf.ndim == 2:
                var_out = varNcf[lower_bound:upper_bound, ...]
            elif split_dim == 'lon' and varNcf.ndim > 2:
                var_out = varNcf[..., lower_bound:upper_bound]
            elif split_dim == 'lon' and varNcf.ndim == 2:
                var_out = varNcf[lower_bound:upper_bound, ...]
            longname_txt, units_txt = get_longname_units(fNcdf, ivar)
            Log.log_variable(ivar, var_out, varNcf.dimensions,
                             longname_txt, units_txt)
        else:
            # ivar is ivar OR ivar is not a dim of ivar
            if (ivar in ['pfull', 'lat', 'lon', 'phalf', 'pk', 'bk',
                         'pstd', 'zstd', 'zagl', 'time'] and 
                ivar != split_dim):
                # ivar is a dimension
                print(f'{Cyan}Copying axis: {ivar}...{Nclr}')
                Log.copy_Ncaxis_with_content(fNcdf.variables[ivar])
            elif ivar != split_dim:
                # ivar is not itself and not a dimension of ivar
                print(f'{Cyan}Copying variable: {ivar}...{Nclr}')
                Log.copy_Ncvar(fNcdf.variables[ivar])
    Log.close()
    fNcdf.close()
    return

# ==================================================================
#                   Time-Shifting Implementation
#                            Victoria H.
# ==================================================================

def time_shift(file_list):
    """
    This function converts the data in diurn files with a time_of_day_XX
    dimension to universal local time.

    :param file_list: list of file names
    :type file_list: list
    """
    if parser.parse_args().tshift == 999:
        # Target local times requested by user
        target_list = None 
    else:
        target_list = np.fromstring(parser.parse_args().tshift, 
                                    dtype = float, sep=" ")

    for file in file_list:
        # Add path unless full path is provided
        if not ("/" in file):
            input_file_name = f"{data_dir}/{file}"
        else:
            input_file_name = file
        output_file_name = f"{input_file_name[:-3]}_T.nc"

        if parser.parse_args().ext:
            # Append extension, if any:
            output_file_name = (f"{output_file_name[:-3]}_"
                                f"{parser.parse_args().ext}.nc")

        fdiurn = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")
        # Define a netcdf object from the netcdf wrapper module
        fnew = Ncdf(output_file_name)
        # Copy some dimensions from the old file to the new file
        fnew.copy_all_dims_from_Ncfile(fdiurn)

        # Find the "time of day" variable name
        tod_name_in = find_tod_in_diurn(fdiurn)
        _, zaxis = FV3_file_type(fdiurn)

        # Copy some variables from the old file to the new file
        fnew.copy_Ncaxis_with_content(fdiurn.variables["lon"])
        fnew.copy_Ncaxis_with_content(fdiurn.variables["lat"])
        fnew.copy_Ncaxis_with_content(fdiurn.variables["time"])
        fnew.copy_Ncaxis_with_content(fdiurn.variables["scalar_axis"])

        # Only create a vertical axis if orig. file has 3D fields
        if zaxis in fdiurn.dimensions.keys():
            fnew.copy_Ncaxis_with_content(fdiurn.variables[zaxis])

        # Take care of TOD dimension in new file
        tod_orig = np.array(fdiurn.variables[tod_name_in])

        if target_list is None:
            # If user does not specify which TOD(s) to do, do all 24
            tod_name_out = tod_name_in
            fnew.copy_Ncaxis_with_content(fdiurn.variables[tod_name_in])
            # Only copy "areo" if it exists in the original file
            if "areo" in fdiurn.variables.keys():
                fnew.copy_Ncvar(fdiurn.variables["areo"])
        else:
            # If user requests specific local times, update the old
            # axis as necessary
            tod_name_out = f"time_of_day_{(len(target_list)):02}"
            fnew.add_dim_with_content(tod_name_out, target_list,
                                      longname_txt = "time of day",
                                      units_txt = ("[hours since 0000-00-00 "
                                                   "00:00:00]"),
                                      cart_txt = "")
            # Create areo variable with the new size
            areo_shape = fdiurn.variables["areo"].shape
            areo_dims = fdiurn.variables["areo"].dimensions

            # Update shape with new time_of_day
            areo_shape = (areo_shape[0], len(target_list), areo_shape[2])
            areo_dims = (areo_dims[0], tod_name_out, areo_dims[2])

            areo_out = np.zeros(areo_shape)

            for i in range(len(target_list)):
                # For new target_list, e.g [3,15]
                # Get the closest "time_of_day" index
                j = np.argmin(np.abs(target_list[i] - tod_orig))
                areo_out[:, i, 0] = fdiurn.variables["areo"][:, j, 0]

            fnew.add_dim_with_content("scalar_axis", [0], longname_txt = "none",
                                      units_txt = "none")
            fnew.log_variable("areo", areo_out, areo_dims, "areo", "degrees")

        # Read in 4D field(s) and do the time shift. Exclude vars
        # not listed after --include in var_list
        lons = np.array(fdiurn.variables["lon"])
        var_list = filter_vars(fdiurn, parser.parse_args().include)

        for var in var_list:
            print(f"{Cyan}Processing: {var}...{Nclr}")
            value = fdiurn.variables[var][:]
            dims = fdiurn.variables[var].dimensions
            longname_txt, units_txt = get_longname_units(fdiurn, var)
            
            if (len(dims) >= 4):
                y = dims.index("lat")
                x = dims.index("lon")
                t = dims.index("time")
                tod = dims.index(tod_name_in)

            if (len(dims) == 4):
                # time, tod, lat, lon
                var_val_tmp = np.transpose(value, (x, y, t, tod))
                var_val_T = tshift(var_val_tmp, lons, tod_orig,
                                   timex = target_list)
                var_out = np.transpose(var_val_T, (2, 3, 1, 0))
                fnew.log_variable(var, var_out,
                                  ["time", tod_name_out, "lat", "lon"],
                                  longname_txt, units_txt)
            if (len(dims) == 5):
                # time, tod, Z, lat, lon
                z = dims.index(zaxis)
                var_val_tmp = np.transpose(value, (x, y, z, t, tod))
                var_val_T = tshift(var_val_tmp, lons, tod_orig,
                                   timex = target_list)
                var_out = np.transpose(var_val_T, (3, 4, 2, 1, 0))
                fnew.log_variable(var, var_out,
                                  ["time", tod_name_out, zaxis, "lat", "lon"],
                                  longname_txt, units_txt)
        fnew.close()
        fdiurn.close()
    return

# ======================================================================
#                           MAIN PROGRAM
# ======================================================================

def main():
    global data_dir
    file_list = parser.parse_args().input_file
    data_dir = os.getcwd()

    # Make a list of input files including the full path to the dir
    full_file_list = []
    for file in file_list:
        if not ("/" in file):
            full_file_list.append(f"{data_dir}/{file}")
        else:
            full_file_list.append(f"{file}")
            
    if parser.parse_args().fv3 and parser.parse_args().combine:
        print(f"{Red}Use --fv3 and --combine separately to avoid ambiguity")
        exit()

    # ==================================================================
    #               Conversion Legacy -> MGCM Format
    #                    Richard U. and Alex. K.
    # ==================================================================

    # Convert to MGCM Output Format
    if parser.parse_args().fv3:
        for req_file in parser.parse_args().fv3:
            if req_file not in ["fixed", "average", "daily", "diurn"]:
                print(f"{Red}{req_file} is invalid. Select ``fixed``, "
                      f"``average``, ``daily``, or ``diurn``{Nclr}")

        # lsmin = None
        # lsmax = None

        if full_file_list[0][-3:] == ".nc":
            print("Processing Legacy MGCM netCDF files")
            for f in full_file_list:
                # file_name = os.path.basename(f)
                # ls_l = file_name[-12:-9]
                # ls_r = file_name[-6:-3]

                # if lsmin is None:
                #     lsmin = ls_l
                # else:
                #     lsmin = str(min(int(lsmin), int(ls_l))).zfill(3)
                # if lsmax is None:
                #     lsmax = ls_r
                # else:
                #     lsmax = str(max(int(lsmax), int(ls_r))).zfill(3)
                make_FV3_files(f, parser.parse_args().fv3, True)
        else:
            print("Processing fort.11 files")
            for f in full_file_list:
                file_name = Fort(f)
                if "fixed" in parser.parse_args().fv3:
                    file_name.write_to_fixed()
                if "average" in parser.parse_args().fv3:
                    file_name.write_to_average()
                if "daily" in parser.parse_args().fv3:
                    file_name.write_to_daily()
                if "diurn" in parser.parse_args().fv3:
                    file_name.write_to_diurn()
                    
    elif parser.parse_args().combine:
        # Combine files along the time dimension
        combine_files(file_list, full_file_list)
    
    elif parser.parse_args().split:
        # Split file along the specified dimension. If none specified,
        # default to time dimension
        split_files(file_list, parser.parse_args().dim)
        
    elif parser.parse_args().tshift:
        # Time-shift files
        time_shift(file_list)

    # ==================================================================
    #               Bin a daily file as an average file
    #                               Alex K.
    # ==================================================================
    elif (parser.parse_args().bin_average and not
          parser.parse_args().bin_diurn):

        # Generate output file name
        bin_period = parser.parse_args().bin_average
        for file in file_list:
            # Add path unless full path is provided
            if not ("/" in file):
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}_to_average.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (
                    f"{output_file_name[:-3]}_{parser.parse_args().ext}.nc"
                )

            fdaily = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, parser.parse_args().include)

            time = fdaily.variables["time"][:]


            time_increment = time[1] - time[0]
            dt_per_day = int(np.round(1/time_increment))
            dt_total = int(dt_per_day * bin_period)

            bins = len(time) / (dt_per_day * bin_period)
            bins_even = len(time) // dt_total
            bins_left = len(time) % dt_total

            if bins_left != 0:
                print(f"{Yellow}*** Warning *** Requested a {bin_period}-sol "
                         f"bin period but the file has a total of {len(time)}"
                         f"timesteps ({dt_per_day} per sol) and {len(time)}/"
                         f"({bin_period}x{dt_per_day})={bins} is not a round "
                         f"number.\nWill use {bins_even} bins with "
                         f"{bin_period}x{dt_per_day}={dt_total} timesteps per "
                         f"bin ({bins_even*dt_total} timsteps total) and "
                         f"discard {bins_left} timesteps.{Nclr}")

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["time"])

            # Calculate and log the new time array
            fnew.add_dimension("time", None)
            time_out = daily_to_average(time[:], time_increment, bin_period)
            fnew.log_axis1D("time", time_out, "time", 
                            longname_txt = "sol number",
                            units_txt = "days since 0000-00-00 00:00:00", 
                            cart_txt = "T")

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]

                if "time" in varNcf.dimensions:
                    print(f"{Cyan}Processing: {ivar}{Nclr}")
                    var_out = daily_to_average(varNcf[:], time_increment, 
                                               bin_period)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(ivar, var_out, varNcf.dimensions,
                                      longname_txt, units_txt)

                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ==================================================================
    #               Bin a daily file as a diurn file
    #                               Alex K.
    # ==================================================================
    elif parser.parse_args().bin_diurn:
        # Use defaut binning period of 5 days (like average files)
        if parser.parse_args().bin_average is None:
            bin_period = 5
        else:
            bin_period = parser.parse_args().bin_average

        for file in file_list:
            # Add path unless full path is provided
            if not ("/" in file):
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}_to_diurn.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (f"{output_file_name[:-3]}_"
                                    f"{parser.parse_args().ext}.nc")

            fdaily = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, parser.parse_args().include)

            time = fdaily.variables["time"][:]

            time_increment = time[1] - time[0]
            dt_per_day = int(np.round(1/time_increment))

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but "time" from the old file to the
            # new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["time"])

            # If no binning is requested, copy time axis as-is
            fnew.add_dimension("time", None)
            time_out = daily_to_average(time[:], time_increment, bin_period)
            fnew.add_dim_with_content("time", time_out, 
                                      longname_txt = "sol number",
                                      units_txt = ("days since 0000-00-00 "
                                                   "00:00:00"), 
                                      cart_txt = "T")

            # Create a new time_of_day dimension
            tod_name = f"time_of_day_{dt_per_day:02}"
            time_tod = np.squeeze(daily_to_diurn(time[0:dt_per_day],
                                                 time[0:dt_per_day]))
            tod = np.mod(time_tod*24, 24)
            fnew.add_dim_with_content(tod_name, tod, 
                                      longname_txt = "time of day",
                                      units_txt = ("hours since 0000-00-00 "
                                                   "00:00:00"), 
                                      cart_txt = "N")

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]
                
                if "time" in varNcf.dimensions and ivar != "time":
                    # If time is the dimension (not just an array)
                    print(f"{Cyan}Processing: {ivar}{Nclr}")
                    dims_in = varNcf.dimensions
                    dims_out = (dims_in[0],)+(tod_name,)+dims_in[1:]
                    var_out = daily_to_diurn(varNcf[:], time[0:dt_per_day])
                    if bin_period != 1:
                        # dt is 1 sol between two diurn timesteps
                        var_out = daily_to_average(var_out, 1., bin_period)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(ivar, var_out, dims_out,
                                      longname_txt, units_txt)

                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    elif ivar != "time":
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ==================================================================
    #                       Transient Wave Analysis
    #                       Alex K. & R. J. Wilson
    # ==================================================================

    elif (parser.parse_args().high_pass_filter or
          parser.parse_args().low_pass_filter or
          parser.parse_args().band_pass_filter):

        # This functions requires scipy > 1.2.0. Import package here.
        from amescap.Spectral_utils import zeroPhi_filter

        if parser.parse_args().high_pass_filter:
            btype = "high"
            out_ext = "_hpf"
            nsol = np.asarray(
                parser.parse_args().high_pass_filter
                ).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                print(f"{Red}***Error*** sol_min accepts only one value")
                exit()
        if parser.parse_args().low_pass_filter:
            btype = "low"
            out_ext = "_lpf"
            nsol = np.asarray(
                parser.parse_args().low_pass_filter
                ).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                print(f"{Red}sol_max accepts only one value")
                exit()
        if parser.parse_args().band_pass_filter:
            btype = "band"
            out_ext = "_bpf"
            nsol = np.asarray(
                parser.parse_args().band_pass_filter
                ).astype(float)
            if len(np.atleast_1d(nsol)) != 2:
                print(f"{Red}Requires two values: sol_min sol_max")
                exit()
        if parser.parse_args().no_trend:
            out_ext = f"{out_ext}_no_trend"

        for file in file_list:
            if not ("/" in file):
                # Add path unless full path is provided
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}{out_ext}.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (
                    f"{output_file_name[:-3]}_{parser.parse_args().ext}.nc"
                )

            fdaily = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")

            var_list = filter_vars(fdaily, parser.parse_args().include)

            time = fdaily.variables["time"][:]

            dt = time[1]-time[0]

            # Check if the frequency domain is allowed
            if any(nn <= 2*dt for nn in nsol):
                print(f"{Red}***Error***  minimum cut-off cannot be smaller "
                      f"than the Nyquist period of 2xdt={2*dt} sol{Nclr}")
                exit()

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the
            # new file
            fnew.copy_all_dims_from_Ncfile(fdaily)

            if btype == "low":
                fnew.add_constant("sol_max", nsol,
                                  "Low-pass filter cut-off period ",
                                  "sol")
            elif btype == "high":
                fnew.add_constant("sol_min", nsol,
                                  "High-pass filter cut-off period ",
                                  "sol")
            elif btype == "band":
                fnew.add_constant("sol_min", nsol[0],
                                  "High-pass filter low cut-off period ",
                                  "sol")
                fnew.add_constant("sol_max", nsol[1],
                                  "High-pass filter high cut-off period ",
                                  "sol")
            dt = time[1]-time[0]

            fs = 1/(dt) # Frequency in sol-1
            if btype == "band":
                # Flip the sols so that the low frequency comes first
                low_highcut = 1/nsol[::-1]
            else:
                low_highcut = 1./nsol

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]

                if ("time" in varNcf.dimensions and
                    ivar not in ["time", "areo"]):
                    print(f"{Cyan}Processing: {ivar}{Nclr}")
                    var_out = zeroPhi_filter(
                        varNcf[:], btype, low_highcut, fs, axis = 0, order = 4,
                        no_trend = parser.parse_args().no_trend)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(ivar, var_out, varNcf.dimensions,
                                      longname_txt, units_txt)
                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ==================================================================
    #                       Zonal Decomposition Analysis
    #                               Alex K.
    # ==================================================================
    # elif (parser.parse_args().high_pass_zonal or 
    #       parser.parse_args().low_pass_zonal or 
    #       parser.parse_args().band_pass_zonal):
    #     # This function requires scipy > 1.2.0. Import the package here
    #     from amescap.Spectral_utils import zonal_decomposition, zonal_construct
    #     # Load the module
    #     # init_shtools()
    #     if parser.parse_args().high_pass_zonal:
    #         btype = "high"
    #         out_ext = "_hpk"
    #         nk = np.asarray(parser.parse_args().high_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk)) != 1:
    #             print(f"{Red}***Error*** kmin accepts only one value")
    #             exit()
    #     if parser.parse_args().low_pass_zonal:
    #         btype = "low"
    #         out_ext = "_lpk"
    #         nk = np.asarray(parser.parse_args().low_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk)) != 1:
    #             print(f"{Red}kmax accepts only one value")
    #             exit()
    #     if parser.parse_args().band_pass_zonal:
    #         btype = "band"
    #         out_ext = "_bpk"
    #         nk = np.asarray(parser.parse_args().band_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk)) != 2:
    #             print(f"{Red}Requires two values: kmin kmax")
    #             exit()
    
    #     if parser.parse_args().no_trend:
    #         out_ext = f"{out_ext}_no_trend"
    
    #     for file in file_list:
    #         # Add path unless full path is provided
    #         if not ("/" in file):
    #             input_file_name = f"{data_dir}/{file}"
    #         else:
    #             input_file_name=file
    #         output_file_name = f"{input_file_name[:-3]}{out_ext}.nc"
    
    #         # Append extension, if any:
    #         if parser.parse_args().ext:
    #             output_file_name = (f"{output_file_name[:-3]}_"
    #                                 f"{parser.parse_args().ext}.nc")
    
    #         fname = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")
    #         # Get all variables
    #         var_list = filter_vars(fname,parser.parse_args().include) 
    #         lon = fname.variables["lon"][:]
    #         lat = fname.variables["lat"][:]
    #         LON, LAT = np.meshgrid(lon,lat)
    
    #         dlat = lat[1] - lat[0]
    #         dx = 2*np.pi*3400
    
    #         # Check if the frequency domain is allowed and display some 
    #         # information
    #         if any(nn > len(lat)/2 for nn in nk):
    #             print(f"{Red}***Warning***  maximum wavenumber cut-off cannot "
    #                   f"be larger than the Nyquist criteria of "
    #                   f"``nlat/2 = ``{len(lat)/2)} sol{Nclr}")
    #         elif btype == "low":
    #             L_max = (1./nk) * dx
    #             print(f"{Yellow}Low pass filter, allowing only wavelength > "
    #                   f"{L_max} km{Nclr}")
    #         elif btype == "high":
    #             L_min = (1./nk) * dx
    #             print(f"{Yellow}High pass filter, allowing only wavelength < "
    #                   f"{L_min} km{Nclr}")
    #         elif btype == "band":
    #             L_min = (1. / nk[1]) * dx
    #             L_max = 1. / max(nk[0], 1.e-20) * dx
    #             if L_max > 1.e20:
    #                 L_max = np.inf
    #             print(f"{Yellow}Band pass filter, allowing only {L_min} km < "
    #                   f"wavelength < {L_max} km{Nclr}")
                
    #         # Define a netcdf object from the netcdf wrapper module
    #         fnew = Ncdf(output_file_name) 
    #         # Copy all dimensions but "time" from old -> new file
    #         fnew.copy_all_dims_from_Ncfile(fname)
    
    #         if btype == "low":
    #             fnew.add_constant("kmax", nk, 
    #                               "Low-pass filter zonal wavenumber ", 
    #                               "wavenumber")
    #         elif btype == "high":
    #             fnew.add_constant("kmin", nk, 
    #                               "High-pass filter zonal wavenumber ", 
    #                               "wavenumber")
    #         elif btype == "band":
    #             fnew.add_constant("kmin", nk[0], 
    #                               "Band-pass filter low zonal wavenumber ", 
    #                               "wavenumber")
    #             fnew.add_constant("kmax", nk[1], 
    #                               "Band-pass filter high zonal wavenumber ", 
    #                               "wavenumber")
    #         low_highcut = nk
                
    #         for ivar in var_list:
    #             # Loop over all variables in the file
    #             varNcf = fname.variables[ivar]
    
    #             if ("lat" in varNcf.dimensions and 
    #                 "lon" in varNcf.dimensions):
    #                 print(f"{Cyan}Processing: {ivar}...{Nclr}")
    #                 # Step 1 : Detrend the data
    #                 TREND = get_trend_2D(varNcf[:], LON, LAT,  "wmean")
    #                 # Step 2 : Calculate spherical harmonic coeffs
    #                 COEFF, PSD = zonal_decomposition(varNcf[:] - TREND)
    #                 # Step 3 : Recompose the variable out of the coeffs
    #                 VAR_filtered=zonal_construct(COEFF, varNcf[:].shape, 
    #                                              btype = btype, 
    #                                              low_highcut = low_highcut)
    #                 #Step 4: Add the trend, if requested
    #                 if parser.parse_args().no_trend:
    #                     var_out = VAR_filtered
    #                 else:
    #                     var_out = VAR_filtered + TREND
    
    #                 fnew.log_variable(ivar, var_out, varNcf.dimensions, 
    #                                   varNcf.long_name, varNcf.units)
    #             else:
    #                 if  ivar in ["pfull", "lat", "lon", "phalf", "pk", "bk", 
    #                              "pstd", "zstd", "zagl", "time"]:
    #                     print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
    #                     fnew.copy_Ncaxis_with_content(fname.variables[ivar])
    #                 else:
    #                     print(f"{Cyan}Copying variable: {ivar}...{Nclr}")
    #                     fnew.copy_Ncvar(fname.variables[ivar])
    #         fnew.close()

    # ==================================================================
    #                           Tidal Analysis
    #                           Alex K. & R. J. Wilson
    # ==================================================================

    elif parser.parse_args().tidal:
        from amescap.Spectral_utils import diurn_extract, reconstruct_diurn
        N = parser.parse_args().tidal[0]
        if len(np.atleast_1d(N)) != 1:
            print(f"{Red}***Error*** N accepts only one value")
            exit()
        out_ext = "_tidal"
        if parser.parse_args().reconstruct:
            out_ext = f"{out_ext}_reconstruct"
        if parser.parse_args().normalize:
            out_ext = f"{out_ext}_norm"

        for file in file_list:
            # Add path unless full path is provided
            if not ("/" in file):
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}{out_ext}.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (f"{output_file_name[:-3]}_"
                                    f"{parser.parse_args().ext}.nc")

            fdiurn = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")

            var_list = filter_vars(fdiurn, parser.parse_args().include)

            # Find time_of_day variable name
            tod_name = find_tod_in_diurn(fdiurn)

            target_tod = fdiurn.variables[tod_name][:]
            lon = fdiurn.variables["lon"][:]
            areo = fdiurn.variables["areo"][:]

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dims but time_of_day from the old file to the
            # new file

            # Harmonics to reconstruct the signal. We use the original
            # time_of_day array.
            if parser.parse_args().reconstruct:
                fnew.copy_all_dims_from_Ncfile(fdiurn)
                # Copy time_of_day axis from initial files
                fnew.copy_Ncaxis_with_content(fdiurn.variables[tod_name])

            else:
                fnew.copy_all_dims_from_Ncfile(fdiurn, 
                                               exclude_dim = [tod_name])
                # Create new dimension holding the harmonics. We reuse
                # the time_of_day name to facilitate. Compatible with
                # other routines, but keep in mind this is the harmonic
                # number
                fnew.add_dim_with_content(
                    f"time_of_day_{N}",
                    np.arange(1, N+1),
                    longname_txt = "tidal harmonics",
                    units_txt = "Diurnal harmonic number",
                    cart_txt = "N"
                )

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdiurn.variables[ivar]
                varIN = varNcf[:]
                longname_txt, units_txt = get_longname_units(fdiurn, ivar)
                var_unit = getattr(varNcf, "units", "")

                if (tod_name in varNcf.dimensions and
                    ivar not in [tod_name, "areo"] and
                    len(varNcf.shape) > 2):
                    print(f"{Cyan}Processing: {ivar}{Nclr}")

                    # Normalize the data
                    if parser.parse_args().normalize:
                        # Normalize and reshape the array along the
                        # time_of_day dimension
                        norm = np.mean(varIN, axis = 1)[:, np.newaxis, ...]
                        varIN = 100*(varIN-norm)/norm
                        #units_txt = f"% of diurnal mean"
                        var_unit = f"% of diurnal mean"

                    amp, phas = diurn_extract(varIN.swapaxes(0, 1), N, 
                                              target_tod, lon)
                    if parser.parse_args().reconstruct:
                        VARN = reconstruct_diurn(amp, phas, target_tod, lon, 
                                                 sumList=[])
                        for nn in range(N):
                            fnew.log_variable(f"{ivar}_N{nn+1}",
                                              VARN[nn, ...].swapaxes(0, 1),
                                              varNcf.dimensions,
                                              (f"harmonic N={nn+1} for "
                                               f"{longname_txt}"),
                                              units_txt)

                    else:
                        # Update the dimensions
                        new_dim = list(varNcf.dimensions)
                        new_dim[1] = f"time_of_day_{N}"
                        fnew.log_variable(f"{ivar}_amp", amp.swapaxes(0,1),
                                          new_dim, 
                                          f"tidal amplitude for {longname_txt}",
                                          units_txt)
                        fnew.log_variable(f"{ivar}_phas", phas.swapaxes(0,1),
                                          new_dim,
                                          f"tidal phase for {longname_txt}",
                                          "hr")

                elif  ivar in ["pfull", "lat", "lon", "phalf", "pk",
                               "bk", "pstd", "zstd", "zagl", "time"]:
                        print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdiurn.variables[ivar])
                elif  ivar in ["areo"]:
                        if parser.parse_args().reconstruct:
                            #time_of_day is the same size as the
                            # original file
                            print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                            fnew.copy_Ncvar(fdiurn.variables["areo"])
                        else:
                            print(f"{Cyan}Processing: {ivar}...{Nclr}")
                            #Create areo variable reflecting the
                            # new shape
                            areo_new=np.zeros((areo.shape[0], N, 1))
                            #Copy areo
                            for xx in range(N):
                                areo_new[:, xx, :] = areo[:, 0, :]
                            #Update the dimensions
                            new_dim = list(varNcf.dimensions)
                            new_dim[1] = f"time_of_day_{N}"
                            #fnew.log_variable(ivar, bareo_new, new_dim,
                            # longname_txt, units_txt)
                            fnew.log_variable(ivar, areo_new, new_dim, 
                                              longname_txt, var_unit)
            fnew.close()

    # ==================================================================
    #                           Regridding Routine
    #                                 Alex K.
    # ==================================================================

    elif parser.parse_args().regrid_source:
        out_ext = "_regrid"
        name_target = parser.parse_args().regrid_source[0]

        # Add path unless full path is provided
        if not ("/" in name_target):
            name_target = f"{data_dir}/{name_target}"
        fNcdf_t = Dataset(name_target, "r")

        for file in file_list:
            # Add path unless full path is provided
            if not ("/" in file):
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}{out_ext}.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (f"{output_file_name[:-3]}_"
                                    f"{parser.parse_args().ext}.nc")

            f_in = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")

            var_list = filter_vars(
                f_in, parser.parse_args().include)  # Get all variables

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)

            # Copy all dims from the target file to the new file
            fnew.copy_all_dims_from_Ncfile(fNcdf_t)

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf     = f_in.variables[ivar]
                longname_txt,units_txt = get_longname_units(f_in, ivar)

                if  ivar in ["pfull", "lat", "lon", "phalf", "pk",
                             "bk", "pstd", "zstd", "zagl", "time", "areo"]:
                        print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                        fnew.copy_Ncaxis_with_content(fNcdf_t.variables[ivar])
                elif varNcf.dimensions[-2:]==("lat", "lon"):
                    #Ignore variables like time_bounds, scalar_axis
                    # or grid_xt_bnds...
                    print(f"{Cyan}Regridding: {ivar}...{Nclr}")
                    var_OUT = regrid_Ncfile(varNcf, f_in, fNcdf_t)
                    fnew.log_variable(ivar, var_OUT, varNcf.dimensions,
                                      longname_txt, units_txt)
            fnew.close()
            fNcdf_t.close()

    # ==================================================================
    #                           Zonal Averaging
    #                              Alex K.
    # ==================================================================
    elif parser.parse_args().zonal_avg:
        for file in file_list:
            if not ("/" in file):
                # Add path unless full path is provided
                input_file_name = f"{data_dir}/{file}"
            else:
                input_file_name = file
            output_file_name = f"{input_file_name[:-3]}_zonal_avg.nc"

            # Append extension, if any:
            if parser.parse_args().ext:
                output_file_name = (f"{output_file_name[:-3]}_"
                                    f"{parser.parse_args().ext}.nc")

            fdaily = Dataset(input_file_name, "r", format = "NETCDF4_CLASSIC")
            var_list = filter_vars(
                fdaily, parser.parse_args().include)  # Get all variables

            lon_in = fdaily.variables["lon"][:]

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the
            # new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["lon"])

            # Add a new dimension for the longitude, size = 1
            fnew.add_dim_with_content("lon", [lon_in.mean()], 
                                      longname_txt = "longitude", 
                                      units_txt = "degrees_E", cart_txt = "X")

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]
                longname_txt,units_txt = get_longname_units(fdaily,ivar)
                if ("lon" in varNcf.dimensions and
                    ivar not in ["lon", "grid_xt_bnds", "grid_yt_bnds"]):
                    print(f"{Cyan}Processing: {ivar}...{Nclr}")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore",
                                              category = RuntimeWarning)
                        var_out = np.nanmean(varNcf[:], 
                                             axis = -1)[..., np.newaxis]
                        fnew.log_variable(ivar, var_out, varNcf.dimensions,
                                          longname_txt, units_txt)
                else:
                    if ivar in ["pfull", "lat", "phalf", "pk", "bk", "pstd", 
                                "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    elif ivar in ["grid_xt_bnds", "grid_yt_bnds", "lon"]:
                        pass
                    else:
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()
    else:
        print(f"{Red}Error: no action requested: use ``MarsFiles *nc --fv3 "
              "--combine, --tshift, --bin_average, --bin_diurn etc ...``")
# END of script

# ======================================================
#                  DEFINITIONS
# ======================================================

def make_FV3_files(fpath, typelistfv3, renameFV3=True):
    """
    Make MGCM-like ``average``, ``daily``, and ``diurn`` files.
    Used if call to [``-fv3 --fv3``] is made AND Legacy files are in
    netCDFformat (not fort.11).

    :param fpath: Full path to the Legacy netcdf files
    :type fpath: str
    :param typelistfv3: MGCM-like file type: ``average``, ``daily``,
        or ``diurn``
    :type typelistfv3: list
    :param renameFV3: Rename the files from Legacy_LsXXX_LsYYY.nc to
        ``XXXXX.atmos_average.nc`` following MGCM output conventions
    :type renameFV3: bool

    :return: The MGCM-like files: ``XXXXX.atmos_average.nc``,
        ``XXXXX.atmos_daily.nc``, ``XXXXX.atmos_diurn.nc``.
    """
    historyDir = os.getcwd()
    histfile = Dataset(fpath, "r", format = "NETCDF4_CLASSIC")
    histdims = histfile.dimensions.keys()
    
    if renameFV3:
        # Convert the first Ls in file to a sol number
        fdate = f"{(ls2sol_1year(histfile.variables['ls'][0])):05}"

    def proccess_file(newf, typefv3):
        """
            Creates required variables and inputs them into the new
            files. Required variables include ``latitude``,
            ``longitude``, ``time``, ``time-of-day`` (if diurn file),
            and vertical layers (``phalf`` and ``pfull``).

            :param newf: path to target file
            :type newf: str
            :param typefv3: identifies type of file: ``average``,
                ``daily``, or ``diurn``
            :type typefv3: str

            :return: netCDF file with minimum required variables
        """
        for dname in histdims:
            if dname == "nlon":
                var = histfile.variables["longitude"]
                npvar = var[:]
                newf.add_dim_with_content("lon", npvar, "longitude",
                                          getattr(var, "units"), "X")
            elif dname == "nlat":
                var = histfile.variables["latitude"]
                npvar = var[:]
                newf.add_dim_with_content("lat", npvar, "latitude",
                                          getattr(var, "units"), "Y")

            elif dname == "time":
                newf.add_dimension("time", None)
            elif dname == "ntod" and typefv3 == "diurn":
                dim = histfile.dimensions[dname]
                newf.add_dimension("time_of_day_16", dim.size)
            elif dname == "nlay":
                nlay = histfile.dimensions[dname]
                num = nlay.size
                nump = num + 1
                pref = 7.01 * 100 # [Pa]
                pk = np.zeros(nump)
                bk = np.zeros(nump)
                pfull = np.zeros(num)
                phalf = np.zeros(nump)

                sgm = histfile.variables["sgm"]
                # Changed pk[0] = .08 to pk[0] = .08/2, otherwise
                # phalf[0] > phalf[1]
                pk[0] = 0.08/2
                # .. NOTE:: pk in amesCAP/mars_data/Legacy.fixed.nc
                # is also updated
                for z in range(num):
                    bk[z+1] = sgm[2*z + 2]
                phalf[:] = pk[:] + pref*bk[:] # [Pa]

                # DEPRECATED:
                # pfull[:] = (phalf[1:]-phalf[:num])/(np.log(phalf[1:])
                #             - np.log(phalf[:num]))

                # First layer:
                if pk[0] == 0 and bk[0] == 0:
                    pfull[0] = 0.5 * (phalf[0]+phalf[1])
                else:
                    pfull[0] = ((phalf[1]-phalf[0])
                                / (np.log(phalf[1]) - np.log(phalf[0])))
                # All other layers:
                pfull[1:] = ((phalf[2:]-phalf[1:-1])
                             / (np.log(phalf[2:]) - np.log(phalf[1:-1])))

                newf.add_dim_with_content("pfull", pfull,
                                          "ref full pressure level", "Pa")
                newf.add_dim_with_content("phalf", phalf,
                                          "ref half pressure level", "Pa")
                newf.log_axis1D("pk", pk, ("phalf"),
                                longname_txt = ("pressure part of the hybrid "
                                              "coordinate"),
                                units_txt = "Pa", cart_txt = "")
                newf.log_axis1D("bk", bk, ("phalf"),
                                longname_txt = ("sigma part of the hybrid "
                                              "coordinate"),
                                units_txt = "Pa", cart_txt = "")
            else:
                dim = histfile.dimensions[dname]
                newf.add_dimension(dname, dim.size)

        # =========== END function ===========

    if "average" in typelistfv3:
        # 5-sol average over "time_of_day" and "time"
        newfname_avg = f"{fdate}.atmos_average.nc"
        newfpath_avg = os.path.join(historyDir, newfname_avg)
        newfavg = Ncdf(newfpath_avg)
        proccess_file(newfavg, "average")
        do_avg_vars(histfile, newfavg, True, True)
        newfavg.close()

    if "daily" in typelistfv3:
        # Daily snapshot of the output
        newfname_daily = f"{fdate}.atmos_daily.nc"
        newfpath_daily = os.path.join(historyDir, newfname_daily)
        newfdaily = Ncdf(newfpath_daily)
        proccess_file(newfdaily, "daily")
        do_avg_vars(histfile, newfdaily, False, False)
        newfdaily.close()

    if "diurn" in typelistfv3:
        # 5-sol average over "time" only
        newfname_diurn = f"{fdate}.atmos_diurn.nc"
        newfpath_diurn = os.path.join(historyDir, newfname_diurn)
        newfdiurn = Ncdf(newfpath_diurn)
        proccess_file(newfdiurn, "diurn")
        do_avg_vars(histfile, newfdiurn, True, False)
        newfdiurn.close()

    if "fixed" in typelistfv3:
        # Copy Legacy.fixed to current directory
        cmd_txt = f"cp {sys.prefix}/mars_data/Legacy.fixed.nc {fdate}.fixed.nc"
        p = subprocess.run(cmd_txt, universal_newlines = True, shell = True)
        print(f"{os.getcwd()}/{fdate}.fixed.nc was copied locally")

def do_avg_vars(histfile, newf, avgtime, avgtod, bin_period=5):
    """
    Performs a time average over all fields in a file.

    :param histfile: file to perform time average on
    :type histfile: str
    :param newf: path to target file
    :type newf: str
    :param avgtime: whether ``histfile`` has averaged fields
        (e.g., ``atmos_average``)
    :type avgtime: bool
    :param avgtod: whether ``histfile`` has a diurnal time dimenion
        (e.g., ``atmos_diurn``)
    :type avgtod: bool
    :param bin_period: the time binning period if `histfile` has
        averaged fields (i.e., if ``avgtime==True``), defaults to 5
    :type bin_period: int, optional

    :return: a time-averaged file
    """
    histvars = histfile.variables.keys()
    for vname in histvars:
        var = histfile.variables[vname]
        npvar = var[:]
        dims = var.dimensions
        ndims = npvar.ndim
        vshape = npvar.shape
        ntod = histfile.dimensions["ntod"]

        # longname_txt, units_txt = get_longname_units(histfile, vname)
        longname_txt = getattr(histfile.variables[vname], "long_name", "")

        if longname_txt == "":
            # On some files, like the LegacyGCM_Ls*** on the NAS data
            # portal, the attribute long_name may be mispelled longname
            longname_txt = getattr(histfile.variables[vname], "longname", "")

        units_txt = getattr(histfile.variables[vname], "units", "")

        if avgtod:
            newdims = replace_dims(dims, True)
        elif avgtime:
            newdims = replace_dims(dims, False)
        else:
            newdims = replace_dims(dims, True)

        if "time" in dims:
            tind = dims.index("time")
            tind_new = newdims.index("time")
            numt = histfile.dimensions["time"].size
        # TODO fix time !!
        # Now do time averages and write to files
        if ndims == 1:
            if vname == "ls":
                if not np.all(npvar[1:] >= npvar[:-1]):
                    # If Ls crosses over into a new year
                    year = 0.
                    for x in range(1, npvar.size):
                        if (350. < npvar[x-1] < 360. and npvar[x] < 10.):
                            year += 1.
                        npvar[x] += 360. * year

                # Create a time array
                time0 = (ls2sol_1year(npvar[0])
                         + np.linspace(0, 10., len(npvar)))

                if avgtime:
                    varnew = np.mean(npvar.reshape(-1, bin_period), axis = 1)
                    time0 = np.mean(time0.reshape(-1, bin_period), axis = 1)

                if not avgtime and not avgtod:
                    # Daily file
                    # Solar longitude
                    ls_start = npvar[0]
                    ls_end = npvar[-1]
                    step = ((ls_end-ls_start) / np.float32(((numt-1) 
                                                            * ntod.size)))
                    varnew = np.arange(0, numt * ntod.size, dtype = np.float32)
                    varnew[:] = varnew[:]*step + ls_start

                    # Time
                    step = ((ls2sol_1year(ls_end) - ls2sol_1year(ls_start))
                            / np.float32((numt * ntod.size)))
                    time0 = np.arange(0, numt * ntod.size, dtype = np.float32)
                    time0[:] = time0[:]*step + ls2sol_1year(ls_start)

                newf.log_axis1D("areo", varnew, dims,
                                longname_txt = "solar longitude",
                                units_txt = "degree",
                                cart_txt = "T")
                newf.log_axis1D("time", time0, dims,
                                longname_txt = "sol number",
                                units_txt = "days since 0000-00-00 00:00:00",
                                cart_txt = "T")
            else:
                continue

        elif ndims == 4:
            varnew = npvar
            if avgtime:
                varnew = np.mean(npvar.reshape(-1, bin_period, vshape[1],
                                               vshape[2], vshape[3]), axis = 1)
            if avgtod:
                varnew = varnew.mean(axis = 1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3])
            # Rename variable
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt)
            # Convert surface pressure from mbar -> Pa
            if vname2 == "ps":
                varnew *= 100.
            newf.log_variable(vname2, varnew, newdims,longname_txt2,
                              units_txt2)
        elif ndims == 5:
            varnew = npvar
            if avgtime:
                varnew = np.mean(npvar.reshape(-1, bin_period, vshape[1], 
                                               vshape[2], vshape[3], 
                                               vshape[4]), axis = 1)
            if avgtod:
                varnew = varnew.mean(axis = 1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3], vshape[4])
            # Rename variables
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt)
            newf.log_variable(vname2, varnew, newdims,longname_txt2,
                              units_txt2)
        elif vname == "tloc":
            if avgtime and not avgtod:
                vname2 = "time_of_day_16"
                longname_txt2 = "time_of_day"
                units_txt2 = "hours since 0000-00-00 00:00:00"
                # Overwrite ``time_of_day`` from 
                # [``time_of_day_16``, ``lon``] -> ``time_of_day_16``
                newdims = ("time_of_day_16")
                # Every 1.5 hours, centered at half timestep
                npvar = np.arange(0.75, 24, 1.5)
                newf.log_variable(vname2, npvar, newdims, longname_txt2,
                                  units_txt2)
    return 0

def change_vname_longname_unit(vname, longname_txt, units_txt):
    """
    Update variable ``name``, ``longname``, and ``units``. This is
    designed to work specifically with LegacyCGM.nc files.

    :param vname: variable name
    :type vname: str
    :param longname_txt: variable description
    :type longname_txt: str
    :param units_txt: variable units
    :type units_txt: str

    :return: variable name and corresponding description and unit
    """
    if vname == "psurf":
        vname = "ps"
        longname_txt = "surface pressure"
        units_txt = "Pa"
    elif vname == "tsurf":
        vname = "ts"
        longname_txt = "surface temperature"
        units_txt = "K"
    elif vname == "dst_core_mass":
        vname = "cor_mass"
        longname_txt = "dust core mass for the water ice aerosol"
        units_txt = "kg/kg"
    elif vname == "h2o_vap_mass":
        vname = "vap_mass"
        longname_txt = "water vapor mixing ratio"
        units_txt = "kg/kg"
    elif vname == "h2o_ice_mass":
        vname = "ice_mass"
        longname_txt = "water ice aerosol mass mixing ratio"
        units_txt = "kg/kg"
    elif vname == "dst_mass":
        vname = "dst_mass"
        longname_txt = "dust aerosol mass mixing ratio"
        units_txt = "kg/kg"
    elif vname == "dst_numb":
        vname = "dst_num"
        longname_txt = "dust aerosol number"
        units_txt = "number/kg"
    elif vname == "h2o_ice_numb":
        vname = "ice_num"
        longname_txt = "water ice aerosol number"
        units_txt = "number/kg"
    elif vname == "temp":
        longname_txt = "temperature"
        units_txt = "K"
    elif vname == "ucomp":
        longname_txt = "zonal wind"
        units_txt = "m/s"
    elif vname == "vcomp":
        longname_txt = "meridional wind"
        units_txt = "m/s"
    else:
        # Return original values
        pass
    return vname, longname_txt, units_txt

def replace_dims(dims, todflag):
    """
    Replaces dimensions with MGCM-like names. Removes ``time_of_day``.
    This is designed to work specifically with LegacyCGM.nc files.

    :param dims: dimensions of the variable
    :type dims: str
    :param todflag: indicates whether there exists a ``time_of_day``
        dimension
    :type todflag: bool

    :return: new dimension names for the variable
    """
    newdims = dims
    if "nlat" in dims:
        newdims = replace_at_index(newdims, newdims.index("nlat"), "lat")
    if "nlon" in dims:
        newdims = replace_at_index(newdims, newdims.index("nlon"), "lon")
    if "nlay" in dims:
        newdims = replace_at_index(newdims, newdims.index("nlay"), "pfull")
    if "ntod" in dims:
        if todflag:
            newdims = replace_at_index(newdims, newdims.index("ntod"), None)
        else:
            newdims = replace_at_index(newdims, newdims.index("ntod"),
                                       "time_of_day_16")
    return newdims

def replace_at_index(tuple_dims, idx, new_name):
    """
    Updates variable dimensions.

    :param tuple_dims: the dimensions as tuples e.g. (``pfull``,
        ``nlat``, ``nlon``)
    :type tuple_dims: tuple
    :param idx: index indicating axis with the dimensions to update
        (e.g. ``idx = 1``  for ``nlat``)
    :type idx: int
    :param new_name: new dimension name (e.g. ``latitude``)
    :type new_name: str

    :return: updated dimensions
    """
    if new_name is None:
        return tuple_dims[:idx] + tuple_dims[idx+1:]
    else:
        return tuple_dims[:idx] + (new_name,) + tuple_dims[idx+1:]

def ls2sol_1year(Ls_deg, offset=True, round10=True):
    """
    Returns a sol number from the solar longitude.

    :param Ls_deg: solar longitude [°]
    :type Ls_deg: float
    :param offset: if True, force year to start at Ls 0
    :type offset: bool
    :param round10: if True, round to the nearest 10 sols
    :type round10: bool

    :returns: ``Ds`` the sol number

    .. NOTE:: For the moment, this is consistent with 0 <= Ls <=
        359.99, but not for monotically increasing Ls.
    """
    Lsp = 250.99    # Ls at perihelion
    tperi = 485.35  # Time (in sols) at perihelion
    Ns = 668.6      # Number of sols in 1 MY
    e = 0.093379    # From MGCM: modules.f90
    nu = (Ls_deg-Lsp)*np.pi/180
    E = 2 * np.arctan(np.tan(nu/2) * np.sqrt((1-e)/(1+e)))
    M = E - e*np.sin(E)
    Ds = M/(2*np.pi)*Ns + tperi

    if offset:
        # Offset correction:
        if len(np.atleast_1d(Ds)) == 1:
            # Ds is a float
            Ds -= Ns
            if Ds < 0:
                Ds += Ns
        else:
            # Ds is an array
            Ds -= Ns
            Ds[Ds < 0] = Ds[Ds < 0] + Ns
    if round:
        # -1 means round to the nearest 10
        Ds = np.round(Ds, -1)
    return Ds

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    main()
