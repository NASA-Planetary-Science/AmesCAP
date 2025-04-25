#!/usr/bin/env python3
"""
The MarsFiles executable has functions for manipulating entire files.
The capabilities include time-shifting, binning, and regridding data,
as well as band pass filtering, tide analysis, zonal averaging, and
extracting variables from files.

The executable requires:

    * ``[input_file]``                  The file for manipulation

and optionally accepts:

    * ``[-bin, --bin_files]``             Produce MGCM 'fixed', 'diurn', 'average' and 'daily' files from Legacy output
    * ``[-c, --concatenate]``             Combine sequential files of the same type into one file
    * ``[-split, --split]``               Split file along a specified dimension or extracts slice at one point along the dim
    * ``[-t, --time_shift]``              Apply a time-shift to 'diurn' files
    * ``[-ba, --bin_average]``            Bin MGCM 'daily' files like 'average' files
    * ``[-bd, --bin_diurn]``              Bin MGCM 'daily' files like 'diurn' files
    * ``[-hpt, --high_pass_temporal]``    Temporal filter: high-pass
    * ``[-lpt, --low_pass_temporal]``     Temporal filter: low-pass
    * ``[-bpt, --band_pass_temporal]``    Temporal filter: band-pass
    * ``[-trend, --add_trend]``           Return amplitudes only (use with temporal filters)
    * ``[-tide, --tide_decomp]``          Extract diurnal tide and its harmonics
    * ``[-recon, --reconstruct]``         Reconstruct the first N harmonics
    * ``[-norm, --normalize]``            Provide ``-tide`` result in % amplitude
    * ``[-regrid, --regrid_XY_to_match]`` Regrid a target file to match a source file
    * ``[-zavg, --zonal_average]``        Zonally average all variables in a file
    * ``[-incl, --include]``              Only include specific variables in a calculation
    * ``[-ext, --extension]``             Create a new file with a unique extension instead of modifying the current file

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
from amescap.Script_utils import (
    Yellow, Cyan, Red, Blue, Yellow, Nclr, Green
)

# Load generic Python Modules
import sys          # System commands
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import subprocess   # Run command-line commands
import warnings     # Suppress errors triggered by NaNs
import re           # Regular expressions
import numpy as np
from netCDF4 import Dataset
import tempfile     # Platform-independent temporary file creation

# Load amesCAP modules
from amescap.Ncdf_wrapper import (Ncdf, Fort)
from amescap.FV3_utils import (
    time_shift_calc, daily_to_average, daily_to_diurn, get_trend_2D
)
from amescap.Script_utils import (
    find_tod_in_diurn, FV3_file_type, filter_vars, regrid_Ncfile,
    get_longname_unit, check_bounds
)

# ------------------------------------------------------
#                  EXTENSION FUNCTION
# ------------------------------------------------------
class ExtAction(argparse.Action):
    def __init__(self, *args, ext_content=None, parser=None, **kwargs):
        self.parser = parser
        # Store the extension content that's specific to this argument
        self.ext_content = ext_content
        # For flags, we need to handle nargs=0 and default=False
        if kwargs.get('nargs') == 0:
            kwargs['default'] = False
            kwargs['const'] = True
            kwargs['nargs'] = 0
        super().__init__(*args, **kwargs)

        # Store this action in the parser's list of actions
        # We'll use this later to set up default info values
        if self.parser:
            if not hasattr(self.parser, '_ext_actions'):
                self.parser._ext_actions = []
            self.parser._ext_actions.append(self)

    def __call__(self, parser, namespace, values, option_string=None):
        # Handle flags (store_true type arguments)
        if self.nargs == 0:
            setattr(namespace, self.dest, True)
        # Handle other types
        elif isinstance(values, list):
            setattr(namespace, self.dest, values)
        elif isinstance(values, str) and ',' in values:
            # Handle comma-separated lists
            setattr(namespace, self.dest, values.split(','))
        else:
            try:
                # Try to convert to int if it's an integer argument
                setattr(namespace, self.dest, int(values))
            except ValueError:
                setattr(namespace, self.dest, values)

        # Set the extension content using the argument name
        ext_attr = f"{self.dest}_ext"
        setattr(namespace, ext_attr, self.ext_content)

class ExtArgumentParser(argparse.ArgumentParser):
    def parse_args(self, *args, **kwargs):
        namespace = super().parse_args(*args, **kwargs)

        # Then set info attributes for any unused arguments
        if hasattr(self, '_ext_actions'):
            for action in self._ext_actions:
                ext_attr = f"{action.dest}_ext"
                if not hasattr(namespace, ext_attr):
                    setattr(namespace, ext_attr, "")

        return namespace

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

# ------------------------------------------------------
#                  ARGUMENT PARSER
# ------------------------------------------------------

parser = ExtArgumentParser(
    prog=('MarsFiles'),
    description=(
        f"{Yellow}MarsFiles is a file manager. Use it to modify a "
        f"netCDF file format.{Nclr}\n\n"
    ),
    formatter_class = argparse.RawTextHelpFormatter
)

parser.add_argument('input_file', nargs='+',
    type=argparse.FileType('rb'),
    help=(f"A netCDF or fort.11 file or list of files.\n\n"))

parser.add_argument('-bin', '--bin_files', nargs='+', type=str,
    choices=['fixed', 'diurn', 'average', 'daily'],
    help=(
        f"Produce MGCM 'fixed', 'diurn', 'average' and 'daily' "
        f"files from Legacy output:\n"
        f"  - ``fixed``  : static fields (e.g., topography)\n"
        f"  - ``daily``  : instantaneous data\n"
        f"  - ``diurn``  : 5-sol averages binned by time of day\n"
        f"  - ``average``: 5-sol averages\n"
        f"Works on 'fort.11' files only.\n"
        f"{Green}Example:\n"
        f"> MarsFiles fort.11_* -bin fixed daily diurn average\n"
        f"> MarsFiles fort.11_* -bin fixed diurn"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-c', '--concatenate', action='store_true',
    help=(
        f"Combine sequential files of the same type into one file.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> ls\n"
        f"00334.atmos_average.nc 00668.atmos_average.nc\n"
        f"> MarsFiles *.atmos_average.nc -c\n"
        f"{Blue}Overwrites 00334.atmos_average.nc with concatenated "
        f"files:{Green}\n"
        f"> ls\n"
        f"   00334.atmos_average.nc\n"
        f"{Yellow}To preserve original files, use [-ext --extension]:"
        f"{Green}\n"
        f"> MarsFiles *.atmos_average.nc -c -ext _concatenated\n"
        f"{Blue}Produces 00334.atmos_average_concatenated.nc and "
        f"preserves all other files:{Green}\n"
        f"> ls\n"
        f"   00334.atmos_average.nc 00668.atmos_average.nc "
        f"00334.atmos_average_concatenated.nc\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-split', '--split', nargs='+', type=float,
    help=(
        f"Extract one value or a range of values along a dimension.\n"
        f"Defaults to ``areo`` (Ls) unless otherwise specified using "
        f"-dim.\nIf a file contains multiple Mars Years of data, the "
        f"function splits the file in the first Mars Year.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_average.nc --split 0 90\n"
        f"> MarsFiles 01336.atmos_average.nc --split 270\n"
        f"{Yellow}Use -dim to specify the dimension:{Green}\n"
        f"> MarsFiles 01336.atmos_average.nc --split 0 90 -dim lat"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-t', '--time_shift', action=ExtAction, type=str,
    ext_content='_T',
    parser=parser,
    nargs="?", const=999,
    help=(
        f"Convert hourly binned ('diurn') files to universal local "
        f"time.\nUseful for comparing data at a specific time of day "
        f"across all longitudes.\nWorks on vertically interpolated "
        f"files.\n"
        f"Works on 'diurn' files only.\n"
        f"{Yellow}Generates a new file ending in ``_T.nc``{Nclr}\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_diurn.nc -t"
        f" {Blue}(outputs data for all 24 local times){Green}\n"
        f"> MarsFiles 01336.atmos_diurn.nc -t '3 15'"
        f" {Blue}(outputs data for target local times only)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ba', '--bin_average', action=ExtAction, type=int,
    ext_content='_to_average',
    parser=parser,
    nargs='?', const=5,
    help=(
        f"Calculate 5-day averages from instantaneous data files\n"
        f"(i.e., convert 'daily' files into 'average' files.\n"
        f"Requires input file to be in MGCM 'daily' format.\n"
        f"{Yellow}Generates a new file ending in ``_to_average.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -ba {Blue}(5-sol bin){Green}\n"
        f"> MarsFiles 01336.atmos_daily.nc -ba 10 {Blue}(10-sol bin)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-bd', '--bin_diurn', action=ExtAction, type=int,
    ext_content='_to_diurn',
    parser=parser,
    nargs=0,
    help=(
        f"Calculate 5-day averages binned by hour from instantaneous "
        f"data files\n(i.e., convert 'daily' files into 'diurn' "
        f"files.\nMay be used jointly with --bin_average.\n"
        f"Requires input file to be in MGCM 'daily' format.\n"
        f"{Yellow}Generates a new file ending in ``_to_diurn.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -bd {Blue}(5-sol bin){Green}\n"
        f"> MarsFiles 01336.atmos_daily.nc -bd -ba 10 "
        f"{Blue}(10-sol bin){Green}\n"
        f"> MarsFiles 01336.atmos_daily.nc -bd -ba 1 "
        f"{Blue}(no binning, mimics raw Legacy output)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-hpt', '--high_pass_temporal', action=ExtAction,
    ext_content='_hpt',
    parser=parser,
    nargs='+', type=float,
    help=(
        f"Temporal high-pass filtering: removes low-frequencies \n"
        f"Only works with 'daily' or 'average' files. Requires a cutoff "
        f"frequency in Sols.\n"
        f"Data is detrended before filtering. Use ``-add_trend`` \n"
        f"to add the original linear trend to the amplitudes \n"
        f"\n{Yellow}Generates a new file ending in ``_hpt.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -hpt 10.\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-lpt', '--low_pass_temporal', action=ExtAction,
    ext_content='_lpt',
    parser=parser,
    nargs='+', type=float,
    help=(
        f"Temporal low-pass filtering: removes high-frequencies \n"
        f"Only works with 'daily' or 'average' files. Requires a cutoff "
        f"frequency in Sols.\n"
        f"Data is detrended before filtering. Use ``-add_trend`` \n"
        f"to add the original linear trend to the amplitudes \n"
        f"\n{Yellow}Generates a new file ending in ``_lpt.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -lpt 0.75\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-bpt', '--band_pass_temporal', action=ExtAction,
    ext_content='_bpt',
    parser=parser,
    nargs="+", type=float,
    help=(
        f"Temporal band-pass filtering"
        f"specified by user.\nOnly works with 'daily' or 'average' "
        f"files. Requires two cutoff frequencies in Sols.\n"
        f"Data is detrended before filtering. Use ``-add_trend`` \n"
        f"to add the original linear trend to the amplitudes \n"
        f"\n{Yellow}Generates a new file ending in ``bpt.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -bpt 0.75 10.\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-tide', '--tide_decomp', action=ExtAction,
    ext_content='_tide_decomp',
    parser=parser,
    nargs='+', type=int,
    help=(
        f"Use fourier decomposition to break down the signal into N "
        f"harmonics.\nOnly works with 'diurn' files.\nReturns the phase "
        f"and amplitude of the variable.\n"
        f"N = 1 diurnal tide, N = 2 semi-diurnal, etc.\n"
        f"Works on 'diurn' files only.\n"
        f"{Yellow}Generates a new file ending in ``_tide_decomp.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_diurn.nc -tide 2 -incl ps temp\n"
        f"{Blue}(extracts semi-diurnal tide component of ps and\ntemp "
        f"variables; 2 harmonics)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-regrid', '--regrid_XY_to_match', action=ExtAction,
    ext_content='_regrid',
    parser=parser,
    nargs='+',
    help=(
        f"Regrid the X and Y dimensions of a target file to match a "
        f"user-provided source file.\nBoth files must have the same "
        f"vertical dimensions (i.e., should be vertically\ninterpolated "
        f"to the same standard grid [zstd, zagl, pstd, etc.].\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Yellow}Generates a new file ending in ``_regrid.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_average_pstd.nc -regrid "
        f"01336.atmos_average_pstd_c48.nc\n"
        f"{Yellow}NOTE: regridded file name does not matter:\n"
        f"{Green}> MarsFiles sim1/01336.atmos_average_pstd.nc -regrid "
        f"sim2/01336.atmos_average_pstd.nc"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-zavg', '--zonal_average', action=ExtAction,
    ext_content='_zavg',
    parser=parser,
    nargs=0,
    help=(
        f"Zonally average the entire file over the longitudinal "
        f"dimension.\nDoes not work if the longitude dimension has "
        f"length <= 1.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Yellow}Generates a new file ending in ``_zavg.nc``\n"
        f"{Green}Example:\n"
        "> MarsFiles 01336.atmos_diurn.nc -zavg"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-dim', '--dim_select', type=str, default=None,
    help=(
        f"Must be used with [-split --split]. Flag indicates the "
        f"dimension on which to trim the file.\nAcceptable values are "
        f"'areo', 'lev', 'lat', 'lon'. Default = 'areo'.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_average.nc --split 0 90 -dim areo\n"
        f"> MarsFiles 01336.atmos_average.nc --split -70 -dim lat"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-norm', '--normalize', action=ExtAction,
    ext_content='_norm',
    parser=parser,
    nargs=0,
    help=(
        f"For use with ``-tide``: Returns the result in percent "
        f"amplitude.\n"
        f"N = 1 diurnal tide, N = 2 semi-diurnal, etc.\n"
        f"Works on 'diurn' files only.\n"
        f"{Yellow}Generates a new file ending in ``_norm.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_diurn.nc -tide 6 -incl ps "
        f"-norm"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-add_trend', '--add_trend', action=ExtAction,
    ext_content='_trended',
    parser=parser,
    nargs=0,
    help=(
        f"Return filtered oscillation amplitudes with the linear trend "
        f"added. Works with 'daily' and 'average' files.\n"
        f"For use with temporal filtering utilities (``-lpt``, "
        f"``-hpt``, ``-bpt``).\n"
        f"{Yellow}Generates a new file ending in ``_trended.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -hpt 10. -add_trend\n"
        f"> MarsFiles 01336.atmos_daily.nc -lpt 0.75 -add_trend\n"
        f"> MarsFiles 01336.atmos_daily.nc -bpt 0.75 10. -add_trend"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-recon', '--reconstruct', action=ExtAction,
    ext_content='_reconstruct',
    parser=parser,
    nargs=0,
    help=(
        f"Reconstructs the first N harmonics.\n"
        f"Works on 'diurn' files only.\n"
        f"{Yellow}Generates a new file ending in ``_reconstruct.nc``\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_diurn.nc -tide 6 -incl ps temp "
        f"-recon"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-incl', '--include', nargs='+',
    help=(
        f"Flag to use only the variables specified in a calculation.\n"
        f"All dimensional and 1D variables are ported to the new file "
        f"automatically.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_daily.nc -ba -incl ps temp ts"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ext', '--extension', type=str, default=None,
    help=(
        f"Must be paired with an argument listed above.\nInstead of "
        f"overwriting a file to perform a function, ``-ext``\ntells "
        f"CAP to create a new file with the extension name specified "
        f"here.\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_average.nc -c -ext _my_concat\n"
        f"{Blue}(Produces 01336.atmos_average_my_concat.nc and "
        f"preserves all other files)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsFiles 01336.atmos_average.nc -t '3 15' --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()
debug = args.debug

if args.input_file:
    for file in args.input_file:
        if not (re.search(".nc", file.name) or 
                re.search("fort.11", file.name)):
            parser.error(
                f"{Red}{file.name} is not a netCDF or fort.11 file{Nclr}"
                )
            exit()

if args.dim_select and not args.split:
    parser.error(f"{Red}[-dim --dim_select] must be used with [-split "
                 f"--split]{Nclr}")
    exit()

if args.normalize and not args.tide_decomp:
    parser.error(f"{Red}[-norm --normalize] must be used with [-tide "
                 f"--tide_decomp]{Nclr}")
    exit()

if args.add_trend and not (args.high_pass_temporal or args.low_pass_temporal
                           or args.band_pass_temporal):
    parser.error(f"{Red}[-add_trend --add_trend] must be used with [-lpt "
                 f"--low_pass_temporal], [-hpt --high_pass_temporal], or "
                 f"[-bpt --band_pass_temporal]{Nclr}")
    exit()

if args.reconstruct and not args.tide_decomp:
    parser.error(f"{Red}[-recon --reconstruct] must be used with [-tide "
                 f"--tide_decomp]{Nclr}")
    exit()

all_args = [args.bin_files, args.concatenate, args.split, args.time_shift,
            args.bin_average, args.bin_diurn, args.high_pass_temporal,
            args.low_pass_temporal, args.band_pass_temporal,
            args.tide_decomp, args.normalize,
            args.regrid_XY_to_match, args.zonal_average]

if (all(v is None or v is False for v in all_args) 
    and args.include is not None):
    parser.error(f"{Red}[-incl --include] must be used with another "
                 f"argument{Nclr}")
    exit()

if (all(v is None or v is False for v in all_args) and 
    args.extension is not None):
    parser.error(f"{Red}[-ext --extension] must be used with another "
                 f"argument{Nclr}")
    exit()

# ------------------------------------------------------
#                  EXTENSION FUNCTION
# ------------------------------------------------------
# Concatenates extensions to append to the file name depending on
# user-provided arguments.
"""
This is the documentation for out_ext
:no-index:
"""
out_ext = (f"{args.time_shift_ext}"
            f"{args.bin_average_ext}"
            f"{args.bin_diurn_ext}"
            f"{args.high_pass_temporal_ext}"
            f"{args.low_pass_temporal_ext}"
            f"{args.band_pass_temporal_ext}"
            f"{args.add_trend_ext}"
            f"{args.tide_decomp_ext}"
            f"{args.reconstruct_ext}"
            f"{args.normalize_ext}"
            f"{args.regrid_XY_to_match_ext}"
            f"{args.zonal_average_ext}"
            )

if args.extension:
    # Append extension, if any:
    out_ext = (f"{out_ext}_{args.extension}")

# ------------------------------------------------------
#                  DEFINITIONS
# ------------------------------------------------------
def concatenate_files(file_list, full_file_list):
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
        for i in range(1, num_files):
            try:
                os.remove(full_file_list[i])
            except OSError as e:
                print(f"{Yellow}Warning: Could not remove {full_file_list[i]}: {e}{Nclr}")
        print(f"{Cyan}Cleaned all but {file_list[0]}{Nclr}")
        exit()

    print(
        f"{Cyan}Merging {num_files} files starting with {file_list[0]}..."
          f"{Nclr}"
          )

    if args.include:
        # Exclude variables NOT listed after --include
        f = Dataset(file_list[0], "r")
        exclude_list = filter_vars(f, args.include, giveExclude = True)
        f.close()
    else:
        exclude_list = []

    # Create a temporary file ending in _tmp.nc to work in
    file_base = os.path.splitext(full_file_list[0])[0]
    tmp_file = f"{file_base}_tmp.nc"
    Log = Ncdf(tmp_file, "Merged file")
    Log.merge_files_from_list(full_file_list, exclude_var=exclude_list)
    Log.close()

    # ----- Delete the files that were used for concatenate -----

    # First, rename temporary file for the final merged file
    #   For Legacy netCDF files, rename using initial and end Ls
    #   For MGCM netCDF files, rename to the first file in the list
    filename_base = os.path.basename(file_list[0])
    if filename_base[:12] == "LegacyGCM_Ls":
        ls_ini = file_list[0][12:15]
        ls_end = file_list[-1][18:21]
        merged_file = f"LegacyGCM_Ls{ls_ini}_Ls{ls_end}.nc"
    else:
        merged_file = full_file_list[0]

    # Delete the files that were concatenated.
    # Apply the new name created above
    for file in full_file_list:
        try:
            os.remove(file)
        except OSError as e:
            print(f"{Yellow}Warning: Could not remove {file}: {e}{Nclr}")
            
    try:
        os.rename(tmp_file, merged_file)
    except OSError as e:
        print(f"{Red}Error: Could not rename {tmp_file} to {merged_file}: {e}{Nclr}")
        # Fallback option using shutil for cross-device moves
        import shutil
        shutil.move(tmp_file, merged_file)
        
    print(f"{Cyan}{merged_file} was created from a merge{Nclr}")

    return

def split_files(file_list, split_dim):
    """
    Extracts variables in the file along the time dimension, unless
    other dimension is specified (lev, lat, or lon).

    :param file_list: list of file names
    :type split_dim: dimension along which to perform extraction

    :returns: new file with sliced dimensions

    """
    if split_dim not in ['time','areo', 'lev', 'lat', 'lon']:
        print(f"{Red}Split dimension must be one of the following:"
              f"    time, areo, lev, lat, lon{Nclr}")

    bounds = np.asarray(args.split).astype(float)

    if len(np.atleast_1d(bounds)) > 2 or len(np.atleast_1d(bounds)) < 1:
        print(
            f"{Red}Accepts only ONE or TWO values: [bound] to reduce one "
            f"dimension to a single value, [lower_bound] [upper_bound] to "
            f"reduce one dimension to a range{Nclr}"
            )
        exit()

    # Use os.path.join for platform-independent path handling
    for file in file_list:
        # Get the file name from the file object
        if hasattr(file, 'name'):
            input_file_name = file.name  # Get the name attribute of the file object
        else:
            input_file_name = file  # In case it's already a string
            
    print(f"{Yellow}split input_file_name = {input_file_name}{Nclr}\n\n")

    file_base = os.path.splitext(input_file_name)[0]
    output_file_name = f"{file_base}{out_ext}.nc"
    print(f"{Cyan}split output_file_name = {output_file_name}{Nclr}")
            
    original_date = os.path.basename(os.path.splitext(input_file_name)[0])

    fNcdf = Dataset(str(input_file_name), 'r', format='NETCDF4_CLASSIC')
    var_list = filter_vars(fNcdf, args.include)

    # Get file type (diurn, average, daily, etc.)
    f_type, interp_type = FV3_file_type(fNcdf)

    if split_dim == 'lev':
        split_dim = interp_type
        if interp_type == 'pstd':
            unt_txt = 'Pa'
            lnm_txt = 'standard pressure'
        elif interp_type == 'zagl':
            unt_txt = 'm'
            lnm_txt = 'altitude above ground level'
        elif interp_type == 'zstd':
            unt_txt = 'm'
            lnm_txt = 'standard altitude'
        else:
            unt_txt = 'mb'
            lnm_txt = 'ref full pressure level'

    # Remove all single dimensions from areo (scalar_axis)
    if f_type == 'diurn':
        if split_dim == 'areo':
            # size areo = (time, tod, scalar_axis)
            reducing_dim = np.squeeze(fNcdf.variables[split_dim][:, 0, :])
        else:
            reducing_dim = np.squeeze(fNcdf.variables[split_dim][:])
    else:
        if split_dim == 'areo':
            # size areo = (time, scalar_axis)
            reducing_dim = np.squeeze(fNcdf.variables[split_dim][:])
        else:
            reducing_dim = np.squeeze(fNcdf.variables[split_dim][:])

    if split_dim == 'areo':
        print(f"\n{Yellow}Areo MOD(360):\n{reducing_dim%360.}\n")
    
    print(f"\n{Yellow}All values in dimension:\n{reducing_dim}\n")

    dx = np.abs(reducing_dim[1] - reducing_dim[0])

    bounds_in = bounds.copy()
    if split_dim == 'areo':
        while (bounds[0] < (reducing_dim[0] - dx)):
            bounds += 360.
        if len(np.atleast_1d(bounds)) == 2:
            while (bounds[1] < bounds[0]):
                bounds[1] += 360.

    if len(np.atleast_1d(bounds)) < 2:
        check_bounds(bounds[0], reducing_dim[0], reducing_dim[-1], dx)
        indices = [(np.abs(reducing_dim - bounds[0])).argmin()]
        dim_out = reducing_dim[indices]
        print(f"Requested value = {bounds[0]}\n"
              f"Nearest value = {dim_out[0]}\n")
    else:
        check_bounds(bounds, reducing_dim[0], reducing_dim[-1], dx)
        if ((split_dim == 'lon') & (bounds[1] < bounds[0])):
            indices = np.where(
                (reducing_dim <= bounds[1]) | 
                (reducing_dim >= bounds[0])
                )[0]
        else:
            indices = np.where(
                (reducing_dim >= bounds[0]) & 
                (reducing_dim <= bounds[1])
                )[0]
        dim_out = reducing_dim[indices]
        print(
            f"Requested range = {bounds_in[0]} - {bounds_in[1]}\n"
            f"Corresponding values = {dim_out}\n"
            )
        
        if len(indices) == 0:
            print(
                f"{Red}Warning, no values were found in the range {split_dim} "
                f"{bounds_in[0]}, {bounds_in[1]}) ({split_dim} values range "
                f"from {reducing_dim[0]:.1f} to {reducing_dim[-1]:.1f})"
                )
            exit()

    if split_dim in ('time', 'areo'):
        time_dim = (np.squeeze(fNcdf.variables['time'][:]))[indices]
        print(f"time_dim = {time_dim}")
        
        # If time_dim is an array with multiple values, use the first one for naming
        if hasattr(time_dim, '__len__') and len(time_dim) > 1:
            time_value = time_dim[0]  # Use the first time value for the filename
        else:
            time_value = time_dim  # Use as is if it's a scalar

    # For maintaining the 00000.atmos_average.nc format:
    base_path = os.path.dirname(input_file_name)
    file_name = os.path.basename(input_file_name)
    # Split at the first dot to preserve the numeric prefix
    file_date, file_base = file_name.split('.', 1)
    if split_dim == 'time':
        if len(np.atleast_1d(bounds)) < 2:
            out_name = (f"{int(time_value):05d}.{file_base}_nearest_sol{int(bounds_in[0]):03d}.nc")
        else:
            out_name = (f"{int(time_value):05d}.{file_base}_sol{int(bounds_in[0]):05d}_{int(bounds_in[1]):05d}.nc")
    elif split_dim =='areo':
        if len(np.atleast_1d(bounds)) < 2:
            out_name = (f"{int(time_value):05d}.{file_base}_nearest_Ls{int(bounds_in[0]):03d}.nc")
        else:
            out_name = (f"{int(time_value):05d}.{file_base}_Ls{int(bounds_in[0]):03d}_{int(bounds_in[1]):03d}.nc")
        split_dim = 'time'
    elif split_dim == 'lat':
        new_bounds = [f"{abs(int(b))}S" if b < 0 
                      else f"{int(b)}N" for b in bounds]
        if len(np.atleast_1d(bounds)) < 2:
            out_name = (f"{original_date}.{file_base}_nearest_{split_dim}_{new_bounds[0]}.nc")
        else:
            print(f"{Yellow}bounds = {bounds[0]} {bounds[1]}")
            print(f"{Yellow}new_bounds = {new_bounds[0]} {new_bounds[1]}")
            out_name = (f"{original_date}.{file_base}_{split_dim}_{new_bounds[0]}_{new_bounds[1]}.nc")
    elif split_dim == interp_type:
        if interp_type == 'pfull':
            new_bounds = [f"{abs(int(b*100))}Pa" if b < 1 
                          else f"{int(b)}{unt_txt} " for b in bounds]
        elif interp_type == 'pstd':
            new_bounds = [f"{abs(int(b*100))}hPa" if b < 1 
                          else f"{int(b)}{unt_txt}" for b in bounds]
        else:
            new_bounds = [f"{int(b)}{unt_txt}" for b in bounds]
            
        if len(np.atleast_1d(bounds)) < 2:
            print(f"{Yellow}bounds = {bounds[0]}")
            print(f"{Yellow}new_bounds = {new_bounds[0]}")
            out_name = (f"{original_date}.{file_base}_{split_dim}_nearest_{new_bounds[0]}.nc")
        else:
            print(f"{Yellow}bounds = {bounds[0]} {bounds[1]}")
            print(f"{Yellow}new_bounds = {new_bounds[0]} {new_bounds[1]}")
            out_name = (f"{original_date}.{file_base}_{split_dim}_{new_bounds[0]}_{new_bounds[1]}.nc")
    else:
        if len(np.atleast_1d(bounds)) < 2:
            out_name = (f"{original_date}.{file_base}_nearest_{split_dim}_{int(bounds[0]):03d}.nc")
        else:
            out_name = (f"{original_date}.{file_base}_{split_dim}_{int(bounds[0]):03d}_{int(bounds[1]):03d}.nc")

    # Append extension, if any:
    output_file_name = os.path.join(base_path, out_name)
    print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

    Log = Ncdf(output_file_name)

    Log.copy_all_dims_from_Ncfile(fNcdf, exclude_dim=[split_dim])

    if split_dim == 'time':
        Log.add_dimension(split_dim, None)
    else:
        Log.add_dimension(split_dim, len(dim_out))

    if split_dim == 'time':
        Log.log_axis1D(
            variable_name = 'time',
            DATAin = dim_out,
            dim_name = 'time',
            longname_txt = 'sol number',
            units_txt = 'days since 0000-00-00 00:00:00',
            cart_txt = 'T'
            )
    elif split_dim == 'lat':
        Log.log_axis1D(
            variable_name = 'lat',
            DATAin = dim_out,
            dim_name = 'lat',
            longname_txt = 'latitude',
            units_txt = 'degrees_N',
            cart_txt = 'T'
            )
    elif split_dim == 'lon':
        Log.log_axis1D(
            variable_name = 'lon',
            DATAin = dim_out,
            dim_name = 'lon',
            longname_txt = 'longitude',
            units_txt = 'degrees_E',
            cart_txt = 'T'
            )
    elif split_dim == interp_type:
        Log.log_axis1D(
            variable_name = split_dim,
            DATAin = dim_out,
            dim_name = split_dim,
            longname_txt = lnm_txt,
            units_txt = unt_txt,
            cart_txt = 'T'
            )

    # Loop over all variables in the file
    for ivar in var_list:
        varNcf = fNcdf.variables[ivar]
        if split_dim in varNcf.dimensions and ivar != split_dim:
            # ivar is a dim of ivar but ivar is not ivar
            print(f'{Cyan}Processing: {ivar}...{Nclr}')
            if split_dim == 'time':
                var_out = varNcf[indices, ...]
            elif split_dim == 'lat' and varNcf.ndim == 5:
                var_out = varNcf[:, :, :, indices, :]
            elif split_dim == 'lat' and varNcf.ndim == 4:
                var_out = varNcf[:, :, indices, :]
            elif split_dim == 'lat' and varNcf.ndim == 3:
                var_out = varNcf[:, indices, :]
            elif split_dim == 'lat' and varNcf.ndim == 2:
                var_out = varNcf[indices, ...]
            elif split_dim == 'lon' and varNcf.ndim > 2:
                var_out = varNcf[..., indices]
            elif split_dim == 'lon' and varNcf.ndim == 2:
                var_out = varNcf[indices, ...]
            elif split_dim == interp_type and varNcf.ndim == 5:
                var_out = varNcf[:, :, indices, :, :]
            elif split_dim == interp_type and varNcf.ndim == 4:
                var_out = varNcf[:, indices, :, :]
            longname_txt, units_txt = get_longname_unit(fNcdf, ivar)
            Log.log_variable(
                ivar, 
                var_out, 
                varNcf.dimensions,
                longname_txt, 
                units_txt
                )
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

# ------------------------------------------------------
#            Time-Shifting Implementation
#                   Victoria H.
# ------------------------------------------------------
def process_time_shift(file_list):
    """
    This function converts the data in diurn files with a time_of_day_XX
    dimension to universal local time.

    :param file_list: list of file names
    :type file_list: list

    """
    if args.time_shift == 999:
        # Target local times requested by user
        target_list = None
    else:
        target_list = np.fromstring(args.time_shift, dtype=float, sep=" ")

    for file in file_list:
        # Get the file name from the file object
        if hasattr(file, 'name'):
            input_file_name = file.name  # Get the name attribute of the file object
        else:
            input_file_name = file  # In case it's already a string
        
        file_base = os.path.splitext(input_file_name)[0]
        output_file_name = f"{file_base}{out_ext}.nc"
        print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")
        
        file_name = os.path.splitext(input_file_name)[0]
        base_path = os.path.dirname(input_file_name)
        output_file_name2 = os.path.join(base_path, f"{file_name}{out_ext}")
        print(f"{Cyan}output_file_name2 = {output_file_name2}{Nclr}")

        fdiurn = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")
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
        try:
            fnew.copy_Ncaxis_with_content(fdiurn.variables["scalar_axis"])
        except:
            print(f'{Red}Could not find scalar axis')

        # Only create a vertical axis if orig. file has 3D fields
        if zaxis in fdiurn.dimensions.keys():
            fnew.copy_Ncaxis_with_content(fdiurn.variables[zaxis])

        # Take care of TOD dimension in new file
        tod_orig = np.array(fdiurn.variables[tod_name_in][:])

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
            fnew.add_dim_with_content(
                tod_name_out, 
                target_list,
                longname_txt = "time of day",
                units_txt = ("[hours since 0000-00-00 00:00:00]"),
                cart_txt = ""
                )
            
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

            fnew.add_dim_with_content(
                dimension_name = "scalar_axis", 
                DATAin = [0], 
                longname_txt = "none",
                units_txt = "none"
                )
            fnew.log_variable("areo", areo_out, areo_dims, "areo", "degrees")

        # Read in 4D field(s) and do the time shift. Exclude vars
        # not listed after --include in var_list
        lons = np.array(fdiurn.variables["lon"][:])
        var_list = filter_vars(fdiurn, args.include)

        for var in var_list:
            print(f"{Cyan}Processing: {var}...{Nclr}")
            value = fdiurn.variables[var][:]
            dims = fdiurn.variables[var].dimensions
            longname_txt, units_txt = get_longname_unit(fdiurn, var)

            if (len(dims) >= 4):
                y = dims.index("lat")
                x = dims.index("lon")
                t = dims.index("time")
                tod = dims.index(tod_name_in)

            if (len(dims) == 4):
                # time, tod, lat, lon
                var_val_tmp = np.transpose(value, (x, y, t, tod))
                var_val_T = time_shift_calc(
                    var_in = var_val_tmp, 
                    lon = lons, 
                    tod = tod_orig,
                    target_times = target_list
                    )
                var_out = np.transpose(var_val_T, (2, 3, 1, 0))
                fnew.log_variable(
                    var, 
                    var_out,
                    ["time", tod_name_out, "lat", "lon"],
                    longname_txt, 
                    units_txt
                    )
            if (len(dims) == 5):
                # time, tod, Z, lat, lon
                z = dims.index(zaxis)
                var_val_tmp = np.transpose(value, (x, y, z, t, tod))
                var_val_T = time_shift_calc(
                    var_in = var_val_tmp, 
                    lon = lons, 
                    tod = tod_orig,
                    target_times = target_list
                    )
                var_out = np.transpose(var_val_T, (3, 4, 2, 1, 0))
                fnew.log_variable(
                    var, 
                    var_out,
                    ["time", tod_name_out, zaxis, "lat", "lon"],
                    longname_txt,
                    units_txt
                    )
        fnew.close()
        fdiurn.close()
    return

# ------------------------------------------------------
#                  MAIN PROGRAM
# ------------------------------------------------------
@debug_wrapper
def main():
    global data_dir
    file_list = [f.name for f in args.input_file]
    data_dir = os.getcwd()

    # Make a list of input files including the full path to the dir
    full_file_list = []
    for file in file_list:
        if not ("/" in file):
            full_file_list.append(f"{data_dir}/{file}")
        else:
            full_file_list.append(f"{file}")

    if args.bin_files and args.concatenate:
        print(
            f"{Red}Use --bin_files and --concatenate separately to avoid "
            f"ambiguity"
            )
        exit()

    # ------------------------------------------------------------------
    #               Conversion Legacy -> MGCM Format
    #                    Richard U. and Alex. K.
    # ------------------------------------------------------------------
    # Convert to MGCM Output Format
    if args.bin_files:
        for req_file in args.bin_files:
            if req_file not in ["fixed", "average", "daily", "diurn"]:
                print(
                    f"{Red}{req_file} is invalid. Select 'fixed', 'average', "
                    f"'daily', or 'diurn'{Nclr}"
                    )

        if full_file_list[0][-3:] == ".nc":
            print("Processing Legacy MGCM netCDF files")
            for f in full_file_list:
                make_FV3_files(f, args.bin_files, True)
        elif "fort.11" in full_file_list[0]:
            print("Processing fort.11 files")
            for f in full_file_list:
                file_name = Fort(f)
                if "fixed" in args.bin_files:
                    file_name.write_to_fixed()
                if "average" in args.bin_files:
                    file_name.write_to_average()
                if "daily" in args.bin_files:
                    file_name.write_to_daily()
                if "diurn" in args.bin_files:
                    file_name.write_to_diurn()
        else:
            print(f"Cannot bin {full_file_list[0]} like MGCM file.")

    elif args.concatenate:
        # Combine files along the time dimension
        concatenate_files(file_list, full_file_list)

    elif args.split:
        # Split file along the specified dimension. If none specified,
        # default to time dimension
        split_dim = 'areo' if args.dim_select == None else args.dim_select
        split_files(file_list, split_dim)

    elif args.time_shift:
        # Time-shift files
        process_time_shift(file_list)

    # ------------------------------------------------------------------
    #               Bin a daily file as an average file
    #                               Alex K.
    # ------------------------------------------------------------------
    elif (args.bin_average and 
          not args.bin_diurn):

        # Generate output file name
        bin_period = args.bin_average
        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            fdaily = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, args.include)

            time = fdaily.variables["time"][:]

            time_increment = time[1] - time[0]
            dt_per_day = int(np.round(1 / time_increment))
            dt_total = int(dt_per_day * bin_period)

            bins = len(time) / (dt_per_day*bin_period)
            bins_even = len(time) // dt_total
            bins_left = len(time) % dt_total

            if bins_left != 0:
                print(
                    f"{Yellow}*** Warning *** Requested a {bin_period}-sol "
                    f"bin period but the file has a total of {len(time)}"
                    f"timesteps ({dt_per_day} per sol) and {len(time)}/"
                    f"({bin_period}x{dt_per_day})={bins} is not a round "
                    f"number.\nWill use {bins_even} bins with {bin_period}x"
                    f"{dt_per_day}={dt_total} timesteps per bin "
                    f"({bins_even*dt_total} timsteps total) and discard "
                    f"{bins_left} timesteps.{Nclr}"
                    )

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["time"])

            # Calculate and log the new time array
            fnew.add_dimension("time", None)
            time_out = daily_to_average(time[:], time_increment, bin_period)
            fnew.log_axis1D(
                "time", 
                time_out, 
                "time",
                longname_txt = "sol number",
                units_txt = "days since 0000-00-00 00:00:00",
                cart_txt = "T"
                )

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]

                if "time" in varNcf.dimensions:
                    print(f"{Cyan}Processing: {ivar}{Nclr}")
                    var_out = daily_to_average(
                        varNcf[:], 
                        time_increment, 
                        bin_period
                        )
                    longname_txt, units_txt = get_longname_unit(fdaily, ivar)
                    fnew.log_variable(
                        ivar, 
                        var_out, 
                        varNcf.dimensions,
                        longname_txt, 
                        units_txt
                        )
                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ------------------------------------------------------------------
    #               Bin a daily file as a diurn file
    #                               Alex K.
    # ------------------------------------------------------------------
    elif args.bin_diurn:
        # Use defaut binning period of 5 days (like average files)
        if args.bin_average is None:
            bin_period = 5
        else:
            bin_period = args.bin_average

        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            fdaily = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, args.include)

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
            fnew.add_dim_with_content(
                "time",
                time_out,
                longname_txt = "sol number",
                units_txt = ("days since 0000-00-00 00:00:00"),
                cart_txt = "T"
                )

            # Create a new time_of_day dimension
            tod_name = f"time_of_day_{dt_per_day:02}"
            time_tod = np.squeeze(daily_to_diurn(time[0:dt_per_day],
                                                 time[0:dt_per_day]))
            tod = np.mod(time_tod*24, 24)
            fnew.add_dim_with_content(
                tod_name,
                tod,
                longname_txt = "time of day",
                units_txt = ("hours since 0000-00-00 00:00:00"),
                cart_txt = "N"
                )

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
                    longname_txt, units_txt = get_longname_unit(fdaily, ivar)
                    fnew.log_variable(
                        ivar, 
                        var_out, 
                        dims_out,
                        longname_txt,
                        units_txt
                        )
                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    elif ivar != "time":
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ------------------------------------------------------------------
    #                       Temporal Filtering
    #                       Alex K. & R. J. Wilson
    # ------------------------------------------------------------------
    elif (args.high_pass_temporal or
          args.low_pass_temporal or
          args.band_pass_temporal):
        # This function requires scipy > 1.2.0. Import package here.
        from amescap.Spectral_utils import zeroPhi_filter

        if args.high_pass_temporal:
            btype = "high"
            nsol = np.asarray(args.high_pass_temporal).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                print(f"{Red}***Error*** sol_min accepts only one value")
                exit()
        if args.low_pass_temporal:
            btype = "low"
            nsol = np.asarray(args.low_pass_temporal).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                print(f"{Red}sol_max accepts only one value")
                exit()
        if args.band_pass_temporal:
            btype = "band"
            nsol = np.asarray(args.band_pass_temporal).astype(float)
            if len(np.atleast_1d(nsol)) != 2:
                print(f"{Red}Requires two values: sol_min sol_max")
                exit()

        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            fdaily = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, args.include)
            time = fdaily.variables["time"][:]
            dt = time[1] - time[0]

            # Check if the frequency domain is allowed
            if any(nn <= 2*dt for nn in nsol):
                print(
                    f"{Red}***Error***  minimum cut-off cannot be smaller "
                    f"than the Nyquist period of 2xdt={2*dt} sol{Nclr}"
                    )
                exit()

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the
            # new file
            fnew.copy_all_dims_from_Ncfile(fdaily)

            if btype == "low":
                fnew.add_constant(
                    "sol_max", 
                    nsol,
                    "Low-pass filter cut-off period ",
                    "sol"
                    )
            elif btype == "high":
                fnew.add_constant(
                    "sol_min", 
                    nsol,
                    "High-pass filter cut-off period ",
                    "sol"
                    )
            elif btype == "band":
                fnew.add_constant(
                    "sol_min",
                    nsol[0],
                    "High-pass filter low cut-off period ",
                    "sol"
                    )
                fnew.add_constant(
                    "sol_max", 
                    nsol[1],
                    "High-pass filter high cut-off period ",
                    "sol"
                    )
            dt = time[1] - time[0]
            if dt == 0:
                print(f"{Red}***Error*** dt = 0, time dimension is not increasing")
                exit()
                
            fs = 1/(dt) # frequency in sol-1
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
                        varNcf[:], 
                        btype, 
                        low_highcut, 
                        fs, 
                        axis = 0, 
                        order = 4,
                        add_trend = args.add_trend
                        )
                    longname_txt, units_txt = get_longname_unit(fdaily, ivar)
                    fnew.log_variable(
                        ivar, 
                        var_out, 
                        varNcf.dimensions,
                        longname_txt, 
                        units_txt
                        )
                else:
                    if ivar in ["pfull", "lat", "lon", "phalf", "pk",
                                "bk", "pstd", "zstd", "zagl"]:
                        print(f"{Cyan}Copying axis: {ivar}{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        print(f"{Cyan}Copying variable: {ivar}{Nclr}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()
            
    # ------------------------------------------------------------------
    #                           Tidal Analysis
    #                           Alex K. & R. J. Wilson
    # ------------------------------------------------------------------
    elif args.tide_decomp:
        from amescap.Spectral_utils import diurn_extract, reconstruct_diurn
        N = args.tide_decomp[0]
        if len(np.atleast_1d(N)) != 1:
            print(f"{Red}***Error*** N accepts only one value")
            exit()

        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            fdiurn = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")

            var_list = filter_vars(fdiurn, args.include)

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
            if args.reconstruct:
                fnew.copy_all_dims_from_Ncfile(fdiurn)
                # Copy time_of_day axis from initial files
                fnew.copy_Ncaxis_with_content(fdiurn.variables[tod_name])

            else:
                fnew.copy_all_dims_from_Ncfile(
                    fdiurn, exclude_dim = [tod_name]
                    )
                # Create new dimension holding the harmonics. We reuse
                # the time_of_day name to facilitate. Compatible with
                # other routines, but keep in mind this is the harmonic
                # number
                fnew.add_dim_with_content(
                    dimension_name = f"time_of_day_{N}",
                    DATAin = np.arange(1, N+1),
                    longname_txt = "tidal harmonics",
                    units_txt = "Diurnal harmonic number",
                    cart_txt = "N"
                )

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdiurn.variables[ivar]
                varIN = varNcf[:]
                longname_txt, units_txt = get_longname_unit(fdiurn, ivar)

                if (tod_name in varNcf.dimensions and
                    ivar not in [tod_name, "areo"] and
                    len(varNcf.shape) > 2):
                    print(f"{Cyan}Processing: {ivar}{Nclr}")

                    # Normalize the data
                    if args.normalize:
                        # Normalize and reshape the array along the
                        # time_of_day dimension
                        norm = np.mean(varIN, axis = 1)[:, np.newaxis, ...]
                        varIN = 100*(varIN-norm)/norm
                        units_txt = f"% of diurnal mean"

                    amp, phas = diurn_extract(
                        varIN.swapaxes(0, 1), 
                        N,
                        target_tod, 
                        lon
                        )
                    if args.reconstruct:
                        VARN = reconstruct_diurn(
                            amp, 
                            phas, 
                            target_tod, 
                            lon,
                            sumList=[]
                            )
                        for nn in range(N):
                            fnew.log_variable(
                                f"{ivar}_N{nn+1}",
                                VARN[nn, ...].swapaxes(0, 1),
                                varNcf.dimensions,
                                (f"harmonic N={nn+1} for {longname_txt}"),
                                units_txt
                                )

                    else:
                        # Update the dimensions
                        new_dim = list(varNcf.dimensions)
                        new_dim[1] = f"time_of_day_{N}"
                        fnew.log_variable(
                            f"{ivar}_amp", 
                            amp.swapaxes(0,1),
                            new_dim,
                            f"tidal amplitude for {longname_txt}",
                            units_txt
                            )
                        fnew.log_variable(
                            f"{ivar}_phas", 
                            phas.swapaxes(0,1),
                            new_dim,
                            f"tidal phase for {longname_txt}",
                            "hr"
                            )

                elif  ivar in ["pfull", "lat", "lon", "phalf", "pk",
                               "bk", "pstd", "zstd", "zagl", "time"]:
                        print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                        fnew.copy_Ncaxis_with_content(fdiurn.variables[ivar])
                elif  ivar in ["areo"]:
                        if args.reconstruct:
                            # time_of_day is the same size as the
                            # original file
                            print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                            fnew.copy_Ncvar(fdiurn.variables["areo"])
                        else:
                            print(f"{Cyan}Processing: {ivar}...{Nclr}")
                            # Create areo variable reflecting the
                            # new shape
                            areo_new = np.zeros((areo.shape[0], N, 1))
                            # Copy areo
                            for xx in range(N):
                                areo_new[:, xx, :] = areo[:, 0, :]
                            # Update the dimensions
                            new_dim = list(varNcf.dimensions)
                            new_dim[1] = f"time_of_day_{N}"
                            # fnew.log_variable(ivar, bareo_new, new_dim,
                            # longname_txt, units_txt)
                            fnew.log_variable(
                                ivar, 
                                areo_new, 
                                new_dim,
                                longname_txt, 
                                units_txt
                                )
            fnew.close()

    # ------------------------------------------------------------------
    #                           Regridding Routine
    #                                 Alex K.
    # ------------------------------------------------------------------
    elif args.regrid_XY_to_match:
        name_target = args.regrid_XY_to_match[0]

        # Add path unless full path is provided
        if not ("/" in name_target):
            name_target = f"{data_dir}/{name_target}"
        fNcdf_t = Dataset(name_target, "r")

        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            f_in = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")

            var_list = filter_vars(f_in, args.include)  # Get all variables

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)

            # Copy all dims from the target file to the new file
            fnew.copy_all_dims_from_Ncfile(fNcdf_t)

            # Loop over all variables in the file
            print(var_list)
            for ivar in var_list:
                varNcf = f_in.variables[ivar]
                longname_txt, units_txt = get_longname_unit(f_in, ivar)

                if  ivar in ["pfull", "lat", "lon", "phalf", "pk",
                             "bk", "pstd", "zstd", "zagl", "time", "areo"]:
                        if ivar in fNcdf_t.variables.keys():
                            print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                            fnew.copy_Ncaxis_with_content(
                                fNcdf_t.variables[ivar]
                                )
                elif varNcf.dimensions[-2:] == ("lat", "lon"):
                    # Ignore variables like time_bounds, scalar_axis
                    # or grid_xt_bnds...
                    print(f"{Cyan}Regridding: {ivar}...{Nclr}")
                    var_OUT = regrid_Ncfile(varNcf, f_in, fNcdf_t)
                    fnew.log_variable(
                        ivar, 
                        var_OUT, 
                        varNcf.dimensions,
                        longname_txt, 
                        units_txt
                        )
            fnew.close()
            fNcdf_t.close()

    # ------------------------------------------------------------------
    #                           Zonal Averaging
    #                              Alex K.
    # ------------------------------------------------------------------
    elif args.zonal_average:
        for file in file_list:
            # Use os.path.join for platform-independent path handling
            print(f"\n\n{Yellow}dir = {os.path.dirname(file)}{Nclr}")
            if os.path.dirname(file) == "":
                input_file_name = file
            else:
                input_file_name = os.path.join(data_dir, file)
            print(f"{Yellow}input_file_name = {input_file_name}{Nclr}\n\n")

            file_base = os.path.splitext(input_file_name)[0]
            output_file_name = f"{file_base}{out_ext}.nc"
            print(f"{Cyan}output_file_name = {output_file_name}{Nclr}")

            fdaily = Dataset(input_file_name, "r", format="NETCDF4_CLASSIC")
            var_list = filter_vars(fdaily, args.include) # Get all variables

            lon_in = fdaily.variables["lon"][:]

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(output_file_name)
            # Copy all dimensions but time from the old file to the
            # new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["lon"])

            # Add a new dimension for the longitude, size = 1
            fnew.add_dim_with_content(
                "lon", 
                [lon_in.mean()],
                longname_txt = "longitude",
                units_txt = "degrees_E", 
                cart_txt = "X"
                )

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]
                longname_txt,units_txt = get_longname_unit(fdaily,ivar)
                if ("lon" in varNcf.dimensions and
                    ivar not in ["lon", "grid_xt_bnds", "grid_yt_bnds"]):
                    print(f"{Cyan}Processing: {ivar}...{Nclr}")
                    with warnings.catch_warnings():
                        warnings.simplefilter(
                            "ignore", category = RuntimeWarning
                            )
                        var_out = np.nanmean(
                            varNcf[:], axis = -1
                            )[..., np.newaxis]
                        fnew.log_variable(
                            ivar, 
                            var_out, 
                            varNcf.dimensions,
                            longname_txt, 
                            units_txt
                            )
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
        print(
            f"{Red}Error: no action requested: use ``MarsFiles *nc "
            f"--bin_files --concatenate, --time_shift, --bin_average, "
            f"--bin_diurn etc ...``{Nclr}")

# ----------------------------------------------------------------------
#                  DEFINITIONS
# ----------------------------------------------------------------------
def make_FV3_files(fpath, typelistfv3, renameFV3=True):
    """
    Make MGCM-like ``average``, ``daily``, and ``diurn`` files.
    Used if call to [``-bin --bin_files``] is made AND Legacy files are in
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
    histfile = Dataset(fpath, "r", format="NETCDF4_CLASSIC")
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
                newf.add_dim_with_content(
                    "lon", 
                    npvar, 
                    "longitude",
                    getattr(var, "units"), 
                    "X"
                    )
            elif dname == "nlat":
                var = histfile.variables["latitude"]
                npvar = var[:]
                newf.add_dim_with_content(
                    "lat", 
                    npvar, 
                    "latitude",
                    getattr(var, "units"), 
                    "Y"
                    )

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
                # ..note:: pk in amesCAP/mars_data/Legacy.fixed.nc
                # is also updated
                for z in range(num):
                    bk[z + 1] = sgm[2*z + 2]
                phalf[:] = pk[:] + pref*bk[:] # [Pa]

                # First layer:
                if pk[0] == 0 and bk[0] == 0:
                    pfull[0] = 0.5*(phalf[0] + phalf[1])
                else:
                    pfull[0] = (
                        (phalf[1]-phalf[0]) 
                        / (np.log(phalf[1]) - np.log(phalf[0]))
                        )
                # All other layers:
                pfull[1:] = (
                    (phalf[2:]-phalf[1:-1]) 
                    / (np.log(phalf[2:]) - np.log(phalf[1:-1]))
                    )

                newf.add_dim_with_content(
                    "pfull", 
                    pfull,
                    "ref full pressure level", 
                    "Pa"
                    )
                newf.add_dim_with_content(
                    "phalf", 
                    phalf,
                    "ref half pressure level", 
                    "Pa"
                    )
                newf.log_axis1D(
                    "pk", 
                    pk, 
                    ("phalf"),
                    longname_txt = ("pressure part of the hybrid coordinate"),
                    units_txt = "Pa", 
                    cart_txt = ""
                    )
                newf.log_axis1D(
                    "bk", 
                    bk, 
                    ("phalf"),
                    longname_txt = ("sigma part of the hybrid coordinate"),
                    units_txt = "Pa", 
                    cart_txt = ""
                    )
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
        source_file = os.path.join(sys.prefix, "mars_data", "Legacy.fixed.nc")
        dest_file = f"{fdate}.fixed.nc"
        import shutil
        try:
            shutil.copy2(source_file, dest_file)
            print(f"{Cyan}{os.path.join(os.getcwd(), dest_file)} was copied locally{Nclr}")
        except OSError as e:
            print(f"{Red}Error copying {source_file} to {dest_file}: {e}{Nclr}")
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

        # longname_txt, units_txt = get_longname_unit(histfile, vname)
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
                time0 = (
                    ls2sol_1year(npvar[0]) + np.linspace(0, 10., len(npvar))
                    )

                if avgtime:
                    varnew = np.mean(npvar.reshape(-1, bin_period), axis=1)
                    time0 = np.mean(time0.reshape(-1, bin_period), axis=1)

                if not avgtime and not avgtod:
                    # Daily file
                    # Solar longitude
                    ls_start = npvar[0]
                    ls_end = npvar[-1]
                    step = (
                        (ls_end-ls_start) / np.float32(((numt-1) * ntod.size))
                        )
                    varnew = np.arange(0, numt * ntod.size, dtype=np.float32)
                    varnew[:] = varnew[:]*step + ls_start

                    # Time
                    step = (
                        (ls2sol_1year(ls_end) - ls2sol_1year(ls_start))
                        / np.float32((numt * ntod.size))
                        )
                    time0 = np.arange(0, numt * ntod.size, dtype=np.float32)
                    time0[:] = time0[:]*step + ls2sol_1year(ls_start)

                newf.log_axis1D(
                    "areo", 
                    varnew, 
                    dims,
                    longname_txt = "solar longitude",
                    units_txt = "degree",
                    cart_txt = "T"
                    )
                newf.log_axis1D(
                    "time", 
                    time0, 
                    dims,
                    longname_txt = "sol number",
                    units_txt = "days since 0000-00-00 00:00:00",
                    cart_txt = "T"
                    )
            else:
                continue

        elif ndims == 4:
            varnew = npvar
            if avgtime:
                varnew = np.mean(
                    npvar.reshape(
                        -1, bin_period, vshape[1], vshape[2], vshape[3]
                        ), 
                    axis = 1
                    )
            if avgtod:
                varnew = varnew.mean(axis = 1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3])
            # Rename variable
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt
                )
            # Convert surface pressure from mbar -> Pa
            if vname2 == "ps":
                varnew *= 100.
            newf.log_variable(
                vname2, 
                varnew, 
                newdims,
                longname_txt2,
                units_txt2
                )
        elif ndims == 5:
            varnew = npvar
            if avgtime:
                varnew = np.mean(
                    npvar.reshape(
                        -1, bin_period, vshape[1], vshape[2], vshape[3], vshape[4]
                        ), 
                    axis = 1
                    )
            if avgtod:
                varnew = varnew.mean(axis = 1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3], vshape[4])
            # Rename variables
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt
                )
            newf.log_variable(
                vname2, 
                varnew, 
                newdims,
                longname_txt2,
                units_txt2
                )
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
                newf.log_variable(
                    vname2, 
                    npvar, 
                    newdims, 
                    longname_txt2,
                    units_txt2
                    )
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
            newdims = replace_at_index(
                newdims, newdims.index("ntod"), "time_of_day_16"
                )
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

    ..note::
        For the moment, this is consistent with 0 <= Ls <= 359.99, but
        not for monotically increasing Ls.

    """
    Ls_perihelion = 250.99    # Ls at perihelion
    tperi = 485.35  # Time (in sols) at perihelion
    Ns = 668.6      # Number of sols in 1 MY
    e = 0.093379    # From MGCM: modules.f90
    nu = (Ls_deg - Ls_perihelion)*np.pi/180
    if nu == np.pi:
        nu = nu + 1e-10  # Adding epsilon of 10^-10
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

# ------------------------------------------------------
#                  END OF PROGRAM
# ------------------------------------------------------
if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
