#!/usr/bin/env python3
"""
The MarsPlot executable is for generating plots from Custom.in template
files. It sources variables from netCDF files in a specified directory.

The executable requires:

    * ``[-template --generate_template]`` Generates a Custom.in template
    * ``[-i --inspect]``         Triggers ncdump-like text to console
    * ``[Custom.in]``            To create plots in Custom.in template

Third-party Requirements:

    * ``numpy``
    * ``netCDF4``
    * ``sys``
    * ``argparse``
    * ``os``
    * ``warnings``
    * ``subprocess``
    * ``matplotlib``

"""

# Make print statements appear in color
from amescap.Script_utils import (
    Yellow, Red, Purple, Nclr, Blue, Green
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import subprocess   # Run command-line commands
import warnings     # Suppress errors triggered by NaNs
import matplotlib
import re           # Regular expressions
import numpy as np
from pypdf import PdfReader, PdfWriter
from netCDF4 import Dataset, MFDataset
from warnings import filterwarnings
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# Force matplotlib NOT to load Xwindows backend
matplotlib.use("Agg")

# Allows operations in square brackets in Custom.in
from numpy import (abs, sqrt, log, exp, min, max, mean)

from matplotlib.ticker import (
    LogFormatter, NullFormatter, LogFormatterSciNotation,
    MultipleLocator
)

# Load amesCAP modules
from amescap.Script_utils import (
    check_file_tape, section_content_amescap_profile, print_fileContent,
    print_varContent, FV3_file_type, find_tod_in_diurn, wbr_cmap,
    rjw_cmap, dkass_temp_cmap,dkass_dust_cmap,hot_cold_cmap
)
from amescap.FV3_utils import (
    lon360_to_180, lon180_to_360, UT_LTtxt, area_weights_deg,
    shiftgrid_180_to_360, shiftgrid_360_to_180, add_cyclic,
    azimuth2cart, mollweide2cart, robin2cart, ortho2cart
)

# Ignore deprecation warnings
filterwarnings("ignore", category = DeprecationWarning)

degr = u"\N{DEGREE SIGN}"
global current_version
current_version = 3.5

import functools
import traceback

import platform
import shutil

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

# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
    prog=('MarsPlot'),
    description=(
        f"{Yellow}MarsPlot V{current_version} is the plotting routine "
        f"for CAP.\nTo get started, use the -template flag to generate "
        f"a Custom.in template file.\nThen, use the Custom.in file to "
        f"create plots.\n"
        f"{Green}Example:\n"
        f"> MarsPlot -template {Blue}generates Custom.in file\n"
        f"modify the Custom.in file to generate desired plots{Green}\n"
        f"> MarsPlot Custom.in {Blue}generates pdf of plots from the "
        f"template"
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('template_file', nargs='?',
    type=argparse.FileType('r'),
    help=(
        f"Pass a template file to MarsPlot to create figures.\n"
        f"Must be a '.in' file.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in\n"
        f"> MarsPlot my_template.in"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-i', '--inspect_file', nargs='?',
    type=argparse.FileType('rb'),
    help=(
        f"Print the content of a netCDF file to the screen.\nVariables "
        f"are sorted by dimension.\n"
        f"Works on ANY netCDF file, including 'daily', diurn', "
        f"'average', and 'fixed'\n"
        f"{Green}Example:\n"
        f"> MarsPlot -i 01336.atmos_daily.nc"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-template', '--generate_template', default=False,
    action='store_true',
    help=(
        f"Generate a file called Custom.in that provides templates "
        f"for making plots with CAP.\n"
        f"{Green}Example:\n"
        f"> MarsPlot -template\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-d', '--date', nargs=1, default=None,
    help=(
        f"Specify the file to use. Default is the last file created.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -d 01336"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-sy', '--stack_years', action='store_true',
    default=False,
    help=(
        f"Plot consecutive years of data over the same axes range (e.g."
        f"Ls=0–360). For 1D time series plots only.\nRequires ADD LINE "
        f"in the Custom.in template (see template for instructions).\n"
        f"Default action is to plot in monotonically increasing "
        f"format.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -sy"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ftype', '--figure_filetype', default=None,
    type=str, choices=['pdf', 'eps', 'png'],
    help=(
        f"Output file format.\n Default is PDF else PNG.\n"
        f"Supported formats: PDF, EPS, PNG.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -ftype png"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-portrait', '--portrait_mode', action='store_true',
    default=False,
    help=(
        f"Output figures in portrait instead of landscape format.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -portrait"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-pw', '--pixel_width', default=2000, type=float,
    help=(
        f"Pixel width of the output figure. Default is 2000.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -pw 1000"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-dir', '--directory', default=os.getcwd(),
    help=(
        f"Target directory if input files are not in current "
        f"directory.\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in -dir path/to/directory"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

# to be used jointly with --generate_template
parser.add_argument('-trim', '--trim_text', action='store_true',
    default=False,
    help=(
        f"Generate a file called Custom.in that provides templates "
        f"for making plots\nwith CAP without the instructions "
        f"at top of the file.\n"
        f"{Green}Example:\n"
        f"> MarsPlot -template -trim\n"
        f"{Nclr}\n\n"
    )
)

# to be used jointly with --inspect
parser.add_argument('-values', '--print_values', nargs='+',
    default=None,
    help=(
        f"For use with ``-i --inspect``:\nPrint the values of the "
        f"specified variable to the screen.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsPlot -i 01336.atmos_daily.nc -values temp\n"
        f"{Blue}(quotes '' req. for browsing dimensions){Green}\n"
        f"> MarsPlot -i 01336.atmos_daily.nc -values 'temp[6,:,30,10]'"
        f"{Nclr}\n\n"
    )
)

# to be used jointly with --inspect
parser.add_argument('-stats', '--statistics', nargs='+', default=None,
    help=(
        f"For use with ``-i --inspect``:\nPrint the min, mean, and max "
        f"values of the specified variable to the screen.\n"
        f"Works on 'daily', 'diurn', and 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsPlot -i 01336.atmos_daily.nc -stats temp\n"
        f"{Blue}(quotes '' req. for browsing dimensions){Green}\n"
        f"> MarsPlot -i 01336.atmos_daily.nc -stats 'temp[6,:,30,10]'"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to print all Python errors and "
        f"status messages to the screen\n"
        f"{Green}Example:\n"
        f"> MarsPlot Custom.in --debug"
        f"{Nclr}\n\n"
    )
 )

# Handle mutually in/exclusive arguments (e.g., -sy requires Custom.in)
args = parser.parse_args()
debug = args.debug

if args.template_file:
    if not (re.search(".in", args.template_file.name) or re.search(".nc", args.template_file.name)):
        parser.error(f"{Red}Template file is not a '.in' or a netCDF file{Nclr}")
        exit()

if args.inspect_file:
    if not re.search(".nc", args.inspect_file.name):
        parser.error(f"{Red}{args.inspect_file.name} is not a netCDF "
                     f"file{Nclr}")
        exit()

if args.date is not None and (args.template_file is None and
                              args.generate_template is False and
                              args.inspect_file is None):
    parser.error(f"{Red}The -d argument requires a template file "
                 f"like Custom.in (e.g., MarsPlot Custom.in -d 01336)"
                 f"{Nclr}")
    exit()

if args.figure_filetype is not None and (
    args.template_file is None and
    args.generate_template is False and
    args.inspect_file is None
    ):
    parser.error(f"{Red}The -f argument requires a template file "
                 f"like Custom.in (e.g., MarsPlot Custom.in -ftype png)"
                 f"{Nclr}")
    exit()

if args.stack_years and (args.template_file is None and
                         args.generate_template is False and
                         args.inspect_file is None):
    parser.error(f"{Red}The -sy argument requires a template file "
                 f"like Custom.in (e.g., MarsPlot Custom.in -sy)"
                 f"{Nclr}")
    exit()

if args.portrait_mode and (args.template_file is None and
                           args.generate_template is False and
                           args.inspect_file is None):
    parser.error(f"{Red}The -portrait argument requires a template "
                 f"file like Custom.in (e.g., MarsPlot Custom.in "
                 f"-portrait){Nclr}")
    exit()

if args.statistics is not None and args.inspect_file is None:
    parser.error(f"{Red}The -stat argument requires a template "
                 f"file like Custom.in (e.g., MarsPlot -i "
                 f"01336.atmos_daily.nc -stats temp{Nclr}")
    exit()

if args.print_values is not None and args.inspect_file is None:
    parser.error(f"{Red}The -values argument requires a template "
                 f"file like Custom.in (e.g., MarsPlot -i "
                 f"01336.atmos_daily.nc -values temp{Nclr}")
    exit()

if args.trim_text and (args.generate_template is False):
    parser.error(f"{Red}The -trim argument requires -template (e.g., "
                 f"MarsPlot -template -trim{Nclr}")
    exit()

# ======================================================================
#                           MAIN PROGRAM
# ======================================================================
@debug_wrapper
def main():
    global output_path, input_paths, out_format, debug
    output_path = os.getcwd()
    out_format = 'pdf' if args.figure_filetype is None else args.figure_filetype
    debug = args.debug
    input_paths = []
    input_paths.append(args.directory)

    global Ncdf_num         # Hosts the simulation timestamps
    global objectList       # Contains all figure objects
    global customFileIN     # The Custom.in template name
    global levels, my_dpi, label_size, title_size, label_factor
    global tick_factor, title_factor

    levels = 21             # Number of contours for 2D plots
    my_dpi = 96.            # Pixels per inch for figure output
    label_size = 18         # Label size for title, xlabel, and ylabel
    title_size = 24         # Label size for title, xlabel, and ylabel
    label_factor = 3/10     # Reduce font size as # of panels increases
    tick_factor = 1/2
    title_factor = 10/12

    global width_inch       # Pixel width for saving figure
    global height_inch      # Pixel width for saving figure
    global vertical_page

    # Set portrait format for outout figures
    vertical_page = args.portrait_mode

    # Directory (dir) containing shared templates
    global shared_dir
    shared_dir = "/path_to_shared_templates"

    # Set figure dimensions
    pixel_width = args.pixel_width
    if vertical_page:
        width_inch = pixel_width / 1.4 / my_dpi
        height_inch = pixel_width / my_dpi
    else:
        width_inch = pixel_width / my_dpi
        height_inch = pixel_width / 1.4 / my_dpi

    objectList = [Fig_2D_lon_lat("fixed.zsurf", True),
                  Fig_2D_lat_lev("atmos_average.ucomp", True),
                  Fig_2D_time_lat("atmos_average.taudust_IR", False),
                  Fig_2D_lon_lev("atmos_average_pstd.temp", False),
                  Fig_2D_time_lev("atmos_average_pstd.temp", False),
                  Fig_2D_lon_time("atmos_average.temp", False),
                  Fig_1D("atmos_average.temp", False)]

    # Group together the first two figures
    objectList[0].subID = 1 # 1st object in a 2-panel figure
    objectList[0].nPan = 2
    objectList[1].subID = 2 # 2nd object in a 2-panel figure
    objectList[1].nPan = 2

    if args.inspect_file:
        # --inspect: Inspect content of netcdf file
        print("Attempting to access file:", args.inspect_file)
        # NAS-specific, check if the file is on tape (Lou only)
        check_file_tape(args.inspect_file)
        if args.print_values:
            # Print variable content to screen
            print_varContent(args.inspect_file,
                             args.print_values, False)
        elif args.statistics:
            # Print variable stats (max, min, mean) to screen
            print_varContent(args.inspect_file,
                             args.statistics, True)
        else:
            # Show information for all variables
            print("doing inspect")
            print_fileContent(args.inspect_file)

    elif args.generate_template:
        make_template()

    elif args.template_file:
        # Case A: Use local Custom.in (most common option)
        print(f"Reading {args.template_file.name}")
        namelist_parser(args.template_file.name)

        if args.date:
            # If optional --date provided, use files matching date(s)
            try:
                # Confirm that input date type = float
                bound = np.asarray(args.date).astype(float)
            except Exception as e:
                print(f"{Red}*** Syntax Error***\nPlease use: ``MarsPlot "
                      f"Custom.in -d XXXX [YYYY] -o out``{Nclr}")
                exit()
        else:
            # If NO --date, default to date of most recent fixed file
            bound = get_Ncdf_num()
            if bound is not None:
                bound = bound[-1]

        # Extract all timestamps in dir
        Ncdf_num = get_Ncdf_num()

        if Ncdf_num is not None:
            # Apply bounds to desired date
            Ncdf_num = select_range(Ncdf_num, bound)

        # Make folder "plots" in cwd
        dir_plot_present = os.path.exists(os.path.join(output_path,"plots"))
        if not dir_plot_present:
            os.makedirs(os.path.join(output_path,"plots"))

        # ============ Update Progress Bar ============
        global i_list
        # Create list of figures
        fig_list = list()
        for i_list in range(0, len(objectList)):

            # Display status of figure in progress
            status = (f"{objectList[i_list].plot_type} :"
                      f"{objectList[i_list].varfull}")
            progress(i_list, len(objectList), status, None)

            objectList[i_list].do_plot()

            if (objectList[i_list].success and
                out_format == "pdf" and not
                debug):
                sys.stdout.write("\033[F")
                # Flush previous output
                sys.stdout.write("\033[K")

            status = (f"{objectList[i_list].plot_type}:"
                      f"{objectList[i_list].varfull}"
                      f"{objectList[i_list].fdim_txt}")
            progress(i_list, len(objectList), status,
                     objectList[i_list].success)

            if objectList[i_list].subID == objectList[i_list].nPan:
                if (i_list < len(objectList)-1 and not
                    objectList[i_list + 1].addLine):
                    fig_list.append(objectList[i_list].fig_name)
                if i_list == len(objectList)-1:
                    fig_list.append(objectList[i_list].fig_name)

        progress(100, 100, "Done")

        # ============ For Multipage PDF ============
        # Make multipage PDF out of figures in /plots. Remove individual 
        # plots. Debug files when complete.
        if out_format == "pdf" and len(fig_list) > 0:
            print("Merging figures...")
            # Construct string of figure names separated by spaces
            all_fig = " "
            for figID in fig_list:
                # Place outer quotes around figID to handle whitespaces
                # in Windows paths, i.e., "/Users/myfolder/plots.pdf"
                figID = (f'"{figID}"')
                all_fig += (f"{figID} ")

            try:
                # If template file = "Custom", use default
                # PDF basename "Diagnostics":
                # e.g., Custom.in -> Diagnostics.pdf, or
                #       Custom_01.in -> Diagnostics_01.pdf
                input_file = (os.path.join(output_path,
                                f"{args.template_file.name}"))

                if platform.system() == "Windows":
                    basename = input_file.split("\\")[-1].split(".")[0].strip()
                else:
                    basename = input_file.split("/")[-1].split(".")[0].strip()
            except:
                # Use default PDF basename "Diagnostics".
                basename = "Custom"

            # Generate PDF name
            if basename == "Custom":
                # If template name = Custom.in -> Diagnostics.pdf
                output_pdf = os.path.join(output_path,"Diagnostics.pdf")
            elif basename[0:7] == "Custom_":
                # If template name = Custom_XX.in -> Diagnostics_XX.pdf
                output_pdf = os.path.join(output_path,f"Diagnostics_{basename[7:9]}.pdf")
            else:
                # If template name is NOT Custom.in, use prefix to
                # generate PDF name
                output_pdf = os.path.join(output_path,f"{basename}.pdf")

            # Add quotes around PDF name (name -> "name")
            output_pdfq = f'"{output_pdf}"'

            # Direct gs output to file instead of printing to screen
            debug_filename = os.path.join(output_path,f".debug_MCMC_plots.txt")
            fdump = open(debug_filename, "w")

            writer = PdfWriter()

            try:
                for pdf_file in fig_list:
                    reader = PdfReader(pdf_file)
                    for page in reader.pages:
                        writer.add_page(page)

                with open(output_pdf, "wb") as f:
                    writer.write(f)
                give_permission(output_pdf)
                print(f"{output_pdfq} was generated")
            except subprocess.CalledProcessError:
                # If gs command fails again, prompt user to try
                # generating PNG instead
                print("ERROR with merging PDF, please "
                      "try a different format, such as PNG.")
                if debug:
                    raise
    else:
        parser.error(f"{Red}No valid argument was passed. Pass a "
                     f"Custom.in template file or use -template or -i "
                     f"to use MarsPlot. Type 'MarsPlot -h' if you need "
                     f"more assistance.{Nclr}")
        exit()
# ======================================================================
#                       DATA OPERATION UTILITIES
# ======================================================================

# User Preferences from amescap_profile
global add_sol_time_axis, lon_coord_type, include_NaNs

# Create a namespace with numpy available
namespace = {'np': np}
# Load preferences in Settings section of amescap_profile
exec(section_content_amescap_profile("MarsPlot Settings"), namespace)

# Determine whether to include sol number in addition to Ls on
# time axis. Default FALSE (= Ls only)
add_sol_time_axis = eval("np.array(add_sol_to_time_axis)", namespace)

# Define longitude coordinates to use. Default = 360 (i.e., 0-360).
# Alt. = 180 (i.e., -180-180)
lon_coord_type = eval("np.array(lon_coordinate)", namespace)

# Determine whether to include or ignore NaNs when computing means.
# Default FALSE = exclude NaNs (use np.nanmean)
# Alt. TRUE = include NaNs (use np.mean)
include_NaNs = eval("np.array(show_NaN_in_slice)", namespace)


def mean_func(arr, axis):
    """
    This function calculates a mean over the selected axis, ignoring or
    including NaN values as specified by ``show_NaN_in_slice`` in
    ``amescap_profile``.

    :param arr: the array to be averaged
    :type arr: array
    
    :param axis: the axis over which to average the array
    :type axis: int

    :return: the mean over the time axis

    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category = RuntimeWarning)
        if include_NaNs:
            return np.mean(arr, axis = axis)
        else:
            return np.nanmean(arr, axis = axis)


def shift_data(lon, data):
    """
    Shifts the longitude data from 0-360 to -180/180 and vice versa.

    :param lon: 1D array of longitude
    :type lon: array [lon]
    
    :param data: 2D array with last dimension = longitude
    :type data: array [1,lon]
    
    :raises ValueError: Longitude coordinate type is not recognized.
    
    :return: longitude (-180/180)
    :rtype: array [lon]
    
    :return: shifted data
    :rtype: array [1,lon]
    
    .. note::
        Use ``np.ma.hstack`` instead of ``np.hstack`` to keep the
        masked array properties

    """
    nlon = len(lon)
    # If 1D plot, reshape array
    if len(data.shape) <= 1:
        data = data.reshape(1, nlon)
    if lon_coord_type == 180:
        lon_out=lon360_to_180(lon)
        data = shiftgrid_360_to_180(lon, data)
    elif lon_coord_type == 360:
        lon_out=lon180_to_360(lon)
        data = shiftgrid_180_to_360(lon, data)
    else:
        raise ValueError("Longitude coordinate type invalid. Please "
                         "specify ``180`` or ``360`` after "
                         "``lon_coordinate`` in ``amescap_profile``")
    # If 1D plot, squeeze array
    if data.shape[0] == 1:
        data = np.squeeze(data)
    return lon_out, data


def MY_func(Ls_cont):
    """
    Returns the Mars Year

    :param Ls_cont: solar longitude (``areo``; continuous)
    :type Ls_cont: array [areo]
    
    :return: the Mars year
    :rtype: int

    """
    return (Ls_cont)//(360.)+1


def get_lon_index(lon_query_180, lons):
    """
    Returns the indices that will extract data from the netCDF file
    according to a range of *longitudes*.

    :param lon_query_180: longitudes in -180/180: value,
        ``[min, max]``, or `None`
    :type lon_query_180: list
    
    :param lons: longitude in 0-360
    :type lons: array [lon]
    
    :return: 1D array of file indices
    :rtype: array
    
    :return: text descriptor for the extracted longitudes
    :rtype: str
    
    .. note::
        The keyword ``all`` passed as ``-99999`` by the rT() functions

    """
    Nlon = len(lons)
    lon_query_180 = np.array(lon_query_180)

    if None in lon_query_180:
        # Zonal average
        lon_query_180 = np.array(-99999)

    if lons.max() > 180:
        # ============== FV3 format ==============
        # If lon = 0-360, convert to -180/180
        # ========================================
        if lon_query_180.size == 1:
            # If one longitude provided
            if lon_query_180 == -99999:
                # Request zonal average
                loni = np.arange(0, Nlon)
                txt_lon = ", zonal avg"
            else:
                # Get closest value
                lon_query_360 = lon180_to_360(lon_query_180)
                loni = np.argmin(abs(lon_query_360-lons))
                txt_lon = f", lon={lon360_to_180(lons[loni]):.1f}"

        elif lon_query_180.size == 2:
            # If range of longitudes provided
            lon_query_360 = lon180_to_360(lon_query_180)
            loni_bounds = np.array([np.argmin(abs(lon_query_360[0]-lons)),
                                    np.argmin(abs(lon_query_360[1]-lons))])
            # Longitude should be increasing for extraction # TODO
            # Normal case (e.g., -45°W > 45°E)
            if loni_bounds[0] < loni_bounds[1]:
                loni = np.arange(loni_bounds[0], loni_bounds[1] + 1)
            else:
                # Loop around (e.g., 160°E > -40°W)
                loni = np.append(np.arange(loni_bounds[0], len(lons)),
                                 np.arange(0, loni_bounds[1] + 1))
                print(f"{Purple}lon360_to_180(lons[loni]){Nclr}")

            lon_bounds_180 = lon360_to_180([lons[loni_bounds[0]],
                                            lons[loni_bounds[1]]])
            # Longitude should be increasing for display
            txt_lon = (f", lon=avg[{lon_bounds_180[0]:.1f}"
                       f"<->{lon_bounds_180[1]:.1f}]")
    else:
        # =========== Legacy Format ===========
        # Lon = -180/180
        # =====================================
        if lon_query_180.size == 1:
            # If one longitude provided
            if lon_query_180 == -99999:
                # Zonal average
                loni = np.arange(0, Nlon)
                txt_lon = ", zonal avg"
            else:
                # Get closest value
                loni = np.argmin(abs(lon_query_180-lons))
                txt_lon = f", lon={lons[loni]:.1f}"

        elif lon_query_180.size == 2:
            # If range of longitudes provided
            loni_bounds = np.array([np.argmin(abs(lon_query_180[0]-lons)),
                                    np.argmin(abs(lon_query_180[1]-lons))])
            if loni_bounds[0] < loni_bounds[1]:
                # Normal case (e.g., -45 °W > 45 °E)
                loni = np.arange(loni_bounds[0], loni_bounds[1]+1)
            else:
                # Loop around (e.g., 160°E > -40°W)
                loni = np.append(np.arange(loni_bounds[0], len(lons)),
                                 np.arange(0, loni_bounds[1]+1))
            txt_lon = (f", lon=avg[{lons[loni_bounds[0]]:.1f}"
                       f"<->{lons[loni_bounds[1]]:.1f}]")
    #.. note:: is lon dimension is degenerate, e.g. (time,lev,lat,1)
    #loni must be a scalar, otherwise f.variables['var'][time,lev,lat,loni] returns an error
    if len(np.atleast_1d(loni))==1 and not np.isscalar(loni):
        loni=loni[0]
    return loni, txt_lon


def get_lat_index(lat_query, lats):
    """
    Returns the indices that will extract data from the netCDF file
    according to a range of *latitudes*.

    :param lat_query: requested latitudes (-90/+90)
    :type lat_query: list
    
    :param lats: latitude
    :type lats: array [lat]
    
    :return: 1d array of file indices
    :rtype: text descriptor for the extracted longitudes
    
    .. note::T
        The keyword ``all`` passed as ``-99999`` by the ``rt()``
        function

    """
    Nlat = len(lats)
    lat_query = np.array(lat_query)

    if None in lat_query:
        # Default to equator
        lat_query = np.array(0.)
    if lat_query.size == 1:
        # If one latitude provided
        if lat_query == -99999:
            # Meridional average
            lati = np.arange(0, Nlat)
            txt_lat = ", merid. avg"
        else:
            # Get closest value
            lati = np.argmin(abs(lat_query-lats))
            txt_lat = f", lat={lats[lati]:g}"
    elif lat_query.size == 2:
        # If range of latitudes provided
        lat_bounds = np.array([np.argmin(abs(lat_query[0] - lats)),
                               np.argmin(abs(lat_query[1] - lats))])
        if lat_bounds[0] > lat_bounds[1]:
            # Latitude should be increasing for extraction
            lat_bounds = np.flipud(lat_bounds)
        lati = np.arange(lat_bounds[0], lat_bounds[1]+1)
        txt_lat = f", lat=avg[{lats[lati[0]]:g}<->{lats[lati[-1]]:g}]"
    return lati, txt_lat


def get_tod_index(tod_query, tods):
    """
    Returns the indices that will extract data from the netCDF file
    according to a range of *times of day*.

    :param tod_query: requested time of day (0-24)
    :type tod_query: list
    
    :param tods: times of day
    :type tods: array [tod]
    
    :return: file indices
    :rtype: array [tod]
    
    :return: descriptor for the extracted time of day
    :rtype: str
    
    .. note::
        The keyword ``all`` is passed as ``-99999`` by the ``rT()``
        function

    """
    Ntod = len(tods)
    tod_query = np.array(tod_query)

    if None in tod_query:
        # Default to 3 pm
        tod_query = np.array(15)

    if tod_query.size == 1:
        # If one time of day provided
        if tod_query == -99999:
            # Diurnal average
            todi = np.arange(0, Ntod)
            txt_tod = ", tod avg"
        else:
            # Get closest value
            todi = np.argmin(abs(tod_query-tods))
            txt_tmp = UT_LTtxt(tods[todi]/24., lon_180 = 0., roundmin = 1)
            txt_tod = f", tod= {txt_tmp}"
    elif tod_query.size == 2:
        # If range of times of day provided
        tod_bounds = np.array([np.argmin(abs(tod_query[0] - tods)),
                               np.argmin(abs(tod_query[1] - tods))])
        if tod_bounds[0] < tod_bounds[1]:
            # Normal case (e.g., 4 am > 10am)
            todi = np.arange(tod_bounds[0], tod_bounds[1]+1)
        else:
            # Loop around (e.g., 18 pm > 6 am)
            todi = np.append(np.arange(tod_bounds[0], len(tods)),
                             np.arange(0, tod_bounds[1]+1))
        txt_tmp = UT_LTtxt(tods[todi[0]]/24., lon_180 = 0., roundmin = 1)
        txt_tmp2 = UT_LTtxt(tods[todi[-1]]/24., lon_180 = 0., roundmin = 1)
        txt_tod = f", tod=avg[{txt_tmp}<->{txt_tmp2}]"
    return todi, txt_tod


def get_level_index(level_query, levs):
    """
    Returns the indices that will extract data from the netCDF file
    according to a range of *pressures* (resp. depth for ``zgrid``).

    :param level_query: requested pressure [Pa] (depth [m])
    :type level_query: float
    
    :param levs: levels (in the native coordinates)
    :type levs: array [lev]
    
    :return: file indices
    :rtype: array
    
    :return: descriptor for the extracted pressure (depth)
    :rtype: str
    
    .. note::
        The keyword ``all`` is passed as ``-99999`` by the ``rT()``
        functions

    """
    level_query = np.array(level_query)
    Nz = len(levs)

    if None in level_query :
        # Default to surface
        # If level_query >>> Psfc (even for a 10-bar Early Mars sim)
        level_query = np.array(2*10**7)
    if level_query.size == 1:
        # If one level provided
        if level_query == -99999:
            # Average
            levi = np.arange(0, Nz)
            txt_level = ", column avg"
        else:
            # Specific level
            levi = np.argmin(abs(level_query-levs))
            if level_query > 10.**7:
                # Provide smart labeling
                # None (i.e.surface was requested)
                txt_level = ", at sfc"
            else:
                txt_level = f", lev={levs[levi]:1.2e} Pa/m"
    elif level_query.size == 2:
        # Bounds are provided
        levi_bounds = np.array([np.argmin(abs(level_query[0] - levs)),
                                np.argmin(abs(level_query[1] - levs))])
        if levi_bounds[0] > levi_bounds[1]:
            # Level should be increasing for extraction
            levi_bounds = np.flipud(levi_bounds)
        levi = np.arange(levi_bounds[0], levi_bounds[1]+1)
        lev_bounds = [levs[levi[0]], levs[levi[-1]]]
        if lev_bounds[0] < lev_bounds[1]:
            # Level should be decreasing for display
            lev_bounds = np.flipud(lev_bounds)
        txt_level = (f", lev=avg[{lev_bounds[0]:1.2e}"
                     f"<->{lev_bounds[1]:1.2e}] Pa/m")
    return levi, txt_level


def get_time_index(Ls_query_360, LsDay):
    """
    Returns the indices that will extract data from the netCDF file
    according to a range of solar longitudes [0-360].

    First try the Mars Year of the last timestep, then try the year
    before that. Use whichever Ls period is closest to the requested
    date.

    :param Ls_query_360: requested solar longitudes
    :type Ls_query_360: list
    
    :param LsDay: continuous solar longitudes
    :type LsDay: array [areo]
    
    :return: file indices
    :rtype: array
    
    :return: descriptor for the extracted solar longitudes
    :rtype: str
    
    .. note::
        The keyword ``all`` is passed as ``-99999`` by the ``rT()``
        function

    """
    if len(np.atleast_1d(LsDay)) == 1:
        # Special case: file has 1 timestep, transform LsDay -> array
        LsDay = np.array([LsDay])

    Nt = len(LsDay)
    Ls_query_360 = np.array(Ls_query_360)
    if None in Ls_query_360:
        # Defaultto last timestep
        Ls_query_360 = np.mod(LsDay[-1], 360.)
    if Ls_query_360.size == 1:
        # If one time provided
        if Ls_query_360 == -99999:
            # Time average
            ti = np.arange(0, Nt)
            txt_time = ", time avg"
        else:
            # Get Mars Year (MY) of last timestep in file
            MY_end = MY_func(LsDay[-1])
            if MY_end >= 1:
                # Check if desired Ls available in this MY
                Ls_query = Ls_query_360 + (MY_end - 1)*360.
                # MY starts at 1, not 0
            else:
                Ls_query = Ls_query_360

            if Ls_query > LsDay[-1] and MY_end > 1:
                # Lok one year back
                MY_end -= 1
                Ls_query = Ls_query_360 + (MY_end - 1)*360.
            ti = np.argmin(abs(Ls_query - LsDay))
            txt_time = f", Ls= (MY{MY_end:02}) {np.mod(LsDay[ti], 360.):.2f}"
    elif Ls_query_360.size == 2:
        # If a range of times provided
        MY_last = MY_func(LsDay[-1])
        if MY_last >= 1:
            # Get MY of last timestep
            Ls_query_last = Ls_query_360[1] + (MY_last-1)*360.
        else:
            Ls_query_last = Ls_query_360[1]
        if Ls_query_last > LsDay[-1] and MY_last > 1:
            # Look one MY back
            MY_last -= 1
            Ls_query_last = Ls_query_360[1] + (MY_last-1)*360.
        ti_last = np.argmin(abs(Ls_query_last - LsDay))
        # Then get first value for that MY
        MY_beg = MY_last.copy()

        # Try MY of last timestep
        Ls_query_beg = Ls_query_360[0] + (MY_beg-1)*360.
        ti_beg = np.argmin(abs(Ls_query_beg - LsDay))

        if ti_beg >= ti_last:
            # Search year before for ti_beg
            MY_beg -= 1
            Ls_query_beg = Ls_query_360[0] + (MY_beg-1)*360.
            ti_beg = np.argmin(abs(Ls_query_beg - LsDay))

        ti = np.arange(ti_beg, ti_last + 1)
        Ls_bounds = [LsDay[ti[0]], LsDay[ti[-1]]]
        txt_time = (f", Ls= avg [(MY{MY_beg:02}) "
                    f"{np.mod(Ls_bounds[0], 360.):.2f} <-> (MY{MY_last:02}) "
                    f"{np.mod(Ls_bounds[1], 360.):.2f}]")
    return ti, txt_time

# ======================================================================
#                          TEMPLATE UTILITIES
# ======================================================================
def filter_input(txt, typeIn="char"):
    """
    Read template for the type of data expected

    :param txt: text input into ``Custom.in`` to the right of an equal
        sign
    :type txt: str
    
    :param typeIn: type of data expected: ``char``, ``float``, ``int``,
        ``bool``, defaults to ``char``
    :type typeIn: str, optional
    
    :return: text input reformatted to ``[val1, val2]``
    :rtype: float or array

    """
    if txt == "None" or not txt:
        # If None or empty string
        return None

    if "," in txt:
        # If 2 values provided
        answ = []
        for i in range(0, len(txt.split(","))):
            # For a char, read all text as one
            if typeIn == "char":
                answ = txt
            if typeIn == "float":
                answ.append(float(txt.split(",")[i].strip()))
            if typeIn == "int":
                answ.append(int(txt.split(",")[i].strip()))
            if typeIn == "bool":
                answ.append(txt.split(",")[i].strip() == "True")
        return answ
    else:
        if typeIn == "char":
            answ = txt
        if typeIn == "bool":
            answ = ("True" == txt)
        if typeIn == "float":
            # Pass the all key words as -99999
            if txt == "all":
                answ = -99999.
            elif txt == "AXIS":
                answ = -88888.
            else:
                answ = float(txt)
        if typeIn == "int":
            # Pass the all key words as -99999
            if txt == "all":
                answ = -99999
            else:
                # True if text matches
                answ = int(txt)
        return answ


def rT(typeIn="char"):
    """
    Read template for the type of data expected. Returns value to
    ``filter_input()``.

    :param typeIn: type of data expected: ``char``, ``float``, ``int``,
        ``bool``, defaults to ``char``
    :type typeIn: str, optional
    
    :return: text input reformatted to ``[val1, val2]``
    :rtype: float or array

    """
    global customFileIN
    raw_input = customFileIN.readline()

    if len(raw_input.split("=")) == 2:
        # Get text on right side of equal sign if 1
        # equal sign in string (e.g., 02400.atmos_average2{lat=20})
        txt = raw_input.split("=")[1].strip()

    elif len(raw_input.split("=")) > 2:
        # Read string manually if 1+ equal signs
        # (e.g., 02400.atmos_average2{lat=20,tod=4})
        current_varfull = ""
        record = False
        for i in range(0, len(raw_input)):
            if record:
                current_varfull += raw_input[i]
            if raw_input[i] == "=":
                record = True
        txt = current_varfull.strip()
    return filter_input(txt, typeIn)


def read_axis_options(axis_options_txt):
    """
    Return axis customization options.

    :param axis_options_txt: a copy of the last line ``Axis Options``
        in ``Custom.in`` templates
    :type axis_options_txt: str
    
    :return: X-axis bounds as a numpy array or ``None`` if undedefined
    :rtype: array or None
    
    :return: Y-axis bounds as a numpy array or ``None`` if undedefined
    :rtype: array or None
    
    :return: colormap (e.g., ``jet``, ``nipy_spectral``) or line
        options (e.g., ``--r`` for dashed red)
    :rtype: str
    
    :return: linear (``lin``) or logarithmic (``log``) color scale
    :rtype: str
    
    :return: projection (e.g., ``ortho -125,45``)
    :rtype: str

    """
    list_txt = axis_options_txt.split(":")[1].split("|")

    # Xaxis: get bounds
    txt = list_txt[0].split("=")[1].replace("[", "").replace("]", "")
    Xaxis = []
    for i in range(0, len(txt.split(","))):
        if txt.split(",")[i].strip() == "None":
            Xaxis = None
            break
        else:
            Xaxis.append(float(txt.split(",")[i].strip()))

    # Yaxis: get bounds
    txt = list_txt[1].split("=")[1].replace("[", "").replace("]", "")
    Yaxis = []
    for i in range(0, len(txt.split(","))):
        if txt.split(",")[i].strip() == "None":
            Yaxis = None
            break
        else:
            Yaxis.append(float(txt.split(",")[i].strip()))

    # Line or colormap
    custom_line1 = list_txt[2].split("=")[1].strip()
    custom_line2 = None
    custom_line3 = None

    # Scale: lin or log (2D plots only)
    if len(list_txt) == 4:
        custom_line2 = list_txt[3].split("=")[1].strip()
        if custom_line2.strip() == "None":
            custom_line2 = None
    if len(list_txt) == 5:
        custom_line2 = list_txt[3].split("=")[1].strip()
        custom_line3 = list_txt[4].split("=")[1].strip()
        if custom_line2.strip() == "None":
            custom_line2 = None
        if custom_line3.strip() == "None":
            custom_line3 = None
    return Xaxis, Yaxis, custom_line1, custom_line2, custom_line3


def split_varfull(varfull):
    """
    Split ``varfull`` object into its component parts

    :param varfull: a ``varfull`` object (e.g,
        ``atmos_average@2.zsurf``, ``02400.atmos_average@2.zsurf``)
    :type varfull: str
    
    :return: (sol_array) a sol number or ``None`` (if none provided)
    :rtype: int or None
    
    :return: (filetype) file type (e.g, ``atmos_average``)
    :rtype: str
    
    :return: (var) variable of interest (e.g, ``zsurf``)
    :rtype: str
    
    :return: (``simuID``) simulation ID (Python indexing starts at 0)
    :rtype: int

    """
    if varfull.count(".") == 1:
        # Default: no sol number provided (e.g., atmos_average2.zsurf).
        # Extract variables and file from varfull
        sol_array = np.array([None])
        # File and ID
        filetypeID = varfull.split(".")[0].strip()
        # Variable name
        var = varfull.split(".")[1].strip()

    # Case 2: sol number provided (e.g., 02400.atmos_average2.zsurf)
    elif varfull.count(".") == 2:
        # Sol number
        sol_array = np.array([int(varfull.split(".")[0].strip())])
        # File and ID
        filetypeID = varfull.split(".")[1].strip()
        # Variable name
        var = varfull.split(".")[2].strip()
    # Split filename and simulation ID

    if "@" in filetypeID:
        filetype = filetypeID.split("@")[0].strip()
        # Simulation ID starts at 0 in the code
        simuID = int(filetypeID.split("@")[1].strip()) - 1
    else:
        # No digit (i.e., reference simulation)
        simuID = 0
        filetype = filetypeID
    return sol_array, filetype, var, simuID


def remove_whitespace(raw_input):
    """
    Remove whitespace inside an expression.

    This is different from the ``.strip()`` method, which only removes
    whitespaces at the edges of a string.

    :param raw_input: user input for variable, (e.g.,
        ``[atmos_average.temp] + 2)``
    :type raw_input: str
    
    :return: raw_input without whitespaces (e.g.,
        ``[atmos_average.temp]+2)``
    :rtype: str

    """
    processed_input = ""
    for i in range(0, len(raw_input)):
        if raw_input[i] != " ":
            processed_input += raw_input[i]
    return processed_input


def clean_comma_whitespace(raw_input):
    """
    Remove commas and whitespaces inside an expression.

    :param raw_input: dimensions specified by user input to Variable
        (e.g., ``lat=3. , lon=2 , lev = 10.``)
    :type raw_input: str
    
    :return: raw_input without whitespaces (e.g.,
        ``lat=3.,lon=2,lev=10.``)
    :rtype: str

    """
    processed_input = ""
    for i in range(0, len(raw_input)):
        if raw_input[i] != ",":
            processed_input += raw_input[i]
    return remove_whitespace(processed_input)


def get_list_varfull(raw_input):
    """
    Return requested variable from a complex ``varfull`` object with ``[]``.

    :param raw_input: complex user input to Variable (e.g.,
        ``2*[atmos_average.temp]+[atmos_average2.ucomp]*1000``)
    :type raw_input: str
    
    :return: list required variables (e.g., [``atmos_average.temp``,
        ``atmos_average2.ucomp``])
    :rtype: str

    """
    var_list = []
    record = False
    current_name = ""
    for i in range(0, len(raw_input)):
        if raw_input[i] == "]":
            record = False
            var_list.append(current_name.strip())
            current_name = ""
        if record:
            current_name += raw_input[i]
        if raw_input[i] == "[":
            record = True
    return var_list


def get_overwrite_dim_2D(varfull_bracket, plot_type, fdim1, fdim2, ftod):
    """
    Return new dimensions that will overwrite default dimensions for a
    varfull object with ``{}`` on a 2D plot.

    ``2D_lon_lat:  fdim1 = ls,  fdim2 = lev``
    ``2D_lat_lev:  fdim1 = ls,  fdim2 = lon``
    ``2D_time_lat: fdim1 = lon, fdim2 = lev``
    ``2D_lon_lev:  fdim1 = ls,  fdim2 = lat``
    ``2D_time_lev: fdim1 = lat, fdim2 = lon``
    ``2D_lon_time: fdim1 = lat, fdim2 = lev``

    :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
        ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
    :type varfull_bracket: str
    
    :param plot_type: the type of the plot template
    :type plot_type: str
    
    :param fdim1: X axis dimension for plot
    :type fdim1: str
    
    :param fdim2: Y axis dimension for plot
    :type fdim2: str
    
    :return: (varfull) required file and variable (e.g.,
        ``atmos_average.temp``);
        (fdim_out1) X axis dimension for plot;
        (fdim_out2) Y axis dimension for plot;
        (ftod_out) if X or Y axis dimension is time of day

    """
    # Initialization: use dimension provided in template
    fdim_out1 = fdim1
    fdim_out2 = fdim2

    # Left of "{":
    varfull_no_bracket = varfull_bracket.split(
        "{")[0].strip()

    # Right of "{", with last "}" removed:
    overwrite_txt = remove_whitespace(varfull_bracket.split("{")[1][:-1])

    # Count number of "=" in string
    ndim_update = overwrite_txt.count("=")

    # Split to different blocks (e.g., lat = 3. and lon = 20)
    split_dim = overwrite_txt.split(";")
    if overwrite_txt.count(";") < (overwrite_txt.count("=")-1):
        print(f"{Yellow}*** Error: use semicolon ';' to separate dimensions "
              f"'{{}}'{Nclr}")

    for i in range(0, ndim_update):
        # Check if the requested dimension exists:
        if (split_dim[i].split("=")[0] not in
            ["ls", "lev", "lon", "lat", "tod"]):
            print(f"{Yellow}*** Warning*** Ignoring dimension: "
                  f"{split_dim[i].split('=')[0]} because it is not recognized."
                  f"Valid dimensions = ls, lev, lon, lat, or tod{Nclr}")
        if plot_type == "2D_lon_lat":
            if split_dim[i].split("=")[0] == "ls":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lev":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")
        if plot_type == "2D_lat_lev":
            if split_dim[i].split("=")[0] == "ls":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lon":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")
        if plot_type == "2D_time_lat":
            if split_dim[i].split("=")[0] == "lon":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lev":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")
        if plot_type == "2D_lon_lev":
            if split_dim[i].split("=")[0] == "ls":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lat":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")
        if plot_type == "2D_time_lev":
            if split_dim[i].split("=")[0] == "lat":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lon":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")
        if plot_type == "2D_lon_time":
            if split_dim[i].split("=")[0] == "lat":
                fdim_out1 = filter_input(split_dim[i].split("=")[1], "float")
            if split_dim[i].split("=")[0] == "lev":
                fdim_out2 = filter_input(split_dim[i].split("=")[1], "float")

        ftod_out = None
        if split_dim[i].split("=")[0] == "tod":
            # Always get time of day
            ftod_out = filter_input(split_dim[i].split("=")[1], "float")
    # .. note:: filter_input() converts text (3 or 4, 5) to variable:
    # (e.g., numpy.array([3.]) or numpy.array([4., 5.]))
    return varfull_no_bracket, fdim_out1, fdim_out2, ftod_out


def get_overwrite_dim_1D(varfull_bracket, t_in, lat_in, lon_in, lev_in,
                         ftod_in):
    """
    Return new dimensions that will overwrite default dimensions for a
    varfull object with ``{}`` for a 1D plot.

    :param varfull_bracket: a ``varfull`` object with ``{}`` (e.g.,
        ``atmos_average.temp{lev=10;ls=350;lon=155;lat=25}``)
    :type varfull_bracket: str
    
    :param t_in: self.t variable
    :type t_in: array [time]
    
    :param lat_in: self.lat variable
    :type lat_in: array [lat]
    
    :param lon_in: self.lon variable
    :type lon_in: array [lon]
    
    :param lev_in: self.lev variable
    :type lev_in: array [lev]
    
    :param ftod_in: self.ftod variable
    :type ftod_in: array [tod]
    
    :return: ``varfull`` object without brackets (e.g.,
        ``atmos_average.temp``);
        :return: (t_out) dimension to update;
        :return: (lat_out) dimension to update;
        :return: (lon_out) dimension to update;
        :return: (lev_out) dimension to update;
        :return: (ftod_out) dimension to update;

    """
    # Initialization: Use dimension provided in template
    t_out = t_in
    lat_out = lat_in
    lon_out = lon_in
    lev_out = lev_in

    # Left of "{":
    varfull_no_bracket = varfull_bracket.split("{")[0].strip()

    # Right of "{", with last "}" removed:
    overwrite_txt = remove_whitespace(varfull_bracket.split("{")[1][:-1])

    # Count number of "=" in string
    ndim_update = overwrite_txt.count("=")

    # Split to different blocks (e.g., lat = 3. and lon = 20)
    split_dim = overwrite_txt.split(";")
    for i in range(0, ndim_update):
        # Check if requested dimension exists:
        if split_dim[i].split("=")[0] not in ["time", "lev", "lon", "lat",
                                              "tod"]:
            print(f"{Yellow}*** Warning*** ignoring dimension: "
                  f"{split_dim[i].split('=')[0]} because it is not recognized."
                  f"Valid dimensions = time, lev, lon, lat, or tod{Nclr}")
        if split_dim[i].split("=")[0] == "ls":
            t_out = filter_input(split_dim[i].split("=")[1], "float")
        if split_dim[i].split("=")[0] == "lat":
            lat_out = filter_input(split_dim[i].split("=")[1], "float")
        if split_dim[i].split("=")[0] == "lon":
            lon_out = filter_input(split_dim[i].split("=")[1], "float")
        if split_dim[i].split("=")[0] == "lev":
            lev_out = filter_input(split_dim[i].split("=")[1], "float")

        # Always get time of day
        ftod_out = None
        if split_dim[i].split("=")[0] == "tod":
            ftod_out = filter_input(split_dim[i].split("=")[1], "float")
    # .. note:: filter_input() converts text ("3" or "4,5") to variable:
    # (e.g., numpy.array([3.]) or numpy.array([4.,5.]))
    return varfull_no_bracket, t_out, lat_out, lon_out, lev_out, ftod_out


def create_exec(raw_input, varfull_list):
    expression_exec = raw_input
    for i in range(0, len(varfull_list)):
        swap_txt = f"[{varfull_list[i]}]"
        expression_exec = expression_exec.replace(swap_txt, f"VAR[{i:0}]")
    return expression_exec


def fig_layout(subID, nPan, vertical_page=False):
    """
    Return figure layout.

    :param subID: current subplot number
    :type subID: int
    
    :param nPan: number of panels desired on page (max = 64, 8x8)
    :type nPan: int
    
    :param vertical_page: reverse the tuple for portrait format if
        ``True``
    :type vertical_page: bool
    
    :return: plot layout (e.g., ``plt.subplot(nrows = out[0], ncols =
        out[1], plot_number = out[2])``)
    :rtype: tuple

    """
    out = list((0, 0, 0))

    if nPan == 1:
        # nrow, ncol
        layout = (1, 1)
    if nPan == 2:
        layout = (1, 2)
    if nPan == 3 or nPan == 4:
        layout = (2, 2)
    if nPan == 5 or nPan == 6:
        layout = (2, 3)
    if nPan == 7 or nPan == 8:
        layout = (2, 4)
    if nPan == 9:
        layout = (3, 3)
    if 10 <= nPan <= 12:
        layout = (3, 4)
    if 13 <= nPan <= 16:
        layout = (4, 4)
    if 17 <= nPan <= 20:
        layout = (4, 5)
    if 21 <= nPan <= 25:
        layout = (5, 5)
    if 26 <= nPan <= 30:
        layout = (5, 6)
    if 30 <= nPan <= 36:
        layout = (6, 6)
    if 36 <= nPan <= 42:
        layout = (7, 6)
    if 42 <= nPan <= 49:
        layout = (7, 7)
    if 49 <= nPan <= 56:
        layout = (8, 7)
    if 56 <= nPan <= 64:
        layout = (8, 8)
    if vertical_page:
        layout = layout[::-1]

    # Finally the current plot
    out[0:2] = layout
    out[2] = subID

    return out


def make_template():
    """
    Generate the ``Custom.in`` template file.
    
    :return: Custom.in blank template

    """
    global customFileIN  # Will be modified
    global current_version
    newname = os.path.join(output_path,"Custom.in")
    newname = create_name(newname)

    customFileIN = open(newname, "w")

    # Add line header. Primary use: change text color in VIM
    lh = "# "

    # Create header with instructions. Add version number to title.
    customFileIN.write(
        f"===================== |MarsPlot V{str(current_version)}| ===================\n")
    if args.trim_text is not None:
        # Additional instructions if requested
        customFileIN.write(
            "# ================================================= INSTRUCTIONS =================================================\n")
        customFileIN.write("# - Copy/paste template for the desired plot type. - Do not modify text left of an equal ``=`` sign. \n")
        customFileIN.write("# - Add comments using ``#``                         - Skip plots by setting <<<< Plot = False >>>> \n")
        customFileIN.write("# - Capitalize ``True``, ``False``, and ``None``.        - Do not use quotes ("") anywhere in this file. \n")
        customFileIN.write("# \n")
        customFileIN.write("# Group figures onto pages using``HOLD ON`` and ``HOLD OFF``. \n")
        customFileIN.write("# Optionally, use ``row,col`` to specify the layout: HOLD ON 2,3``. \n")
        customFileIN.write("# Use ``ADD LINE`` between 1D plots to overplot on the same figure. \n")
        customFileIN.write("# Figures templates must appear after ``START`` and before ``STOP``. \n")
        customFileIN.write("# Set the colorbar range with ``Cmin, Cmax``. Scientific notation (e.g., 1e-6, 2e3) is supported. \n")
        customFileIN.write("# Set the colorbar intervals directly by providing a list (e.g., 1e-6, 1e-4, 1e-2, 1e-0). \n")
        customFileIN.write("# Set the contour intervals for ``2nd Variable`` in a list (e.g., 150, 200, 250, 300, 350). \n")
        customFileIN.write("# The vertical grid of the *.nc file used in the plot determines what ``Level`` refers to.\n")
        customFileIN.write("#    ``Level`` can be: ``level``, ``pfull``, ``pstd``, ``plevs`` [Pa] or ``zstd``, ``zagl``, or ``zgrid`` [m].\n")
        customFileIN.write("# \n")
        customFileIN.write("# ============================================ ALGEBRA ============================================\n")
        customFileIN.write("# Use square brackets ``[]`` for element-wise operations: \n")
        customFileIN.write("#    ``[fixed.zsurf]/(10.**3)``            Convert between units ([m] to [km], in this case).\n")
        customFileIN.write("#    ``[file.var1]/[file.var2]*610``       Multiply variables together.\n")
        customFileIN.write("#    ``[file.var]-[file@2.var]``           Difference plot of ``var`` from 2 simulations.\n")
        customFileIN.write("#    ``[file.var]-[file.var{lev=10}]``     Difference plot of ``var`` at two levels.\n")
        customFileIN.write("# Square brackets support the following expressions: sqrt, log, exp, abs, min, max, & mean.\n")
        customFileIN.write("# \n")
        customFileIN.write("# ========================================= FREE DIMENSIONS =========================================\n")
        customFileIN.write("# Dimensions can be ``time``, ``lev``, ``lat``, ``lon``, or ``tod``.\n")
        customFileIN.write("# Dimensions default to None when a value or range is not specified. None corresponds to: \n")
        customFileIN.write("#    time = -1      The last (most recent) timestep (Nt).\n")
        customFileIN.write("#    lev  = sfc     Nz for *.nc files, 0 for *_pstd.nc files.\n")
        customFileIN.write("#    lat  = 0       Equator\n")
        customFileIN.write("#    lon  = ``all``   Zonal average over all longitudes\n")
        customFileIN.write("#    tod  = ``15``    3 PM UT \n")
        customFileIN.write("# Setting a dimension equal to a number finds the value closest to that number. \n")
        customFileIN.write("# Setting a dimension equal to ``all`` averages the dimension over all values. \n")
        customFileIN.write("# Setting a dimension equal to a range averages the dimension over the values in the range. \n")
        customFileIN.write("# You can also overwrite a dimension in the Main Variable input using curvy brackets ``{}`` and the\n")
        customFileIN.write("#    dimension name. Separate the arguments with semi-colons ``;`` \n")
        customFileIN.write("#        e.g., Main Variable = atmos_average.temp{ls = 90; lev= 5.,10; lon= all; lat=45} \n")
        customFileIN.write("#    Values must correspond to the units of the variable in the file: \n")
        customFileIN.write("#        time [Ls], lev [Pa/m], lon [+/-180°], and lat [°]. \n")
        customFileIN.write("# * You can only select a time of day (tod) in diurn files using this syntax: \n")
        customFileIN.write("#        e.g., Main Variable = atmos_diurn.ps{tod = 20} \n")
        customFileIN.write("# You can also specify the fontsize in Title using curvy brackets and ``size``:\n")
        customFileIN.write("#        e.g., Title = Temperature [K] {size = 20}.\n")
        customFileIN.write("# \n")
        customFileIN.write("# ==================================== TIME SERIES AND 1D PLOTS ====================================\n")
        customFileIN.write("# Set the X axis variable by indicating AXIS after the appropriate dimension: \n")
        customFileIN.write("#        e.g., Ls = AXIS \n")
        customFileIN.write("# The other dimensions remain FREE DIMENSIONS and accept values as described above. \n")
        customFileIN.write("# The ``Diurnal [hr]`` dimension only accepts ``AXIS`` or ``None``. Indicate time of day only using the``\n")
        customFileIN.write("#    ``tod`` syntax as described in FREE DIMENSIONS. \n")
        customFileIN.write("# \n")
        customFileIN.write("# ================================== AXIS OPTIONS AND PROJECTIONS ==================================\n")
        customFileIN.write("# Set the X and Y axis limits, map projection, colormap, and linestyle under Axis Options. \n")
        customFileIN.write("# All Matplolib styles are supported. \n")
        customFileIN.write("#    ``cmap``  colormap    ``jet`` (winds), ``nipy_spectral`` (temperature), ``bwr`` (diff plot), etc. \n")
        customFileIN.write("#    ``scale`` gradient    ``lin`` (linear), ``log`` (logarithmic; Cmin, Cmax is typically expected. \n")
        customFileIN.write("#    ``line``  linestyle   ``-r`` (solid red), ``--g`` (dashed green), ``-ob`` (solid blue + markers). \n")
        customFileIN.write("#    ``proj``  projection  Cylindrical: ``cart`` (Cartesian), ``robin`` (Robinson), ``moll`` (Mollweide), \n")
        customFileIN.write("#                        Azithumal: ``Npole lat`` (North Pole), ``Spole lat`` (South Pole),\n")
        customFileIN.write("#                        ``ortho lon,lat`` (Orthographic). \n")
        customFileIN.write("# \n")
        customFileIN.write("# ===================== FILES FROM MULTIPLE SIMULATIONS =====================\n")
        customFileIN.write("# Under <<< Simulations >>>, there are numbered lines (``N>``) for you to use to indicate the \n")
        customFileIN.write("#    path to the *.nc file you want to reference. Empty fields are ignored. \n")
        customFileIN.write("# Provide the FULL PATH on the line, e.g., ``2> /u/User/FV3/path/to/history``. \n")
        customFileIN.write("# Specify the *.nc file from which to plot using the ``@`` symbol + the simulation number:\n")
        customFileIN.write("#    in the call to Main Variable, e.g., Main Variable = atmos_average@2.temp \n")
        customFileIN.write("# \n")

    customFileIN.write(
        "<<<<<<<<<<<<<<<<<<<<<< Simulations >>>>>>>>>>>>>>>>>>>>>\n")
    customFileIN.write("ref> None\n")
    customFileIN.write("2> \n")
    customFileIN.write("3>\n")
    customFileIN.write(
        "=======================================================\n")
    customFileIN.write("START\n\n")

    # For the default list of figures in main(), create a template.
    for i in range(0, len(objectList)):
        if objectList[i].subID == 1 and objectList[i].nPan > 1:
            customFileIN.write("HOLD ON\n")
        objectList[i].make_template()
        customFileIN.write("\n")
        if (objectList[i].nPan > 1 and
            objectList[i].subID == objectList[i].nPan):
            customFileIN.write("HOLD OFF\n")

        # Separate empty templates
        if i == 1:
            customFileIN.write("#=========================================================================\n")
            customFileIN.write("#================== Empty Templates (set to False)========================\n")
            customFileIN.write("#========================================================================= \n")
            customFileIN.write(" \n")
    customFileIN.close()

    # NAS system only: set group permissions & print confirmation
    give_permission(newname)
    print(f"{newname} was created")


def give_permission(filename):
    """
    Sets group permissions for files created on NAS.

    :param filename: name of the file
    :type filename: str

    """
    # NAS system only: set group permissions to file
    try:
        # Catch error and standard output
        subprocess.check_call(["setfacl -v"], shell = True,
                              stdout = open(os.devnull, "w"),
                              stderr = open(os.devnull, "w"))
        cmd_txt = f"setfacl -R -m g:s0846:r {filename}"
        subprocess.call(cmd_txt, shell = True)
    except subprocess.CalledProcessError:
        pass


def namelist_parser(Custom_file):
    """
    Parse a ``Custom.in`` template.

    :param Custom_file: full path to ``Custom.in`` file
    :type Custom_file: str

    :return: updated global variables, ``FigLayout``, ``objectList``

    """
    global objectList
    global customFileIN
    global input_paths

    # Custom.in file provided, flush default figures in main()
    objectList = []     # All individual plots
    panelList = []      # List of panels
    subplotList = []    # Layout of figures
    addLineList = []    # Add several lines to plot on same graph
    layoutList = []
    nobj = 0            # Number for object. "nobj = 1,[2,3],4" means
                        # plots 2 & 3 are 2-panel plot.
    npanel = 1          # Number of panels plotted along this object.
                        # e.g., npanel = 1 = object 1,
                        # e.g., npanel = 2 = objects 2 & 3
    subplotID = 1       # Subplot ID for each object. = 1 = object 1
    holding = False
    addLine = False
    addedLines = 0     # Line plots
    # Plot number at start of a new page (HOLD ON). Used if layout
    # provided with HOLD ON (e.g., HOLD ON 2,3)
    npage = 0
    layout = None
    customFileIN = open(Custom_file, "r")

    # Get version number in header
    version = float(
        customFileIN.readline().split("|")[1].strip().split("V")[1].strip()
    )

    if int(version) != int(current_version):
        # Check if main versions are compatible
        # 1.1 and 1.2 = compatible
        # 1.0 and 2.0 = NOT compatible
        print(f"{Yellow}*** Warning ***\nUsing MarsPlot V{current_version} "
              f"but Custom.in template is deprecated (V{version})"
              f"\n***************{Nclr}")

    while (customFileIN.readline()[0] != "<"):
        # Skip the header
        pass

    while True:
        # Read paths under ``<<<<<<<<<< Simulations >>>>>>>>>>>``
        line = customFileIN.readline()
        if line[0] == "#":
            # Skip comments
            pass
        else:
            if line[0] == "=":
                # Finished reading
                break
            if line.split(">")[0] == "ref":
                # Special case: use reference simulation
                if line.split(">")[1].strip() != "None":
                    # If it is different from default, overwrite it
                    input_paths[0] = line.split(">")[1].strip()
            else:
                if ">" in line:
                    # Line contains ">" symbol
                    if line.split(">")[1].strip():
                        # Line exists and is not blank
                        input_paths.append(line.split(">")[1].strip())

    # Skip lines until keyword START found. Initialize backup counter.
    nsafe = 0
    while True and nsafe < 2000:
        line = customFileIN.readline()
        if line.strip() == "START":
            break
        nsafe += 1
    if nsafe == 2000:
        print(f"{Red}Custom.in is missing a 'START' keyword after the '====='"
              f"simulation block{Nclr}")

    # Start reading figure templates
    while True:
        line = customFileIN.readline()
        if not line or line.strip() == "STOP":
            # Reached end of file
            break
        if line.strip()[0:7] == "HOLD ON":
            holding = True
            subplotID = 1
            # Get layout info
            if "," in line:
                # Layout provided (e.g., HOLD ON 2,3), return as string
                tmp = line.split("ON")[-1].strip()
                # Returns [2,3]
                layout = [int(tmp.split(",")[0]), int(tmp.split(",")[1])]
            else:
                layout = None
        if line.strip() == "ADD LINE":
            # Overplot 1D plot
            addLine = True
        if line[0] == "<":
            # If new figure
            figtype, boolPlot = get_figure_header(line)
            if boolPlot:
                # Only if we want to plot field
                # Add object to list
                if figtype == "Plot 2D lon X lat":
                    objectList.append(Fig_2D_lon_lat())
                if figtype == "Plot 2D time X lat":
                    objectList.append(Fig_2D_time_lat())
                if figtype == "Plot 2D lat X lev":
                    objectList.append(Fig_2D_lat_lev())
                if figtype == "Plot 2D lon X lev":
                    objectList.append(Fig_2D_lon_lev())
                if figtype == "Plot 2D time X lev":
                    objectList.append(Fig_2D_time_lev())
                if figtype == "Plot 2D lon X time":
                    objectList.append(Fig_2D_lon_time())
                if figtype == "Plot 1D":
                    objectList.append(Fig_1D())
                objectList[nobj].read_template()
                nobj += 1
                if holding and not addLine:
                    subplotList.append(subplotID)
                    panelList.append(subplotID)
                    subplotID += 1
                    for iobj in range(npage, nobj-1):
                        # Add +1 panel to all plots on current page
                        panelList[iobj] += 1
                elif holding and addLine:
                    # Do not update subplot ID if adding lines
                    subplotList.append(subplotID-1)
                    panelList.append(subplotID-1)
                else:
                    # One plot per page. Reset page counter.
                    panelList.append(1)
                    subplotList.append(1)
                    npage = nobj
                    layout = None
                if layout:
                    layoutList.append(layout)
                else:
                    layoutList.append(None)
                if addLine:
                    addedLines += 1
                    addLineList.append(addedLines)
                else:
                    # No added lines
                    addLineList.append(0)
                    # Reset line counter
                    addedLines = 0

            # Reset after reading each block
            addLine = False
        if line.strip() == "HOLD OFF":
            holding = False
            subplotID = 1
            npage = nobj

    if holding:
        print(f"{Red}*** Error ***\nMissing ``HOLD OFF`` statement in "
              f"{Custom_file}{Nclr}")
        exit()
    if addLine:
        print(f"{Red}*** Error ***\nCannot have ``ADD LINE`` after the last "
              f"figure in {Custom_file}{Nclr}")
        exit()

    for i in range(0, nobj):
        # Distribute number of figures and panels for each plot
        objectList[i].subID = subplotList[i]
        objectList[i].nPan = panelList[i]
        objectList[i].addLine = addLineList[i]
        objectList[i].layout = layoutList[i]
    customFileIN.close()


def get_figure_header(line_txt):
    """
    Returns the plot type by confirming that template = ``True``.

    :param line_txt: template header from Custom.in (e.g.,
        ``<<<<<<<<<| Plot 2D lon X lat = True |>>>>>>>>``)
    :type line_txt: str
    
    :return: (figtype) figure type (e.g., ``Plot 2D lon X lat``)
    :rtype: str
    
    :return: (boolPlot) whether to plot (``True``) or skip (``False``)
        figure
    :rtype: bool

    """
    # Plot 2D lon X lat = True
    line_cmd = line_txt.split("|")[1].strip()
    # Plot 2D lon X lat
    figtype = line_cmd.split("=")[0].strip()
    # Return True
    boolPlot = line_cmd.split("=")[1].strip() == "True"
    return figtype, boolPlot


def format_lon_lat(lon_lat, type):
    """
    Format latitude and longitude as labels (e.g., 30°S, 30°N, 45°W,
    45°E)

    :param lon_lat: latitude or longitude (+180/-180)
    :type lon_lat: float
    
    :param type: ``lat`` or ``lon``
    :type type: str
    
    :return: formatted label
    :rtype: str

    """
    letter = ""
    if type == "lon":
        if lon_lat < 0:
            letter = "W"
        if lon_lat > 0:
            letter = "E"
    elif type == "lat":
        if lon_lat < 0:
            letter = "S"
        if lon_lat > 0:
            letter = "N"

    # Remove minus sign, if any
    lon_lat = abs(lon_lat)
    return f"{lon_lat}{letter}"

# ======================================================================
#                       FILE SYSTEM UTILITIES
# ======================================================================
def get_Ncdf_num():
    """
    Return the prefix numbers for the netCDF files in the directory.
    Requires at least one ``fixed`` file in the directory.

    :return: a sorted array of sols
    :rtype: array

    """
    # e.g., 00350.fixed.nc
    list_dir = os.listdir(input_paths[0])
    avail_fixed = [k for k in list_dir if ".fixed.nc" in k]
    # Remove .fixed.nc (returning 00350 or 00000)
    list_num = [item[0:5] for item in avail_fixed]
    # Transform to array (returning [0, 350])
    Ncdf_num = np.sort(np.asarray(list_num).astype(float))
    if Ncdf_num.size == 0:
        Ncdf_num = None
    return Ncdf_num


def select_range(Ncdf_num, bound):
    """
    Return the prefix numbers for the netCDF files in the directory
    within the user-defined range.

    :param Ncdf_num: a sorted array of sols
    :type Ncdf_num: array
    
    :param bound: a sol (e.g., 0350) or range of sols ``[min max]``
    :type bound: int or array
    
    :return: a sorted array of sols within the bounds
    :rtype: array

    """
    bound = np.array(bound)
    if bound.size == 1:
        Ncdf_num = Ncdf_num[Ncdf_num == bound]
        if Ncdf_num.size == 0:
            print(f"{Red}*** Error *** \n"
                  f"File {int(bound):05}.fixed.nc not found{Nclr}")
            exit()
    elif bound.size == 2:
        Ncdf_num = Ncdf_num[Ncdf_num >= bound[0]]
        Ncdf_num = Ncdf_num[Ncdf_num <= bound[1]]
        if Ncdf_num.size == 0:
            print(f"{Red}*** Error ***\nNo fixed file with date between "
                  f"[{int(bound[0]):05}-{int(bound[1]):05}] detected. Please "
                  f"double check the range.{Nclr}")
            exit()
    return Ncdf_num


def create_name(root_name):
    """
    Modify file name if a file with that name already exists.

    :param root_name: path + default name for the file type (e.g.,
        ``/path/custom.in`` or ``/path/figure.png``)
    :type root_name: str
    
    :return: the modified name if the file already exists
        (e.g., ``/path/custom_01.in`` or ``/path/figure_01.png``)
    :rtype: str

    """
    n = 1
    # Get extension length (e.g., 2 for *.nc, 3 for *.png)
    len_ext = len(root_name.split(".")[-1])
    ext = root_name[-len_ext:]
    new_name = root_name

    if os.path.isfile(new_name):
        # If example.png already exists, create example_01.png
        new_name = f"{root_name[0:-(len_ext + 1)]}_{n:02}.{ext}"

    while os.path.isfile(f"{root_name[0:-(len_ext + 1)]}_{n:02}.{ext}"):
        # If example_01.png already exists, create example_02.png etc
        n = n + 1
        new_name = f"{root_name[0:-(len_ext + 1)]}_{n:02}.{ext}"
    return new_name

def progress(k, Nmax, txt="", success=True):
    """
    Display a progress bar when performing heavy calculations.

    :param k: current iteration of the outer loop
    :type k: float
    
    :param Nmax: max iteration of the outer loop
    :type Nmax: float
    
    :return: progress bar (EX: ``Running... [#---------] 10.64 %``)

    """
    progress = float(k)/Nmax
    barLength = 10
    block = int(round(barLength*progress))
    bar = f"[{('#'*block ) + ('-'*(barLength-block))}]"
    if success == True:
        status = f"{int(100*progress):>3} % {Green}({txt}){Nclr}"
    elif success == False:
        status = f"{int(100*progress):>3} % {Red}({txt}){Nclr}"
    elif success == None:
        status = f"{int(100*progress):>3} % ({txt})"
    text = (f"\r{bar}{status}\n")
    sys.stdout.write(text)
    if not debug:
        sys.stdout.flush()


def prep_file(var_name, file_type, simuID, sol_array):
    """
    Open the file as a Dataset or MFDataset object depending on its
        status on Lou. Note that the input arguments are typically
        extracted from a ``varfull`` object (e.g.,
        ``03340.atmos_average.ucomp``) and not from a file whose disk
        status is known beforehand.

    :param var_name: variable to extract (e.g., ``ucomp``)
    :type var_name: str
    
    :param file_type: MGCM output file type (e.g., ``average``)
    :type file_name: str
    
    :param simuID: simulation ID number (e.g., 2 for 2nd simulation)
    :type simuID: int
    
    :param sol_array: date in file name (e.g., [3340,4008])
    :type sol_array: list

    :return: Dataset or MFDataset object;
        (var_info) longname and units;
        (dim_info) dimensions e.g., (``time``, ``lat``,``lon``);
        (dims) shape of the array e.g., [133,48,96]

    """
    global input_paths
    # Holds sol numbers (e.g., [1500,2400])
    global Ncdf_num
    # Specific sol requested (e.g., [2400])
    Sol_num_current = [0]

    if os.path.isfile(os.path.join(f"{input_paths[simuID]}",f"{file_type}.nc")):
        # First check if file on tape without a sol number
        # (e.g., Luca_dust_MY24_dust.nc exists on disk)
        file_has_sol_number = False
    else:
        # If file does NOT exist, append sol number in MarsPlot
        # (e.g., Custom.in -d XXXX) or of last file in dir.
        file_has_sol_number = True
        if sol_array != [None]:
            # File number explicitly provided in varfull
            # (e.g., 01336.atmos_average.nc)
            Sol_num_current = sol_array
        elif Ncdf_num is not None:
            # File number NOT provided in varfull
            Sol_num_current = Ncdf_num

    # Create list of files (even if only one file provided)
    nfiles = len(Sol_num_current)
    file_list = [None]*nfiles

    # Loop over requested timesteps
    for i in range(0, nfiles):
        if file_has_sol_number:
            # Sol number
            file_list[i] = (os.path.join(input_paths[simuID],
                            f"{int(Sol_num_current[i]):05}.{file_type}.nc"))
        else:
            # No sol number
            file_list[i] = os.path.join(input_paths[simuID],f"{file_type}.nc")

        check_file_tape(file_list[i])

    try:
        # Files on tape, open with MFDataset if aggregate dim detected
        f = MFDataset(file_list, "r")
    except IOError:
        # IOError should be: "master dataset ***.nc does not have
        # an aggregation dim". Use Dataset otherwise
        f = Dataset(file_list[0], "r")

    var_info = (f"{getattr(f.variables[var_name], 'long_name', '')} "
                f"[{getattr(f.variables[var_name], 'units', '')}]")
    dim_info = f.variables[var_name].dimensions
    dims = f.variables[var_name].shape
    return f, var_info, dim_info, dims

class CustomTicker(LogFormatterSciNotation):
    def __call__(self, x, pos=None):
        if x < 0:
            return LogFormatterSciNotation.__call__(self, x, pos = None)
        else:
            return "{x:g}".format(x = x)

# ======================================================================
#                           FIGURE DEFINITIONS
# ======================================================================
class Fig_2D(object):
    def __init__(self, varfull="fileYYY.XXX", doPlot=False, varfull2=None):

        self.title = None
        self.varfull = varfull
        self.range = None
        self.fdim1 = None
        self.fdim2 = None
        self.ftod = None # Time of day
        self.varfull2 = varfull2
        self.contour2 = None
        # Logic
        self.doPlot = doPlot
        self.plot_type = self.__class__.__name__[4:]

        # Extract filetype, variable, and simulation ID (initialized
        # only for the default plots). Note that varfull objects
        # for default plots are simple (e.g., atmos_average.ucomp)
        (self.sol_array,
         self.filetype,
         self.var,
         self.simuID) = split_varfull(self.varfull)

        if self.varfull2:
            (self.sol_array2,
             self.filetype2,
             self.var2,
             self.simuID2) = split_varfull(self.varfull2)

        # Multipanel
        self.nPan = 1
        self.subID = 1
        self.layout = None  # e.g., [2,3], used only if HOLD ON 2,3

        # Annotation for free dimensions
        self.fdim_txt = ""
        self.success = False
        self.addLine = False
        self.vert_unit = "" # m or Pa

        # Axis options
        self.Xlim = None
        self.Ylim = None
        self.axis_opt1 = "jet"
        self.axis_opt2 = "lin" # Linear or logscale
        self.axis_opt3 = None  # Placeholder for projection type


    def make_template(self, plot_txt, fdim1_txt, fdim2_txt, Xaxis_txt,
                      Yaxis_txt):
        customFileIN.write(
            f"<<<<<<<<<<<<<<| {plot_txt:<15} = {self.doPlot} |>>>>>>>>>>>>>\n")
        customFileIN.write(f"Title          = {self.title}\n")          # 1
        customFileIN.write(f"Main Variable  = {self.varfull}\n")        # 2
        customFileIN.write(f"Cmin, Cmax     = {self.range}\n")          # 3
        customFileIN.write(f"{fdim1_txt:<15}= {self.fdim1}\n")          # 4
        customFileIN.write(f"{fdim2_txt:<15}= {self.fdim2}\n")          # 4
        customFileIN.write(f"2nd Variable   = {self.varfull2}\n")       # 6
        customFileIN.write(f"Contours Var 2 = {self.contour2}\n")       # 7

        # Write colormap AND projection if plot is 2D_lon_lat (Line # 8)
        if self.plot_type == "2D_lon_lat":
            customFileIN.write(
                f"Axis Options  : {Xaxis_txt} = [None,None] | {Yaxis_txt} = "
                f"[None,None] | cmap = jet | scale = lin | proj = cart \n")
        else:
            #Special case of Xaxis_txt is Ls, the axis are set using the sol array
            #This is useful in the case multiple years are displayed
            if Xaxis_txt =='Ls': Xaxis_txt='Sol'
            customFileIN.write(
                f"Axis Options  : {Xaxis_txt} = [None,None] | {Yaxis_txt} = "
                f"[None,None] | cmap = jet |scale = lin \n")


    def read_template(self):
        self.title = rT("char")  # 1
        self.varfull= rT("char")  # 2
        self.range = rT("float") # 3
        self.fdim1 = rT("float") # 4
        self.fdim2 = rT("float") # 5
        self.varfull2 = rT("char")  # 6
        self.contour2 = rT("float") # 7
        (self.Xlim,
         self.Ylim,
         self.axis_opt1,
         self.axis_opt2,
         self.axis_opt3) = read_axis_options(customFileIN.readline()) # 8

        # Various sanity checks
        if self.range and len(np.atleast_1d(self.range)) == 1:
            print(f"{Yellow}*** Warning *** In plot {self.varfull}, Cmin, "
                  f"Cmax must be two values. Resetting to default{Nclr}")
            self.range = None


    def data_loader_2D(self, varfull, plot_type):
        if not "[" in varfull:
            # Plot 1 of the variables in the file
            if "{" in varfull:
                # If overwriting dim, get new dim and trim varfull from
                # {lev=5.}
                (varfull, fdim1_extract,
                 fdim2_extract, ftod_extract) = get_overwrite_dim_2D(varfull,
                    plot_type, self.fdim1, self.fdim2, self.ftod)
            else:
                # If no "{}" in varfull, do not overwrite dims, use
                # plot defaults
                fdim1_extract = self.fdim1
                fdim2_extract = self.fdim2
                ftod_extract = self.ftod

            sol_array, filetype, var, simuID = split_varfull(varfull)
            xdata, ydata, var, var_info = self.read_NCDF_2D(var, filetype,
                                                            simuID, sol_array,
                                                            plot_type,
                                                            fdim1_extract,
                                                            fdim2_extract,
                                                            ftod_extract)
        else:
            # Recognize operation on variables
            VAR = []
            # Extract individual variables and prepare for execution
            varfull = remove_whitespace(varfull)
            varfull_list = get_list_varfull(varfull)
            # Initialize list of requested dims
            fdim1_list = [None]*len(varfull_list)
            fdim2_list = [None]*len(varfull_list)
            ftod_list = [None]*len(varfull_list)
            expression_exec = create_exec(varfull, varfull_list)

            for i in range(0, len(varfull_list)):
                if "{" in varfull_list[i]:
                    # If overwriting dim, get new dim and trim varfull
                    # from {lev=5.}
                    (varfull_list[i], fdim1_list[i],
                     fdim2_list[i],ftod_list[i]) = get_overwrite_dim_2D(
                         varfull_list[i],plot_type, self.fdim1,
                         self.fdim2, self.ftod)
                else:
                    # If no "{}" in varfull, do not overwrite dims, use
                    # plot defaults
                    fdim1_list[i] = self.fdim1
                    fdim2_list[i] = self.fdim2
                    ftod_list[i] = self.ftod

                sol_array, filetype, var, simuID = split_varfull(
                    varfull_list[i])
                xdata, ydata, temp, var_info = self.read_NCDF_2D(
                    var, filetype, simuID, sol_array, plot_type,
                    fdim1_list[i], fdim2_list[i], ftod_list[i]
                )

                VAR.append(temp)
            var_info = varfull
            var = eval(expression_exec)#TODO removed ,namespace

        return xdata, ydata, var, var_info


    def read_NCDF_2D(self, var_name, file_type, simuID, sol_array, plot_type,
                     fdim1, fdim2, ftod):
        f, var_info, dim_info, dims = prep_file(var_name, file_type,
                                                simuID, sol_array)

        # Get file type (fixed, diurn, average, daily) and interp type
        # (pfull, zstd, etc.)
        f_type, interp_type = FV3_file_type(f)

        # Initialize dims (in all .nc files)
        lat = f.variables["lat"][:]
        lati = np.arange(0, len(lat))
        lon = f.variables["lon"][:]
        loni = np.arange(0, len(lon))

        # If self.fdim is empty, add variable name (do only once)
        add_fdim = False
        if not self.fdim_txt.strip():
            add_fdim = True

        # ------------------------ Time of Day -------------------------
        # For diurn files, select data on time of day axis and update
        # dims so that resulting variable is same as in average and
        # daily files. Time of day always 2nd dim (dim_info[1])

        if f_type == "diurn" and dim_info[1][:11] == "time_of_day":
            tod = f.variables[dim_info[1]][:]
            todi, temp_txt = get_tod_index(ftod, tod)
            # Update dim_info
            # time, time_of_day_XX, lat, lon -> time, lat, lon
            # OR
            # time, time_of_day_XX, pfull, lat, lon -> time, pfull, lat, lon
            dim_info = (dim_info[0],) + dim_info[2:]

            if add_fdim:
                self.fdim_txt += temp_txt
        # --------------------------------------------------------------

        # Load variable depending on requested free dims
        # ====== static ======= ignore level and time dim
        if dim_info == ("lat", "lon"):
            var = f.variables[var_name][lati, loni]
            f.close()
            return lon, lat, var, var_info

        # ====== time, lat, lon =======
        if dim_info == ("time", "lat", "lon"):
            # Initialize dimension
            t = f.variables["time"][:]
            LsDay = np.squeeze(f.variables["areo"][:])
            ti = np.arange(0, len(t))
            # For diurn file, change time_of_day[time, 24, 1] ->
            # time_of_day[time] at midnight UT
            if f_type == "diurn" and len(LsDay.shape) > 1:
                LsDay = np.squeeze(LsDay[:, 0])
            # Stack time and areo array as one variable
            t_stack = np.vstack((t, LsDay))

            if plot_type == "2D_lon_lat":
                ti, temp_txt = get_time_index(fdim1, LsDay)
            if plot_type == "2D_time_lat":
                loni, temp_txt = get_lon_index(fdim1, lon)
            if plot_type == "2D_lon_time":
                lati, temp_txt = get_lat_index(fdim1, lat)

            if add_fdim:
                self.fdim_txt += temp_txt

            # Extract data and close file
            if f_type == "diurn":
                # Do time of day average first
                var = f.variables[var_name][ti, todi, lati, loni].reshape(
                    len(np.atleast_1d(ti)),
                    len(np.atleast_1d(todi)),
                    len(np.atleast_1d(lati)),
                    len(np.atleast_1d(loni))
                )
                var = mean_func(var, axis = 1)
            else:
                var = f.variables[var_name][ti, lati, loni].reshape(
                    len(np.atleast_1d(ti)),
                    len(np.atleast_1d(lati)),
                    len(np.atleast_1d(loni))
                )
            f.close()
            w = area_weights_deg(var.shape, lat[lati])

            # Return data
            if plot_type == "2D_lon_lat":
                # Time average
                return lon, lat, mean_func(var, axis = 0), var_info
            if plot_type == "2D_time_lat":
                # Transpose, X dim must be in last column of variable
                return t_stack, lat, mean_func(var, axis = 2).T, var_info
            if plot_type == "2D_lon_time":
                return (lon, t_stack, np.average(var, weights = w, axis = 1),
                        var_info)

        # ====== [time, lev, lat, lon] =======
        if (dim_info == ("time", "pfull", "lat", "lon")
           or dim_info == ("time", "level", "lat", "lon")
           or dim_info == ("time", "pstd",  "lat", "lon")
           or dim_info == ("time", "zstd",  "lat", "lon")
           or dim_info == ("time", "zagl",  "lat", "lon")
           or dim_info == ("time", "zgrid", "lat", "lon")
           or dim_info == ("zgrid", "lat",  "lon")):

            if dim_info[1] in ["pfull", "level", "pstd"]:
                self.vert_unit = "Pa"
            if dim_info[1] in ["zagl", "zstd", "zgrid"]:
                self.vert_unit = "m"
            if dim_info[0] in ["zgrid"]:
                # Thermal inertia = special case
                self.vert_unit = "m"

            # Initialize dims
            levs = f.variables[dim_info[1]][:]
            zi = np.arange(0, len(levs))
            t = f.variables["time"][:]
            LsDay = np.squeeze(f.variables["areo"][:])
            ti = np.arange(0, len(t))
            # For diurn file, change time_of_day[time, 24, 1] ->
            # time_of_day[time] at midnight UT
            if f_type == "diurn" and len(LsDay.shape) > 1:
                LsDay = np.squeeze(LsDay[:, 0])
            # Stack time and areo arrays as 1 variable
            t_stack = np.vstack((t, LsDay))

            if plot_type == "2D_lon_lat":
                ti, temp_txt = get_time_index(fdim1, LsDay)
                if add_fdim:
                    self.fdim_txt += temp_txt
                zi, temp_txt = get_level_index(fdim2, levs)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if plot_type == "2D_time_lat":
                loni, temp_txt = get_lon_index(fdim1, lon)
                if add_fdim:
                    self.fdim_txt += temp_txt
                zi, temp_txt = get_level_index(fdim2, levs)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if plot_type == "2D_lat_lev":
                ti, temp_txt = get_time_index(fdim1, LsDay)
                if add_fdim:
                    self.fdim_txt += temp_txt
                loni, temp_txt = get_lon_index(fdim2, lon)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if plot_type == "2D_lon_lev":
                ti, temp_txt = get_time_index(fdim1, LsDay)
                if add_fdim:
                    self.fdim_txt += temp_txt
                lati, temp_txt = get_lat_index(fdim2, lat)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if plot_type == "2D_time_lev":
                lati, temp_txt = get_lat_index(fdim1, lat)
                if add_fdim:
                    self.fdim_txt += temp_txt
                loni, temp_txt = get_lon_index(fdim2, lon)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if plot_type == "2D_lon_time":
                lati, temp_txt = get_lat_index(fdim1, lat)
                if add_fdim:
                    self.fdim_txt += temp_txt
                zi, temp_txt = get_level_index(fdim2, levs)
                if add_fdim:
                    self.fdim_txt += temp_txt

            if f_type == "diurn":
                # time of day average
                var = f.variables[var_name][ti, todi, zi, lati, loni].reshape(
                    len(np.atleast_1d(ti)),
                    len(np.atleast_1d(todi)),
                    len(np.atleast_1d(zi)),
                    len(np.atleast_1d(lati)),
                    len(np.atleast_1d(loni))
                )
                var = mean_func(var, axis = 1)

            else:
                var = f.variables[var_name][ti, zi, lati, loni].reshape(
                    len(np.atleast_1d(ti)),
                    len(np.atleast_1d(zi)),
                    len(np.atleast_1d(lati)),
                    len(np.atleast_1d(loni))
                )

            f.close()
            w = area_weights_deg(var.shape, lat[lati])

            if plot_type == "2D_lon_lat":
                return (lon, lat,
                        mean_func(mean_func(var, axis = 1), axis = 0),
                        var_info)
            if plot_type == "2D_time_lat":
                # Transpose
                return (t_stack, lat,
                        mean_func(mean_func(var, axis = 1), axis = 2).T,
                        var_info)
            if plot_type == "2D_lat_lev":
                return (lat, levs,
                        mean_func(mean_func(var, axis = 3), axis = 0),
                        var_info)
            if plot_type == "2D_lon_lev":
                return (lon, levs,
                        mean_func(np.average(var,weights = w, axis = 2),
                                    axis = 0),
                        var_info)
            if plot_type == "2D_time_lev":
                # Transpose
                return (t_stack, levs,
                        mean_func(np.average(var, weights = w, axis = 2),
                                    axis = 2).T,
                        var_info)
            if plot_type == "2D_lon_time":
                return (lon,t_stack,
                        mean_func(np.average(var, weights = w, axis = 2),
                                    axis = 1),
                        var_info)


    def plot_dimensions(self):
        print(f"{Yellow}{self.ax.get_position()}{Nclr}")


    def make_title(self, var_info, xlabel, ylabel):
        if self.title:
            # Title provided
            if "{fontsize = " in self.title:
                # Fontsize specified
                fs = int(remove_whitespace(
                    (self.title).split("{fontsize = ")[1].split("}")[0]))
                title_text = ((self.title).split("{fontsize = ")[0])
                plt.title(title_text,
                          fontsize = (fs - self.nPan*title_factor),
                          wrap = False)
            else:
                # Fontsize not specified
                plt.title(self.title,
                          fontsize = title_size - self.nPan*title_factor)
        else:
            # Title NOT provided
            plt.title(f"{var_info}\n{self.fdim_txt[1:]}",
                      fontsize = (title_size - self.nPan*title_factor),
                      wrap = False)

        plt.xlabel(xlabel, fontsize = (label_size - self.nPan*label_factor))
        plt.ylabel(ylabel, fontsize = (label_size - self.nPan*label_factor))


    def make_colorbar(self, levs):
        if self.axis_opt2 == "log":
            formatter = LogFormatter(10, labelOnlyBase = False)
            if self.range:
                cbar = plt.colorbar(ticks = levs,
                                    orientation = "horizontal",
                                    aspect = 30,
                                    format = formatter)
            else:
                cbar = plt.colorbar(orientation = "horizontal",
                                    aspect = 30,
                                    format = formatter)

        else:
            cbar = plt.colorbar(orientation = "horizontal", aspect = 30)

        # Shrink colorbar label as number of subplots increases
        cbar.ax.tick_params(labelsize=(label_size - self.nPan*label_factor))


    def return_norm_levs(self):
        norm = None
        levs = None
        if self.axis_opt2 == "log":
            # Logarithmic colormap
            norm = LogNorm()
        else:
            # Linear colormap (default)
            self.axis_opt2 = "lin"
            norm = None
        if self.range:
            if self.axis_opt2 == "lin":
                # If 2 numbers provided (e.g., Cmin,Cmax)
                if len(self.range) == 2:
                    levs = np.linspace(self.range[0], self.range[1], levels)
                # If list provided, set intervals explicitly
                else:
                    levs = self.range

            if self.axis_opt2 == "log":
                if self.range[0] <= 0 or self.range[1] <= 0:
                    print(f"{Red}*** Error using log scale, bounds cannot be "
                          f"zero or negative{Nclr}")
                levs = np.logspace(
                    np.log10(self.range[0]), np.log10(self.range[1]), levels)
        return norm, levs


    def exception_handler(self, e, ax):
        if debug:
            raise
        sys.stdout.write("\033[F")
        # Cursor up one line, then clear lines previous output
        sys.stdout.write("\033[K")
        print(f"{Yellow}*** Warning *** {e}{Nclr}")
        ax.text(0.5, 0.5, f"ERROR: {e}",
            horizontalalignment = "center",
            verticalalignment = "center",
            bbox = dict(boxstyle="round",
                        ec = (1., 0.5, 0.5),
                        fc = (1., 0.8, 0.8),),
            transform = ax.transAxes, wrap = True, fontsize = 16)


    def fig_init(self):
        # Create figure
        if self.layout is None:
            # No layout specified
            out = fig_layout(self.subID, self.nPan, vertical_page)
        else:
            # Layout specified
            out = np.append(self.layout, self.subID)
        if self.subID == 1:
            # Create figure if 1st panel
            # 1.4 is ratio (16:9 screen would be 1.77)
            fig = plt.figure(facecolor="white", figsize = (width_inch,
                                                           height_inch))
        ax = plt.subplot(out[0], out[1], out[2]) # nrow, ncol, subID
        # Nans = grey
        ax.patch.set_color(".1")
        return ax


    def fig_save(self):
        # Save figure
        if self.subID == self.nPan:
            # Last subplot
            if self.subID == 1:
                # 1 plot
                if not "[" in self.varfull:
                    # Add split "{" in case varfull contains layer.
                    # Does not do anything else.
                    sensitive_name = self.varfull.split("{")[0].strip()
                    # If varfull = complex expression
                else:
                    expr = (get_list_varfull(
                        self.varfull)[0].split('{')[0].strip())
                    sensitive_name = (f"expression_{expr}")
            else:
                # Multipanel
                sensitive_name = "multi_panel"

            plt.tight_layout()
            self.fig_name = (os.path.join(output_path,"plots",
                             f"{sensitive_name}.{out_format}"))
            self.fig_name = create_name(self.fig_name)
            plt.savefig(self.fig_name, dpi=my_dpi)
            if out_format != "pdf":
                print(f"Saved:{self.fig_name}")


    def filled_contour(self, xdata, ydata, var):
        cmap = self.axis_opt1
        if cmap == "wbr":
            cmap = wbr_cmap()
        if cmap == "rjw":
            cmap = rjw_cmap()
        if cmap == "dkass_temp":
            cmap = dkass_temp_cmap()
        if cmap == "dkass_dust":
            cmap = dkass_dust_cmap()
        if cmap == "hot_cold":
            cmap = hot_cold_cmap()
        norm, levs = self.return_norm_levs()

        if self.range:
            plt.contourf(xdata, ydata, var, levs,
                         extend = "both", cmap = cmap, norm = norm)
        else:
            plt.contourf(xdata, ydata, var, levels, cmap = cmap, norm = norm)

        self.make_colorbar(levs)


    def solid_contour(self, xdata, ydata, var, contours):
        # Prevent error message when drawing contours
        np.seterr(divide="ignore", invalid="ignore")
        if contours is None:
            CS = plt.contour(xdata, ydata, var, 11, colors = "k",
                             linewidths = 2)
        else:
            # If one contour provided (as float), convert to array
            if type(contours) == float:
                contours = [contours]
            CS = plt.contour(xdata, ydata, var, contours,
                             colors = "k", linewidths = 2)
        plt.clabel(CS, inline = 1, fontsize = 14, fmt = "%g")


# ===============================

class Fig_2D_lon_lat(Fig_2D):


    # Make_template calls method from parent class
    def make_template(self):
        super(Fig_2D_lon_lat, self).make_template(
            "Plot 2D lon X lat", "Ls 0-360", "Level Pa/m", "Lon", "Lat")


    def get_topo_2D(self, varfull, plot_type):
        """
        This function returns the longitude, latitude, and topography
        to overlay as contours in a ``2D_lon_lat`` plot. Because the
        main variable requested may be complex
        (e.g., ``[01336.atmos_average_psdt2.temp]/1000.``), we will
        ensure to load the matching topography (here ``01336.fixed.nc``
        from the 2nd simulation). This function essentially does a
        simple task in a complicated way. Note that a great deal of
        the code is borrowed from the ``data_loader_2D()`` function.

        :param varfull: variable input to main_variable in Custom.in
            (e.g., ``03340.atmos_average.ucomp``)
        :type varfull: str
        
        :param plot_type: plot type (e.g.,
            ``Plot 2D lon X time``)
        :type plot_type: str
        
        :return: topography or ``None`` if no matching ``fixed`` file

        """

        if not "[" in varfull:
            # If overwriting dim, get new dimension and trim varfull
            # from {lev=5.}
            if "{" in varfull:
                varfull, _, _, _ = get_overwrite_dim_2D(
                    varfull, plot_type, self.fdim1, self.fdim2, self.ftod)
            sol_array, filetype, var, simuID = split_varfull(varfull)
        # Recognize an operation on the variables
        else:
            # Extract individual variables and prepare for execution
            varfull = remove_whitespace(varfull)
            varfull_list = get_list_varfull(varfull)
            f = get_list_varfull(varfull)
            sol_array, filetype, var, simuID = split_varfull(varfull_list[0])

        # If requesting a lat-lon plot for 01336.atmos_average.nc,
        # try to find matching 01336.fixed.nc
        try:
            f, var_info, dim_info, dims = prep_file(
                "zsurf", "fixed", simuID, sol_array)
            # Get file type (fixed, diurn, average, daily)
            # and interp type (pfull, zstd, etc.)
            zsurf = f.variables["zsurf"][:, :]
            f.close()
        except:
            # No corresponding fixed file, return None
            zsurf = None
        return zsurf


    def do_plot(self):
        # Create figure
        ax = super(Fig_2D_lon_lat, self).fig_init()
        try:
            # Try to create figure, else return error
            lon, lat, var, var_info = super(
                Fig_2D_lon_lat, self).data_loader_2D(self.varfull,
                                                     self.plot_type)
            lon_shift, var = shift_data(lon, var)
            try:
                # Try to get topography if a matching fixed file exists
                _, zsurf = shift_data(lon, zsurf)
                add_topo = True
            except:
                add_topo = False

            projfull = self.axis_opt3

            # ----------------------------------------------------------
            # If proj = cart, use generic contours from Fig_2D() class
            # ----------------------------------------------------------
            if projfull == "cart":

                super(Fig_2D_lon_lat, self).filled_contour(lon_shift,
                                                           lat, var)
                if add_topo:
                    # Add topography contour
                    plt.contour(lon_shift, lat, zsurf, 11, colors = "k",
                                linewidths = 0.5, linestyles = "solid")

                if self.varfull2:
                    (_, _, var2,var_info2) = super(
                        Fig_2D_lon_lat, self).data_loader_2D(self.varfull2,
                                                             self.plot_type)
                    lon_shift, var2 = shift_data(lon, var2)
                    super(Fig_2D_lon_lat, self).solid_contour(lon_shift,
                                                              lat, var2,
                                                              self.contour2)
                    var_info += f" (& {var_info2})"

                if self.Xlim:
                    plt.xlim(self.Xlim[0], self.Xlim[1])
                if self.Ylim:
                    plt.ylim(self.Ylim[0], self.Ylim[1])

                super(Fig_2D_lon_lat, self).make_title(var_info, "Longitude",
                                                       "Latitude")
             # --- Annotation---
                ax.xaxis.set_major_locator(MultipleLocator(30))
                ax.xaxis.set_minor_locator(MultipleLocator(10))
                ax.yaxis.set_major_locator(MultipleLocator(15))
                ax.yaxis.set_minor_locator(MultipleLocator(5))
                plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                           rotation = 0)
                plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                           rotation = 0)

            # ----------------------------------------------------------
            #                      Special Projections
            # ----------------------------------------------------------
            else:
                cmap = self.axis_opt1
                if cmap == "wbr":
                    cmap = wbr_cmap()
                if cmap == "rjw":
                    cmap = rjw_cmap()
                norm, levs = super(Fig_2D_lon_lat, self).return_norm_levs()

                ax.axis("off")
                # Nans = white for projections
                ax.patch.set_color("1")
                # ------------------------------------------------------

                if projfull == "robin":
                    LON, LAT = np.meshgrid(lon_shift, lat)
                    X, Y = robin2cart(LAT, LON)

                    for mer in np.arange(-180, 180, 30):
                        # Add meridans and parallels
                        xg, yg = robin2cart(lat, lat*0+mer)
                        plt.plot(xg, yg, ":k", lw = 0.5)

                    for mer in np.arange(-180, 181, 90):
                        # Label every other meridian
                        xl, yl = robin2cart(lat.min(), mer)
                        lab_txt = format_lon_lat(mer, "lon")
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size
                                             - self.nPan*label_factor),
                                 verticalalignment = "top",
                                 horizontalalignment = "center")

                    for par in np.arange(-60, 90, 30):
                        xg, yg = robin2cart(lon_shift*0+par, lon_shift)
                        plt.plot(xg, yg, ":k", lw = 0.5)
                        xl, yl = robin2cart(par, 180)
                        lab_txt = format_lon_lat(par, "lat")
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size
                                             - self.nPan*label_factor))
                # ------------------------------------------------------

                if projfull == "moll":
                    LON, LAT = np.meshgrid(lon_shift, lat)
                    X, Y = mollweide2cart(LAT, LON)

                    for mer in np.arange(-180, 180, 30):
                        # Add meridans
                        xg, yg = mollweide2cart(lat, lat*0+mer)
                        plt.plot(xg, yg, ":k", lw = 0.5)

                    for mer in [-180, 0, 180]:
                        # Label every other meridian
                        xl, yl = mollweide2cart(lat.min(), mer)
                        lab_txt = format_lon_lat(mer, "lon")
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size
                                             - self.nPan*label_factor),
                                 verticalalignment = "top",
                                 horizontalalignment = "center")

                    for par in np.arange(-60, 90, 30):
                        # Add parallels
                        xg, yg = mollweide2cart(lon_shift*0+par, lon_shift)
                        xl, yl = mollweide2cart(par, 180)
                        lab_txt = format_lon_lat(par, "lat")
                        plt.plot(xg, yg, ":k", lw = 0.5)
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size
                                             - self.nPan*label_factor))

                if projfull[0:5] in ["Npole", "Spole", "ortho"]:
                    # Common to all azimuthal projections
                    ax.set_aspect("equal")
                    lon180_original = lon_shift.copy()
                    var, lon_shift = add_cyclic(var, lon_shift)
                    if add_topo:
                        zsurf, _ = add_cyclic(zsurf, lon180_original)
                    lon_lat_custom = None
                    lat_b = None

                    if len(projfull) > 5:
                        # Get custom lat-lon, if any
                        lon_lat_custom = filter_input(projfull[5:], "float")

                if projfull[0:5] == "Npole":
                    # Reduce data
                    lat_b = 60
                    if not(lon_lat_custom is None):
                        # Bounding latitude
                        lat_b = lon_lat_custom
                    lat_bi, _ = get_lat_index(lat_b, lat)
                    lat = lat[lat_bi:]
                    var = var[lat_bi:, :]
                    if add_topo:
                        zsurf = zsurf[lat_bi:, :]
                    LON, LAT = np.meshgrid(lon_shift, lat)
                    X, Y = azimuth2cart(LAT, LON, 90, 0)

                    for mer in np.arange(-180, 180, 30):
                        # Add meridans and parallels
                        xg, yg = azimuth2cart(lat, lat*0+mer, 90)
                        plt.plot(xg, yg, ":k", lw = 0.5)

                    for mer in np.arange(-150, 180, 30):
                        # Skip 190°W to leave room for title
                        # Place label 3° S of bounding latitude
                        xl, yl = azimuth2cart(lat.min()-3, mer, 90)
                        lab_txt = format_lon_lat(mer, "lon")
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size - self.nPan*label_factor),
                                 verticalalignment = "top",
                                 horizontalalignment = "center")

                    for par in np.arange(80, lat.min(), -10):
                        # Parallels start from 80°N, every 10°
                        xg, yg = azimuth2cart(lon_shift*0+par, lon_shift, 90)
                        plt.plot(xg, yg, ":k", lw = 0.5)
                        xl, yl = azimuth2cart(par, 180, 90)
                        lab_txt = format_lon_lat(par, "lat")
                        plt.text(xl, yl, lab_txt, fontsize = 5)

                if projfull[0:5] == "Spole":
                    lat_b = -60
                    if not(lon_lat_custom is None):
                        # Bounding latitude
                        lat_b = lon_lat_custom
                    lat_bi, _ = get_lat_index(lat_b, lat)
                    lat = lat[:lat_bi]
                    var = var[:lat_bi, :]
                    if add_topo:
                        zsurf = zsurf[:lat_bi, :]
                    LON, LAT = np.meshgrid(lon_shift, lat)
                    X, Y = azimuth2cart(LAT, LON, -90, 0)

                    for mer in np.arange(-180, 180, 30):
                        # Add meridans and parallels
                        xg, yg = azimuth2cart(lat, lat*0+mer, -90)
                        plt.plot(xg, yg, ":k", lw = 0.5)

                    for mer in np.append(np.arange(-180, 0, 30),
                                         np.arange(30, 180, 30)):
                        # Skip 0 to leave room for title
                        # Place label 3°N of bounding latitude
                        xl, yl = azimuth2cart(lat.max()+3, mer, -90)
                        lab_txt = format_lon_lat(mer, "lon")
                        plt.text(xl, yl, lab_txt,
                                 fontsize = (label_size
                                             - self.nPan*label_factor),
                                 verticalalignment = "top",
                                 horizontalalignment = "center")

                    for par in np.arange(-80, lat.max(), 10):
                        # Parallels start from 80°S, every 10°
                        xg, yg = azimuth2cart(lon_shift*0+par, lon_shift, -90)
                        plt.plot(xg, yg, ":k", lw = 0.5)
                        xl, yl = azimuth2cart(par, 180, -90)
                        lab_txt = format_lon_lat(par, "lat")
                        plt.text(xl, yl, lab_txt, fontsize = 5)

                if projfull[0:5] == "ortho":
                    # Initialization
                    lon_p, lat_p = -120, 20
                    if not(lon_lat_custom is None):
                        lon_p = lon_lat_custom[0]
                        # Bounding latitude
                        lat_p = lon_lat_custom[1]
                    LON, LAT = np.meshgrid(lon_shift, lat)
                    # Mask opposite side of planet
                    X, Y, MASK = ortho2cart(LAT, LON, lat_p, lon_p)
                    var = var*MASK
                    if add_topo:
                        zsurf = zsurf*MASK

                    for mer in np.arange(-180, 180, 30):
                        # Add meridans and parallels
                        xg, yg, maskg = ortho2cart(
                            lat, lat*0+mer, lat_p, lon_p)
                        plt.plot(xg*maskg, yg, ":k", lw = 0.5)
                    for par in np.arange(-60, 90, 30):
                        xg, yg, maskg = ortho2cart(
                            lon_shift*0+par, lon_shift, lat_p, lon_p)
                        plt.plot(xg*maskg, yg, ":k", lw = 0.5)

                if self.range:
                    plt.contourf(X, Y, var, levs, extend = "both",
                                 cmap = cmap, norm = norm)
                else:
                    plt.contourf(X, Y, var, levels, cmap = cmap, norm = norm)

                super(Fig_2D_lon_lat, self).make_colorbar(levs)

                if add_topo:
                    # Add topography contours
                    plt.contour(X, Y, zsurf, 11, colors = "k",
                                linewidths = 0.5, linestyles = "solid")

                # ======================================================
                # =========== Solid Contour 2nd Variable ===============
                # ======================================================
                if self.varfull2:
                    lon, lat, var2, var_info2 = super(
                        Fig_2D_lon_lat, self).data_loader_2D(self.varfull2,
                                                             self.plot_type)
                    lon_shift, var2 = shift_data(lon, var2)

                    if projfull == "robin":
                        LON, LAT = np.meshgrid(lon_shift, lat)
                        X, Y = robin2cart(LAT, LON)

                    if projfull == "moll":
                        LON, LAT = np.meshgrid(lon_shift, lat)
                        X, Y = mollweide2cart(LAT, LON)

                    if projfull[0:5] in ["Npole", "Spole", "ortho"]:
                        # Common to all azithumal projections
                        var2, lon_shift = add_cyclic(var2, lon_shift)
                        lon_lat_custom = None
                        lat_b = None

                        if len(projfull) > 5:
                            # Get custom lat-lon, if any
                            lon_lat_custom = filter_input(
                                projfull[5:], "float")

                    if projfull[0:5] == "Npole":
                        # Reduce data
                        lat_b = 60
                        if not(lon_lat_custom is None):
                            # Bounding latitude
                            lat_b = lon_lat_custom
                        lat_bi, _ = get_lat_index(lat_b, lat)
                        lat = lat[lat_bi:]
                        var2 = var2[lat_bi:, :]
                        LON, LAT = np.meshgrid(lon_shift, lat)
                        X, Y = azimuth2cart(LAT, LON, 90, 0)
                    if projfull[0:5] == "Spole":
                        lat_b = -60
                        if not(lon_lat_custom is None):
                            # Bounding latitude
                            lat_b = lon_lat_custom
                        lat_bi, _ = get_lat_index(lat_b, lat)
                        lat = lat[:lat_bi]
                        var2 = var2[:lat_bi, :]
                        LON, LAT = np.meshgrid(lon_shift, lat)
                        X, Y = azimuth2cart(LAT, LON, -90, 0)

                    if projfull[0:5] == "ortho":
                        lon_p, lat_p = -120, 20
                        if not(lon_lat_custom is None):
                            lon_p = lon_lat_custom[0]
                            # Bounding latitude
                            lat_p = lon_lat_custom[1]
                        LON, LAT = np.meshgrid(lon_shift, lat)
                        # Mask opposite side of planet
                        X, Y, MASK = ortho2cart(LAT, LON, lat_p, lon_p)
                        var2 = var2 * MASK

                    # Prevent error message for "contours not found"
                    np.seterr(divide="ignore", invalid="ignore")

                    if self.contour2 is None:
                        CS = plt.contour(
                            X, Y, var2, 11, colors = "k", linewidths = 2)
                    else:
                        # 1 contour provided (float), convert to array
                        if type(self.contour2) == float:
                            self.contour2 = [self.contour2]
                        CS = plt.contour(X, Y, var2, self.contour2,
                                         colors = "k", linewidths = 2)
                    plt.clabel(CS, inline = 1, fontsize = 14, fmt = "%g")

                    var_info += f" (& {var_info2})"

                if self.title:
                    plt.title((self.title),
                              fontsize = (title_size - self.nPan*title_factor))
                else:
                    plt.title(f"{var_info}\n{self.fdim_txt[1:]}",
                              fontsize = (title_size - self.nPan*title_factor),
                              wrap = False)

            self.success = True

        except Exception as e:
            super(Fig_2D_lon_lat, self).exception_handler(e, ax)

        super(Fig_2D_lon_lat, self).fig_save()


class Fig_2D_time_lat(Fig_2D):

    def make_template(self):
        # Calls method from parent class
        super(Fig_2D_time_lat, self).make_template("Plot 2D time X lat",
                                                   "Lon +/-180",
                                                   "Level [Pa/m]",
                                                   "Ls", "Lat")


    def do_plot(self):
        # Create figure
        ax = super(Fig_2D_time_lat, self).fig_init()
        try:
            # Try to create figure, else return error
            t_stack, lat, var, var_info = super(
                Fig_2D_time_lat, self).data_loader_2D(self.varfull,
                                                      self.plot_type)
            SolDay = t_stack[0, :]
            LsDay = t_stack[1, :]

            super(Fig_2D_time_lat, self).filled_contour(LsDay, lat, var)

            if self.varfull2:
                _, _, var2, var_info2 = super(
                    Fig_2D_time_lat, self).data_loader_2D(self.varfull2,
                                                          self.plot_type)
                super(Fig_2D_time_lat, self).solid_contour(LsDay, lat,
                                                           var2, self.contour2)
                var_info += f" (& {var_info2})"

            # Axis formatting
            if self.Xlim:
                idmin = np.argmin(abs(SolDay - self.Xlim[0]))
                idmax = np.argmin(abs(SolDay - self.Xlim[1]))
                plt.xlim([LsDay[idmin], LsDay[idmax]])

            if self.Ylim:
                plt.ylim(self.Ylim[0], self.Ylim[1])

            Ls_ticks = [item for item in ax.get_xticks()]
            labels = [item for item in ax.get_xticklabels()]

            for i in range(0, len(Ls_ticks)):
                # Find timestep closest to this tick
                id = np.argmin(abs(LsDay-Ls_ticks[i]))
                if add_sol_time_axis:
                    labels[i] = (f"{np.mod(Ls_ticks[i], 360.):g}{degr}"
                                 f"\nsol {SolDay[id]}")
                else:
                    labels[i] = (f"{np.mod(Ls_ticks[i], 360.):g}{degr}")
            #Clean-up Ls labels at the edges.
            labels[0]='';labels[-1]=''
            ax.set_xticks(Ls_ticks)
            ax.set_xticklabels(labels,
                               fontsize = (label_size - self.nPan*tick_factor),
                               rotation = 0)

            super(Fig_2D_time_lat, self).make_title(var_info, "L$_s$",
                                                    "Latitude")

            ax.yaxis.set_major_locator(MultipleLocator(15))
            ax.yaxis.set_minor_locator(MultipleLocator(5))
            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)

            self.success = True

        except Exception as e:
            super(Fig_2D_time_lat, self).exception_handler(e, ax)

        super(Fig_2D_time_lat, self).fig_save()


class Fig_2D_lat_lev(Fig_2D):


    def make_template(self):
        # Calls method from parent class
        super(Fig_2D_lat_lev, self).make_template("Plot 2D lat X lev",
                                                  "Ls 0-360 ", "Lon +/-180",
                                                  "Lat", "Level[Pa/m]")


    def do_plot(self):
        # Create figure
        ax = super(Fig_2D_lat_lev, self).fig_init()
        try:
            # Try to create figure, else return error
            lat, pfull, var, var_info = super(
                Fig_2D_lat_lev, self).data_loader_2D(self.varfull,
                                                     self.plot_type)
            super(Fig_2D_lat_lev, self).filled_contour(lat, pfull, var)

            if self.varfull2:
                _, _, var2, var_info2 = super(
                    Fig_2D_lat_lev, self).data_loader_2D(self.varfull2,
                                                         self.plot_type)
                super(Fig_2D_lat_lev, self).solid_contour(
                    lat, pfull, var2, self.contour2)
                var_info += f" (& {var_info2})"

            if self.vert_unit == "Pa":
                ax.set_yscale("log")
                ax.invert_yaxis()
                ax.yaxis.set_major_formatter(CustomTicker())
                ax.yaxis.set_minor_formatter(NullFormatter())
                ylabel_txt = "Pressure [Pa]"
            else:
                ylabel_txt = "Altitude [m]"

            if self.Xlim:
                plt.xlim(self.Xlim)
            if self.Ylim:
                plt.ylim(self.Ylim)

            super(Fig_2D_lat_lev, self).make_title(var_info, "Latitude",
                                                   ylabel_txt)

            ax.xaxis.set_major_locator(MultipleLocator(15))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)

            self.success = True
        except Exception as e:
            super(Fig_2D_lat_lev, self).exception_handler(e, ax)

        super(Fig_2D_lat_lev, self).fig_save()


class Fig_2D_lon_lev(Fig_2D):


    def make_template(self):
        """
        Calls method from parent class

        """
        super(Fig_2D_lon_lev, self).make_template("Plot 2D lon X lev",
                                                  "Ls 0-360 ", "Latitude",
                                                  "Lon +/-180", "Level[Pa/m]")


    def do_plot(self):
        """
        Create figure

        """
        ax = super(Fig_2D_lon_lev, self).fig_init()
        try:
            # Try to create figure, else return error
            lon, pfull, var, var_info = super(
                Fig_2D_lon_lev, self).data_loader_2D(self.varfull,
                                                     self.plot_type)
            lon_shift, var = shift_data(lon, var)

            super(Fig_2D_lon_lev, self).filled_contour(lon_shift, pfull, var)

            if self.varfull2:
                _, _, var2, var_info2 = super(
                    Fig_2D_lon_lev, self).data_loader_2D(self.varfull2,
                                                         self.plot_type)
                _, var2 = shift_data(lon, var2)
                super(Fig_2D_lon_lev, self).solid_contour(
                    lon_shift, pfull, var2, self.contour2)
                var_info += f" (& {var_info2})"

            if self.vert_unit == "Pa":
                ax.set_yscale("log")
                ax.invert_yaxis()
                ax.yaxis.set_major_formatter(CustomTicker())
                ax.yaxis.set_minor_formatter(NullFormatter())
                ylabel_txt = "Pressure [Pa]"
            else:
                ylabel_txt = "Altitude [m]"

            if self.Xlim:
                plt.xlim(self.Xlim)
            if self.Ylim:
                plt.ylim(self.Ylim)

            super(Fig_2D_lon_lev, self).make_title(
                var_info, "Longitude", ylabel_txt)

            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(10))
            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)

            self.success = True
        except Exception as e:
            super(Fig_2D_lon_lev, self).exception_handler(e, ax)

        super(Fig_2D_lon_lev, self).fig_save()


class Fig_2D_time_lev(Fig_2D):


    def make_template(self):
        # Calls method from parent class
        super(Fig_2D_time_lev, self).make_template("Plot 2D time X lev",
                                                   "Latitude", "Lon +/-180",
                                                   "Ls", "Level[Pa/m]")


    def do_plot(self):
        # Create figure
        ax = super(Fig_2D_time_lev, self).fig_init()
        try:
            # Try to create figure, else return error
            t_stack, pfull, var, var_info = super(
                Fig_2D_time_lev, self).data_loader_2D(self.varfull,
                                                      self.plot_type)
            SolDay = t_stack[0, :]
            LsDay = t_stack[1, :]
            super(Fig_2D_time_lev, self).filled_contour(LsDay, pfull, var)

            if self.varfull2:
                _, _, var2, var_info2 = super(
                    Fig_2D_time_lev, self).data_loader_2D(self.varfull2,
                                                          self.plot_type)
                super(Fig_2D_time_lev, self).solid_contour(LsDay, pfull,
                                                           var2, self.contour2)
                var_info += f" (& {var_info2})"

            # Axis formatting
            if self.Xlim:
                idmin = np.argmin(abs(SolDay - self.Xlim[0]))
                idmax = np.argmin(abs(SolDay - self.Xlim[1]))
                plt.xlim([LsDay[idmin], LsDay[idmax]])
            if self.Ylim:
                plt.ylim(self.Ylim)

            Ls_ticks = [item for item in ax.get_xticks()]
            labels = [item for item in ax.get_xticklabels()]

            for i in range(0, len(Ls_ticks)):
                # Find timestep closest to this tick
                id = np.argmin(abs(LsDay-Ls_ticks[i]))
                if add_sol_time_axis:
                    labels[i] = (f"{np.mod(Ls_ticks[i], 360.)}{degr}"
                                 f"\nsol {SolDay[id]}")
                else:
                    labels[i] = f"{np.mod(Ls_ticks[i], 360.)}{degr}"
            #Clean-up Ls labels at the edges.
            labels[0]='';labels[-1]=''
            ax.set_xticks(Ls_ticks)
            ax.set_xticklabels(labels,
                               fontsize = label_size - self.nPan*tick_factor,
                               rotation = 0)

            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)

            if self.vert_unit == "Pa":
                ax.set_yscale("log")
                ax.invert_yaxis()
                ax.yaxis.set_major_formatter(CustomTicker())
                ax.yaxis.set_minor_formatter(NullFormatter())
                ylabel_txt = "Pressure [Pa]"
            else:
                ylabel_txt = "Altitude [m]"

            super(Fig_2D_time_lev, self).make_title(
                var_info, "L$_s$", ylabel_txt)

            self.success = True

        except Exception as e:
            super(Fig_2D_time_lev, self).exception_handler(e, ax)

        super(Fig_2D_time_lev, self).fig_save()


class Fig_2D_lon_time(Fig_2D):


    def make_template(self):
        # Calls method from parent class
        super(Fig_2D_lon_time, self).make_template("Plot 2D lon X time",
                                                   "Latitude", "Level [Pa/m]",
                                                   "Lon +/-180", "Ls")


    def do_plot(self):
        # Create figure
        ax = super(Fig_2D_lon_time, self).fig_init()
        try:
            # Try to create figure, else return error
            lon, t_stack, var, var_info = super(
                Fig_2D_lon_time, self).data_loader_2D(self.varfull,
                                                      self.plot_type)
            lon_shift, var = shift_data(lon, var)

            SolDay = t_stack[0, :]
            LsDay = t_stack[1, :]
            super(Fig_2D_lon_time, self).filled_contour(lon_shift, LsDay, var)

            if self.varfull2:
                _, _, var2, var_info2 = super(
                    Fig_2D_lon_time, self).data_loader_2D(self.varfull2,
                                                          self.plot_type)
                _, var2 = shift_data(lon, var2)
                super(Fig_2D_lon_time, self).solid_contour(lon_shift, LsDay,
                                                           var2, self.contour2)
                var_info += (f" (& {var_info2})")

            # Axis formatting
            if self.Xlim:
                plt.xlim(self.Xlim)

            if self.Ylim:
                idmin = np.argmin(abs(SolDay - self.Ylim[0]))
                idmax = np.argmin(abs(SolDay - self.Ylim[1]))
                plt.ylim([LsDay[idmin], LsDay[idmax]])

            Ls_ticks = [item for item in ax.get_yticks()]
            labels = [item for item in ax.get_yticklabels()]

            for i in range(0, len(Ls_ticks)):
                # Find timestep closest to this tick
                id = np.argmin(abs(LsDay-Ls_ticks[i]))
                if add_sol_time_axis:
                    labels[i] = (f"{np.mod(Ls_ticks[i], 360.):g}{degr}"
                                 f"\nsol {SolDay[id]}")
                else:
                    labels[i] = (f"{np.mod(Ls_ticks[i], 360.):g}{degr}")
            ax.set_yticklabels(labels,
                               fontsize = label_size - self.nPan*tick_factor,
                               rotation = 0)

            ax.xaxis.set_major_locator(MultipleLocator(30))
            ax.xaxis.set_minor_locator(MultipleLocator(10))

            super(Fig_2D_lon_time, self).make_title(
                var_info, "Longitude", "L$_s$")
            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)

            self.success = True

        except Exception as e:
            super(Fig_2D_lon_time, self).exception_handler(e, ax)

        super(Fig_2D_lon_time, self).fig_save()


class Fig_1D(object):
    # Parent class for 1D figure
    def __init__(self, varfull="atmos_average.ts", doPlot=True):

        self.title = None
        self.legend = None
        self.varfull = varfull
        self.t = "AXIS" # Default value for AXIS
        self.lat = None
        self.lon = None
        self.lev = None
        self.ftod = None # Time of day, requested input
        self.hour = None # Hour of day, bool, for diurn plots only
        # Logic
        self.doPlot = doPlot
        self.plot_type = "1D_time"

        # Extract filetype, variable, and simulation ID
        # (initialization only)
        (self.sol_array, self.filetype,
         self.var, self.simuID) = split_varfull(self.varfull)

        # Multipanel
        self.nPan = 1
        self.subID = 1
        self.addLine = False
        # Page layout, e.g., [2,3] if HOLD ON 2,3
        self.layout = None
        # Annotation for free dims
        self.fdim_txt = ""
        self.success = False
        # Vertical unit is m or Pa
        self.vert_unit = ""
        # Axis options

        # Dim limit
        self.Dlim = None
        # Variable limit
        self.Vlim = None
        self.axis_opt1 = "-"


    def make_template(self):
        customFileIN.write(
            f"<<<<<<<<<<<<<<| Plot 1D = {self.doPlot} |>>>>>>>>>>>>>\n")
        customFileIN.write(f"Title          = {self.title}\n")   # 1
        customFileIN.write(f"Legend         = {self.legend}\n")  # 2
        customFileIN.write(f"Main Variable  = {self.varfull}\n") # 3
        customFileIN.write(f"Ls 0-360       = {self.t}\n")       # 4
        customFileIN.write(f"Latitude       = {self.lat}\n")     # 5
        customFileIN.write(f"Lon +/-180     = {self.lon}\n")     # 6
        customFileIN.write(f"Level [Pa/m]   = {self.lev}\n")     # 7
        customFileIN.write(f"Diurnal  [hr]  = {self.hour}\n")    # 8
        customFileIN.write(
            f"Axis Options  : lat,lon+/-180,[Pa/m],Ls = [None,None] | "
            f"var = [None,None] | linestyle = - | axlabel = None \n") # 9


    def read_template(self):
        self.title = rT("char")     # 1
        self.legend = rT("char")    # 2
        self.varfull = rT("char")   # 3
        self.t = rT("float")        # 4
        self.lat = rT("float")      # 5
        self.lon = rT("float")      # 6
        self.lev = rT("float")      # 7
        self.hour = rT("float")     # 8
        (self.Dlim, self.Vlim,
         self.axis_opt1, self.axis_opt2, _) = read_axis_options(
             customFileIN.readline()) # 7

        self.plot_type = self.get_plot_type()


    def get_plot_type(self):
        """
        Note that the ``self.t == "AXIS" test`` and the
        ``self.t = -88888`` assignment are only used when MarsPlot is
        not passed a template.

        :return: type of 1D plot to create (1D_time, 1D_lat, etc.)

        """
        ncheck = 0
        graph_type = "Error"
        if self.t == -88888 or self.t == "AXIS":
            self.t = -88888
            graph_type = "1D_time"
            ncheck += 1
        if self.lat == -88888 or self.lat == "AXIS":
            self.lat = -88888
            graph_type = "1D_lat"
            ncheck += 1
        if self.lon == -88888 or self.lon == "AXIS":
            self.lon = -88888
            graph_type = "1D_lon"
            ncheck += 1
        if self.lev == -88888 or self.lev == "AXIS":
            self.lev = -88888
            graph_type = "1D_lev"
            ncheck += 1
        if self.hour == -88888 or self.hour == "AXIS":
            self.hour = -88888
            graph_type = "1D_diurn"
            ncheck += 1
        if ncheck == 0:
            print(f"{Yellow}*** Warning *** In 1D plot, {self.varfull}: use "
                  f"``AXIS`` to set the varying dimension{Nclr}")
        if ncheck > 1:
            print(f"{Yellow}*** Warning *** In 1D plot, {self.varfull}: "
                  f"``AXIS`` keyword can only be used once{Nclr}")
        return graph_type


    def data_loader_1D(self, varfull, plot_type):

        if not "[" in varfull:
            if "{" in varfull:
                (varfull, t_req, lat_req,
                 lon_req, lev_req,
                 ftod_req) = get_overwrite_dim_1D(varfull, self.t,
                                                  self.lat,self.lon,
                                                  self.lev, self.ftod)
                # t_req, lat_req, lon_req, lev_req contain dims to
                # overwrite if "{}" provided, else default to self.t,
                # self.lat, self.lon,
                # self.lev
            else:
                # No "{}" to overwrite dims, copy plot defaults
                t_req = self.t
                lat_req = self.lat
                lon_req = self.lon
                lev_req = self.lev
                ftod_req = self.ftod

            sol_array, filetype, var, simuID = split_varfull(varfull)
            xdata, var, var_info = self.read_NCDF_1D(var, filetype, simuID,
                                                     sol_array, plot_type,
                                                     t_req, lat_req, lon_req,
                                                     lev_req, ftod_req)
            leg_text = f"{var_info}"
            varlabel = f"{var_info}"

        else:
            VAR = []
            # Extract individual variables and prepare for execution
            varfull = remove_whitespace(varfull)
            varfull_list = get_list_varfull(varfull)
            expression_exec = create_exec(varfull, varfull_list)

            # Initialize list of requested dims
            t_list = [None] * len(varfull_list)
            lat_list = [None] * len(varfull_list)
            lon_list = [None] * len(varfull_list)
            lev_list = [None] * len(varfull_list)
            ftod_list = [None] * len(varfull_list)
            expression_exec = create_exec(varfull, varfull_list)

            for i in range(0, len(varfull_list)):
                # If overwriting dim, get new dim and trim varfull from
                # {lev=5.}
                if "{" in varfull_list[i]:
                    (varfull_list[i], t_list[i],
                     lat_list[i], lon_list[i],
                     lev_list[i], ftod_list[i]) = get_overwrite_dim_1D(
                        varfull_list[i], self.t, self.lat,
                        self.lon, self.lev, self.ftod)
                else:
                    # No "{}" to overwrite dims, copy plot defaults
                    t_list[i] = self.t
                    lat_list[i] = self.lat
                    lon_list[i] = self.lon
                    lev_list[i] = self.lev
                    ftod_list[i] = self.ftod

                (sol_array, filetype,
                 var, simuID) = split_varfull(varfull_list[i])
                xdata, temp, var_info = self.read_NCDF_1D(var,
                                                          filetype,
                                                          simuID,
                                                          sol_array,
                                                          plot_type,
                                                          t_list[i],
                                                          lat_list[i],
                                                          lon_list[i],
                                                          lev_list[i],
                                                          ftod_list[i])
                VAR.append(temp)
            leg_text = (f"{var} {var_info.split(' ')[-1]}"
                        f"{expression_exec.split(']')[-1]}")
            varlabel = f"{var}"
            var_info = varfull
            var = eval(expression_exec) #TODO removed ,namespace

        return xdata, var, var_info, leg_text, varlabel


    def read_NCDF_1D(self, var_name, file_type, simuID, sol_array,
                     plot_type, t_req, lat_req, lon_req, lev_req, ftod_req):
        """
        Parse a Main Variable expression object that includes a square
        bracket [] (for variable calculations) for the variable to
        plot.

        :param var_name: variable name (e.g., ``temp``)
        :type var_name: str
        
        :param file_type: MGCM output file type. Must be ``fixed`` or
            ``average``
        :type file_type: str
        
        :param simuID: number identifier for netCDF file directory
        :type simuID: str
        
        :param sol_array: sol if different from default
            (e.g., ``02400``)
        :type sol_array:  str
        
        :param plot_type: ``1D_lon``, ``1D_lat``, ``1D_lev``, or
            ``1D_time``
        :type plot_type: str
        
        :param t_req: Ls requested
        :type t_req: str
        
        :param lat_req: lat requested
        :type lat_req: str
        
        :param lon_req: lon requested
        :type lon_req: str
        
        :param lev_req: level [Pa/m] requested
        :type lev_req: str
        
        :param ftod_req: time of day requested
        :type ftod_req: str

        :return: (dim_array) the axis (e.g., an array of longitudes),
                 (var_array) the variable extracted

        """

        f, var_info, dim_info, dims = prep_file(var_name, file_type,
                                                simuID, sol_array)

        # Get file type (fixed, diurn, average, daily) and interp type
        # (pfull, zstd, etc.)
        f_type, interp_type = FV3_file_type(f)

        add_fdim = False
        if not self.fdim_txt.strip():
            # If self.fdim is empty, add variable (do only once)
            add_fdim = True

        # Initialize dims in all .nc files
        lat = f.variables["lat"][:]
        lati = np.arange(0, len(lat))
        lon = f.variables["lon"][:]
        loni = np.arange(0, len(lon))

        # ------------------------Time of Day --------------------------
        # *** Performed only for 1D_lat, 1D_lev, or 1D_time plots ***
        #                     from a diurn file
        # For plotting 1D_lat, 1D_lev, or 1D_time figures from diurn
        # files, select data on time of day axis and update dims so that
        # resulting variable is in format of average and daily files.
        # This simplifies logic so that all daily, average, and diurn
        # files are treated the same. Naturally, plot type 1D_diurn is
        # an exeception & following lines are skipped.

        if ((f_type == "diurn" and
             dim_info[1][:11] == "time_of_day") and not
            plot_type == "1D_diurn"):
            # Time of day is always 2nd dim (dim_info[1])
            tod = f.variables[dim_info[1]][:]
            todi, temp_txt = get_tod_index(ftod_req, tod)
            # Update dim_info from
            # time, time_of_day_XX, lat, lon -> time, lat, lon
            # OR
            # time, time_of_day_XX, pfull, lat, lon -> time, pfull, lat, lon
            dim_info = (dim_info[0],) + dim_info[2:]
            if add_fdim:
                self.fdim_txt += temp_txt

        # Static: Ignore level and time dims
        if dim_info == (u"lat", u"lon"):
            if plot_type == "1D_lat":
                loni, temp_txt = get_lon_index(lon_req, lon)
            elif plot_type == "1D_lon":
                lati, temp_txt = get_lat_index(lat_req, lat)

            if add_fdim:
                self.fdim_txt += temp_txt
            var = f.variables[var_name][lati, loni].reshape(
                len(np.atleast_1d(lati)), len(np.atleast_1d(loni))
            )
            f.close()
            w = area_weights_deg(var.shape, lat[lati])

            if plot_type == "1D_lat":
                return lat, mean_func(var, axis = 1), var_info
            if plot_type == "1D_lon":
                return lon, np.average(var, weights = w, axis = 0), var_info

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        #       ~~ For 1D_time, 1D_lat, 1D_lon, and 1D_lev only ~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if not plot_type == "1D_diurn":
            # ====== time, lat, lon =======
            if dim_info == (u"time", u"lat", u"lon"):

                # Initialize dimension
                t = f.variables["time"][:]
                LsDay = np.squeeze(f.variables["areo"][:])
                ti = np.arange(0, len(t))

                if f_type == "diurn" and len(LsDay.shape) > 1:
                    # For diurn file, change time_of_day[time, 24, 1] to
                    # time_of_day[time] at midnight UT
                    LsDay = np.squeeze(LsDay[:, 0])

                # Stack time and areo arrays as 1 variable
                t_stack = np.vstack((t, LsDay))

                if plot_type == "1D_lat":
                    ti, temp_txt = get_time_index(t_req, LsDay)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    loni, temp_txt = get_lon_index(lon_req, lon)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                if plot_type == "1D_lon":
                    lati, temp_txt = get_lat_index(lat_req, lat)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    ti, temp_txt = get_time_index(t_req, LsDay)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                if plot_type == "1D_time":
                    loni, temp_txt = get_lon_index(lon_req, lon)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    lati, temp_txt = get_lat_index(lat_req, lat)
                    if add_fdim:
                        self.fdim_txt += temp_txt

                if f_type == "diurn":
                    var = f.variables[var_name][ti, todi, lati, loni].reshape(
                        len(np.atleast_1d(ti)),
                        len(np.atleast_1d(todi)),
                        len(np.atleast_1d(lati)),
                        len(np.atleast_1d(loni)))
                    var = mean_func(var, axis = 1)
                else:
                    var = f.variables[var_name][ti, lati, loni].reshape(
                        len(np.atleast_1d(ti)),
                        len(np.atleast_1d(lati)),
                        len(np.atleast_1d(loni)))
                f.close()

                w = area_weights_deg(var.shape, lat[lati])

                if plot_type == "1D_lat":
                    return (lat,
                            mean_func(mean_func(var, axis = 2), axis = 0),
                            var_info)
                if plot_type == "1D_lon":
                    return (lon,
                            mean_func(np.average(var, weights = w, axis = 1),
                                      axis = 0),
                            var_info)
                if plot_type == "1D_time":
                    return (t_stack,
                            mean_func(np.average(var, weights = w, axis = 1),
                                      axis = 1),
                            var_info)

            # ====== [time, lev, lat, lon] =======
            if (dim_info == (u"time", u"pfull", u"lat", u"lon")
                or dim_info == (u"time", u"level", u"lat", u"lon")
                or dim_info == (u"time", u"pstd", u"lat", u"lon")
                or dim_info == (u"time", u"zstd", u"lat", u"lon")
                or dim_info == (u"time", u"zagl", u"lat", u"lon")
                    or dim_info == (u"time", u"zgrid", u"lat", u"lon")):

                if dim_info[1] in ["pfull", "level", "pstd"]:
                    self.vert_unit = "Pa"
                if dim_info[1] in ["zagl", "zstd", "zgrid"]:
                    self.vert_unit = "m"

                # Initialize dims
                levs = f.variables[dim_info[1]][:]
                zi = np.arange(0, len(levs))
                t = f.variables["time"][:]
                LsDay = np.squeeze(f.variables["areo"][:])
                ti = np.arange(0, len(t))

                if f_type == "diurn" and len(LsDay.shape) > 1:
                    # For diurn file, change time_of_day[time, 24, 1] ->
                    # time_of_day[time] at midnight UT
                    LsDay = np.squeeze(LsDay[:, 0])

                # Stack time and areo arrays as 1 variable
                t_stack = np.vstack((t, LsDay))

                if plot_type == "1D_lat":
                    ti, temp_txt = get_time_index(t_req, LsDay)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    loni, temp_txt = get_lon_index(lon_req, lon)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    zi, temp_txt = get_level_index(lev_req, levs)
                    if add_fdim:
                        self.fdim_txt += temp_txt

                if plot_type == "1D_lon":
                    lati, temp_txt = get_lat_index(lat_req, lat)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    ti, temp_txt = get_time_index(t_req, LsDay)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    zi, temp_txt = get_level_index(lev_req, levs)
                    if add_fdim:
                        self.fdim_txt += temp_txt

                if plot_type == "1D_time":
                    loni, temp_txt = get_lon_index(lon_req, lon)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    lati, temp_txt = get_lat_index(lat_req, lat)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    zi, temp_txt = get_level_index(lev_req, levs)
                    if add_fdim:
                        self.fdim_txt += temp_txt

                if plot_type == "1D_lev":
                    ti, temp_txt = get_time_index(t_req, LsDay)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    lati, temp_txt = get_lat_index(lat_req, lat)
                    if add_fdim:
                        self.fdim_txt += temp_txt
                    loni, temp_txt = get_lon_index(lon_req, lon)
                    if add_fdim:
                        self.fdim_txt += temp_txt

                # Fix for new netCDF4 version: Get array elements
                # instead of manipulating variable
                # It used to be that var = f.variables[var_name]

                if f_type == "diurn":
                    # Do time of day average first
                    var0 = f.variables[var_name][ti, todi, zi, lati, loni]
                    var = var0.reshape(
                        len(np.atleast_1d(ti)),
                        len(np.atleast_1d(todi)),
                        len(np.atleast_1d(zi)),
                        len(np.atleast_1d(lati)),
                        len(np.atleast_1d(loni))
                    )
                    var = mean_func(var, axis = 1)
                else:
                    reshape_shape = [len(np.atleast_1d(ti)),
                                     len(np.atleast_1d(zi)),
                                     len(np.atleast_1d(lati)),
                                     len(np.atleast_1d(loni))]
                    var = f.variables[var_name][ti, zi,lati, loni].reshape(
                        reshape_shape)
                f.close()

                w = area_weights_deg(var.shape, lat[lati])

                if plot_type == "1D_lat":
                    return (lat,
                            mean_func(mean_func(mean_func(var, axis = 3), axis = 1), axis = 0),
                            var_info)
                if plot_type == "1D_lon":
                    return (lon,
                            mean_func(mean_func(np.average(var, weights = w, axis = 2), axis = 1), axis = 0),
                            var_info)
                if plot_type == "1D_time":
                    return (t_stack,
                            mean_func(mean_func(np.average(var, weights = w, axis = 2), axis = 2), axis = 1),
                            var_info)
                if plot_type == "1D_lev":
                    return (levs,
                            mean_func(mean_func(np.average(var, weights = w, axis = 2), axis = 2), axis = 0),
                            var_info)

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~ This Section is for 1D_diurn only ~~~~~~~~~~~~~~~~~
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        else:
            # Find name of time of day variable
            # (i.e., time_of_day_16 or time_of_day_24)
            tod_dim_name = find_tod_in_diurn(f)
            tod = f.variables[tod_dim_name][:]
            todi = np.arange(0, len(tod))

            # ====== time, lat, lon =======
            if f.variables[var_name].dimensions == ("time", tod_dim_name,
                                                    "lat", "lon"):

                # Initialize dim
                t = f.variables["time"][:]
                LsDay = np.squeeze(f.variables["areo"][:])
                ti = np.arange(0, len(t))

                if f_type == "diurn" and len(LsDay.shape) > 1:
                    # For diurn file, change time_of_day[time, 24, 1] ->
                    # time_of_day[time] at midnight UT
                    LsDay = np.squeeze(LsDay[:, 0])

                # Stack time and areo arrays as 1 variable
                t_stack = np.vstack((t, LsDay))

                loni, temp_txt = get_lon_index(lon_req, lon)
                if add_fdim:
                    self.fdim_txt += temp_txt
                lati, temp_txt = get_lat_index(lat_req, lat)
                if add_fdim:
                    self.fdim_txt += temp_txt
                ti, temp_txt = get_time_index(t_req, LsDay)
                if add_fdim:
                    self.fdim_txt += temp_txt

                reshape_shape = [len(np.atleast_1d(ti)),
                                 len(np.atleast_1d(tod)),
                                 len(np.atleast_1d(lati)),
                                 len(np.atleast_1d(loni))]

                # Broadcast dims before extraction. New req. for numpy
                var0 = f.variables[var_name][ti, :, lati, loni]
                var = var0.reshape(reshape_shape)
                f.close()

                w = area_weights_deg(var.shape, lat[lati])
                return (tod,
                        mean_func(mean_func(np.average(var, weights = w, axis = 2), axis = 2), axis = 0),
                        var_info)

            # ====== [time, lev, lat, lon] =======
            if (dim_info == ("time", tod_dim_name, "pfull", "lat", "lon") or
                dim_info == ("time", tod_dim_name, "level", "lat", "lon") or
                dim_info == ("time", tod_dim_name, "pstd", "lat", "lon") or
                dim_info == ("time", tod_dim_name, "zstd", "lat", "lon") or
                dim_info == ("time", tod_dim_name, "zagl", "lat", "lon") or
                dim_info == ("time", tod_dim_name, "zgrid", "lat", "lon")):

                if dim_info[1] in ["pfull", "level", "pstd"]:
                    self.vert_unit = "Pa"
                if dim_info[1] in ["zagl", "zstd", "zgrid"]:
                    self.vert_unit = "m"

                # Initialize dims
                levs = f.variables[dim_info[2]][:]

                t = f.variables["time"][:]
                LsDay = np.squeeze(f.variables["areo"][:])
                ti = np.arange(0, len(t))

                if f_type == "diurn" and len(LsDay.shape) > 1:
                    # For diurn file, change time_of_day[time, 24, 1] ->
                    # time_of_day[time] at midnight UT
                    LsDay = np.squeeze(LsDay[:, 0])

                # Stack time and areo arrays as 1 variable
                t_stack = np.vstack((t, LsDay))

                ti, temp_txt = get_time_index(t_req, LsDay)
                if add_fdim:
                    self.fdim_txt += temp_txt
                lati, temp_txt = get_lat_index(lat_req, lat)
                if add_fdim:
                    self.fdim_txt += temp_txt
                loni, temp_txt = get_lon_index(lon_req, lon)
                if add_fdim:
                    self.fdim_txt += temp_txt
                zi, temp_txt = get_level_index(lev_req, levs)
                if add_fdim:
                    self.fdim_txt += temp_txt

                reshape_shape = [len(np.atleast_1d(ti)),
                                 len(np.atleast_1d(tod)),
                                 len(np.atleast_1d(zi)),
                                 len(np.atleast_1d(lati)),
                                 len(np.atleast_1d(loni))]

                var = f.variables[var_name][ti, :, zi, lati, loni].reshape(
                    reshape_shape)
                f.close()

                w = area_weights_deg(var.shape, lat[lati])

                return (tod,
                        mean_func(mean_func(mean_func(np.average(var, weights = w, axis = 3), axis = 3), axis = 2), axis = 0),
                        var_info)


    def exception_handler(self, e, ax):
        if debug:
            raise

        sys.stdout.write("\033[F")
        sys.stdout.write("\033[K")
        print(f"{Yellow}*** Warning *** Attempting {self.plot_type} profile "
              f"for {self.varfull}: {str(e)}{Nclr}")
        ax.text(0.5, 0.5, f"ERROR:{str(e)}",
                horizontalalignment = "center",
                verticalalignment = "center",
                bbox = dict(boxstyle = "round",
                          ec = (1., 0.5, 0.5),
                          fc = (1., 0.8, 0.8),),
                transform = ax.transAxes, wrap = True, fontsize = 16)


    def fig_init(self):
        # Create figure
        if self.layout is None:
            # No layout specified
            out = fig_layout(self.subID, self.nPan, vertical_page)
        else:
            out = np.append(self.layout, self.subID)

        if self.subID == 1 and not self.addLine:
            # Create figure if first panel
            fig = plt.figure(facecolor="white",
                             figsize = (width_inch, height_inch))
        if not self.addLine:
            # nrow, ncol, subID
            ax = plt.subplot(out[0], out[1], out[2])
        else:

            ax = plt.gca()

        return ax


    def fig_save(self):
        # Save the figure
        if self.subID == self.nPan:
            # Last subplot
            if self.subID == 1:
                # If 1 plot
                if not "[" in self.varfull:
                    # Add split "{" if varfull contains layer.
                    # Does not do anything otherwise
                    sensitive_name = self.varfull.split("{")[0].strip()
                else:
                    sensitive_name = ("expression_" + get_list_varfull(
                        self.varfull)[0].split("{")[0].strip())
            else:
                # Multipanel
                sensitive_name = "multi_panel"

            self.fig_name = (
                os.path.join(output_path,"plots",f"{sensitive_name}.{out_format}")
            )
            self.fig_name = create_name(self.fig_name)

            if (i_list < len(objectList)-1 and not
                objectList[i_list+1].addLine):
                plt.savefig(self.fig_name, dpi=my_dpi)
                if out_format != "pdf":
                    print(f"Saved: {self.fig_name}")
            # Last subplot
            if i_list == len(objectList)-1:
                plt.savefig(self.fig_name, dpi=my_dpi)
                if out_format != "pdf":
                    print(f"Saved: {self.fig_name}")


    def do_plot(self):
        # Create figure
        ax = self.fig_init()
        try:
            # Try to create figure, else return error
            (xdata, var, var_info,
             leg_text, varlabel) = self.data_loader_1D(self.varfull,
                                                       self.plot_type)

            if self.legend:
                txt_label = self.legend
            else:
                # Remove 1st comma in fdim_txt to print to new line
                if self.nPan > 1:
                    txt_label = leg_text
                else:
                    # Remove 1st comma in fdim_txt to print to new line
                    txt_label = f"{var_info}\n{self.fdim_txt[1:]}"

            if self.title:
                if "{" in self.title:
                    fs = int(remove_whitespace(
                        (self.title).split("=")[1].split("}")[0]))
                    title_text = ((self.title).split("{")[0])
                    plt.title(title_text,
                              fontsize = (fs - self.nPan*title_factor),
                              wrap = False)
                else:
                    plt.title((self.title),
                              fontsize = (title_size - self.nPan*title_factor))
            else:
                plt.title(f"{var_info}\n{self.fdim_txt[1:]}",
                          fontsize = (title_size - self.nPan*title_factor),
                          wrap = False)

            if self.plot_type == "1D_lat":
                plt.plot(var, xdata, self.axis_opt1, lw = 3,
                         ms = 7, label = txt_label)
                plt.ylabel("Latitude",
                           fontsize = (label_size - self.nPan*label_factor))

                # Label provided
                if self.axis_opt2:
                    plt.xlabel(self.axis_opt2,
                               fontsize = (label_size
                                           - self.nPan*label_factor))
                else:
                    plt.xlabel(varlabel,fontsize = (label_size
                                                    - self.nPan*label_factor))

                ax.yaxis.set_major_locator(MultipleLocator(15))
                ax.yaxis.set_minor_locator(MultipleLocator(5))
                if self.Dlim:
                    plt.ylim(self.Dlim)
                if self.Vlim:
                    plt.xlim(self.Vlim)

            if self.plot_type == "1D_lon":
                lon_shift, var = shift_data(xdata, var)

                plt.plot(lon_shift, var, self.axis_opt1, lw = 3, ms = 7,
                         label = txt_label)
                plt.xlabel("Longitude",
                           fontsize = (label_size - self.nPan*label_factor))
                # Label provided
                if self.axis_opt2:
                    plt.ylabel(self.axis_opt2,
                               fontsize = (label_size
                                           - self.nPan*label_factor))
                else:
                    plt.ylabel(varlabel, fontsize = (label_size
                                                     - self.nPan*label_factor))

                ax.xaxis.set_major_locator(MultipleLocator(30))
                ax.xaxis.set_minor_locator(MultipleLocator(10))
                if self.Dlim:
                    plt.xlim(self.Dlim)
                if self.Vlim:
                    plt.ylim(self.Vlim)

            if self.plot_type == "1D_time":
                SolDay = xdata[0, :]
                LsDay = xdata[1, :]

                if args.stack_years:
                    # If simulations span different years, stack (overplot)
                    LsDay = np.mod(LsDay, 360)

                plt.plot(LsDay, var, self.axis_opt1, lw = 3, ms = 7,
                         label = txt_label)
                plt.xlabel("L$_s$",
                           fontsize = (label_size - self.nPan*label_factor))
                # Label provided
                if self.axis_opt2:
                    plt.ylabel(self.axis_opt2,
                               fontsize = (label_size
                                           - self.nPan*label_factor))
                else:
                    plt.ylabel(varlabel,fontsize = (label_size
                                                    - self.nPan*label_factor))

                # Axis formatting
                if self.Vlim:
                    plt.ylim(self.Vlim)

                if self.Dlim:
                    plt.xlim(self.Dlim) # TODO

                Ls_ticks = [item for item in ax.get_xticks()]
                labels = [item for item in ax.get_xticklabels()]

                for i in range(0, len(Ls_ticks)):
                    # Find timestep closest to this tick
                    id = np.argmin(abs(LsDay-Ls_ticks[i]))
                    if add_sol_time_axis:
                        labels[i] = (f"{np.mod(Ls_ticks[i], 360.)}{degr}"
                                     f"\nsol {SolDay[id]}")
                    else:
                        labels[i] = (f"{np.mod(Ls_ticks[i], 360.)}{degr}")
                #Clean-up Ls labels at the edges.
                labels[0]='';labels[-1]=''
                ax.set_xticks(Ls_ticks)
                ax.set_xticklabels(labels,
                                   fontsize = (label_size
                                               - self.nPan*tick_factor),
                                   rotation = 0)

            if self.plot_type == "1D_lev":
                plt.plot(var, xdata, self.axis_opt1, lw = 3, ms = 7,
                         label = txt_label)

                # Label provided
                if self.axis_opt2:
                    plt.xlabel(self.axis_opt2,
                               fontsize = (label_size
                                           - self.nPan*label_factor))
                else:
                    plt.xlabel(varlabel, fontsize = (label_size
                                                     - self.nPan*label_factor))

                if self.vert_unit == "Pa":
                    ax.set_yscale("log")
                    ax.invert_yaxis()
                    ax.yaxis.set_major_formatter(CustomTicker())
                    ax.yaxis.set_minor_formatter(NullFormatter())
                    ylabel_txt = "Pressure [Pa]"
                else:
                    ylabel_txt = "Altitude [m]"

                plt.ylabel(ylabel_txt,
                           fontsize = (label_size - self.nPan*label_factor))

                if self.Dlim:
                    plt.ylim(self.Dlim)
                if self.Vlim:
                    plt.xlim(self.Vlim)

            if self.plot_type == "1D_diurn":
                plt.plot(xdata, var, self.axis_opt1,
                         lw = 3, ms = 7, label = txt_label)
                plt.xlabel("Time [hr]",
                           fontsize = (label_size - self.nPan*label_factor))

                # Label provided
                if self.axis_opt2:
                    plt.ylabel(self.axis_opt2,
                               fontsize = (label_size
                                           - self.nPan*label_factor))
                else:
                    plt.ylabel(varlabel, fontsize = (label_size
                                                     - self.nPan*label_factor))

                ax.xaxis.set_major_locator(MultipleLocator(4))
                ax.xaxis.set_minor_locator(MultipleLocator(1))
                # Default: set X dim to 0-24. Can be overwritten.
                plt.xlim([0, 24])

                # Axis formatting
                if self.Dlim:
                    plt.xlim(self.Dlim)
                if self.Vlim:
                    plt.ylim(self.Vlim)

            # ==== Common labeling ====
            plt.xticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.yticks(fontsize = (label_size - self.nPan*tick_factor),
                       rotation = 0)
            plt.legend(fontsize = (title_size - self.nPan*title_factor))
            plt.grid(True)

            self.success = True

        except Exception as e:
            self.exception_handler(e, ax)

        self.fig_save()

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
