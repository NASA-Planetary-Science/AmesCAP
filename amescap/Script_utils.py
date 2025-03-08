# !/usr/bin/env python3
"""
Script_utils contains internal functions for processing netCDF files.
These functions can be used on their own outside of CAP if they are
imported as a module::

    from /u/path/Script_utils import MY_func

Third-party Requirements:
    * ``numpy``
    * ``netCDF4``
    * ``re``
    * ``os``
    * ``subprocess``
"""

# Load generic Python modules
import sys
import os           # Access operating system functions
import subprocess   # Run command-line commands
import numpy as np
from netCDF4 import Dataset, MFDataset
import re

# ======================================================================
#                           DEFINITIONS
# ======================================================================

global Cyan, Blue, Yellow, Nclr, Red, Green, Purple

# Variables for colored text
Cyan = "\033[96m"
Blue = "\033[94m"
Yellow = "\033[93m"
Nclr = "\033[00m"
Red = "\033[91m"
Green = "\033[92m"
Purple = "\033[95m"

# Old definitions for colored text, delete after all autodocs are merged
def prRed(skk): print("\033[91m{}\033[00m".format(skk))
def prGreen(skk): print("\033[92m{}\033[00m".format(skk))
def prCyan(skk): print("\033[96m{}\033[00m".format(skk))
def prYellow(skk): print("\033[93m{}\033[00m".format(skk))
def prPurple(skk): print("\033[95m{}\033[00m".format(skk))
def prLightPurple(skk): print("\033[94m{}\033[00m".format(skk))

def MY_func(Ls_cont):
    """
    Returns the Mars Year

    :param Ls_cont: solar longitude (continuous)
    :type Ls_cont: array

    :return: the Mars year
    :rtype: int
    """
    return (Ls_cont)//(360.) + 1

def find_tod_in_diurn(fNcdf):
    """
    Returns the variable for the local time axis in diurn files
    (e.g., time_of_day_24).
    Original implementation by Victoria H.

    :param fNcdf: a netCDF file
    :type fNcdf: netCDF file object

    :return: the name of the time of day dimension
    :rtype: str
    """
    regex = re.compile("time_of_day.")
    varset = fNcdf.variables.keys()
    # Extract the 1st element of the list
    return [string for string in varset if re.match(regex, string)][0]



def print_fileContent(fileNcdf):
    """
    Prints the contents of a netCDF file to the screen. Variables sorted
    by dimension.

    :param fileNcdf: full path to the netCDF file
    :type fileNcdf: str

    :return: None
    """
    if not os.path.isfile(fileNcdf.name):
        print(f"{fileNcdf.name} not found")
    else:
        f = Dataset(fileNcdf.name, "r")
        print("==================== DIMENSIONS ====================")
        print(list(f.dimensions.keys()))
        print(str(f.dimensions))
        print("===================== CONTENT =====================")
        # Get all variables
        all_var = f.variables.keys()
        # Initialize empty list
        all_dims = list()
        for ivar in all_var:
            # Get all the variable dimensions
            all_dims.append(f.variables[ivar].dimensions)

        # Filter duplicates. An object of type set() is an unordered
        # collection of distinct objects
        all_dims = set(all_dims)
        # Sort dimensions by length (e.g., (lat) will come before
        # (lat, lon))
        all_dims = sorted(all_dims, key = len)

        for idim in all_dims:
            for ivar in all_var:
                if f.variables[ivar].dimensions == idim:
                    txt_dim = getattr(f.variables[ivar], "dimensions", "")
                    txt_shape = getattr(f.variables[ivar], "shape", "")
                    txt_long_name = getattr(f.variables[ivar], "long_name", "")
                    txt_units = getattr(f.variables[ivar], "units", "")
                    print(
                        f"{Green}{ivar.ljust(15)}: {Purple}{txt_dim}= "
                        f"{Cyan}{txt_shape}, {Yellow}{txt_long_name}  "
                        f"[{txt_units}]{Nclr}")
        try:
            # Skipped if the netCDF file does not contain time
            t_ini = f.variables["time"][0]
            t_end = f.variables["time"][-1]
            Ls_ini = np.squeeze(f.variables["areo"])[0]
            Ls_end = np.squeeze(f.variables["areo"])[-1]
            MY_ini = MY_func(Ls_ini)
            MY_end = MY_func(Ls_end)
            print(f"\nLs ranging from {(np.mod(Ls_ini, 360.)):.2f} to "
                  f"{(np.mod(Ls_end, 360.)):.2f}: {(t_end-t_ini):.2f} days"
                  f"               (MY {MY_ini:02})   (MY {MY_end:02})")
        except:
            pass
        f.close()
        print("=====================================================")

def print_varContent(fileNcdf, list_varfull, print_stat=False):
    """
    Print variable contents from a variable in a netCDF file. Requires
    a XXXXX.fixed.nc file in the current directory.

    :param fileNcdf: full path to a netcdf file
    :type fileNcdf: str
    :param list_varfull: list of variable names and optional slices
        (e.g., ``["lon", "ps[:, 10, 20]"]``)
    :type list_varfull: list
    :param print_stat: If True, print min, mean, and max. If False,
        print values. Defaults to False
    :type print_stat: bool, optional

    :return: None
    """
    if not os.path.isfile(fileNcdf.name):
        print(f"{fileNcdf.name} not found")
    else:
        if print_stat:
            print(f"{Cyan}____________________________________________________"
                  f"______________________\n"
                  f"           VAR            |      MIN      |      MEAN     "
                  f"|      MAX      |\n"
                  f"__________________________|_______________|_______________"
                  f"|_______________|{Nclr}")
        for varfull in list_varfull:
            try:
                slice = "[:]"
                if "[" in varfull:
                    varname, slice = varfull.strip().split("[")
                    slice = f"[{slice}"
                else:
                    varname = varfull.strip()
                cmd_txt = f"f.variables['{varname}']{slice}"
                f = Dataset(fileNcdf.name, "r")
                var = eval(cmd_txt)

                if print_stat:
                    Min = np.nanmin(var)
                    Mean = np.nanmean(var)
                    Max = np.nanmax(var)
                    print(f"{Cyan}{varfull:>26s}|{Min:>15g}|{Mean:>15g}|"
                          f"{Max:>15g}|{Nclr}")
                    if varname == "areo":
                        # If variable is areo then print modulo
                        print(f"{Cyan}{varfull:>17s}(mod 360)|"
                              f"({(np.nanmin(var%360)):>13g})|"
                              f"({(np.nanmean(var%360)):>13g})|"
                              f"({(np.nanmax(var%360)):>13g})|{Nclr}")
                else:
                    if varname != "areo":
                        print(f"{Cyan}{varfull}= {Nclr}")
                        print(f"{Cyan}{var}{Nclr}")
                    else:
                        # Special case for areo then print modulo
                        print(f"{Cyan}areo (areo mod 360)={Nclr}")
                        for ii in var:
                            print(ii, ii%360)

                    print(f"{Cyan}____________________________________________"
                          f"__________________________{Nclr}")
            except:
                if print_stat:
                    print(f"{Red}{varfull:>26s}|{'':>15s}|{'':>15s}|"
                          f"{'':>15s}|")
                else:
                    print(f"{Red}{varfull}")
        if print_stat:
            # Last line for the table
            print(f"{Cyan}__________________________|_______________|_________"
                  f"______|_______________|{Nclr}")
        f.close()

def give_permission(filename):
    """ Sets group file permissions for the NAS system """
    try:
        # catch error and standard output
        subprocess.check_call(["setfacl -v"],
                              shell = True,
                              stdout = open(os.devnull, "w"),
                              stderr = open(os.devnull, "w"))
        cmd_txt = f"setfacl -R -m g:s0846:r {filename}"
        subprocess.call(cmd_txt, shell = True)
    except subprocess.CalledProcessError:
        pass

def check_file_tape(fileNcdf, abort=False):
    """
    For use in the NAS environnment only.
    Checks whether a file is exists on the disk by running the command
    ``dmls -l`` on NAS. This prevents the program from stalling if the
    files need to be migrated from the disk to the tape.

    :param fileNcdf: full path to a netcdf file or a file object with a name attribute
    :type fileNcdf: str or file object
    :param abort: If True, exit the program. Defaults to False
    :type abort: bool, optional

    :return: None
    """
    # Get the filename, whether passed as string or as file object
    filename = fileNcdf if isinstance(fileNcdf, str) else fileNcdf.name
    
    # If filename is not a netCDF file, exit program
    if not re.search(".nc", filename):
        print(f"{Red}{filename} is not a netCDF file{Nclr}")
        exit()
        
    try:
        # Check if the file exists on the disk, exit otherwise. If it
        # exists, copy it over from the disk.
        # Check if dmls command is available
        subprocess.check_call(["dmls"], shell = True,
                              stdout = open(os.devnull, "w"),
                              stderr = open(os.devnull, "w"))
        # Get the last columns of the ls command (filename and status)
        cmd_txt = f"dmls -l {filename}| awk '{{print $8,$9}}'"
        # Get 3-letter identifier from dmls -l command, convert byte to
        # string for Python3
        dmls_out = subprocess.check_output(cmd_txt,
                                           shell = True).decode("utf-8")
        if dmls_out[1:4] not in ["DUL", "REG", "MIG"]:
            # If file is OFFLINE, UNMIGRATING etc...
            if abort :
                print(f"{Red}*** Error ***\n{dmls_out}\n{dmls_out[6:-1]} "
                      f"is not available on disk, status is: {dmls_out[0:5]}"
                      f"CHECK file status with ``dmls -l *.nc`` and run "
                      f"``dmget *.nc`` to migrate the files.\nExiting now..."
                      f"\n{Nclr}")
                exit()
            else:
                print(f"{Yellow}*** Warning ***\n{dmls_out[6:-1]} is not  "
                      f"available on the disk. Status: {dmls_out[0:5]}.\n"
                      f"Consider checking file status with ``dmls -l *.nc`` "
                      f"and run ``dmget *.nc`` to migrate the files.\n"
                      f"Waiting for file to be migrated to the disk, this "
                      f"may take awhile...")
    except subprocess.CalledProcessError:
        # Return an eror message
        if abort:
             exit()
        else:
            pass


def get_Ncdf_path(fNcdf):
    """
    Returns the full path for a netCDF file object.

    .. NOTE:: ``Dataset`` and multi-file dataset (``MFDataset``) have
    different attributes for the path, hence the need for this function.

    :param fNcdf: Dataset or MFDataset object
    :type fNcdf: netCDF file object
    :return: string list for the Dataset (MFDataset)
    :rtype: str(list)
    """
    # Only MFDataset has the_files attribute
    fname_out = getattr(fNcdf, "_files", False)
    # Regular Dataset
    if not fname_out:
        fname_out = getattr(fNcdf, "filepath")()
    return fname_out

def extract_path_basename(filename):
    """
    Returns the path and basename of a file. If only the filename is
    provided, assume it is in the current directory.

    :param filename: name of the netCDF file (may include full path)
    :type filename: str

    :return: full file path & name of file

    .. NOTE:: This routine does not confirm that the file exists.
        It operates on the provided input string.
    """
    # Get the filename without the path
    if ("/" in filename or
        "\\" in filename):
        filepath, basename = os.path.split(filename)
    else:
        filepath = os.getcwd()
        basename = filename

    if "~" in filepath:
        # If the home ("~") symbol is included, expand the user path
        filepath = os.path.expanduser(filepath)
    return filepath, basename

def FV3_file_type(fNcdf):
    """
    Return the type of the netCDF file (i.e., ``fixed``, ``diurn``,
    ``average``, ``daily``) and the format of the Ls array ``areo``
    (i.e., ``fixed``, ``continuous``, or ``diurn``).

    :param fNcdf: an open Netcdf file
    :type fNcdf: Netcdf file object

    :return: The Ls array type (string, ``fixed``, ``continuous``, or
        ``diurn``) and the netCDF file type (string ``fixed``,
        ``diurn``, ``average``, or ``daily``)
    """
    # Get the full path from the file
    fullpath = get_Ncdf_path(fNcdf)

    if type(fullpath) == list:
        # If MFDataset, get the 1st file in the list
        fullpath = fullpath[0]

    # Get the filename without the path
    _, filename = os.path.split(fullpath)

    # Initialize, assume the file is continuous
    f_type = "continuous"
    interp_type = "unknown"
    tod_name = "n/a"

    # model=read_variable_dict_amescap_profile(fNcdf)

    if "time" not in fNcdf.dimensions.keys():
        # If ``time`` is not a dimension, assume it is a fixed file
        f_type = "fixed"
    try:
        tod_name = find_tod_in_diurn(fNcdf)
        if tod_name in fNcdf.dimensions.keys():
            # If ``tod_name_XX`` is a dimension, it is a diurn file
            f_type = "diurn"
    except:
        pass

    dims = fNcdf.dimensions.keys()
    if "pfull" in dims:
        interp_type = "pfull"
    if "pstd"  in dims:
        interp_type = "pstd"
    if "zstd"  in dims:
        interp_type = "zstd"
    if "zagl"  in dims:
        interp_type = "zagl"
    if "zgrid" in dims:
        interp_type = "zgrid"

    return f_type, interp_type

def alt_FV3path(fullpaths, alt, test_exist=True):
    """
    Returns the original or fixed file given an interpolated daily,
    diurn or average file.

    :param fullpaths: full path to a file or a list of full paths to
        more than one file
    :type fullpaths: str
    :param alt: type of file to return (i.e., original or fixed)
    :type alt: str
    :param test_exist: Whether file exists on the disk, defaults to True
    :type test_exist: bool, optional

    :raises ValueError: _description_

    :return: path to original or fixed file
        (e.g., "/u/path/00010.atmos_average.nc" or
        "/u/path/00010.fixed.nc")
    :rtype: str
    """
    out_list = []
    one_element = False

    if type(fullpaths) == str:
        # Convert to a list for generality
        one_element = True
        fullpaths = [fullpaths]

    for fullpath in fullpaths:
        path, filename = os.path.split(fullpath)
        # Get the date
        DDDDD = filename.split(".")[0]
        # Get the extension
        ext = filename[-8:]

        if alt == "raw":
            # This is an interpolated file
            if ext in ["_pstd.nc", "_zstd.nc", "_zagl.nc", "plevs.nc"]:
                if ext == "plevs.nc":
                    file_raw = f"{filename[0:-9]}.nc"
                else:
                    file_raw = f"{filename[0:-8]}.nc"
            else:
                raise ValueError(f"In alt_FV3path(), FV3 file {filename} "
                                 f"not recognized")
            new_full_path = f"{path}/{file_raw}"
        if alt == "fixed":
            new_full_path = f"{path}/{DDDDD}.fixed.nc"
        if (test_exist and not os.path.exists(new_full_path)):
            raise ValueError(f"In alt_FV3path(), {new_full_path} does not "
                             f"exist")

        out_list.append(new_full_path)

    if one_element:
        out_list = out_list[0]

    return out_list

def smart_reader(fNcdf, var_list, suppress_warning=False):
    """
    Alternative to ``var = fNcdf.variables["var"][:]`` for handling
    *processed* files that also checks for a matching average or daily
    and XXXXX.fixed.nc file.

    :param fNcdf: an open netCDF file
    :type fNcdf: netCDF file object
    :param var_list: a variable or list of variables (e.g., ``areo`` or
        [``pk``, ``bk``, ``areo``])
    :type var_list: _type_
    :param suppress_warning: suppress debug statement. Useful if a
        variable is not expected to be in the file anyway. Defaults to
        False
    :type suppress_warning: bool, optional

    :return: variable content (single or values to unpack)
    :rtype: list

    Example::

        from netCDF4 import Dataset

        fNcdf = Dataset("/u/akling/FV3/00668.atmos_average_pstd.nc", "r")

        # Approach using var = fNcdf.variables["var"][:]
        ucomp = fNcdf.variables["ucomp"][:]
        # New approach that checks for matching average/daily & fixed
        vcomp = smart_reader(fNcdf, "vcomp")

        # This will pull "areo" from an original file if it is not
        # available in the interpolated file. If pk and bk are also not
        # in the average file, it will check for them in the fixed file.
        pk, bk, areo = smart_reader(fNcdf, ["pk", "bk", "areo"])

    .. NOTE:: Only the variable content is returned, not attributes
    """
    # This out_list is for the variable
    out_list = []
    one_element = False
    file_is_MF = False

    # Return string (Dataset) or list (MFDataset)
    Ncdf_path = get_Ncdf_path(fNcdf)
    if type(Ncdf_path) == list:
        file_is_MF = True

    # Convert to list for generality if only one variable is provided
    # (e.g., areo -> [areo])
    if type(var_list) == str:
        one_element = True
        var_list = [var_list]

    for ivar in var_list:
        if ivar in fNcdf.variables.keys():
            # Try to read in the original file
            out_list.append(fNcdf.variables[ivar][:])
        else:
            full_path_try = alt_FV3path(Ncdf_path, alt = "raw",
                                        test_exist = True)

            if file_is_MF:
                f_tmp = MFDataset(full_path_try, "r")
            else:
                f_tmp = Dataset(full_path_try, "r")

            if ivar in f_tmp.variables.keys():
                out_list.append(f_tmp.variables[ivar][:])
                if not suppress_warning:
                    print(f"**Warning*** Using variable {ivar} in "
                          f"{full_path_try} instead of original file(s)")
                f_tmp.close()
            else:
                f_tmp.close()
                full_path_try = alt_FV3path(Ncdf_path,
                                            alt = "fixed",
                                            test_exist = True)

                if file_is_MF:
                    full_path_try = full_path_try[0]

                f_tmp = Dataset(full_path_try, "r")

                if ivar in f_tmp.variables.keys():
                    out_list.append(f_tmp.variables[ivar][:])
                    f_tmp.close()
                    if not suppress_warning:
                        print(f"**Warning*** Using variable {ivar} in "
                              f"{full_path_try} instead of original file(s)")
                else:
                    print(f"***ERROR*** Variable {ivar} not found in "
                          f"{full_path_try}, NOR in raw output or fixed file\n"
                          f"            >>> Assigning {ivar} to NaN")
                    f_tmp.close()
                    out_list.append(np.nan)

    if one_element:
        out_list = out_list[0]
    return out_list

def regrid_Ncfile(VAR_Ncdf, file_Nc_in, file_Nc_target):
    """
    Regrid a netCDF variable from one file structure to another.
    Requires a file with the desired file structure to mimic.
    [Alex Kling, May 2021]

    :param VAR_Ncdf: a netCDF variable object to regrid
        (e.g., ``f_in.variable["temp"]``)
    :type VAR_Ncdf: netCDF file variable
    :param file_Nc_in: an open netCDF file to source for the variable
        (e.g., ``f_in = Dataset("filename", "r")``)
    :type file_Nc_in: netCDF file object
    :param file_Nc_target: an open netCDF file with the desired file
        structure (e.g., ``f_out = Dataset("filename", "r")``)
    :type file_Nc_target: netCDF file object

    :return: the values of the variable interpolated to the target file
        grid.
    :rtype: array

    .. NOTE:: While the KDTree interpolation can handle a 3D dataset
    (lon/lat/lev instead of just 2D lon/lat), the grid points in the
    vertical are just a few (10--100s) meters in the PBL vs a few
    (10-100s) kilometers in the horizontal. This results in excessive
    weighting in the vertical, which is why the vertical dimension is
    handled separately.
    """
    from amescap.FV3_utils import interp_KDTree, axis_interp
    ftype_in, zaxis_in = FV3_file_type(file_Nc_in)
    ftype_t, zaxis_t = FV3_file_type(file_Nc_target)

    # Sanity check
    if ftype_in != ftype_t:
        print(f"*** Warning*** in regrid_Ncfile, input file {ftype_in} "
              f"and target file {ftype_t} must have the same type")

    if zaxis_in != zaxis_t:
        print(f"*** Warning*** in regrid_Ncfile, input file {zaxis_in} "
              f"and target file {zaxis_t} must have the same vertical grid")

    if zaxis_in == "pfull" or zaxis_t == "pfull":
        print(f"*** Warning*** in regrid_Ncfile, input file {zaxis_in} "
              f"and target file {zaxis_t} must be vertically interpolated")

    # Get target dimensions
    lon_t = file_Nc_target.variables["lon"][:]
    lat_t = file_Nc_target.variables["lat"][:]
    if "time" in VAR_Ncdf.dimensions:
        areo_t = file_Nc_target.variables["areo"][:]
        time_t = file_Nc_target.variables["time"][:]

    # Get input dimensions
    lon_in = file_Nc_in.variables["lon"][:]
    lat_in = file_Nc_in.variables["lat"][:]
    if "time" in VAR_Ncdf.dimensions:
        areo_in = file_Nc_in.variables["areo"][:]
        time_in = file_Nc_in.variables["time"][:]

    # Get array elements
    var_OUT = VAR_Ncdf[:]


    if not (np.array_equal(lat_in,lat_t) and np.array_equal(lon_in,lon_t)):
        # STEP 1: lat/lon interpolation are always performed unless
        # target lon and lat are identical
        if len(np.atleast_1d(lon_in)) == 1:
            # Special case if input len(lon)=1 (slice or zonal avg),
            # only interpolate on the latitude axis
            var_OUT = axis_interp(var_OUT, lat_in, lat_t, axis = -2,
                                  reverse_input = False, type_int = "lin")
        elif len(np.atleast_1d(lat_in)) == 1:
            # Special case if input len(lat)=1 (slice or medidional avg),
            # only interpolate on the longitude axis
            var_OUT = axis_interp(var_OUT, lon_in, lon_t, axis = -1,
                                  reverse_input = False, type_int = "lin")
        else:
            # Bi-directional interpolation
            var_OUT = interp_KDTree(var_OUT, lat_in, lon_in, lat_t, lon_t)

    if zaxis_in in VAR_Ncdf.dimensions:
        # STEP 2: linear or log interpolation. If there is a vertical
        # axis, get position: pstd is 1 in (time, pstd, lat, lon)
        pos_axis = VAR_Ncdf.dimensions.index(zaxis_in)
        lev_in = file_Nc_in.variables[zaxis_in][:]
        lev_t = file_Nc_target.variables[zaxis_t][:]

        # Check if the input needs to be reversed. Reuses find_n(),
        # which was designed for pressure interpolation so the values
        # are reversed if up = increasing
        if lev_in[0] < lev_in[-1]:
            reverse_input = True
        else:
            reverse_input = False
        if zaxis_in in ["zagl", "zstd"]:
            intType = "lin"
        elif zaxis_in == "pstd":
            intType = "log"
        var_OUT = axis_interp(var_OUT, lev_in,lev_t, pos_axis,
                              reverse_input = reverse_input,
                              type_int = intType)

    if "time" in VAR_Ncdf.dimensions:
        # STEP 3: Linear interpolation in Ls
        pos_axis = 0
        var_OUT = axis_interp(var_OUT, np.squeeze(areo_in)%360,
                              np.squeeze(areo_t)%360, pos_axis,
                              reverse_input = False, type_int = "lin")

    if ftype_in == "diurn":
        # STEP 4: Linear interpolation in time of day
        # TODO: (the interpolation scheme is not cyclic. If available
        # diurn times are 04 10 16 22 and requested time is 23, value
        # is set to zero, not interpololated from 22 and 04
        pos_axis = 1
        tod_name_in = find_tod_in_diurn(file_Nc_in)
        tod_name_t = find_tod_in_diurn(file_Nc_target)
        tod_in = file_Nc_in.variables[tod_name_in][:]
        tod_t = file_Nc_target.variables[tod_name_t][:]
        var_OUT = axis_interp(var_OUT, tod_in, tod_t, pos_axis,
                              reverse_input = False,
                              type_int = "lin")
    return var_OUT

def progress(k, Nmax):
    """
    Displays a progress bar to monitor heavy calculations.

    :param k: current iteration of the outer loop
    :type k: int
    :param Nmax: max iteration of the outer loop
    :type Nmax: int
    """
    # For rounding to the 2nd digit
    from math import ceil
    progress = float(k)/Nmax
    # Modify barLength to change length of progress bar
    barLength = 10
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "Error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength * progress))
    text = (f"Running... [{'#'*block + '-'*(barLength-block)}] "
            f"{ceil(progress*100*100)/100} {status}%")
    sys.stdout.write(text)
    sys.stdout.flush()

def section_content_amescap_profile(section_ID):
    """
    Executes first code section in ``~/.amescap_profile`` to read in
    user-defined plot & interpolation settings.

    :param section_ID: the section to load (e.g., Pressure definitions
        for pstd)
    :type section_ID: str

    :return: the relevant line with Python syntax
    """
    import os
    import numpy as np
    input_file = os.environ["HOME"]+"/.amescap_profile"
    try:
        f = open(input_file, "r")
        contents = ""
        rec = False
        while True:
            line = f.readline()
            if not line:
                # End of File
                break
            if line[0] == "<":
                rec = False
                if line.split("|")[1].strip() == section_ID:
                    rec = True
                    line = f.readline()
            if rec:
                contents += line
        f.close()
        if contents == "":
            print(f"{Red}No content found for <<< {section_ID} >>> block")
        return contents

    except FileNotFoundError:
        print(f"{Red}Error: {input_file} config file not found.\n"
              f"{Yellow}To use this feature, create a hidden config "
              f"file from the template in your home directory with:\n"
              f"{Cyan}    ``cp AmesCAP/mars_templates/amescap_profile  "
              f"~/.amescap_profile``")
        exit()
    except Exception as exception:
        # Return the error
        print(f"{Red}Error")
        print(exception)

def filter_vars(fNcdf, include_list=None, giveExclude=False):
    """
    Filters the variable names in a netCDF file for processing. Returns
    all dimensions (``lon``, ``lat``, etc.), the ``areo`` variable, and
    any other variable listed in ``include_list``.

    :param fNcdf: an open netCDF object for a diurn, daily, or average
        file
    :type fNcdf: netCDF file object
    :param include_list:list of variables to include (e.g., [``ucomp``,
        ``vcomp``], defaults to None
    :type include_list: list or None, optional
    :param giveExclude: if True, returns variables to be excluded from
        the file, defaults to False
    :type giveExclude: bool, optional
    :return: list of variable names to include in the processed file
    """
    var_list = fNcdf.variables.keys()
    if include_list is None:
        # If no list is provided, return all variables:
        return var_list

    input_list_filtered = []
    for ivar in include_list:
        if ivar in var_list:
            # Make sure the requested variables are present in file
            input_list_filtered.append(ivar)
        else:
            print(f"{Yellow}***Warning*** In filter_vars(), variable(s)"
                  f"{ivar} not found in file")

    baseline_var = []
    for ivar in  var_list:
        # Compute baseline variables, i.e. all dimensions, axis etc...
        if (ivar == "areo" or
            len(fNcdf.variables[ivar].dimensions) <= 2):
            baseline_var.append(ivar)

    out_list = baseline_var + input_list_filtered

    if giveExclude:
        # Return the two lists
        exclude_list = list(var_list)
        for ivar in out_list:
            exclude_list.remove(ivar)
        out_list = exclude_list
    return out_list

def find_fixedfile(filename):
    """
    Finds the relevant fixed file for a given average, daily, or diurn
    file.
    [Batterson, Updated by Alex Nov 29 2022]

    :param filename: an average, daily, or diurn netCDF file
    :type filename: str
    :return: full path to the correspnding fixed file
    :rtype: str

    Compatible file types::

            DDDDD.atmos_average.nc                  -> DDDDD.fixed.nc
            atmos_average.tileX.nc                  -> fixed.tileX.nc
            DDDDD.atmos_average_plevs.nc            -> DDDDD.fixed.nc
            DDDDD.atmos_average_plevs_custom.nc     -> DDDDD.fixed.nc
            atmos_average.tileX_plevs.nc            -> fixed.tileX.nc
            atmos_average.tileX_plevs_custom.nc     -> fixed.tileX.nc
            atmos_average_custom.tileX_plevs.nc     -> fixed.tileX.nc

    """
    filepath, fname = extract_path_basename(filename)
    if "tile" in fname:
        # Try the tile or standard version of the fixed files
        name_fixed = f"{filepath}/fixed.tile{fname.split('tile')[1][0]}.nc"
    else:
        name_fixed = f"{filepath}/{fname.split('.')[0]}.fixed.nc"

    if not  os.path.exists(name_fixed):
        # If neither is found, set-up a default name
        name_fixed = "FixedFileNotFound"
    return name_fixed

def get_longname_unit(fNcdf, varname):
    """
    Returns the longname and unit attributes of a variable in a netCDF
    file. If the attributes are unavailable, returns blank strings to
    avoid an error.

    :param fNcdf: an open netCDF file
    :type fNcdf: netCDF file object
    :param varname: variable to extract attribute from
    :type varname: str

    :return: longname and unit attributes
    :rtype: str

    .. NOTE:: Some functions in MarsVars edit the units
    (e.g., [kg] -> [kg/m]), therefore the empty string is 4 characters
    in length ("    " instead of "") to allow for editing by
    ``editing units_txt[:-2]``, for example.
    """
    return (getattr(fNcdf.variables[varname], "long_name", "    "),
            getattr(fNcdf.variables[varname], "units", "    "))

def wbr_cmap():
    """
    Returns a color map that goes from
    white -> blue -> green -> yellow -> red
    """
    from matplotlib.colors import ListedColormap
    tmp_cmap = np.zeros((254, 4))
    tmp_cmap[:, 3] = 1.
    tmp_cmap[:, 0:3] = np.array([
        [255,255,255], [252,254,255], [250,253,255], [247,252,254],
        [244,251,254], [242,250,254], [239,249,254], [236,248,253],
        [234,247,253], [231,246,253], [229,245,253], [226,244,253],
        [223,243,252], [221,242,252], [218,241,252], [215,240,252],
        [213,239,252], [210,238,251], [207,237,251], [205,236,251],
        [202,235,251], [199,234,250], [197,233,250], [194,232,250],
        [191,231,250], [189,230,250], [186,229,249], [183,228,249],
        [181,227,249], [178,226,249], [176,225,249], [173,224,248],
        [170,223,248], [168,222,248], [165,221,248], [162,220,247],
        [157,218,247], [155,216,246], [152,214,245], [150,212,243],
        [148,210,242], [146,208,241], [143,206,240], [141,204,238],
        [139,202,237], [136,200,236], [134,197,235], [132,195,234],
        [129,193,232], [127,191,231], [125,189,230], [123,187,229],
        [120,185,228], [118,183,226], [116,181,225], [113,179,224],
        [111,177,223], [109,175,221], [106,173,220], [104,171,219],
        [102,169,218], [100,167,217], [ 97,165,215], [ 95,163,214],
        [ 93,160,213], [ 90,158,212], [ 88,156,211], [ 86,154,209],
        [ 83,152,208], [ 81,150,207], [ 79,148,206], [ 77,146,204],
        [ 72,142,202], [ 72,143,198], [ 72,144,195], [ 72,145,191],
        [ 72,146,188], [ 72,147,184], [ 72,148,181], [ 72,149,177],
        [ 72,150,173], [ 72,151,170], [ 72,153,166], [ 72,154,163],
        [ 72,155,159], [ 72,156,156], [ 72,157,152], [ 72,158,148],
        [ 72,159,145], [ 72,160,141], [ 72,161,138], [ 73,162,134],
        [ 73,163,131], [ 73,164,127], [ 73,165,124], [ 73,166,120],
        [ 73,167,116], [ 73,168,113], [ 73,169,109], [ 73,170,106],
        [ 73,172,102], [ 73,173, 99], [ 73,174, 95], [ 73,175, 91],
        [ 73,176, 88], [ 73,177, 84], [ 73,178, 81], [ 73,179, 77],
        [ 73,181, 70], [ 78,182, 71], [ 83,184, 71], [ 87,185, 72],
        [ 92,187, 72], [ 97,188, 73], [102,189, 74], [106,191, 74],
        [111,192, 75], [116,193, 75], [121,195, 76], [126,196, 77],
        [130,198, 77], [135,199, 78], [140,200, 78], [145,202, 79],
        [150,203, 80], [154,204, 80], [159,206, 81], [164,207, 81],
        [169,209, 82], [173,210, 82], [178,211, 83], [183,213, 84],
        [188,214, 84], [193,215, 85], [197,217, 85], [202,218, 86],
        [207,220, 87], [212,221, 87], [217,222, 88], [221,224, 88],
        [226,225, 89], [231,226, 90], [236,228, 90], [240,229, 91],
        [245,231, 91], [250,232, 92], [250,229, 91], [250,225, 89],
        [250,222, 88], [249,218, 86], [249,215, 85], [249,212, 84],
        [249,208, 82], [249,205, 81], [249,201, 80], [249,198, 78],
        [249,195, 77], [248,191, 75], [248,188, 74], [248,184, 73],
        [248,181, 71], [248,178, 70], [248,174, 69], [248,171, 67],
        [247,167, 66], [247,164, 64], [247,160, 63], [247,157, 62],
        [247,154, 60], [247,150, 59], [247,147, 58], [246,143, 56],
        [246,140, 55], [246,137, 53], [246,133, 52], [246,130, 51],
        [246,126, 49], [246,123, 48], [246,120, 47], [245,116, 45],
        [245,113, 44], [245,106, 41], [244,104, 41], [243,102, 41],
        [242,100, 41], [241, 98, 41], [240, 96, 41], [239, 94, 41],
        [239, 92, 41], [238, 90, 41], [237, 88, 41], [236, 86, 41],
        [235, 84, 41], [234, 82, 41], [233, 80, 41], [232, 78, 41],
        [231, 76, 41], [230, 74, 41], [229, 72, 41], [228, 70, 41],
        [228, 67, 40], [227, 65, 40], [226, 63, 40], [225, 61, 40],
        [224, 59, 40], [223, 57, 40], [222, 55, 40], [221, 53, 40],
        [220, 51, 40], [219, 49, 40], [218, 47, 40], [217, 45, 40],
        [217, 43, 40], [216, 41, 40], [215, 39, 40], [214, 37, 40],
        [213, 35, 40], [211, 31, 40], [209, 31, 40], [207, 30, 39],
        [206, 30, 39], [204, 30, 38], [202, 30, 38], [200, 29, 38],
        [199, 29, 37], [197, 29, 37], [195, 29, 36], [193, 28, 36],
        [192, 28, 36], [190, 28, 35], [188, 27, 35], [186, 27, 34],
        [185, 27, 34], [183, 27, 34], [181, 26, 33], [179, 26, 33],
        [178, 26, 32], [176, 26, 32], [174, 25, 31], [172, 25, 31],
        [171, 25, 31], [169, 25, 30], [167, 24, 30], [165, 24, 29],
        [164, 24, 29], [162, 23, 29], [160, 23, 28], [158, 23, 28],
        [157, 23, 27], [155, 22, 27], [153, 22, 27], [151, 22, 26],
        [150, 22, 26], [146, 21, 25]])/255.
    return ListedColormap(tmp_cmap)

def rjw_cmap():
    """
    Returns John Wilson's preferred color map
    (red -> jade -> wisteria)
    """
    from matplotlib.colors import ListedColormap
    tmp_cmap = np.zeros((55, 4))
    tmp_cmap[:, 3] = 1.
    tmp_cmap[:, 0:3] = np.array([
        [255,   0, 244], [248,  40, 244], [241,  79, 244], [234, 119, 244],
        [228, 158, 244], [221, 197, 245], [214, 190, 245], [208, 182, 245],
        [201, 175, 245], [194, 167, 245], [188, 160, 245], [181, 152, 246],
        [175, 145, 246], [140, 140, 247], [105, 134, 249], [ 70, 129, 251],
        [ 35, 124, 253], [  0, 119, 255], [  0, 146, 250], [  0, 173, 245],
        [  0, 200, 241], [  0, 227, 236], [  0, 255, 231], [  0, 255, 185],
        [  0, 255, 139], [  0, 255,  92], [  0, 255,  46], [  0, 255,   0],
        [ 63, 247,  43], [127, 240,  87], [191, 232, 130], [255, 225, 174],
        [255, 231, 139], [255, 237, 104], [255, 243,  69], [255, 249,  34],
        [255, 255,   0], [255, 241,  11], [255, 227,  23], [255, 213,  35],
        [255, 199,  47], [255, 186,  59], [255, 172,  71], [255, 160,  69],
        [255, 148,  67], [255, 136,  64], [255, 124,  62], [255, 112,  60],
        [255, 100,  58], [255,  80,  46], [255,  60,  34], [255,  40,  23],
        [255,  20,  11], [255,   0,   0], [237,  17,   0]])/255.
    return ListedColormap(tmp_cmap)

def hot_cold_cmap():
    """
    Returns Dark blue > light blue>white>yellow>red colormap
    Based on Matlab's bipolar colormap
    """
    from matplotlib.colors import ListedColormap
    tmp_cmap = np.zeros((128,4))
    tmp_cmap [:,3]=1. #set alpha
    tmp_cmap[:,0:3]=np.array([
    [0,0,255],[0,7,255],[0,15,255],[0,23,255],
    [0,30,255],[1,38,255],[2,45,255],[3,52,255],
    [4,60,255],[5,67,255],[6,73,255],[7,80,255],
    [9,87,255],[10,93,255],[12,99,255],[14,105,255],
    [16,112,255],[18,118,255],[20,124,255],[22,129,255],
    [25,135,255],[27,140,255],[30,145,255],[33,151,255],
    [36,156,255],[39,161,255],[42,166,255],[46,170,255],
    [49,175,255],[53,179,255],[56,183,255],[60,188,255],
    [64,192,255],[68,196,255],[73,199,255],[77,203,255],
    [81,207,255],[86,210,255],[91,213,255],[96,216,255],
    [101,220,255],[106,223,255],[111,225,255],[116,228,255],
    [122,230,255],[127,233,255],[133,235,255],[139,237,255],
    [145,239,255],[151,241,255],[158,243,255],[164,245,255],
    [170,246,255],[177,248,255],[184,249,255],[191,250,255],
    [198,251,255],[205,252,255],[212,253,255],[220,253,255],
    [227,254,255],[235,254,255],[242,254,255],[250,254,255],
    [255,254,250],[255,254,242],[255,254,235],[255,254,227],
    [255,253,220],[255,253,212],[255,252,205],[255,251,198],
    [255,250,191],[255,249,184],[255,248,177],[255,246,170],
    [255,245,164],[255,243,158],[255,241,151],[255,239,145],
    [255,237,139],[255,235,133],[255,233,127],[255,230,122],
    [255,228,116],[255,225,111],[255,223,106],[255,220,101],
    [255,216,96],[255,213,91],[255,210,86],[255,207,81],
    [255,203,77],[255,199,73],[255,196,68],[255,192,64],
    [255,188,60],[255,183,56],[255,179,53],[255,175,49],
    [255,170,46],[255,166,42],[255,161,39],[255,156,36],
    [255,151,33],[255,145,30],[255,140,27],[255,135,25],
    [255,129,22],[255,124,20],[255,118,18],[255,112,16],
    [255,105,14],[255,99,12],[255,93,10],[255,87,9],
    [255,80,7],[255,73,6],[255,67,5],[255,60,4],
    [255,52,3],[255,45,2],[255,38,1],[255,30,0],
    [255,23,0],[255,15,0],[255,7,0],[255,0,0]])/255

    return ListedColormap(tmp_cmap)

def dkass_dust_cmap():
    """
    Returns a color map useful for dust cross-sections.
    (yellow -> orange -> red -> purple)
    Provided by Courtney Batterson.
    """
    from matplotlib.colors import ListedColormap, hex2color
    tmp_cmap = np.zeros((256, 4))
    tmp_cmap[:, 3] = 1.
    dkass_cmap = [
        "#ffffa3", "#fffea1", "#fffc9f", "#fffa9d", "#fff99b", "#fff799",
        "#fff597", "#fef395", "#fef293", "#fef091", "#feee8f", "#feec8d",
        "#fdea8b", "#fde989", "#fde787", "#fde584", "#fce382", "#fce180",
        "#fce07e", "#fcde7c", "#fcdc7a", "#fbda78", "#fbd976", "#fbd774",
        "#fbd572", "#fbd370", "#fad16e", "#fad06c", "#face6a", "#facc68",
        "#faca66", "#f9c964", "#f9c762", "#f9c560", "#f9c35e", "#f9c15c",
        "#f8c05a", "#f8be58", "#f8bc56", "#f8ba54", "#f8b952", "#f7b750",
        "#f7b54e", "#f7b34b", "#f7b149", "#f7b047", "#f6ae45", "#f6ac43",
        "#f6aa41", "#f6a83f", "#f5a73d", "#f5a53b", "#f5a339", "#f5a137",
        "#f5a035", "#f49e33", "#f49c31", "#f49a2f", "#f4982d", "#f4972b",
        "#f39529", "#f39327", "#f39125", "#f39023", "#f38e21", "#f28b22",
        "#f28923", "#f18724", "#f18524", "#f18225", "#f08026", "#f07e27",
        "#f07c27", "#ef7a28", "#ef7729", "#ee752a", "#ee732a", "#ee712b",
        "#ed6e2c", "#ed6c2d", "#ed6a2d", "#ec682e", "#ec652f", "#eb6330",
        "#eb6130", "#eb5f31", "#ea5d32", "#ea5a33", "#ea5833", "#e95634",
        "#e95435", "#e85136", "#e84f36", "#e84d37", "#e74b38", "#e74839",
        "#e74639", "#e6443a", "#e6423b", "#e5403c", "#e53d3c", "#e53b3d",
        "#e4393e", "#e4373f", "#e4343f", "#e33240", "#e33041", "#e22e42",
        "#e22b42", "#e22943", "#e12744", "#e12545", "#e12345", "#e02046",
        "#e01e47", "#df1c48", "#df1a48", "#df1749", "#de154a", "#de134b",
        "#de114c", "#dd0e4c", "#dd0c4d", "#dc0a4e", "#dc084f", "#dc064f",
        "#db0350", "#db0151", "#da0052", "#d90153", "#d70154", "#d60256",
        "#d40257", "#d30258", "#d2035a", "#d0035b", "#cf045c", "#cd045e",
        "#cc055f", "#cb0560", "#c90562", "#c80663", "#c60664", "#c50766",
        "#c30767", "#c20868", "#c1086a", "#bf096b", "#be096c", "#bc096e",
        "#bb0a6f", "#ba0a70", "#b80b72", "#b70b73", "#b50c74", "#b40c75",
        "#b30c77", "#b10d78", "#b00d79", "#ae0e7b", "#ad0e7c", "#ac0f7d",
        "#aa0f7f", "#a90f80", "#a71081", "#a61083", "#a51184", "#a31185",
        "#a21287", "#a01288", "#9f1389", "#9e138b", "#9c138c", "#9b148d",
        "#99148f", "#981590", "#961591", "#951693", "#941694", "#921695",
        "#911797", "#8f1798", "#8e1899", "#8d189a", "#8b199c", "#8a199d",
        "#881a9e", "#871aa0", "#861aa1", "#841ba2", "#831ba4", "#811ca5",
        "#801ca4", "#7f1ba2", "#7f1ba0", "#7e1b9e", "#7d1a9b", "#7c1a99",
        "#7b1a97", "#7a1995", "#791993", "#781991", "#77198f", "#76188d",
        "#75188b", "#751889", "#741786", "#731784", "#721782", "#711680",
        "#70167e", "#6f167c", "#6e167a", "#6d1578", "#6c1576", "#6b1574",
        "#6b1471", "#6a146f", "#69146d", "#68136b", "#671369", "#661367",
        "#651265", "#641263", "#631261", "#62125f", "#61115c", "#61115a",
        "#601158", "#5f1056", "#5e1054", "#5d1052", "#5c0f50", "#5b0f4e",
        "#5a0f4c", "#590f4a", "#580e48", "#570e45", "#570e43", "#560d41",
        "#550d3f", "#540d3d", "#530c3b", "#520c39", "#510c37", "#500c35",
        "#4f0b33", "#4e0b30", "#4d0b2e", "#4d0a2c", "#4c0a2a", "#4b0a28",
        "#4a0926", "#490924", "#480922", "#470820"]
    RGB_T = np.array([hex2color(x) for x in dkass_cmap])
    tmp_cmap[:, 0:3] = RGB_T
    return ListedColormap(tmp_cmap)

def dkass_temp_cmap():
    """
    Returns a color map that highlights the 200K temperatures.
    (black -> purple -> blue -> green -> yellow -> orange -> red)
    Provided by Courtney Batterson.
    """
    from matplotlib.colors import ListedColormap, hex2color
    tmp_cmap = np.zeros((256, 4))
    tmp_cmap[:, 3] = 1.
    dkass_cmap = [
        "#200000", "#230104", "#250208", "#27040c", "#290510", "#2c0614",
        "#2e0718", "#30081c", "#320a20", "#350b24", "#370c28", "#390d2c",
        "#3c0f30", "#3e1034", "#401138", "#42123c", "#451340", "#471544",
        "#491648", "#4b174c", "#4e1850", "#501954", "#521b58", "#541c5c",
        "#571d60", "#591e64", "#5b2068", "#5d216c", "#602270", "#622374",
        "#642478", "#66267c", "#692780", "#6b2884", "#6d2988", "#6f2a8c",
        "#722c90", "#742d94", "#762e98", "#782f9c", "#7b30a0", "#7d32a4",
        "#7f33a9", "#8134ad", "#8435b1", "#8637b5", "#8838b9", "#8a39bd",
        "#8d3ac1", "#8f3bc5", "#913dc9", "#933ecd", "#963fd1", "#9840d5",
        "#9a41d9", "#9c43dd", "#9f44e1", "#a145e5", "#a346e9", "#a548ed",
        "#a849f1", "#aa4af5", "#ac4bf9", "#ae4cfd", "#af4eff", "#ad50ff",
        "#aa53ff", "#a755ff", "#a458ff", "#a25aff", "#9f5cff", "#9c5fff",
        "#9961ff", "#9764ff", "#9466ff", "#9168ff", "#8e6bff", "#8c6dff",
        "#8970ff", "#8672ff", "#8374ff", "#8177ff", "#7e79ff", "#7b7cff",
        "#787eff", "#7581ff", "#7383ff", "#7085ff", "#6d88ff", "#6a8aff",
        "#688dff", "#658fff", "#6291ff", "#5f94ff", "#5d96ff", "#5a99ff",
        "#579bff", "#549dff", "#52a0ff", "#4fa2ff", "#4ca5ff", "#49a7ff",
        "#46aaff", "#44acff", "#41aeff", "#3eb1ff", "#3bb3ff", "#39b6ff",
        "#36b8ff", "#33baff", "#30bdff", "#2ebfff", "#2bc2ff", "#28c4ff",
        "#25c6ff", "#23c9ff", "#20cbff", "#1dceff", "#1ad0ff", "#17d3ff",
        "#15d5ff", "#12d7ff", "#0fdaff", "#0cdcff", "#0adfff", "#07e1ff",
        "#04e3ff", "#01e6ff", "#02e7fe", "#06e8fa", "#0ae8f6", "#0ee8f2",
        "#12e9ee", "#16e9ea", "#1ae9e6", "#1eeae2", "#22eade", "#26ebda",
        "#2aebd6", "#2eebd2", "#32ecce", "#36ecca", "#3aecc6", "#3eedc2",
        "#42edbe", "#46eeba", "#4aeeb6", "#4eeeb2", "#52efae", "#55efaa",
        "#59f0a6", "#5df0a1", "#61f09d", "#65f199", "#69f195", "#6df191",
        "#71f28d", "#75f289", "#79f385", "#7df381", "#81f37d", "#85f479",
        "#89f475", "#8df471", "#91f56d", "#95f569", "#99f665", "#9df661",
        "#a1f65d", "#a5f759", "#a9f755", "#adf751", "#b1f84d", "#b5f849",
        "#b9f945", "#bdf941", "#c1f93d", "#c5fa39", "#c9fa35", "#cdfa31",
        "#d1fb2d", "#d5fb29", "#d9fc25", "#ddfc21", "#e1fc1d", "#e5fd19",
        "#e9fd15", "#edfd11", "#f1fe0d", "#f5fe09", "#f8ff05", "#fcff01",
        "#fdfc00", "#fdf800", "#fef400", "#fef000", "#feec00", "#fee800",
        "#fee400", "#fee000", "#fedc00", "#fed800", "#fed400", "#fed000",
        "#fecc00", "#fec800", "#fec400", "#fec000", "#febc00", "#feb800",
        "#feb400", "#feb000", "#feac00", "#fea800", "#fea400", "#fea000",
        "#fe9c00", "#fe9800", "#fe9400", "#fe9000", "#fe8c00", "#fe8800",
        "#fe8400", "#fe8000", "#fe7c00", "#fe7800", "#fe7400", "#fe7000",
        "#fe6c00", "#fe6800", "#fe6400", "#fe6000", "#fe5c00", "#fe5800",
        "#fe5400", "#ff5000", "#ff4c00", "#ff4800", "#ff4400", "#ff4000",
        "#ff3c00", "#ff3800", "#ff3400", "#ff3000", "#ff2c00", "#ff2800",
        "#ff2400", "#ff2000", "#ff1c00", "#ff1800", "#ff1400", "#ff1000",
        "#ff0c00", "#ff0800", "#ff0400", "#ff0000"]
    RGB_T = np.array([hex2color(x) for x in dkass_cmap])
    tmp_cmap[:, 0:3] = RGB_T
    return ListedColormap(tmp_cmap)

def pretty_print_to_fv_eta(var, varname, nperline=6):
    """
    Print the ``ak`` or ``bk`` coefficients for copying to
    ``fv_eta.f90``.

    :param var: ak or bk data
    :type var: array
    :param varname: the variable name ("a" or "b")
    :type varname: str
    :param nperline: the number of elements per line, defaults to 6
    :type nperline: int, optional

    :return: a print statement for copying into ``fv_eta.f90``
    """
    NLAY = len(var) - 1
    ni = 0
    # Print the piece of code to copy/paste in ``fv_eta.f90``

    if varname == "a":
        # If == a, print the variable definitions before the content
        print(f"\n      real a{int(NLAY)}({int(NLAY+1)}),b{int(NLAY)}"
              f"({int(NLAY+1)})\n")

    # Initialize the first line
    print(f"data {varname}{NLAY} /      &")
    sys.stdout.write("    ")

    for i in range(0, len(var)-1):
         # Loop over all elements
        sys.stdout.write(f"{var[i]:16.10e}")
        ni += 1
        if ni == nperline:
            sys.stdout.write(" &\n    ")
            ni = 0
    sys.stdout.write(f"{var[NLAY]:16.10e} /\n")

    if varname == "b":
        # If b, print the code snippet after the displaying the variable
        ks = 0
        while var[ks] == 0.:
            ks += 1
        print("")

        # Remove 1 because it takes 2 boundary points to form 1 layer
        print(f"        case ({int(NLAY)})")
        print(f"          ks = {int(ks-1)}")
        print(f"          do k=1,km+1")
        print(f"            ak(k) = a{int(NLAY)}(k)")
        print(f"            bk(k) = b{int(NLAY)}(k)")
        print(f"          enddo  ")
        print(f" ")


def replace_dims(Ncvar_dim, vert_dim_name=None):
    """
    Updates the name of the variable dimension to match the format of
    the new NASA Ames Mars GCM output files.

    :param Ncvar_dim: netCDF variable dimensions
        (e.g., ``f_Ncdf.variables["temp"].dimensions``)
    :type Ncvar_dim: str
    :param vert_dim_name: the vertical dimension if it is ambiguous
        (``pstd``, ``zstd``, or ``zagl``). Defaults to None
    :type vert_dim_name: str, optional
    :return: updated dimensions
    :rtype: str
    """
    # Set input dictionary options recognizable as MGCM variables
    lat_dic = ["lat", "lats", "latitudes", "latitude"]
    lon_dic = ["lon", "lon", "longitude", "longitudes"]
    lev_dic = ["pressure", "altitude"]
    areo_dic = ["ls"]

    # Set the desired output names
    dims_out = list(Ncvar_dim).copy()
    for ii, idim in enumerate(Ncvar_dim):
        # Rename axes
        if idim in lat_dic: dims_out[ii] = "lat"
        if idim in lon_dic: dims_out[ii] = "lon"
        if idim in lev_dic:
            if vert_dim_name is None:
                # Vertical coordinate: If no input provided, assume it
                # is standard pressure ``pstd``
                dims_out[ii] = "pstd"
            else:
                # Use provided dimension
                dims_out[ii] = vert_dim_name
    return tuple(dims_out)

def ak_bk_loader(fNcdf):
    """
    Return ``ak`` and ``bk`` arrays from the current netCDF file. If
    these are not found in the current file, search the fixed file in
    the same directory. If not there, then search the tiled fixed files.

    :param fNcdf: an open netCDF file
    :type fNcdf: a netCDF file object

    :return: the ``ak`` and ``bk`` arrays

    .. NOTE:: This routine will look for both ``ak`` and ``bk``. There
    are cases when it is convenient to load the ``ak``, ``bk`` once when
    the files are first opened in ``MarsVars``, but the ``ak`` and
    ``bk`` arrays may not be necessary for in the calculation as is the
    case for ``MarsVars XXXXX.atmos_average_psd.nc --add msf``, which
    operates on a pressure interpolated (``_pstd.nc``) file.
    """
    # First try to read ak and bk in the current netCDF file:
    allvars = fNcdf.variables.keys()

    # Get netCDF file and path (for debugging)
    Ncdf_name = get_Ncdf_path(fNcdf)
    filepath,fname = extract_path_basename(Ncdf_name)
    fullpath_name = os.path.join(filepath, fname)

    if ("pk" in allvars or "ak" in allvars) and "bk" in allvars:
        # Check for ak first, then pk (pk for backwards compatibility)
        if "ak" in allvars:
            ak = np.array(fNcdf.variables["ak"])
        else:
            ak = np.array(fNcdf.variables["pk"])
        bk = np.array(fNcdf.variables["bk"])
        print("``ak`` and ``bk`` are in the file")

    else:
        try:
            name_fixed = find_fixedfile(fullpath_name)
            f_fixed = Dataset(name_fixed, "r",
                              format = "NETCDF4_CLASSIC")
            allvars = f_fixed.variables.keys()

            if "ak" in allvars:
                # Check for ``ak`` first, then ``pk``
                ak = np.array(f_fixed.variables["ak"])
            else:
                ak = np.array(f_fixed.variables["pk"])
            bk = np.array(f_fixed.variables["bk"])
            f_fixed.close()
            print("``pk`` and ``bk`` are in the fixed file")
        except:
            print(f"{Red}Fixed file does not exist in {filepath}. "
                  f"Make sure the fixed file you are referencing "
                  f"matches the FV3 filetype (e.g., ``fixed.tileX.nc`` "
                  f"for operations on tile X)")
            exit()
    return ak, bk

def read_variable_dict_amescap_profile(f_Ncdf=None):
    '''
    Inspect a Netcdf file and return the name of the variables and dimensions based on the content of ~/.amescap_profile.
    Calling this function allows to remove hard-coded calls in CAP.
    For example, to f.variables['ucomp'] is replaced by f.variables["ucomp"], with "ucomp" taking the values of'ucomp', 'U'
    Args:
        f_Ncdf: An opened Netcdf file object
    Returns:
        model: a dictionary with the dimensions and variables, e.g. "ucomp"='U' or "dim_lat"='latitudes'

    ***NOTE***
    The defaut names for variables are defined in () parenthesis in  ~/.amescap_profile :
    'X direction wind        [m/s]                   (ucomp)>'

    The defaut names for dimensions are defined in {} parenthesis in  ~/.amescap_profile :
    Ncdf Y latitude dimension    [integer]          {lat}>lats

    The dimensions (lon,lat,pfull,pstd) are loaded in the dictionary as "dim_lon", "dim_lat"
    '''

    if f_Ncdf is not None:
        var_list_Ncdf=list(f_Ncdf.variables.keys())
        dim_list_Ncdf=list(f_Ncdf.dimensions.keys())
    else:
        var_list_Ncdf=[]
        dim_list_Ncdf=[]

    all_lines=section_content_amescap_profile('Variable dictionary')
    lines=all_lines.split('\n')
    #Remove empty lines:
    while("" in lines):lines.remove("")

    #Initialize model
    class model(object):
        pass
    MOD=model()

    #Read through all lines in the Variable dictionary section of amesgcm_profile:
    for il in lines:
        var_list=[]
        #e.g. 'X direction wind [m/s]                          (ucomp)>U,u'
        left,right=il.split('>') #Split on either side of '>'

        #If using {var}, current entry is a dimension. If using (var), it is a variable
        if '{' in left:
            sep1='{';sep2='}';type_input='dimension'
        elif '(' in left:
            sep1='(';sep2=')';type_input='variable'

        # First get 'ucomp' from  'X direction wind [m/s]      (ucomp)
        _,tmp=FV3_var=left.split(sep1)
        FV3_var=tmp.replace(sep2,'').strip() #THIS IS THE FV3 NAME OF THE CURRENT VARIABLE
        #Then, get the list of variable on the righ-hand side, e.g.  'U,u'
        all_vars=right.split(',')
        for ii in all_vars:
            var_list.append(ii.strip())  #var_list IS A LIST OF POTENTIAL CORRESPONDING VARIABLES
        #Set the attribute to the dictionary
        #If the list is empty, e.g just [''], use the default FV3 variable presents in () or {}
        if len(var_list)==1 and var_list[0]=='':var_list[0]=FV3_var #var_list IS A LIST OF POTENTIAL CORRESPONDING VARIABLES
        found_list=[]

        #Place the input in the appropriate varialbe () or dimension {} dictionary
        #print('var_list>>>',var_list)
        if type_input=='variable':
            for ivar in var_list:
                if ivar in var_list_Ncdf:found_list.append(ivar)

            if len(found_list)==0:
                setattr(MOD,FV3_var,FV3_var)
            elif  len(found_list)==1:
                setattr(MOD,FV3_var,found_list[0])
            else:
                setattr(MOD,FV3_var,found_list[0])
                prYellow('''***Warning*** more than one possible variable '%s' found in file: %s'''%(FV3_var,found_list))
        if type_input=='dimension':
            for ivar in var_list:
                if ivar in dim_list_Ncdf:found_list.append(ivar)
            if len(found_list)==0:
                setattr(MOD,'dim_'+FV3_var,FV3_var)
            elif  len(found_list)==1:
                setattr(MOD,'dim_'+FV3_var,found_list[0])
            else:
                setattr(MOD,'dim_'+FV3_var,found_list[0])
                prYellow('''***Warning*** more than one possible dimension '%s' found in file: %s'''%(FV3_var,found_list))

    return MOD

def reset_FV3_names(MOD):
    '''
    This  function reset the model dictionary to the native FV3's variables, e.g.
    model.dim_lat = 'latitude' > model.dim_lat = 'lat'
    model.ucomp   = 'U'        > model.ucomp = 'ucomp'
    etc...

    Args:
        model: a class object generated with  read_variable_dict_amescap_profile()
    Returns:
        model: same object with updated names for the dimensions and variables.
    '''
    atts_list=dir(MOD) #Get all attributes
    vars_list=[k for k in atts_list if '__' not in k] #do not touch all the __init__ etc..

    for ivar in vars_list:
        name=ivar #get the native name, e.g ucomp
        if 'dim_' in ivar:
            name=ivar[4:] #if attribute is dim_lat, just get the 'lat' part
        setattr(MOD,ivar,name) #reset the original names
    return MOD