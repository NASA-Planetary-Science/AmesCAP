:py:mod:`amescap.Script_utils`
==============================

.. py:module:: amescap.Script_utils

.. autoapi-nested-parse::

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
       * ``sys``
       * ``math``
       * ``matplotlib.colors``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.Script_utils.FV3_file_type
   amescap.Script_utils.MY_func
   amescap.Script_utils.ak_bk_loader
   amescap.Script_utils.alt_FV3path
   amescap.Script_utils.check_bounds
   amescap.Script_utils.check_file_tape
   amescap.Script_utils.dkass_dust_cmap
   amescap.Script_utils.dkass_temp_cmap
   amescap.Script_utils.except_message
   amescap.Script_utils.extract_path_basename
   amescap.Script_utils.filter_vars
   amescap.Script_utils.find_fixedfile
   amescap.Script_utils.find_tod_in_diurn
   amescap.Script_utils.get_Ncdf_path
   amescap.Script_utils.get_longname_unit
   amescap.Script_utils.give_permission
   amescap.Script_utils.hot_cold_cmap
   amescap.Script_utils.pretty_print_to_fv_eta
   amescap.Script_utils.print_fileContent
   amescap.Script_utils.print_varContent
   amescap.Script_utils.progress
   amescap.Script_utils.read_variable_dict_amescap_profile
   amescap.Script_utils.regrid_Ncfile
   amescap.Script_utils.replace_dims
   amescap.Script_utils.reset_FV3_names
   amescap.Script_utils.rjw_cmap
   amescap.Script_utils.section_content_amescap_profile
   amescap.Script_utils.smart_reader
   amescap.Script_utils.wbr_cmap



Attributes
~~~~~~~~~~

.. autoapisummary::

   amescap.Script_utils.Blue
   amescap.Script_utils.Cyan
   amescap.Script_utils.Green
   amescap.Script_utils.Nclr
   amescap.Script_utils.Purple
   amescap.Script_utils.Red
   amescap.Script_utils.Yellow


.. py:function:: FV3_file_type(fNcdf)

   Return the type of the netCDF file.

   Returns netCDF file type (i.e., ``fixed``, ``diurn``, ``average``,
   ``daily``) and the format of the Ls array ``areo`` (i.e., ``fixed``,
   ``continuous``, or ``diurn``).

   :param fNcdf: an open Netcdf file
   :type  fNcdf: Netcdf file object
   :return: The Ls array type (string, ``fixed``, ``continuous``, or
       ``diurn``) and the netCDF file type (string ``fixed``,
       ``diurn``, ``average``, or ``daily``)
   :rtype:  tuple
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute


.. py:function:: MY_func(Ls_cont)

   Returns the Mars Year.

   :param Ls_cont: solar longitude (continuous)
   :type  Ls_cont: array

   :return: the Mars year
   :rtype:  int

   :raises ValueError: if Ls_cont is not in the range [0, 360)


.. py:function:: ak_bk_loader(fNcdf)

   Loads the ak and bk variables from a netCDF file.

   This function will first check the current netCDF file for the
   ``ak`` and ``bk`` variables. If they are not found, it will
   search the fixed file in the same directory. If they are still
   not found, it will search the tiled fixed files. The function
   will return the ``ak`` and ``bk`` arrays.

   :param fNcdf: an open netCDF file
   :type  fNcdf: a netCDF file object
   :return: the ``ak`` and ``bk`` arrays
   :rtype:  tuple
   :raises ValueError: if the ``ak`` and ``bk`` variables are not
       found in the netCDF file, the fixed file, or the tiled fixed files

   .. note::
       This routine will look for both ``ak`` and ``bk``. There
       are cases when it is convenient to load the ``ak``, ``bk`` once
       when the files are first opened in ``MarsVars``, but the ``ak``
       and ``bk`` arrays may not be necessary for in the calculation
       as is the case for ``MarsVars XXXXX.atmos_average_psd.nc
       --add msf``, which operates on a pressure interpolated
       (``_pstd.nc``) file.


.. py:function:: alt_FV3path(fullpaths, alt, test_exist=True)

   Returns the original or fixed file for a given path.

   :param fullpaths: full path to a file or a list of full paths to
       more than one file
   :type  fullpaths: str
   :param alt: type of file to return (i.e., original or fixed)
   :type  alt: str
   :param test_exist: Whether file exists on the disk, defaults to True
   :type  test_exist: bool, optional
   :return: path to original or fixed file (e.g.,
       /u/path/00010.atmos_average.nc or /u/path/00010.fixed.nc)
   :rtype:  str
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path


.. py:function:: check_bounds(values, min_val, max_val, dx)

   Checks the bounds of a variable in a netCDF file.

   This function checks if the values in a netCDF file are within
   the specified bounds. If any value is out of bounds, it will
   print an error message and exit the program.
   The function can handle both single values and arrays.

   Parameters:
   :param values: Single value or array of values to check
   :type   values: array-like
   :param min_val: Minimum allowed value
   :type   min_val: float
   :param max_val: Maximum allowed value
   :type   max_val: float
   :return values: The validated value(s)
   :rtype:  array or float
   :raises ValueError: if values is out of bounds or if values is not
       a number, array, or list


.. py:function:: check_file_tape(fileNcdf)

   Checks whether a file exists on the disk.

   If on a NAS system, also checks if the file needs to be migrated
   from tape.

   :param fileNcdf: full path to a netcdf file or a file object with a name attribute
   :type  fileNcdf: str or file object

   :return: None

   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises subprocess.CalledProcessError: if the dmls command fails


.. py:function:: dkass_dust_cmap()

   Color map (From Courtney Batterson).

   Returns a color map useful for dust cross-sections that highlight
   dust mixing ratios > 4 ppm. The color map goes from
   white -> yellow -> orange -> red -> purple.
   [Courtney Batterson, Nov 2022]

   :return: color map
   :rtype:  array


.. py:function:: dkass_temp_cmap()

   Color map (From Courtney Batterson).

   Returns a color map useful for highlighting the 200 K temperature
   level. The color map goes from
   black -> purple -> blue -> green -> yellow -> orange -> red.
   [Courtney Batterson, Nov 2022]

   :return: color map
   :rtype:  array


.. py:function:: except_message(debug, exception, varname, ifile, pre='', ext='')

   Prints an error message if a variable is not found.

   It also contains a special error in the case of a pre-existing
   variable.

   :param debug: Flag for debug mode
   :type  debug: logical
   :param exception: Exception from try statement
   :type  exception: class object
   :param varname: Name of variable causing exception
   :type  varname: string
   :param ifile: Name of input file
   :type  ifile: string
   :param pre: Prefix to new variable
   :type  pre: string
   :param ext: Extension to new variable
   :type  ext: string
   :return: None
   :rtype:  None
   :raises ValueError: if debug is True, exception is not a class
       object or string, varname is not a string, ifile is not a
       string, pre is not a string, or ext is not a string


.. py:function:: extract_path_basename(filename)

   Returns the path and basename of a file.

   If only the filename is provided, assume it is in the current
   directory.

   :param filename: name of the netCDF file (may include full path)
   :type  filename: str
   :return: full file path & name of file
   :rtype:  tuple
   :raises ValueError: if the filename is not a string
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the filename is not a string
   :raises OSError: if the filename is not a valid path

   .. note::
       This routine does not confirm that the file exists. It operates
       on the provided input string.


.. py:function:: filter_vars(fNcdf, include_list=None, giveExclude=False)

   Filters the variable names in a netCDF file for processing.

   Returns all dimensions (``lon``, ``lat``, etc.), the ``areo``
   variable, and any other variable listed in ``include_list``.

   :param fNcdf: an open netCDF object for a diurn, daily, or average
       file
   :type  fNcdf: netCDF file object
   :param include_list:list of variables to include (e.g., [``ucomp``,
       ``vcomp``], defaults to None
   :type  include_list: list or None, optional
   :param giveExclude: if True, returns variables to be excluded from
       the file, defaults to False
   :type  giveExclude: bool, optional
   :return: list of variable names to include in the processed file
   :rtype:  list
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path
   :raises KeyError: if the variable is not found in the file


.. py:function:: find_fixedfile(filename)

   Finds the relevant fixed file for an average, daily, or diurn file.

   [Courtney Batterson, updated by Alex Nov 29 2022]

   :param filename: an average, daily, or diurn netCDF file
   :type  filename: str
   :return: full path to the correspnding fixed file
   :rtype:  str
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path

   Compatible file types::

       DDDDD.atmos_average.nc              -> DDDDD.fixed.nc
       atmos_average.tileX.nc              -> fixed.tileX.nc
       DDDDD.atmos_average_plevs.nc        -> DDDDD.fixed.nc
       DDDDD.atmos_average_plevs_custom.nc -> DDDDD.fixed.nc
       atmos_average.tileX_plevs.nc        -> fixed.tileX.nc
       atmos_average.tileX_plevs_custom.nc -> fixed.tileX.nc
       atmos_average_custom.tileX_plevs.nc -> fixed.tileX.nc


.. py:function:: find_tod_in_diurn(fNcdf)

   Returns the variable for the local time axis in diurn files.

   (e.g., time_of_day_24). Original implementation by Victoria H.

   :param fNcdf: a netCDF file
   :type  fNcdf: netCDF file object

   :return: the name of the time of day dimension
   :rtype:  str
   :raises ValueError: if the time of day variable is not found


.. py:function:: get_Ncdf_path(fNcdf)

   Returns the full path for a netCDF file object.

   ``Dataset`` and multi-file dataset (``MFDataset``) have different
   attributes for the path, hence the need for this function.

   :param fNcdf: Dataset or MFDataset object
   :type  fNcdf: netCDF file object

   :return: string list for the Dataset (MFDataset)
   :rtype:  str(list)

   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist


.. py:function:: get_longname_unit(fNcdf, varname)

   Returns the longname and unit of a variable.

   If the attributes are unavailable, returns blank strings to avoid
   an error.

   :param fNcdf: an open netCDF file
   :type  fNcdf: netCDF file object
   :param varname: variable to extract attribute from
   :type  varname: str
   :return: longname and unit attributes
   :rtype:  str
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path
   :raises KeyError: if the variable is not found in the file

   .. note::
       Some functions in MarsVars edit the units
       (e.g., [kg] -> [kg/m]), therefore the empty string is 4
       characters in length ("    " instead of "") to allow for
       editing by ``editing units_txt[:-2]``, for example.


.. py:function:: give_permission(filename)

   Sets group file permissions for the NAS system

   :param filename: full path to the netCDF file
   :type  filename: str

   :return: None

   :raises subprocess.CalledProcessError: if the setfacl command fails
   :raises FileNotFoundError: if the file is not found


.. py:function:: hot_cold_cmap()

   Returns a color map (From Alex Kling, based on bipolar cmap).

   Color map goes from dark blue -> light blue -> white -> yellow -> red.
   Based on Matlab's bipolar colormap.
   [Alex Kling, Nov 2022]

   :return: color map
   :rtype:  array


.. py:function:: pretty_print_to_fv_eta(var, varname, nperline=6)

   Print the ``ak`` or ``bk`` coefficients for copying to ``fv_eta.f90``.

   The ``ak`` and ``bk`` coefficients are used to calculate the
   vertical coordinate transformation.

   :param var: ak or bk data
   :type  var: array
   :param varname: the variable name ("a" or "b")
   :type  varname: str
   :param nperline: the number of elements per line, defaults to 6
   :type  nperline: int, optional
   :return: a print statement for copying into ``fv_eta.f90``
   :rtype:  None
   :raises ValueError: if varname is not "a" or "b"
   :raises ValueError: if nperline is not a positive integer
   :raises ValueError: if var is not a 1D array of length NLAY+1


.. py:function:: print_fileContent(fileNcdf)

   Prints the contents of a netCDF file to the screen.

   Variables sorted by dimension.

   :param fileNcdf: full path to the netCDF file
   :type  fileNcdf: str

   :return: None


.. py:function:: print_varContent(fileNcdf, list_varfull, print_stat=False)

   Print variable contents from a variable in a netCDF file.

   Requires a XXXXX.fixed.nc file in the current directory.

   :param fileNcdf: full path to a netcdf file
   :type  fileNcdf: str
   :param list_varfull: list of variable names and optional slices
       (e.g., ``["lon", "ps[:, 10, 20]"]``)
   :type  list_varfull: list
   :param print_stat: If True, print min, mean, and max. If False,
       print values. Defaults to False
   :type  print_stat: bool, optional

   :return: None

   :raises ValueError: if the variable is not found in the file
   :raises FileNotFoundError: if the file is not found
   :raises Exception: if the variable is not found in the file
   :raises Exception: if the file is not found


.. py:function:: progress(k, Nmax)

   Displays a progress bar to monitor heavy calculations.

   :param k: current iteration of the outer loop
   :type  k: int
   :param Nmax: max iteration of the outer loop
   :type  Nmax: int
   :return: None
   :raises ValueError: if k or Nmax are not integers, k > Nmax, or k < 0


.. py:function:: read_variable_dict_amescap_profile(f_Ncdf=None)

   Reads a variable dictionary from the ``amescap_profile`` file.

   This function will read the variable dictionary from the
   ``amescap_profile`` file and return a dictionary with the
   variable names and dimensions. The function will also check
   the opened netCDF file for the variable names and dimensions.

   :param f_Ncdf: An opened Netcdf file object
   :type  f_Ncdf: File object
   :return: MOD, a class object with the variable names and dimensions
       (e.g., ``MOD.ucomp`` = 'U' or ``MOD.dim_lat`` = 'latitudes')
   :rtype  MOD: class object
   :raises ValueError: if the ``amescap_profile`` file is not found
       or if the variable dictionary is not found
   :raises ValueError: if the variable or dimension name is not found
       in the netCDF file
   :raises ValueError: if the variable or dimension name is not found
       in the``amescap_profile``


.. py:function:: regrid_Ncfile(VAR_Ncdf, file_Nc_in, file_Nc_target)

   Regrid a netCDF variable from one file structure to another.

   Requires a file with the desired file structure to mimic.
   [Alex Kling, May 2021]

   :param VAR_Ncdf: a netCDF variable object to regrid
       (e.g., ``f_in.variable["temp"]``)
   :type  VAR_Ncdf: netCDF file variable
   :param file_Nc_in: an open netCDF file to source for the variable
       (e.g., ``f_in = Dataset("filename", "r")``)
   :type  file_Nc_in: netCDF file object
   :param file_Nc_target: an open netCDF file with the desired file
       structure (e.g., ``f_out = Dataset("filename", "r")``)
   :type  file_Nc_target: netCDF file object
   :return: the values of the variable interpolated to the target file
       grid.
   :rtype:  array
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path

   .. note::
       While the KDTree interpolation can handle a 3D dataset
       (lon/lat/lev instead of just 2D lon/lat), the grid points in
       the vertical are just a few (10--100s) meters in the PBL vs a
       few (10-100s) kilometers in the horizontal. This results in
       excessive weighting in the vertical, which is why the vertical
       dimension is handled separately.


.. py:function:: replace_dims(Ncvar_dim, vert_dim_name=None)

   Replaces the dimensions of a variable in a netCDF file.

   Updates the name of the variable dimension to match the format of
   the new NASA Ames Mars GCM output files.

   :param Ncvar_dim: netCDF variable dimensions
       (e.g., ``f_Ncdf.variables["temp"].dimensions``)
   :type  Ncvar_dim: str
   :param vert_dim_name: the vertical dimension if it is ambiguous
       (``pstd``, ``zstd``, or ``zagl``). Defaults to None
   :type  vert_dim_name: str, optional
   :return: updated dimensions
   :rtype:  str
   :raises ValueError: if Ncvar_dim is not a list
   :raises ValueError: if vert_dim_name is not a string
   :raises ValueError: if vert_dim_name is not in the list of
       recognized vertical dimensions


.. py:function:: reset_FV3_names(MOD)

   Resets the FV3 variable names in a netCDF file.

   This function resets the model dictionary to the native FV3
   variable names, e.g.::

       model.dim_lat = 'latitude' > model.dim_lat = 'lat'
       model.ucomp   = 'U'        > model.ucomp = 'ucomp'

   :param MOD: Generated with read_variable_dict_amescap_profile()
   :type  MOD: class object
   :return: same object with updated names for the dimensions and
   variables
   :rtype:  class object
   :raises ValueError: if the MOD object is not a class object or does
       not contain the expected attributes


.. py:function:: rjw_cmap()

   Returns a color map (From R. John Wilson).

   Color map goes from red -> jade -> wisteria.
   [R. John Wilson, Nov 2022]

   :return: color map
   :rtype:  array


.. py:function:: section_content_amescap_profile(section_ID)

   Executes first code section in ``~/.amescap_profile``.

   Reads in user-defined plot & interpolation settings.
   [Alex Kling, Nov 2022]

   :param section_ID: the section to load (e.g., Pressure definitions
       for pstd)
   :type  section_ID: str
   :return: the relevant line with Python syntax
   :rtype:  str
   :raises FileNotFoundError: if the file is not found


.. py:function:: smart_reader(fNcdf, var_list, suppress_warning=False)

   Reads a variable from a netCDF file.

   If the variable is not found in the file, it checks for the variable
   in the original file (e.g., atmos_average.nc) or the fixed file
   (e.g., 00010.fixed.nc).

   Alternative to ``var = fNcdf.variables["var"][:]`` for handling
   *processed* files.

   :param fNcdf: an open netCDF file
   :type  fNcdf: netCDF file object
   :param var_list: a variable or list of variables (e.g., ``areo`` or
       [``pk``, ``bk``, ``areo``])
   :type  var_list: _type_
   :param suppress_warning: suppress debug statement. Useful if a
       variable is not expected to be in the file anyway. Defaults to
       False
   :type  suppress_warning: bool, optional
   :return: variable content (single or values to unpack)
   :rtype:  list
   :raises ValueError: if the file is not a netCDF file
   :raises FileNotFoundError: if the file does not exist
   :raises TypeError: if the file is not a Dataset or MFDataset
   :raises AttributeError: if the file does not have the _files or
       filepath attribute
   :raises OSError: if the file is not a valid path

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

   .. note::
       Only the variable content is returned, not attributes


.. py:function:: wbr_cmap()

   Returns a color map (From R. John Wilson).

   Color map goes from white -> blue -> green -> yellow -> red
   [R. John Wilson, Nov 2022]

   :return: color map
   :rtype:  array


.. py:data:: Blue
   :value: '\x1b[94m'

   

.. py:data:: Cyan
   :value: '\x1b[96m'

   

.. py:data:: Green
   :value: '\x1b[92m'

   

.. py:data:: Nclr
   :value: '\x1b[00m'

   

.. py:data:: Purple
   :value: '\x1b[95m'

   

.. py:data:: Red
   :value: '\x1b[91m'

   

.. py:data:: Yellow
   :value: '\x1b[93m'

   

