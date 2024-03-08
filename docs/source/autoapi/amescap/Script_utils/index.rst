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



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   amescap.Script_utils.MY_func
   amescap.Script_utils.find_tod_in_diurn
   amescap.Script_utils.print_fileContent
   amescap.Script_utils.print_varContent
   amescap.Script_utils.give_permission
   amescap.Script_utils.check_file_tape
   amescap.Script_utils.get_Ncdf_path
   amescap.Script_utils.extract_path_basename
   amescap.Script_utils.FV3_file_type
   amescap.Script_utils.alt_FV3path
   amescap.Script_utils.smart_reader
   amescap.Script_utils.regrid_Ncfile
   amescap.Script_utils.progress
   amescap.Script_utils.section_content_amescap_profile
   amescap.Script_utils.filter_vars
   amescap.Script_utils.find_fixedfile
   amescap.Script_utils.get_longname_units
   amescap.Script_utils.wbr_cmap
   amescap.Script_utils.rjw_cmap
   amescap.Script_utils.dkass_dust_cmap
   amescap.Script_utils.dkass_temp_cmap
   amescap.Script_utils.pretty_print_to_fv_eta
   amescap.Script_utils.replace_dims
   amescap.Script_utils.ak_bk_loader



.. py:function:: MY_func(Ls_cont)

   Returns the Mars Year

   :param Ls_cont: solar longitude (continuous)
   :type Ls_cont: array

   :return: the Mars year
   :rtype: int


.. py:function:: find_tod_in_diurn(fNcdf)

   Returns the variable for the local time axis in diurn files
   (e.g., time_of_day_24).
   Original implementation by Victoria H.

   :param fNcdf: a netCDF file
   :type fNcdf: netCDF file object

   :return: the name of the time of day dimension
   :rtype: str


.. py:function:: print_fileContent(fileNcdf)

   Prints the contents of a netCDF file to the screen. Variables sorted
   by dimension.

   :param fileNcdf: full path to the netCDF file
   :type fileNcdf: str

   :return: None


.. py:function:: print_varContent(fileNcdf, list_varfull, print_stat=False)

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


.. py:function:: give_permission(filename)

   Sets group file permissions for the NAS system 


.. py:function:: check_file_tape(fileNcdf, abort=False)

   For use in the NAS environnment only.
   Checks whether a file is exists on the disk by running the command
   ``dmls -l`` on NAS. This prevents the program from stalling if the
   files need to be migrated from the disk to the tape.

   :param fileNcdf: full path to a netcdf file
   :type fileNcdf: _type_
   :param abort: If True, exit the program. Defaults to False
   :type abort: bool, optional

   :return: None


.. py:function:: get_Ncdf_path(fNcdf)

   Returns the full path for a netCDF file object.

   .. NOTE:: ``Dataset`` and multi-file dataset (``MFDataset``) have
   different attributes for the path, hence the need for this function.

   :param fNcdf: Dataset or MFDataset object
   :type fNcdf: netCDF file object
   :return: string list for the Dataset (MFDataset)
   :rtype: str(list)


.. py:function:: extract_path_basename(filename)

   Returns the path and basename of a file. If only the filename is
   provided, assume it is in the current directory.

   :param filename: name of the netCDF file (may include full path)
   :type filename: str

   :return: full file path & name of file

   .. NOTE:: This routine does not confirm that the file exists.
       It operates on the provided input string.


.. py:function:: FV3_file_type(fNcdf)

   Return the type of the netCDF file (i.e., ``fixed``, ``diurn``,
   ``average``, ``daily``) and the format of the Ls array ``areo``
   (i.e., ``fixed``, ``continuous``, or ``diurn``).

   :param fNcdf: an open Netcdf file
   :type fNcdf: Netcdf file object

   :return: The Ls array type (string, ``fixed``, ``continuous``, or
       ``diurn``) and the netCDF file type (string ``fixed``,
       ``diurn``, ``average``, or ``daily``)


.. py:function:: alt_FV3path(fullpaths, alt, test_exist=True)

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


.. py:function:: smart_reader(fNcdf, var_list, suppress_warning=False)

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


.. py:function:: regrid_Ncfile(VAR_Ncdf, file_Nc_in, file_Nc_target)

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


.. py:function:: progress(k, Nmax)

   Displays a progress bar to monitor heavy calculations.

   :param k: current iteration of the outer loop
   :type k: int
   :param Nmax: max iteration of the outer loop
   :type Nmax: int


.. py:function:: section_content_amescap_profile(section_ID)

   Executes first code section in ``~/.amescap_profile`` to read in
   user-defined plot & interpolation settings.

   :param section_ID: the section to load (e.g., Pressure definitions
       for pstd)
   :type section_ID: str

   :return: the relevant line with Python syntax


.. py:function:: filter_vars(fNcdf, include_list=None, giveExclude=False)

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


.. py:function:: find_fixedfile(filename)

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



.. py:function:: get_longname_units(fNcdf, varname)

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


.. py:function:: wbr_cmap()

   Returns a color map that goes from
   white -> blue -> green -> yellow -> red


.. py:function:: rjw_cmap()

   Returns John Wilson's preferred color map
   (red -> jade -> wisteria)


.. py:function:: dkass_dust_cmap()

   Returns a color map useful for dust cross-sections.
   (yellow -> orange -> red -> purple)
   Provided by Courtney Batterson.


.. py:function:: dkass_temp_cmap()

   Returns a color map that highlights the 200K temperatures.
   (black -> purple -> blue -> green -> yellow -> orange -> red)
   Provided by Courtney Batterson.


.. py:function:: pretty_print_to_fv_eta(var, varname, nperline=6)

   Print the ``ak`` or ``bk`` coefficients for copying to
   ``fv_eta.f90``.

   :param var: ak or bk data
   :type var: array
   :param varname: the variable name ("a" or "b")
   :type varname: str
   :param nperline: the number of elements per line, defaults to 6
   :type nperline: int, optional

   :return: a print statement for copying into ``fv_eta.f90``


.. py:function:: replace_dims(Ncvar_dim, vert_dim_name=None)

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


.. py:function:: ak_bk_loader(fNcdf)

   Return ``ak`` and ``bk`` arrays from the current netCDF file. If
   these are not found in the current file, search the fixed file in
   the same directory. If not there, then search the tiled fixed files.

   :param fNcdf: an open netCDF file
   :type fNcdf: a netCDF file object

   :return: the ``ak`` and ``bk`` arrays

   .. NOTE:: This routine will look for both ``ak`` and ``bk``. There
   are cases when it is convenient to load the ``ak``, ``bk`` once when
   the files are first opened in ``MarsVars.py``, but the ``ak`` and
   ``bk`` arrays may not be necessary for in the calculation as is the
   case for ``MarsVars.py XXXXX.atmos_average_psd.nc --add msf``, which
   operates on a pressure interpolated (``_pstd.nc``) file.


