:orphan:

:py:mod:`amescap.Script_utils`
==============================

.. py:module:: amescap.Script_utils


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

   This function return the Mars Year
   Args:
       Ls_cont: solar longitude, contineuous
   Returns:
       MY : int the Mars year


.. py:function:: find_tod_in_diurn(fNcdf)

   Return the variable for the local time axis in diurn files.
   Original implementation by Victoria H.
   Args:
       fNcdf: an (open) Netcdf file object
   Return:
       tod (string): 'time_of_day_16'or 'time_of_day_24'


.. py:function:: print_fileContent(fileNcdf)

   Print the content of a Netcdf file in a compact format. Variables are sorted by dimensions.
   Args:
       fileNcdf: full path to netcdf file
   Returns:
       None (print in the terminal)


.. py:function:: print_varContent(fileNcdf, list_varfull, print_stat=False)

   Print the content of a variable inside a Netcdf file
   This test is based on the existence of a least one  00XXX.fixed.nc in the current directory.
   Args:
       fileNcdf:      full path to netcdf file
       list_varfull:  list of variable names and optional slices, e.g ['lon' ,'ps[:,10,20]']
       print_stat:  if true, print min, mean and max instead of values
   Returns:
       None (print in the terminal)


.. py:function:: give_permission(filename)

   # NAS system only: set group permission to the file


.. py:function:: check_file_tape(fileNcdf, abort=False)

   Relevant for use on the NASA Advanced Supercomputing (NAS) environnment only
   Check if a file is present on the disk by running the NAS dmls -l data migration command.
   This avoid the program to stall if the files need to be migrated from the disk to the tape
   Args:
       fileNcdf: full path to netcdf file
       exit: boolean. If True, exit the program (avoid stalling the program if file is not on disk)
   Returns:
       None (print status and abort program)


.. py:function:: get_Ncdf_path(fNcdf)

   Return the full path of a Netcdf object.
   Note that 'Dataset' and multi-files dataset (i.e. 'MFDataset') have different
   attributes for the path, hence the need for this function.
   Args:
       fNcdf : Dataset or  MFDataset object
   Returns :
       Path: string (list) for Dataset (MFDataset)


.. py:function:: extract_path_basename(filename)

   Return the path and basename of a file. If only the filename is provided, assume it is the current directory
   Args:
       filename: e.g. 'XXXXX.fixed.nc', './history/XXXXX.fixed.nc' or '/home/user/history/XXXXX.fixed.nc'
   Returns:
       filepath : '/home/user/history/XXXXX.fixed.nc' in all the cases above
       basename:   XXXXX.fixed.nc in all the cases above

   ***NOTE***
   This routine does not check for file existence and only operates on the provided input string.


.. py:function:: FV3_file_type(fNcdf)

   Return the type of output files:
   Args:
       fNcdf: an (open) Netcdf file object
   Return:
      f_type (string): 'fixed', 'contineous', or 'diurn'
      interp_type (string): 'pfull','pstd','zstd','zagl'


.. py:function:: alt_FV3path(fullpaths, alt, test_exist=True)

   Internal function. given an interpolated daily, diurn or average file
   return the raw or fixed file. Accept string or list as input
   Args:
       fullpaths : e.g '/u/path/00010.atmos_average_pstd.nc' or LIST
           alt: alternative type 'raw' or 'fixed'
           test_exist=True test file existence on disk
   Returns :
       Alternative path to raw or fixed file, e.g.
                   '/u/path/00010.atmos_average.nc'
                   '/u/path/00010.fixed.nc'


.. py:function:: smart_reader(fNcdf, var_list, suppress_warning=False)

   Smarter alternative to using var=fNcdf.variables['var'][:] when handling PROCESSED files that also check
   matching XXXXX.atmos_average.nc (or daily...) and XXXXX.fixed.nc files

   Args:
       fNcdf: Netcdf file object (i.e. already opened with Dataset or MFDataset)
       var_list: variable or list of variables, e.g 'areo' or ['pk','bk','areo']
       suppress_warning: Suppress debug statement, useful if variable is not expected to be found in the file anyway
   Returns:
       out_list: variables content as singleton or values to unpack

   -------
   Example:

   from netCDF4 import Dataset

   fNcdf=Dataset('/u/akling/FV3/00668.atmos_average_pstd.nc','r')

   ucomp= fNcdf.variables['ucomp'][:]   # << this is the regular way
   vcomp= smart_reader(fNcdf,'vcomp')   # << this is exacly equivalent
   pk,bk,areo= smart_reader(fNcdf,['pk','bk','areo'])  # this will get 'areo' from 00668.atmos_average.nc is not available in the original _pstd.nc file
                                                       # if pk and bk are absent from 0668.atmos_average.nc, it will also check 00668.fixed.n
   *** NOTE ***
       -Only the variables' content is returned, not the attributes


.. py:function:: regrid_Ncfile(VAR_Ncdf, file_Nc_in, file_Nc_target)

   Regrid a Ncdf variable from one file's structure to match another file  [Alex Kling , May 2021]
   Args:
       VAR_Ncdf: A netCDF4 variable OBJECT, e.g. 'f_in.variables['temp']' from the source file
       file_Nc_in: The opened netcdf file object  for that input variable, e.g f_in=Dataset('fname','r')
       file_Nc_target: An opened netcdf file object  for the target grid t e.g f_out=Dataset('fname','r')
   Returns:
       VAR_OUT: the VALUES of VAR_Ncdf[:], interpolated on the grid for the target file.

   *** Note***
   While the KDTree interpolation can handle a 3D dataset (lon/lat/lev instead of just 2D lon/lat) , the grid points in the vertical are just a few 10's -100's meter in the PBL vs few 10'-100's km in the horizontal. This would results in excessive weighting in the vertical, which is why the vertical dimension is handled separately.


.. py:function:: progress(k, Nmax)

   Display a progress bar to monitor heavy calculations.
   Args:
       k: current iteration of the outer loop
       Nmax: max iteration of the outer loop
   Returns:
       Running... [#---------] 10.64 %


.. py:function:: section_content_amescap_profile(section_ID)

   Execude code section in /home/user/.amescap_profile
   Args:
       section_ID: string defining the section to load, e.g 'Pressure definitions for pstd'
   Returns
       return line in that section as python code


.. py:function:: filter_vars(fNcdf, include_list=None, giveExclude=False)

   Filter variable names in netcdf file for processing.
   Will return all dimensions (e.g. 'lon', 'lat'...), the 'areo' variable, and any variable included in include_list
   Args:
       fNcdf: an open netcdf4 object pointing to a diurn, daily or average file
       include_list: a list of variables to include, e.g. ['ucomp','vcomp']
       giveExclude: if True, instead return the variables that must be excluded from the file, i.e.
                    exclude_var = [all the variables] - [axis & dimensions] - [include_list]
   Return:
       var_list


.. py:function:: find_fixedfile(filename)

   Batterson, Updated by Alex Nov 29 2022
   Args:
       filename   =  name of FV3 data file in use, i.e.
                  'atmos_average.tile6.nc'
   Returns:
       name_fixed: fullpath to correspnding fixed file

           DDDDD.atmos_average.nc  -> DDDDD.fixed.nc
           atmos_average.tileX.nc  -> fixed.tileX.nc

           *variations of these work too*

           DDDDD.atmos_average_plevs.nc            -> DDDDD.fixed.nc
           DDDDD.atmos_average_plevs_custom.nc     -> DDDDD.fixed.nc
           atmos_average.tileX_plevs.nc            -> fixed.tileX.nc
           atmos_average.tileX_plevs_custom.nc     -> fixed.tileX.nc
           atmos_average_custom.tileX_plevs.nc     -> fixed.tileX.nc



.. py:function:: get_longname_units(fNcdf, varname)

   Return the 'long_name' and 'units'  attributes of a netcdf variable.
   If the attributes are not present, this function will return blank strings instead of raising an error
   Args:
       fNcdf: an opened netcdf file
       varname:  A variable to extract the attribute from (e.g. 'ucomp')
   Return:
       longname_txt : long_name attribute, e.g. 'zonal winds'
       units_txt    : units attribute, e.g. [m/s]

   *** NOTE***
   Some functions in MarsVars edit the units, e.g. turn [kg]  to [kg/m], therefore the empty string is made 4
   characters in length ('    ' instead of '') to allow for editing by editing units_txt[:-2] for example


.. py:function:: wbr_cmap()

   Returns a color map that goes from white>blue>green>yellow>red or 'wbr'


.. py:function:: rjw_cmap()

   Returns a color map that goes from red<jade<wisteria or 'rjw'


.. py:function:: dkass_dust_cmap()

   Returns a color map that goes from yellow>orange>red>purple
   Provided by Courtney B.


.. py:function:: dkass_temp_cmap()

   Returns a color map that goes from black>purple>blue>green>yellow>orange>red
   Provided by Courtney B.


.. py:function:: pretty_print_to_fv_eta(var, varname, nperline=6)

   Print the ak or bk coefficients to copy paste in fv_eta.f90
   Args:
       data: ak or bk data
       varname: the variable name, 'a' or 'b'
       nperline:the number of elements per line
   Returns:
        The print statement ready to copy-paste in fv_eta.f90



.. py:function:: replace_dims(Ncvar_dim, vert_dim_name=None)

   Update the name for the variables dimension to match FV3's
   Args:
       Ncvar_dim: Netcdf variable dimensions, e.g f_Ncdf.variables['temp'].dimensions
       vert_dim_name(optional): 'pstd', 'zstd', or 'zagl' if the vertical dimensions is ambigeous
   Return:
       dims_out: updated dimensions matching FV3's naming convention



.. py:function:: ak_bk_loader(fNcdf)

   Return the ak and bk. First look in the current netcdf file.
   If not found, this routine will check the XXXXX.fixed.nc in the same directory and then in the XXXXX.fixed.tileX.nc files if present
   Args:
       fNcdf: an opened netcdf file
   Returns:
       ak,bk : the ak, bk values
   ***NOTE***

   This routine will look for both 'ak' and 'pk'.

   There are cases when it is convenient to load the  pk, bk once at the begining of the files in MarsVars.py,
   However the pk, bk may not be used at all in the calculation. This is the case with MarsVars.py XXXXX.atmos_average_psd.nc --add msf (which operates on the _pstd.nc file)




