:py:mod:`MarsFiles`
===================

.. py:module:: MarsFiles

.. autoapi-nested-parse::

   The MarsFiles executable is for ...

   The executable requires x arguments:
       * [-x --x]      define

   Third-party Requirements:
       * numpy
       * argparse
       * requests

   List of Functions:
       * x



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsFiles.make_FV3_files
   MarsFiles.change_vname_longname_unit
   MarsFiles.replace_dims
   MarsFiles.replace_at_index
   MarsFiles.ls2sol_1year



.. py:function:: make_FV3_files(fpath, typelistfv3, renameFV3=True)

   Make MGCM-like 'average', 'daily', and 'diurn' files.

   Used if call to -fv3 --fv3 is made AND Legacy files are in netCDF
   format (not fort.11).

   Parameters
   ----------
   fpath : str
       Full path to the Legacy netcdf files
   typelistfv3 : list
       MGCM-like file type: 'average', 'daily', or 'diurn'
   renameFV3 : bool
       Rename the files from Legacy_LsXXX_LsYYY.nc to             XXXXX.atmos_average.nc following MGCM output conventions

   Returns
   -------
   The MGCM-like files: XXXXX.atmos_average.nc, XXXXX.atmos_daily.nc,         XXXXX.atmos_diurn.nc


.. py:function:: change_vname_longname_unit(vname, longname_txt, units_txt)

   Update variable name, longname, and units.
   This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_dims(dims, todflag)

   Function to replace dimensions with MGCM-like names and remove 'time_of_day'.
   This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_at_index(tuple_dims, idx, new_name)

   Function to update dimensions.
   Args:
       tup      : the dimensions as tuples e.g. ('pfull', 'nlat', 'nlon')
       idx      : index indicating axis with the dimensions to update (e.g. idx = 1  for 'nlat')
       new_name : new dimension name (e.g. 'latitude')


.. py:function:: ls2sol_1year(Ls_deg, offset=True, round10=True)

   Returns a sol number from the solar longitude.
   Args:
       Ls_deg  : solar longitude in degrees
       offset  : if True, force year to start at Ls 0
       round10 : if True, round to the nearest 10 sols
   Returns:
       Ds: sol number
   ***NOTE***
   For the moment, this is consistent with Ls 0 -> 359.99, but not for monotically increasing Ls.


