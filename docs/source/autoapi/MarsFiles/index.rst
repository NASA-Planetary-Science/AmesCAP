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



.. py:function:: make_FV3_files(fpath, typelistfv3, renameFV3=True, cwd=None)

   Make MGCM-like 'average', 'daily', and 'diurn' files.
   Args:
       fpath       : full path to the Legacy netcdf files
       typelistfv3 : MGCM-like file type: 'average', 'daily', or 'diurn'
       renameFV3   : rename the files from Legacy_Lsxxx_Lsyyy.nc to XXXXX.atmos_average.nc following MGCM output conventions
       cwd         : the output path
   Returns:
       atmos_average, atmos_daily, atmos_diurn


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


