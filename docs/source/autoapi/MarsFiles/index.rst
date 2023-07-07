:orphan:

:py:mod:`MarsFiles`
===================

.. py:module:: MarsFiles


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

   Make FV3-type atmos_average,atmos_daily,atmos_diurn
   Args:
       fpath     : full path to Legacy .nc files
       typelistfv3: e.g['average', 'daily', 'diurn']
       renameFV3 : rename files from Legacy_Lsxxx_Lsyyy.nc to XXXXX.atmos_average.nc folllowing FV3's convention
       cwd       : output path
   Returns:
       atmos_average,atmos_daily,atmos_diurn


.. py:function:: change_vname_longname_unit(vname, longname_txt, units_txt)

   Update variables names, longname and units. This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_dims(dims, todflag)

   Function to replace dimensions with fv3 names and remove tod.
   This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_at_index(tuple_dims, idx, new_name)

   Function to update dimensions.
   Args:
       tup: dimensions as tuples, e.g. ('pfull','nlat','nlon')
       idx      :  index:  axis of dimensions to update (e.g. ix=1  for 'nlat')
       new_name:    new dimension's name (e.g. 'latitude')


.. py:function:: ls2sol_1year(Ls_deg, offset=True, round10=True)

   Returns a sol number from the solar longitude.
   Args:
       Ls_deg: solar longitude in degree
       offset : if True, make year starts at Ls 0
       round10 : if True, round to the nearest 10 sols
   Returns:
       Ds :sol number
   ***NOTE***
   For the moment this is consistent with Ls 0->359.99, not for monotically increasing Ls


