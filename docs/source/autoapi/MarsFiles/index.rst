:py:mod:`MarsFiles`
===================

.. py:module:: MarsFiles

.. autoapi-nested-parse::

   The MarsFiles executable has functions for manipulating entire files. The capabilities include time-shifting, binning, and regridding data, as well as band pass filtering, tide analysis, zonal averaging, and extracting variables from files. 

   The executable requires:
       * ``[input_file]``      The file for manipulation

   and optionally accepts:
       * ``[-fv3, --fv3]``                 Produce MGCM ``fixed``,         ``diurn``, ``average`` and ``daily`` files from Legacy output
       * ``[-c, --combine]``               Combine sequential files of         the same type into one file
       * ``[-t, --tshift]``                Apply a time-shift to         ``diurn`` files
       * ``[-ba, --bin_average]``          Bin MGCM ``daily`` files like         ``average`` files
       * ``[-bd, --bin_diurn]``            Bin MGCM ``daily`` files like         ``diurn`` files
       * ``[-hp, --high_pass_filter]``     Temporal filtering: high-pass
       * ``[-lp, --low_pass_filter]``      Temporal filtering: low-pass
       * ``[-bp, --band_pass_filter]``     Temporal filtering: band-pass
       * ``[-no_trend, --no_trend]``       Filter and compute amplitudes         only (use with filtering)
       * ``[-hpk, --high_pass_zonal]``     Spatial filtering: high-pass
       * ``[-lpk, --low_pass_zonal]``      Spatial filtering: low-pass
       * ``[-bpk, --band_pass_zonal]``     Spatial filtering: band-pass
       * ``[-tidal, --tidal]``             Extracts diurnal tide and its         harmonics
       * ``[-reconstruct, --reconstruct]`` Reconstructs the first N         harmonics
       * ``[-norm, --normalize]``          Provides ``-tidal`` result in         percent amplitude
       * ``[-rs, --regrid_source]``        Regrid a target file to match         a source file
       * ``[-za, --zonal_avg]``            Zonally average all variables         in a file
       * ``[-include, --include]``         Only include specific         variables from the target file
       * ``[-e, --ext]``                   Create a new file with a         unique extension instead of overwriting current file
       
   Third-party Requirements:
       * ``numpy``
       * ``netCDF4``
       * ``sys``
       * ``argparse``
       * ``os``
       * ``subprocess``
       * ``warnings``



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

   Make MGCM-like ``average``, ``daily``, and ``diurn`` files.
   Used if call to [``-fv3 --fv3``] is made AND Legacy files are in     netCDFformat (not fort.11).

   :param fpath: Full path to the Legacy netcdf files
   :type fpath: str
   :param typelistfv3: MGCM-like file type: ``average``, ``daily``,         or ``diurn``
   :type typelistfv3: list
   :param renameFV3: Rename the files from Legacy_LsXXX_LsYYY.nc to         ``XXXXX.atmos_average.nc`` following MGCM output conventions
   :type renameFV3: bool

   :return: The MGCM-like files: ``XXXXX.atmos_average.nc``,         ``XXXXX.atmos_daily.nc``, ``XXXXX.atmos_diurn.nc``.


.. py:function:: change_vname_longname_unit(vname, longname_txt, units_txt)

   Update variable ``name``, ``longname``, and ``units``.
   This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_dims(dims, todflag)

   Function replaces dimensions with MGCM-like names and remove 
   time_of_day. This was designed specifically for LegacyCGM.nc files.


.. py:function:: replace_at_index(tuple_dims, idx, new_name)

   Update dimensions.

   :param tuple_dims: the dimensions as tuples e.g. (``pfull``,         ``nlat``, ``nlon``)
   :type tuple_dims: tuple
   :param idx: index indicating axis with the dimensions to update         (e.g. ``idx = 1``  for ``nlat``)
   :type idx: int
   :param new_name: new dimension name (e.g. ``latitude``)
   :type new_name: str

   :return: updated dimensions


.. py:function:: ls2sol_1year(Ls_deg, offset=True, round10=True)

   Returns a sol number from the solar longitude.

   :param Ls_deg: solar longitude in degrees
   :type Ls_deg: float
   :param offset: if True, force year to start at Ls 0
   :type offset: bool
   :param round10: if True, round to the nearest 10 sols
   :type round10: bool

   :returns: ``Ds`` sol number

   .. NOTE:: For the moment, this is consistent with 0 <= Ls <=         359.99, but not for monotically increasing Ls.


