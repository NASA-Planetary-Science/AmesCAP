:py:mod:`MarsFiles`
===================

.. py:module:: MarsFiles

.. autoapi-nested-parse::

   The MarsFiles executable has functions for manipulating entire files.
   The capabilities include time-shifting, binning, and regridding data,
   as well as band pass filtering, tide analysis, zonal averaging, and
   extracting variables from files.

   The executable requires:
       * ``[input_file]``                  the file for manipulation

   and optionally accepts:
       * ``[-fv3, --fv3]``                 produce MGCM ``fixed``,
           ``diurn``, ``average`` and ``daily`` files from Legacy output
       * ``[-c, --combine]``               Combine sequential files of
           the same type into one file
       * ``[-t, --tshift]``                apply a time-shift to
           ``diurn`` files
       * ``[-ba, --bin_average]``          bin MGCM ``daily`` files like
           ``average`` files
       * ``[-bd, --bin_diurn]``            bin MGCM ``daily`` files like
           ``diurn`` files
       * ``[-hp, --high_pass_filter]``     temporal filtering: high-pass
       * ``[-lp, --low_pass_filter]``      temporal filtering: low-pass
       * ``[-bp, --band_pass_filter]``     temporal filtering: band-pass
       * ``[-no_trend, --no_trend]``       filter and compute amplitudes
           only (use with filtering)
       * ``[-hpk, --high_pass_zonal]``     spatial filtering: high-pass
       * ``[-lpk, --low_pass_zonal]``      spatial filtering: low-pass
       * ``[-bpk, --band_pass_zonal]``     spatial filtering: band-pass
       * ``[-tidal, --tidal]``             extracts diurnal tide and its
           harmonics
       * ``[-reconstruct, --reconstruct]`` reconstructs the first N
           harmonics
       * ``[-norm, --normalize]``          provides ``-tidal`` result in
           percent amplitude
       * ``[-rs, --regrid_source]``        regrid a target file to match
           a source file
       * ``[-za, --zonal_avg]``            zonally average all variables
           in a file
       * ``[-include, --include]``         only include specific
           variables from the target file
       * ``[-e, --ext]``                   create a new file with a
           unique extension instead of overwriting current file

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

   MarsFiles.combine_files
   MarsFiles.time_shift
   MarsFiles.make_FV3_files
   MarsFiles.do_avg_vars
   MarsFiles.change_vname_longname_unit
   MarsFiles.replace_dims
   MarsFiles.replace_at_index
   MarsFiles.ls2sol_1year



.. py:function:: combine_files(file_list, full_file_list)

   Concatenates sequential output files in chronological order.

   :param file_list: list of file names
   :type file_list: list
   :param full_file_list: list of file names and full paths
   :type full_file_list: list


.. py:function:: time_shift(file_list)

   This function converts the data in diurn files with a time_of_day_XX
   dimension to universal local time.

   :param file_list: list of file names
   :type file_list: list


.. py:function:: make_FV3_files(fpath, typelistfv3, renameFV3=True)

   Make MGCM-like ``average``, ``daily``, and ``diurn`` files.
   Used if call to [``-fv3 --fv3``] is made AND Legacy files are in
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


.. py:function:: do_avg_vars(histfile, newf, avgtime, avgtod, bin_period=5)

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


.. py:function:: change_vname_longname_unit(vname, longname_txt, units_txt)

   Update variable ``name``, ``longname``, and ``units``. This is
   designed to work specifically with LegacyCGM.nc files.

   :param vname: variable name
   :type vname: str
   :param longname_txt: variable description
   :type longname_txt: str
   :param units_txt: variable units
   :type units_txt: str

   :return: variable name and corresponding description and unit


.. py:function:: replace_dims(dims, todflag)

   Replaces dimensions with MGCM-like names. Removes ``time_of_day``.
   This is designed to work specifically with LegacyCGM.nc files.

   :param dims: dimensions of the variable
   :type dims: str
   :param todflag: indicates whether there exists a ``time_of_day``
       dimension
   :type todflag: bool

   :return: new dimension names for the variable


.. py:function:: replace_at_index(tuple_dims, idx, new_name)

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


.. py:function:: ls2sol_1year(Ls_deg, offset=True, round10=True)

   Returns a sol number from the solar longitude.

   :param Ls_deg: solar longitude [°]
   :type Ls_deg: float
   :param offset: if True, force year to start at Ls 0
   :type offset: bool
   :param round10: if True, round to the nearest 10 sols
   :type round10: bool

   :returns: ``Ds`` the sol number

   .. NOTE:: For the moment, this is consistent with 0 <= Ls <=
       359.99, but not for monotically increasing Ls.


