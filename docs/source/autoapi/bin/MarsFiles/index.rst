:py:mod:`bin.MarsFiles`
=======================

.. py:module:: bin.MarsFiles

.. autoapi-nested-parse::

   The MarsFiles executable has functions for manipulating entire files.
   The capabilities include time-shifting, binning, and regridding data,
   as well as band pass filtering, tide analysis, zonal averaging, and
   extracting variables from files.

   The executable requires:

       * ``[input_file]``                  The file for manipulation

   and optionally accepts:

       * ``[-bin, --bin_files]``             Produce MGCM 'fixed', 'diurn', 'average' and 'daily' files from Legacy output
       * ``[-c, --concatenate]``             Combine sequential files of the same type into one file
       * ``[-split, --split]``               Split file along a specified dimension or extracts slice at one point along the dim
       * ``[-t, --time_shift]``              Apply a time-shift to 'diurn' files
       * ``[-ba, --bin_average]``            Bin MGCM 'daily' files like 'average' files
       * ``[-bd, --bin_diurn]``              Bin MGCM 'daily' files like 'diurn' files
       * ``[-hpt, --high_pass_temporal]``    Temporal filter: high-pass
       * ``[-lpt, --low_pass_temporal]``     Temporal filter: low-pass
       * ``[-bpt, --band_pass_temporal]``    Temporal filter: band-pass
       * ``[-trend, --add_trend]``           Return amplitudes only (use with temporal filters)
       * ``[-hps, --high_pass_spatial]``     Spatial filter: high-pass
       * ``[-lps, --low_pass_spatial]``      Spatial filter: low-pass
       * ``[-bps, --band_pass_spatial]``     Spatial filter: band-pass
       * ``[-tide, --tide_decomp]``          Extract diurnal tide and its harmonics
       * ``[-recon, --reconstruct]``         Reconstruct the first N harmonics
       * ``[-norm, --normalize]``            Provide ``-tide`` result in percent amplitude
       * ``[-prop, --prop_tides]``           Extract propagating tide harmonics
       * ``[-regrid, --regrid_XY_to_match]`` Regrid a target file to match a source file
       * ``[-zavg, --zonal_average]``        Zonally average all variables in a file
       * ``[-incl, --include]``              Only include specific variables in a calculation
       * ``[-ext, --extension]``             Create a new file with a unique extension instead of modifying the current file

   Third-party Requirements:

       * ``sys``
       * ``argparse``
       * ``os``
       * ``warnings``
       * ``re``
       * ``numpy``
       * ``netCDF4``
       * ``shutil``
       * ``functools``
       * ``traceback``



Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   bin.MarsFiles.ExtAction
   bin.MarsFiles.ExtArgumentParser



Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsFiles.change_vname_longname_unit
   bin.MarsFiles.concatenate_files
   bin.MarsFiles.debug_wrapper
   bin.MarsFiles.do_avg_vars
   bin.MarsFiles.ls2sol_1year
   bin.MarsFiles.main
   bin.MarsFiles.make_FV3_files
   bin.MarsFiles.process_time_shift
   bin.MarsFiles.replace_at_index
   bin.MarsFiles.replace_dims
   bin.MarsFiles.split_files



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsFiles.all_args
   bin.MarsFiles.args
   bin.MarsFiles.debug
   bin.MarsFiles.exit_code
   bin.MarsFiles.out_ext
   bin.MarsFiles.out_ext
   bin.MarsFiles.parser


.. py:class:: ExtAction(*args, ext_content=None, parser=None, **kwargs)


   Bases: :py:obj:`argparse.Action`

   Custom action for argparse to handle file extensions.

   This action is used to add an extension to the output file.

   :param ext_content: The content to be added to the file name.
   :type  ext_content: str
   :param parser: The parser instance to which this action belongs.
   :type  parser: argparse.ArgumentParser
   :param args: Additional positional arguments.
   :type  args: tuple
   :param kwargs: Additional keyword arguments.
   :type  kwargs: dict
   :param ext_content: The content to be added to the file name
   :type  ext_content: str
   :return: None
   :rtype:  None
   :raises ValueError: If ext_content is not provided.
   :raises TypeError: If parser is not an instance of argparse.ArgumentParser.
   :raises Exception: If an error occurs during the action.
   :raises AttributeError: If the parser does not have the specified attribute.
   :raises ImportError: If the parser cannot be imported.
   :raises RuntimeError: If the parser cannot be run.
   :raises KeyError: If the parser does not have the specified key.
   :raises IndexError: If the parser does not have the specified index.
   :raises IOError: If the parser cannot be opened.
   :raises OSError: If the parser cannot be accessed.
   :raises EOFError: If the parser cannot be read.
   :raises MemoryError: If the parser cannot be allocated.
   :raises OverflowError: If the parser cannot be overflowed.

   .. py:method:: __call__(parser, namespace, values, option_string=None)


   .. py:method:: __repr__()

      Return repr(self).


   .. py:method:: format_usage()



.. py:class:: ExtArgumentParser(prog=None, usage=None, description=None, epilog=None, parents=[], formatter_class=HelpFormatter, prefix_chars='-', fromfile_prefix_chars=None, argument_default=None, conflict_handler='error', add_help=True, allow_abbrev=True, exit_on_error=True)


   Bases: :py:obj:`argparse.ArgumentParser`

   Custom ArgumentParser that handles file extensions for output files.

   This class extends the argparse.ArgumentParser to add functionality
   for handling file extensions in the command-line arguments.

   :param args: Additional positional arguments.
   :type  args: tuple
   :param kwargs: Additional keyword arguments.
   :type  kwargs: dict
   :param ext_content: The content to be added to the file name.
   :type  ext_content: str
   :param parser: The parser instance to which this action belongs.
   :type  parser: argparse.ArgumentParser
   :return: None
   :rtype:  None
   :raises ValueError: If ext_content is not provided.
   :raises TypeError: If parser is not an instance of argparse.ArgumentParser.
   :raises Exception: If an error occurs during the action.
   :raises AttributeError: If the parser does not have the specified attribute.
   :raises ImportError: If the parser cannot be imported.
   :raises RuntimeError: If the parser cannot be run.
   :raises KeyError: If the parser does not have the specified key.
   :raises IndexError: If the parser does not have the specified index.
   :raises IOError: If the parser cannot be opened.
   :raises OSError: If the parser cannot be accessed.
   :raises EOFError: If the parser cannot be read.
   :raises MemoryError: If the parser cannot be allocated.
   :raises OverflowError: If the parser cannot be overflowed.

   .. py:method:: __repr__()

      Return repr(self).


   .. py:method:: add_argument(*args, **kwargs)

      add_argument(dest, ..., name=value, ...)
      add_argument(option_string, option_string, ..., name=value, ...)


   .. py:method:: add_argument_group(*args, **kwargs)


   .. py:method:: add_mutually_exclusive_group(**kwargs)


   .. py:method:: add_subparsers(**kwargs)


   .. py:method:: convert_arg_line_to_args(arg_line)


   .. py:method:: error(message)

      error(message: string)

      Prints a usage message incorporating the message to stderr and
      exits.

      If you override this in a subclass, it should not return -- it
      should either exit or raise an exception.


   .. py:method:: exit(status=0, message=None)


   .. py:method:: format_help()


   .. py:method:: format_usage()


   .. py:method:: get_default(dest)


   .. py:method:: parse_args(*args, **kwargs)


   .. py:method:: parse_intermixed_args(args=None, namespace=None)


   .. py:method:: parse_known_args(args=None, namespace=None)


   .. py:method:: parse_known_intermixed_args(args=None, namespace=None)


   .. py:method:: print_help(file=None)


   .. py:method:: print_usage(file=None)


   .. py:method:: register(registry_name, value, object)


   .. py:method:: set_defaults(**kwargs)



.. py:function:: change_vname_longname_unit(vname, longname_txt, units_txt)

   Update variable ``name``, ``longname``, and ``units``. This is
   designed to work specifically with LegacyCGM.nc files.

   :param vname: variable name
   :type  vname: str
   :param longname_txt: variable description
   :type  longname_txt: str
   :param units_txt: variable units
   :type  units_txt: str
   :return: variable name and corresponding description and unit
       (e.g. ``vname = "ps"``)
   :rtype:  tuple
   :raises KeyError: if the required variables are not found
   :raises ValueError: if the required dimensions are not found
   :raises AttributeError: if the required attributes are not found

   :note::
       The ``diurn`` file is created by binning the Legacy files.
       The ``average`` and ``daily`` files are created by
       averaging over the ``diurn`` file.


.. py:function:: concatenate_files(file_list, full_file_list)

   Concatenates sequential output files in chronological order.

   :param file_list: list of file names
   :type  file_list: list
   :param full_file_list: list of file names and full paths
   :type  full_file_list: list
   :returns: None
   :rtype:  None
   :raises OSError: If the file cannot be removed.
   :raises IOError: If the file cannot be moved.
   :raises Exception: If the file cannot be opened.
   :raises ValueError: If the file cannot be accessed.
   :raises TypeError: If the file is not of the correct type.
   :raises IndexError: If the file does not have the correct index.
   :raises KeyError: If the file does not have the correct key.


.. py:function:: debug_wrapper(func)

   A decorator that wraps a function with error handling
   based on the --debug flag.

   If the --debug flag is set, it prints the full traceback
   of any exception that occurs. Otherwise, it prints a
   simplified error message.

   :param func: The function to wrap.
   :type   func: function
   :return: The wrapped function.
   :rtype:  function
   :raises Exception: If an error occurs during the function call.
   :raises TypeError: If the function is not callable.
   :raises ValueError: If the function is not found.
   :raises NameError: If the function is not defined.
   :raises AttributeError: If the function does not have the
       specified attribute.
   :raises ImportError: If the function cannot be imported.
   :raises RuntimeError: If the function cannot be run.
   :raises KeyError: If the function does not have the
       specified key.
   :raises IndexError: If the function does not have the
       specified index.
   :raises IOError: If the function cannot be opened.
   :raises OSError: If the function cannot be accessed.
   :raises EOFError: If the function cannot be read.
   :raises MemoryError: If the function cannot be allocated.
   :raises OverflowError: If the function cannot be overflowed.
   :raises ZeroDivisionError: If the function cannot be divided by zero.
   :raises StopIteration: If the function cannot be stopped.
   :raises KeyboardInterrupt: If the function cannot be interrupted.
   :raises SystemExit: If the function cannot be exited.
   :raises AssertionError: If the function cannot be asserted.


.. py:function:: do_avg_vars(histfile, newf, avgtime, avgtod, bin_period=5)

   Performs a time average over all fields in a file.

   :param histfile: file to perform time average on
   :type  histfile: str
   :param newf: path to target file
   :type  newf: str
   :param avgtime: whether ``histfile`` has averaged fields
       (e.g., ``atmos_average``)
   :type  avgtime: bool
   :param avgtod: whether ``histfile`` has a diurnal time dimenion
       (e.g., ``atmos_diurn``)
   :type  avgtod: bool
   :param bin_period: the time binning period if `histfile` has
       averaged fields (i.e., if ``avgtime==True``), defaults to 5
   :type  bin_period: int, optional
   :return: a time-averaged file
   :rtype:  None
   :raises KeyError: if the required variables are not found
   :raises ValueError: if the required dimensions are not found
   :raises AttributeError: if the required attributes are not found

   :note::
       The ``diurn`` file is created by binning the Legacy files.
       The ``average`` and ``daily`` files are created by
       averaging over the ``diurn`` file.


.. py:function:: ls2sol_1year(Ls_deg, offset=True, round10=True)

   Returns a sol number from the solar longitude.

   This is consistent with the MGCM model. The Ls is the solar
   longitude in degrees. The sol number is the number of sols since
   the perihelion (Ls = 250.99 degrees).

   :param Ls_deg: solar longitude [°]
   :type  Ls_deg: float
   :param offset: if True, force year to start at Ls 0
   :type  offset: bool
   :param round10: if True, round to the nearest 10 sols
   :type  round10: bool
   :returns: ``Ds`` the sol number
   :rtype:  float
   :raises ValueError: if the required variables are not found
   :raises KeyError: if the required variables are not found
   :raises AttributeError: if the required attributes are not found

   ..note::
       This is consistent with 0 <= Ls <= 359.99, but not for
       monotically increasing Ls.


.. py:function:: main()

   Main entry point for MarsFiles data processing utility.

   This function processes input NetCDF or legacy MGCM files according
   to command-line arguments. It supports a variety of operations,
   including:

   - Conversion of legacy MGCM files to FV3 format.
   - Concatenation and splitting of NetCDF files along specified
     dimensions.
   - Temporal binning of daily files into average or diurnal files.
   - Temporal filtering (high-pass, low-pass, band-pass) using spectral
     methods.
   - Spatial (zonal) filtering and decomposition using spherical
     harmonics.
   - Tidal analysis and harmonic decomposition of diurnal files.
   - Regridding of data to match a target NetCDF file's grid.
   - Zonal averaging of variables over longitude.

   The function handles file path resolution, argument validation, and
   orchestrates the appropriate processing routines based on
   user-specified options. Output files are written in NetCDF format,
   with new variables and dimensions created as needed.

   Global Variables:
       data_dir (str): The working directory for input/output files.
   Arguments:
       None directly. Uses global 'args' parsed from command-line.
   Returns:
       None. Outputs are written to disk.
   Raises:
       SystemExit: For invalid argument combinations or processing
       errors.


.. py:function:: make_FV3_files(fpath, typelistfv3, renameFV3=True)

   Make MGCM-like ``average``, ``daily``, and ``diurn`` files.

   Used if call to [``-bin --bin_files``] is made AND Legacy files are
   in netCDF format (not fort.11).

   :param fpath: Full path to the Legacy netcdf files
   :type  fpath: str
   :param typelistfv3: MGCM-like file type: ``average``, ``daily``,
       or ``diurn``
   :type  typelistfv3: list
   :param renameFV3: Rename the files from Legacy_LsXXX_LsYYY.nc to
       ``XXXXX.atmos_average.nc`` following MGCM output conventions
   :type  renameFV3: bool
   :return: The MGCM-like files: ``XXXXX.atmos_average.nc``,
       ``XXXXX.atmos_daily.nc``, ``XXXXX.atmos_diurn.nc``.
   :rtype:  None

   :note::
       The ``average`` and ``daily`` files are created by
       averaging over the ``diurn`` file. The ``diurn`` file is
       created by binning the Legacy files.

   :note::
       The ``diurn`` file is created by binning the Legacy files.


.. py:function:: process_time_shift(file_list)

   Converts diurn files to local time.

   This function is used to convert diurn files to local time.

   :param file_list: list of file names
   :type  file_list: list
   :returns: None
   :rtype:  None
   :raises OSError: If the file cannot be removed.
   :raises IOError: If the file cannot be moved.
   :raises Exception: If the file cannot be opened.
   :raises ValueError: If the file cannot be accessed.
   :raises TypeError: If the file is not of the correct type.
   :raises IndexError: If the file does not have the correct index.


.. py:function:: replace_at_index(tuple_dims, idx, new_name)

   Replaces the dimension at the given index with a new name.

   If ``new_name`` is None, the dimension is removed.
   This is designed to work specifically with LegacyCGM.nc files.

   :param tuple_dims: the dimensions as tuples e.g. (``pfull``,
       ``nlat``, ``nlon``)
   :type  tuple_dims: tuple
   :param idx: index indicating axis with the dimensions to update
       (e.g. ``idx = 1``  for ``nlat``)
   :type  idx: int
   :param new_name: new dimension name (e.g. ``latitude``)
   :type  new_name: str
   :return: updated dimensions
   :rtype:  tuple
   :raises KeyError: if the required variables are not found
   :raises ValueError: if the required dimensions are not found
   :raises AttributeError: if the required attributes are not found


.. py:function:: replace_dims(dims, todflag)

   Replaces dimensions with MGCM-like names. Removes ``time_of_day``.
   This is designed to work specifically with LegacyCGM.nc files.

   :param dims: dimensions of the variable
   :type  dims: str
   :param todflag: indicates whether there exists a ``time_of_day``
       dimension
   :type  todflag: bool
   :return: new dimension names for the variable
   :rtype:  tuple
   :raises KeyError: if the required variables are not found
   :raises ValueError: if the required dimensions are not found
   :raises AttributeError: if the required attributes are not found


.. py:function:: split_files(file_list, split_dim)

   Extracts variables in the file along the specified dimension.

   The function creates a new file with the same name as the input
   file, but with the specified dimension sliced to the requested
   value or range of values. The new file is saved in the same
   directory as the input file.

   The function also checks the bounds of the requested dimension
   against the actual values in the file. If the requested value
   is outside the bounds of the dimension, an error message is
   printed and the function exits.

   The function also checks if the requested dimension is valid
   (i.e., time, lev, lat, or lon). If the requested dimension is
   invalid, an error message is printed and the function exits.

   The function also checks if the requested dimension is a
   single dimension (i.e., areo). If the requested dimension is
   a single dimension, the function removes all single dimensions
   from the areo dimension (i.e., scalar_axis) before performing
   the extraction.

   :param file_list: list of file names
   :type  split_dim: dimension along which to perform extraction
   :returns: new file with sliced dimensions
   :rtype:  None
   :raises OSError: If the file cannot be removed.
   :raises IOError: If the file cannot be moved.
   :raises Exception: If the file cannot be opened.
   :raises ValueError: If the file cannot be accessed.
   :raises TypeError: If the file is not of the correct type.
   :raises IndexError: If the file does not have the correct index.


.. py:data:: all_args

   

.. py:data:: args

   

.. py:data:: debug

   

.. py:data:: exit_code

   

.. py:data:: out_ext

   

.. py:data:: out_ext

   

.. py:data:: parser

   

