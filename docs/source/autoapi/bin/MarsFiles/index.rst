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
       * ``[-norm, --normalize]``            Provide ``-tide`` result in % amplitude
       * ``[-regrid, --regrid_XY_to_match]`` Regrid a target file to match a source file
       * ``[-zavg, --zonal_average]``        Zonally average all variables in a file
       * ``[-incl, --include]``              Only include specific variables in a calculation
       * ``[-ext, --extension]``             Create a new file with a unique extension instead of modifying the current file

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
   bin.MarsFiles.out_ext
   bin.MarsFiles.out_ext
   bin.MarsFiles.parser


.. py:class:: ExtAction(*args, ext_content=None, parser=None, **kwargs)


   Bases: :py:obj:`argparse.Action`

   Information about how to convert command line strings to Python objects.

   Action objects are used by an ArgumentParser to represent the information
   needed to parse a single argument from one or more strings from the
   command line. The keyword arguments to the Action constructor are also
   all attributes of Action instances.

   Keyword Arguments:

       - option_strings -- A list of command-line option strings which
           should be associated with this action.

       - dest -- The name of the attribute to hold the created object(s)

       - nargs -- The number of command-line arguments that should be
           consumed. By default, one argument will be consumed and a single
           value will be produced.  Other values include:
               - N (an integer) consumes N arguments (and produces a list)
               - '?' consumes zero or one arguments
               - '*' consumes zero or more arguments (and produces a list)
               - '+' consumes one or more arguments (and produces a list)
           Note that the difference between the default and nargs=1 is that
           with the default, a single value will be produced, while with
           nargs=1, a list containing a single value will be produced.

       - const -- The value to be produced if the option is specified and the
           option uses an action that takes no values.

       - default -- The value to be produced if the option is not specified.

       - type -- A callable that accepts a single string argument, and
           returns the converted value.  The standard Python types str, int,
           float, and complex are useful examples of such callables.  If None,
           str is used.

       - choices -- A container of values that should be allowed. If not None,
           after a command-line argument has been converted to the appropriate
           type, an exception will be raised if it is not a member of this
           collection.

       - required -- True if the action must always be specified at the
           command line. This is only meaningful for optional command-line
           arguments.

       - help -- The help string describing the argument.

       - metavar -- The name to be used for the option's argument with the
           help string. If None, the 'dest' value will be used as the name.

   .. py:method:: format_usage()



.. py:class:: ExtArgumentParser(prog=None, usage=None, description=None, epilog=None, parents=[], formatter_class=HelpFormatter, prefix_chars='-', fromfile_prefix_chars=None, argument_default=None, conflict_handler='error', add_help=True, allow_abbrev=True, exit_on_error=True)


   Bases: :py:obj:`argparse.ArgumentParser`

   Object for parsing command line strings into Python objects.

   Keyword Arguments:
       - prog -- The name of the program (default:
           ``os.path.basename(sys.argv[0])``)
       - usage -- A usage message (default: auto-generated from arguments)
       - description -- A description of what the program does
       - epilog -- Text following the argument descriptions
       - parents -- Parsers whose arguments should be copied into this one
       - formatter_class -- HelpFormatter class for printing help messages
       - prefix_chars -- Characters that prefix optional arguments
       - fromfile_prefix_chars -- Characters that prefix files containing
           additional arguments
       - argument_default -- The default value for all arguments
       - conflict_handler -- String indicating how to handle conflicts
       - add_help -- Add a -h/-help option
       - allow_abbrev -- Allow long options to be abbreviated unambiguously
       - exit_on_error -- Determines whether or not ArgumentParser exits with
           error info when an error occurs

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
   :type vname: str
   :param longname_txt: variable description
   :type longname_txt: str
   :param units_txt: variable units
   :type units_txt: str

   :return: variable name and corresponding description and unit


.. py:function:: concatenate_files(file_list, full_file_list)

   Concatenates sequential output files in chronological order.

   :param file_list: list of file names
   :type file_list: list
   :param full_file_list: list of file names and full paths
   :type full_file_list: list


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


.. py:function:: main()


.. py:function:: make_FV3_files(fpath, typelistfv3, renameFV3=True)

   Make MGCM-like ``average``, ``daily``, and ``diurn`` files.
   Used if call to [``-bin --bin_files``] is made AND Legacy files are in
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


.. py:function:: process_time_shift(file_list)

   This function converts the data in diurn files with a time_of_day_XX
   dimension to universal local time.

   :param file_list: list of file names
   :type file_list: list


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


.. py:function:: replace_dims(dims, todflag)

   Replaces dimensions with MGCM-like names. Removes ``time_of_day``.
   This is designed to work specifically with LegacyCGM.nc files.

   :param dims: dimensions of the variable
   :type dims: str
   :param todflag: indicates whether there exists a ``time_of_day``
       dimension
   :type todflag: bool

   :return: new dimension names for the variable


.. py:function:: split_files(file_list, split_dim)

   Extracts variables in the file along the time dimension, unless
   other dimension is specified (lev, lat, or lon).

   :param file_list: list of file names
   :type split_dim: dimension along which to perform extraction
   :returns: new file with sliced dimensions


.. py:data:: all_args

   

.. py:data:: args

   

.. py:data:: out_ext

   

.. py:data:: out_ext

   

.. py:data:: parser

   

