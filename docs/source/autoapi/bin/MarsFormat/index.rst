bin.MarsFormat
==============

.. py:module:: bin.MarsFormat

.. autoapi-nested-parse::

   The MarsFormat executable is for converting non-MGCM data, such as that
   from EMARS, OpenMARS, PCM, and MarsWRF, into MGCM-like netCDF data
   products. The MGCM is the NASA Ames Mars Global Climate Model developed
   and maintained by the Mars Climate Modeling Center (MCMC). The MGCM
   data repository is available at data.nas.nasa.gov/mcmc.

   The executable requires two arguments:

       * ``[input_file]``         The file to be transformed
       * ``[-gcm --gcm_name]``    The GCM from which the file originates

   and optionally accepts:

       * ``[-rn --retain_names]`` Preserve original variable and dimension names
       * ``[-ba, --bin_average]`` Bin non-MGCM files like 'average' files
       * ``[-bd, --bin_diurn]``   Bin non-MGCM files like 'diurn' files

   Third-party requirements:

       * ``numpy``
       * ``netCDF4``
       * ``sys``
       * ``argparse``
       * ``os``
       * ``re``
       * ``functools``
       * ``traceback``
       * ``xarray``
       * ``amescap``



Attributes
----------

.. autoapisummary::

   bin.MarsFormat.args
   bin.MarsFormat.debug
   bin.MarsFormat.exit_code
   bin.MarsFormat.parser
   bin.MarsFormat.path2data
   bin.MarsFormat.ref_press


Functions
---------

.. autoapisummary::

   bin.MarsFormat.debug_wrapper
   bin.MarsFormat.get_time_dimension_name
   bin.MarsFormat.main
   bin.MarsFormat.marswrf_eta_mid_to_interfaces
   bin.MarsFormat.marswrf_fit_eta
   bin.MarsFormat.marswrf_time_axis


Module Contents
---------------

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


.. py:function:: get_time_dimension_name(DS, model)

   Find the time dimension name in the dataset.

   Updates the model object with the correct dimension name.

   :param DS: The xarray Dataset
   :type  DS: xarray.Dataset
   :param model: Model object with dimension information
   :type  model: object
   :return: The actual time dimension name found
   :rtype:  str
   :raises KeyError: If no time dimension is found
   :raises ValueError: If the model object is not defined
   :raises TypeError: If the dataset is not an xarray Dataset
   :raises AttributeError: If the model object does not have the
       specified attribute
   :raises ImportError: If the xarray module cannot be imported


.. py:function:: main()

   Main processing function for MarsFormat.

   This function processes NetCDF files from various Mars General
   Circulation Models (GCMs)
   including MarsWRF, OpenMars, PCM, and EMARS, and reformats them for
   use in the AmesCAP
   framework.

   It performs the following operations:
       - Validates the selected GCM type and input files.
       - Loads NetCDF files and reads model-specific variable and
       dimension mappings.
       - Applies model-specific post-processing, including:
           - Unstaggering variables (for MarsWRF and EMARS).
           - Creating and orienting pressure coordinates (pfull, phalf,
           ak, bk).
           - Standardizing variable and dimension names.
           - Converting longitude ranges to 0-360 degrees east.
           - Adding scalar axes where required.
           - Handling vertical dimension orientation, especially for
           PCM files.
       - Optionally performs time binning:
           - Daily, average (over N sols), or diurnal binning.
           - Ensures correct time units and bin sizes.
           - Preserves or corrects vertical orientation after binning.
       - Writes processed datasets to new NetCDF files with appropriate
       naming conventions.

   Args:
       None. Uses global `args` for configuration and file selection.

   Raises:
       KeyError: If required dimensions or variables are missing in
       the input files.
       ValueError: If dimension swapping fails for PCM files.
       SystemExit: If no valid GCM type is specified.

   Outputs:
       Writes processed NetCDF files to disk, with suffixes indicating
       the type of processing
       (e.g., _daily, _average, _diurn, _nat).

   Note:
       This function assumes the presence of several helper functions
       and global variables,
       such as `read_variable_dict_amescap_profile`,
       `get_time_dimension_name`,  `reset_FV3_names`, and color
       constants for printing.
       


.. py:function:: marswrf_eta_mid_to_interfaces(eta_mid)

   Interface (w-level) eta from mid (mass-level) eta, surface first:
   eta_w[0] = 1 and eta_w[k+1] = 2*eta_u[k] - eta_w[k]. The last
   value is clipped to >= 0. Input may be in either vertical order;
   output follows the input order.


.. py:function:: marswrf_fit_eta(p3d, ps, tdim_axis=0)

   Recover the WRF eta coordinate from a 3D mid-level pressure field
   and surface pressure, for reduced files that did not keep
   ZNU/ZNW/P_TOP. WRF: p = p_top + eta*(ps - p_top), so per level a
   least-squares fit p = a_k + b_k*ps over all columns and times
   gives eta_k = b_k and p_top = a_k/(1 - b_k).

   p3d: array [time, lev, lat, lon] (any vertical order)
   ps : array [time, lat, lon]
   Returns (eta_mid, p_top).


.. py:function:: marswrf_time_axis(DS, model)

   Build the time coordinate [sols] for a planetWRF file from the
   best available clock and return (DS, source_description).

   Priority:
     1. 'Times' character array 'YYYY-DDDDD_HH:MM:SS' (absolute,
        carries the time of sol)
     2. 'XTIME' [minutes since SIMULATION_START_DATE] plus the sol of
        SIMULATION_START_DATE
     3. 'JULIAN' (fractional sol of year) plus MODEL_MARS_YEAR
     4. START_DATE global attribute plus a uniform spacing inferred
        from the L_S slope (reduced files that kept only L_S)

   The result is stored under model.time (normally 'XTIME') so that
   the rest of MarsFormat sees one time variable with units in days.
   Ls is left in model.areo ('L_S'); if it is missing it is computed
   from the sol axis with sol2ls().


.. py:data:: args
   :value: None


.. py:data:: debug

.. py:data:: exit_code
   :value: None


.. py:data:: parser

.. py:data:: path2data

.. py:data:: ref_press
   :value: 725


