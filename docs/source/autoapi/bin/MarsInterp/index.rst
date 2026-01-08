:py:mod:`bin.MarsInterp`
========================

.. py:module:: bin.MarsInterp

.. autoapi-nested-parse::

   The MarsInterp executable is for interpolating files to pressure or
   altitude coordinates. Options include interpolation to standard
   pressure (``pstd``), standard altitude (``zstd``), altitude above
   ground level (``zagl``), or a custom vertical grid.

   The executable requires:

       * ``[input_file]``          The file to be transformed

   and optionally accepts:

       * ``[-t --interp_type]``    Type of interpolation to perform (altitude, pressure, etc.)
       * ``[-v --vertical_grid]``  Specific vertical grid to interpolate to
       * ``[-incl --include]``     Variables to include in the new interpolated file
       * ``[-ext --extension]``    Custom extension for the new file
       * ``[-print --print_grid]`` Print the vertical grid to the screen

   Third-party Requirements:

       * ``numpy``
       * ``netCDF4``
       * ``argparse``
       * ``os``
       * ``time``
       * ``matplotlib``
       * ``re``
       * ``functools``
       * ``traceback``
       * ``sys``
       * ``amescap``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsInterp.debug_wrapper
   bin.MarsInterp.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsInterp.Cp
   bin.MarsInterp.M_co2
   bin.MarsInterp.R
   bin.MarsInterp.args
   bin.MarsInterp.debug
   bin.MarsInterp.exit_code
   bin.MarsInterp.filepath
   bin.MarsInterp.fill_value
   bin.MarsInterp.g
   bin.MarsInterp.parser
   bin.MarsInterp.rgas


.. py:function:: debug_wrapper(func)

   A decorator that wraps a function with error handling
   based on the --debug flag.
   If the --debug flag is set, it prints the full traceback
   of any exception that occurs. Otherwise, it prints a
   simplified error message.

   :param func: The function to wrap.
   :type  func: function
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


.. py:function:: main()

   Main function for performing vertical interpolation on Mars
   atmospheric model NetCDF files.

   This function processes one or more input NetCDF files,
   interpolating variables from their native vertical coordinate
   (e.g., model pressure levels) to a user-specified standard vertical
   grid (pressure, altitude, or altitude above ground level).
   The interpolation type and grid can be customized via command-line
   arguments.

   Workflow:
       1. Parses command-line arguments for input files, interpolation
       type, custom vertical grid, and other options.
       2. Loads standard vertical grid definitions (pressure, altitude,
       or altitude above ground level) or uses a custom grid.
       3. Optionally prints the vertical grid and exits if requested.
       4. For each input file:
           - Checks file existence.
           - Loads necessary variables (e.g., pk, bk, ps, temperature).
           - Computes the 3D vertical coordinate field for
             interpolation.
           - Creates a new NetCDF output file with updated vertical
             dimension.
           - Interpolates selected variables to the new vertical grid.
           - Copies or interpolates other variables as appropriate.
       5. Handles both regular and diurnal-cycle files, as well as
       FV3-tiled and lat/lon grids.

   Command-line arguments (via `args`):
       - input_file: List of input NetCDF files to process.
       - interp_type: Type of vertical interpolation ('pstd', 'zstd',
         or 'zagl').
       - vertical_grid: Custom vertical grid definition (optional).
       - print_grid: If True, prints the vertical grid and exits.
       - extension: Optional string to append to output filenames.
       - include: List of variable names to include in interpolation.
       - debug: Enable debug output.

   Notes:
       - Requires several helper functions and classes (e.g.,
         section_content_amescap_profile, find_fixedfile, Dataset,
         Ncdf, vinterp).
       - Handles both FV3-tiled and regular lat/lon NetCDF files.
       - Exits with an error message if required files or variables are
         missing.


.. py:data:: Cp
   :value: 735.0

   

.. py:data:: M_co2
   :value: 0.044

   

.. py:data:: R
   :value: 8.314

   

.. py:data:: args

   

.. py:data:: debug

   

.. py:data:: exit_code

   

.. py:data:: filepath

   

.. py:data:: fill_value
   :value: 0.0

   

.. py:data:: g
   :value: 3.72

   

.. py:data:: parser

   

.. py:data:: rgas
   :value: 189.0

   

