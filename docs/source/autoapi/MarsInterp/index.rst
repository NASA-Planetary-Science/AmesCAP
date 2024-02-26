:py:mod:`MarsInterp`
====================

.. py:module:: MarsInterp

.. autoapi-nested-parse::

   The MarsInterp executable is for interpolating files to pressure or altitude coordinates. Options include interpolation to standard pressure (``pstd``), standard altitude (``zstd``), altitude above ground level (``zagl``), or a custom vertical grid.

   The executable requires:
       * ``[input_file]``          the file to be transformed

   and optionally accepts:
       * ``[-t --type]``           type of interpolation to perform         (altitude, pressure, etc.)
       * ``[-l --level]``          specific vertical grid to interpolate to
       * ``[-include --include]``  variables to include in the new         interpolated file
       * ``[-e --ext]``            custom extension for the new file
       * ``[-g --grid]``           print the vertical grid chosen by         [-l --level] to the screen


   Third-party Requirements:
       * ``numpy``
       * ``netCDF4``
       * ``argparse``
       * ``os``
       * ``time``
       * ``matplotlib``



