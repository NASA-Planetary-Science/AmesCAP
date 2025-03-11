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



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsInterp.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsInterp.Cp
   bin.MarsInterp.M_co2
   bin.MarsInterp.R
   bin.MarsInterp.args
   bin.MarsInterp.filepath
   bin.MarsInterp.fill_value
   bin.MarsInterp.g
   bin.MarsInterp.parser
   bin.MarsInterp.rgas


.. py:function:: main()


.. py:data:: Cp
   :value: 735.0

   

.. py:data:: M_co2
   :value: 0.044

   

.. py:data:: R
   :value: 8.314

   

.. py:data:: args

   

.. py:data:: filepath

   

.. py:data:: fill_value
   :value: 0.0

   

.. py:data:: g
   :value: 3.72

   

.. py:data:: parser

   

.. py:data:: rgas
   :value: 189.0

   

