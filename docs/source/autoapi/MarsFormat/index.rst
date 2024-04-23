:py:mod:`MarsFormat`
====================

.. py:module:: MarsFormat

.. autoapi-nested-parse::

   The MarsFormat executable is a routine that transforms non-MGCM model
   output into MGCM-like model output for compatibility with CAP.

   MarsFormat changes variable names, dimension names, dimension order,
   and units to the configuration expected by CAP. In some cases, such as
   for MarsWRF, variables are derived and regridded onto a standard grid.

   The executable requires:
       * ``[input_file]``              the file to be transformed

   and optionally accepts:
       * ``[-openmars --openmars]``    convert openMars data to MGCM format
       * ``[-marswrf --marswrf]``      convert MarsWRF data to MGCM format

   Third-party Requirements:
       * ``numpy``
       * ``os``
       * ``argparse``
       * ``xarray``



