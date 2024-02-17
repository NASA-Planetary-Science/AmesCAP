:py:mod:`MarsFormat`
====================

.. py:module:: MarsFormat

.. autoapi-nested-parse::

   The MarsFormat executable is a routine that transforms non-MGCM model
   output into MGCM-like model output for compatibility with CAP.

   MarsFormat changes variable names, dimension names, dimension order,
   and units to the configuration expected by CAP. In some cases, such as
   for MarsWRF, variables are derived and regridded onto a standard grid.

   The executable requires 1 argument:
       * ``[input_file]``      the file to be transformed

   and optionally accepts 2 arguments:
       * ``[-openmars --openmars]``    convert openMars data to MGCM format
       * ``[-marswrf --marswrf]``      convert MarsWRF data to MGCM format

   Third-party Requirements:
       * ``numpy``
       * ``os``
       * ``argparse``
       * ``xarray``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsFormat.marswrf_to_mgcm
   MarsFormat.openmars_to_mgcm



.. py:function:: marswrf_to_mgcm(DS)

   Converts variables in MarsWRF output files to MGCM-like format and
   derives variables from their perturbations.

   WRF data is output to dimensions: [time, pfull, lat, lon] or
   [t,z,y,x], just like MGCM data. Some WRF variables are on staggered
   grids, referred to in the comments using ' (prime; like y').

   The dimensions of the native WRF variables can be:

   ======== ==================== ========================
   Variable MarsWRF Dimensions   MGCM Equivalent Variable
   ======== ==================== ========================
   ``t``    ``time``             ``time``
   ``z``    ``bottom_top``       ``pfull``
   ``z'``   ``bottom_top_stag``  ``phalf``
   ``y``    ``south_north``      ``lat``
   ``y'``   ``south_north_stag``
   ``x``    ``west_east``        ``lon``
   ``x'``   ``west_east_stag``
   ======== ==================== ========================

   The variables transferred or derived from WRF output and piped to
   the MGCM daily file are listed below.

   ========== =============== ========== ================================
   MarsWRF    MGCM Equiv.     Units      Notes
   ========== =============== ========== ================================
   ``XTIME``  ``time``        days       converted from minutes to                                           days since simulation start
   ``L_S``    ``areo``        degree     
   ``PSFC``   ``ps``          Pa         
   ``XLONG``  ``lon``         degree E   
   ``XLAT``   ``lat``         degree N   
   ``HGT``    ``zsurf``       meters     
   ``U``      ``ucomp``       m/s        Requires interpolation to a                                           regular grid
   ``V``      ``vcomp``       m/s        Requires interpolation to a                                           regular grid
   ``W``      ``w``           m/s        Requires interpolation to a                                           regular grid
   ``H2OICE`` ``h2o_ice_sfc`` kg/m2      
   ``CO2ICE`` ``co2_ice_sfc`` kg/m2      
   ``ZNW``    ``bk``                     
   ``TSK``    ``ts``          K          
   ``P_TOP``  ``pk[0]``       Pa         model top pressure
   ========== =============== ========== ================================

   :param DS: The dataset created by xarray when it opens the         user-supplied input file.
   :type DS: xarray dataset

   Returns
   -------
   :return: ``var_dict`` Dictionary with variable names as keys and a        list of attributes[values, dimensions, longname, units] as         values.

       ``time`` (array) Minutes since simulation start

       ``lat`` (array) Latitude on a regular grid

       ``lon`` (array) Longitude on a regular grid

       ``phalf`` (array) Half pressure levels

       ``pfull`` (array) Full pressure levels



.. py:function:: openmars_to_mgcm(DS)

   Converts variables in openMars output files to MGCM-like format.

   openMars data is similar to MGCM data already. This function derives
   pfull and phalf but otherwise only needs to rename variables and
   update units, longnames, and dimensions to match MGCM output.

   :param DS: The dataset created by xarray when it opens the         user-supplied input file.
   :type DS: xarray dataset

   Returns
   -------
   :return: ``var_dict`` Dictionary with variable names as keys and a        list of attributes[values, dimensions, longname, units] as         values.

       ``time`` (array) Minutes since simulation start

       ``lat`` (array) Latitude on a regular grid

       ``lon`` (array) Longitude on a regular grid

       ``phalf`` (array) Half pressure levels

       ``pfull`` (array) Full pressure levels



