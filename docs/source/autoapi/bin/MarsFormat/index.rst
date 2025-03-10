:py:mod:`bin.MarsFormat`
========================

.. py:module:: bin.MarsFormat

.. autoapi-nested-parse::

   The MarsFormat executable is for converting non-MGCM data, such as thatfrom EMARS, OpenMARS, PCM, and MarsWRF, into MGCM-like netCDF data products. The MGCM is the NASA Ames Mars Global Climate Model developedand maintained by the Mars Climate Modeling Center (MCMC). The MGCM data repository is available at data.nas.nasa.gov/mcmc.

   The executable requires 2 arguments:
       * ``[input_file]``         The file to be transformed
       * ``[-gcm --gcm_name]``    The GCM from which the file originates
       
   and optionally accepts:
       * ``[-rn --retain_names]`` Preserve original variable and dimension names
       * ``[-ba, --bin_average]`` Bin non-MGCM files like 'average' files
       * ``[-bd, --bin_diurn]``   Bin non-MGCM files like 'diurn' files
       
   Third-party Requirements:
       * ``numpy``
       * ``netCDF4``
       * ``sys``
       * ``argparse``
       * ``os``

   List of Functions:
       * download - Queries the requested file from the NAS Data Portal.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsFormat.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsFormat.args
   bin.MarsFormat.parser
   bin.MarsFormat.path2data
   bin.MarsFormat.ref_press


.. py:function:: main()


.. py:data:: args

   

.. py:data:: parser

   

.. py:data:: path2data

   

.. py:data:: ref_press
   :value: 725

   

