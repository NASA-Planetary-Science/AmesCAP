:py:mod:`MarsPull`
==================

.. py:module:: MarsPull

.. autoapi-nested-parse::

   The MarsPull executable is for querying data from the Mars Climate Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

   The executable requires 2 arguments:
       * ``[-id --id]``      The simulation identifier, AND
       * ``[-ls --ls]``      the desired solar longitude(s), OR
       * ``[-f --filename]`` the name(s) of the desired file(s).

   Third-party Requirements:
       * ``numpy``
       * ``sys``
       * ``argparse``
       * ``os``
       * ``requests``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsPull.download



.. py:function:: download(file_name, simulation_id)

   Downloads a file from the MCMC Legacy GCM directory on the NAS Data
   Portal (data.nas.nasa.gov).

   This function specifies the file to download by appending to the     URL to the subdirectory, indicated by the user-specified     simulation identifier [-id --id], and the name of the file. The     file name is either provided by the user directly using     [-f --filename] or determined based on the user-specified solar         longitude [-ls --ls].

   Parameters
   ----------
   simulation_id : str
       The simulation identifier, i.e., the name of the directory to         query from: https://data.nas.nasa.gov/mcmc/data_legacygcm.php

   file_name : str
       The name of the file to download.

   Raises
   ------
   A file-not-found error.

   Returns
   -------
   The requested file(s), downloaded & saved to the current directory.


