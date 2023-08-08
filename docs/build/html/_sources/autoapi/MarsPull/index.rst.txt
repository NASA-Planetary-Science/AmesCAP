:py:mod:`MarsPull`
==================

.. py:module:: MarsPull

.. autoapi-nested-parse::

   MarsPull

   The MarsPull executable is for querying data from the Mars Climate
   Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on
   the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

   The executable requires two arguments:
       * [-id --id]      the simulation identifier, AND
       * [-ls --ls]      the desired solar longitude(s), OR
       * [-f -filename]  the name(s) of the desired file(s)

   Third-party Requirements:
       * numpy
       * argparse
       * requests

   List of Functions:
       * download - queries the requested file from the NAS Data Portal



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsPull.download



.. py:function:: download(fName, simID)

   Downloads a file from the MCMC Legacy GCM directory on the NAS Data
   Portal (data.nas.nasa.gov).

   This function specifies the file to download by appending to the 
   URL the subdirectory, indicated by the user-specified 
   simulation identifier [-id --id], and the name of the file. The 
   file name is either provided by the user directly using 
   [-f --filename] or determined based on the user-specified solar 
   longitude [-ls --ls].

   Parameters
   ----------
   simID : str
       The simulation identifier, i.e., the name of the directory
       to query at:
       https://data.nas.nasa.gov/mcmc/data_legacygcm.php

   fName : str
       The name of the file to download

   Raises
   ------
   rsp.status_code
       A file-not-found error

   Returns
   -------
   The requested file(s), downloaded and saved to the current directory


