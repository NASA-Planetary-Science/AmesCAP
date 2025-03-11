:py:mod:`bin.MarsPull`
======================

.. py:module:: bin.MarsPull

.. autoapi-nested-parse::

   The MarsPull executable is for querying data from the Mars Climate Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

   The executable requires 2 arguments:

       * The directory from which to pull data from, AND
       * ``[-ls --ls]``      The desired solar longitude(s), OR
       * ``[-f --filename]`` The name(s) of the desired file(s)

   Third-party Requirements:

       * ``numpy``
       * ``argparse``
       * ``requests``

   List of Functions:

       * download - Queries the requested file from the NAS Data Portal.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsPull.download
   bin.MarsPull.file_list
   bin.MarsPull.main



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsPull.Ls_end
   bin.MarsPull.Ls_ini
   bin.MarsPull.args
   bin.MarsPull.parser
   bin.MarsPull.saveDir


.. py:function:: download(url, filename)

   Downloads a file from the NAS Data Portal (data.nas.nasa.gov).

   :param url: The url to download, e.g   'https://data.nas.nasa.gov/legacygcm/download_data.php?file=/legacygcmdata/LegacyGCM_Ls000_Ls004.nc'
   :type url: str

   :param filename: The local filename e.g  '/lou/la4/akling/Data/LegacyGCM_Ls000_Ls004.nc'
   :type filename: str

   :return: The requested file(s), downloaded and saved to the current     directory.

   :raises FileNotFoundError: A file-not-found error.



.. py:function:: file_list(list_of_files)


.. py:function:: main()


.. py:data:: Ls_end

   

.. py:data:: Ls_ini

   

.. py:data:: args

   

.. py:data:: parser

   

.. py:data:: saveDir

   

