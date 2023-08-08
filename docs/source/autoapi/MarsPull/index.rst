:orphan:

:py:mod:`MarsPull`
==================

.. py:module:: MarsPull


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsPull.download



.. py:function:: download(url, filename)

   Downloads a file from  https://data.nas.nasa.gov.

   The file to download is specified by appending the above URL with
   the legacy gcm subdirectory + the filename. The filename can be 
   provided by the user directly or determined based on the 
   user-requested solar longitude (Ls). The simulation identifier (ID)
   must always be provided.

   Parameters
   ----------
   URL: str
       The URL to download from. This is built from:
       https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/
       by appending the simulation ID to the end of the URL.
   filename: str
       The name of the file to download

   Raises
   ------
   rsp.status_code
       A file-not-found error

   Returns
   -------
   downloaded file


