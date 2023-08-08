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

   Downloads a file from data.nas.nasa.gov

   This function specifies the file to download by appending to the 
   URL the subdirectory, indicated by the user-specified 
   simulation identifier [-id --id], and the name of the file. The 
   filename is either provided by the user directly using 
   [-f --filename] or determined based on the user-specified solar 
   longitude [-ls --ls].

   Parameters
   ----------
   URL: str
       The URL of the file to download. This is built from:
       https://data.nas.nasa.gov/mcmc/legacygcmdata
       by appending the simulation ID to the end of the URL.

   filename: str
       The name of the file to download

   Raises
   ------
   rsp.status_code
       A file-not-found error

   Returns
   -------
   The requested file(s), downloaded and saved to the current directory


