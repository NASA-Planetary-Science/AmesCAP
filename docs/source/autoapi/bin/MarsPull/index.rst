:py:mod:`bin.MarsPull`
======================

.. py:module:: bin.MarsPull

.. autoapi-nested-parse::

   The MarsPull executable is for querying data from the Mars Climate
   Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on
   the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

   The executable requires 2 arguments:

       * The directory from which to pull data from, AND
       * ``[-ls --ls]``      The desired solar longitude(s), OR
       * ``[-f --filename]`` The name(s) of the desired file(s)

   Third-party Requirements:

       * ``sys``
       * ``argparse``
       * ``os``
       * ``re``
       * ``numpy``
       * ``functools``
       * ``traceback``
       * ``requests``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsPull.debug_wrapper
   bin.MarsPull.download
   bin.MarsPull.main
   bin.MarsPull.print_file_list



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsPull.Ls_end
   bin.MarsPull.Ls_ini
   bin.MarsPull.args
   bin.MarsPull.debug
   bin.MarsPull.exit_code
   bin.MarsPull.parser
   bin.MarsPull.save_dir


.. py:function:: debug_wrapper(func)

   A decorator that wraps a function with error handling
   based on the --debug flag.
   If the --debug flag is set, it prints the full traceback
   of any exception that occurs. Otherwise, it prints a
   simplified error message.

   :param func: The function to wrap.
   :type  func: function
   :return: The wrapped function.
   :rtype:  function
   :raises Exception: If an error occurs during the function call.
   :raises TypeError: If the function is not callable.
   :raises ValueError: If the function is not found.
   :raises AttributeError: If the function does not have the
       specified attribute.
   :raises IndexError: If the function does not have the
       specified index.


.. py:function:: download(url, file_name)

   Downloads a file from the NAS Data Portal (data.nas.nasa.gov).
   The function takes a URL and a file name as input, and downloads the
   file from the URL, saving it to the specified file name. It also
   provides a progress bar to show the download progress if the file
   size is known. If the file size is unknown, it simply downloads the
   file without showing a progress bar.
   The function handles errors during the download process and prints
   appropriate messages to the console.

   :param url: The url to download from, e.g.,
       'https://data.nas.nasa.gov/legacygcm/fv3betaout1data/03340.fixed.nc'
   :type  url: str
   :param file_name: The local file_name e.g.,
       '/files/Data/LegacyGCM_Ls000_Ls004.nc'
   :type  file_name: str
   :return: The requested file(s), downloaded and saved to the current
       directory.
   :rtype:  None
   :raises FileNotFoundError: A file-not-found error.
   :raises PermissionError: A permission error.
   :raises OSError: An operating system error.
   :raises ValueError: A value error.
   :raises TypeError: A type error.
   :raises requests.exceptions.RequestException: A request error.
   :raises requests.exceptions.HTTPError: An HTTP error.
   :raises requests.exceptions.ConnectionError: A connection error.
   :raises requests.exceptions.Timeout: A timeout error.
   :raises requests.exceptions.TooManyRedirects: A too many redirects
       error.
   :raises requests.exceptions.URLRequired: A URL required error.
   :raises requests.exceptions.InvalidURL: An invalid URL error.
   :raises requests.exceptions.InvalidSchema: An invalid schema error.
   :raises requests.exceptions.MissingSchema: A missing schema error.
   :raises requests.exceptions.InvalidHeader: An invalid header error.
   :raises requests.exceptions.InvalidProxyURL: An invalid proxy URL
       error.
   :raises requests.exceptions.InvalidRequest: An invalid request error.
   :raises requests.exceptions.InvalidResponse: An invalid response
       error.


.. py:function:: main()

   The main function that handles the command-line arguments

   Handles the command-line arguments and coordinates the download
   process. It checks for the presence of the required arguments,
   validates the input, and calls the appropriate functions to download
   the requested files. It also handles the logic for listing available
   directories and files, as well as downloading files based on
   specified solar longitudes (Ls) or file names.

   :return: 0 if successful, 1 if an error occurred.
   :rtype:  int
   :raises SystemExit: If an error occurs during the execution of the
       program, the program will exit with a non-zero status code.


.. py:function:: print_file_list(list_of_files)

   Prints a list of files.

   :param list_of_files: The list of files to print.
   :type  list_of_files: list
   :return: None
   :rtype:  None
   :raises TypeError: If list_of_files is not a list.
   :raises ValueError: If list_of_files is empty.
   :raises IndexError: If list_of_files is out of range.
   :raises KeyError: If list_of_files is not found.
   :raises OSError: If list_of_files is not accessible.
   :raises IOError: If list_of_files is not open.


.. py:data:: Ls_end

   

.. py:data:: Ls_ini

   

.. py:data:: args

   

.. py:data:: debug

   

.. py:data:: exit_code

   

.. py:data:: parser

   

.. py:data:: save_dir

   

