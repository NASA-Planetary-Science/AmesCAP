:py:mod:`MarsCalendar`
======================

.. py:module:: MarsCalendar

.. autoapi-nested-parse::

   The MarsCalendar executable accepts an input Ls or day-of-year (sol) and returns the corresponding sol or Ls, respectively.

   The executable requires 1 of the following arguments:
       * [-sol --sol]          The sol to convert to Ls, OR
       * [-ls --ls]            the Ls to convert to sol.

   and optionally accepts 2 arguments:
       * [-my --marsyear]      The Mars Year of the simulation to compute                            sol or Ls from, AND/OR
       * [-c --cumulative]     Returns Ls in cumulative form.

   Third-party Requirements:
       * numpy
       * argparse

   List of Functions:
       * parse_array - Formats the requested sol/Ls for conversion to                     Ls/sol. Computes arrays from [start, stop, step]                     if necessary.



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   MarsCalendar.parse_array



.. py:function:: parse_array(numIN)

   Formats the input array for conversion.

   Confirms that either [-ls --ls] or [-sol --sol] was passed as an     argument. Creates an array that ls2sol or sol2ls can read for the     conversion from sol -> Ls or Ls -> sol.

   Parameters
   ----------
   numIN : float
       The input Ls or sol to convert. Can either be one number:

           300

       or start stop step:

           300 310 2

   Raises
   ------
   Error if neither [-ls --ls] or [-sol --sol] is provided.

   Returns
   -------
   arrayIN
       An array formatted for input into ls2sol or sol2ls. 

           If numIN = 300, arrayIN=[300]

           If numIN = 300 310 2, arrayIN = [300, 302, 304, 306, 308]    


