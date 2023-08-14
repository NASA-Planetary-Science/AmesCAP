:py:mod:`MarsCalendar`
======================

.. py:module:: MarsCalendar

.. autoapi-nested-parse::

   The MarsCalendar executable accepts an input Ls or day-of-year (sol) and
   returns the corresponding sol or Ls, respectively.

   The executable requires 1 of the following arguments:
       * [-sol --sol]          The sol to convert to Ls, OR
       * [-ls --ls]            the Ls to convert to sol.

   and optionally accepts 2 arguments:
       * [-my --marsyear]      The Mars Year of the simulation to compute 
                               sol or Ls from, AND/OR
       * [-c --cumulative]     Returns Ls in cumulative form.

   Third-party Requirements:
       * numpy
       * argparse

   List of Functions:
       * parse_array - Formats the requested sol/Ls for conversion to
                       Ls/sol. Computes arrays from [start, stop, step] if
                       necessary.



