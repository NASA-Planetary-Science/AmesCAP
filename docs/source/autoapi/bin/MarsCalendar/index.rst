:py:mod:`bin.MarsCalendar`
==========================

.. py:module:: bin.MarsCalendar

.. autoapi-nested-parse::

   The MarsCalendar executable accepts an input Ls or day-of-year (sol)
   and returns the corresponding sol or Ls, respectively.

   The executable requires 1 of the following arguments:

       * ``[-sol --sol]``      The sol to convert to Ls, OR
       * ``[-ls --ls]``        The Ls to convert to sol

   and optionally accepts:

       * ``[-my --marsyear]``  The Mars Year of the simulation to compute sol or Ls from, AND/OR
       * ``[-c --continuous]`` Returns Ls in continuous form

   Third-party Requirements:

       * ``numpy``
       * ``argparse``
       * ``functools``
       * ``traceback``
       * ``sys``
       * ``amescap``



Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   bin.MarsCalendar.debug_wrapper
   bin.MarsCalendar.main
   bin.MarsCalendar.parse_array



Attributes
~~~~~~~~~~

.. autoapisummary::

   bin.MarsCalendar.args
   bin.MarsCalendar.debug
   bin.MarsCalendar.exclusive_group
   bin.MarsCalendar.exit_code
   bin.MarsCalendar.group
   bin.MarsCalendar.parser


.. py:function:: debug_wrapper(func)

   A decorator that wraps a function with error handling
   based on the --debug flag.

   If the --debug flag is set, it prints the full traceback
   of any exception that occurs. Otherwise, it prints a
   simplified error message.

   :param func: The function to wrap.
   :type   func: function
   :return: The wrapped function.
   :rtype:  function
   :raises Exception: If an error occurs during the function call.
   :raises TypeError: If the function is not callable.
   :raises ValueError: If the function is not found.
   :raises NameError: If the function is not defined.
   :raises AttributeError: If the function does not have the
       specified attribute.
   :raises ImportError: If the function cannot be imported.
   :raises RuntimeError: If the function cannot be run.
   :raises KeyError: If the function does not have the
       specified key.
   :raises IndexError: If the function does not have the
       specified index.
   :raises IOError: If the function cannot be opened.
   :raises OSError: If the function cannot be accessed.
   :raises EOFError: If the function cannot be read.
   :raises MemoryError: If the function cannot be allocated.
   :raises OverflowError: If the function cannot be overflowed.
   :raises ZeroDivisionError: If the function cannot be divided by zero.
   :raises StopIteration: If the function cannot be stopped.
   :raises KeyboardInterrupt: If the function cannot be interrupted.
   :raises SystemExit: If the function cannot be exited.
   :raises AssertionError: If the function cannot be asserted.


.. py:function:: main()

   Main function for MarsCalendar command-line tool.

   This function processes user-specified arguments to convert between 
   Mars solar longitude (Ls) and sol (Martian day) values for a given 
   Mars year. It supports both continuous and discrete sol 
   calculations, and can handle input in either Ls or sol, returning 
   the corresponding converted values. The results are printed in a 
   formatted table.

   Arguments are expected to be provided via the `args` namespace:

       - args.marsyear: Mars year (default is 0)
       - args.continuous: If set, enables continuous sol calculation
       - args.ls: List of Ls values to convert to sol
       - args.sol: List of sol values to convert to Ls
       
   :param args: Command-line arguments parsed using argparse.
   :type  args: argparse.Namespace
   :raises ValueError: If the input is not a valid number or
       range.
   :returns: 0 if the program runs successfully, 1 if an error occurs.
   :rtype:  int
   :raises RuntimeError: If the input is not a valid runtime.
   :raises TypeError: If the input is not a valid type.
   :raises IndexError: If the input is not a valid index.
   :raises KeyError: If the input is not a valid key.
   :raises AttributeError: If the input is not a valid attribute.
   :raises ImportError: If the input is not a valid import.


.. py:function:: parse_array(len_input)

   Formats the input array for conversion.

   Confirms that either ``[-ls --ls]`` or ``[-sol --sol]`` was passed
   as an argument. Creates an array that ls2sol or sol2ls can read
   for the conversion from sol -> Ls or Ls -> sol.

   :param len_input: The input Ls or sol to convert. Can either be
       one number (e.g., 300) or start stop step (e.g., 300 310 2).
   :type  len_input: float or list of floats
   :raises: Error if neither ``[-ls --ls]`` or ``[-sol --sol]`` are
       provided.
   :return: ``input_as_arr`` An array formatted for input into
       ``ls2sol`` or ``sol2ls``. If ``len_input = 300``, then
       ``input_as_arr=[300]``. If ``len_input = 300 310 2``, then
       ``input_as_arr = [300, 302, 304, 306, 308]``.
   :rtype:  list of floats
   :raises ValueError: If the input is not a valid number or
       range.
   :raises TypeError: If the input is not a valid type.
   :raises IndexError: If the input is not a valid index.
   :raises KeyError: If the input is not a valid key.
   :raises AttributeError: If the input is not a valid attribute.
   :raises ImportError: If the input is not a valid import.
   :raises RuntimeError: If the input is not a valid runtime.
   :raises OverflowError: If the input is not a valid overflow.
   :raises MemoryError: If the input is not a valid memory.


.. py:data:: args

   

.. py:data:: debug

   

.. py:data:: exclusive_group

   

.. py:data:: exit_code

   

.. py:data:: group

   

.. py:data:: parser

   

