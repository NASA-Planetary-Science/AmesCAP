#!/usr/bin/env python
"""
The MarsCalendar executable accepts an input Ls or day-of-year (sol)
and returns the corresponding sol or Ls, respectively.

The executable requires 1 of the following arguments:

    * ``[-sol --sol]``          The sol to convert to Ls, OR
    * ``[-ls --ls]``            The Ls to convert to sol

and optionally accepts:

    * ``[-my --marsyear]``      The Mars Year of the simulation to compute sol or Ls from, AND/OR
    * ``[-c --continuous]``     Returns Ls in continuous form

Third-party Requirements:

    * ``numpy``
    * ``argparse``

"""

# Make print statements appear in color
from amescap.Script_utils import (
    Yellow, Nclr, Green, Red
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import numpy as np

# Load amesCAP modules
from amescap.FV3_utils import (sol2ls, ls2sol)

import functools
import traceback

def debug_wrapper(func):
    """
    A decorator that wraps a function with error handling based on the 
    --debug flag.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global debug
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if debug:
                # In debug mode, show the full traceback
                print(f"{Red}ERROR in {func.__name__}: {str(e)}{Nclr}")
                traceback.print_exc()
            else:
                # In normal mode, show a clean error message
                print(f"{Red}ERROR in {func.__name__}: {str(e)}\nUse "
                      f"--debug for more information.{Nclr}")
            return 1  # Error exit code
    return wrapper

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    prog=('MarsCalendar'),
    description=(
        f"{Yellow}Returns the solar longitude (Ls) corresponding to a "
        f"sol or vice-versa. Adapted from areols.py by Tanguy Bertrand."
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

group = parser.add_argument_group(
    "Required Arguments:",
    f"{Yellow}MarsCalendar requires either -ls or -sol.{Nclr}\n")
exclusive_group = parser.add_mutually_exclusive_group(required=True)

exclusive_group.add_argument('-sol', '--sol', nargs='+', type=float,
    help=(
        f"Input sol number. Required. Can either be one sol or a"
        f"range with an increment ``[start stop step]``.\n"
        f"{Green}Example:\n"
        f"> MarsCalendar -sol 10\n"
        f"> MarsCalendar -sol 0 100 20"
        f"{Nclr}\n\n"
    )
)

exclusive_group.add_argument('-ls', '--ls', nargs='+', type=float,
    help=(
        f"Return the sol number corresponding to this Ls.\n"
        f"{Green}Example:\n"
        f"> MarsCalendar -ls 180\n"
        f"> MarsCalendar -ls 90 180 10"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-my', '--marsyear', nargs='+', type=float,
    default = 0.,
    help=(
        f"Return the sol or Ls corresponding to the Ls or sol of a "
        f"particular year of the simulation. \n"
        f"Req. ``[-ls --ls]`` or ``[-sol --sol]``. \n"
        f"``MY=0`` for sol=0-667, ``MY=1`` for sol=668-1335 etc.\n"
        f"{Green}Example:\n"
        f"> Usage: MarsCalendar -ls 350 -my 2"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-c', '--continuous', action='store_true',
    help=(
        f"Return Ls from sol in continuous form. Req. ``[-sol --sol]``."
        f"\nEX: Returns Ls=360-720 instead of Ls=0-360 for "
        f"sol=669-1336\n"
        f"{Green}Example:\n"
        f"> MarsCalendar -sol 700 -c"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsCalendar -sol 700 --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()

if args.sol is None and args.ls is None:
    parser.error(f"{Red}MarsCalendar requires either ``[-sol --sol]`` or "
                 f"``[-ls --ls]``. See ``MarsCalendar -h`` for additional "
                 f"help.{Nclr}")
    exit()

# ======================================================================
#                               DEFINITIONS
# ======================================================================

def parse_array(len_input):
    """
    Formats the input array for conversion.

    Confirms that either ``[-ls --ls]`` or ``[-sol --sol]`` was passed
    as an argument. Creates an array that ls2sol or sol2ls can read
    for the conversion from sol -> Ls or Ls -> sol.

    :param len_input: The input Ls or sol to convert. Can either be
        one number (e.g., 300) or start stop step (e.g., 300 310 2).
    :type len_input: float

    :raises: Error if neither ``[-ls --ls]`` or ``[-sol --sol]`` are
        provided.

    :return: ``input_as_arr`` An array formatted for input into
        ``ls2sol`` or ``sol2ls``. If ``len_input = 300``, then
        ``input_as_arr=[300]``. If ``len_input = 300 310 2``, then
        ``input_as_arr = [300, 302, 304, 306, 308]``.\n

    """
    if len(len_input) == 1:
        input_as_arr = len_input

    elif len(len_input) == 3:
        start, stop, step = len_input[0], len_input[1], len_input[2]
        while (stop < start):
            stop += 360.
        input_as_arr = np.arange(start, stop, step)

    else:
        print(f"{Red}ERROR either ``[-ls --ls]`` or ``[-sol --sol]`` are "
              f"required. See ``MarsCalendar -h`` for additional "
              f"help.{Nclr}")
        exit()
    return(input_as_arr)

# ======================================================================
#                           MAIN PROGRAM
# ======================================================================

@debug_wrapper
def main():
    # Load in user-specified Mars year, if any. Default = 0
    MY = np.squeeze(args.marsyear)
    print(f"MARS YEAR = {MY}")

    if args.continuous:
        # Set Ls to continuous, if requested
        accumulate = True
    else:
        accumulate = False

    if args.ls:
        # If [-Ls --Ls] is input, return sol
        input_num = np.asarray(args.ls).astype(float)
        head_text = "\n   Ls    |    Sol    \n-----------------------"
        input_arr = parse_array(input_num)
        output_arr = ls2sol(input_arr) +MY*668.
        input_arr = input_arr%360

    elif args.sol:
        # If [-sol --sol] is input, return Ls
        input_num = np.asarray(args.sol).astype(float)
        head_text = "\n    SOL  |    Ls    \n-----------------------"
        input_arr = parse_array(input_num)
        output_arr = sol2ls(input_arr, cumulative=accumulate)
        input_arr = input_arr+MY*668.

    # If scalar, return as float
    output_arr = np.atleast_1d(output_arr)

    print(head_text)
    for i in range(0, len(input_arr)):
        # Print input_arr and corresponding output_arr
        print(f" {input_arr[i]:<6.2f}  |  {(output_arr[i]):.2f}")

    print("\n")

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
