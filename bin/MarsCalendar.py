#!/usr/bin/env python
"""
The MarsCalendar executable accepts an input Ls or day-of-year (sol) \
and returns the corresponding sol or Ls, respectively.

The executable requires 1 of the following arguments:
    * [-sol --sol]          The sol to convert to Ls, OR
    * [-ls --ls]            the Ls to convert to sol.

and optionally accepts 2 arguments:
    * [-my --marsyear]      The Mars Year of the simulation to compute\
                            sol or Ls from, AND/OR
    * [-c --cumulative]     Returns Ls in cumulative form.

Third-party Requirements:
    * numpy
    * argparse

List of Functions:
    * parse_array - Formats the requested sol/Ls for conversion to \
                    Ls/sol. Computes arrays from [start, stop, step] \
                    if necessary.
"""

# make print statements appear in color
from amescap.Script_utils import prRed, Yellow, NoColor, Green

# load generic Python modules
import argparse     # parse arguments
import numpy as np

# load amesCAP modules
from amescap.FV3_utils import sol2ls, ls2sol

# ======================================================
#                  ARGUMENT PARSER
# ======================================================
parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}Returns the solar longitude (Ls) corresponding to a "
        f"sol or vice-versa. Adapted from areols.py.{NoColor}\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    '-sol', '--sol', nargs='+', type=float,
    help=(
        f"Input sol number. Required. Can either be one sol or a"
        f"range with an increment (start stop step).\n"
        f"{Green}Usage:\n"
        f"> MarsCalendar.py -sol 750\n"
        f"> MarsCalendar.py -sol 750 800 5"
        f"{NoColor}\n\n"
    )
)

parser.add_argument(
    '-ls', '--ls', nargs='+', type=float,
    help=(
        f"Return the sol number corresponding to this Ls.\n"
        f"{Green}Usage:\n"
        f"> MarsCalendar.py -ls 350\n"
        f"> MarsCalendar.py -ls 350 360 5"
        f"{NoColor}\n\n"
    )
)

parser.add_argument(
    '-my', '--marsyear', nargs='+', type=float, default = 0.,
    help=(
        f"Return the sol or Ls corresponding to the Ls or sol of a "
        f"particular year of the simulation. \n"
        f"Req. [-ls --ls] or [-sol --sol]. \n"
        f"MY=0 for sol=0-667, MY=1 for sol=668-1335 etc.\n"
        f"{Green}Usage:\n"
        f"> Usage: MarsCalendar.py -ls 350 -my 2"
        f"{NoColor}\n\n"
    )
)

parser.add_argument(
    '-c', '--cumulative', action='store_true',
    help=(
        f"Return Ls from sol in cumulative form. Req. [-sol --sol].\n"
        f"EX: Returns Ls=360-720 instead of Ls=0-360 for input "
        f"sol=669-1336 \n"
        f"{Green}Usage:\n"
        f"> MarsCalendar.py -sol 700 -c"
        f"{NoColor}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Debug flag: release the exceptions.\n\n"
    )
)

# ======================================================
#                  DEFINITIONS
# ======================================================

def parse_array(len_input):
    """
    Formats the input array for conversion.

    Confirms that either [-ls --ls] or [-sol --sol] was passed as an \
    argument. Creates an array that ls2sol or sol2ls can read for the \
    conversion from sol -> Ls or Ls -> sol.

    Parameters
    ----------
    len_input : float
        The input Ls or sol to convert. Can either be one number:\n
            300\n
        or start stop step:\n
            300 310 2

    Raises
    ------
    Error if neither [-ls --ls] or [-sol --sol] is provided.

    Returns
    -------
    input_as_arr
        An array formatted for input into ls2sol or sol2ls.\n
            If len_input = 300:\n
                input_as_arr=[300]\n
            If len_input = 300 310 2:\n
                input_as_arr = [300, 302, 304, 306, 308]\n
    """

    if len(len_input) == 1:
        input_as_arr = len_input

    elif len(len_input) == 3:
        start, stop, step = len_input[0], len_input[1], len_input[2]
        input_as_arr = np.arange(start, stop, step)

    else:
        prRed("ERROR either [-ls --ls] or [-sol --sol] are required. "
              "See 'MarsCalendar.py -h' for additional help.")
        exit()
    return(input_as_arr)

# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
    # load in user-specified Mars year, if any. Default = 0
    MY = np.squeeze(parser.parse_args().marsyear)
    print(f"MARS YEAR = {MY}")

    if parser.parse_args().cumulative:
        # set Ls to cumulative, if requested
        accumulate = True
    else:
        accumulate = False

    if parser.parse_args().ls:
        # if [-Ls --Ls] is input, return sol
        input_num = np.asarray(parser.parse_args().ls).astype(float)
        head_text = '\n   Ls    |    Sol    \n-----------------------'
        input_arr = parse_array(input_num)
        output_arr = ls2sol(input_arr)

    elif parser.parse_args().sol:
        # if [-sol --sol] is input, return Ls
        input_num = np.asarray(parser.parse_args().sol).astype(float)
        head_text = '\n    SOL  |    Ls    \n-----------------------'
        input_arr = parse_array(input_num)
        output_arr = sol2ls(input_arr, cumulative=accumulate)

    # if scalar, return as float
    output_arr = np.atleast_1d(output_arr)

    print(head_text)
    for i in range(0, len(input_arr)):
        # print input_arr and corresponding output_arr
        print(f" {input_arr[i]:.2f}  |  {(output_arr[i]+MY*668.):.2f}")

    print("\n")

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
