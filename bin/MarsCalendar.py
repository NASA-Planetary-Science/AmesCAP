#!/usr/bin/env python
"""
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
"""

# make print statements appear in color
from amescap.Script_utils import prRed
Cyan = "\033[96m"
Yellow = "\033[93m"
Default = "\033[00m"

# load generic Python modules
import argparse     # parse arguments
import numpy as np

# load amesCAP modules
from amescap.FV3_utils import sol2ls, ls2sol

# ======================================================
#                  ARGUMENT PARSER
# ======================================================
parser = argparse.ArgumentParser(
    description=(f"{Yellow}Returns the solar longitude (Ls) "
                 f"corresponding to a sol or vice-versa. Adapted "
                 f"from areols.py {Default}"),
    formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument('-sol', '--sol', nargs = '+', type = float,
                    help = (f"Input sol number. Required. Can either "
                    f"be one sol or a range with an increment "
                    f"(start stop step).\n"
                    f"> Usage: MarsCalendar.py -sol 750\n"
                    f">        MarsCalendar.py -sol 750 800 5\n\n"))

parser.add_argument('-ls', '--ls', nargs = '+', type = float,
                    help = (f"Return the sol number corresponding "
                    f"to this Ls.\n"
                    f"> Usage: MarsCalendar.py -ls 350\n"
                    f">        MarsCalendar.py -ls 350 360 5\n\n"))


parser.add_argument('-my', '--marsyear', nargs = '+', type = float, 
                    default = 0.,
                    help = (f"Return the sol or Ls corresponding "
                    f"to the Ls or sol of a particular year of the "
                    f"simulation. Req. [-ls --ls] or [-sol --sol]. \n"
                    f"MY=0 for sol=0-667, MY=1 for sol=668-1335 etc.\n"
                    f"> Usage: MarsCalendar.py -ls 350 -my 2\n\n"))

parser.add_argument('-c', '--cumulative', action='store_true',
                    help = (f"Return Ls from sol in cumulative form. "
                    f"Req. [-sol --sol]. \n"
                    f"EX: Returns Ls=360-720 instead of Ls=0-360 "
                    f"for input sol=669-1336 \n"
                    f"> Usage: MarsCalendar.py -sol 700 -c\n\n"))

# ======================================================
#                  DEFINITIONS
# ======================================================

def parse_array(numIN):
    """
    Formats the input array for conversion.
    
    Confirms that either [-ls --ls] or [-sol --sol] was passed as an
    argument. Creates an array that ls2sol or sol2ls can read for the
    conversion from sol -> Ls or Ls -> sol.
    
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
    """
    
    if len(numIN) == 1:
        arrayIN = numIN

    elif len(numIN) == 3:
        start, stop, step = numIN[0], numIN[1], numIN[2]
        arrayIN = np.arange(start, stop, step)

    else:
        prRed("ERROR either [-ls --ls] or [-sol --sol] are required. "
                "See 'MarsCalendar.py -h' for additional help.")
        exit()
    return(arrayIN)

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
        inputNum = np.asarray(parser.parse_args().ls).astype(float)
        headerText = '\n   Ls    |    Sol    \n-----------------------'
        inputArr = parse_array(inputNum)
        outputArr = ls2sol(inputArr)
        
    elif parser.parse_args().sol:
        # if [-sol --sol] is input, return Ls
        inputNum = np.asarray(parser.parse_args().sol).astype(float)
        headerText = '\n    SOL  |    Ls    \n-----------------------'
        inputArr = parse_array(inputNum)
        outputArr = sol2ls(inputArr, cumulative=accumulate)

    # if scalar, return as float
    outputArr = np.atleast_1d(outputArr)
    
    print(headerText)
    for i in range(0, len(inputArr)):
        # print inputArr and corresponding outputArr
        print(f" {inputArr[i]:.2f}  |  {(outputArr[i]+MY*668.):.2f}  ")
    
    print(" ")

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
