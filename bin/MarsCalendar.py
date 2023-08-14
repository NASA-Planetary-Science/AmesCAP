#!/usr/bin/env python
"""
The MarsCalendar executable accepts an input Ls or day-of-year (sol) and
returns the corresponding sol or Ls, respectively.

The executable requires x arguments:
    * [-x --x]      define

Third-party Requirements:
    * numpy
    * argparse

List of Functions:
    * x
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
                    help = (f"Return the sol number corresponding "
                    f"to this Ls of a particular year of the "
                    f"simulation. Requires [-ls --ls]. \n"
                    f"MY=0 for sol=0-667, MY=1 for sol=668-1335 etc.\n"
                    f"> Usage: MarsCalendar.py -ls 350 -my 2\n\n"))

parser.add_argument('-c', '--cumulative', action='store_true',
                    help = (f"Return Ls from sol in cumulative form. "
                    f"Requires [-sol --sol]. \n"
                    f"EX: Returns Ls=360-720 instead of Ls=0-360 "
                    f"for input sol=669-1336 \n"
                    f"> Usage: MarsCalendar.py -sol 700 -c\n\n"))

# ======================================================
#                  DEFINITIONS
# ======================================================

def parse_array(numIN):
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
    # MY = np.asarray(parser.parse_args().marsyear).astype(float)
    MY = parser.parse_args().marsyear
    # If MY is a float(e.g. psfc=700.), make it a 1-element array(e.g. psfc=[700])
    MY = np.squeeze(MY)
    print(f"MARS YEAR = {MY}")
    
    if parser.parse_args().cumulative:
        # set Ls to cumulative, if requested
        accumulate = True
    else:
        accumulate = False

    if parser.parse_args().ls:
        # if [-Ls --Ls] is input, return sol
        numIN = np.asarray(parser.parse_args().ls).astype(float)
        headerText = '\n   Ls    |    Sol    \n-----------------------'
        arrayIN = parse_array(numIN)
        arrayOUT = ls2sol(arrayIN)
        
    elif parser.parse_args().sol:
        # if [-sol --sol] is input, return Ls
        numIN = np.asarray(parser.parse_args().sol).astype(float)
        headerText = '\n    SOL  |    Ls    \n-----------------------'
        arrayIN = parse_array(numIN)
        arrayOUT = sol2ls(arrayIN, cumulative=accumulate)

    # if scalar, return as float
    arrayOUT = np.atleast_1d(arrayOUT)
    
    print(headerText)
    for i in range(0, len(arrayIN)):
        # print arrayIN and corresponding arrayOUT
        print(f" {arrayIN[i]:.2f}  |  {(arrayOUT[i]+MY*668.):.2f}  ")
    
    print(" ")

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
