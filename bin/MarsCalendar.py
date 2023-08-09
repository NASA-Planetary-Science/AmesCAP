#!/usr/bin/env python
"""
The MarsCalendar executable is for ...

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
                    help = (f"Input sol number. Required. Can either \n"
                    f"be one sol or a range with an increment \n"
                    f"(start stop step)\n"
                    f"> Usage: MarsCalendar.py -sol 750\n"
                    f">        MarsCalendar.py -sol 750 800 5\n\n"))

parser.add_argument('-ls', '--ls', nargs = '+', type = float,
                    help = (f"Return the sol number corresponding \n"
                    f"to this Ls\n"
                    f"> Usage: MarsCalendar.py -ls 350\n"
                    f">        MarsCalendar.py -ls 350 360 5\n\n"))


parser.add_argument('-my', '--marsyear', nargs = '+', type = float, 
                    default = 0.,
                    help = (f"The year of the simulation. \n"
                    f"MY=0 for sol=0-667, MY=1 for sol=668-1335 etc.\n"
                    f"> Usage: MarsCalendar.py -ls 350 -my 2\n\n"))

parser.add_argument('-cum', '--cum', action = 'store_true',
                    help = (f"Return Ls in a cummulative form. \n"
                    f"Example: instead of Ls=0-360, Ls can be 0-720\n"
                    f"> Usage: MarsCalendar.py -sol 670 -cum\n\n"))

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
    MY = np.asarray(parser.parse_args().marsyear).astype(float)
    
    if parser.parse_args().cum:
        # set Ls to cumulative, if requested
        cum = True
    else:
        cum = False

    if parser.parse_args().ls:
        # if [-Ls --Ls] is input, return sol
        numIN = np.asarray(parser.parse_args().ls).astype(float)
        txt_multi='\n   Ls    |    Sol    '
        arrayIN = parse_array(numIN)
        arrayOUT = ls2sol(arrayIN)
        
    elif parser.parse_args().sol:
        # if [-sol --sol] is input, return Ls
        numIN = np.asarray(parser.parse_args().sol).astype(float)
        txt_multi='\n    SOL  |    Ls    '
        arrayIN = parse_array(numIN)
        arrayOUT = sol2ls(arrayIN, cummulative = cum)

    # if scalar, turn as float
    if len(np.atleast_1d(arrayOUT)) == 1:
        arrayOUT = [arrayOUT]
    
    # print arrayIN and corresponding arrayOUT
    print(txt_multi)
    print('-----------------------')
    for i in range(0, len(arrayIN)):
        print(f" {arrayIN[i]:.2f}  |  {(arrayOUT[i]+MY*668.):.2f}  ")
    
    print(" ")

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
