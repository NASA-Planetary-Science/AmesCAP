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
                    help = ("Input sol number. Required. Can either be "
                    "one sol or a range with an increment "
                    "(start stop step)\n"
                    "> Usage: MarsCalendar.py 750.\n"
                    ">        MarsCalendar.py 750 800 5\n\n"))

parser.add_argument('-ls', '--ls', action = 'store_true',
                    help = ("Return the sol number corresponding"
                    "to this Ls\n"
                    "> Usage: MarsCalendar.py start stop step -ls\n\n"))


parser.add_argument('-my', nargs = '+', type = float, default = 0.,
                    help = ("The year of the simulation. MY=0 for "
                    "sol=0-667, MY=1 for sol=668-1335 etc..\n"
                    "> Usage: MarsCalendar.py 350 -ls -my 2\n\n"))

parser.add_argument('-cum', action = 'store_true',
                    help = ("Return Ls in a cummulative form. Example: instead of Ls=0-360, Ls can be 0-720\n"
                    "> Usage: MarsCalendar.py 670 -cum\n\n"))

# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
    # load in user-specified Mars year, if any. Default = 0
    my = np.asarray(parser.parse_args().my).astype(float)
    cum = False
    # set Ls to cumulative, if requested
    if parser.parse_args().cum:
        cum = True

    # load in the user-specified sols to be processed
    # data_input = np.asarray(parser.parse_args().sol).astype(float)
    data_input = parser.parse_args().sol

    if len(data_input) == 1:
        in_array = data_input
        
    elif len(data_input) == 3:
        start, stop, step = data_input[0], data_input[1], data_input[2]
        in_array = np.arange(start, stop, step)
    
    else:
        prRed(f"ERROR a sol number or range (with an increment) is required. See 'MarsCalendar.py -h' for additional help.")
        exit()

    if parser.parse_args().ls:
        # if [-Ls --Ls] is input
        txt_multi='   Ls    |    Sol    '
        result = ls2sol(in_array)
        
    else:
        # if [-Ls --Ls] is excluded
        txt_multi='    SOL    |    Ls    '
        result = sol2ls(in_array, cummulative = cum)

    # if scalar, turn as float
    if len(np.atleast_1d(result)) == 1:
        result = float(result)
    
    #Display data
    print(txt_multi)
    print('-----------------------')
    for i in range(0, len(in_array)):
        print(' %7.2f   |    %7.3f  '%(in_array[i],result[i]+my*668.))

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
