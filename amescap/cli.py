#!/usr/bin/env python3
import argparse
import sys
from amescap.Script_utils import Yellow, Nclr, Green, Cyan

def main():
    # Custom help formatter to override the default help message
    class CustomHelpFormatter(argparse.HelpFormatter):
        def format_help(self):
            help_message = f"""
{Cyan}
Welcome to the NASA Ames Community Analysis Pipeline (CAP)!
-----------------------------------------------------------{Nclr}
The Community Analysis Pipeline (CAP) is a Python-based command-line tool that performs analysis and creates plots from netCDF files output by the Mars Global Climate Model (MGCM). The offical user guide for CAP is available on readthedocs at:
{Yellow}https://amescap.readthedocs.io/en/latest/index.html{Nclr}

{Cyan}How do I use CAP?
-----------------{Nclr}
Below is a list of the executables in CAP. Use this list to find the executable that performs the operation you desire. To see the arguments for each executable, use:
{Green}<command> -h
 Example: MarsVars.py -h{Nclr}

Then, change to the directory hosting your netCDF output files and pass an argument to a CAP executable to perform the operation. A good place to start is to use the example command shown below your desired operation.{Nclr}

{Cyan}Available Commands
------------------{Nclr}
{Yellow}MarsCalendar.py - Converts Ls into day-of-year (sol) and vice versa.
{Cyan}MarsFiles.py    - Manipulates entire files (e.g., time-shift, regrid, filter, etc.)
{Yellow}MarsFormat.py   - Transforms non-MGCM model output for compatibility with CAP.
{Cyan}MarsInterp.py   - Interpolates files to pressure or altitude coordinates.
{Yellow}MarsPlot.py     - Generates plots from Custom.in template files.
{Cyan}MarsPull.py     - Queries data from the MGCM repository at data.nas.nasa.gov/mcmc
{Yellow}MarsVars.py     - Performs variable manipulations (e.g., deriving secondary variables, column-integration, etc.)

{Green}For detailed help on each command, use:
 <command> -h
 Example: MarsVars.py -h{Nclr}

{Cyan}Additional Information
----------------------{Nclr}
CAP is currently compatible with output from:
- The MCMC Legacy MGCM
- The MCMC FV3-based MGCM (which are publicly available on GitHub)
- Mars Weather Research and Forecasting Model (MarsWRF)
- OpenMars

CAP is developed and maintained by the **Mars Climate Modeling Center (MCMC)** at NASA's Ames Research Center in Mountain View, CA. For more information, visit our website at:
{Yellow}https://www.nasa.gov/space-science-and-astrobiology-at-ames/division-overview/planetary-systems-branch-overview-stt/mars-climate-modeling-center/
{Nclr}
"""
            return help_message

    parser = argparse.ArgumentParser(
        formatter_class=CustomHelpFormatter,
        add_help=False  # This prevents the default help message
    )
    
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                       help=argparse.SUPPRESS)
    parser.add_argument('command', nargs='?', default='help',
                       help=argparse.SUPPRESS)

    args = parser.parse_args()

    # Print help for all cases
    parser.print_help()
    return 0

if __name__ == '__main__':
    sys.exit(main())