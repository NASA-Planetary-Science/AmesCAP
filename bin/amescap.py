#!/usr/bin/env python3
import argparse
import sys


def create_help_message():
    return """
Welcome to the NASA Ames Community Analysis Pipeline (CAP)!
-----------------------------------------------------------

The Community Analysis Pipeline (CAP) is a Python-based command-line tool 
that performs analysis and creates plots from netCDF files output by the Mars 
Global Climate Model (MGCM). The offical user guide for CAP is available on 
readthedocs at:
https://amescap.readthedocs.io/en/latest/index.html

How do I use CAP?
-----------------
Below is a list of the executables in CAP. Use this list to find the 
executable that performs the operation you desire. To see the arguments for 
each executable, use:
  <command> -h
  Example: MarsVars.py -h 

Then, change to the directory hosting your netCDF output files and pass an 
argument to a CAP executable to perform the operation. A good place to start 
is to use the example command shown below your desired operation.

Available commands:
  MarsCalendar.py  - Converts Ls into day-of-year (sol) and vice versa.
  MarsFiles.py     - Functions for manipulating entire files 
                     (e.g., time-shifting, binning, regridding, temporal and 
                     spatial filtering, tide analysis, zonal averaging, etc.)
  MarsFormat.py    - Transforms non-MGCM model output into MGCM-like model 
                     output for compatibility with CAP.
  MarsInterp.py    - Interpolates files to pressure or altitude coordinates.
  MarsPlot.py      - Generates plots from Custom.in template files.
  MarsPull.py      - Query data from the MGCM repository on the NASA NAS Data 
                     Portal (data.nas.nasa.gov/mcmc).
  MarsVars.py      - Performs variable manipulations (e.g., deriving secondary 
                     variables, column-integration, etc.)

For detailed help on each command, use:
  <command> -h
  Example: MarsVars.py -h

CAP is currently compatible with output from the MCMC’s Legacy and FV3-based 
MGCMs, which are publicly available on GitHub, as well as output from the Mars 
Weather Research and Forecasting Model (MarsWRF) and OpenMars.

CAP is developed and maintained by the Mars Climate Modeling Center (MCMC) at 
NASA’s Ames Research Center in Mountain View, CA. For more information, visit 
our website at:
https://www.nasa.gov/space-science-and-astrobiology-at-ames/division-overview/planetary-systems-branch-overview-stt/mars-climate-modeling-center/
    """
        
def main():
    # Create custom formatter that uses our help message for all help requests
    class CustomFormatter(argparse.HelpFormatter):
        def format_help(self):
            return create_help_message()

    parser = argparse.ArgumentParser(
        description='Welcome to AMESCAP!',
        formatter_class=CustomFormatter
    )
    
    parser.add_argument('command', nargs='?', default='help',
                       help='Command to execute (use "help" for more information)')

    args = parser.parse_args()

    # Print help message for both 'help' command and no command
    if len(sys.argv) == 1 or args.command == 'help':
        print(create_help_message())
        return 0

    return 0

if __name__ == '__main__':
    sys.exit(main())