#!/usr/bin/env python3
import argparse
import sys
from amescap.Script_utils import Yellow, Nclr, Green, Cyan

def main():
    parser = argparse.ArgumentParser(
        description='Welcome to AMESCAP!',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS
    )
    
    parser.add_argument('command', nargs='?', default='help',
                       help=argparse.SUPPRESS)

    args = parser.parse_args()

    help_message = f"""
{Green}Welcome to the NASA Ames Community Analysis Pipeline (CAP)!
-----------------------------------------------------------{Nclr}

The Community Analysis Pipeline (CAP) is a Python-based command-line tool that performs analysis and creates plots from netCDF files output by the Mars Global Climate Model (MGCM). The offical user guide for CAP is available on readthedocs at:
{Yellow}https://amescap.readthedocs.io/en/latest/index.html{Nclr}

{Cyan}How do I use CAP?
-----------------
Below is a list of the executables in CAP. Use this list to find the executable that performs the operation you desire. To see the arguments for each executable, use:
  {Green}<command> -h
  Example: MarsVars.py -h{Cyan}

Then, change to the directory hosting your netCDF output files and pass an argument to a CAP executable to perform the operation. A good place to start is to use the example command shown below your desired operation.{Nclr}

Available commands:
  {Green}MarsCalendar.py{Nclr}  - Converts Ls into day-of-year (sol) and vice versa.
  {Green}MarsFiles.py{Nclr}     - Functions for manipulating entire files (e.g., time-shifting, binning, regridding, temporal and spatial filtering, tide analysis, zonal averaging, etc.)
  {Green}MarsFormat.py{Nclr}    - Transforms non-MGCM model output into MGCM-like model output for compatibility with CAP.
  {Green}MarsInterp.py{Nclr}    - Interpolates files to pressure or altitude coordinates.
  {Green}MarsPlot.py{Nclr}      - Generates plots from Custom.in template files.
  {Green}MarsPull.py{Nclr}      - Query data from the MGCM repository on the NASA NAS Data Portal (data.nas.nasa.gov/mcmc).
  {Green}MarsVars.py{Nclr}      - Performs variable manipulations (e.g., deriving secondary variables, column-integration, etc.)

{Green}For detailed help on each command, use:
  <command> -h
  Example: MarsVars.py -h{Nclr}

CAP is currently compatible with output from the MCMC’s Legacy and FV3-based MGCMs, which are publicly available on GitHub, as well as output from the Mars Weather Research and Forecasting Model (MarsWRF) and OpenMars.

CAP is developed and maintained by the Mars Climate Modeling Center (MCMC) at NASA’s Ames Research Center in Mountain View, CA. For more information, visit our website at:
https://www.nasa.gov/space-science-and-astrobiology-at-ames/division-overview/planetary-systems-branch-overview-stt/mars-climate-modeling-center/
    """
    print(help_message)
    return 0

if __name__ == '__main__':
    sys.exit(main())