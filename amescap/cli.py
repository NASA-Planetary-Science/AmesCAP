#!/usr/bin/env python3
import argparse
import sys
import os
import time
from amescap.Script_utils import Yellow, Nclr, Green, Cyan

def get_install_info():
    # Get the location and timestamp of cli.py 
    cli_path = os.path.abspath(__file__)
    install_time = time.ctime(os.path.getctime(cli_path))
    return f"""
{Cyan}CAP Installation Information
----------------------------{Nclr}
{Cyan}Version:{Nclr} 0.3
{Cyan}Install Date:{Nclr} {install_time}
{Cyan}Install Location:{Nclr} {os.path.dirname(os.path.dirname(cli_path))}
"""

def main():
    # Custom help formatter to override the default help message
    class CustomHelpFormatter(argparse.HelpFormatter):
        def format_help(self):
            help_message = f"""
            
{Cyan}
Welcome to the NASA Ames Community Analysis Pipeline (CAP)!
-----------------------------------------------------------{Nclr}
The Community Analysis Pipeline (CAP) is a Python-based command-line tool that performs analysis and creates plots from netCDF files output by the Mars Global Climate Model (MGCM). The offical user guide for CAP is available on readthedocs:
{Yellow}https://amescap.readthedocs.io/en/latest/index.html{Nclr}

{Cyan}How do I use CAP?
-----------------{Nclr}
Below is a list of the executables in CAP. Use this list to find the executable that performs the operation you desire. 
To see the arguments for each executable, use:
{Green}<command> -h
 Example: MarsVars.py -h{Nclr}

Then, change to the directory hosting your netCDF output files and pass an argument to a CAP executable to perform the operation. A good place to start is to use the example command shown below your desired operation.{Nclr}

{Cyan}Available Commands
------------------{Nclr}
{Green}MarsCalendar.py {Nclr}- Converts Ls into day-of-year (sol) and vice versa.
{Green}MarsFiles.py    {Nclr}- Manipulates entire files (e.g., time-shift, regrid, filter, etc.)
{Green}MarsFormat.py   {Nclr}- Transforms non-MGCM model output for compatibility with CAP.
{Green}MarsInterp.py   {Nclr}- Interpolates files to pressure or altitude coordinates.
{Green}MarsPlot.py     {Nclr}- Generates plots from Custom.in template files.
{Green}MarsPull.py     {Nclr}- Queries data from the MGCM repository at data.nas.nasa.gov/mcmc
{Green}MarsVars.py     {Nclr}- Performs variable manipulations (e.g., deriving secondary variables, column-integration, etc.)

{Cyan}Model Compatibility
-------------------{Nclr}
CAP is currently compatible with output from the:
- MCMC Legacy MGCM
- MCMC FV3-based MGCM
- Mars Weather Research and Forecasting Model (MarsWRF)
- OpenMars

{Cyan}Additional Information
----------------------{Nclr}
CAP is developed and maintained by the **Mars Climate Modeling Center (MCMC) at NASA's Ames Research Center** in Mountain View, CA. For more information, visit our website:
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

    # Check for version/info command before printing help
    if args.command in ['version', 'info']:
        print(get_install_info())
        return 0
    
    # Print help for all cases
    parser.print_help()
    return 0

if __name__ == '__main__':
    sys.exit(main())