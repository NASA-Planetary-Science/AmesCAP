#!/usr/bin/env python3
"""
The MarsPull executable is for querying data from the Mars Climate \
Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on \
the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

The executable requires 2 arguments:

    * The directory from which to pull data from, AND
    * ``[-ls --ls]``      The desired solar longitude(s), OR
    * ``[-f --filename]`` The name(s) of the desired file(s)

Third-party Requirements:

    * ``sys``
    * ``argparse``
    * ``os``
    * ``re``
    * ``numpy``
    * ``functools``
    * ``traceback``
    * ``requests``

List of Functions:

    * download - Queries the requested file from the NAS Data Portal.
"""

# make print statements appear in color
from amescap.Script_utils import (
    prYellow, prCyan, Green, Yellow, Nclr, Cyan, Blue, Red
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import re           # Regular expressions
import numpy as np
import functools    # For function decorators
import traceback    # For printing stack traces
import requests     # Download data from website


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


# ------------------------------------------------------
#                  ARGUMENT PARSER
# ------------------------------------------------------

parser = argparse.ArgumentParser(
    prog=('MarsPull'),
    description=(
        f"{Yellow}Uility for downloading NASA Ames Mars Global Climate "
        f"Model output files from the NAS Data Portal at:\n"
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/\n{Nclr}"
        f"Requires ``-f`` or ``-ls``."
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('directory_name', type=str, nargs='?',
    choices=[
        'FV3BETAOUT1', 'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS',
        'ACTIVECLDS_NCDF'],
    help=(
        f"Selects the simulation directory from the "
        f"NAS data portal ("
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/){Nclr}\n"
        f"Current directory options are:\n{Yellow}FV3BETAOUT1, ACTIVECLDS, "
        f"ACTIVECLDS, INERTCLDS, NEWBASE_ACTIVECLDS, ACTIVECLDS_NCDF\n"
        f"{Red}MUST be used with either ``-f`` or ``-ls``\n"
        f"{Green}Example:\n"
        f"> MarsPull INERTCLDS -f fort.11_0690\n"
        f"{Blue}OR{Green}\n"
        f"> MarsPull INERTCLDS -ls 90\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-list', '--list_files', action='store_true',
    help=(
        f"Return a list of all the files available for download from "
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/{Nclr}\n"
        f"{Green}Example:\n"
        f"> MarsPull -list"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-f', '--filename', nargs='+', type=str,
    help=(
        f"The name(s) of the file(s) to download.\n"
        f"{Green}Example:\n"
        f"> MarsPull INERTCLDS -f fort.11_0690"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ls', '--ls', nargs='+', type=float,
    help=(
        f"Selects the file(s) to download based on a range of solar "
        f"longitudes (Ls).\n"
        f"This only works on data in the {Yellow}ACTIVECLDS{Nclr} and "
        f"{Yellow}INERTCLDS{Nclr} folders.\n"
        f"{Green}Example:\n"
        f"> MarsPull INERTCLDS -ls 90\n"
        f"> MarsPull INERTCLDS -ls 90 180"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsPull INERTCLDS -ls 90 --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()
debug = args.debug


# ------------------------------------------------------
#                  DEFINITIONS
# ------------------------------------------------------
save_dir = (f"{os.getcwd()}/")

# available files by Ls:
Ls_ini = np.array([
    0, 5, 10, 15, 19, 24, 29, 34, 38, 43, 48, 52, 57, 61, 66, 70, 75,
    79, 84, 88, 93, 97, 102, 106, 111, 116, 121, 125, 130, 135, 140,
    146, 151, 156, 162, 167, 173, 179, 184, 190, 196, 202, 209, 215,
    221, 228, 234, 241, 247, 254, 260, 266, 273, 279, 286, 292, 298,
    304, 310, 316,322, 328, 333, 339, 344, 350, 355
])

Ls_end = np.array([
    4, 9, 14, 19, 24, 29, 33, 38, 42, 47, 52, 56, 61, 65, 70, 74, 79,
    83, 88, 92, 97, 101, 106, 111, 115, 120, 125, 130, 135, 140, 145,
    150, 156, 161, 167, 172, 178, 184, 190, 196, 202, 208, 214, 221,
    227, 233, 240, 246, 253, 259, 266, 272, 279, 285, 291, 297, 304,
    310, 316, 321, 327, 333, 338, 344, 349, 354, 0
])


def download(url, file_name):
    """
    Downloads a file from the NAS Data Portal (data.nas.nasa.gov).

    :param url: The url to download from, e.g 'https://data.nas.nasa.\
        gov/legacygcm/download_data.php?file=/legacygcmdata/LegacyGCM_\
            Ls000_Ls004.nc'
    :type url: str
    :param file_name: The local file_name e.g  '/lou/la4/akling/Data/L\
        egacyGCM_Ls000_Ls004.nc'
    :type file_name: str
    :return: The requested file(s), downloaded and saved to the \
        current directory.
    :raises FileNotFoundError: A file-not-found error.
    """

    _, fname = os.path.split(file_name)
    response = requests.get(url, stream=True)
    total = response.headers.get('content-length')

    if response.status_code == 404:
        print(f'Error during download, error code is: {response.status_code}')
    else:
        if total is not None:
            # If file size is known, return progress bar
            with open(file_name, 'wb') as f:
                downloaded = 0
                if total:
                    total = int(total)
                for data in response.iter_content(
                    chunk_size = max(int(total/1000), 1024*1024)
                    ):
                    downloaded += len(data)
                    f.write(data)
                    status = int(50*downloaded/total)
                    sys.stdout.write(f"\r[{'#'*status}{'.'*(50 - status)}]")
                    sys.stdout.flush()
            sys.stdout.write('\n')
        else:
            # If file size is unknown, skip progress bar
            print(f"Downloading {fname}...")
            with open(file_name, 'wb')as f:
                f.write(response.content)
            print(f"{fname} Done")


def print_file_list(list_of_files):
    for file in list_of_files:
        print(file)


# ------------------------------------------------------
#                  MAIN FUNCTION
# ------------------------------------------------------

@debug_wrapper
def main():
    # validation
    if not args.list_files and not args.directory_name:
        print("Error: You must specify either -list or a directory.")
        sys.exit(1)

    if args.list_files:
        # Send an HTTP GET request to the URL and store the response.
        legacy_home = requests.get(
            'https://data.nas.nasa.gov/mcmcref/legacygcm/'
            )
        fv3_home = requests.get(
            'https://data.nas.nasa.gov/mcmcref/fv3betaout1/'
            )

        # Access the text content of the response, which contains the
        # webpage's HTML.
        legacy_dir_text = legacy_home.text
        fv3_file_text = fv3_home.text
                
        # Search for the URLs beginning with the below string
        legacy_subdir_search = (
            "https://data\.nas\.nasa\.gov/legacygcm/legacygcmdata/"
            )

        legacy_urls = re.findall(
            fr"{legacy_subdir_search}[a-zA-Z0-9_\-\.~:/?#\[\]@!$&'()*+,;=]+",
            legacy_dir_text
            )

        print(f"\nAvailable directories:")
        print(f"---------------------")
        for url in legacy_urls:
            legacy_dir_option = url.split('legacygcmdata/')[1]
            print(f"(Legacy MGCM) {legacy_dir_option} {Cyan}({url}){Nclr}")
        print(f"(FV3-based MGCM) FV3BETAOUT1 {Cyan}({fv3_home}){Nclr}")

        if args.directory_name:
            # If a directory is provided, list the files in that directory
            portal_dir=args.directory_name
            if portal_dir == 'FV3BETAOUT1':
                # FV3-based MGCM
                print(f"\nAvailable files from the FV3-based MGCM's FV3BETAOUT1 directory:")
                print(f"-------------------------------------------------")
                print_file_list(re.findall(r'href="[^"]*\/([^"\/]+\.nc)"',
                                            fv3_file_text))
                print(f"\n{Cyan}({fv3_home}){Nclr}")
            
            elif portal_dir in [
                'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS',
                'ACTIVECLDS_NCDF'
                ]:
                # Legacy MGCM
                print(f"\nAvailable files from the Legacy MGCM's {portal_dir} "
                      f"directory:")
                print(f"-------------------------------------------------")
                legacy_dir_url = legacy_subdir_search + portal_dir + r'\/'
                legacy_file_text = (requests.get(legacy_dir_url)).text
                print_file_list(re.findall(r'download="(fort\.11_[0-9]+)"',
                                                legacy_file_text))
                print(f"\n{Cyan}({legacy_dir_url}){Nclr}")

    if args.directory_name and not args.list_files:
        portal_dir=args.directory_name
        if portal_dir in [
            'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS', 'ACTIVECLDS_NCDF'
            ]:
            requested_url = (
                "https://data.nas.nasa.gov/legacygcm/legacygcmdata/"
                + portal_dir
                + "/"
                )
        elif portal_dir in ['FV3BETAOUT1']:
            requested_url = (
                "https://data.nas.nasa.gov/legacygcm/fv3betaout1data/"
                )

        if not (args.ls or args.filename):
            prYellow("ERROR No file requested. Use [-ls --ls] or "
                     "[-f --filename] to specify a file to download.")
            sys.exit(1)  # Return a non-zero exit code

        if args.ls:
            data_input = np.asarray(args.ls)
            if len(data_input) == 1:
                # Query the file that contains this Ls
                closest_index = np.argmin(np.abs(Ls_ini - data_input))
                if data_input < Ls_ini[closest_index]:
                    closest_index -= 1
                file_list = np.arange(closest_index, closest_index + 1)

            elif len(data_input) == 2:
                # Query files within this range of Ls
                i_start = np.argmin(np.abs(Ls_ini - data_input[0]))
                if data_input[0] < Ls_ini[i_start]:
                    i_start -= 1

                i_end = np.argmin(np.abs(Ls_end - data_input[1]))
                if data_input[1] > Ls_end[i_end]:
                    i_end += 1

                file_list = np.arange(i_start, i_end + 1)

            prCyan(f"Saving {len(file_list)} file(s) to {save_dir}")

            for ii in file_list:
                if portal_dir == 'ACTIVECLDS_NCDF':
                    # Legacy .nc files
                    file_name = (
                        f"LegacyGCM_Ls{Ls_ini[ii]:03d}_Ls{Ls_end[ii]:03d}.nc"
                        )
                else:
                    # fort.11 files
                    file_name = f"fort.11_{670+ii:04d}"

                url = requested_url + file_name
                file_name = save_dir + file_name
                print(f"Downloading {url}...")
                download(url, file_name)

        elif args.filename:
            requested_files = np.asarray(args.filename)
            for f in requested_files:
                url = requested_url + f
                file_name = save_dir + f
                print(f"Downloading {url}...")
                download(url, file_name)

    elif not args.list_files:
        # If no directory is provided and its not a -list request
        prYellow("ERROR: A directory must be specified unless using -list.")
        sys.exit(1)

# ------------------------------------------------------
#                  END OF PROGRAM
# ------------------------------------------------------

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
