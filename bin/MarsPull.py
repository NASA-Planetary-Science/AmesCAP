#!/usr/bin/env python3
"""
The MarsPull executable is for querying data from the Mars Climate \
Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on \
the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

The executable requires:
    * ``[-id --id]``      The simulation identifier, AND
    * ``[-ls --ls]``      the desired solar longitude(s), OR
    * ``[-f --filename]`` the name(s) of the desired file(s).

Third-party Requirements:
    * ``numpy``
    * ``sys``
    * ``argparse``
    * ``os``
    * ``requests``
"""

# Make print statements appear in color
from amescap.Script_utils import (
    prYellow, prCyan, Yellow, NoColor, Green, Cyan
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import requests     # Download data from website
import numpy as np

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}Uility for querying files on the MCMC NAS Data "
        f"Portal."
        f"{NoColor}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-id", "--id", type=str,
    help=(
        f"Query data by simulation identifier corresponding to \n"
        f"a subdirectory of legacygcmdata/:\n"
        f"{Cyan}https://data.nas.nasa.gov/legacygcm/data_legacygcm.php?"
        f"dir=/legacygcmdata{NoColor}\n"
        f"Current options include: ``{Yellow}ACTIVECLDS{NoColor}``, "
        f"``{Yellow}INERTCLDS{NoColor}``, and "
        f"``{Yellow}ACTIVECLDS_NCDF``\n"
        f"{Green}Usage:\n"
        f"> MarsPull.py -id  INERTCLDS..."
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-ls", "--ls", nargs="+", type=float,
    help=(
        f"Query data by solar longitude (Ls). Requires a simulation "
        f"identifier (``--id``)\n"
        f"{Green}Usage:\n"
        f"> MarsPull.py -id ACTIVECLDS -ls 90.\n"
        f"> MarsPull.py -id ACTIVECLDS -ls [start] [stop]"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("-f", "--filename", nargs="+", type=str,
    help=(
        f"Query data by file name. Requires a simulation identifier "
        f"(``--id``)\n"
        f"{Green}Usage:\n"
        f"> MarsPull.py -id ACTIVECLDS -f fort.11_0730 fort.11_0731"
        f"{NoColor}\n\n"
    )
)

parser.add_argument("--debug", action="store_true",
    help = (f"Debug flag: do not bypass errors.\n\n"))


# ======================================================
#                  DEFINITIONS
# ======================================================


saveDir = (f"{os.getcwd()}/")

global SAVEDIR, ls_0, ls_N

SAVEDIR = (f"{os.getcwd()}/")

# available files by Ls:
Ls_ini = np.array([
    0, 5, 10, 15, 19, 24, 29, 34, 38, 43, 48, 52, 57, 61, 66, 70, 75,
    79, 84, 88, 93, 97, 102, 106, 111, 116, 121, 125, 130, 135, 140,
    146, 151, 156, 162, 167, 173, 179, 184, 190, 196, 202, 209, 215,
    221, 228, 234, 241, 247, 254, 260, 266, 273, 279, 286, 292, 298,
    304, 310, 316,322, 328, 333, 339, 344, 350, 355
])
# available files by Ls:
ls_0 = np.array([
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

ls_N = np.array([
    4, 9, 14, 19, 24, 29, 33, 38, 42, 47, 52, 56, 61, 65, 70, 74, 79,
    83, 88, 92, 97, 101, 106, 111, 115, 120, 125, 130, 135, 140, 145,
    150, 156, 161, 167, 172, 178, 184, 190, 196, 202, 208, 214, 221,
    227, 233, 240, 246, 253, 259, 266, 272, 279, 285, 291, 297, 304,
    310, 316, 321, 327, 333, 338, 344, 349, 354, 0
])


def download(file_name, simulation_id):
    """
    Downloads a file from the MCMC Legacy GCM directory on the NAS Data
    Portal (data.nas.nasa.gov).

    This function specifies the file to download by appending to the \
    URL to the subdirectory, indicated by the user-specified \
    simulation identifier ``[-id --id]``, and the name of the file. The\
    file name is either provided by the user directly using \
    ``[-f --filename]`` or determined based on the user-specified solar\
    longitude ``[-ls --ls]``.

    Parameters
    ----------
    :param simulation_id: The simulation identifier, i.e., the name of \
        the directory to query from: \
        https://data.nas.nasa.gov/mcmc/data_legacygcm.php
    :type simulation_id: str
    :param file_name: The name of the file to download.
    :type file_name: str

    :raises: file-not-found error.

    :return: downloads & saves the requested file(s) to the current \
        directory.
    """

    baseURL = (
        f"https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.\
        php?file=/legacygcmdata/{simulation_id}/")

    URL = baseURL + file_name

    filename = SAVEDIR + file_name

    # Use a context manager to make an HTTP request and file
    rsp = requests.get(URL, stream=True)

    # Get the total size, in bytes, from the response header
    total_size = rsp.headers.get("content-length")

    if rsp.status_code == 404:
        prYellow(f"ERROR File not found! Error code: {rsp.status_code}")
        exit()

    else:
        # Download the data and show a progress bar
        prCyan(f"Downloading {file_name}...")

        with open(filename, "wb") as f:
            downloaded = 0
            if total_size:
                # Convert total_size from str to int
                total_size = int(total_size)
                # Define the size of the chunk to iterate over (Mb)
                chunk_size = max(int(total_size/1000), 1024*1024)

            # Iterate over every chunk and calculate % of total_size
            for chunk in rsp.iter_content(chunk_size=chunk_size):
                downloaded += len(chunk)
                f.write(chunk)

                # Calculate current %
                status = int(50*downloaded/total_size)

                # Print progress to console then flush console
                sys.stdout.write(f"\r[{'#'*status}{'.'*(50-status)}]")
                sys.stdout.flush()
        sys.stdout.write("\n")


# ======================================================================
#                           MAIN PROGRAM
# ======================================================================

def main():
    # parse out the name of the file
    simulation_id = parser.parse_args().id

    if simulation_id is None:
        prYellow("ERROR [-id, --id] is required. Use 'MarsPull.py -h' \
            for additional help.")
        exit()

    if parser.parse_args().ls:
        # If the user input an ID and specified 1+ Ls
        ls_input = np.asarray(parser.parse_args().ls)

        if len(ls_input) == 1:
            # If the user input 1 Ls
            i_start = np.argmin(np.abs(ls_0 - ls_input))
            if ls_input < ls_0[i_start]:
                i_start -= 1
            i_end = i_start

        elif len(ls_input) == 2:
            # If the user input a range of Ls
            i_start = np.argmin(np.abs(ls_0 - ls_input[0]))
            if ls_input[0] < ls_0[i_start]:
                i_start -= 1
            i_end = np.argmin(np.abs(ls_N - ls_input[1]))
            if ls_input[1] > ls_N[i_end]:
                i_end += 1

        num_files = np.arange(i_start, i_end+1)
        prCyan(f"Saving {len(num_files)} file(s) to {SAVEDIR}")

        for n in num_files:
            # For netCDF files
            if simulation_id == "ACTIVECLDS_NCDF":
                file_name = (
                    f"LegacyGCM_Ls{ls_0[n]:03d}_Ls{ls_N[n]:03d}.nc"
                )
            # For fort.11 files
            else:
                file_name = (f"fort.11_{(670 + n):04d}")
            # Trigger the file download
            download(file_name, simulation_id)

    elif parser.parse_args().filename:
        # If the user input an ID and a file name
        file_name_in = np.asarray(parser.parse_args().filename)
        prCyan(f"Saving {len(file_name_in)} file(s) to {SAVEDIR}")
        for file_name in file_name_in:
            # trigger the file download
            download(file_name, simulation_id)
    else:
        # If the user did not specify Ls or a file name
        prYellow("ERROR No file requested. Use ``[-ls --ls]`` or \
            ``[-f --filename]`` with ``[-id --id]`` to specify a file \
            to download.")
        exit()

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == "__main__":
    main()
