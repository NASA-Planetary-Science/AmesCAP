#!/usr/bin/env python3
"""
The MarsPull executable is for querying data from the Mars Climate \
Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on \
the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

The executable requires 2 arguments:
    * ``[-id --id]``      The simulation identifier, AND
    * ``[-ls --ls]``      The desired solar longitude(s), OR
    * ``[-f --filename]`` The name(s) of the desired file(s)

Third-party Requirements:
    * ``numpy``
    * ``argparse``
    * ``requests``

List of Functions:
    * download - Queries the requested file from the NAS Data Portal.
"""

# make print statements appear in color
from amescap.Script_utils import (
    prYellow, prCyan, Green, Yellow, Nclr, Cyan
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import requests     # Download data from website
import numpy as np

# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}Uility for querying files on the MCMC NAS Data "
        f"Portal.{Nclr}"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-id", "--id", type=str,
    help=(
        f"Query data by simulation identifier corresponding to \n"
        f"a subdirectory of :\n"
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/ \n"
        f"Current options include: '{Yellow}FV3BETAOUT1{Nclr}' '{Yellow}ACTIVECLDS{Nclr}', "
        f"'{Yellow}INERTCLDS{Nclr}', {Yellow}NEWBASE_ACTIVECLDS{Nclr}  and '{Yellow}ACTIVECLDS_NCDF\n"
        f"{Green}Usage:\n"
        f"> MarsPull -id  INERTCLDS..."
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-f", "--filename", nargs="+", type=str,
    help=(
        f"Query data by file name. Requires a simulation identifier "
        f"(--id)\n"
        f"{Green}Usage:\n"
        f"> MarsPull -id ACTIVECLDS -f fort.11_0730 fort.11_0731"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("-ls", "--ls", nargs="+", type=float,
    help=(
        f"Legacy GCM only: Query data by solar longitude (Ls). Requires a simulation "
        f"identifier (--id)\n"
        f"{Green}Usage:\n"
        f"> MarsPull -id ACTIVECLDS -ls 90.\n"
        f"> MarsPull -id ACTIVECLDS -ls [start] [stop]"
        f"{Nclr}\n\n"
    )
)

parser.add_argument("--debug", action="store_true",
    help=(
        f"More verbosity in status and error messages when running CAP."
        f"\n\n"
    )
 )

# ======================================================
#                  DEFINITIONS
# ======================================================


saveDir = (f"{os.getcwd()}/")

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



def download(url, filename):
    """
    Downloads a file from the NAS Data Portal (data.nas.nasa.gov).

    Parameters
    ----------
    url : str
        The url to download, e.g   'https://data.nas.nasa.gov/legacygcm/download_data.php?file=/legacygcmdata/LegacyGCM_Ls000_Ls004.nc'
    filename : str
        The local filename e.g  '/lou/la4/akling/Data/LegacyGCM_Ls000_Ls004.nc'

    Returns
    -------
    The requested file(s), downloaded and saved to the current \
    directory.


    Raises
    ------
    A file-not-found error.
    """

    _ , fname=os.path.split(filename)
    response = requests.get(url, stream=True)
    total = response.headers.get('content-length')

    if response.status_code == 404:
        print('Error during download, error code is: ',response.status_code)
    else:

        #If we have access to the size of the file, return progress bar
        if total is not None:
            with open(filename, 'wb') as f:
                downloaded = 0
                if total :total = int(total)
                for data in response.iter_content(chunk_size=max(int(total/1000), 1024*1024)):
                    downloaded += len(data)
                    f.write(data)
                    status = int(50*downloaded/total)
                    sys.stdout.write(f"\r[{'#'*status}{'.'*(50-status)}]")
                    sys.stdout.flush()
            sys.stdout.write('\n')
        else:

            #Header is unknown yet, skip the use of a progress bar
            print('Downloading %s ...'%(fname))
            with open(filename, 'wb')as f:
                f.write(response.content)
            print('%s Done'%(fname))

# ======================================================
#                  MAIN PROGRAM
# ======================================================
def main():

    #Original
    #URLbase="https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/"
    simu_ID=parser.parse_args().id

    if simu_ID is None :
        prYellow("***Error*** simulation ID [-id --id] is required. See 'MarsPull -h' for help")
        exit()

    #URLbase='https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/'+simu_ID+'/'
    print('new URL base')
    if simu_ID in ['ACTIVECLDS','INERTCLDS', 'NEWBASE_ACTIVECLDS','ACTIVECLDS_NCDF']:
        URLbase='https://data.nas.nasa.gov/legacygcm/legacygcmdata/'+simu_ID+'/'
    elif simu_ID in ['FV3BETAOUT1']:
        URLbase='https://data.nas.nasa.gov/legacygcm/fv3betaout1data/'

    if parser.parse_args().ls :
        data_input=np.asarray(parser.parse_args().ls)
        if len(data_input)==1: #query only  the file that contains this Ls
            i_start=np.argmin(np.abs(Ls_ini-data_input))
            if data_input<Ls_ini[i_start]:i_start-=1
            num_files=np.arange(i_start,i_start+1)

        elif len(data_input)==2: #start stop  is provided
            i_start=np.argmin(np.abs(Ls_ini-data_input[0]))
            if data_input[0]<Ls_ini[i_start]:i_start-=1

            i_end=np.argmin(np.abs(Ls_end-data_input[1]))
            if data_input[1]>Ls_end[i_end]:i_end+=1

            num_files=np.arange(i_start,i_end+1)
        prCyan(f"Saving {len(num_files)} file(s) to {saveDir}")
        for ii in num_files:
            #Legacy .nc files
            if simu_ID=='ACTIVECLDS_NCDF':
                fName='LegacyGCM_Ls%03d_Ls%03d.nc'%(Ls_ini[ii],Ls_end[ii])
            #fort.11 files
            else:
                fName='fort.11_%04d'%(670+ii)

            url = URLbase+fName
            filename=saveDir+fName
            #print('Downloading '+ fName+ '...')
            print('Downloading '+ url+ '...')
            download(url,filename)

    elif parser.parse_args().filename:
        f_input=np.asarray(parser.parse_args().filename)
        for ff in f_input :
            url = URLbase+ff
            filename=saveDir+ff
            print('Downloading '+ url+ '...')#ff
            download(url,filename)
    else:
        prYellow("ERROR No file requested. Use [-ls --ls] or "
                 "[-f --filename] with [-id --id] to specify a file to "
                 "download.")
        exit()

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
