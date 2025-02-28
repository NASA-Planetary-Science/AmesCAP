#!/usr/bin/env python3
"""
The MarsPull executable is for querying data from the Mars Climate \
Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on \
the NASA NAS Data Portal at data.nas.nasa.gov/mcmc.

The executable requires 2 arguments:
    * The directory from which to pull data from 
    (https://data.nas.nasa.gov/mcmcref/), AND
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
    prYellow, prCyan, Green, Yellow, Nclr, Cyan, Blue, Red
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import requests     # Download data from website
import re           # Regular expressions
import numpy as np
import requests     # make HTTP requests
# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
    prog=('MarsPull'),
    description=(
        f"{Yellow}Uility for downloading NASA Ames Mars Global Climate "
        f"Model output files from the NAS Data Portal at:"
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/\n{Nclr}\n"
        f"Requires the ``-id`` argument AND EITHER ``-f`` or ``-ls``."
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('-list', '--list_files', action='store_true',
    help=(
        f"Return a list of all the files available for download from:\n"
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/\n{Nclr}\n"
        f"{Green}Example:\n"
        f"> MarsPull -list"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('directory_name', type=str,
    choices=[
        'FV3BETAOUT1', 'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS',
        'ACTIVECLDS_NCDF'],
    help=(
        f"Selects the simulation directory from the "
        f"NAS data portal:\n"
        f"{Cyan}https://data.nas.nasa.gov/mcmcref/\n{Nclr}\n"
        f"Current options are:\n{Yellow}FV3BETAOUT1\nACTIVECLDS\n"
        f"INERTCLDS\nNEWBASE_ACTIVECLDS\nACTIVECLDS_NCDF\n"
        f"{Red}MUST be used with either ``-f`` or ``-ls``.\n"
        f"{Green}Example:\n"
        f"> MarsPull ACTIVECLDS -f fort.11_0730\n"
        f"{Blue}OR{Green}\n"
        f"> MarsPull ACTIVECLDS -ls 90\n"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-f', '--filename', nargs='+', type=str,
    help=(
        f"The name(s) of the file(s) to download.\n"
        f"{Green}Example:\n"
        f"> MarsPull ACTIVECLDS -f fort.11_0730 fort.11_0731"
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
        f"> MarsPull ACTIVECLDS -ls 90\n"
        f"> MarsPull ACTIVECLDS -ls 180 360"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsPull ACTIVECLDS -ls 90 --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()

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

def file_list(list_of_files):
    print("Available files:")
    for file in list_of_files:
        print(file)

# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
    if args.list_files:
        # Send an HTTP GET request to the URL and store the response.
        legacy_data = requests.get('https://data.nas.nasa.gov/mcmcref/legacygcm/')
        fv3_data = requests.get('https://data.nas.nasa.gov/mcmcref/fv3betaout1/')

        # Access the text content of the response, which contains the webpage's HTML.
        legacy_dir_text = legacy_data.text
        fv3_dir_text = fv3_data.text

        # Search for the URLs beginning with the below string
        legacy_dir_search = "https://data\.nas\.nasa\.gov/legacygcm/legacygcmdata/"

        legacy_urls = re.findall(
            fr"{legacy_dir_search}[a-zA-Z0-9_\-\.~:/?#\[\]@!$&'()*+,;=]+", legacy_dir_text
            )

        print("Available directories:")
        for url in legacy_urls:
            dir_option = url.split('legacygcmdata/')[1]
            print(dir_option)

        legacy_data = requests.get(legacy_urls[0])
        legacy_dir_text = legacy_data.text

        legacy_files_available = re.findall(r'download="(fort\.11_[0-9]+)"', legacy_dir_text)
        fv3_files_available = re.findall(r'href="[^"]*\/([^"\/]+\.nc)"', fv3_dir_text)


        file_list(legacy_files_available)
        file_list(fv3_files_available)
    
    portal_dir=args.directory_name
    if portal_dir in ['ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS', 'ACTIVECLDS_NCDF']:
        url_requested="https://data.nas.nasa.gov/legacygcm/legacygcmdata/"+portal_dir+'/'
    elif portal_dir in ['FV3BETAOUT1']:
        url_requested="https://data.nas.nasa.gov/legacygcm/fv3betaout1data/"

    if args.ls :
        data_input=np.asarray(args.ls)
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
            if portal_dir=='ACTIVECLDS_NCDF':
                file_name='LegacyGCM_Ls%03d_Ls%03d.nc'%(Ls_ini[ii],Ls_end[ii])
            #fort.11 files
            else:
                file_name='fort.11_%04d'%(670+ii)

            url = url_requested+file_name
            filename=saveDir+file_name
            #print('Downloading '+ file_name+ '...')
            print('Downloading '+ url+ '...')
            download(url,filename)

    elif args.filename:
        f_input=np.asarray(args.filename)
        for ff in f_input :
            url = url_requested+ff
            filename=saveDir+ff
            print('Downloading '+ url+ '...')#ff
            download(url,filename)
    else:
        prYellow("ERROR No file requested. Use [-ls --ls] or "
                 "[-f --filename] to specify a file to download.")
        exit()

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
