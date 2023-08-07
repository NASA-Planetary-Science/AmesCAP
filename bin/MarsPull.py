#!/usr/bin/env python3

# Load generic Python modules
import sys        # system command
import os         # access operating systems function

# make print statements appear in color
def prCyan(skk): print("\033[96m{}\033[00m".format(skk))
def prYellow(skk): print("\033[93m{}\033[00m".format(skk))

# try to import specific scientic modules
try:
    import numpy as np
    import argparse     # parse arguments
    import requests
    import time

except ImportError as error_msg:
    prYellow("Error while importing modules")
    print("Error was: " + error_msg.message)
    exit()

except Exception as exception:
    # output unexpected Exceptions
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)
    exit()

# ======================================================
#                  ARGUMENT PARSER
# ======================================================
parser = argparse.ArgumentParser(
    description=("""\033[93mUilities for accessing files on the MCMC 
                NAS Data Portal \033[00m """), 
                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-id', '--id', type=str,
                    help=("Query data by simulation identifier "
                    "corresponding to \n"
                    "a subdirectory of legacygcmdata/:\n"
                    "\033[96mhttps://data.nas.nasa.gov/legacygcm/"
                    "data_legacygcm.php?dir=/legacygcmdata\033[00m\n"
                    "Current options include: "
                    "'\033[93mACTIVECLDS\033[00m', "
                    "'\033[93mINERTCLDS\033[00m', "
                    "and '\033[93mACTIVECLDS_NCDF\033[00m'\n"
                    "> Usage: MarsPull.py -id  INERTCLDS \n\n"))

parser.add_argument('-ls', '--ls', nargs='+', type=float,
                    help="Query data by solar longitude (Ls)\n"
                    "> Usage: MarsPull.py -ls 90.\n"
                    ">        MarsPull.py -ls [start] [stop] \n\n")

parser.add_argument('-f', '--filename', nargs='+', type=str,
                    help=("Query data by filename\n"
                    "> Usage: MarsPull.py -id ACTIVECLDS_NCDF -f "
                    "fort.11_0730 fort.11_0731"))

# ======================================================
#                  DEFINITIONS
# ======================================================
saveDir = os.getcwd() + '/'

# available files by Ls:
Ls_ini = np.array([  0,   5,  10,  15,  19,  24,  29,  34,  38,  43,
                    48,  52,  57,  61,  66,  70,  75,  79,  84,  88,
                    93,  97, 102, 106, 111, 116, 121, 125, 130, 135,
                   140, 146, 151, 156, 162, 167, 173, 179, 184, 190, 
                   196, 202, 209, 215, 221, 228, 234, 241, 247, 254, 
                   260, 266, 273, 279, 286, 292, 298, 304, 310, 316,
                   322, 328, 333, 339, 344, 350, 355])

Ls_end = np.array([  4,   9,  14,  19,  24,  29,  33,  38,  42,  47,
                    52,  56,  61,  65,  70,  74,  79,  83,  88,  92,
                    97, 101, 106, 111, 115, 120, 125, 130, 135, 140,
                   145, 150, 156, 161, 167, 172, 178, 184, 190, 196,
                   202, 208, 214, 221, 227, 233, 240, 246, 253, 259,
                   266, 272, 279, 285, 291, 297, 304, 310, 316, 321,
                   327, 333, 338, 344, 349, 354,   0])


def download(url, filename):
    """
    Downloads a file from a URL to the local computer.
    
    Parameters
    ----------
    URL: str
        The URL to download from. This is built from:
        https://data.nas.nasa.gov/legacygcm/download_data_legacygcm.php?file=/legacygcmdata/
        by appending the simulation ID to the end of the URL.
    filename: str
        The name of the file to download
    
    Raises
    ------
    response.status_code
        A file-not-found error
    
    """

    _, fname = os.path.split(filename)
    # use a context manager to make an HTTP request and file
    response = requests.get(url, stream=True)
    # get the total size, in bytes, from the response header
    total_size = response.headers.get('content-length')

    if response.status_code == 404:
        print('File not found! Error code: ', response.status_code)
    
    else:
        # if the header is found the file size is known. Return a
        # progress bar
        if total_size is not None:
            with open(filename, 'wb') as f:
                downloaded = 0
                if total_size:
                    total_size = int(total_size)
                    # define the size of the chunk to iterate over (Mb)
                    chunk_size = max(int(total_size/1000), 1024*1024)
                # iterate over every chunk and calculate % of total_size
                for i, data in enumerate(response.iter_content(chunk_size=chunk_size)):
                    downloaded += len(data)
                    f.write(data)
                    # calculate current %
                    done = int(50*downloaded/total_size)
                    c = i * chunk_size / total_size * 100
                    # print progress to console then flush console
                    sys.stdout.write('\r[{}{}]'.format(
                        '#' * done, '.' * (50-done)))
                    sys.stdout.flush()
                    print(f"\r{round(c, 4)}%")
                    print(f"\r{round(done, 4)}%")
                    time.sleep(.1)
                    sys.stdout.flush()
            sys.stdout.write('\n')
        
        else:
            # If the header is not found, skip the progressbar
            print('Downloading %s ...' % (fname))
            with open(local_file, 'wb')as f:
                f.write(data.content)
            print('%s Done' % (fname))


# ======================================================
#                  MAIN PROGRAM
# ======================================================
def main():
    simu_ID = parser.parse_args().id

    if simu_ID is None:
        prYellow(
            "***Error*** simulation ID [-id, --id] is required. See 'MarsPull.py -h' for help.")
        exit()

    URLbase = ("https://data.nas.nasa.gov/legacygcm/"
               "download_data_legacygcm.php?file=/legacygcmdata/"
               + simu_ID + "/")

    if parser.parse_args().ls:
        data_input = np.asarray(parser.parse_args().ls)
        if len(data_input) == 1:  # Wuery only the file containing this Ls
            i_start = np.argmin(np.abs(Ls_ini-data_input))
            if data_input < Ls_ini[i_start]:
                i_start -= 1
            i_request = np.arange(i_start, i_start+1)

        elif len(data_input) == 2:  # start, stop is provided
            i_start = np.argmin(np.abs(Ls_ini-data_input[0]))
            if data_input[0] < Ls_ini[i_start]:
                i_start -= 1

            i_end = np.argmin(np.abs(Ls_end-data_input[1]))
            if data_input[1] > Ls_end[i_end]:
                i_end += 1

            i_request = np.arange(i_start, i_end+1)

        print('Saving  %i files in %s ' % (len(i_request), saveDir))
        for ii in i_request:
            # Legacy .nc files
            if simu_ID == 'ACTIVECLDS_NCDF':
                fName = 'LegacyGCM_Ls%03d_Ls%03d.nc' % (Ls_ini[ii], Ls_end[ii])
            # fort.11 files
            else:
                fName = 'fort.11_%04d' % (670+ii)

            url = URLbase+fName
            filename = saveDir+fName
            #print('Downloading '+ fName+ '...')
            print('Downloading ' + url + '...')
            download(url, filename)

    elif parser.parse_args().filename:
        f_input = np.asarray(parser.parse_args().filename)
        for ff in f_input:
            url = URLbase+ff
            filename = saveDir+ff
            print('Downloading ' + url + '...')  # ff
            download(url, filename)
    else:
        prYellow('No data requested. Use -ls or -f to specify data to download.')

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == '__main__':
    main()
