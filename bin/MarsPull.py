#!/usr/bin/env python3
"""
The MarsPull executable is for querying data from the Mars Climate
Modeling Center (MCMC) Mars Global Climate Model (MGCM) repository on
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

    * debug_wrapper - A decorator that wraps a function with error handling
    * download - Downloads a file from the NAS Data Portal (data.nas.nasa.gov).
    * print_file_list - Prints a list of files.
    * main - The main function that handles the command-line arguments
        and coordinates the download process.
        It checks for the presence of the required arguments, validates
        the input, and calls the appropriate functions to download the
        requested files. It also handles the logic for listing available
        directories and files, as well as downloading files based on
        specified solar longitudes (Ls) or file names.
"""

# make print statements appear in color
from amescap.Script_utils import (
    Green, Yellow, Nclr, Cyan, Blue, Red
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
    A decorator that wraps a function with error handling
    based on the --debug flag.
    If the --debug flag is set, it prints the full traceback
    of any exception that occurs. Otherwise, it prints a
    simplified error message.

    :param func: The function to wrap.
    :type  func: function
    :return: The wrapped function.
    :rtype:  function
    :raises Exception: If an error occurs during the function call.
    :raises TypeError: If the function is not callable.
    :raises ValueError: If the function is not found.
    :raises NameError: If the function is not defined.
    :raises AttributeError: If the function does not have the
        specified attribute.
    :raises ImportError: If the function cannot be imported.
    :raises RuntimeError: If the function cannot be run.
    :raises KeyError: If the function does not have the
        specified key.
    :raises IndexError: If the function does not have the
        specified index.
    :raises IOError: If the function cannot be opened.
    :raises OSError: If the function cannot be accessed.
    :raises EOFError: If the function cannot be read.
    :raises MemoryError: If the function cannot be allocated.
    :raises OverflowError: If the function cannot be overflowed.
    :raises ZeroDivisionError: If the function cannot be divided by zero.
    :raises StopIteration: If the function cannot be stopped.
    :raises KeyboardInterrupt: If the function cannot be interrupted.
    :raises SystemExit: If the function cannot be exited.
    :raises AssertionError: If the function cannot be asserted.
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global debug
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if debug:
                # In debug mode, show the full traceback
                print(f'{Red}ERROR in {func.__name__}: {str(e)}{Nclr}')
                traceback.print_exc()
            else:
                # In normal mode, show a clean error message
                print(f'{Red}ERROR in {func.__name__}: {str(e)}\nUse '
                      f'--debug for more information.{Nclr}')
            return 1  # Error exit code
    return wrapper


# ------------------------------------------------------
#                  ARGUMENT PARSER
# ------------------------------------------------------

parser = argparse.ArgumentParser(
    prog=('MarsPull'),
    description=(
        f'{Yellow}Uility for downloading NASA Ames Mars Global Climate '
        f'Model output files from the NAS Data Portal at:\n'
        f'{Cyan}https://data.nas.nasa.gov/mcmcref/\n{Nclr}'
        f'Requires ``-f`` or ``-ls``.'
        f'{Nclr}\n\n'
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument('directory_name', type=str, nargs='?',
    choices=[
        'FV3BETAOUT1', 'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS',
        'ACTIVECLDS_NCDF'],
    help=(
        f'Selects the simulation directory from the '
        f'NAS data portal ('
        f'{Cyan}https://data.nas.nasa.gov/mcmcref/){Nclr}\n'
        f'Current directory options are:\n{Yellow}FV3BETAOUT1, ACTIVECLDS, '
        f'ACTIVECLDS, INERTCLDS, NEWBASE_ACTIVECLDS, ACTIVECLDS_NCDF\n'
        f'{Red}MUST be used with either ``-f`` or ``-ls``\n'
        f'{Green}Example:\n'
        f'> MarsPull INERTCLDS -f fort.11_0690\n'
        f'{Blue}OR{Green}\n'
        f'> MarsPull INERTCLDS -ls 90\n'
        f'{Nclr}\n\n'
    )
)

parser.add_argument('-list', '--list_files', action='store_true',
    help=(
        f'Return a list of the directories and files available for download '
        f'from {Cyan}https://data.nas.nasa.gov/mcmcref/{Nclr}\n'
        f'{Green}Example:\n'
        f'> MarsPull -list {Blue}# lists all directories{Green}\n'
        f'> MarsPull -list ACTIVECLDS {Blue}# lists files under ACTIVECLDS '
        f'{Nclr}\n\n'
    )
)

parser.add_argument('-f', '--filename', nargs='+', type=str,
    help=(
        f'The name(s) of the file(s) to download.\n'
        f'{Green}Example:\n'
        f'> MarsPull INERTCLDS -f fort.11_0690'
        f'{Nclr}\n\n'
    )
)

parser.add_argument('-ls', '--ls', nargs='+', type=float,
    help=(
        f'Selects the file(s) to download based on a range of solar '
        f'longitudes (Ls).\n'
        f'This only works on data in the {Yellow}ACTIVECLDS{Nclr} and '
        f'{Yellow}INERTCLDS{Nclr} folders.\n'
        f'{Green}Example:\n'
        f'> MarsPull INERTCLDS -ls 90\n'
        f'> MarsPull INERTCLDS -ls 90 180'
        f'{Nclr}\n\n'
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('--debug', action='store_true',
    help=(
        f'Use with any other argument to pass all Python errors and\n'
        f'status messages to the screen when running CAP.\n'
        f'{Green}Example:\n'
        f'> MarsPull INERTCLDS -ls 90 --debug'
        f'{Nclr}\n\n'
    )
 )

args = parser.parse_args()
debug = args.debug


# ------------------------------------------------------
#                  DEFINITIONS
# ------------------------------------------------------
save_dir = (f'{os.getcwd()}/')

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
    The function takes a URL and a file name as input, and downloads the
    file from the URL, saving it to the specified file name. It also
    provides a progress bar to show the download progress if the file
    size is known. If the file size is unknown, it simply downloads the
    file without showing a progress bar.
    The function handles errors during the download process and prints
    appropriate messages to the console.

    :param url: The url to download from, e.g.,
        'https://data.nas.nasa.gov/legacygcm/fv3betaout1data/03340.fixed.nc'
    :type  url: str
    :param file_name: The local file_name e.g.,
        '/lou/la4/akling/Data/LegacyGCM_Ls000_Ls004.nc'
    :type  file_name: str
    :return: The requested file(s), downloaded and saved to the current
        directory.
    :rtype:  None
    :raises FileNotFoundError: A file-not-found error.
    :raises PermissionError: A permission error.
    :raises OSError: An operating system error.
    :raises ValueError: A value error.
    :raises TypeError: A type error.
    :raises requests.exceptions.RequestException: A request error.
    :raises requests.exceptions.HTTPError: An HTTP error.
    :raises requests.exceptions.ConnectionError: A connection error.
    :raises requests.exceptions.Timeout: A timeout error.
    :raises requests.exceptions.TooManyRedirects: A too many redirects
        error.
    :raises requests.exceptions.URLRequired: A URL required error.
    :raises requests.exceptions.InvalidURL: An invalid URL error.
    :raises requests.exceptions.InvalidSchema: An invalid schema error.
    :raises requests.exceptions.MissingSchema: A missing schema error.
    :raises requests.exceptions.InvalidHeader: An invalid header error.
    :raises requests.exceptions.InvalidProxyURL: An invalid proxy URL
        error.
    :raises requests.exceptions.InvalidRequest: An invalid request error.
    :raises requests.exceptions.InvalidResponse: An invalid response
        error.
    """

    _, fname = os.path.split(file_name)
    response = requests.get(url, stream=True)
    total = response.headers.get('content-length')

    if response.status_code == 404:
        print(f'{Red}Error during download, error code is: '
              f'{response.status_code}{Nclr}')
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
                    sys.stdout.write(
                        f'\rProgress: '
                        f'[{"#"*status}{"."*(50 - status)}] {status}%'
                        )
                    sys.stdout.flush()
            sys.stdout.write('\n\n')
        else:
            # If file size is unknown, skip progress bar
            print(f'Downloading {fname}...')
            with open(file_name, 'wb')as f:
                f.write(response.content)
            print(f'{fname} Done')


def print_file_list(list_of_files):
    """
    Prints a list of files.

    :param list_of_files: The list of files to print.
    :type  list_of_files: list
    :return: None
    :rtype:  None
    :raises TypeError: If list_of_files is not a list.
    :raises ValueError: If list_of_files is empty.
    :raises IndexError: If list_of_files is out of range.
    :raises KeyError: If list_of_files is not found.
    :raises OSError: If list_of_files is not accessible.
    :raises IOError: If list_of_files is not open.
    """
    for file in list_of_files:
        print(file)


# ------------------------------------------------------
#                  MAIN FUNCTION
# ------------------------------------------------------

@debug_wrapper
def main():
    """
    The main function that handles the command-line arguments

    Handles the command-line arguments and coordinates the download
    process. It checks for the presence of the required arguments,
    validates the input, and calls the appropriate functions to download
    the requested files. It also handles the logic for listing available
    directories and files, as well as downloading files based on
    specified solar longitudes (Ls) or file names.

    :return: 0 if successful, 1 if an error occurred.
    :rtype:  int
    :raises SystemExit: If an error occurs during the execution of the
        program, the program will exit with a non-zero status code.
    """
    global debug

    if not args.list_files and not args.directory_name:
        print('Error: You must specify either -list or a directory.')
        sys.exit(1)

    base_dir = 'https://data.nas.nasa.gov'
    legacy_home_url = f'{base_dir}/mcmcref/legacygcm/'
    legacy_data_url = f'{base_dir}/legacygcm/legacygcmdata/'
    fv3_home_url = f'{base_dir}/mcmcref/fv3betaout1/'
    fv3_data_url = f'{base_dir}/legacygcm/fv3betaout1data/'

    if args.list_files:
        # Send an HTTP GET request to the URL and store the response.
        legacy_home_html = requests.get(f'{legacy_home_url}')
        fv3_home_html = requests.get(f'{fv3_home_url}')

        # Access the text content of the response, which contains the
        # webpage's HTML.
        legacy_dir_text = legacy_home_html.text
        fv3_dir_text = fv3_home_html.text

        # Search for the URLs beginning with the below string
        legacy_subdir_search = (
            'https://data\.nas\.nasa\.gov/legacygcm/legacygcmdata/'
            )
        fv3_subdir_search = (
            'https://data\.nas\.nasa\.gov/legacygcm/fv3betaout1data/'
            )

        legacy_urls = re.findall(
            fr'{legacy_subdir_search}[a-zA-Z0-9_\-\.~:/?#\[\]@!$&"()*+,;=]+',
            legacy_dir_text
            )

        # NOTE: The FV3-based MGCM data only has one directory and it is
        #       not listed in the FV3BETAOUT1 directory. The URL is
        #       hardcoded below. The regex below is commented out, but
        #       left in place in case the FV3BETAOUT1 directory is
        #       updated with subdirectories in the future.
        # fv3_urls = re.findall(
        #     fr'{fv3_subdir_search}[a-zA-Z0-9_\-\.~:/?#\[\]@!$&"()*+,;=]+',
        #     fv3_dir_text
        #     )
        fv3_urls = [f'{fv3_data_url}']

        print(f'\nSearching for available directories...')
        print(f'--------------------------------------')
        for url in legacy_urls:
            legacy_dir_option = url.split('legacygcmdata/')[1]
            print(f'{"(Legacy MGCM)":<17} {legacy_dir_option:<20} '
                  f'{Cyan}{url}{Nclr}')

        # NOTE: See above comment for the FV3-based MGCM data note
        # for url in fv3_urls:
        #     fv3_dir_option = url.split('fv3betaout1data/')[1]
        #     print(f'{"(FV3-based MGCM)":<17} {fv3_dir_option:<17} '
        #           f'{Cyan}{url}{Nclr}')
        print(f'{"(FV3-based MGCM)":<17} {"FV3BETAOUT1":<20} '
              f'{Cyan}{fv3_home_url}{Nclr}')
        print(f'---------------------\n')

        if args.directory_name:
            # If a directory is provided, list the files in that directory
            portal_dir = args.directory_name
            if portal_dir == 'FV3BETAOUT1':
                # FV3-based MGCM
                print(f'\n{Green}Selected: (FV3-based MGCM) FV3BETAOUT1{Nclr}')
                print(f'\nSearching for available files...')
                print(f'--------------------------------')
                fv3_dir_url = f'{fv3_home_url}'
                fv3_data = requests.get(fv3_dir_url)
                fv3_file_text = fv3_data.text

                # This looks for download attributes or href links
                # ending with the .nc pattern
                fv3_files_available = []

                # Try multiple patterns to find .nc files
                download_files = re.findall(
                    r'download="([^"]+\.nc)"',
                    fv3_file_text
                    )
                if download_files:
                    fv3_files_available = download_files
                else:
                    # Look for href attributes with .nc files
                    href_files = re.findall(
                        r'href="[^"]*\/([^"\/]+\.nc)"',
                        fv3_file_text
                        )
                    if href_files:
                        fv3_files_available = href_files
                    else:
                        # Look for links with .nc text
                        link_files = re.findall(
                            r'<a[^>]*>([^<]+\.nc)</a>',
                            fv3_file_text
                            )
                        if link_files:
                            fv3_files_available = link_files

                # Filter out any potential HTML or Javascript that might
                # match the pattern
                fv3_files_available = [f for f in fv3_files_available if (
                    not f.startswith('<') and
                    not f.startswith('function') and
                    not f.startswith('var') and
                    '.nc' in f
                )]

                # Sort the files
                fv3_files_available.sort()

                # Print the files
                if fv3_files_available:
                    print_file_list(fv3_files_available)
                else:
                    print('No .nc files found. Run with --debug for more info')
                    if debug:
                        # Try a different approach for debugging
                        table_rows = re.findall(
                            r'<tr>.*?</tr>',
                            fv3_file_text,
                            re.DOTALL
                            )
                        for row in table_rows:
                            if '.nc' in row:
                                print(f'Debug - Found row with .nc: {row}')

                print(f'---------------')
                # The download URL differs from the listing URL
                print(f'{Cyan}({fv3_dir_url}){Nclr}\n')

                print(f'{Yellow}You can download files using the -f '
                      f'option with the directory name, e.g.\n'
                      f'> MarsPull FV3BETAOUT1 -f 03340.fixed.nc\n'
                      f'> MarsPull FV3BETAOUT1 -f 03340.fixed.nc '
                      f'03340.atmos_average.nc{Nclr}\n')

            elif portal_dir in [
                'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS',
                'ACTIVECLDS_NCDF'
                ]:
                # Legacy MGCM
                print(f'\n{Green}Selected: (Legacy MGCM) {portal_dir}{Nclr}')
                print(f'\nAvailable files:')
                print(f'---------------')
                legacy_dir_url = (f'{legacy_data_url}' + portal_dir + r'/')
                legacy_data = requests.get(legacy_dir_url)
                legacy_file_text = legacy_data.text

                # This looks for download attributes or href links
                # ending with the fort.11_ pattern
                legacy_files_available = []

                # First try to find download attributes which are more reliable
                download_files = re.findall(
                    r'download="(fort\.11_[0-9]+)"',
                    legacy_file_text
                    )
                if download_files:
                    legacy_files_available = download_files
                else:
                    # Fallback to looking for href links with fort.11_ pattern
                    href_files = re.findall(
                        r'href="[^"]*\/?(fort\.11_[0-9]+)"',
                        legacy_file_text
                        )
                    if href_files:
                        legacy_files_available = href_files
                    # If still empty, try another pattern to match links
                    if not legacy_files_available:
                        href_files = re.findall(
                            r'<a href="[^"]*"[^>]*>(fort\.11_[0-9]+)</a>',
                            legacy_file_text
                            )
                        legacy_files_available = href_files

                print_file_list(legacy_files_available)
                print(f'---------------')
                print(f'{Cyan}({legacy_dir_url}){Nclr}\n')

                print(f'{Yellow}You can download these files using the '
                      f'-f or -ls options with the directory name, e.g.\n'
                      f'> MarsPull ACTIVECLDS -f fort.11_0690\n'
                      f'> MarsPull ACTIVECLDS -f fort.11_0700 fort.11_0701 \n'
                      f'> MarsPull ACTIVECLDS -ls 90\n'
                      f'> MarsPull ACTIVECLDS -ls 90 180{Nclr}\n')

            else:
                print(f'Error: Directory {portal_dir} does not exist.')
                sys.exit(1)
            sys.exit(0)

        else:
            # If no directory is provided, exit with an error
            print(f'{Yellow}You can list the files in a directory by using '
                  f'the -list option with a directory name, e.g.\n'
                  f'> MarsPull -list ACTIVECLDS{Nclr}\n')

    if args.directory_name and not args.list_files:
        portal_dir = args.directory_name
        if portal_dir in [
            'ACTIVECLDS', 'INERTCLDS', 'NEWBASE_ACTIVECLDS', 'ACTIVECLDS_NCDF'
            ]:
            requested_url = (f'{legacy_data_url}' + portal_dir + '/')
        elif portal_dir in ['FV3BETAOUT1']:
            requested_url = (f'{fv3_data_url}')

        if not (args.ls or args.filename):
            print(f'{Yellow}ERROR No file requested. Use [-ls --ls] or '
                  f'[-f --filename] to specify a file to download.{Nclr}')
            sys.exit(1)  # Return a non-zero exit code
        portal_dir = args.directory_name

        if portal_dir == 'FV3BETAOUT1' and args.ls:
            print(f'{Red}ERROR: The FV3BETAOUT1 directory does not support '
                f'[-ls --ls] queries. Please query by file name(s) '
                f'[-f --filename], e.g.\n'
                f'> MarsPull FV3BETAOUT1 -f 03340.fixed.nc{Nclr}')
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

            for ii in file_list:
                if portal_dir == 'ACTIVECLDS_NCDF':
                    # Legacy .nc files
                    file_name = (
                        f'LegacyGCM_Ls{Ls_ini[ii]:03d}_Ls{Ls_end[ii]:03d}.nc'
                        )
                else:
                    # fort.11 files
                    file_name = f'fort.11_{670+ii:04d}'

                url = requested_url + file_name
                file_name = save_dir + file_name
                print(f'\nDownloading {Cyan}{url}{Nclr}...')
                print(f'Saving {Cyan}{len(file_list)}{Nclr} file(s) to '
                      f'{Cyan}{save_dir}{Nclr}')
                download(url, file_name)

        elif args.filename:
            requested_files = np.asarray(args.filename)
            for f in requested_files:
                url = requested_url + f
                file_name = save_dir + f
                print(f'\nDownloading {url}...')
                download(url, file_name)

    elif not args.list_files:
        # If no directory is provided and its not a -list request
        print(f'{Yellow}ERROR: A directory must be specified unless using '
              f'-list.{Nclr}')
        sys.exit(1)

# ------------------------------------------------------
#                  END OF PROGRAM
# ------------------------------------------------------

if __name__ == '__main__':
    exit_code = main()
    sys.exit(exit_code)
