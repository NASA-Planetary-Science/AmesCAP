#!/usr/bin/env python3
"""
The MarsFiles executable is for ...

The executable requires x arguments:
    * [-x --x]      define

Third-party Requirements:
    * numpy
    * argparse
    * requests

List of Functions:
    * x
"""

# make print statements appear in color
from amescap.Script_utils import prYellow, prCyan, prRed
Cyan = "\033[96m"
Blue = '\033[94m'
Yellow = "\033[93m"
nocolor = "\033[00m"
Red = "\033[91m"
Green = "\033[92m"
Purple = "\033[95m"

# load generic Python Modules
import argparse     # parse arguments
import sys          # system commands
import os           # access operating system functions
import subprocess   # run command-line commands
import numpy as np
from netCDF4 import Dataset
import warnings     # suppress errors triggered by NaNs

# load amesCAP modules
from amescap.Ncdf_wrapper import Ncdf, Fort
from amescap.FV3_utils import tshift, daily_to_average, daily_to_diurn
from amescap.Script_utils import (find_tod_in_diurn, FV3_file_type, 
                                  filter_vars, regrid_Ncfile, 
                                  get_longname_units)


# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
            description=(f"{Yellow}MarsFiles is a file manager. Use it "
                 f"to modify a netCDF file format.{nocolor} \n\n"),
            formatter_class = argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs = '+',
            help = (f"A netCDF file or list of netCDF files. \n\n"))

parser.add_argument('-fv3', '--fv3', nargs = '+',
            help = (f"Produce MGCM 'diurn', 'average' and 'daily' "
            f"files from Legacy output. \n"
            f"Available options are: \n"
            f"  - 'fixed'  : static fields (e.g., topography) \n"
            f"  - 'average': 5-sol averages \n"
            f"  - 'daily'  : 5-sol continuous \n"
            f"  - 'diurn'  : 5-sol averages for each time of day \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py filename.nc -fv3 fixed \n"
            f"> MarsFiles.py filename.nc -fv3 fixed diurn "
            f"{nocolor}\n\n"))

parser.add_argument('-c', '--combine', action='store_true',
            help = (f"Combine sequential files of the same type into "
            f"one file. \n"
            f"Works with all file types ('fixed', 'average', "
            f"'daily' and 'diurn'). \n"
            f"{Yellow}Overwrites the first file in the series. "
            f"To override, use --ext.{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_average.nc --combine "
            f"{nocolor}\n\n"))

parser.add_argument('-t', '--tshift', nargs = '?', const = 999, 
            type = str,
            help = (f"Apply a time-shift to {Yellow}'diurn'{nocolor} "
            f"files. \n"
            "Vertically interpolated 'diurn' files OK. \n"
            f"{Yellow}Generates a new file ending in '_T.nc'{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_diurn.nc --tshift \n"
            f"  {Blue}(outputs data for all 24 local times){Green} \n"
            f"> MarsFiles.py *.atmos_diurn.nc --tshift '3 15'"
            f"\n"
            f"  {Blue}(outputs data for target local times only)"
            f"{nocolor}\n\n"))

parser.add_argument('-ba', '--bin_average', nargs = '?', const = 5, 
            type = int,
            help = (f"Bin MGCM 'daily' files like 'average' files. \n"
            f"{Yellow}Generates a new file ending in '_to_average.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -ba \n"
            f"  {Blue}(nocolor, bin 5 days){Green} \n"
            f"> MarsFiles.py *.atmos_daily_pstd.nc -ba 10 \n"
            f"  {Blue}(bin 10 days)"
            f"{nocolor}\n\n"))

parser.add_argument('-bd', '--bin_diurn', action = 'store_true',
            help = (f"Bin MGCM 'daily' files like 'diurn' files. \n"
            f"May be used jointly with --bin_average. \n"
            f"{Yellow}Generates a new file ending in '_to_diurn.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -bd \n"
            f"  {Blue}(default 5-day bin){Green} \n"
            f"> MarsFiles.py *.atmos_daily_pstd.nc -bd -ba 10 \n"
            f"  {Blue}(10-day bin){Green} \n"
            f"> MarsFiles.py *.atmos_daily_pstd.nc -bd -ba 1 \n"
            f"  {Blue}(No binning. Mimics raw Legacy output)"
            f"{nocolor}\n\n"))


parser.add_argument('-hpf', '--high_pass_filter', nargs = '+', 
            type = float,
            help = (f"Temporal filtering utilities: low-, high-, and "
            f"band-pass filters. \n"
            f"Use '--no_trend' to compute amplitudes only. \n"
            f"Data detrended before filtering. \n"
            f"{Yellow}Generates a new file ending in '_hpf.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -hpf 10. \n"
            f"  {Blue}(-hpf) --high_pass_filter sol_min "
            f"{nocolor}\n\n"))

parser.add_argument('-lpf', '--low_pass_filter', nargs = '+', 
            type = float,
            help = (f"Temporal filtering utilities: low-, high-, and "
            f"band-pass filters. \n"
            f"Use '--no_trend' to compute amplitudes only. \n"
            f"Data detrended before filtering. \n"
            f"{Yellow}Generates a new file ending in '_lpf.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -lpf 0.5 \n"
            f"  {Blue}(-lpf) --low_pass_filter sol_max "
            f"{nocolor}\n\n"))

parser.add_argument('-bpf', '--band_pass_filter', nargs = '+',
            help = (f"Temporal filtering utilities: low-, high-, and "
            f"band-pass filters. \n"
            f"Use '--no_trend' to compute amplitudes only. \n"
            f"Data detrended before filtering. \n"
            f"{Yellow}Generates a new file ending in 'bpf.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -hpf 0.5 10. \n"
            f"  {Blue}(-bpf) --band_pass_filter sol_min sol max "
            f"{nocolor}\n\n"))

parser.add_argument('-no_trend', '--no_trend', action='store_true',
            help = (f"Filter and compute amplitudes only. \n"
            f"For use with temporal filtering utilities (-lpf, -hpf, "
            f"-bpf). \n"
            f"{Yellow}Generates a new file ending in '_no_trend.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -hpf 10. "
            f"--no_trend \n"
            f"> MarsFiles.py *.atmos_daily.nc -lpf 0.5 "
            f"--no_trend \n"
            f"> MarsFiles.py *.atmos_daily.nc -hpf 0.5 10. "
            f"--no_trend "
            f"{nocolor}\n\n"))

# Decomposition in zonal harmonics, disabled for initial CAP release:
# parser.add_argument('-hpk', '--high_pass_zonal', nargs = '+', 
#             type = int,
#             help = (f"Spatial filtering utilities: low-, high-, and "
#             f"band pass filters. \n"
#             f"Use '--no_trend' to compute amplitudes only. \n"
#             f"Data detrended before filtering. \n"
#             f"{Yellow}Generates a new file ending in '_hpk.nc'"
#             f"{nocolor} \n"
#             f"{Green}Usage: \n"
#             f"    > MarsFiles.py *.atmos_daily.nc -hpk 10 "
#             f"--no_trend \n"
#             f"      {Blue}(-hpk)  --high_pass_zonal kmin "
#             f"{nocolor}\n\n"))

# parser.add_argument('-lpk', '--low_pass_zonal', nargs = '+', type = int,
#             help = (f"Spatial filtering utilities: low-, high-, and "
#             f"band pass filters. \n"
#             f"Use '--no_trend' to compute amplitudes only. \n"
#             f"Data detrended before filtering. \n"
#             f"{Yellow}Generates a new file ending in '_lpk.nc'"
#             f"{nocolor} \n"
#             f"{Green}Usage: \n"
#             f"    > MarsFiles.py *.atmos_daily.nc -lpk 20 "
#             f"--no_trend\n"
#             f"      {Blue}(-lpk)  --low_pass_zonal  kmax "
#             f"{nocolor}\n\n"))

# parser.add_argument('-bpk', '--band_pass_zonal', nargs = '+',
#             help = (f"Spatial filtering utilities: low-, high-, and "
#             f"band pass filters. \n"
#             f"Use '--no_trend' to compute amplitudes only. \n"
#             f"Data detrended before filtering. \n"
#             f"{Yellow}Generates a new file ending in '_bpk.nc'"
#             f"{nocolor} \n"
#             f"{Green}Usage: \n"
#             f"    > MarsFiles.py *.atmos_daily.nc -bpk 10 20 "
#             f"--no_trend\n"
#             f"      {Blue}(-bpk)  --band_pass_zonal kmin kmax "
#             f"{nocolor}\n\n"))

parser.add_argument('-tidal', '--tidal', nargs = '+', type = int,
            help = (f"Performs a tidal analyis on 'diurn' files. \n"
            f"Extracts diurnal tide and its harmonics. \n"
            f"N = 1 diurnal, N = 2 semi-diurnal etc. \n"
            f"{Yellow}Generates a new file ending in '_tidal.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_diurn.nc -tidal 4 \n"
            f"  {Blue}(extracts 4 harmonics) "
            f"{nocolor}\n\n"))

parser.add_argument('-reconstruct', '--reconstruct', 
            action = 'store_true',
            help = (f"Reconstructs the first N harmonics. \n"
            f"{Yellow}Generates a new file ending in '_reconstruct.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_diurn.nc -tidal 6 "
            f"--include ps temp --reconstruct "
            f"{nocolor}\n\n"))

parser.add_argument('-norm', '--normalize', action = 'store_true',
            help = (f"Provides result in percent amplitude. \n"
            f"{Yellow}Generates a new file ending in '_norm.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_diurn.nc -tidal 6 "
            f"--include ps --normalize "
            f"{nocolor}\n\n"))

parser.add_argument('-rs', '--regrid_source', nargs = '+',
            help = (f"Regrid a target file to match a source file. \n"
            f"Both source and target files should be vertically \n"
            f"interpolated to the same standard grid \n"
            f"(e.g. zstd, zagl, pstd, etc.). \n"
            f"{Yellow}Generates a new file ending in '_regrid.nc'"
            f"{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsInterp.py *.atmos.average_pstd.nc -rs "
            f"simu2/00668.atmos_average_pstd.nc "
            f"{nocolor}\n\n"))

parser.add_argument('-za', '--zonal_avg', action = 'store_true',
            help = (f"Zonally average all variables in a file. \n"
            f"{Yellow}Generates a new file ending in '_zonal_avg.nc'"
            f"{nocolor} \n"
            "   > MarsFiles.py *.atmos_diurn.nc -za "
            f"{nocolor}\n\n"))

parser.add_argument('-include', '--include', nargs = '+',
            help = (f"Flag to include only the variables listed \n"
            f"after -include in the target file. \n"
            f"All dimensional and 1D variables are always included. \n"
            f"{Yellow}Overwrites existing target file. To override, "
            f"use --ext.{nocolor} \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos_daily.nc -ba --include ps "
            f"ts ucomp "
            f"{nocolor}\n\n"))

parser.add_argument('-e', '--ext', type = str, default = None,
            help = (f"Do not overwrite file. Append the extension \n"
            f"provided after --ext to the new file. \n"
            f"{Green}Usage: \n"
            f"> MarsFiles.py *.atmos.average.nc --combine "
            f"--ext _combined \n"
            f"  {Blue}(produces *.atmos.average_combined.nc) "
            f"{nocolor}\n\n"))

parser.add_argument('--debug', action='store_true',
            help = (f"Debug flag: release the exceptions. \n\n"))

# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
    file_list = parser.parse_args().input_file
    cwd       = os.getcwd()
    path2data = os.getcwd()

    if parser.parse_args().fv3 and parser.parse_args().combine:
        prRed("Use --fv3 and --combine sequentially to avoid ambiguity")
        exit()

    # "===========================================================================
    # ==========  Conversion Legacy -> FV3 by Richard U. and Alex. K. ===========
    # ===========================================================================

    # ======= Convert to MGCM Output Format =======
    if parser.parse_args().fv3:
        for irequest in parser.parse_args().fv3:
            if irequest not in ['fixed', 'average', 'daily', 'diurn']:
                prRed(f"{irequest} is not available, select 'fixed', 'average', 'daily', or 'diurn'")

    # Argument Definitions:
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        # Get files to process
        histlist = []
        for filei in file_list:
            if not ('/' in filei):
                histlist.append(path2data+'/'+filei)
            else:
                histlist.append(filei)
        fnum = len(histlist)

        lsmin = None
        lsmax = None

        if histlist[0][-3:] == '.nc':
            print("Processing LegacyGCM_*.nc files")
            for f in histlist:
                histname = os.path.basename(f)
                ls_l = histname[-12:-9]
                ls_r = histname[-6:-3]
                if lsmin is None:
                    lsmin = ls_l
                else:
                    lsmin = str(min(int(lsmin), int(ls_l))).zfill(3)
                if lsmax is None:
                    lsmax = ls_r
                else:
                    lsmax = str(max(int(lsmax), int(ls_r))).zfill(3)
                a = make_FV3_files(f, parser.parse_args().fv3, True, cwd)

        else:
            print("Processing fort.11 files")
            for fname in histlist:
                f = Fort(fname)
                if 'fixed' in parser.parse_args().fv3:
                    f.write_to_fixed()
                if 'average' in parser.parse_args().fv3:
                    f.write_to_average()
                if 'daily' in parser.parse_args().fv3:
                    f.write_to_daily()
                if 'diurn' in parser.parse_args().fv3:
                    f.write_to_diurn()

    # ===========================================================================
    # =============  Append netcdf files along the 'time' dimension =============
    # ===========================================================================
    elif parser.parse_args().combine:
        prYellow("Using internal method for concatenation")

        # Get files to process
        histlist = []
        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                histlist.append(path2data+'/'+filei)
            else:
                histlist.append(filei)

        fnum = len(histlist)
        # Easy case: merging **.fixed.nc means deleting all but the first file:
        if file_list[0][5:] == '.fixed.nc' and fnum >= 2:
            rm_cmd = 'rm -f '
            for i in range(1, fnum):
                rm_cmd += ' '+histlist[i]
            p = subprocess.run(rm_cmd, universal_newlines=True, shell=True)
            prCyan(f"Cleaned all but {file_list[0]}")
            exit()

        # =========
        fnum = len(histlist)
        prCyan(f"Merging {fnum} files, starting with {file_list[0]} ...")

        # This section iexcludes any variable not listed after --include
        if parser.parse_args().include:
            f = Dataset(file_list[0], 'r')
            exclude_list = filter_vars(f, parser.parse_args(
            ).include, giveExclude=True)  # variable to exclude
            f.close()
        else:
            exclude_list = []

        # This creates a temporaty file ***_tmp.nc to work in
        file_tmp = histlist[0][:-3]+'_tmp'+'.nc'
        Log = Ncdf(file_tmp, 'Merged file')
        Log.merge_files_from_list(histlist, exclude_var=exclude_list)
        Log.close()

        # ===== Delete the files that were combined ====

        # Rename merged file LegacyGCM_LsINI_LsEND.nc or first files of the list (e.g 00010.atmos_average.nc)
        if file_list[0][:12] == 'LegacyGCM_Ls':
            ls_ini = file_list[0][12:15]
            ls_end = file_list[-1][18:21]
            fileout = 'LegacyGCM_Ls%s_Ls%s.nc' % (ls_ini, ls_end)
        else:
            fileout = histlist[0]

        # Assemble 'remove' and 'move' commands to execute
        rm_cmd = 'rm -f '
        for ifile in histlist:
            rm_cmd += ' '+ifile
        cmd_txt = 'mv '+file_tmp+' '+fileout
        p = subprocess.run(rm_cmd, universal_newlines=True, shell=True)
        p = subprocess.run(cmd_txt, universal_newlines=True, shell=True)
        prCyan(f"{fileout} was merged")

# ===============================================================================
# ============= Time-Shifting Implementation by Victoria H. =====================
# ===============================================================================

    elif parser.parse_args().tshift:
        # target_list holds the target local times
        if parser.parse_args().tshift == 999:
            target_list = None
        else:
            target_list = np.fromstring(
                parser.parse_args().tshift, dtype=float, sep=' ')

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+'_T'+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdiurn = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy some dimensions from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdiurn)

            # Find the "time of day" variable name
            tod_name_in = find_tod_in_diurn(fdiurn)
            _, zaxis = FV3_file_type(fdiurn)

            # Copy some variables from the old file to the new file
            fnew.copy_Ncaxis_with_content(fdiurn.variables['lon'])
            fnew.copy_Ncaxis_with_content(fdiurn.variables['lat'])
            fnew.copy_Ncaxis_with_content(fdiurn.variables['time'])
            fnew.copy_Ncaxis_with_content(fdiurn.variables['scalar_axis'])

            # Only create a vertical axis if the original file contains 3D fields
            if zaxis in fdiurn.dimensions.keys():
                fnew.copy_Ncaxis_with_content(fdiurn.variables[zaxis])

            # Copy some dimensions from the old file to the new file
            if target_list is None:
                # Same input local times are used as target local times, use the old axis as-is
                tod_orig = np.array(fdiurn.variables[tod_name_in])
                tod_name_out = tod_name_in
                fnew.copy_Ncaxis_with_content(fdiurn.variables[tod_name_in])
                # tod_in=np.array(fdiurn.variables[tod_name_in])
                tod_in = None
                # Only copy 'areo' if it exists in the original file
                if 'areo' in fdiurn.variables.keys():
                    fnew.copy_Ncvar(fdiurn.variables['areo'])
            else:

                tod_orig = np.array(fdiurn.variables[tod_name_in])
                # Copy all dimensions but time_of_day. Update time_of_day array.
                # fnew.copy_all_dims_from_Ncfile(fdiurn,exclude_dim=tod_name_in)
                tod_in = target_list
                tod_name_out = 'time_of_day_%02i' % (len(tod_in))
                fnew.add_dim_with_content(tod_name_out, tod_in, longname_txt="time of day",
                                          units_txt='[hours since 0000-00-00 00:00:00]', cart_txt='')

                # Create 'areo' variable with the new size
                areo_in = fdiurn.variables['areo'][:]
                areo_shape = areo_in.shape
                dims_out = fdiurn.variables['areo'].dimensions

                # Update shape with new time_of_day
                areo_shape = (areo_shape[0], len(tod_in), areo_shape[2])
                dims_out = (dims_out[0], tod_name_out, dims_out[2])
                areo_out = np.zeros(areo_shape)
                # For new tod_in, e.g [3,15]
                for ii in range(len(tod_in)):
                    # Get the closest 'time_of_day' index in the input array
                    it = np.argmin(np.abs(tod_in[ii]-tod_orig))
                    areo_out[:, ii, 0] = areo_in[:, it, 0]

                fnew.add_dim_with_content(
                    'scalar_axis', [0], longname_txt="none", units_txt='none')
                fnew.log_variable('areo', areo_out, dims_out,
                                  'areo', 'degrees')

            # Read 4D field and do the time shift
            longitude = np.array(fdiurn.variables['lon'])
            var_list = filter_vars(
                fdiurn, parser.parse_args().include)  # Get all variables

            for ivar in var_list:
                prCyan(f"Processing: {ivar} ...")
                varNcf = fdiurn.variables[ivar]
                varIN = varNcf[:]
                vkeys = varNcf.dimensions
                longname_txt, units_txt = get_longname_units(fdiurn, ivar)
                if (len(vkeys) == 4):
                    ilat = vkeys.index('lat')
                    ilon = vkeys.index('lon')
                    itime = vkeys.index('time')
                    itod = vkeys.index(tod_name_in)
                    newvar = np.transpose(varIN, (ilon, ilat, itime, itod))
                    newvarOUT = tshift(newvar, longitude,
                                       tod_orig, timex=tod_in)
                    varOUT = np.transpose(newvarOUT, (2, 3, 1, 0))
                    fnew.log_variable(
                        ivar, varOUT, ['time', tod_name_out, 'lat', 'lon'], longname_txt, units_txt)
                if (len(vkeys) == 5):
                    ilat = vkeys.index('lat')
                    ilon = vkeys.index('lon')
                    iz = vkeys.index(zaxis)
                    itime = vkeys.index('time')
                    itod = vkeys.index(tod_name_in)
                    newvar = np.transpose(varIN, (ilon, ilat, iz, itime, itod))
                    newvarOUT = tshift(newvar, longitude,
                                       tod_orig, timex=tod_in)
                    varOUT = np.transpose(newvarOUT, (3, 4, 2, 1, 0))
                    fnew.log_variable(ivar, varOUT, [
                                      'time', tod_name_out, zaxis, 'lat', 'lon'], longname_txt, units_txt)
            fnew.close()
            fdiurn.close()

    # ===========================================================================
    # ===============  Bin a 'daily' file to an 'average' file ==================
    # ===========================================================================
    elif parser.parse_args().bin_average and not parser.parse_args().bin_diurn:
        nday = parser.parse_args().bin_average
        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+'_to_average'+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdaily = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
            var_list = filter_vars(
                fdaily, parser.parse_args().include)  # Get all variables

            time_in = fdaily.variables['time'][:]
            Nin = len(time_in)

            dt_in = time_in[1]-time_in[0]
            iperday = int(np.round(1/dt_in))
            combinedN = int(iperday*nday)

            N_even = Nin//combinedN
            N_left = Nin % combinedN

            if N_left != 0:
                prYellow(f"***Warning*** requested {nday} sols bin period. File has {iperday} timestep/sols and {Nin}/({nday} x {iperday}) is not a round number")
                prYellow(f"    Will use {N_even}  bins of ({nday} x {iperday})={combinedN} timesteps ({N_even*combinedN}) and discard {N_left} timesteps")

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy all dimensions but 'time' from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim=['time'])

            # Calculate and log the new time array
            fnew.add_dimension('time', None)
            time_out = daily_to_average(time_in[:], dt_in, nday)
            fnew.log_axis1D('time', time_out, 'time', longname_txt="sol number",
                            units_txt='days since 0000-00-00 00:00:00', cart_txt='T')

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]

                if 'time' in varNcf.dimensions:
                    prCyan(f"Processing: {ivar}")
                    var_out = daily_to_average(varNcf[:], dt_in, nday)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(
                        ivar, var_out, varNcf.dimensions, longname_txt, units_txt)

                else:
                    if ivar in ['pfull', 'lat', 'lon', 'phalf', 'pk', 'bk', 'pstd', 'zstd', 'zagl']:
                        prCyan(f"Copying axis: {ivar}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        prCyan(f"Copying variable: {ivar}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ===========================================================================
    # ===============  Bin a 'daily' file to a 'diurn' file =====================
    # ===========================================================================
    elif parser.parse_args().bin_diurn:
        # Use defaut binning period of 5 days (like 'average' files)
        if parser.parse_args().bin_average is None:
            nday = 5
        else:
            nday = parser.parse_args().bin_average

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+'_to_diurn'+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdaily = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
            var_list = filter_vars(
                fdaily, parser.parse_args().include)  # Get all variables

            time_in = fdaily.variables['time'][:]
            Nin = len(time_in)

            dt_in = time_in[1]-time_in[0]
            iperday = int(np.round(1/dt_in))

            # define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy all dimensions but 'time' from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim=['time'])

            # If no binning is requested, copy time axis as-is
            fnew.add_dimension('time', None)
            time_out = daily_to_average(time_in[:], dt_in, nday)
            fnew.add_dim_with_content('time', time_out, longname_txt="sol number",
                                      units_txt='days since 0000-00-00 00:00:00', cart_txt='T')

            # Create a new 'time_of_day' dimension
            tod_name = 'time_of_day_%02d' % (iperday)
            time_tod = np.squeeze(daily_to_diurn(
                time_in[0:iperday], time_in[0:iperday]))
            tod = np.mod(time_tod*24, 24)
            fnew.add_dim_with_content(tod_name, tod, longname_txt="time of day",
                                      units_txt="hours since 0000-00-00 00:00:00", cart_txt='N')

            # Loop over all variables in the file
            for ivar in var_list:

                varNcf = fdaily.variables[ivar]

                # If 'time' is the dimension (not just a 'time' array)
                if 'time' in varNcf.dimensions and ivar != 'time':
                    prCyan(f"Processing: {ivar}")
                    dims_in = varNcf.dimensions
                    dims_out = (dims_in[0],)+(tod_name,)+dims_in[1:]
                    var_out = daily_to_diurn(varNcf[:], time_in[0:iperday])
                    if nday != 1:
                        # dt is 1 sol between two 'diurn' timesteps
                        var_out = daily_to_average(var_out, 1., nday)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(ivar, var_out, dims_out,
                                      longname_txt, units_txt)

                else:
                    if ivar in ['pfull', 'lat', 'lon', 'phalf', 'pk', 'bk', 'pstd', 'zstd', 'zagl']:
                        prCyan(f"Copying axis: {ivar}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    elif ivar != 'time':
                        prCyan(f"Copying variable: {ivar}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ===========================================================================
    # ========================  Transient wave analysis =========================
    # ===========================================================================

    elif parser.parse_args().high_pass_filter or parser.parse_args().low_pass_filter or parser.parse_args().band_pass_filter:

        # This functions requires scipy > 1.2.0. We import the package here.
        from amescap.Spectral_utils import zeroPhi_filter

        if parser.parse_args().high_pass_filter:
            btype = 'high'
            out_ext = '_hpf'
            nsol = np.asarray(
                parser.parse_args().high_pass_filter).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                prRed("***Error*** sol_min accepts only one value")
                exit()
        if parser.parse_args().low_pass_filter:
            btype = 'low'
            out_ext = '_lpf'
            nsol = np.asarray(
                parser.parse_args().low_pass_filter).astype(float)
            if len(np.atleast_1d(nsol)) != 1:
                prRed("sol_max accepts only one value")
                exit()
        if parser.parse_args().band_pass_filter:
            btype = 'band'
            out_ext = '_bpf'
            nsol = np.asarray(
                parser.parse_args().band_pass_filter).astype(float)
            if len(np.atleast_1d(nsol)) != 2:
                prRed("Requires two values: sol_min sol_max")
                exit()
        if parser.parse_args().no_trend:
            out_ext = out_ext+'_no_trend'

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+out_ext+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdaily = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')

            var_list = filter_vars(
                fdaily, parser.parse_args().include)  # Get all variables

            time_in = fdaily.variables['time'][:]

            dt = time_in[1]-time_in[0]

            # Check if the frequency domain is allowed
            if any(nn <= 2*dt for nn in nsol):
                prRed(f"***Error***  minimum cut-off cannot be smaller than the Nyquist period of 2xdt={2*dt} sol")
                exit()

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy all dimensions but 'time' from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily)

            if btype == 'low':
                fnew.add_constant(
                    'sol_max', nsol, "Low-pass filter cut-off period ", "sol")
            elif btype == 'high':
                fnew.add_constant(
                    'sol_min', nsol, "High-pass filter cut-off period ", "sol")
            elif btype == 'band':
                fnew.add_constant(
                    'sol_min', nsol[0], "High-pass filter low cut-off period ", "sol")
                fnew.add_constant(
                    'sol_max', nsol[1], "High-pass filter high cut-off period ", "sol")
            dt = time_in[1]-time_in[0]

            fs = 1/(dt)  # Frequency in sol-1
            if btype == 'band':
                # Flip the sols so that the low frequency comes first
                low_highcut = 1/nsol[::-1]
            else:
                low_highcut = 1./nsol

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdaily.variables[ivar]

                if 'time' in varNcf.dimensions and ivar not in ['time', 'areo']:
                    prCyan(f"Processing: {ivar}")
                    var_out = zeroPhi_filter(
                        varNcf[:], btype, low_highcut, fs, axis=0, order=4, no_trend=parser.parse_args().no_trend)
                    longname_txt, units_txt = get_longname_units(fdaily, ivar)
                    fnew.log_variable(
                        ivar, var_out, varNcf.dimensions, longname_txt, units_txt)
                else:
                    if ivar in ['pfull', 'lat', 'lon', 'phalf', 'pk', 'bk', 'pstd', 'zstd', 'zagl']:
                        prCyan(f"Copying axis: {ivar}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    else:
                        prCyan(f"Copying variable: {ivar}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()

    # ===========================================================================
    # ========================  Zonal Decomposition Analysis ====================
    # ===========================================================================

    # elif parser.parse_args().high_pass_zonal or parser.parse_args().low_pass_zonal or parser.parse_args().band_pass_zonal:
    #
    #     # This function requires scipy > 1.2.0. We import the package here
    #     from amescap.Spectral_utils import zonal_decomposition, zonal_construct
    #     #Load the module
    #     #init_shtools()
    #
    #     if parser.parse_args().high_pass_zonal:
    #         btype='high';out_ext='_hpk';nk=np.asarray(parser.parse_args().high_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk))!=1:
    #             prRed("***Error*** kmin accepts only one value")
    #             exit()
    #     if parser.parse_args().low_pass_zonal:
    #         btype='low';out_ext='_lpk';nk=np.asarray(parser.parse_args().low_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk))!=1:
    #             prRed("kmax accepts only one value")
    #             exit()
    #     if parser.parse_args().band_pass_zonal:
    #         btype='band';out_ext='_bpk';nk=np.asarray(parser.parse_args().band_pass_zonal).astype(int)
    #         if len(np.atleast_1d(nk))!=2:
    #             prRed("Requires two values: kmin kmax")
    #             exit()
    #
    #     if parser.parse_args().no_trend:out_ext =f"{out_ext}_no_trend"
    #
    #     for filei in file_list:
    #         # Add path unless full path is provided
    #         if not ('/' in filei):
    #             fullnameIN = path2data + '/' + filei
    #         else:
    #             fullnameIN=filei
    #         fullnameOUT = fullnameIN[:-3]+out_ext+'.nc'
    #
    #         # Append extension, if any:
    #         if parser.parse_args().ext:fullnameOUT=fullnameOUT[:-3]+'_'+parser.parse_args().ext+'.nc'
    #
    #         fname = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
    #
    #         var_list = filter_vars(fname,parser.parse_args().include) # Get all variables
    #
    #         lon=fname.variables['lon'][:]
    #         lat=fname.variables['lat'][:]
    #         LON,LAT=np.meshgrid(lon,lat)
    #
    #         dlat=lat[1]-lat[0]
    #         dx=2*np.pi*3400
    #
    #         # Check if the frequency domain is allowed and display some information
    #
    #         if any(nn > len(lat)/2 for nn in nk):
    #             prRed(f"***Warning***  maximum wavenumber cut-off cannot be larger than the Nyquist criteria of nlat/2= {len(lat)/2)} sol")
    #         elif btype=='low':
    #             L_max=(1./nk)*dx
    #             prYellow("Low pass filter, allowing only wavelength > {L_max} km")
    #         elif btype=='high':
    #             L_min=(1./nk)*dx
    #             prYellow("High pass filter, allowing only wavelength < {L_min} km")
    #         elif btype=='band':
    #             L_min=(1./nk[1])*dx
    #             L_max=1./max(nk[0],1.e-20)*dx
    #             if L_max>1.e20:L_max=np.inf
    #             prYellow("Band pass filter, allowing only {L_min} km < wavelength < {L_max} km")

   ##
    #         fnew = Ncdf(fullnameOUT) # Define a netcdf object from the netcdf wrapper module
    #         # Copy all dimensions but 'time' from the old file to the new file
    #         fnew.copy_all_dims_from_Ncfile(fname)
    #
    #         if btype=='low':
    #             fnew.add_constant('kmax',nk,"Low-pass filter zonal wavenumber ","wavenumber")
    #         elif btype=='high':
    #             fnew.add_constant('kmin',nk,"High-pass filter zonal wavenumber ","wavenumber")
    #         elif btype=='band':
    #             fnew.add_constant('kmin',nk[0],"Band-pass filter low zonal wavenumber ","wavenumber")
    #             fnew.add_constant('kmax',nk[1],"Band-pass filter high zonal wavenumber ","wavenumber")
    #
    #         low_highcut = nk
    #
    #         #Loop over all variables in the file
    #         for ivar in var_list:
    #             varNcf = fname.variables[ivar]
    #
    #             if ('lat' in varNcf.dimensions) and ('lon' in varNcf.dimensions):
    #                 prCyan(f"Processing: {ivar} ...")
    #
    #                 # Step 1 : Detrend the data
    #                 TREND=get_trend_2D(varNcf[:],LON,LAT,'wmean')
    #                 # Step 2 : Calculate spherical harmonic coefficients
    #                 COEFF,PSD=zonal_decomposition(varNcf[:]-TREND)
    #                 # Step 3 : Recompose the variable out of the coefficients
    #                 VAR_filtered=zonal_construct(COEFF,varNcf[:].shape,btype=btype,low_highcut=low_highcut)
    #                 #Step 4: Add the trend, if requested
    #                 if parser.parse_args().no_trend:
    #                     var_out=VAR_filtered
    #                 else:
    #                     var_out=VAR_filtered+TREND
    #
    #                 fnew.log_variable(ivar,var_out,varNcf.dimensions,varNcf.long_name,varNcf.units)
    #             else:
    #                 if  ivar in ['pfull', 'lat', 'lon','phalf','pk','bk','pstd','zstd','zagl','time']:
    #                     prCyan(f"Copying axis: {ivar} ...")
    #                     fnew.copy_Ncaxis_with_content(fname.variables[ivar])
    #                 else:
    #                     prCyan(f"Copying variable: {ivar} ...")
    #                     fnew.copy_Ncvar(fname.variables[ivar])
    #         fnew.close()

    # ===========================================================================
    # ============================  Tidal Analysis ==============================
    # ===========================================================================

    elif parser.parse_args().tidal:
        from amescap.Spectral_utils import diurn_extract, reconstruct_diurn
        N = parser.parse_args().tidal[0]
        if len(np.atleast_1d(N)) != 1:
            prRed("***Error*** N accepts only one value")
            exit()
        out_ext = '_tidal'
        if parser.parse_args().reconstruct:
            out_ext = out_ext+'_reconstruct'
        if parser.parse_args().normalize:
            out_ext = out_ext+'_norm'

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+out_ext+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdiurn = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')

            var_list = filter_vars(
                fdiurn, parser.parse_args().include)  # Get all variables

            # Find 'time_of_day' variable name
            tod_name = find_tod_in_diurn(fdiurn)

            tod_in = fdiurn.variables[tod_name][:]
            lon = fdiurn.variables['lon'][:]
            areo = fdiurn.variables['areo'][:]

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy all dims but 'time_of_day' from the old file to the new file

            # Harmonics to reconstruct the signal. We use the original time_of_day array.
            if parser.parse_args().reconstruct:
                fnew.copy_all_dims_from_Ncfile(fdiurn)
                # Copy time_of_day axis from initial files
                fnew.copy_Ncaxis_with_content(fdiurn.variables[tod_name])

            else:
                fnew.copy_all_dims_from_Ncfile(fdiurn, exclude_dim=[tod_name])
                # Create new dimension holding the harmonics. We reuse the 'time_of_day' name to facilitate
                # Compatible with other routines, but keep in mind this is the harmonic number
                fnew.add_dim_with_content('time_of_day_%i' % (N), np.arange(
                    1, N+1), longname_txt="tidal harmonics", units_txt="Diurnal harmonic number", cart_txt='N')

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf = fdiurn.variables[ivar]
                varIN = varNcf[:]
                longname_txt, units_txt = get_longname_units(fdiurn, ivar)
                var_unit = getattr(varNcf, 'units', '')

                if tod_name in varNcf.dimensions and ivar not in [tod_name, 'areo'] and len(varNcf.shape) > 2:
                    prCyan(f"Processing: {ivar}")

                    # Normalize the data
                    if parser.parse_args().normalize:
                        # Normalize and reshape the array along the time_of_day dimension
                        norm = np.mean(varIN, axis=1)[:, np.newaxis, ...]
                        varIN = 100*(varIN-norm)/norm
                        #units_txt='% of diurnal mean'
                        var_unit = '% of diurnal mean'

                    amp, phas = diurn_extract(
                        varIN.swapaxes(0, 1), N, tod_in, lon)
                    if parser.parse_args().reconstruct:
                        VARN = reconstruct_diurn(
                            amp, phas, tod_in, lon, sumList=[])
                        for nn in range(N):
                            fnew.log_variable(f"{ivar}_N{nn+1}", VARN[nn, ...].swapaxes(
                                0, 1), varNcf.dimensions, f"harmonic N={nn+1} for {longname_txt}", units_txt)

                    else:
                        #Update the dimensions
                        new_dim=list(varNcf.dimensions)
                        new_dim[1]=f"time_of_day_{N}"
                        fnew.log_variable(f"{ivar}_amp",amp.swapaxes(0,1),new_dim,f"tidal amplitude for {longname_txt}",units_txt)
                        fnew.log_variable(f"{ivar}_phas",phas.swapaxes(0,1),new_dim,f"tidal phase for {longname_txt}","hr")

                elif  ivar in ['pfull', 'lat', 'lon','phalf','pk','bk','pstd','zstd','zagl','time']:
                        prCyan(f"Copying axis: {ivar}...")
                        fnew.copy_Ncaxis_with_content(fdiurn.variables[ivar])
                elif  ivar in ['areo']:
                        if parser.parse_args().reconstruct:
                            #time_of_day is the same size as the original file
                            prCyan(f"Copying axis: {ivar}...")
                            fnew.copy_Ncvar(fdiurn.variables['areo'])
                        else:
                            prCyan(f"Processing: {ivar} ...")
                            #Create areo variable reflecting the new shape
                            areo_new=np.zeros((areo.shape[0],N,1))
                            #Copy areo
                            for xx in range(N):areo_new[:,xx,:]=areo[:,0,:]
                            #Update the dimensions
                            new_dim=list(varNcf.dimensions)
                            new_dim[1]=f"time_of_day_{N}"
                            #fnew.log_variable(ivar,areo_new,new_dim,longname_txt,units_txt)
                            fnew.log_variable(ivar,areo_new,new_dim,longname_txt,var_unit)

            fnew.close()

    # ===========================================================================
    # =============================  Regrid  files ==============================
    # ===========================================================================

    elif parser.parse_args().regrid_source:
        out_ext = '_regrid'
        name_target = parser.parse_args().regrid_source[0]

        # Add path unless full path is provided
        if not ('/' in name_target):
            name_target = path2data + '/' + name_target
        fNcdf_t = Dataset(name_target, 'r')

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+out_ext+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            f_in = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')

            var_list = filter_vars(
                f_in, parser.parse_args().include)  # Get all variables

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)

            # Copy all dims from the target file to the new file
            fnew.copy_all_dims_from_Ncfile(fNcdf_t)

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf     = f_in.variables[ivar]
                longname_txt,units_txt=get_longname_units(f_in,ivar)

                if  ivar in ['pfull', 'lat', 'lon','phalf','pk','bk','pstd','zstd','zagl','time','areo']:
                        prCyan(f"Copying axis: {ivar}...")
                        fnew.copy_Ncaxis_with_content(fNcdf_t.variables[ivar])
                elif varNcf.dimensions[-2:]==('lat', 'lon'): #Ignore variables like  'time_bounds', 'scalar_axis' or 'grid_xt_bnds'...
                    prCyan(f"Regridding: {ivar} ...")
                    var_OUT=regrid_Ncfile(varNcf,f_in,fNcdf_t)
                    fnew.log_variable(ivar,var_OUT,varNcf.dimensions,longname_txt,units_txt)

            fnew.close()
            fNcdf_t.close()

    # ===========================================================================
    # =======================  Zonal averaging    ===============================
    # ===========================================================================

    elif parser.parse_args().zonal_avg:

        for filei in file_list:
            # Add path unless full path is provided
            if not ('/' in filei):
                fullnameIN = path2data + '/' + filei
            else:
                fullnameIN = filei
            fullnameOUT = fullnameIN[:-3]+'_zonal_avg'+'.nc'

            # Append extension, if any:
            if parser.parse_args().ext:
                fullnameOUT = fullnameOUT[:-3] + \
                    '_'+parser.parse_args().ext+'.nc'

            fdaily = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
            var_list = filter_vars(
                fdaily, parser.parse_args().include)  # Get all variables

            lon_in = fdaily.variables['lon'][:]

            # Define a netcdf object from the netcdf wrapper module
            fnew = Ncdf(fullnameOUT)
            # Copy all dimensions but 'time' from the old file to the new file
            fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim=['lon'])

            # Add a new dimension for the longitude, size = 1
            fnew.add_dim_with_content('lon', [lon_in.mean(
            )], longname_txt="longitude", units_txt="degrees_E", cart_txt='X')

            # Loop over all variables in the file
            for ivar in var_list:
                varNcf     = fdaily.variables[ivar]
                longname_txt,units_txt=get_longname_units(fdaily,ivar)
                if 'lon' in varNcf.dimensions and ivar not in ['lon','grid_xt_bnds','grid_yt_bnds']:
                    prCyan(f"Processing: {ivar} ...")
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=RuntimeWarning)
                        var_out=np.nanmean(varNcf[:],axis=-1)[...,np.newaxis]
                        fnew.log_variable(ivar,var_out,varNcf.dimensions,longname_txt,units_txt)
                else:
                    if ivar in ['pfull', 'lat', 'phalf', 'pk', 'bk', 'pstd', 'zstd', 'zagl']:
                        prCyan("Copying axis: {ivar}")
                        fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
                    elif ivar in ['grid_xt_bnds', 'grid_yt_bnds', 'lon']:
                        pass

                    else:
                        prCyan(f"Copying variable: {ivar}")
                        fnew.copy_Ncvar(fdaily.variables[ivar])
            fnew.close()
    else:
        prRed("Error: no action requested: use 'MarsFiles *nc --fv3 --combine, --tshift, --bin_average, --bin_diurn etc ...'")

# END of script

# *******************************************************************************
# ************ Definitions for the functions used in this script ****************
# *******************************************************************************


def make_FV3_files(fpath, typelistfv3, renameFV3=True, cwd=None):
    """
    Make MGCM-like 'average', 'daily', and 'diurn' files.
    
    Parameters
    ----------
    fpath : str
        Full path to the Legacy netcdf files
    typelistfv3 : list
        MGCM-like file type: 'average', 'daily', or 'diurn'
    renameFV3 : bool
        Rename the files from Legacy_LsXXX_LsYYY.nc to \
            XXXXX.atmos_average.nc following MGCM output conventions
    cwd : str
        Sets the current working directory
    
    Returns
    -------
    The MGCM-like files: XXXXX.atmos_average.nc, XXXXX.atmos_daily.nc, \
        XXXXX.atmos_diurn.nc
    """

    histname = os.path.basename(fpath)
    if cwd is None:
        histdir = os.path.dirname(fpath)
        print(f"using fpath {histdir}")
    else:
        histdir = cwd
        print(f"using cwd {histdir}")

    histfile = Dataset(fpath, 'r', format='NETCDF4_CLASSIC')
    histvars = histfile.variables.keys()
    histdims = histfile.dimensions.keys()

    # Convert the first Ls in file to a sol number
    if renameFV3:
        fdate = '%05i' % (ls2sol_1year(histfile.variables['ls'][0]))

    def proccess_file(newf, typefv3):
        """
            Summary line - fits on one line

            Extended description of function, what it outputs

            Parameters
            ----------
            newf : str
                Path to target file
            typefv3 : str
                Identifies type of file: 'average', 'daily', or 'diurn'

            Raises
            ------
            SomeError
                Error description if X

            Returns
            -------
            returned_val
                define returned_val
        """
        for dname in histdims:
            if dname == 'nlon':
                var = histfile.variables['longitude']
                npvar = var[:]
                newf.add_dim_with_content('lon', npvar, 'longitude', getattr(var, 'units'), 'X')
            elif dname == 'nlat':
                var = histfile.variables['latitude']
                npvar = var[:]
                newf.add_dim_with_content('lat', npvar, 'latitude', getattr(var, 'units'), 'Y')

            elif dname == 'time':
                newf.add_dimension('time', None)
            elif dname == 'ntod' and typefv3 == 'diurn':
                dim = histfile.dimensions[dname]
                newf.add_dimension('time_of_day_16', dim.size)
            elif dname == 'nlay':
                nlay = histfile.dimensions[dname]
                num = nlay.size
                nump = num+1
                pref = 7.01*100  # in Pa
                pk = np.zeros(nump)
                bk = np.zeros(nump)
                pfull = np.zeros(num)
                phalf = np.zeros(nump)

                sgm = histfile.variables['sgm']
                # [AK] changed pk[0]=.08 to pk[0]=.08/2, otherwise phalf[0] would be greater than phalf[1]
                pk[0] = 0.08/2
                # *** NOTE that pk in amesCAP/mars_data/Legacy.fixed.nc is also updated***
                for z in range(num):
                    bk[z+1] = sgm[2*z+2]
                phalf[:] = pk[:]+pref*bk[:]  # Output in Pa

                # DEPRECATED: pfull[:] = (phalf[1:]-phalf[:num])/(np.log(phalf[1:])-np.log(phalf[:num]))
                # First layer:
                if pk[0] == 0 and bk[0] == 0:
                    pfull[0] = 0.5*(phalf[0]+phalf[1])
                else:
                    pfull[0] = (phalf[1]-phalf[0]) / \
                        (np.log(phalf[1])-np.log(phalf[0]))
                # All other layers:
                pfull[1:] = (phalf[2:]-phalf[1:-1]) / \
                    (np.log(phalf[2:])-np.log(phalf[1:-1]))

                newf.add_dim_with_content(
                    'pfull', pfull, 'ref full pressure level', 'Pa')
                newf.add_dim_with_content(
                    'phalf', phalf, 'ref half pressure level', 'Pa')
                newf.log_axis1D(
                    'pk', pk, ('phalf'), longname_txt='pressure part of the hybrid coordinate', units_txt='Pa', cart_txt='')
                newf.log_axis1D(
                    'bk', bk, ('phalf'), longname_txt='sigma part of the hybrid coordinate', units_txt='Pa', cart_txt='')
            else:
                dim = histfile.dimensions[dname]
                newf.add_dimension(dname, dim.size)

        # ===========END function========

    if 'average' in typelistfv3:
        newfname_avg = fdate+'.atmos_average.nc'  # 5 sol average over 'time_of_day' and 'time'
        newfpath_avg = os.path.join(histdir, newfname_avg)
        newfavg = Ncdf(newfpath_avg)
        proccess_file(newfavg, 'average')
        do_avg_vars(histfile, newfavg, True, True)
        newfavg.close()

    if 'daily' in typelistfv3:
        # Daily snapshot of the output
        newfname_daily = fdate+'.atmos_daily.nc'
        newfpath_daily = os.path.join(histdir, newfname_daily)
        newfdaily = Ncdf(newfpath_daily)
        proccess_file(newfdaily, 'daily')
        do_avg_vars(histfile, newfdaily, False, False)
        newfdaily.close()

    if 'diurn' in typelistfv3:
        newfname_diurn = fdate+'.atmos_diurn.nc'  # 5 sol average over 'time' only
        newfpath_diurn = os.path.join(histdir, newfname_diurn)
        newfdiurn = Ncdf(newfpath_diurn)
        proccess_file(newfdiurn, 'diurn')
        do_avg_vars(histfile, newfdiurn, True, False)
        newfdiurn.close()

    if 'fixed' in typelistfv3:
        # Copy Legacy.fixed to current directory
        cmd_txt = f"{cp} {sys.prefix}/mars_data/Legacy.fixed.nc {fdate}.fixed.nc"
        p = subprocess.run(cmd_txt, universal_newlines=True, shell=True)
        print(f"{cwd}/{fdate}.fixed.nc was copied locally")


# Function to perform time averages over all fields
def do_avg_vars(histfile, newf, avgtime, avgtod, Nday=5):
    histvars = histfile.variables.keys()
    for vname in histvars:
        var = histfile.variables[vname]
        npvar = var[:]
        dims = var.dimensions
        ndims = npvar.ndim
        vshape = npvar.shape
        ntod = histfile.dimensions['ntod']

        # longname_txt, units_txt = get_longname_units(histfile, vname)
        longname_txt = getattr(histfile.variables[vname], 'long_name', '')

        # On some files, like the LegacyGCM_Ls*** on the NAS data portal, the attribute 'long_name' may be mispelled 'longname'
        if longname_txt == '':
            longname_txt = getattr(histfile.variables[vname], 'longname', '')

        units_txt = getattr(histfile.variables[vname], 'units', '')

        if avgtod:
            newdims = replace_dims(dims, True)
        elif avgtime:
            newdims = replace_dims(dims, False)
        else:
            newdims = replace_dims(dims, True)

        if 'time' in dims:
            tind = dims.index('time')
            tind_new = newdims.index('time')
            numt = histfile.dimensions['time'].size
        # TODO fix time !!
        # now do various time averaging and write to files
        if ndims == 1:
            if vname == 'ls':

                # first check if ls crosses over to a new year
                if not np.all(npvar[1:] >= npvar[:-1]):
                    year = 0.
                    for x in range(1, npvar.size):
                        if 350. < npvar[x-1] < 360. and npvar[x] < 10.:
                            year += 1.
                        npvar[x] += 360.*year

                # Create a 'time' array
                time0 = ls2sol_1year(npvar[0])+np.linspace(0, 10., len(npvar))

                if avgtime:
                    varnew = np.mean(npvar.reshape(-1, Nday), axis=1)
                    time0 = np.mean(time0.reshape(-1, Nday), axis=1)

                if not avgtime and not avgtod:  # i.e 'daily' file
                    # Solar longitude
                    ls_start = npvar[0]
                    ls_end = npvar[-1]
                    step = (ls_end-ls_start)/np.float32(((numt-1)*ntod.size))
                    varnew = np.arange(0, numt*ntod.size, dtype=np.float32)
                    varnew[:] = varnew[:]*step+ls_start

                    # Time
                    step = (ls2sol_1year(ls_end)-ls2sol_1year(ls_start)
                            )/np.float32((numt*ntod.size))
                    time0 = np.arange(0, numt*ntod.size, dtype=np.float32)
                    time0[:] = time0[:]*step+ls2sol_1year(ls_start)

                newf.log_axis1D(
                    'areo', varnew, dims, longname_txt='solar longitude', units_txt='degree', cart_txt='T')
                newf.log_axis1D('time', time0, dims, longname_txt='sol number',
                                units_txt='days since 0000-00-00 00:00:00', cart_txt='T')  # added AK
            else:
                continue
        elif ndims == 4:
            varnew = npvar
            if avgtime:
                varnew = np.mean(
                    npvar.reshape(-1, Nday, vshape[1], vshape[2], vshape[3]), axis=1)
            if avgtod:
                varnew = varnew.mean(axis=1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3])
            # Rename variable
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt)
            # AK convert surface pressure from mbar to Pa
            if vname2 == 'ps':
                varnew *= 100.
            newf.log_variable(vname2, varnew, newdims,
                              longname_txt2, units_txt2)
        elif ndims == 5:
            varnew = npvar
            if avgtime:
                varnew = np.mean(
                    npvar.reshape(-1, Nday, vshape[1], vshape[2], vshape[3], vshape[4]), axis=1)
            if avgtod:
                varnew = varnew.mean(axis=1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1, vshape[2], vshape[3], vshape[4])
            # Rename variables
            vname2, longname_txt2, units_txt2 = change_vname_longname_unit(
                vname, longname_txt, units_txt)
            newf.log_variable(vname2, varnew, newdims,
                              longname_txt2, units_txt2)
        elif vname == 'tloc':
            if avgtime and not avgtod:
                vname2 = 'time_of_day_16'
                longname_txt2 = 'time_of_day'
                units_txt2 = 'hours since 0000-00-00 00:00:00'
                # Overwrite 'time_of_day' from ('time_of_day_16', 'lon') to 'time_of_day_16'
                newdims = ('time_of_day_16')
                # Every 1.5 hours, centered at half timestep ? AK
                npvar = np.arange(0.75, 24, 1.5)
                newf.log_variable(vname2, npvar, newdims,
                                  longname_txt2, units_txt2)

    return 0


def change_vname_longname_unit(vname, longname_txt, units_txt):
    '''
    Update variable name, longname, and units.
    This was designed specifically for LegacyCGM.nc files.
    '''

    if vname == 'psurf':
        vname = 'ps'
        longname_txt = 'surface pressure'
        units_txt = 'Pa'
    elif vname == 'tsurf':
        vname = 'ts'
        longname_txt = 'surface temperature'
        units_txt = 'K'
    elif vname == 'dst_core_mass':
        vname = 'cor_mass'
        longname_txt = 'dust core mass for the water ice aerosol'
        units_txt = 'kg/kg'

    elif vname == 'h2o_vap_mass':
        vname = 'vap_mass'
        longname_txt = 'water vapor mixing ratio'
        units_txt = 'kg/kg'

    elif vname == 'h2o_ice_mass':
        vname = 'ice_mass'
        longname_txt = 'water ice aerosol mass mixing ratio'
        units_txt = 'kg/kg'

    elif vname == 'dst_mass':
        vname = 'dst_mass'
        longname_txt = 'dust aerosol mass mixing ratio'
        units_txt = 'kg/kg'

    elif vname == 'dst_numb':
        vname = 'dst_num'
        longname_txt = 'dust aerosol number'
        units_txt = 'number/kg'

    elif vname == 'h2o_ice_numb':
        vname = 'ice_num'
        longname_txt = 'water ice aerosol number'
        units_txt = 'number/kg'
    elif vname == 'temp':
        longname_txt = 'temperature'
        units_txt = 'K'
    elif vname == 'ucomp':
        longname_txt = 'zonal wind'
        units_txt = 'm/s'
    elif vname == 'vcomp':
        longname_txt = 'meridional wind'
        units_txt = 'm/s'
    else:
        # Return original values
        pass
    return vname, longname_txt, units_txt


def replace_dims(dims, todflag):
    '''
    Function to replace dimensions with MGCM-like names and remove 'time_of_day'.
    This was designed specifically for LegacyCGM.nc files.
    '''
    newdims = dims
    if 'nlat' in dims:
        newdims = replace_at_index(newdims, newdims.index('nlat'), 'lat')
    if 'nlon' in dims:
        newdims = replace_at_index(newdims, newdims.index('nlon'), 'lon')
    if 'nlay' in dims:
        newdims = replace_at_index(newdims, newdims.index('nlay'), 'pfull')
    if 'ntod' in dims:
        if todflag:
            newdims = replace_at_index(newdims, newdims.index('ntod'), None)
        else:
            newdims = replace_at_index(
                newdims, newdims.index('ntod'), 'time_of_day_16')
    return newdims


def replace_at_index(tuple_dims, idx, new_name):
    '''
    Function to update dimensions.
    Args:
        tup      : the dimensions as tuples e.g. ('pfull', 'nlat', 'nlon')
        idx      : index indicating axis with the dimensions to update (e.g. idx = 1  for 'nlat')
        new_name : new dimension name (e.g. 'latitude')
    '''
    if new_name is None:
        return tuple_dims[:idx]+tuple_dims[idx+1:]
    else:
        return tuple_dims[:idx] + (new_name,) + tuple_dims[idx+1:]


def ls2sol_1year(Ls_deg, offset=True, round10=True):
    '''
    Returns a sol number from the solar longitude.
    Args:
        Ls_deg  : solar longitude in degrees
        offset  : if True, force year to start at Ls 0
        round10 : if True, round to the nearest 10 sols
    Returns:
        Ds: sol number
    ***NOTE***
    For the moment, this is consistent with Ls 0 -> 359.99, but not for monotically increasing Ls.
    '''
    Lsp = 250.99    # Ls at perihelion
    tperi = 485.35  # Time (in sols) at perihelion
    Ns = 668.6      # Number of sols in 1 MY
    e = 0.093379    # From MGCM: modules.f90
    nu = (Ls_deg-Lsp)*np.pi/180
    E = 2*np.arctan(np.tan(nu/2)*np.sqrt((1-e)/(1+e)))
    M = E-e*np.sin(E)
    Ds = M/(2*np.pi)*Ns+tperi
    # Offset correction:
    if offset:
        # Ds is a float
        if len(np.atleast_1d(Ds)) == 1:
            Ds -= Ns
            if Ds < 0:
                Ds += Ns
        # Ds is an array
        else:
            Ds -= Ns
            Ds[Ds < 0] = Ds[Ds < 0]+Ns
    if round:
        Ds = np.round(Ds, -1)  # -1 means round to the nearest 10
    return Ds

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == "__main__":
    main()
