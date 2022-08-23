#!/usr/bin/env python3

#Load generic Python Modules
import argparse   # parse arguments
import os         # access operating systems function
import subprocess # run command
import sys        # system command
import time       # monitor interpolation time
import re         # string matching module to handle time_of_day_XX

from amesgcm.FV3_utils import fms_press_calc,fms_Z_calc,vinterp,find_n,polar2XYZ,interp_KDTree,axis_interp
from amesgcm.Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
from amesgcm.Script_utils import section_content_amesgcm_profile,find_tod_in_diurn,filter_vars,find_fixedfile
from amesgcm.Ncdf_wrapper import Ncdf

#=====Attempt to import specific scientic modules one may not find in the default python on NAS ====
try:
    import matplotlib
    matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
    import numpy as np
    from netCDF4 import Dataset, MFDataset

except ImportError as error_msg:
    prYellow("Error while importing modules")
    prYellow('Your are using python '+str(sys.version_info[0:3]))
    prYellow('Please, source your virtual environment');prCyan('    source envPython3.7/bin/activate.csh \n')
    print("Error was: "+ error_msg.message)
    exit()
except Exception as exception:
    # Output unexpected Exceptions.
    print(exception, False)
    print(exception.__class__.__name__ + ": " + exception.message)
    exit()


#======================================================
#                  ARGUMENTS PARSER
#======================================================
parser = argparse.ArgumentParser(description="""\033[93m MarsInterp, pressure interpolation on fixed layers\n \033[00m""",
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+', #sys.stdin
                             help='***.nc file or list of ***.nc files ')
parser.add_argument('-t','--type',type=str,default='pstd',
                 help=""">  --type may be 'pstd', 'zstd' or 'zagl' [DEFAULT is pstd, 36 levels] \n"""
                      """>  Usage: MarsInterp.py ****.atmos.average.nc \n"""
                      """          MarsInterp.py ****.atmos.average.nc -t zstd \n""")

parser.add_argument('-l','--level',type=str,default=None,
                 help=""">  Layers ID as defined in your personal ~/.amesgcm_profile hidden file \n"""
                      """"(For 1st time set-up, copy \033[96mcp ~/amesGCM3/mars_templates/amesgcm_profile ~/.amesgcm_profile\033[00m )   \n"""
                      """>  Usage: MarsInterp.py ****.atmos.average.nc -t pstd -l p44 \n"""
                      """          MarsInterp.py ****.atmos.average.nc -t zstd -l phalf_mb  \n""")

parser.add_argument('-include','--include',nargs='+',
                     help="""Only include listed variables. Dimensions and 1D variables are always included \n"""
                         """> Usage: MarsInterp.py *.atmos_daily.nc --include ps ts temp     \n"""
                         """\033[00m""")

parser.add_argument('-e','--ext',type=str,default=None,
                 help="""> Append an extension _ext.nc to the output file instead of replacing any existing file \n"""
                      """>  Usage: MarsInterp.py ****.atmos.average.nc -ext B \n"""
                      """   This will produce   ****.atmos.average_pstd_B.nc files     \n""")
parser.add_argument('--debug',  action='store_true', help='Debug flag: release the exceptions')


#=====================================================================
#=====================================================================
#=====================================================================
#TODO : if only one time step, reshape from (lev,lat, lon) to (time.lev,lat, lon)

#Fill values for NaN Do not use np.NaN as this would raise issue when running runpinterp
fill_value=0.
#===Define constants=========
rgas = 189.            # J/(kg-K) => m2/(s2 K)
g    = 3.72            # m/s2
R=8.314 #J/ mol. K
Cp   = 735.0 #J/K
M_co2 =0.044 # kg/mol
#===========================
filepath=os.getcwd()



def main():
    start_time = time.time()
    debug =parser.parse_args().debug
    #load all the .nc files
    file_list=parser.parse_args().input_file
    interp_type=parser.parse_args().type   #e.g.  'pstd'
    custom_level=parser.parse_args().level #e.g.  'p44'

    #The fixed file is needed if pk, bk are not available in the requested file, or
    # to load the topography is zstd output is requested
    name_fixed=find_fixedfile(filepath,file_list[0])

    # PRELIMINARY DEFINITIONS
    #===========================pstd============================================
    if interp_type=='pstd':
        longname_txt= 'standard pressure'
        units_txt= 'Pa'
        need_to_reverse=False
        interp_technic='log'
        if custom_level:
            content_txt=section_content_amesgcm_profile('Pressure definitions for pstd')
            #print(content_txt)
            exec(content_txt) #load all variables in that section
            lev_in=eval('np.array('+custom_level+')') #copy requested variable
        else:
            #Default levels, this is size 36
            lev_in=np.array([1.0e+03, 9.5e+02, 9.0e+02, 8.5e+02, 8.0e+02, 7.5e+02, 7.0e+02,
                6.5e+02, 6.0e+02, 5.5e+02, 5.0e+02, 4.5e+02, 4.0e+02, 3.5e+02,
                3.0e+02, 2.5e+02, 2.0e+02, 1.5e+02, 1.0e+02, 7.0e+01, 5.0e+01,
                3.0e+01, 2.0e+01, 1.0e+01, 7.0e+00, 5.0e+00, 3.0e+00, 2.0e+00,
                1.0e+00, 5.0e-01, 3.0e-01, 2.0e-01, 1.0e-01, 5.0e-02, 3.0e-02,
                1.0e-02])
    #===========================zstd============================================
    elif interp_type=='zstd':
        longname_txt= 'standard altitude'
        units_txt= 'm'
        need_to_reverse=True
        interp_technic='lin'
        if custom_level:
            content_txt=section_content_amesgcm_profile('Altitude definitions for zstd')
            exec(content_txt) #load all variables in that section
            lev_in=eval('np.array('+custom_level+')') #copy requested variable
        else:
            #Default levels, this is size 45
            lev_in=np.array([-7000,-6000,-5000,-4500,-4000,-3500,-3000,-2500,-2000,-1500,-1000,
                    -500,0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,
                    6000,7000,8000,9000,10000,12000,14000,16000,18000,
                    20000,25000,30000,35000,40000,45000,50000,55000,
                    60000,70000,80000,90000,100000])
        try:
            f_fixed=Dataset(name_fixed,'r')
            zsurf=f_fixed.variables['zsurf'][:]
            f_fixed.close()
        except FileNotFoundError:
            prRed('***Error*** Topography is needed for zstd interpolation, however')
            prRed('file %s not found'%(name_fixed))
            exit()
    #===========================zagl============================================
    elif interp_type=='zagl':
        longname_txt= 'altitude above ground level'
        units_txt= 'm'
        need_to_reverse=True
        interp_technic='lin'
        if custom_level:
            content_txt=section_content_amesgcm_profile('Altitude definitions for zagl')
            #print(content_txt)
            exec(content_txt) #load all variables in that section
            lev_in=eval('np.array('+custom_level+')') #copy requested variable
        else:
            #Default levels, this is size 45
            lev_in=np.array([0,500,1000,1500,2000,2500,3000,3500,4000,4500,5000,
                    6000,7000,8000,9000,10000,12000,14000,16000,18000,
                    20000,25000,30000,35000,40000,45000,50000,55000,
                    60000,70000,80000,90000,100000,110000])
    else:
        prRed("Interpolation type '%s' is not supported, use  'pstd','zstd' or 'zagl'"%(interp_type))
        exit()

    #For all the files
    for ifile in file_list:
        #First check if file is present on the disk (Lou only)
        check_file_tape(ifile)

        #Append extension, in any
        if parser.parse_args().ext:
            newname=filepath+'/'+ifile[:-3]+'_'+interp_type+'_'+parser.parse_args().ext+'.nc'
        else:
            newname=filepath+'/'+ifile[:-3]+'_'+interp_type+'.nc'

        #=================================================================
        #=======================Interpolate action========================
        #=================================================================

        fNcdf=Dataset(ifile, 'r', format='NETCDF4_CLASSIC')
        # Load pk and bk and ps for 3D pressure field calculation.
        # We will read the pk and bk for each file in case the vertical resolution is changed.

        try:
            #First try to read pk and bk in the file
            pk=np.array(fNcdf.variables['pk'])
            bk=np.array(fNcdf.variables['bk'])
        except:
            #If pk and bk are not available in the file, try the matching XXXXX.fixed.nc
            name_fixed=find_fixedfile(filepath,ifile)
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4_CLASSIC')
            pk=np.array(f_fixed.variables['pk'])
            bk=np.array(f_fixed.variables['bk'])
            f_fixed.close()

        ps=np.array(fNcdf.variables['ps'])

        if len(ps.shape)==3:
            do_diurn=False
            tod_name='not_used'
            permut=[1,0,2,3] # Put vertical axis first for 4D variable, e.g (time,lev,lat,lon) >>> (lev,time,lat,lon)
                             #                              ( 0    1   2   3 ) >>> ( 1   0    2   3 )
        elif len(ps.shape)==4:
            do_diurn=True
            #find time of day variable name
            tod_name=find_tod_in_diurn(fNcdf)
            permut=[2,1,0,3,4]  #Same for diun files, e.g (time,time_of_day_XX,lev,lat,lon) >>> (lev,time_of_day_XX,time,lat,lon)
                                #                         (  0        1         2   3   4)  >>> ( 2       1          0    3   4 )
        #== Compute levels in the file, these are permutted arrays

        # Suppress divided by zero error ==
        with np.errstate(divide='ignore', invalid='ignore'):
            if interp_type=='pstd':
                L_3D_P= fms_press_calc(ps,pk,bk,lev_type='full') #permuted by default, e.g lev is first

            elif interp_type=='zagl':
                temp=fNcdf.variables['temp'][:]
                L_3D_P= fms_Z_calc(ps,pk,bk,temp.transpose(permut),topo=0.,lev_type='full')

            elif interp_type=='zstd':
                temp=fNcdf.variables['temp'][:]
                #Expend the zsurf array to the time dimension
                zflat=np.repeat(zsurf[np.newaxis,:],ps.shape[0],axis=0)
                if do_diurn:
                    zflat=np.repeat(zflat[:,np.newaxis,:,:],ps.shape[1],axis=1)

                L_3D_P= fms_Z_calc(ps,pk,bk,temp.transpose(permut),topo=zflat,lev_type='full')


        fnew = Ncdf(newname,'Pressure interpolation using MarsInterp.py')
        #===========      Replicate existing DIMENSIONS but pfull  =================
        #get all variables in file:
        ###var_list=fNcdf.variables.keys()
        var_list = filter_vars(fNcdf,parser.parse_args().include) # get the variables

        fnew.copy_all_dims_from_Ncfile(fNcdf,exclude_dim=['pfull'])
        fnew.add_dim_with_content(interp_type,lev_in,longname_txt,units_txt) #Add new vertical dimension


        if 'tile' in ifile:
            fnew.copy_Ncaxis_with_content(fNcdf.variables['grid_xt'])
            fnew.copy_Ncaxis_with_content(fNcdf.variables['grid_yt'])
        else:
            fnew.copy_Ncaxis_with_content(fNcdf.variables['lon'])
            fnew.copy_Ncaxis_with_content(fNcdf.variables['lat'])

        fnew.copy_Ncaxis_with_content(fNcdf.variables['time'])

        if do_diurn:fnew.copy_Ncaxis_with_content(fNcdf.variables[tod_name])

        #We will re-use the indices for each files, this speeds-up the calculation
        compute_indices=True
        for ivar in var_list:
            if (fNcdf.variables[ivar].dimensions==('time','pfull', 'lat', 'lon') or
             fNcdf.variables[ivar].dimensions==('time',tod_name,'pfull', 'lat', 'lon') or
             fNcdf.variables[ivar].dimensions==('time','pfull', 'grid_yt', 'grid_xt')):
                if compute_indices:
                    prCyan("Computing indices ...")
                    index=find_n(L_3D_P,lev_in,reverse_input=need_to_reverse)
                    compute_indices=False

                prCyan("Interpolating: %s ..."%(ivar))
                varIN=fNcdf.variables[ivar][:]
                #==This with loop suppresses divided by zero errors==
                with np.errstate(divide='ignore', invalid='ignore'):
                    varOUT=vinterp(varIN.transpose(permut),L_3D_P,
                                   lev_in,type_int=interp_technic,reverse_input=need_to_reverse,
                                   masktop=True,index=index).transpose(permut)
                if not do_diurn:
                    if 'tile' in ifile:
                        fnew.log_variable(ivar,varOUT,('time',interp_type, 'grid_yt', 'grid_xt'),
                                        fNcdf.variables[ivar].long_name,fNcdf.variables[ivar].units)
                    else:
                        fnew.log_variable(ivar,varOUT,('time',interp_type, 'lat', 'lon'),
                                        fNcdf.variables[ivar].long_name,fNcdf.variables[ivar].units)
                else:
                    if 'tile' in ifile:
                        fnew.log_variable(ivar,varOUT,('time',tod_name,interp_type, 'grid_yt', 'grid_xt'),
                                      fNcdf.variables[ivar].long_name,fNcdf.variables[ivar].units)
                    else:
                        fnew.log_variable(ivar,varOUT,('time',tod_name,interp_type, 'lat', 'lon'),
                                      fNcdf.variables[ivar].long_name,fNcdf.variables[ivar].units)
            else:

                if  ivar not in ['time','pfull', 'lat', 'lon','phalf','pk','bk','pstd','zstd','zagl',tod_name,'grid_xt','grid_yt']:
                    #print("\r Copying over: %s..."%(ivar), end='')
                    prCyan("Copying over: %s..."%(ivar))
                    fnew.copy_Ncvar(fNcdf.variables[ivar])


        print('\r ', end='')
        fNcdf.close()
        fnew.close()
        print("Completed in %.3f sec" % (time.time() - start_time))
if __name__ == '__main__':
    main()
