#!/usr/bin/env python3

#Load generic Python Modules 
import argparse #parse arguments
import os       #access operating systems function
import subprocess #run command
import sys       #system command
import time

#TODO delete when done testing
'''
sys.path.append('/Users/akling/amesgcm/amesgcm/')
from Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
from FV3_utils import fms_press_calc,fms_Z_calc,pinterp, find_n
from Ncdf_wrapper import Ncdf
'''
from amesgcm.FV3_utils import fms_press_calc,fms_Z_calc,pinterp,find_n
from amesgcm.Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
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
                 help="""> Usage: MarsInterp ****.atmos.average.nc -t pstd \n"""
                      """Default is 'pstd' \n"""  # TODO, 'zstd', 'zagl'
                      """ \n""")      

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

#Default levels
pstd_in=np.array([1.0e+03, 9.5e+02, 9.0e+02, 8.5e+02, 8.0e+02, 7.5e+02, 7.0e+02,
       6.5e+02, 6.0e+02, 5.5e+02, 5.0e+02, 4.5e+02, 4.0e+02, 3.5e+02,
       3.0e+02, 2.5e+02, 2.0e+02, 1.5e+02, 1.0e+02, 7.0e+01, 5.0e+01,
       3.0e+01, 2.0e+01, 1.0e+01, 7.0e+00, 5.0e+00, 3.0e+00, 2.0e+00,
       1.0e+00, 5.0e-01, 3.0e-01, 2.0e-01, 1.0e-01, 5.0e-02, 3.0e-02,
       1.0e-02])
       
def compute_p_3D(ps,ak,bk,temp):
    dim_out=temp.shape
    p_3D= fms_press_calc(ps,ak,bk,lev_type='full')
    if len(p_3D.shape)==4:p_3D=p_3D.transpose([0,3,1,2])# p_3D [tim,lat,lon,lev] ->[tim, lev, lat, lon]
    if len(p_3D.shape)==3:p_3D=p_3D.transpose([2,0,1]) #p_3D [lat,lon,lev] ->    [lev, lat, lon]
    return p_3D.reshape(dim_out)


def compute_zfull(ps,ak,bk,temp):
    """
    Compute the altitude AGL in m
    """
    dim_out=temp.shape
    zfull=fms_Z_calc(ps,ak,bk,temp.transpose([0,2,3,1]),topo=0.,lev_type='full') #transpose temp like PRESS_h
    if len(zfull.shape)==4:
        zfull=fms_Z_calc(ps,ak,bk,temp.transpose([0,2,3,1]),topo=0.,lev_type='full') #transpose temp like PRESS_h
        zfull=zfull.transpose([0,3,1,2])# p_3D [tim,lat,lon,lev] ->[tim, lev, lat, lon]
    if len(zfull.shape)==3:
        zfull=fms_Z_calc(ps,ak,bk,temp.transpose([1,2,0]),topo=0.,lev_type='full') #transpose temp like PRESS_h
        zfull=zfull.transpose([2,0,1]) #p_3D [lat,lon,lev] ->    [lev, lat, lon]
    return zfull.reshape(dim_out)


def main():
    start_time = time.time()
    debug =parser.parse_args().debug
    #load all the .nc files
    file_list=parser.parse_args().input_file
    interp_type=parser.parse_args().type
    
    
    if interp_type not in ['pstd','zstd ','zagl']: 
        prRed("Interpolation type '%s' is not supported, use  'pstd','zstd' or 'zagl'"%(interp_type))
        exit()
    #TODO delete the following when  zagl and zstd  are supported
    if interp_type not in ['pstd']: 
        prRed("Interpolation type '%s' is not yet supported"%(interp_type))
        exit()
    
    #For all the files
    for ifile in file_list:
        #First check if file is present on the disk (Lou only)
        check_file_tape(ifile)
        newname=filepath+'/'+ifile[:-3]+'_'+interp_type+'.nc'
        
        if interp_type=='pstd':
            longname_txt= 'standard pressure'
            units_txt= 'Pa'

        
        #=================================================================
        #=======================Interpolate action========================
        #=================================================================

        fNcdf=Dataset(ifile, 'r', format='NETCDF4_CLASSIC')
        # Load ak and bk for pressure calculation, do only once.
        if ifile==file_list[0]:
            ak=np.array(fNcdf.variables['pk'])
            bk=np.array(fNcdf.variables['bk'])
            
        ps=np.array(fNcdf.variables['ps'])
        
        if len(ps.shape)==3:
            transposeP=[3,0,1,2] #p_3D [tim,lat,lon,lev] ->[lev,tim,lat, lon]
        elif len(ps.shape)==4:    
            prRed('Interpolation of atmos_diurn not yet supported')
            exit()
        #==This with loop suppresses divided by zero error for the Legacy GCM==
        with np.errstate(divide='ignore', invalid='ignore'):    
            p_3D= fms_press_calc(ps,ak,bk,lev_type='full').transpose(transposeP)
            
        var_list=fNcdf.variables.keys()
        fnew = Ncdf(newname,'Pressure interpolation using MarsInterp.py')
        #===========      Replicate existing DIMENSIONS but pfull  =================
        fnew.copy_all_dims_from_Ncfile(fNcdf,exclude_dim=['pfull'])
        
        fnew.add_dim_with_content('pstd',pstd_in,longname_txt,units_txt) #Add new dimensions
        fnew.copy_Ncaxis_with_content(fNcdf.variables['lon'])
        fnew.copy_Ncaxis_with_content(fNcdf.variables['lat'])
        fnew.copy_Ncaxis_with_content(fNcdf.variables['time'])
        
        #---Time deserve special treatment:
        #fnew.copy_Ncdim_with_content(fNcdf.variables['time'])
        '''
        fnew.add_dimension('time',None)
        t=fNcdf.variables['time'][:]
        fnew.log_axis1D('time',t,('time'),longname_txt='sol number',unit_txt='',cart_txt='')
        '''

        
        #We will re-use the indices for each files, this speeds-up the calculation
        compute_indices=True    
        for ivar in var_list: 
            
            if fNcdf.variables[ivar].dimensions==('time','pfull', 'lat', 'lon'):
                if compute_indices:
                    prCyan("Computing indices ...")
                    index=find_n(p_3D,pstd_in)
                    compute_indices=False
                    
                #print("\r Interpolating: %s ..."%(ivar), end='')
                prCyan("Interpolating: %s ..."%(ivar))
                varIN=fNcdf.variables[ivar][:]
                #==This with loop suppresses divided by zero error for the Legacy GCM==
                with np.errstate(divide='ignore', invalid='ignore'):
                    varOUT=pinterp(varIN.transpose([1,0,2,3]),p_3D ,pstd_in,True,index).transpose([1,0,2,3])
                fnew.log_variable(ivar,varOUT,('time','pstd', 'lat', 'lon'),fNcdf.variables[ivar].long_name,fNcdf.variables[ivar].units)
            else:
                
                if  ivar not in ['time','pfull', 'lat', 'lon','pstd','phalf','pk','bk']:
                    #print("\r Copying over: %s..."%(ivar), end='')
                    prCyan("Copying over: %s..."%(ivar))
                    fnew.copy_Ncvar(fNcdf.variables[ivar]) 
                    
             
            '''
            try:
                
                fileNC=Dataset(ifile, 'a', format='NETCDF4')
                #---temp and ps are always needed---
                dim_out=fileNC.variables['temp'].dimensions #get dimension
                temp=fileNC.variables['temp'][:,:,:,:]
                ps=fileNC.variables['ps'][:,:,:]
                #zfull=compute_zfull(ps,ak,bk,temp)
                OUT=compute_zfull(ps,ak,bk,temp)

                #filter nana
                OUT[np.isnan(OUT)]=fill_value
                #Log the variable
                var_Ncdf = fileNC.createVariable(ivar,'f4',dim_out)
                var_Ncdf.long_name=VAR[ivar][0]
                var_Ncdf.units=    VAR[ivar][1]
                var_Ncdf[:]= OUT
                fileNC.close()

                print('%s: \033[92mDone\033[00m'%(ivar))
            except Exception as exception:
                if debug:raise
                print(exception.__class__.__name__ + ": " + exception.message)
                if exception.message=='NetCDF: String match to name in use':
                    prCyan("""Delete existing var %s with 'MarsVars %s -rm %s'"""%(ivar,ifile,ivar))

            '''
        print('\r ', end='')    
        fNcdf.close()
        fnew.close()
        print("Completed in %.3f sec" % (time.time() - start_time))  
if __name__ == '__main__':
    main()
