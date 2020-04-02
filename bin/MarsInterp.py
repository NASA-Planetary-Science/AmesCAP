#!/usr/bin/env python

#Load generic Python Modules
import argparse #parse arguments
import os       #access operating systems function
import subprocess #run command
import sys       #system command

#TODO delete when done testing
sys.path.append('/Users/akling/amesgcm/amesgcm/')
from Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent
from FV3_utils import fms_press_calc,fms_Z_calc,Ncdf

#from amesgcm.FV3_utils import fms_press_calc,fms_Z_calc,dvar_dh
#from amesgcm.Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent

#=====Attempt to import specific scientic modules one may not find in the default python on NAS ====
try:
    import matplotlib
    matplotlib.use('Agg') # Force matplotlib to not use any Xwindows backend.
    import numpy as np
    from netCDF4 import Dataset, MFDataset

except ImportError as error_msg:
    prYellow("Error while importing modules")
    prYellow('Your are using python '+str(sys.version_info[0:3])+', recommend version is (2, 7, 12) ' )
    prYellow('Please load recommended version with:')
    prCyan('        module load python/2.7.12\n')
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
parser.add_argument('-t','--type',type=str,default='plevs',
                 help="""> Usage: MarsInterp ****.atmos.average.nc -t plevs \n"""
                      """Available type are  'plevs', 'zlevs', 'zagl'\n"""
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

def compute_p_3D(ps,ak,bk,temp):
    dim_out=temp.shape
    p_3D= fms_press_calc(ps,ak,bk,lev_type='full')
    if len(p_3D.shape)==4:p_3D=p_3D.transpose([0,3,1,2])# p_3D [tim,lat,lon,lev] ->[tim, lev, lat, lon]
    if len(p_3D.shape)==3:p_3D=p_3D.transpose([2,0,1]) #p_3D [lat,lon,lev] ->    [lev, lat, lon]
    return p_3D.reshape(dim_out)

def compute_rho(ps,ak,bk,temp):
    """
    Return the density in kg/m3
    """
    return compute_p_3D(ps,ak,bk,temp)/(rgas*temp)

def compute_theta(ps,ak,bk,temp):
    """
    Return the potential temperature in K
    """

    theta_exp= R/(M_co2*Cp)
    p_3D= compute_p_3D(ps,ak,bk,temp)

    ps_shape=ps.shape #(time,lat,lon) is transformed into (time,1,lat,lon) in the next line with np.reshape
    return temp*(np.reshape(ps,(ps_shape[0],1,ps_shape[1],ps_shape[2]))/p_3D)**(theta_exp) #potential temperature

def compute_w(rho,omega):
    return -omega/(rho*g)

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
    debug =parser.parse_args().debug
    #load all the .nc files
    file_list=parser.parse_args().input_file
    interp_type=parser.parse_args().type
    print(type(interp_type))
    if interp_type not in ['plevs','zlevs','zagl']: 
        prRed("Interpolation type '%s' is not supported, choose 'plevs','zlevs' or 'zagl'"%(interp_type))
        exit()
    
    #For all the files
    for ifile in file_list:
        #First check if file is present on the disk (Lou only)
        check_file_tape(ifile)
        newname=filepath+'/'+ifile[:-3]+'_'+interp_type+'.nc'
        prCyan(newname)
        
        #=================================================================
        #=======================Interpolate action========================
        #=================================================================

        fNcdf=Dataset(ifile, 'r', format='NETCDF4_CLASSIC')
        # Load ak and bk for pressure calculation, do only once.
        if ifile==file_list[0]:
            ak=np.array(fNcdf.variables['pk'])
            bk=np.array(fNcdf.variables['bk'])
            
        var_list=fNcdf.variables.keys()
        fnew = Ncdf(newname,'Pressure interpolation using MarsInterp.py')
        #===========      Replicate existing dimensions  =================
        var=histfile.variables['longitude']
        npvar=var[:]
        a=add_dim(histfile,newfavg,'lon',npvar.size,npvar,getattr(var,'units'))
        var=histfile.variables['latitude']
        npvar=var[:]
        a=add_dim(histfile,newfavg,'lat',npvar.size,npvar,getattr(var,'units'))
        newfavg.add_dimension('time',None)

            #----Check if the variable is currently supported---
        for ivar in var_list: 
            print('Processing: %s...'%(ivar))
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
        fNcdf.close()
        fnew.close()
if __name__ == '__main__':
    main()
