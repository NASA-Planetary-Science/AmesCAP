#!/usr/bin/env python

#Load generic Python Modules
import argparse #parse arguments
import os       #access operating systems function
import subprocess #run command
import sys       #system command

sys.path.append(r'/u/mkahre/MCMC/analysis/lib/')
from FV3_utils import fms_press_calc,fms_Z_calc,dvar_dh
from Script_utils import check_file_tape,prYellow,prRed,prCyan,print_fileContent

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
parser = argparse.ArgumentParser(description="""\033[93m MarsVars, variable manager, use to add or remove variables to the diagnostic files\n \033[00m""",
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+', #sys.stdin
                             help='***.nc file or list of ***.nc files ')
parser.add_argument('-add','--add', nargs='+',default=[],
                 help='Add a new variable to file  \n'
                      '> Usage: MarsVars ****.atmos.average.nc -add rho\n'
                      'SUPPORTED:\n\033[96m'
                      'rho   (density)              Req. [ps,temp] \n'
                      'theta (pot. temperature)     Req. [ps,temp] \n'
                      'w     (vertical winds)       Req. [ps,temp,omega] \n'
                      'zfull (altitude AGL)         Req. [ps,temp] \n'
                      'N     (Brunt Vaisala freq)   Req. [ps,temp] \n'
                      'Ri    (Richardson number)    Req. [ps,temp] \n'
                       '\033[00m')

parser.add_argument('-zdiff','--zdiff', nargs='+',default=[],
                 help="""Differentiate a variable 'var' with respect to the z axis\n"""
                      """A new a variable dvar_dz in [Unit/m] will be added o the file\n"""
                      """> Usage: MarsVars ****.atmos.average.nc -zdiff temp\n"""
                      """ \n""")


parser.add_argument('-rm','--remove', nargs='+',default=[],
                 help='Remove any variable from the file  \n'
                      '> Usage: MarsVars ****.atmos.average.nc -rm rho \n')

parser.add_argument('--debug',  action='store_true', help='Debug flag: release the exceptions')

#=====================================================================
# This is the list of supported variable to add: short name, longname, units
#=====================================================================
VAR= {'rho'   :['density (added postprocessing)','kg/m3'],
      'theta' :['potential temperature (added postprocessing)','K'],
      'w'     :['vertical wind (added postprocessing)','m/s'],
      'zfull' :['layer midpoint AGL (added postprocessing)','m'],
      'N'     :['Brunt Vaisala frequency (added postprocessing)','rad/s'],
      'Ri'    :['Richardson number (added postprocessing)','none']    }
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

def compute_N(theta,zfull):
    """
    Compute the Brunt Vaisala freqency in rad/s
    """
    dtheta_dz  = dvar_dh(theta.transpose([1,0,2,3]),zfull.transpose([1,0,2,3])).transpose([1,0,2,3])

    return np.sqrt(g/theta*dtheta_dz)

def main():
    #load all the .nc files
    file_list=parser.parse_args().input_file
    add_list=parser.parse_args().add
    zdiff_list=parser.parse_args().zdiff
    remove_list=parser.parse_args().remove
    debug =parser.parse_args().debug

    #Check if an operation is requested, otherwise print file content.
    if not (add_list or zdiff_list or remove_list):
        print_fileContent(file_list[0])
        prYellow(''' ***Notice***  No operation requested, use '-add var',  '-zdiff var' ''')
        exit() #Exit cleanly

    #For all the files
    for ifile in file_list:
        #First check if file is present on the disk (Lou only)
        check_file_tape(ifile)

        #=================================================================
        #====================Remove action================================
        #=================================================================

        if remove_list:
            cmd_txt='ncks --version'
            try:
                #Test if ncks is available, make recommendation and  exit otherwise--
                subprocess.check_call(cmd_txt,shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
            except subprocess.CalledProcessError:
                prYellow("Enable the -rm option by loading the ncks dependencies:")
                prCyan('    module load comp-intel/2018.0.128')
                prCyan('    module load hdf5/1.8.18_serial')
                prCyan('    module load hdf4/4.2.12')
                prCyan('    module load netcdf/4.4.1.1_serial')
                prCyan('    module load nco')
                exit()


            for ivar in remove_list:
                print('Creating new file %s without %s:'%(ifile,ivar))
                cmd_txt='ncks -C -O -x -v %s %s %s'%(ivar,ifile,ifile)
                try:
                    subprocess.check_call(cmd_txt,shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
                except Exception as exception:
                    print(exception.__class__.__name__ + ": " + exception.message)


        #=================================================================
        #=======================Add action================================
        #=================================================================

        #If the list is not empty, load ak and bk for pressure calculation, those are always needed.
        if add_list:
            name_fixed=ifile[0:5]+'.fixed.nc'
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4')
            variableNames = f_fixed.variables.keys();
            ak=np.array(f_fixed.variables['pk'])
            bk=np.array(f_fixed.variables['bk'])
            f_fixed.close()
        #----
            #----Check if the variable is currently supported---
        for ivar in add_list:
            if ivar not in VAR.keys():
                 prRed("Variable '%s' is not supported"%(ivar))
            else:
                print('Processing: %s...'%(ivar))
                try:
                    fileNC=Dataset(ifile, 'a', format='NETCDF4')
                    #---temp and ps are always needed---
                    dim_out=fileNC.variables['temp'].dimensions #get dimension
                    temp=fileNC.variables['temp'][:,:,:,:]
                    ps=fileNC.variables['ps'][:,:,:]
                    #----

                    if ivar=='rho':OUT=compute_rho(ps,ak,bk,temp)
                    if ivar=='theta':OUT=compute_theta(ps,ak,bk,temp)
                    if ivar=='w':
                        omega=fileNC.variables['omega'][:,:,:,:]
                        rho  = compute_rho(ps,ak,bk,temp)
                        OUT=compute_w(rho,omega)
                    if ivar=='zfull': OUT=compute_zfull(ps,ak,bk,temp)
                    if ivar=='N':
                        theta=compute_theta(ps,ak,bk,temp)
                        zfull=compute_zfull(ps,ak,bk,temp)
                        OUT=compute_N(theta,zfull)
                    if ivar=='Ri':
                        theta=compute_theta(ps,ak,bk,temp)
                        zfull=compute_zfull(ps,ak,bk,temp)
                        N=compute_N(theta,zfull)

                        ucomp=fileNC.variables['ucomp'][:,:,:,:]
                        vcomp=fileNC.variables['vcomp'][:,:,:,:]
                        du_dz=dvar_dh(ucomp.transpose([1,0,2,3]),zfull.transpose([1,0,2,3])).transpose([1,0,2,3])
                        dv_dz=dvar_dh(vcomp.transpose([1,0,2,3]),zfull.transpose([1,0,2,3])).transpose([1,0,2,3])
                        OUT=N**2/(du_dz**2+dv_dz**2)
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

        #=================================================================
        #=============Vertical Differentiation action=====================
        #=================================================================

        #ak and bk are needed to derive the distance between layher pfull
        if zdiff_list:
            name_fixed=ifile[0:5]+'.fixed.nc'
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4')
            variableNames = f_fixed.variables.keys();
            ak=np.array(f_fixed.variables['pk'])
            bk=np.array(f_fixed.variables['bk'])
            f_fixed.close()

        for idiff in zdiff_list:
            fileNC=Dataset(ifile, 'a', format='NETCDF4')

            if idiff not in fileNC.variables.keys():
                prRed("zdiff error: variable '%s' is not present in %s"%(idiff, ifile))
                fileNC.close()
            else:
                print('Differentiating: %s...'%(idiff))

                try:
                    var=fileNC.variables[idiff][:,:,:,:]
                    newUnits=fileNC.variables[idiff].units[:-2]+'/m]' #remove the last ']' to update units, e.g turn '[kg]' to '[kg/m]'
                    newLong_name='vertical gradient of '+fileNC.variables[idiff].long_name

                    #---temp and ps are always needed---
                    dim_out=fileNC.variables['temp'].dimensions #get dimension
                    temp=fileNC.variables['temp'][:,:,:,:]
                    ps=fileNC.variables['ps'][:,:,:]
                    zfull=fms_Z_calc(ps,ak,bk,temp,topo=0.,lev_type='full').transpose([3,0,1,2]) #z is first axis

                    #differentiate the variable zith respect to z:
                    darr_dz=dvar_dh(var.transpose([1,0,2,3]),zfull).transpose([1,0,2,3])

                    #Log the variable
                    var_Ncdf = fileNC.createVariable('d_dz_'+idiff,'f4',dim_out)
                    var_Ncdf.long_name=newLong_name
                    var_Ncdf.units=   newUnits
                    var_Ncdf[:]=  darr_dz
                    fileNC.close()

                    print('%s: \033[92mDone\033[00m'%('d_dz_'+idiff))
                except Exception as exception:
                    if debug:raise
                    print(exception.__class__.__name__ + ": " + exception.message)
                    if exception.message=='NetCDF: String match to name in use':
                        prCyan("""Delete existing var %s with 'MarsVars %s -rm %s'"""%(ivar,ifile,ivar))




if __name__ == '__main__':
    main()
#
