#!/usr/bin/env python3

#Load generic Python Modules
import argparse #parse arguments
import os       #access operating systems function
import subprocess #run command
import sys       #system command
import warnings #Suppress certain errors when dealing with NaN arrays


from amesgcm.FV3_utils import fms_press_calc,fms_Z_calc,dvar_dh,cart_to_azimut_TR,mass_stream,zonal_detrend,spherical_div,spherical_curl,frontogenesis
from amesgcm.Script_utils import check_file_tape,prYellow,prRed,prCyan,prGreen,prPurple, print_fileContent,FV3_file_type,filter_vars
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
    prYellow('Please, source your virtual environment');prCyan('    source amesGCM3/bin/activate \n')
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
parser = argparse.ArgumentParser(description="""\033[93m MarsVars, variable manager,  utility to add or remove variables to the diagnostic files\n Use MarsFiles ****.atmos.average.nc to view file content \033[00m""",
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+', #sys.stdin
                             help='***.nc file or list of ***.nc files ')  
parser.add_argument('-add','--add', nargs='+',default=[],
                 help='Add a new variable to file  \n'
                      '> Usage: MarsVars ****.atmos.average.nc -add rho\n'
                      '\033[96mON NATIVE FILES:\n'
                      'rho        (density)                         Req. [ps,temp] \n'
                      'theta      (pot. temperature)                Req. [ps,temp] \n'
                      'pfull3D    (pressure at layer midpoint)      Req. [ps,temp] \n' 
                      'DP         (layer pressure thickness)        Req. [ps,temp] \n' 
                      'zfull      (altitude AGL)                    Req. [ps,temp] \n'       
                      'DZ         (layer altitude thickness)        Req. [ps,temp] \n'                
                      'w          (vertical winds)                  Req. [ps,temp,omega] \n'
                      'wdir       (wind direction)                  Req. [ucomp,vcomp] \n'
                      'wspeed     (wind magnitude)                  Req. [ucomp,vcomp] \n'
                      'N          (Brunt Vaisala freq)              Req. [ps,temp] \n'
                      'Ri         (Richardson number)               Req. [ps,temp] \n'
                      'Tco2       (CO2 condensation temperature)    Req. [ps,temp] \n'
                      'scorer_wl  (Scorer horizontal wavelength)    Req. [ps,temp,ucomp] \n'
                      'div        (divergence)                      Req. [ucomp,vcomp] \n'
                      'curl       (relative vorticity)              Req. [ucomp,vcomp] \n'
                      'fn         (frontogenesis)                   Req. [ucomp,vcomp,theta] \n'
                      ' \nNOTE:                    \n'
                      '     Some support on interpolated files, in particular if pfull3D \n'
                      '         and zfull are added before interpolation to _pstd, _zagl, _zstd. \n'
                      '\033[00m\n' 
                      '\033[93mON INTERPOLATED FILES :                                     \n'
                      'msf        (mass stream function)              Req. [vcomp] \n'
                      'ep         (wave potential energy)             Req. [temp] \n'
                      'ek         (wave kinetic energy)               Req. [ucomp,vcomp] \n'
                      'mx         (vertical flux of zonal momentum)   Req. [ucomp,w] \n' 
                      'my         (vertical flux of merid. momentum)  Req. [vcomp,w]  \n'
                      'ax         (zonal wave-mean flow forcing)      Req. [ucomp,w,rho]  \n'
                      'ay         (merid. wave-mean flow forcing)     Req. [ucomp,w,rho]  \n'
                      'tp_t       (norm. temperature perturbation)    Req. [temp]  \n'
                      '\033[00m')  


parser.add_argument('-zdiff','--zdiff', nargs='+',default=[],
                 help="""Differentiate a variable 'var' with respect to the z axis\n"""
                      """A new a variable d_dz_var in [Unit/m] will be added o the file\n"""
                      """> Usage: MarsVars ****.atmos.average.nc -zdiff temp\n"""  
                      """ \n""")    
parser.add_argument('-col','--col', nargs='+',default=[],
                 help="""Integrate a mixing ratio of a  variable 'var' through the column\n"""
                      """A new a variable  var_col in [kg/m2] will be added o the file\n"""
                      """> Usage: MarsVars ****.atmos.average.nc -col ice_mass \n"""
                      """ \n""")               
                             
parser.add_argument('-zd','--zonal_detrend', nargs='+',default=[],
                 help="""Detrend a variable by substracting the zonal mean value \n"""
                      """A new a variable  var_p (for prime) will be added to the file \n"""
                      """> Usage: MarsVars ****.atmos.average.nc -zd ucomp \n"""
                      """ \n""")       
                      
parser.add_argument('-rm','--remove', nargs='+',default=[],
                 help='Remove any variable from the file  \n'
                      '> Usage: MarsVars ****.atmos.average.nc -rm rho theta \n')      
                      
parser.add_argument('-extract','--extract', nargs='+',default=[],
                 help='Extract variable(s) to a new  _extract.nc file \n'
                      '> Usage: MarsVars ****.atmos.average.nc -extract ps ts \n')                                          

parser.add_argument('--debug',  action='store_true', help='Debug flag: release the exceptions')

#=====================================================================
# This is the list of supported variable to add: short name, longname, units
#=====================================================================
VAR= {'rho'       :['density (added postprocessing)','kg/m3'],
      'theta'     :['potential temperature (added postprocessing)','K'],
      'w'         :['vertical wind (added postprocessing)','m/s'],    
      'pfull3D'   :['pressure at layer midpoint (added postprocessing)','Pa'], 
      'DP'        :['layer thickness (added postprocessing)','Pa'],   
      'zfull'     :['altitude  AGL at layer midpoint (added postprocessing)','m'],
      'DZ'        :['layer thickness (added postprocessing)','m'],      
      'wdir'      :['wind direction (added postprocessing)','deg'],
      'wspeed'    :['wind speed (added postprocessing)','m/s'],
      'N'         :['Brunt Vaisala frequency (added postprocessing)','rad/s'],   
      'Ri'        :['Richardson number (added postprocessing)','none'], 
      'Tco2'      :['condensation temerature of CO2  (added postprocessing)','K'],
      'div'       :['divergence of the wind field  (added postprocessing)','Hz'],
      'curl'      :['relative vorticity for the wind field  (added postprocessing)','Hz'],
      'scorer_wl' :['Scorer horizontal wavelength L=2.pi/sqrt(l**2)   (added postprocessing)','m'],
      'msf'       :['mass stream function  (added postprocessing)','1.e8 x kg/s'],
      'ep'        :['wave potential energy (added postprocessing)','J/kg'],
      'ek'        :['wave kinetic energy (added postprocessing)','J/kg'] ,
      'mx'        :['vertical flux of zonal momentum (added postprocessing)','J/kg'] ,
      'my'        :['vertical flux of merididional momentum(added postprocessing)','J/kg'] ,
      'ax'        :['zonal wave-mean flow forcing (added postprocessing)','m/s/s'] ,
      'ay'        :['meridional wave-mean flow forcing (added postprocessing)','m/s/s'] ,
      'tp_t'      :['normalized temperature perturbation (added postprocessing)','None'],
      'fn'        :['frontogenesis (added postprocessing)','K m-1 s-1'],
          }                                                                                                       
#=====================================================================
#=====================================================================
#=====================================================================
#TODO : if only one time step, reshape from (lev,lat, lon) to (time.lev,lat, lon)

#Fill values for NaN Do not use np.NaN as this would raise issue when running runpinterp
fill_value=0.
#===Define constants=========
rgas = 189.  # J/(kg-K) => m2/(s2 K)
g    = 3.72  # m/s2
R=    8.314  # J/ mol. K
Cp   = 735.0 #J/K
M_co2 =0.044 # kg/mol
N=0.01       #rad/s  This is used for the Ep calculation
#===========================

def compute_p_3D(ps,ak,bk,shape_out):
    """
    Retunr the 3D pressure field at the layer midpoint. 
    *** NOTE***
    The shape_out argument ensures that, when time=1 (one timestep) results are returned as (1,lev,lat,lon), not (lev,lat,lon)
    """
    p_3D= fms_press_calc(ps,ak,bk,lev_type='full')
    p_3D=p_3D.transpose(lev_T)# p_3D [lev,tim,lat,lon] ->[tim, lev, lat, lon]
    return p_3D.reshape(shape_out)  
    
def compute_rho(p_3D,temp):
    """
    Return the density in [kg/m3]
    """
    return p_3D/(rgas*temp)
    
def compute_theta(p_3D,ps,temp,f_type):
    """
    Return the potential temperature in [K]
    """
    theta_exp= R/(M_co2*Cp)
    #Broadcast dimensions
    ps_shape=ps.shape
    if f_type=='diurn':
        ps_shape=[ps_shape[0],ps_shape[1],1,ps_shape[2],ps_shape[3]]#(time,tod,lat,lon) is transformed into (time,tod,1,lat,lon)
    else:   
        ps_shape=[ps_shape[0],1,ps_shape[1],ps_shape[2]] #(time,lat,lon) is transformed into (time,1,lat,lon)

    return temp*(np.reshape(ps,ps_shape)/p_3D)**(theta_exp) 
    
def compute_w(rho,omega):    
    return -omega/(rho*g)

def compute_zfull(ps,ak,bk,temp):
    """
    Compute the altitude AGL in [m]
    """
    dim_out=temp.shape
    zfull=fms_Z_calc(ps,ak,bk,temp.transpose(lev_T),topo=0.,lev_type='full') # (lev, time, tod, lat,lon)
    zfull=zfull.transpose(lev_T_out)# p_3D [lev,tim,lat,lon] ->[tim, lev, lat, lon] # temp: [tim,tod,lev,lat,lon,lev] ->[lev,time, tod,lat, lon]
    return zfull

def compute_zhalf(ps,ak,bk,temp):
    """
    Compute the altitude AGL in [m]
    """
    dim_out=temp.shape
    zhalf=fms_Z_calc(ps,ak,bk,temp.transpose(lev_T),topo=0.,lev_type='half') # temp: [tim,lev,lat,lon,lev] ->[lev,time, lat, lon]
    zhalf=zhalf.transpose(lev_T_out)# p_3D [lev+1,tim,lat,lon] ->[tim, lev+1, lat, lon]  
    return zhalf

def compute_DZ_full_pstd(pstd,temp,ftype='average'):
    '''
    Return the distance between two layers  mid-point from the standard pressure levels

    Args:
        pstd: 1D  array of standard pressure in [Pa] 
        temp : 3D array of temperature
        ftype: 'daily', 'aveage' or 'diurn'
    Returns:
        DZ_full_pstd: 3D array of distancez  between adjacent layers
        
    *** NOTE***
    In this context p_full = p_std, with the half layers boundaries defined somewhere in between successive layers 
        
    --- Nk --- TOP        ========  p_half
    --- Nk-1 ---                    
                         --------  p_full = p_std   ^
                                                    | DZ_full_pstd
                         ========  p_half           |
    --- 1 ---            --------  p_full = p_std   v
    --- 0 --- SFC        ========  p_half 
                        / / / / 
    '''               
    rgas = 189.  # J/(kg-K) => m2/(s2 K)  
    g    = 3.72  # m/s2   
    if ftype=='diurn':
        axis=2
    else:
        axis=1    
        
    temp=np.swapaxes(temp,0,axis)  

    #Create broadcasting array for pstd
    shape_out=temp.shape
    reshape_shape=[1 for i in range(0,len(shape_out))]     
    reshape_shape[0]=len(pstd)  #e.g [28,1,1,1]
    pstd_b=pstd.reshape(reshape_shape)
    
    DZ_full_pstd=np.zeros_like(temp)
    
    # We Use the average temperature for both layers

    DZ_full_pstd[0:-1,...]=-rgas*0.5*(temp[1:,...]+temp[0:-1,...])/g*np.log(pstd_b[1:,...]/pstd_b[0:-1,...])
    
    # There  is nothing to differentiate the last layer with, so we will copy over the value at N-1
    # Note  that unless you fine-tuned the standard pressure  levels to match the model top, there is typically missing data in the 
    # last few layers so this is not be a big issue.
    DZ_full_pstd[-1,...]=DZ_full_pstd[-2,...]   
    return np.swapaxes(DZ_full_pstd,0,axis)


def compute_N(theta,zfull):
    """
    Compute the Brunt Vaisala freqency in [rad/s]
    """
    dtheta_dz  = dvar_dh(theta.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)        
    return np.sqrt(g/theta*dtheta_dz)


def compute_Tco2(P_3D,temp):
    """
    Compute the frost point of CO2 in [K]
    From [Fannale 1982] Mars: The regolit-atmosphere cap system and climate change. Icarus 
    """
    return np.where(P_3D<518000,-3167.8/(np.log(0.01*P_3D)-23.23),684.2-92.3*np.log(P_3D)+4.32*np.log(P_3D)**2)

def compute_scorer(N,ucomp,zfull):
    """
    Compute the Scorer wavelenght in [m]
    """
    dudz=dvar_dh(ucomp.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)
    dudz2=dvar_dh(dudz.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)
    scorer2= N**2/ucomp**2 -1./ucomp*dudz2    
    return 2*np.pi/np.sqrt(scorer2) 

def compute_DP_3D(ps,ak,bk,shape_out):
    """
    Compute the thickness of a layer in [Pa]
    """
    p_half3D= fms_press_calc(ps,ak,bk,lev_type='half') #[lev,tim,lat,lon]
    DP_3D=p_half3D[1:,...,]- p_half3D[0:-1,...] 
    DP_3D=DP_3D.transpose(lev_T)# p_3D [lev,tim,lat,lon] ->[tim, lev, lat, lon]
    out=DP_3D.reshape(shape_out)
    return out

def compute_DZ_3D(ps,ak,bk,temp,shape_out):
    """
    Compute the thickness of a layer in [Pa]
    """
    z_half3D= fms_Z_calc(ps,ak,bk,temp.transpose(lev_T),topo=0.,lev_type='half')
    DZ_3D=z_half3D[0:-1,...]-z_half3D[1:,...,] #Note the reverse order as Z decreases with increasing levels
    DZ_3D=DZ_3D.transpose(lev_T)# DZ_3D [lev,tim,lat,lon] ->[tim, lev, lat, lon]
    out=DZ_3D.reshape(shape_out)
    return out

def compute_Ep(temp):
    """
    Return the wave potential energy: Ep= 1/2 (g/N)**2 (T'/T)**2 in [J/kg]
    """
    return 0.5*g**2*(zonal_detrend(temp)/(temp*N))**2 
 
def compute_Ek(ucomp,vcomp):
    """
    Return the wave kinetic energy: Ek= 1/2 (u'**2+v'**2) in[J/kg]
    """
    return 0.5*(zonal_detrend(ucomp)**2+zonal_detrend(vcomp)**2 )  

def compute_MF(UVcomp,w):
    """
    Return the zonal or meridional momentum fluxes u'w' or v'w' 
    """
    return zonal_detrend(UVcomp)*zonal_detrend(w)

def compute_WMFF(MF,rho,lev,interp_type):
    """
    Return the zonal or meridional wave-mean flow forcing ax= -1/rho d(rho u'w')/dz in [m/s/s]
    ***NOTE***                                            ay= -1/rho d(rho v'w')/dz in [m/s/s]
    For pstd, we have:
        du/dz= (du/dp).(dp/dz) > du/dz=-rho g (du/dp) with dp/dz = -rho g 
    """
    #Differentiate the variable 
    darr_dz=dvar_dh((rho*MF).transpose(lev_T),lev).transpose(lev_T)
    
    if interp_type=='pstd':
        # We just computed du/dp, we need to multiply by (-rho g) to obtain du/dz
        return g* darr_dz
    else: #With zagl and zstd, levs are already in meters so we already computed du/dz 
        return -1/rho*darr_dz

def main():
    #load all the .nc files
    file_list=parser.parse_args().input_file
    add_list=parser.parse_args().add
    zdiff_list=parser.parse_args().zdiff
    zdetrend_list=parser.parse_args().zonal_detrend
    col_list=parser.parse_args().col
    remove_list=parser.parse_args().remove
    extract_list=parser.parse_args().extract
    debug =parser.parse_args().debug
    
    global lev_T #an array to swap vertical axis first and back: [1,0,2,3] for [time,lev,lat,lon], and [2,1,0,3,4]for [tim, tod,lev, lat, lon]
    global lev_T_out #reshape in zfull, zhalf calculation
    
    #Check if an operation is requested, otherwise print file content.
    if not (add_list or zdiff_list or zdetrend_list or remove_list or col_list or extract_list): 
        print_fileContent(file_list[0])
        prYellow(''' ***Notice***  No operation requested, use '-add var',  '-zdiff var','-zd var', '-col var', '-rm var' ''')
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
                #If ncks is available, use it:--
                subprocess.check_call(cmd_txt,shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w")) 
                for ivar in remove_list:
                    print('Creating new file %s without %s:'%(ifile,ivar))
                    cmd_txt='ncks -C -O -x -v %s %s %s'%(ivar,ifile,ifile)
                    try:
                        subprocess.check_call(cmd_txt,shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w")) 
                    except Exception as exception:
                        print(exception.__class__.__name__ + ": " + exception.message)
            #ncks is not available, we use internal method.            
            except subprocess.CalledProcessError:
                f_IN=Dataset(ifile, 'r', format='NETCDF4_CLASSIC')
                ifile_tmp=ifile[:-3]+'_tmp'+'.nc'
                Log=Ncdf(ifile_tmp,'Edited in postprocessing')
                Log.copy_all_dims_from_Ncfile(f_IN)
                Log.copy_all_vars_from_Ncfile(f_IN,remove_list)
                f_IN.close()
                Log.close()
                cmd_txt='mv '+ifile_tmp+' '+ifile
                p = subprocess.run(cmd_txt, universal_newlines=True, shell=True)
                prCyan(ifile+' was updated')
                
        #=================================================================
        #====================Extract action================================
        #=================================================================
        
        if extract_list: 
            f_IN=Dataset(ifile, 'r', format='NETCDF4_CLASSIC')
            exclude_list = filter_vars(f_IN,parser.parse_args().extract,giveExclude=True) # variable to exclude
            print()
            ifile_tmp=ifile[:-3]+'_extract.nc'
            Log=Ncdf(ifile_tmp,'Edited in postprocessing')
            Log.copy_all_dims_from_Ncfile(f_IN)
            Log.copy_all_vars_from_Ncfile(f_IN,exclude_list)
            f_IN.close()
            Log.close()
            prCyan(ifile+' was created')                
 
        #=================================================================
        #=======================Add action================================
        #=================================================================
        
        #If the list is not empty, load ak and bk for pressure calculation, those are always needed.
        if add_list: 
            name_fixed=ifile[0:5]+'.fixed.nc'
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4_CLASSIC')
            variableNames = f_fixed.variables.keys();
            ak=f_fixed.variables['pk'][:]
            bk=f_fixed.variables['bk'][:]
            f_fixed.close()
        #----    
            #----Check if the variable is currently supported---
        for ivar in add_list:
            if ivar not in VAR.keys():
                 prRed("Variable '%s' is not supported"%(ivar))
            else:
                print('Processing: %s...'%(ivar))                
                try:
                    fileNC=Dataset(ifile, 'a', format='NETCDF4_CLASSIC')
                    f_type,interp_type=FV3_file_type(fileNC)
                    #---temp and ps are always needed---
                    dim_out=fileNC.variables['temp'].dimensions #get dimension
                    temp=fileNC.variables['temp'][:]
                    shape_out=temp.shape
                    if f_type=='diurn':
                        lev_T= [2,1,0,3,4] # [tim, tod,lev, lat, lon] >[lev,tod,time,lat,lon] > [time, tod,lev, lat, lon]
                        #                      0    1   2    3    4      2   1    0   3   4       2    1   0    3    4 
                        lev_T_out=[1,2,0,3,4]
                        lev_axis=2 # in atmos_diurn,the levels is the 3rd axis" (time,tod,lev,lat,lon)
                    else:
                        lev_T= [1,0,2,3] # [tim,lev, lat, lon] > [lev,time,lat,lon] > [tim,lev,lat,lon]
                        #                    0   1   2    3        1   0    2    3      1  0    2     3
                        lev_T_out=lev_T
                        lev_axis=1 # in atmos_average, atmos_daily, the levels is the 2nd axis" (time,lev,lat,lon)
                            
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                    #~~~~~~~~~~~~  Non interpolated files ~~~~~~~~~~~~~~~~~
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                     
                    #These are often needed so we will calculate once here
                    if interp_type=='pfull':
                        ps=fileNC.variables['ps'][:]
                        p_3D=compute_p_3D(ps,ak,bk,shape_out)
                    
                    #If using 'pstd', calculating the 3D pressure field is easy    
                    elif interp_type=='pstd':
                        lev=fileNC.variables['pstd'][:]
                        reshape_shape=[1 for i in range(0,len(shape_out))]                 #(  0   1   2   3 )
                        reshape_shape[lev_axis]=len(lev)  #e.g [1,28,1,1]
                        p_3D=lev.reshape(reshape_shape)
                    #If inter_type is 'zstd', or 'zagl', we need the field pfull3D  pre-computed before interpolation.
                    # We use a 'try' statement as some computations (e.g. wind speed) do not require pfull and will work anyway
                    else:
                        try:
                            p_3D=fileNC.variables['pfull3D'][:]
                        except:
                            pass    
                    
                    if ivar=='pfull3D': OUT=p_3D
                    if ivar=='DP':      OUT=compute_DP_3D(ps,ak,bk,shape_out)
                    if ivar=='rho':
                        OUT=compute_rho(p_3D,temp)    
                    if ivar=='theta':
                        OUT=compute_theta(p_3D,ps,temp,f_type)
                    if ivar=='w':
                        omega=fileNC.variables['omega'][:]     
                        rho=compute_rho(p_3D,temp)         
                        OUT=compute_w(rho,omega)
                         
                    if ivar=='zfull': OUT=compute_zfull(ps,ak,bk,temp) #TODO not with _pstd
                    if ivar=='DZ': OUT=compute_DZ_3D(ps,ak,bk,temp,shape_out)
                    
                    if ivar=='wspeed' or ivar=='wdir': 
                        ucomp=fileNC.variables['ucomp'][:] 
                        vcomp=fileNC.variables['vcomp'][:] 
                        theta,mag=cart_to_azimut_TR(ucomp,vcomp,mode='from')
                        if ivar=='wdir':OUT=theta
                        if ivar=='wspeed':OUT=mag
                            
                    if ivar=='N':
                        theta=compute_theta(p_3D,ps,temp,f_type)    
                        zfull=compute_zfull(ps,ak,bk,temp)  #TODO not with _pstd
                        OUT=compute_N(theta,zfull)
                        
                    if ivar=='Ri':
                        theta=compute_theta(p_3D,ps,temp,f_type)   
                        zfull=compute_zfull(ps,ak,bk,temp) #TODO not with _pstd
                        N=compute_N(theta,zfull)  
                        
                        ucomp=fileNC.variables['ucomp'][:] 
                        vcomp=fileNC.variables['vcomp'][:]    
                        du_dz=dvar_dh(ucomp.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)
                        dv_dz=dvar_dh(vcomp.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)
                        OUT=N**2/(du_dz**2+dv_dz**2)

                    if ivar=='Tco2':OUT=compute_Tco2(p_3D,temp)
                    
                    if ivar=='scorer_wl':
                        ucomp=fileNC.variables['ucomp'][:]
                        theta=compute_theta(p_3D,ps,temp,f_type)   
                        zfull=compute_zfull(ps,ak,bk,temp)
                        N=compute_N(theta,zfull)
                        OUT=compute_scorer(N,ucomp,zfull)
                    
                    if ivar in ['div','curl','fn']:
                        lat=fileNC.variables['lat'][:]
                        lon=fileNC.variables['lon'][:]
                        ucomp=fileNC.variables['ucomp'][:]
                        vcomp=fileNC.variables['vcomp'][:]
                    
                    if ivar =='div': OUT=spherical_div(ucomp,vcomp,lon,lat,R=3400*1000.,spacing='regular')
                    if ivar =='curl':OUT=spherical_curl(ucomp,vcomp,lon,lat,R=3400*1000.,spacing='regular')
                    if ivar =='fn': 
                        theta=fileNC.variables['theta'][:]
                        OUT=frontogenesis(ucomp,vcomp,theta,lon,lat,R=3400*1000.,spacing='regular')  
                                          
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
                    #~~~~~~~~~~~~~~~   Interpolated files ~~~~~~~~~~~~~~~~~~~
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
                    
                    #Common to all interpolated files:
                    if interp_type!='pfull':
                        lev=fileNC.variables[interp_type][:]
                    
                    if ivar=='msf':
                        vcomp=fileNC.variables['vcomp'][:]
                        lat=fileNC.variables['lat'][:]
                        if f_type=='diurn':
                            #[time,tod,lev,lat,lon] > [lev,lat,time,tod,lon]  >  [time,tod,lev,lat,lon]
                            #  0    1   2   3   4       2   3    0   1   4         2    3   0   1   4
                            OUT=mass_stream(vcomp.transpose([2,3,0,1,4]),lat,lev,type=interp_type).transpose([2,3,0,1,4]) 
                        else:    
                            #[time,lev,lat,lon] > [lev,lat,lon,time]  >  [time,lev,lat,lon] 
                            #  0    1   2   3       1   2    3   0         3    0   1   2
                            OUT=mass_stream(vcomp.transpose([1,2,3,0]),lat,lev,type=interp_type).transpose([3,0,1,2]) 
                        
                    if ivar=='ep':
                        #TODO    N=fileNC.variables['N'][:]  >>> Replaced by constant N=0.01 
                        OUT=compute_Ep(temp)       
                    if ivar=='ek':
                        ucomp=fileNC.variables['ucomp'][:]
                        vcomp=fileNC.variables['vcomp'][:]
                        OUT=compute_Ek(ucomp,vcomp)     
                    
                    if ivar=='mx':OUT=compute_MF(fileNC.variables['ucomp'][:],fileNC.variables['w'][:])                      

                    if ivar=='my':OUT=compute_MF(fileNC.variables['vcomp'][:],fileNC.variables['w'][:])                      
                    
                    if ivar=='ax':
                        mx=compute_MF(fileNC.variables['ucomp'][:],fileNC.variables['w'][:])  
                        rho=fileNC.variables['rho'][:]
                        OUT=compute_WMFF(mx,rho,lev,interp_type)
                    
                    if ivar=='ay':
                        my=compute_MF(fileNC.variables['vcomp'][:],fileNC.variables['w'][:])  
                        rho=fileNC.variables['rho'][:]
                        OUT=compute_WMFF(my,rho,lev,interp_type)    
                        
                    if ivar=='tp_t':OUT=zonal_detrend(temp)/temp
                    
                    #filter nan for native files 
                    if interp_type=='pfull':
                        OUT[np.isnan(OUT)]=fill_value
                        
                    #Add nan for interpolated file
                    else :
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore", category=RuntimeWarning)
                            OUT[OUT>1.e30]=np.NaN  
                            OUT[OUT<-1.e30]=np.NaN                          
                    
                    #Log the variable
                    var_Ncdf = fileNC.createVariable(ivar,'f4',dim_out) 
                    var_Ncdf.long_name=VAR[ivar][0]
                    var_Ncdf.units=    VAR[ivar][1]
                    var_Ncdf[:]= OUT
                    fileNC.close()
    
                    print('%s: \033[92mDone\033[00m'%(ivar))
                except Exception as exception:
                    if debug:raise
                    if str(exception)=='NetCDF: String match to name in use':
                        prYellow("""***Error*** Variable already exists""")
                        prYellow("""Delete existing variables %s with 'MarsVars.py %s -rm %s'"""%(ivar,ifile,ivar))
                    
        #=================================================================
        #=============Vertical Differentiation action=====================
        #=================================================================
        
        #ak and bk are needed to derive the distance between layer pfull
        if zdiff_list: 
            name_fixed=ifile[0:5]+'.fixed.nc'
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4_CLASSIC')
            variableNames = f_fixed.variables.keys();
            ak=np.array(f_fixed.variables['pk'])
            bk=np.array(f_fixed.variables['bk'])
            f_fixed.close()
        
        for idiff in zdiff_list:
            fileNC=Dataset(ifile, 'a', format='NETCDF4_CLASSIC')
            f_type,interp_type=FV3_file_type(fileNC)
            if idiff not in fileNC.variables.keys():
                prRed("zdiff error: variable '%s' is not present in %s"%(idiff, ifile))
                fileNC.close()
            else:    
                print('Differentiating: %s...'%(idiff))  
                if f_type=='diurn':
                    lev_T= [2,1,0,3,4] # [tim, tod,lev, lat, lon]
                else: #[time,lat,lon] 
                    lev_T= [1,0,2,3] # [tim,lev, lat, lon]
                try:
                    var=fileNC.variables[idiff][:]
                    newUnits=fileNC.variables[idiff].units[:-2]+'/m]' #remove the last ']' to update units, e.g turn '[kg]' to '[kg/m]'
                    newLong_name='vertical gradient of '+fileNC.variables[idiff].long_name
                    #---temp and ps are always needed---
                    dim_out=fileNC.variables['temp'].dimensions #get dimension
                    if interp_type=='pfull':
                        temp=fileNC.variables['temp'][:]
                        ps=fileNC.variables['ps'][:]
                        zfull=fms_Z_calc(ps,ak,bk,temp.transpose(lev_T),topo=0.,lev_type='full') #z is first axis
                        # atmos_average: zfull= (lev, time, lat, lon)
                        # atmos_diurn :   zfull= (lev, tod, time, lat, lon)
                        #differentiate the variable with respect to z:
                        darr_dz=dvar_dh(var.transpose(lev_T),zfull).transpose(lev_T) 
                        
                    elif interp_type=='pstd':
                        #If pstd, we need the zfull variable
                        if 'zfull' in fileNC.variables.keys():
                            zfull=fileNC.variables['zfull'][:]
                            darr_dz=dvar_dh(var.transpose(lev_T),zfull.transpose(lev_T)).transpose(lev_T)
                        else:
                            lev=fileNC.variables[interp_type][:]
                            temp=fileNC.variables['temp'][:]
                            dzfull_pstd=compute_DZ_full_pstd(lev,temp)
                            darr_dz=dvar_dh(var.transpose(lev_T)).transpose(lev_T)/dzfull_pstd                      
                        
                    elif interp_type in ['zagl','zstd']: #zagl, zstd
                        lev=fileNC.variables[interp_type][:]
                        darr_dz=dvar_dh(var.transpose(lev_T),lev).transpose(lev_T) 

                        
                    #Log the variable
                    var_Ncdf = fileNC.createVariable('d_dz_'+idiff,'f4',dim_out) 
                    var_Ncdf.long_name=newLong_name
                    var_Ncdf.units=   newUnits
                    var_Ncdf[:]=  darr_dz
                    fileNC.close()
                    
                    print('%s: \033[92mDone\033[00m'%('d_dz_'+idiff))
                except Exception as exception:
                    if debug:raise
                    if str(exception)=='NetCDF: String match to name in use':
                        prYellow("""***Error*** Variable already exists""")
                        prYellow("""Delete existing variable %s with 'MarsVars %s -rm %s'"""%('d_dz_'+idiff,ifile,'d_dz_'+idiff))
                         

        #=================================================================
        #================ Zonal detrending action=========================
        #=================================================================
        
        
        for izdetrend in zdetrend_list:
            fileNC=Dataset(ifile, 'a', format='NETCDF4_CLASSIC')
            f_type,interp_type=FV3_file_type(fileNC)
            if izdetrend not in fileNC.variables.keys():
                prRed("zdiff error: variable '%s' is not present in %s"%(izdetrend, ifile))
                fileNC.close()
            else:    
                print('Detrending: %s...'%(izdetrend))   
                 
                try:
                    var=fileNC.variables[izdetrend][:]
                    newUnits=fileNC.variables[izdetrend].units 
                    newLong_name='zonal perturbation of '+fileNC.variables[izdetrend].long_name
                    dim_out=fileNC.variables[izdetrend].dimensions #get dimension
                            
                    #Log the variable
                    var_Ncdf = fileNC.createVariable(izdetrend+'_p','f4',dim_out) 
                    var_Ncdf.long_name=newLong_name
                    var_Ncdf.units=   newUnits
                    var_Ncdf[:]=  zonal_detrend(var)
                    fileNC.close()
                    
                    print('%s: \033[92mDone\033[00m'%(izdetrend+'_p'))
                except Exception as exception:
                    if debug:raise
                    if str(exception)=='NetCDF: String match to name in use':
                        prYellow("""***Error*** Variable already exists""")
                        prYellow("""Delete existing variable %s with 'MarsVars %s -rm %s'"""%('d_dz_'+idiff,ifile,'d_dz_'+idiff))
                        

        #=================================================================
        #=============  Column  integration   ============================
        #=================================================================
        '''
                          z_top
                          ⌠
        We have col=      ⌡ var rho dz  with dp/dz=-rho g => rho dz = -dp/g
                          0
                      
                      ___ p_sfc
             >  col = \       
                      /__ (var dp/g)  
                        p_top
        ''' 
        #ak and bk are needed to derive the distance between layer pfull
        if col_list:
            name_fixed=ifile[0:5]+'.fixed.nc'
            f_fixed=Dataset(name_fixed, 'r', format='NETCDF4_CLASSIC')
            variableNames = f_fixed.variables.keys();
            ak=np.array(f_fixed.variables['pk'])
            bk=np.array(f_fixed.variables['bk'])
            f_fixed.close()

        for icol in col_list:
            fileNC=Dataset(ifile, 'a') #, format='NETCDF4_CLASSIC
            f_type,interp_type=FV3_file_type(fileNC)
            if icol not in fileNC.variables.keys():
                prRed("column integration error: variable '%s' is not present in %s"%(icol, ifile))
                fileNC.close()
            else:
                print('Performing colum integration: %s...'%(icol))

                try:
                    var=fileNC.variables[icol][:]
                    #prRed(fileNC.variables[icol].units+'|')
                    newUnits=fileNC.variables[icol].units[:-3]+'/m2' # turn 'kg/kg'> to 'kg/m2'
                    newLong_name='column integration of '+fileNC.variables[icol].long_name

                    #---temp and ps are always needed---
                    dim_in=fileNC.variables['temp'].dimensions #get dimension
                    shape_in=fileNC.variables['temp'].shape
                    #TODO edged cases where time =1
                    if f_type=='diurn':
                        #[time,tod,lat,lon] 
                        lev_T= [2,1,0,3,4] # [tim, tod,lev, lat, lon]
                        dim_out=tuple([dim_in[0],dim_in[1],dim_in[3],dim_in[4]])
                        lev_axis=2 # in atmos_diurn,the levels is the 3rd axis" (time,tod,lev,lat,lon)
                    else: #[time,lat,lon] 
                        lev_T= [1,0,2,3] # [tim,lev, lat, lon]
                        dim_out=tuple([dim_in[0],dim_in[2],dim_in[3]])
                        lev_axis=1

                    ps=fileNC.variables['ps'][:] 
                    DP=compute_DP_3D(ps,ak,bk,shape_in)
                    out=np.sum(var*DP/g,axis=lev_axis) 
                    
                    #Log the variable
                    var_Ncdf = fileNC.createVariable(icol+'_col','f4',dim_out)
                    var_Ncdf.long_name=newLong_name
                    var_Ncdf.units=   newUnits
                    var_Ncdf[:]=  out
                    
                    fileNC.close()

                    print('%s: \033[92mDone\033[00m'%(icol+'_col'))
                except Exception as exception:
                    if debug:raise
                    if str(exception)=='NetCDF: String match to name in use':
                        prYellow("""***Error*** Variable already exists""")
                        prYellow("""Delete existing variable %s with 'MarsVars %s -rm %s'"""%(icol+'_col',ifile,icol+'_col'))            

             
if __name__ == '__main__':
    main()        
#
