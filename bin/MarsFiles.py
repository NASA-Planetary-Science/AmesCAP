#!/usr/bin/env python3

#
#Assumes the companion scripts are in the same directory
#Load generic Python Modules
import argparse #parse arguments
import sys
import getopt
import os
import re
import glob
import shutil
import subprocess
import numpy as np
from netCDF4 import Dataset


#TODO remove this block to use package instead
#==============
#sys.path.append('/Users/akling/amesgcm/amesgcm/')
#from FV3_utils import Ncdf
#from Script_utils import prYellow,prCyan,prRed
#===========
from amesgcm.FV3_utils import Ncdf
from amesgcm.Script_utils import prYellow,prCyan,prRed
#---


#======================================================
#                  ARGUMENTS PARSER
#======================================================
parser = argparse.ArgumentParser(description="""\033[93m MarsFiles files manager. Used to convert Legacy GCM to FV3 format \n \033[00m""",
                                formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+',
                                help='***.nc file or list of ***.nc files ')

parser.add_argument('-fv3','--fv3', action='store_true',
                    help="""Produce FV3's diurn, average and daily files from Legacy outputs  \n"""
                        """> Usage: MarsFiles LegacyGCM*.nc -fv3 \n"""
                        """\033[00m""")

parser.add_argument('-c','--combine', action='store_true',
                    help="""Combine a sequence of similar files as one file\n"""
                        """> Usage: MarsVars *.atmos_average.nc --combine \n"""
                        """ \n""")

parser.add_argument('--debug',  action='store_true', help='Debug flag: release the exceptions')


#======================================================
#======================================================
#======================================================
#TODO needed?
#sys.path.append(os.getcwd())

def main():

        #exit()
    file_list=parser.parse_args().input_file
    #=======Convert to FV3================
    if parser.parse_args().fv3:

    #argument definitions:
        cwd=os.getcwd()
        path2data=os.getcwd()
        do_multi= False
        do_1year= False #Used with LegacyGCM_1year.nc'
    
        #Get files to process
        histlist=[]
        for filei in file_list:
            histlist.append(path2data+'/'+filei)
        fnum = len(histlist)
        if fnum >= 0:do_multi = True #TODO why not 1?
    
        try:
            hist1year= path2data+'/LegacyGCM_1year.nc'
            file1year= Dataset(hist1year,'r',format='NETCDF4_CLASSIC')
            do_1year = True
        except:
            hist1year= None
            do_1year = False
    
        lsmin = None
        lsmax = None
    
        if do_multi:
            for f in histlist:
                histname = os.path.basename(f)
                ls_l = histname[-12:-9]
                ls_r = histname[-6:-3]
                if lsmin is None:
                    lsmin=ls_l
                else:
                    lsmin = str(min(int(lsmin),int(ls_l))).zfill(3)
                if lsmax is None:
                    lsmax=ls_r
                else:
                    lsmax = str(max(int(lsmax),int(ls_r))).zfill(3)
                a=make_FV3_files(f,True,cwd)
                
    
        #TODO is using sys.prefix reliable outside a virtual environment ?


    elif parser.parse_args().combine:
        #Test if ncks is available , make recommendation and  exit otherwise--
        try:
            subprocess.check_call('ncks --version',shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
        except subprocess.CalledProcessError:
            prYellow("***Error*** ncks dependency is  required to concatenate files")
            exit()
        #now cat together the files
        newfavg="Ls"+lsmin+"_Ls"+lsmax+".atmos_average.nc"
        newfdai="Ls"+lsmin+"_Ls"+lsmax+".atmos_daily.nc"
        newfdiu="Ls"+lsmin+"_Ls"+lsmax+".atmos_diurn.nc"
        tempdir=os.path.join(cwd,'temp')
        os.makedirs(tempdir, exist_ok=True)
        os.chdir(tempdir)

        catavg = "ncrcat ../*.atmos_average.nc "+"00000.atmos_average.nc"
        catdai = "ncrcat ../*.atmos_daily.nc "+"00000.atmos_daily.nc"
        catdiu = "ncrcat ../*.atmos_diurn.nc "+"00000.atmos_diurn.nc"
        p = subprocess.Popen(catavg, universal_newlines=True, shell=True)
        p.wait()
        p = subprocess.Popen(catdai, universal_newlines=True, shell=True)
        p.wait()
        p = subprocess.Popen(catdiu, universal_newlines=True, shell=True)
        p.wait()
        os.chdir(cwd)
        p = subprocess.run('rm -f Ls*.nc', universal_newlines=True, shell=True)
        p = subprocess.run('mv temp/*.nc .', universal_newlines=True, shell=True)
        p = subprocess.run('rm -rf temp/', universal_newlines=True, shell=True)
        #TODO
        #MarsVars.beta 00000.atmos_average.nc -colint vap_mass_micro ice_mass_micro dst_mass_micro
        #runpinterp -l 24
        # TO DO MarsVars.beta 00000.atmos_average_plevs.nc -add mass_stream
        #MarsPlot --do legacy
        if do_1year:
            a=make_FV3_files(hist1year,cwd)
            
    else:
        prRed("""Error: no action requested: use 'MarsFiles *nc --fv3' or 'MarsFiles *nc --fv3 --combine '""")    
        
        


def make_FV3_files(fpath,renameFV3=True,cwd=None):
    '''
    Make FV3-type atmos_average,atmos_daily,atmos_diurn
    Args:
        fpath     : full path to Legacy .nc files
        renameFV3 : rename files from Legacy_Lsxxx_Lsyyy.nc to XXXXX.atmos_average.nc folllowing FV3's convention
        cwd       : output path
    Returns:
        atmos_average,atmos_daily,atmos_diurn
    '''
    
    histname = os.path.basename(fpath)
    if cwd is None:
        histdir = os.path.dirname(fpath)
    else:
        histdir = cwd
    result   = re.search('LegacyGCM_(.*).nc',histname)
    fdate    = result.group(1)
    histfile = Dataset(fpath,'r',format='NETCDF4_CLASSIC')
    histvars = histfile.variables.keys()
    histdims = histfile.dimensions.keys()
    
    #Convert the first Ls in file to a sol number
    if renameFV3:fdate= '%05i'%(ls2sol_1year(histfile.variables['ls'][0]))
    
    newfname_avg = fdate+'.atmos_average.nc'  #5 sol averages over tod and time
    newfname_daily = fdate+'.atmos_daily.nc'  #daily snapshot output...this is exactly the current output?
    newfname_diurn = fdate+'.atmos_diurn.nc'  #5 sol averages over time only
    newfpath_avg = os.path.join(histdir,newfname_avg)
    newfpath_daily = os.path.join(histdir,newfname_daily)
    newfpath_diurn = os.path.join(histdir,newfname_diurn)
    newfavg = Ncdf(newfpath_avg)
    newfdaily = Ncdf(newfpath_daily)
    newfdiurn = Ncdf(newfpath_diurn)


    for dname in histdims:
        if dname == 'nlon':
            var=histfile.variables['longitude']
            npvar=var[:]
            newfavg.add_dim_with_content('lon',npvar,'longitudes',getattr(var,'units'))
            newfdaily.add_dim_with_content('lon',npvar,'longitudes',getattr(var,'units'))
            newfdiurn.add_dim_with_content('lon',npvar,'longitudes',getattr(var,'units'))
        if dname == 'nlat':
            var=histfile.variables['latitude']
            npvar=var[:]
            newfavg.add_dim_with_content('lat',npvar,'latitudes',getattr(var,'units'))
            newfdaily.add_dim_with_content('lat',npvar,'latitudes',getattr(var,'units'))
            newfdiurn.add_dim_with_content('lat',npvar,'latitudes',getattr(var,'units'))
        if dname == 'time':
            newfavg.add_dimension('time',None)
            newfdaily.add_dimension('time',None)
            newfdiurn.add_dimension('time',None)
        if dname == 'ntod':
            dim=histfile.dimensions[dname]
            newfdiurn.add_dimension('time_of_day_16',dim.size)
        if dname == 'nlay':
            nlay=histfile.dimensions[dname]
            num =nlay.size
            nump=num+1
            pref=7.01 #TODO  * 100  change this to Pa
            pk=np.zeros(nump)
            bk=np.zeros(nump)
            pfull=np.zeros(num)
            phalf=np.zeros(nump)
            dsgm=histfile.variables['dsgm']
            sgm =histfile.variables['sgm']
            pk[0]=.08
            for z in range(num):
                bk[z+1] = sgm[2*z+2]
            phalf[:]=pk[:]/100.+pref*bk[:] #TODO remove /100 to output in to Pa
            pfull[:] = (phalf[1:]-phalf[:num])/(np.log(phalf[1:])-np.log(phalf[:num]))
            #

            newfavg.add_dim_with_content('pfull',pfull,'ref full pressure level','mb')
            newfavg.add_dim_with_content('phalf',phalf,'ref half pressure level','mb')
            newfavg.log_axis1D('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pa',cart_txt='')
            newfavg.log_axis1D('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='Pa',cart_txt='')
            #
            newfdaily.add_dim_with_content('pfull',pfull,'ref full pressure level','mb')
            newfdaily.add_dim_with_content('phalf',phalf,'ref half pressure level','mb')
            newfdaily.log_axis1D('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pascal',cart_txt='')
            newfdaily.log_axis1D('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='',cart_txt='')
            #
            newfdiurn.add_dim_with_content('pfull',pfull,'ref full pressure level','mb')
            newfdiurn.add_dim_with_content('phalf',phalf,'ref half pressure level','mb')
            newfdiurn.log_axis1D('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pascal',cart_txt='')
            newfdiurn.log_axis1D('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='',cart_txt='')

#perform various averages
    a=do_avg_vars(histfile,newfavg,True,True)
    a=do_avg_vars(histfile,newfdaily,False,False)
    a=do_avg_vars(histfile,newfdiurn,True,False)

    #Copy Legacy.fixed to current directory 
    cmd_txt='cp '+sys.prefix+'/mars_data/Legacy.fixed.nc '+fdate+'.fixed.nc'
    p = subprocess.run(cmd_txt, universal_newlines=True, shell=True)
 
    newfavg.close()
    newfdaily.close()
    newfdiurn.close()
    print(cwd+'/'+fdate+'.fixed.nc was copied locally')


    return 0


#Function to perform time averages over all fields
def do_avg_vars(histfile,newf,avgtime,avgtod):
    histvars = histfile.variables.keys()
    for vname in histvars:
        var     = histfile.variables[vname]
        npvar = var[:]
        dims  = var.dimensions
        ndims = npvar.ndim
        vshape= npvar.shape
        ntod  = histfile.dimensions['ntod']
        if avgtod:
            newdims  = replace_dims(dims,True)
        elif avgtime:
            newdims  = replace_dims(dims,False)
        else:
            newdims  = replace_dims(dims,True)


        if 'time' in dims:
            tind = dims.index('time')
            tind_new= newdims.index('time')
            numt = histfile.dimensions['time'].size
        #TODO fix time !!
        #now do various time averaging and write to files
        if ndims == 1:
            if vname == 'ls':
                #first check if spans new year
                if not np.all(npvar[1:] >= npvar[:-1]):
                    year = 0.
                    for x in range(1,npvar.size):
                        if 350. < npvar[x-1] < 360. and npvar[x] < 10.:
                            year += 1.
                        npvar[x] += 360.*year
                if avgtime:
                    varnew = np.mean(npvar.reshape(-1,5),axis=1)
                if not avgtime and not avgtod:
                    ls_start = npvar[0]
                    ls_end   = npvar[-1]
                    step     = (ls_end-ls_start)/np.float32(((numt-1)*ntod.size))
                    varnew = np.arange(0,numt*ntod.size,dtype=np.float32)
                    varnew[:] = varnew[:]*step+ls_start
#                print(npvar,varnew)
                newf.log_axis1D('areo',varnew,dims,longname_txt='time',unit_txt='',cart_txt='')
                time0=np.arange(0,len(varnew))
                newf.log_axis1D('time',varnew,dims,longname_txt='sol number',unit_txt='',cart_txt='')#added AK
            else:
                continue
        elif ndims == 4:
            varnew = npvar
            if avgtime:
                varnew = np.mean(npvar.reshape(-1,5,vshape[1],vshape[2],vshape[3]),axis=1)
            if avgtod:
                varnew = varnew.mean(axis=1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1,vshape[2],vshape[3])
            vname2=change_vname(vname)
            #AK convert surface pressure from mbar to Pa
            if vname2=='ps':varnew*=100.
            newf.log_variable(vname2,varnew,newdims,longname_txt=vname,unit_txt='')
        elif ndims == 5:
            varnew = npvar
            if avgtime:
                varnew = np.mean(npvar.reshape(-1,5,vshape[1],vshape[2],vshape[3],vshape[4]),axis=1)
            if avgtod:
                varnew = varnew.mean(axis=1)
            if not avgtime and not avgtod:
                varnew = npvar.reshape(-1,vshape[2],vshape[3],vshape[4])
            vname2=change_vname(vname)
            newf.log_variable(vname2,varnew,newdims,longname_txt=vname,unit_txt='')
        elif vname == 'tloc':
            if avgtime and not avgtod:
                newf.log_variable(vname,npvar,newdims,longname_txt=vname,unit_txt='')
    return 0


def change_vname(vname):
    vname2=vname
    if vname == 'psurf':
        vname2 = 'ps'
    if vname == 'tsurf':
        vname2 = 'ts'
    if vname.endswith('_mass'):
        if 'core' in vname:
            vname2 = 'cor_mass_micro'
        else:
            vname2 = vname[-8:-5]+'_mass_micro'
    if vname.endswith('_numb'):
        vname2 = vname[-8:-5]+'_num_micro'
    return vname2

#Function to replace dimensions with fv3 names and remove tod
def replace_dims(dims,todflag):
    newdims = dims
    if 'nlat' in dims:
        newdims = replace_at_index(newdims,newdims.index('nlat'),'lat')
    if 'nlon' in dims:
        newdims = replace_at_index(newdims,newdims.index('nlon'),'lon')
    if 'nlay' in dims:
        newdims = replace_at_index(newdims,newdims.index('nlay'),'pfull')
    if 'ntod' in dims:
        if todflag:
            newdims = replace_at_index(newdims,newdims.index('ntod'),None)
        else:
            newdims = replace_at_index(newdims,newdims.index('ntod'),'time_of_day_16')
    return newdims

def replace_at_index(tup, ix, val):
    if val is None:
        return tup[:ix]+tup[ix+1:]
    else:
        return tup[:ix] + (val,) + tup[ix+1:]


def ls2sol_1year(Ls_deg,offset=True,round10=True):
    '''
    Returns a sol number from the solar longitude.
    Args:
        Ls_deg: solar longitude in degree
        offset : if True, make year starts at Ls 0
        round10 : if True, round to the nearest 10 sols
    Returns:
        Ds :sol number
    ***NOTE***
    For the moment this is consistent with Ls 0->359.99, not for monotically increasing Ls
    '''
    Lsp=250.99   #Ls at perihelion
    tperi=485.35 #Time (in sols) at perihelion
    Ns=668.6     #Number of soils in 1 MY
    e=0.093379   #from GCM: modules.f90
    nu=(Ls_deg-Lsp)*np.pi/180
    E=2*np.arctan(np.tan(nu/2)*np.sqrt((1-e)/(1+e)))
    M=E-e*np.sin(E)
    Ds= M/(2*np.pi)*Ns+tperi
    #====Offset correction======
    if offset:
        #Ds is a float
        if len(np.atleast_1d(Ds))==1:
            Ds-=Ns
            if Ds<0:Ds+=Ns
        #Ds is an array
        else:
            Ds-=Ns
            Ds[Ds<0]=Ds[Ds<0]+Ns
    if round: Ds=np.round(Ds,-1)  #-1 means round to the nearest 10      
    return Ds

if __name__ == "__main__":
    main()
