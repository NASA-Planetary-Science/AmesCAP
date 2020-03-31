#!/usr/bin/env python

#
#Assumes the companion scripts are in the same directory
#Load generic Python Modules
import argparse #parse arguments
import sys
import getopt
import os
import glob
import shutil
import subprocess
import numpy as np
from netCDF4 import Dataset
#from amesgcm.lib.FV3_utils import find_n
sys.path.append(os.getcwd())

#======================================================
#                  ARGUMENTS PARSER
#======================================================
parser = argparse.ArgumentParser(description="""\033[93m MarsFiles files manager. Used to convert Legacy GCM to FV3 format \n \033[00m""",
								formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+', 
								help='***.nc file or list of ***.nc files ')
								
parser.add_argument('-fv3','--fv3', nargs='+',default=[],
					help="""Produce FV3's diurn, average and daily files from Legacy outputs  \n"""
						"""> Usage: MarsFiles LegacyGCM*.nc -fv3 \n"""
						"""\033[00m""")

parser.add_argument('-c','--combine', nargs='+',default=['p'],
					help="""Combine a sequence of similar files as one file\n"""
						"""> Usage: MarsVars *.atmos_average.nc --combine \n"""
						""" \n""")

parser.add_argument('--debug',  action='store_true', help='Debug flag: release the exceptions')


def main(argv):

	#Test if ncks is available , make recommendation and  exit otherwise--	    
	try:
        	subprocess.check_call('ncks --version',shell=True,stdout=open(os.devnull, "w"), stderr=open(os.devnull, "w"))
	except subprocess.CalledProcessError:
		prYellow("ncks dependency is needed to concatenate files. On NAS, use:")
		prCyan('    module load nco')
		#exit()


	print("====== Running postprocess_legacy.py =======")
	
	file_list=parser.parse_args().input_file

#argument definitions:
	cwd=os.getcwd()
	path2data=os.getcwd()
	do_multi= False
	do_1year= False #Used with LegacyGCM_1year.nc'
	
	#Get files to process
	histlist=[]
	for filei in file_list:
		histlist.append(path2data+'/'+filei)
	prCyan(histlist)
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
			a=make_atmos_avg(f,cwd)

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
		a=make_atmos_avg(hist1year,cwd)

	



def make_atmos_avg(fpath,cwd=None):	
	import re
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
			a=add_dim(histfile,newfavg,'lon',npvar.size,npvar,getattr(var,'units'))
			a=add_dim(histfile,newfdaily,'lon',npvar.size,npvar,getattr(var,'units'))
			a=add_dim(histfile,newfdiurn,'lon',npvar.size,npvar,getattr(var,'units'))
		if dname == 'nlat':
			var=histfile.variables['latitude']
			npvar=var[:]
			a=add_dim(histfile,newfavg,'lat',npvar.size,npvar,getattr(var,'units'))
			a=add_dim(histfile,newfdaily,'lat',npvar.size,npvar,getattr(var,'units'))
			a=add_dim(histfile,newfdiurn,'lat',npvar.size,npvar,getattr(var,'units'))
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
			pref=7.01 * 100  #AK changed this to Pa
			pk=np.zeros(nump)
			bk=np.zeros(nump)
			pfull=np.zeros(num)
			phalf=np.zeros(nump)
			dsgm=histfile.variables['dsgm']
			sgm =histfile.variables['sgm']
			pk[0]=.08
			for z in range(num):
				bk[z+1] = sgm[2*z+2]
			phalf[:]=pk[:]+pref*bk[:] #AK changed this to Pa
			pfull[:] = (phalf[1:]-phalf[:num])/(np.log(phalf[1:])-np.log(phalf[:num]))
			a=add_dim(histfile,newfavg,'pfull',num,pfull,'Pa')
			a=add_dim(histfile,newfavg,'phalf',nump,phalf,'Pa')
			newfavg.log_var1d('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pa',cart_txt='')
			newfavg.log_var1d('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='Pa',cart_txt='')
			a=add_dim(histfile,newfdaily,'pfull',num,pfull,'Pa')
			a=add_dim(histfile,newfdaily,'phalf',nump,phalf,'Pa')
			newfdaily.log_var1d('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pascal',cart_txt='')
			newfdaily.log_var1d('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='',cart_txt='')
			a=add_dim(histfile,newfdiurn,'pfull',num,pfull,'Pa')
			a=add_dim(histfile,newfdiurn,'phalf',nump,phalf,'Pa')
			newfdiurn.log_var1d('pk',pk,('phalf'),longname_txt='pressure part of the hybrid coordinate',unit_txt='Pascal',cart_txt='')
			newfdiurn.log_var1d('bk',bk,('phalf'),longname_txt='sigma part of the hybrid coordinate',unit_txt='',cart_txt='')

#perform various averages
	a=do_avg_vars(histfile,newfavg,True,True)
	a=do_avg_vars(histfile,newfdaily,False,False)
	a=do_avg_vars(histfile,newfdiurn,True,False)
	
	
	
	newfavg.close()
	newfdaily.close()
	newfdiurn.close()
	
	
	return 0


#Function to perform time averages over all fields
def do_avg_vars(histfile,newf,avgtime,avgtod):
	histvars = histfile.variables.keys()
	for vname in histvars:
		var	 = histfile.variables[vname]
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
#				print(npvar,varnew)
				newf.log_var1d('areo',varnew,dims,longname_txt='time',unit_txt='',cart_txt='')
				time0=np.arange(0,len(varnew))
				newf.log_var1d('time',varnew,dims,longname_txt='sol number',unit_txt='',cart_txt='')#added AK
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
				
			


#Function to add and define dimension	
def add_dim(ncfile,newf,dname,dsize,dvals,varunit):
	newf.add_dimension(dname,dsize)
	dims=(dname)
	newf.log_var1d(dname,dvals,dims,longname_txt=dname,unit_txt=varunit)
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

	class Ncdf(object):
		'''
		Alex K.
		NetCdf wrapper for quick archiving of data into netcdf format
	
		USAGE:
	
		from netcdf_wrapper import Ncdf
	
		Fgeo= 0.03 #W/m2, a constant
		TG=np.ones((24,8)) #ground temperature
	
		#---create file---
		filename="/lou/s2n/mkahre/MCMC/analysis/working/myfile.nc"
		description="results from new simulation, Alex 01-01-19"
		Log=Ncdf(filename,description)
	
		#---Save the constant to the file---
		Log.add_constant('Fgeo',Fgeo,"geothermal flux","W/m2")
	
		#---Save the TG array to the file---
		Log.add_dimension('Nx',8)
		Log.add_dimension('time',24)
	
		Log.log_variable('TG',TG,('time','Nx'),'soil temperature','K')
	
		Log.close()
	
	
		'''
		def __init__(self,filename=None,description_txt="",action='w'):
			if filename:
				if filename[-3:]!=".nc":
				#assume that only path is provided so make a name for the file
					import datetime;now = datetime.datetime.now()
					filename=filename+\
					'/run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
			else:	#create a default file name	 if path and filename are not provided
				import os #use a default path if not provided
				pathname=os.getcwd()+'/'
				import datetime;now = datetime.datetime.now()
				filename=pathname+\
				'run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
			self.filename=filename
			from netCDF4 import Dataset
			if action=='w':
				self.f_Ncdf = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
				self.f_Ncdf.description = description_txt
			elif action=='a': #append to file
				self.f_Ncdf = Dataset(filename, 'a', format='NETCDF4_CLASSIC')
			#create dictionaries to hold dimensions and variables
			self.dim_dict=dict()
			self.var_dict=dict()
			print(filename+ " was created")
	
		def close(self):
			self.f_Ncdf.close()
			print(self.filename+" was closed")
	
		def add_dimension(self,dimension_name,length):
			self.dim_dict[dimension_name]= self.f_Ncdf.createDimension(dimension_name,length)
	
		def print_dimension(self):
			print(self.dim_dict.items())
		def print_variable(self):
			print(self.var_dict.keys())
	
		def add_constant(self,variable_name,value,longname_txt="",unit_txt=""):
			if not any('constant' in s for s in self.dim_dict.keys()):
				self.add_dimension('constant',1)
			longname_txt =longname_txt+' (%g)'%(value)	 #add the value to the longname
			self.def_variable(variable_name,('constant'),longname_txt,unit_txt)
			self.var_dict[variable_name][:]=value
	
		def def_variable(self,variable_name,dim_array,longname_txt="",unit_txt=""):
			self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f4',dim_array)
			self.var_dict[variable_name].units=unit_txt
			self.var_dict[variable_name].long_name=longname_txt
	#		self.var_dict[variable_name].checksum=chk_txt
		def log_variable(self,variable_name,DATAin,dim_array,longname_txt="",unit_txt="",):
			if not any(variable_name == s for s in self.var_dict.keys()):
				self.def_variable(variable_name,dim_array,longname_txt,unit_txt)
			self.var_dict[variable_name].long_name=longname_txt
			self.var_dict[variable_name].units=unit_txt
			self.var_dict[variable_name][:]=DATAin
	
		def def_var1d(self,variable_name,dim_array,longname_txt="",unit_txt="",cart_txt=""):
			self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f8',dim_array)
			self.var_dict[variable_name].units=unit_txt
			self.var_dict[variable_name].long_name=longname_txt
			self.var_dict[variable_name].cartesian_axis=cart_txt
		def log_var1d(self,variable_name,DATAin,dim_array,longname_txt="",unit_txt="",cart_txt="",):
			if not any(variable_name == s for s in self.var_dict.keys()):
				self.def_var1d(variable_name,dim_array,longname_txt,unit_txt,cart_txt)
			self.var_dict[variable_name].long_name=longname_txt
			self.var_dict[variable_name].units=unit_txt
			self.var_dict[variable_name].cartesian_axis=cart_txt
			self.var_dict[variable_name][:]=DATAin

#AK : Function for output
def prCyan(skk): print("\033[96m{}\033[00m" .format(skk))
def prYellow(skk): print("\033[93m{}\033[00m" .format(skk))

if __name__ == "__main__":
	main(sys.argv[1:])
