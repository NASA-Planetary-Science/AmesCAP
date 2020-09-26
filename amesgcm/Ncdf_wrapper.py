import numpy as np
from netCDF4 import Dataset,MFDataset
#========================================================================= 
#=============Wrapper for creation of netcdf files========================
#=========================================================================

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
    def __init__(self,filename=None,description_txt="",action='w',ncformat='NETCDF4_CLASSIC'):
        if filename: 
            if filename[-3:]!=".nc":
            #assume that only path is provided so make a name for the file
                import datetime;now = datetime.datetime.now()
                filename=filename+\
                '/run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
        else:   #create a default file name  if path and filename are not provided
            import os #use a default path if not provided
            pathname=os.getcwd()+'/'
            import datetime;now = datetime.datetime.now()
            filename=pathname+\
            'run_%02d-%02d-%04d_%i-%i-%i.nc'%(now.day,now.month,now.year,now.hour,now.minute,now.second)
        self.filename=filename
        from netCDF4 import Dataset 
        if action=='w':
            self.f_Ncdf = Dataset(filename, 'w', format=ncformat)
            self.f_Ncdf.description = description_txt
        elif action=='a': #append to file
            self.f_Ncdf = Dataset(filename, 'a', format=ncformat)
        #create dictionaries to hold dimensions and variables
        self.dim_dict=dict()
        self.var_dict=dict()
        #print(filename+ " was created")
        
    def close(self):
        self.f_Ncdf.close()
        print(self.filename+" was created")
        
    def add_dimension(self,dimension_name,length):
        self.dim_dict[dimension_name]= self.f_Ncdf.createDimension(dimension_name,length)
        
    def print_dimensions(self):
        print(self.dim_dict.items())
    def print_variables(self):
        print(self.var_dict.keys())    
        
    def add_constant(self,variable_name,value,longname_txt="",units_txt=""):
        if'constant' not in self.dim_dict.keys():self.add_dimension('constant',1)
        longname_txt =longname_txt+' (%g)'%(value)   #add the value to the longname
        self._def_variable(variable_name,('constant'),longname_txt,units_txt)
        self.var_dict[variable_name][:]=value
    #=====Private definitions=====   
    def _def_variable(self,variable_name,dim_array,longname_txt="",units_txt=""):
        self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f4',dim_array)    
        self.var_dict[variable_name].units=units_txt
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].dim_name=str(dim_array)  

    def _def_axis1D(self,variable_name,dim_array,longname_txt="",units_txt="",cart_txt=""):
        self.var_dict[variable_name]= self.f_Ncdf.createVariable(variable_name,'f4',dim_array)
        self.var_dict[variable_name].units=units_txt
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].cartesian_axis=cart_txt
        
    def _test_var_dimensions(self,Ncvar):
        all_dim_OK=True
        for s in Ncvar.dimensions:
            if s not in self.dim_dict.keys():
                print("""***Warning***, dimension '"""+s+"""' not yet defined, skipping variable '"""+Ncvar._name+"""'""" )
                all_dim_OK=False
        return all_dim_OK     
    #The cartesian axis attribute is replicated if cartesian_axis is present in the original variables and if the dimensions is 1.
    #This will exclude FV3' T-cell latitudes grid_xt_bnds and grid_yt_bnds, which are of size  ('lon', 'bnds')  ('lat', 'bnds') (dim=2)
    def _is_cart_axis(self,Ncvar):
        cart_axis=False
        tmp_cart= getattr(Ncvar,'cartesian_axis',False)
        tmp_size= getattr(Ncvar,'dimensions')
        if tmp_cart and  len(tmp_size)==1: cart_axis=True
        return cart_axis     
    #================================    
    #Example: Log.log_variable('TG',TG,('time','Nx'),'soil temperature','K')
    def log_variable(self,variable_name,DATAin,dim_array,longname_txt="",units_txt=""):
        if variable_name not in self.var_dict.keys():
            self._def_variable(variable_name,dim_array,longname_txt,units_txt)
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].dim_name=str(dim_array)  
        self.var_dict[variable_name].units=units_txt
        self.var_dict[variable_name][:]=DATAin 
        
    #Example: Log.add_dim_with_content('lon',lon_array,'longitudes','degree','X')
    def log_axis1D(self,variable_name,DATAin,dim_name,longname_txt="",units_txt="",cart_txt=""):
        if variable_name not in self.var_dict.keys():
            self._def_axis1D(variable_name,dim_name,longname_txt,units_txt,cart_txt)
        self.var_dict[variable_name].long_name=longname_txt
        self.var_dict[variable_name].units=units_txt
        self.var_dict[variable_name].cartesian_axis=cart_txt
        self.var_dict[variable_name][:]=DATAin
        
    #Function to define a dimension and add a variable with at the same time.
    #Equivalent to add_dimension(), followed by  log_axis1D()
    #lon_array=np.linspace(0,360)
    #Example: Log.add_dim_with_content('lon',lon_array,'longitudes','degree','X')
    def add_dim_with_content(self,dimension_name,DATAin,longname_txt="",units_txt="",cart_txt=''):
        if dimension_name not in self.dim_dict.keys():self.add_dimension(dimension_name,len(DATAin))
        #---If no longname is provided, simply use dimension_name as default longname---
        if longname_txt=="":longname_txt=dimension_name
        if dimension_name not in self.var_dict.keys():
            self._def_axis1D(dimension_name,dimension_name,longname_txt,units_txt,cart_txt)
        self.var_dict[dimension_name].long_name=longname_txt
        self.var_dict[dimension_name].units=units_txt
        self.var_dict[dimension_name].cartesian_axis=cart_txt
        self.var_dict[dimension_name][:]=DATAin 
    
    
    #***Note***
    #The attribute 'name'  was replaced by '_name' to allow compatibility with MFDataset:
    #When using f=MFDataset(fname,'r'), f.variables[var] does not have a 'name' attribute but does have '_name'
    #________     
    #Copy a netcdf DIMENSION variable e.g Ncdim is:  f.variables['lon']
    # if the dimension for that variable does not exist yet, it will be created
    def copy_Ncaxis_with_content(self,Ncdim_var):
        longname_txt=getattr(Ncdim_var,'long_name',Ncdim_var._name)
        units_txt=    getattr(Ncdim_var,'units','')
        cart_txt=    getattr(Ncdim_var,'cartesian_axis','')
        self.add_dim_with_content(Ncdim_var._name,Ncdim_var[:],longname_txt,units_txt,cart_txt)
               
    #Copy a netcdf variable from another file, e.g Ncvar is: f.variables['ucomp']
    #All dimensions must already exist. If swap_array is provided, the original values will be 
    #swapped with this array.
    def copy_Ncvar(self,Ncvar,swap_array=None):
        if Ncvar._name not in self.var_dict.keys():
            dim_array=Ncvar.dimensions
            longname_txt=getattr(Ncvar,'long_name',Ncvar._name)
            units_txt=    getattr(Ncvar,'units','')
            self._def_variable(Ncvar._name,Ncvar.dimensions,longname_txt,units_txt)
            if np.any(swap_array):
                self.log_variable(Ncvar._name,swap_array[:],Ncvar.dimensions,longname_txt,units_txt)
            else:   
                self.log_variable(Ncvar._name,Ncvar[:],Ncvar.dimensions,longname_txt,units_txt) 
        else:    
            print("""***Warning***, '"""+Ncvar._name+"""' is already defined, skipping it"""  )
        
    #Copy all variables, dimensions and attributes from another Netcdfile.   
    def copy_all_dims_from_Ncfile(self,Ncfile_in,exclude_dim=[],time_unlimited=True):
        #----First include dimensions-------
        all_dims=Ncfile_in.dimensions.keys()
        for idim in all_dims:
            if idim not in exclude_dim:
                if idim=='time' and time_unlimited:
                    self.add_dimension(Ncfile_in.dimensions[idim]._name,None) 
                else:
                    self.add_dimension(Ncfile_in.dimensions[idim]._name,Ncfile_in.dimensions[idim].size)

    def copy_all_vars_from_Ncfile(self,Ncfile_in,exclude_var=[]):
        #----First include variables-------             
        all_vars=Ncfile_in.variables.keys()
        for ivar in all_vars:
            if ivar not in exclude_var:
                #Test if all dimensions are availalbe, skip variables otherwise
                if self._test_var_dimensions(Ncfile_in.variables[ivar]):
                    if self._is_cart_axis(Ncfile_in.variables[ivar]):
                        self.copy_Ncaxis_with_content(Ncfile_in.variables[ivar])
                    else:
                        self.copy_Ncvar(Ncfile_in.variables[ivar])
    
    def merge_files_from_list(self,Ncfilename_list,exclude_var=[]):
        Mf_IN=MFDataset(Ncfilename_list,'r')
        self.copy_all_dims_from_Ncfile(Mf_IN)
        self.copy_all_vars_from_Ncfile(Mf_IN,exclude_var=exclude_var)
        Mf_IN.close()
    
##=====TEST ONLY=======
'''
#fname='/Users/akling/test/00010.atmos_average.nc'
#fname='/Users/akling/test/00010.fixed.nc'
fname='/Users/akling/Data/FV3/C24_L36_CG_drag_ON/03335.atmos_average.nc'

#f=Dataset(fname,'r')
f=MFDataset(fname,'r')
varnames=f.variables.keys()
lat=f.variables['lat']
t=f.variables['time']
areo=f.variables['areo']
pfull=f.variables['pfull']
lon=f.variables['lon'][:]  
ucomp=f.variables['ucomp']  
    
f.variables['ucomp']

# for name, variable in f.variables.items():   
#     print(name);#print(variable)
#     for attrname in variable.ncattrs():
#         print("{} -- {}".format(attrname, getattr(variable, attrname)))

Log=Ncdf('/Users/akling')
Log.copy_all_dims_from_Ncfile(f)
#Log.copy_Ncvar(f.variables['ucomp'])
#Log.copy_Ncvar(f.variables['lat'])
Log.copy_all_vars_from_Ncfile(f)
#Fgeo= 0.03 #W/m2, a constant

#Log.add_dimension('Nx',8)
# Log.add_dim_with_content('lon',lon,'longitudes','degree')
# Log.copy_Ncdim_with_content(f.variables['lat'])
# Log.copy_Ncdim_with_content(f.variables['pfull'])
# Log.copy_Ncdim_with_content(f.variables['time'])
# Log.copy_Ncvar(f.variables['ucomp']  )

#Log.close()  
'''
  