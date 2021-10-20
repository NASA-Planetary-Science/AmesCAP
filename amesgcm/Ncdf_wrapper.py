import numpy as np
from netCDF4 import Dataset,MFDataset
from scipy.io import FortranFile
from amesgcm.FV3_utils import daily_to_average, daily_to_diurn
import os

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

#====================================================================================== 
#====Wrapper for creation of netcdf-like object from Legacy GCM Fortran binaries=======
#======================================================================================

    
class Fort(object):
    '''
    Alex K.
    A class that generate an object from fort.11 ... with similar attributes as a Netcdf file, e.g.: 
    >>  f.variables.keys()
    >>  f.variables['var'].long_name
    >>  f.variables['var'].units
    >>  f.variables['var'].dimensions
    
    Create a Fort object using the following:
    f=Fort('/Users/akling/test/fort.11/fort.11_0684')
    
    PUBLIC METHODS:
    >> f.write_to_fixed(), f.write_to_average()  f.write_to_daily()  and f.write_to_diurn() can be used to generate FV3-like netcdf files
    '''
    
    #===Inner class for fortran_variables (Fort_var) that make up the Fort file===
    class Fort_var(np.ndarray):
        '''
        Sub-class that emulate a netcdf-like variable by adding the name, long_name, units, dimensions attribute to a numpy array. [A. Kling]
        *** NOTE***
        A useful resource on subclassing in available at:
        https://numpy.org/devdocs/reference/arrays.classes.html
        
        Note that because we use an existing numpy.ndarray to define the object, we do not use a call to __array_finalize__(self, obj) 
        '''
    
            
        def __new__(cls, input_vals,*args, **kwargs):
            return np.asarray(input_vals).view(cls)    
    
        def __init__(self, input_vals,name_txt,long_name_txt,units_txt,dimensions_tuple):
            self.name = name_txt
            self.long_name = long_name_txt
            self.units= units_txt
            self.dimensions=dimensions_tuple   
            
            
    #==== End of inner class===
    
    def __init__(self,filename=None,description_txt=""):
        from scipy.io import FortranFile
        self.filename=filename
        self.path,self.name=os.path.split(filename) 
        print('Reading '+filename + ' ...')
        self.f = FortranFile(filename)
        if len(filename)==12:
            self.fort_type=filename[-7:-5] #Get output number, e.g. 11 for fort.11_0070, 45 for fort.45_0070 etc..
        else:
            #Case if file is simply named 'fort.11'  which is the case for the first file of a cold  start
            self.fort_type=filename[-2:]
            
        self.nperday=16  # TODO Hard-coded: 16 outputs per day
        self.nsolfile=10 # TODO Hard-coded: 10 sols per output
        #Add time of day dimensions
        self.tod_name=tod_name='time_of_day_%02d'%(self.nperday)    
        self.tod=np.arange(0.5*24/self.nperday,24,24/self.nperday)  # i.e np.arange(0.75,24,1.5) every 1.5 hours, centered at half timestep =0.75
        
        self.dimensions={} #Initialize dictionary        
        self.variables={} #Initialize dictionary

        if self.fort_type=='11':
            self._read_Fort11_header()
            self._read_Fort11_constants()
            self._read_Fort11_static()
            self._create_dims()
            self._read_Fort11_dynamic()
            self._add_axis_as_variables()
            #TODO monotically increasing MY: Get date as FV3 file e.g. 00000
            self.fdate="%05i"%self._ls2sol_1year(self.variables['areo'][0])
            
    #Public methods       
    def write_to_fixed(self):
        '''
        Create 'fixed' file, i.e.  all static variables
        '''
        
        Log=Ncdf(self.path+'/'+self.fdate+'.fixed.nc')
        
        #Define dimensions
        for ivar in ['lat','lon','pfull','phalf','zgrid']:
            if ivar =='lon':cart_ax='X'
            if ivar =='lat':cart_ax='Y'
            if ivar in ['pfull' ,'phalf','zgrid']:cart_ax='Z'
            fort_var=self.variables[ivar]
            Log.add_dim_with_content(dimension_name=ivar,DATAin=fort_var,longname_txt=fort_var.long_name,units_txt=fort_var.units,cart_txt=cart_ax)
        
        #Log static variables
        for ivar in self.variables.keys():
            if 'time' not in self.variables[ivar].dimensions:
                fort_var=self.variables[ivar]
                Log.log_variable(variable_name=ivar,DATAin=fort_var,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        Log.close()
        
    def write_to_daily(self):
        '''
        Create daily file, e.g. contineuous time serie
        '''
        Log=Ncdf(self.path+'/'+self.fdate+'.atmos_daily.nc')
        
        #Define dimensions
        for ivar in ['lat','lon','pfull','phalf']:
            if ivar =='lon':cart_ax='X'
            if ivar =='lat':cart_ax='Y'
            if ivar in ['pfull' ,'phalf','zgrid']:cart_ax='Z'
            fort_var=self.variables[ivar]
            Log.add_dim_with_content(dimension_name=ivar,DATAin=fort_var,longname_txt=fort_var.long_name,units_txt=fort_var.units,cart_txt=cart_ax)
        
        #Add scalar_axis dimension (size 1, only used with areo) 
        Log.add_dimension('scalar_axis',1)    
            
        #Add aggregation dimension (None size for unlimited) 
        Log.add_dimension('time',None)
        fort_var=self.variables['time']
        Log.log_axis1D(variable_name='time',DATAin=fort_var,dim_name='time',longname_txt=fort_var.long_name,units_txt=fort_var.units,cart_txt='T')
        
        #Special case for the solar longitude (areo): needs to be interpolated linearly every 16 timesteps
        ivar='areo';fort_var=self.variables[ivar]
        var_out=self._linInterpLs(np.squeeze(fort_var[:]),16).reshape([len(fort_var),1]) #areo is reshaped as [time,scalar_axis]=[160,1]
        Log.log_variable(variable_name=ivar,DATAin=var_out,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        
        #Log dynamic variables, as well as pk, bk
        for ivar in self.variables.keys():
            if 'time' in self.variables[ivar].dimensions and ivar!='areo' or ivar in ['pk','bk']:
                fort_var=self.variables[ivar]
                Log.log_variable(variable_name=ivar,DATAin=fort_var,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        Log.close()

    def write_to_average(self,day_average=5):
        '''
        Create average file, e.g. N day averages (typically 5)
        '''
        Log=Ncdf(self.path+'/'+self.fdate+'.atmos_average.nc')
        #Define dimensions
        for ivar in ['lat','lon','pfull','phalf']:
            if ivar =='lon':cart_ax='X'
            if ivar =='lat':cart_ax='Y'
            if ivar in ['pfull' ,'phalf','zgrid']:cart_ax='Z'
            fort_var=self.variables[ivar]
            Log.add_dim_with_content(dimension_name=ivar,DATAin=fort_var,longname_txt=fort_var.long_name,units_txt=fort_var.units,cart_txt=cart_ax)
        
        #Add scalar_axis dimension (size 1, only used with areo) 
        Log.add_dimension('scalar_axis',1)
        
        #Add aggregation dimension (None size for unlimited) 
        Log.add_dimension('time',None)
        
        #Perform day average and log new time axis
        time_in=self.variables['time']
        time_out=daily_to_average(varIN=fort_var,dt_in=time_in[1]-time_in[0],nday=day_average,trim=True)
        Log.log_axis1D(variable_name='time',DATAin=time_out,dim_name='time',longname_txt=time_in.long_name,units_txt=time_in.units,cart_txt='T')

        #Log static variables 
        for ivar in ['pk','bk']:
            fort_var=self.variables[ivar]
            Log.log_variable(variable_name=ivar,DATAin=fort_var,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        
        #Log dynamic variables 
        for ivar in self.variables.keys():
            if 'time' in self.variables[ivar].dimensions:
                fort_var=self.variables[ivar]
                var_out=daily_to_average(fort_var,time_in[1]-time_in[0],nday=day_average,trim=True)
                Log.log_variable(variable_name=ivar,DATAin=var_out,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        
        Log.close()
        

            
    def write_to_diurn(self,day_average=5):
        '''
        Create diurn file, e.g.variable are organized by time of day. Additionally, the data is also binned  (typically 5)
        '''
        Log=Ncdf(self.path+'/'+self.fdate+'.atmos_diurn.nc')
        #Define dimensions
        for ivar in ['lat','lon','pfull','phalf']:
            if ivar =='lon':cart_ax='X'
            if ivar =='lat':cart_ax='Y'
            if ivar in ['pfull' ,'phalf','zgrid']:cart_ax='Z'
            fort_var=self.variables[ivar]
            Log.add_dim_with_content(dimension_name=ivar,DATAin=fort_var,longname_txt=fort_var.long_name,units_txt=fort_var.units,cart_txt=cart_ax)
        
        #Add scalar_axis dimension (size 1, only used with areo) 
        Log.add_dimension('scalar_axis',1)
        
        #Add time_of_day dimensions
        Log.add_dim_with_content(dimension_name=self.tod_name,DATAin=self.tod,longname_txt='time of day',units_txt='hours since 0000-00-00 00:00:00',cart_txt='N')
        
        #Add aggregation dimension (None size for unlimited) 
        Log.add_dimension('time',None)
        
        #Perform day average and log new time axis
        time_in=self.variables['time']
        time_out=daily_to_average(varIN=fort_var,dt_in=time_in[1]-time_in[0],nday=day_average,trim=True)
        Log.log_axis1D(variable_name='time',DATAin=time_out,dim_name='time',longname_txt=time_in.long_name,units_txt=time_in.units,cart_txt='T')
        
        #Log static variables 
        for ivar in ['pk','bk']:
            fort_var=self.variables[ivar]
            Log.log_variable(variable_name=ivar,DATAin=fort_var,dim_array=fort_var.dimensions,longname_txt=fort_var.long_name,units_txt=fort_var.units)
        
        #Loop over all variables in file
        for ivar in self.variables.keys():
            if 'time' in self.variables[ivar].dimensions :
                fort_var=self.variables[ivar]
                #If time is the dimension (but not just a time array)
                if 'time' in fort_var.dimensions and ivar!='time':   
                    dims_in=fort_var.dimensions
                    if type(dims_in)==str: #dimensions has 'time' only, it is a string
                        dims_out=(dims_in,)+(self.tod_name,)
                    else: #dimensions is a tuple, e.g. ('time','lat','lon')
                        dims_out=(dims_in[0],)+(self.tod_name,)+dims_in[1:]
                            
                    var_out=daily_to_diurn(fort_var[:],time_in[0:self.nperday])
                    if day_average!=1:var_out=daily_to_average(var_out,1.,day_average) #dt is 1 sol between two diurn timestep
                    Log.log_variable(ivar,var_out,dims_out,fort_var.long_name,fort_var.units)
                   
        Log.close()
        
    
    #Public method       
    def close(self):
        self.f.close()
        print(self.filename+" was closed")
    #Private methods    

    def _read_Fort11_header(self):
        '''
        Return values from fort.11 header:
        Args:
            f: an opened scipy.io.FortranFile object
        Return: RUNNUM, JM,IM,LM,NL,ntrace,version and SM
        ***NOTE***    
        In myhist.f: 
        write(11) RUNNUM (float), JM, IM, LAYERS, NL, NTRACE (ints), version (char= 7)
        
         >> Those are save as attributes, e.g. used f.LAYERS to access the number of layers
        '''
        Rec=self.f.read_record('f4','(1,5)i4','S7')
        self.RUNNUM=Rec[0][0];self.JM=Rec[1][0,0]; self.IM=Rec[1][0,1]; self.LM=Rec[1][0,2];self.NL=Rec[1][0,3]; self.ntrace=Rec[1][0,4];self.version=Rec[2][0]
        
        #Also compute subsurface grid (boundaries)
        self.SM=2*self.NL+1

    def _read_Fort11_constants(self):
        '''
        Return run constants from fort.11 header:
        Args:
            f: an opened scipy.io.FortranFile object
        Return:
        ***NOTE***    
        In myhist.f: 

        write(11) DSIG, DXYP, GRAV, RGAS, cp, stbo, xlhtc, kapa, 
        *          cmk, decmax, eccn, orbinc, vinc, sdepth, alicen,
        *          alices, egoco2n, egoco2s, npcwikg, gidn, gids
        
        >> Those are save as attributes, e.g. use f.rgas to access the gas constant for the simulation
        '''
        Rec=self.f.read_record('(1,{0})f4'.format(self.LM),'(1,{0})f4'.format(self.JM),'f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','(1,{0})f4'.format(self.SM),'f4','f4','f4','f4','f4','f4','f4')
        self.dsig=np.array(Rec[0][0,:]);self.dxyp=np.array(Rec[1][0,:]);self.grav=Rec[2][0]
        self.rgas=Rec[3][0];self.cp=Rec[4][0];self.stbo=Rec[5][0];self.xlhtc=Rec[6][0];self.kapa=Rec[7][0]
        self.cmk=Rec[8][0];self.decmax=Rec[9][0];self.eccn=Rec[10][0];self.orbinc=Rec[11][0];self.vinc=Rec[12][0]
        self.sdepth=np.array(Rec[13][0,:]);self.alicen=Rec[14][0];
        self.alices=Rec[15][0];self.egoco2n=Rec[16][0];self.egoco2s=Rec[17][0]
        
    def _read_Fort11_static(self):
        '''
        Return values from fort.11 header:
        Args:
            f: an opened scipy.io.FortranFile object
        Return:
        ***NOTE***    
        In myhist.f: 

        write(11) TOPOG, ALSP, ZIN, NPCFLAG

        >> Those are save as variables, e.g. used f.variables['zsurf'] to access the topography
        
        '''
        Rec=self.f.read_record('({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM),'({0},{1},{2})f4'.format(self.NL,self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM))
        
        
        #Add static variables to the 'variables' dictionary. 

        self.variables['zsurf']=   self.Fort_var(-np.array(Rec[0].T/self.grav),'zsurf','surface height','m',('lat', 'lon'))
        self.variables['alb']=     self.Fort_var(np.array(Rec[1].T),'alb','Surface Albedo','mks',('lat', 'lon'))     
        self.variables['thin']=    self.Fort_var(np.array(Rec[2].transpose([0,2,1])),'thin','Surface Thermal Inertia','J/m2/K/s1/2',('zgrid','lat', 'lon')) 
        self.variables['npcflag']= self.Fort_var(np.array(Rec[3].T),'npcflag','Polar ice flag','none',('lat', 'lon'))    

    def _create_dims(self):
        '''
        Create dimensions axis from IM, JM after reading the header. Also compute vertical grid structure that includes
        sigma values at the layers' boundaries AND  at layers' midpoints for the radiation code. Total size is therefore 2*LM+2
        '''
        JM=self.JM;IM=self.IM;LM=self.LM;NL=self.NL
        self.lat = -90.0 + (180.0/JM)*np.arange(1,JM+1)
        self.lon=-180.+(360./IM)*np.arange(1,IM+1) 
        
        #compute sigma layer. Temporary arrays:
        sigK=np.zeros(2*LM+3)   #Layer midpoints + boundaries
        sigL=np.zeros(LM)       #Layer midpoints only
        
        #these are edges and midpoints for output
        self.sigm=np.zeros(2*LM+1)
        
        sigK[0:3]=0.
        for l in range(0,LM):
            k=2*(l)+3
            sigK[k] = sigK[k-2]+self.dsig[l]  
        for k in range(4,2*LM+3-1,2):
            sigK[k] = 0.5*(sigK[k+1]+sigK[k-1])
        for l in range(0,LM):
            sigL[l] = sigK[l*2+1]
        sigK[2*LM+2]=1.0
        self.sigm[:]=sigK[2:]
        
        # Subsurface layer
        
        # Assume this is bound |midpoint bound so we take every other point starting with the 2nd
        
        self.zgrid = self.sdepth[1::2]    #TODO check
        
    def _ra_1D(self,new_array,name_txt):
        '''
        _ra stands for 'Return array': Append single timesteps along the first (time) dimensions
        '''
        if type(new_array)!=np.ndarray:new_array=np.array([new_array])
        
        #Add time axis to new data e.g. turn [lat,lon]to as [1,lat,lon]
        new_shape=np.append([1],new_array.shape)
        #First time that varialbe is encountered
        if name_txt not in self.variables.keys():
            return new_array
        else: #Get values from existing array and append to it. Note that np.concatenate((x,y)) takes a tuple as argument.
            return np.append(self.variables[name_txt],new_array)

    def _ra_2D(self,name_txt,Rec=None):
        '''
        _ra stands for 'Return array': Append single timesteps along the first (time) dimensions
        Args:
            name_txt : char, name of variables, e.g. 'temp'
            Rec: Record to archive. if None, read directly in fortran binary.  
        
        '''
        if Rec is None:
            Rec=self.f.read_reals('f4').reshape(self.JM,self.IM, order='F')
        #Set to pole point to value at N-1
        Rec[-1,...]=Rec[-2,...]
        #Add time axis to new data e.g. turn [lat,lon]to as [1,lat,lon]
        new_shape=np.append([1],Rec.shape)
        #First time that varialbe is encountered
        if name_txt not in self.variables.keys():
            return Rec.reshape(new_shape)
        else: #Get values from existing array and append to it. Note that np.concatenate((x,y)) takes a tuple as argument.
            return np.concatenate((self.variables[name_txt],Rec.reshape(new_shape)))
            
    def _ra_3D_atmos(self,name_txt,Rec=None):
        '''
        _ra stands for 'Return array': Append single timesteps along the first (time) dimensions
        '''
        if Rec is None:
            Rec=self.f.read_reals('f4').reshape(self.JM,self.IM, self.LM, order='F')
        #Set to pole point to value at N-1
        Rec[-1,...]=Rec[-2,...]
        Rec=Rec.transpose([2,0,1])
        #Add time axis to new data e.g. turn [lat,lon]to as [1,lat,lon]
        new_shape=np.append([1],Rec.shape)
        #First time that varialbe is encountered
        if name_txt not in self.variables.keys():
            return Rec.reshape(new_shape)
        else: #Get values from existing array and append to it. Note that np.concatenate((x,y)) takes a tuple as argument.
            return np.concatenate((self.variables[name_txt],Rec.reshape(new_shape)))
    
    def _log_var(self,name_txt,long_name,unit_txt,dimensions,Rec=None,scaling=None):
        '''
        
        '''
        #No Record is provided, read from file
        if Rec is None:
            if dimensions==('time','lat','lon'):
                Rec=self.f.read_reals('f4').reshape(self.JM,self.IM, order='F')
            if dimensions==('time','pfull','lat','lon'):
                Rec=self.f.read_reals('f4').reshape(self.JM,self.IM, self.LM, order='F')
            if dimensions==('time','zgrid','lat','lon'):
                Rec=self.f.read_reals('f4').reshape(self.JM,self.IM, self.NL, order='F')
            #If scaling, scale it!    
            if scaling:Rec*=scaling
        #Reorganize 3D vars        
        if dimensions==('time','pfull','lat','lon') or dimensions==('time','zgrid','lat','lon'):Rec=Rec.transpose([2,0,1])
        #Set to pole point to value at N-1
        Rec[...,-1,:]=Rec[...,-2,:]
        
        #Add time axis to new data e.g. turn [lat,lon]to as [1,lat,lon]
        new_shape=np.append([1],Rec.shape)
        #First time that varialbe is encountered
        if name_txt not in self.variables.keys():
            Rec=Rec.reshape(new_shape)
        else: #Get values from existing array and append to it. Note that np.concatenate((x,y)) takes a tuple as argument.
            Rec= np.concatenate((self.variables[name_txt],Rec.reshape(new_shape)))
            
        #Log the variable    
        self.variables[name_txt]=  self.Fort_var(Rec ,name_txt,long_name,unit_txt,dimensions)
        
        
    
    def _read_Fort11_dynamic(self):
        '''
        Read variables from fort.11 files that changes with each timestep.
        
        In mhistv.f :
        
        WRITE(11) TAU, VPOUT, RSDIST, TOFDAY, PSF, PTROP, TAUTOT,
            *          RPTAU, SIND, GASP
            WRITE(11) NC3, NCYCLE
        
            WRITE(11) P
            WRITE(11) T
            WRITE(11) U
            WRITE(11) V
            WRITE(11) GT
            WRITE(11) CO2ICE
            WRITE(11) STRESSX
            WRITE(11) STRESSY
            WRITE(11) TSTRAT
            WRITE(11) TAUSURF
            WRITE(11) SSUN
            WRITE(11) QTRACE
            WRITE(11) QCOND
            write(11) STEMP
            write(11) fuptopv, fdntopv, fupsurfv, fdnsurfv
            write(11) fuptopir, fupsurfir, fdnsurfir
            write(11) surfalb
            write(11) dheat
            write(11) geot
        
        '''   
        nsteps=   self.nperday* self.nsolfile  #typically 16 x 10 =160
        append=False
        for iwsol in range(0,nsteps): 
            Rec=self.f.read_record('f4')
            #TAU=Rec[0];VPOUT=Rec[1]; RSDIST=Rec[2]; TOFDAY=Rec[3]; PSF=Rec[4]; PTROP=Rec[5]; TAUTOT=Rec[6]; RPTAU=Rec[7]; SIND=Rec[8]; GASP2=Rec[9]
            
            self.variables['time']=  self.Fort_var(self._ra_1D(Rec[0]/24,'time')  ,'time','elapsed time from the start of the run','days since 0000-00-00 00:00:00',('time'))
            self.variables['areo']= self.Fort_var(self._ra_1D(Rec[1].reshape([1,1]),'areo')     ,'areo','solar longitude','degree',('time','scalar_axis'))     #TODO monotically increasing ?
            self.variables['rdist']= self.Fort_var(self._ra_1D(Rec[2],'rdist')    ,'rdist','square of the Sun-Mars distance','(AU)**2',('time'))    
            self.variables['tofday']=self.Fort_var(self._ra_1D(Rec[3],'tofday')   ,'npcflag','time of day','hours since 0000-00-00 00:00:00',('time')) #TODO edge or center ?
            self.variables['psf']=   self.Fort_var(self._ra_1D(Rec[4]*100,'psf')  ,'psf','Initial global surface pressure','Pa',('time'))   
            self.variables['ptrop']= self.Fort_var(self._ra_1D(Rec[5],'ptrop')    ,'ptrop','pressure at the tropopause','Pa',('time'))   
            self.variables['tautot']=self.Fort_var(self._ra_1D(Rec[6],'tautot')   ,'tautot','Input (global) dust optical depth at the reference pressure','none',('time'))   
            self.variables['rptau']= self.Fort_var(self._ra_1D(Rec[7]*100,'rptau'),'rptau','reference pressure for dust optical depth','Pa',('time'))   
            self.variables['sind']=  self.Fort_var(self._ra_1D(Rec[8],'sind')     ,'sind','sine of the sub-solar latitude','none',('time'))   
            self.variables['gasp']=  self.Fort_var(self._ra_1D(Rec[9]*100,'gasp') ,'gasp','global average surface pressure','Pa',('time'))   
            
            Rec=self.f.read_record('i4')
            #NC3=Rec[0]; NCYCLE=Rec[1]
            
            self.variables['nc3']=     self.Fort_var(self._ra_1D(Rec[0],'nc3')     ,'nc3','full COMP3 is done every nc3 time steps.','None',('time'))   
            self.variables['ncycle']=  self.Fort_var(self._ra_1D(Rec[0],'ncycle')  ,'ncycle','ncycle','none',('time')) 
            
            self._log_var('ps','surface pressure','Pa',('time','lat','lon'),scaling=100)
            self._log_var('temp','temperature','K',('time','pfull','lat','lon'))
            self._log_var('ucomp','zonal wind','m/sec',('time','pfull','lat','lon'))
            self._log_var('vcomp','meridional wind','m/s',('time','pfull','lat','lon'))
            self._log_var('ts','surface temperature','K',('time','lat','lon'))
            self._log_var('snow','surface amount of CO2 ice on the ground','kg/m2',('time','lat','lon'))
            self._log_var('stressx','zonal component of surface stress','kg/m2',('time','lat','lon'))
            self._log_var('stressy','merdional component of surface stress','kg/m2',('time','lat','lon'))
            self._log_var('tstrat','stratosphere temperature','K',('time','lat','lon'))
            self._log_var('tausurf','visible dust optical depth at the surface.','none',('time','lat','lon'))
            self._log_var('ssun','solar energy absorbed by the atmosphere','W/m2',('time','lat','lon'))

            
            QTRACE=self.f.read_reals('f4').reshape(self.JM,self.IM,self.LM,self.ntrace,order='F') ;QTRACE[-1,:,:,:]=QTRACE[-2,:,:,:]
            QCOND=self.f.read_reals('f4').reshape(self.JM,self.IM,self.ntrace,order='F');QCOND[-1,:,:]=QCOND[-2,:,:]
            STEMP=self.f.read_reals('f4').reshape(self.JM,self.IM,self.NL,order='F');STEMP[-1,:,:]=STEMP[-2,:,:]
            
            #write(11) fuptopv, fdntopv, fupsurfv, fdnsurfv
            Rec=self.f.read_record('({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM))
            
            self._log_var('fuptopv','upward visible flux at the top of the atmosphere','W/m2',('time','lat','lon'),Rec=Rec[0])
            self._log_var('fdntopv','downward visible flux at the top of the atmosphere','W/m2',('time','lat','lon'),Rec=Rec[1])
            self._log_var('fupsurfv','upward visible flux at the surface','W/m2',('time','lat','lon'),Rec=Rec[2])
            self._log_var('fdnsurfv','downward visible flux at the surface','W/m2',('time','lat','lon'),Rec=Rec[3])
       
            #write(11) fuptopir, fupsurfir, fdnsurfir
            Rec=self.f.read_record('({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM),'({0},{1})f4'.format(self.IM,self.JM))

            self._log_var('fuptopir','upward IR flux at the top of the atmosphere','W/m2',('time','lat','lon'),Rec=Rec[0])
            self._log_var('fupsurfir','upward IR flux at the surface','W/m2',('time','lat','lon'),Rec=Rec[1])
            self._log_var('fdnsurfir','downward IR flux at the surface','W/m2',('time','lat','lon'),Rec=Rec[2])
            
            #write(11) surfalb
            self._log_var('surfalb','surface albedo in the visible, soil or H2O, CO2 ices if present','none',('time','lat','lon')) 

            #write(11) dheat
            #write(11) geot
            self._log_var('dheat','diabatic heating rate','K/sol',('time','pfull','lat','lon')) 
            self._log_var('geot','geopotential','m2/s2',('time','pfull','lat','lon'))  
 
    def _add_axis_as_variables(self):
        '''
        Add dimensions to the file as variables 
        '''
        self.variables['lat']=   self.Fort_var(self.lat,'lat','latitude','degree N',('lat'))
        self.variables['lon']=   self.Fort_var(self.lon,'lon','longitude','degrees_E',('lon')) 
        
    
        sgm =self.sigm
        pref=self.variables['psf'][0]# 7.01*100   in Pa
        pk=np.zeros(self.LM+1)
        bk=np.zeros(self.LM+1)
        pfull=np.zeros(self.LM)
        phalf=np.zeros(self.LM+1)
        
        pk[0]=0.08/2 #TODO [AK] changed pk[0]=.08 to pk[0]=.08/2, otherwise phalf[0] would be greater than phalf[1]
        for iz in range(self.LM):
            bk[iz+1] = sgm[2*iz+2]
        phalf[:]=pk[:]+pref*bk[:] # output in  Pa

        #First layer
        if pk[0]==0 and bk[0]==0: 
            pfull[0]=0.5*(phalf[0]+phalf[1])
        else:
            pfull[0]=(phalf[1]-phalf[0])/(np.log(phalf[1])-np.log(phalf[0]))
        #Rest of layers:   
        pfull[1:] = (phalf[2:]-phalf[1:-1])/(np.log(phalf[2:])-np.log(phalf[1:-1]))
        
        self.variables['phalf']= self.Fort_var(phalf,'phalf','ref half pressure level','Pa',('phalf'))
        self.variables['pfull']= self.Fort_var(pfull,'pfull','ref full pressure level','Pa',('pfull'))
        self.variables['bk']=    self.Fort_var(bk,'bk','vertical coordinate sigma value','none',('phalf'))
        self.variables['pk']=    self.Fort_var(pk,'pk','pressure part of the hybrid coordinate','Pa',('phalf'))
        self.variables['zgrid']= self.Fort_var(self.zgrid,'zgrid','depth at the mid-point of each soil layer','m',('zgrid')) 
        
       
    def _ls2sol_1year(self,Ls_deg,offset=True,round10=True):
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
        Ns=668.6     #Number of sols in 1 MY
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
        
    def _linInterpLs(self,Ls,stride=16):
        '''
        Linearly interpolate a step-wise 1D array
        Args:
            Ls     (float):   Input solar longitude
            stride   (int):   Default stride
        Returns
            Ls_out (float):   Ls
        ***NOTE***
        In the Legacy GCM fortran binaries, the solar longitude is only updated once per day, implying that 16 successive timesteps would have the same ls value.
        This routine linearly interpolate the ls between those successive values.   
        '''
        Ls=np.array(Ls);Ls_out=np.zeros_like(Ls)
        Lsdi=Ls[::stride]
        #Add a end point using the last Delta Ls:
        Lsdi=np.append(Lsdi,2*Lsdi[-1]-Lsdi[-2])
    
        for i in range(len(Ls)//stride):
            Ls_out[i*stride:(i+1)*stride]=np.arange(0,stride)/np.float(stride)*(Lsdi[i+1]-Lsdi[i])+Lsdi[i]
        return Ls_out
  