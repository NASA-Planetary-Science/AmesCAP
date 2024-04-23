#!/usr/bin/env python3
"""
The MarsFormat executable is a routine that transforms non-MGCM model
output into MGCM-like model output for compatibility with CAP.

MarsFormat changes variable names, dimension names, dimension order,
and units to the configuration expected by CAP. In some cases, such as
for MarsWRF, variables are derived and regridded onto a standard grid.

The executable requires:
    * ``[input_file]``              the file to be transformed

and optionally accepts:
    * ``[-openmars --openmars]``    convert openMars data to MGCM format
    * ``[-marswrf --marswrf]``      convert MarsWRF data to MGCM format

Third-party Requirements:
    * ``numpy``
    * ``os``
    * ``argparse``
    * ``xarray``
"""

# Make print statements appear in color
from amescap.Script_utils import (Cyan, Yellow, Nclr, Green)

# Load generic Python modules
import argparse     # Parse arguments
import os           # Access operating system functions
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# Load amesCAP modules
from amescap.Script_utils import prPurple,prCyan,prLightPurple,prRed,read_variable_dict_amescap_profile,prYellow,filter_vars,get_longname_units
from amescap.FV3_utils import daily_to_average, daily_to_diurn,layers_mid_point_to_boundary
from amescap.Ncdf_wrapper import Ncdf, Fort
xr.set_options(keep_attrs=True)

#---
# MarsFormat.py
# Routine to Transform Model Input (variable names, dimension names, array order)
# to expected configuration CAP

# ======================================================================
#                           ARGUMENT PARSER
# ======================================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}MarsFormat is for converting non-MGCM output "
        f"to MGCM format.{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument(
    "input_file", nargs="+",
    help=(
        f"A netCDF file or list of netCDF files.\n\n"
    )
)

parser.add_argument('-t', '--type', type=str, default='legacy',
                    help=""">  --type can be 'openmars', 'marswrf' or 'legacy' [DEFAULT is legacy] \n"""
                    """>  Usage: MarsFormat.py ****.nc \n"""
                    """          MarsFormat.py ****.nc -t openmars \n""")


# ===========================
path2data = os.getcwd()

def main():

   if not (parser.parse_args().type):
         prYellow(''' ***Notice***  No operation requested. Use '-type' and specify openmars, marswrf, legacy ''')
         exit()  # Exit cleanly

   path2data = os.getcwd()
   
   # Load all of the netcdf files
   file_list    = parser.parse_args().input_file
   model_type  = parser.parse_args().type  # e.g. 'legacy'
   for filei in file_list:
      #Add path unless full path is provided
      if not ('/' in filei):
         fullnameIN = path2data + '/' + filei
      else:
         fullnameIN=filei
      fullnameOUT = fullnameIN[:-3]+'_atmos_daily.nc'

      print('Processing...')
      #Load model variables,dimensions
      fNcdf=Dataset(fullnameIN,'r')
      model=read_variable_dict_amescap_profile(fNcdf)
      fNcdf.close()
      prCyan('Reading model attributes from ~.amescap_profile:')
      prCyan(vars(model)) #Print attribute
      #dataDIR = path+filename+'.nc'
      DS = xr.open_dataset(fullnameIN, decode_times=False)

      #=================================================================
      # ===================MarsWRF Specific Processing==================
      #=================================================================
      if model_type == 'marswrf':
         #TODO longname is 'description' for MarsWRF
         '''
         print('Input File content (description) and (description) attibutes:')
         print('------')
         for ivar in  DS.keys():
            print(ivar,DS[ivar].attrs['description'],DS[ivar].attrs['units'])
         print('------')
         '''
         #==================================================================
         # First save all variable descriptions in attrs longname
         #================================================================== 
         for var_name in DS.data_vars:
            var = DS[var_name]
            if 'description' in var.attrs:
               var.attrs['long_name'] = var.attrs['description'] 
         
         #==================================================================
         # Reformat Dimension Variables/Coords as Needed
         #================================================================== 
         time        = DS[model.time]/ 60/ 24         # minutes since simulation start [m]
         lat = DS[model.lat][0,:,0]
         lon2D = DS[model.lon][0,:]
         lon = np.squeeze(lon2D[0,:])
         DS[model.lon],DS[model.lat],DS[model.time]=lon,lat,time
         
         # Derive half and full reference pressure levels (Pa)
         pfull = DS.P_TOP[0]+ DS.ZNU[0,:]* DS.P0 
         phalf= DS.P_TOP[0]+ DS.ZNW[0,:]* DS.P0
         DS = DS.assign_coords(phalf=phalf,pfull=pfull)
         DS.phalf.attrs['long_name'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
         DS.phalf.attrs['description'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
         DS.phalf.attrs['units'] = 'Pa'
 
         # Update dimensions
         DS = DS.assign_coords(dimensions='phalf')

         # Update Variable Description, Longname and Unit
         DS[model.lon].attrs['long_name'] = '(MODIFIED IN POST PROCESSING) ' + DS[model.lon].attrs['description']
         DS[model.lat].attrs['long_name'] = '(MODIFIED IN POST PROCESSING) ' + DS[model.lat].attrs['description']
        
         DS[model.time].attrs['long_name'] = '(MODIFIED IN POST PROCESSING) days since simulation start, time/60/24'
         DS[model.lon].attrs['description'] = '(MODIFIED IN POST PROCESSING) ' + DS[model.lon].attrs['description']
         DS[model.lat].attrs['description'] = '(MODIFIED IN POST PROCESSING) ' + DS[model.lat].attrs['description']

         DS[model.time].attrs['description'] = '(MODIFIED IN POST PROCESSING) days since simulation start, time/60/24'
         DS[model.time].attrs['units'] = '(MODIFIED IN POST PROCESSING) days '
        

         #==================================================================
         # Derive ak, bk
         #==================================================================
         ak  = np.zeros(len(DS.phalf))
         bk = np.zeros(len(DS.phalf))
         ak[-1]=DS.P_TOP[0]  #MarsWRF comes with pressure increasing with N
         bk[:]=DS.ZNW[0,:]

         DS['ak'],DS['bk'] = ak,bk

         # Update Variable Description & Longname
         DS['ak'].attrs['long_name'], DS['bk'].attrs['long_name'] = '(ADDED IN POST PROCESSING)', '(ADDED IN POST PROCESSING)'
         DS['ak'].attrs['description'], DS['bk'].attrs['description'] = '(ADDED IN POST PROCESSING)', '(ADDED IN POST PROCESSING)'
        
         #==================================================================
         # Calculate *Level* Heights above the Surface (i.e. above topo)
         #==================================================================
         zagl_lvl = (DS.PH[:,:-1,:,:] + DS.PHB[0,:-1,:,:]) / DS.G - DS.HGT[0,:,:]

         #==================================================================
         # Find Layer Pressures [Pa]
         #=================================================================
         try:
            pfull3D = DS.P_TOP + DS.PB[0,:] # perturb. pressure + base state pressure, time-invariant
         except NameError:
            pfull3D = DS[model.ps][:,:-1,:-1] * DS.ZNU[:,:-1]

         #==================================================================
         # Derive attmospherice temperature [K]
         #==================================================================
         gamma = DS.CP / (DS.CP - DS.R_D)
         temp = (DS.T + DS.T0) * (pfull3D / DS.P0)**((gamma-1.) / gamma)
         DS = DS.assign(temp=temp)
         DS['temp'].attrs['description'] = '(ADDED IN POST PROCESSING) Temperature'
         DS['temp'].attrs['long_name'] = '(ADDED IN POST PROCESSING) Temperature'
         DS['temp'].attrs['units'] = 'K'

         #==================================================================
         # Interpolate U, V, W, Zfull onto Regular Mass Grid (from staggered)
         #==================================================================
         # For variables staggered x (lon) [t,z,y,x'] -> regular mass grid [t,z,y,x]:
         # Step 1: Identify variables with the dimension 'west_east_stag' and _U not in the variable name (these are staggered grid identifiers like the staggered latitude, I think
         variables_with_west_east_stag = [var for var in DS.variables if 'west_east_stag' in DS[var].dims and '_U' not in var]

         prCyan('Interpolating Staggered Variables to Standard Grid')
         # Loop through, and unstagger. the dims_list process finds the dimensoins of the variable and replaces west_east_stag with west_east
         prLightPurple('     From west_east_stag to west_east: ' + ', '.join(variables_with_west_east_stag))
         for var_name in variables_with_west_east_stag:
            var = getattr(DS, var_name)
            dims = var.dims
            # replace west_east_stag in dimensions with west_east
            dims_list = list(dims)
            for i, dim in enumerate(dims_list):
               if dim == 'west_east_stag':
                  dims_list[i]='west_east'
                  break # Stop the loop once the replacement is made
            new_dims = tuple(dims_list)
            
            transformed_var = 0.5 * (var.sel(west_east_stag=slice(None, -1)) + var.sel(west_east_stag=slice(1, None)))
            DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, coords={'XLAT':DS['XLAT']})
           
            DS[var_name].attrs['description'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
            DS[var_name].attrs['long_name'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
            DS[var_name].attrs['stagger'] = 'USTAGGERED IN POST-PROCESSING'
           
        
         # For variables staggered y (lat) [t,z,y',x] -> regular mass grid [t,z,y,x]: 
         variables_with_south_north_stag = [var for var in DS.variables if 'south_north_stag' in DS[var].dims and '_V' not in var]
         prLightPurple('     From south_north_stag to south_north: '  + ', '.join(variables_with_south_north_stag))
         for var_name in variables_with_south_north_stag:
            var = getattr(DS, var_name)
            dims = var.dims
            # replace west_east_stag in dimensions with west_east
            dims_list = list(dims)
            for i, dim in enumerate(dims_list):
               if dim == 'south_north_stag':
                  dims_list[i]='south_north'
                  break # Stop the loop once the replacement is made
            new_dims = tuple(dims_list)
            
            transformed_var = 0.5 * (var.sel(south_north_stag=slice(None, -1)) + var.sel(south_north_stag=slice(1, None)))
            DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, coords={'XLONG':DS['XLONG']})
           
            DS[var_name].attrs['description'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']

            DS[var_name].attrs['long_name'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
            DS[var_name].attrs['stagger'] = 'USTAGGERED IN POST-PROCESSING'
         
         #DS[model.vcomp] = 0.5 * (DS[model.vcomp][:,:,:-1,:] + DS[model.vcomp][:,:,1:,:])

         # For variables staggered p/z (height) [t,z',y,x] -> regular mass grid [t,z,y,x]:
         variables_with_bottom_top_stag = [var for var in DS.variables if 'bottom_top_stag' in DS[var].dims and 'ZNW' not in var and 'phalf' not in var]
         prLightPurple('     From bottom_top_stag to bottom_top: '  + ', '.join(variables_with_bottom_top_stag))
         for var_name in variables_with_bottom_top_stag:
            var = getattr(DS, var_name)
            dims = var.dims
            # replace bottom_top_stag
            dims_list = list(dims)
            for i, dim in enumerate(dims_list):
               if dim == 'bottom_top_stag':
                  dims_list[i]='bottom_top'
                  break # Stop the loop once the replacement is made
            new_dims = tuple(dims_list)
            transformed_var = 0.5 * (var.sel(bottom_top_stag=slice(None, -1)) + var.sel(bottom_top_stag=slice(1, None)))
            DS[var_name] = xr.DataArray(transformed_var, dims=new_dims)
            DS[var_name].attrs['description'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
            DS[var_name].attrs['long_name'] = '(UNSTAGGERED IN POST-PROCESSING) ' + DS[var_name].attrs['description']
            DS[var_name].attrs['stagger'] = 'USTAGGERED IN POST-PROCESSING'
         

         #DS[model.w] = 0.5 * (DS[model.w][:,:-1,:,:] + DS[model.w][:,1:,:,:])

         # ALSO INTERPOLATE TO FIND *LAYER* HEIGHTS ABOVE THE SURFACE (i.e., above topography; m)
         zfull3D = 0.5 * (zagl_lvl[:,:-1,:,:] + zagl_lvl[:,1:,:,:])

      #=================================================================
      # ===================OpenMars Specific Processing==================
      #=================================================================
      elif model_type == 'openmars':
         '''
         print('Input File content (FIELDNAM) and (UNITS) attibutes:')
         print('------')
         for ivar in  DS.keys():
            print(ivar,DS[ivar].attrs['FIELDNAM'],DS[ivar].attrs['UNITS'])
         print('------')
         '''
         #==================================================================
         # First save all variable FIELDNAM in attrs longname
         #==================================================================
         for var_name in DS.data_vars:
            var = DS[var_name]
            if 'FIELDNAM' in var.attrs:
               var.attrs['long_name'] = var.attrs['FIELDNAM']

         #==================================================================
         # Define Coordinates for New DataFrame
         #==================================================================
         ref_press=720 #TODO this is added on to create ak/bk
         time        = DS[model.dim_time]         # minutes since simulation start [m]
         lat = DS[model.dim_lat]  #Replace DS.lat
         lon = DS[model.dim_lon]

         # Derive half and full reference pressure levels (Pa)
         DS = DS.assign(pfull=DS[model.dim_pfull]*ref_press)
         DS['pfull'].attrs['FIELDNAM']='(MODIFIED IN POST-PROCESSING) ' + DS['pfull'].attrs['FIELDNAM']
         DS['pfull'].attrs['long_name']='(MODIFIED IN POST-PROCESSING) ' + DS['pfull'].attrs['FIELDNAM']

         #==================================================================
         # add ak,bk as variables
         # add p_half dimensions as vertical grid coordinate 
         #==================================================================

         #Compute sigma values. Swap the sigma array upside down twice  with [::-1] since the layers_mid_point_to_boundary() needs sigma[0]=0, sigma[-1]=1) and then to reorganize the array in the original openMars format with sigma[0]=1, sigma[-1]=0
         bk = layers_mid_point_to_boundary(DS[model.dim_pfull][::-1],1.)[::-1]
         ak = np.zeros(len(DS[model.dim_pfull]) + 1)

         DS['phalf']= ak + ref_press*bk
         DS.phalf.attrs['long_name'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
         DS.phalf.attrs['description'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
         DS.phalf.attrs['units'] = 'Pa'
          
         DS = DS.assign(bk=(model.dim_phalf, np.array(bk)))
         DS = DS.assign(ak=(model.dim_phalf, np.zeros(len(DS[model.dim_pfull]) + 1)))
 
         # Update Variable Description & Longname
         DS['ak'].attrs['long_name'], DS['bk'].attrs['long_name'] = '(ADDED IN POST PROCESSING)', '(ADDED IN POST PROCESSING)'
         DS['ak'].attrs['FIELDNAM'], DS['bk'].attrs['FIELDNAM'] = '(ADDED IN POST PROCESSING)', '(ADDED IN POST PROCESSING)'
        


      #==================================================================
      # START PROCESSING FOR ALL MODELS
      #==================================================================

      #==================================================================
      # check that vertical grid starts at toa with highest level at surface
      #==================================================================

      if DS[model.dim_pfull][0] != DS[model.dim_pfull].min(): # if toa, lev = 0 is surface then flip
          DS = DS.isel(**{model.dim_pfull: slice(None, None, -1)})
          DS=DS.isel(**{model.dim_phalf: slice(None, None, -1)}) #Also flip phalf,ak, bk
          prRed('NOTE: all variables flipped along vertical dimension, so that the top of the atmosphere is now index 0')

      #==================================================================
      # reorder dimensions
      #==================================================================
      prCyan('Transposing variable dimensions to match order expected in CAP') 
      DS = DS.transpose(model.dim_time, model.dim_pfull ,model.dim_lat,model.dim_lon, ...)

      #==================================================================
      # change longitude from -180-179 to 0-360
      #==================================================================
      if min(DS[model.dim_lon]) < 0:      
            tmp = np.array(DS[model.dim_lon])
            tmp = np.where(tmp<0,tmp+360,tmp)
            DS[model.dim_lon] = tmp
#DS=DS.assign_coords({model.dim_lon:(model.dim_lon,tmp,DS[model.dim_lon].attrs)})
            DS = DS.sortby(model.dim_lon)
            prRed('NOTE: Longitude changed to 0-360E and all variables appropriately reindexed')

      #==================================================================
      # add scalar axis to areo [time, scalar_axis])
      #==================================================================
      inpt_dimlist = DS.dims
      # first check if dimensions are correct and don't need to be modified
      if 'scalar_axis' not in inpt_dimlist:           # first see if scalar axis is a dimension
            scalar_axis = DS.assign_coords(scalar_axis=1)
      if DS[model.areo].dims != (model.time,scalar_axis):
            DS[model.areo] = DS[model.areo].expand_dims('scalar_axis', axis=1)
            DS[model.areo].attrs['long_name'] = '(SCALAR AXIS ADDED IN POST-PROCESSING) ' + DS[model.areo].attrs['long_name']
           
            prRed('NOTE: scalar axis added to aerocentric longitude')


      # Output Processed Data to New **atmos_daily.nc File
      #==================================================================
      DS.to_netcdf(fullnameOUT)
      prCyan(fullnameOUT +' was created')


      # Create **atmos_average.nc file
      #==================================================================
      fullnameOUT = fullnameIN[:-3]+'_atmos_average'+'.nc'

      # Figure out number of timesteps per 5 sol
      dt_in = DS[model.time][1]-DS[model.time][0]
      iperday = int(np.round(1/dt_in))
      combinedN = int(iperday*5)
      time = model.dim_time
      
      # Coarsen the 'time' dimension by a factor of 5 and average over each window
      DS_average = DS.coarsen(**{model.dim_time:combinedN}).mean()

      # Update the time dimension's coordinate values to reflect the new time axis
      # save the middle value of the time interval
      # start indexing from 2 (3 steps in or center of 5 sol segment) and jumping in
      # groups of 5

      # Update the time coordinate attribute
      DS_average[model.time].attrs['long_name'] = 'time averaged over 5 sols'

      # Creat New File
      DS_average.to_netcdf(fullnameOUT)
      prCyan(fullnameOUT +' was created')

      #==================================================================
      # Create **atmos_diurn.nc file
      #==================================================================
      fullnameOUT = fullnameIN[:-3]+'_atmos_diurn'+'.nc'

      # Figure out number of timesteps per 5 sol
      dt_in = DS[model.time][1]-DS[model.time][0]
      iperday = int(np.round(1/dt_in))
      days=len(DS[model.time])/iperday

      # create a new time of day dimension
      tod_name = 'time_of_day_%02d' %(iperday)
      prLightPurple('tod_name'+ tod_name)
      DS_diurn = DS.copy()

      # get the time of day in hours
      tod = np.mod(DS[model.time][0:iperday]*24,24).values
      
      # specify labels for new dimensions
      
      ind = pd.MultIndex.from_product((days,tod),names=('time',tod_name))
      DS_reshaped = DS_dirun.assign(time=ind).unstack('time')
      
      print(DS_reshaped)
      stop
      #reshape methood
      DS_reshaped = xr.Dataset({var_name: var.isel(time=slice(0,None,iperday)).stack(time=[days,tod]) for var_name, var in DS_diurn.data_vars.items() if 'time' in var.dims})


      print(DS_reshaped['tsurf'][0,:,0,0].values)
      stop   
 
      # first make a coarser version (1 time step per sol)
      DS_coarse = DS_diurn.coarsen(time=iperday).mean()
      
      print(DS_diurn['tsurf'][0:11,0,0].values)
      print(DS_coarse['tsurf'][0,0,0].values)
     
      # now add a new time of day dimension if time is a dimensiono in the variable
      DS_coarse.update({var_name: var.expand_dims(**{tod_name:iperday},axis=1) for var_name, var in DS_coarse.data_vars.items() if 'time' in var.dims})
     
      print(DS_coarse['tsurf'][0,:,0,0].values)

      for i in range(iperday):
          print(i, DS_diurn['tsurf'][i,0,0].values)
          #DS_coarse['tsurf'].isel(**{tod_name:i},time=0).values = DS_diurn['tsurf'].isel(time=i).values
          DS_coarse['tsurf'][0,i,:,:].values = DS_diurn['tsurf'][i,:,:].values
          #.isel(**{tod_name:i},time=0).values = DS_diurn['tsurf'].isel(time=i).values
     
      print(DS_coarse['tsurf'][0,:,0,0].values)

      stop
      # add a new time of day dimension if time is a dimension in the variable
      DS_diurn.update({var_name: var.expand_dims(**{tod_name:iperday},axis=1) for var_name, var in DS_diurn.data_vars.items() if 'time' in var.dims})

      # get the time of day in hours
      tod = np.mod(DS[model.time][0:iperday]*24,24).values
      print('tod=', tod)

      # Sort the time of day, e.g. if tod = [6,12,18,0] re-arrange into [0,6,12,18]
      # every element in array must be greater than the one too its left
      if not np.all(tod[1:] >= tod[:-1]):

          # this returns the permutation, e..g if tod= [6,12,18,0] i_sort = [3,0,1,2]
          i_sort = np.argsort(tod)
          print(i_sort)
          # reorder tod
          tod = tod[i_sort]
          DS_diurn.update({var_name:var.isel(**{tod_name:i_sort}) for var_name, var in DS_diurn.data_vars.items() if 'time' in var.dims})

     
      print(DS_diurn['tsurf'][0:11,0,0,0])

      # now squeeze time(i:iperday)==i and time_of-day grabs the rest
      DS_coarse = DS_diurn.coarsen(time=iperday).mean()
      
      print(DS_coarse['tsurf'][0,:,0,0])
      stop


       # get the time of day in hours
      tod = np.mod(DS[model.time][0:iperday]*24,24).values
      DS_coarse.coords[tod_name] = (tod_name,tod)

      # add a new time of day dimension if time is a dimension in the variable
      DS_coarse.update({var_name: var.expand_dims(**{tod_name:iperday},axis=1) for var_name, var in DS_coarse.data_vars.items() if 'time' in var.dims})

      print(DS_coarse['tsurf'][0:11,0,0,0])
      for i in range(iperday):
         DS_coarse['tsurf'].isel(**{tod_name:i},time=0).values=DS['tsurf'].isel(time=i).values

      print(DS.tsurf[0:11,0,0])
      print(DS_coarse.tsurf[0,:,0,0])


      # add a new time of day dimension if time is a dimension in the variable
      DS_diurn.update({var_name: var.expand_dims(**{tod_name:iperday},axis=1) for var_name, var in DS_diurn.data_vars.items() if 'time' in var.dims})

      print(DS_diurn.tsurf[-1,:,24,48])
      # get the time of day in hours
      tod = np.mod(DS[model.time][0:iperday]*24,24).values
      print('tod=', tod)

      # Sort the time of day, e.g. if tod = [6,12,18,0] re-arrange into [0,6,12,18]
      # every element in array must be greater than the one too its left
      if not np.all(tod[1:] >= tod[:-1]):

          # this returns the permutation, e..g if tod= [6,12,18,0] i_sort = [3,0,1,2]
          i_sort = np.argsort(tod)
          print(i_sort)
          # reorder tod
          tod = tod[i_sort]
          DS_diurn.update({var_name:var.isel(**{tod_name:i_sort}) for var_name, var in DS_diurn.data_vars.items() if 'time' in var.dims})

      # now coarsen to 5 sol resolution
      combinedN = int(iperday*5)   # number of 5 sol groups in file
      time = model.dim_time
      
      # Coarsen the 'time' dimension by a factor of 5 and average over each window
      DS_diurnave = DS_diurn.coarsen(**{model.dim_time:combinedN}).mean()

      # Update the time dimension's coordinate values to reflect the new time axis
      # save the middle value of the time interval
      # start indexing from 2 (3 steps in or center of 5 sol segment) and jumping in
      # groups of 5

      # Update the time coordinate attribute
      DS_diurnave[model.time].attrs['long_name'] = 'time averaged over 5 sols'

      print(DS_diurnave)
      # Creat New File
      DS_diurnave.to_netcdf(fullnameOUT)
      prCyan(fullnameOUT +' was created')
 
     


# Creat New File
      #DS_diurn.to_netcdf(fullnameOUT)

# ===========================================================================
# ===============  Bin a 'daily' file to an 'average' file ==================
# ===========================================================================
def bin_average(fullnameIN,model):

    nday = 5

    fullnameOUT = fullnameIN[:-3]+'_to_average'+'.nc'

    fdaily = Dataset(fullnameIN, 'r', format='NETCDF4_CLASSIC')
    var_list = filter_vars(fdaily) # parser.parse_args().include)  # Get all variables
   
    #print(fdaily.dimensions) 
    time_in = fdaily.variables[model.time][:]
    Nin = len(time_in)
    dt_in = time_in[1]-time_in[0]
    iperday = int(np.round(1/dt_in))
    combinedN = int(iperday*nday)

    N_even = Nin//combinedN
    N_left = Nin % combinedN

    if N_left != 0:
       prYellow('***Warning*** requested  %i sols bin period. File has %i timestep/sols and %i/(%i x %i) is not a round number' %
          (nday, iperday, Nin, nday, iperday))
       prYellow('    Will use %i  bins of (%i x %i)=%i timesteps (%i) and discard %i timesteps' % (
          N_even, nday, iperday, combinedN, N_even*combinedN, N_left))

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

       #print(ivar, varNcf.dimensions)
       if model.time in varNcf.dimensions:
          prCyan("Processing: %s ..." % (ivar))
          var_out = daily_to_average(varNcf[:], dt_in, nday)
          longname_txt, units_txt = get_longname_units(fdaily, ivar)
          fnew.log_variable(
             ivar, var_out, varNcf.dimensions, longname_txt, units_txt)
       else:
          if ivar in ['pfull', 'lat', 'lon', 'phalf', 'pk', 'bk', 'pstd', 'zstd', 'zagl']:
             prCyan("Copying axis: %s..." % (ivar))
             fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
          else:
             prCyan("Copying variable: %s..." % (ivar))
             fnew.copy_Ncvar(fdaily.variables[ivar])
       fnew.close()

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == '__main__':
    main()
