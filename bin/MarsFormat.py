#!/usr/bin/env python3
"""
The MarsFormat executable is for converting non-MGCM data, such as that
from EMARS, OpenMARS, PCM, and MarsWRF, into MGCM-like netCDF data 
products. The MGCM is the NASA Ames Mars Global Climate Model developed
and maintained by the Mars Climate Modeling Center (MCMC). The MGCM 
data repository is available at data.nas.nasa.gov/mcmc.

The executable requires two arguments:

    * ``[input_file]``         The file to be transformed
    * ``[-gcm --gcm_name]``    The GCM from which the file originates

and optionally accepts:

    * ``[-rn --retain_names]`` Preserve original variable and dimension names
    * ``[-ba, --bin_average]`` Bin non-MGCM files like 'average' files
    * ``[-bd, --bin_diurn]``   Bin non-MGCM files like 'diurn' files

Third-party requirements:

    * ``numpy``
    * ``netCDF4``
    * ``sys``
    * ``argparse``
    * ``os``

List of Functions:

    * download - Queries the requested file from the NAS Data Portal.

"""

# Make print statements appear in color
from amescap.Script_utils import (
    Green, Yellow, Red, Blue, Cyan, Purple, Nclr,
)

# Load generic Python modules
import sys          # System commands
import argparse     # Parse arguments
import os           # Access operating system functions
import re           # Regular expressions
import numpy as np
import xarray as xr
from netCDF4 import Dataset

# Load amesCAP modules
from amescap.Script_utils import (
    read_variable_dict_amescap_profile, reset_FV3_names
)
from amescap.FV3_utils import layers_mid_point_to_boundary

xr.set_options(keep_attrs=True)

import functools
import traceback

def debug_wrapper(func):
    """
    A decorator that wraps a function with error handling based on the 
    --debug flag.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        global debug
        try:
            return func(*args, **kwargs)
        except Exception as e:
            if debug:
                # In debug mode, show the full traceback
                print(f"{Red}ERROR in {func.__name__}: {str(e)}{Nclr}")
                traceback.print_exc()
            else:
                # In normal mode, show a clean error message
                print(f"{Red}ERROR in {func.__name__}: {str(e)}\nUse "
                      f"--debug for more information.{Nclr}")
            return 1  # Error exit code
    return wrapper

# ======================================================
#                  ARGUMENT PARSER
# ======================================================
parser = argparse.ArgumentParser(
    prog=('MarsFormat'),
    description=(
        f"{Yellow} Converts model output to MGCM-like format for "
        f"compatibility with CAP."
        f"{Nclr}\n\n"
    ),
    formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('input_file', nargs='+',
    type=argparse.FileType('rb'),
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument('-gcm', '--gcm_name', type=str,
    choices=['marswrf', 'openmars', 'pcm', 'emars'],
    help=(
        f"Acceptable types include 'openmars', 'marswrf', 'emars', "
        f"and 'pcm' \n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars\n"
        f"{Blue}(Creates openmars_file_daily.nc)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-ba', '--bin_average', nargs="?", const=5,
    type=int,
    help=(
        f"Calculate 5-day averages from instantaneous data. Generates "
        f"MGCM-like 'average' files.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -ba\n"
        f"{Blue}(Creates openmars_file_average.nc; 5-sol bin){Green}\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -ba 10\n"
        f"{Blue}(10-sol bin)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('-bd', '--bin_diurn', action='store_true',
    default=False,
    help=(
        f"Calculate 5-day averages binned by hour from instantaneous "
        f"data. Generates MGCM-like 'diurn' files.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -bd\n"
        f"{Blue}(Creates openmars_file_daily.nc;  5-sol bin){Green}\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -bd -ba 10\n"
        f"{Blue}(10-sol bin)"
        f"{Nclr}\n\n"
    )
)

# Secondary arguments: Used with some of the arguments above

parser.add_argument('-rn', '--retain_names', action='store_true',
    default=False,
    help=(
        f"Preserves the names of the variables and dimensions in the"
        f"original file.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars -rn\n"
        f"{Blue}(Creates openmars_file_daily.nc)"
        f"{Nclr}\n\n"
    )
)

parser.add_argument('--debug', action='store_true',
    help=(
        f"Use with any other argument to pass all Python errors and\n"
        f"status messages to the screen when running CAP.\n"
        f"{Green}Example:\n"
        f"> MarsFormat openmars_file.nc -gcm openmars --debug"
        f"{Nclr}\n\n"
    )
 )

args = parser.parse_args()
debug = args.debug

if args.input_file:
    for file in args.input_file:
        if not re.search(".nc", file.name):
            parser.error(f"{Red}{file.name} is not a netCDF file{Nclr}")
            exit()

# ===========================
path2data = os.getcwd()
ref_press = 725  # TODO hard-coded reference pressure

@debug_wrapper
def main():
    ext = ''  # Initialize empty extension
    if args.gcm_name not in ['marswrf', 'openmars', 'pcm', 'emars']:
        print(f"{Yellow}***Notice***  No operation requested. Use '-gcm' and specify openmars, marswrf, pcm, emars")
        exit()  # Exit cleanly

    path2data = os.getcwd()

    # Load all of the netcdf files
    file_list = [f.name for f in args.input_file]
    model_type = args.gcm_name  # e.g. 'legacy'
    for filei in file_list:
        # Add path unless full path is provided
        if '/' not in filei:
            fullnameIN = path2data + '/' + filei
        else:
            fullnameIN = filei

        print('Processing...')
        # Load model variables, dimensions
        fNcdf = Dataset(fullnameIN, 'r')
        model = read_variable_dict_amescap_profile(fNcdf)
        fNcdf.close()

        print(f"{Cyan}Reading model attributes from ~.amescap_profile:")
        print(f"{Cyan}{vars(model)}") # Print attribute

        # Open dataset with xarray
        DS = xr.open_dataset(fullnameIN, decode_times=False)

        #=================================================================
        # ===================MarsWRF Specific Processing==================
        #=================================================================
        if model_type == 'marswrf':
            # First save all variable descriptions in attrs longname
            print(f"{Cyan}Current variables at top of marswrf processing:\n{list(DS.variables)}{Nclr}\n")
            for var_name in DS.data_vars:
                var = DS[var_name]
                if 'description' in var.attrs:
                    var.attrs['long_name'] = var.attrs['description']

            # Ensure mandatory dimensions and coordinates exist
            try:
                # Reformat Dimension Variables/Coords as Needed
                # Handle potential missing dimensions or coordinates
                if model.time not in DS:
                    raise KeyError(f"Time dimension {model.time} not found")
                if model.lat not in DS:
                    raise KeyError(f"Latitude dimension {model.lat} not found")
                if model.lon not in DS:
                    raise KeyError(f"Longitude dimension {model.lon} not found")

                # Time conversion (minutes to days)
                time = DS[model.time] / 60 / 24 if 'time' in DS else None

                # Handle latitude and longitude 
                lat = DS[model.lat][0,:,0] if len(DS[model.lat].shape) > 1 else DS[model.lat]
                lon = DS[model.lon][0,0,:] if len(DS[model.lon].shape) > 1 else DS[model.lon]
                
                # Convert longitudes to 0-360 range
                lon360 = (lon + 360) % 360

                # Update coordinates 
                if time is not None:
                    DS[model.time] = time
                DS[model.lon] = lon360
                DS[model.lat] = lat

                # Derive pressure levels 
                if hasattr(DS, 'P_TOP') and hasattr(DS, 'ZNU') and hasattr(DS, 'P0'):
                    pfull = DS.P_TOP[0] + DS.ZNU[0,:] * DS.P0
                    phalf = DS.P_TOP[0] + DS.ZNW[0,:] * DS.P0
                    
                    DS = DS.assign_coords(phalf=phalf, pfull=pfull)
                    DS.phalf.attrs['long_name'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
                    DS.phalf.attrs['description'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
                    DS.phalf.attrs['units'] = 'Pa'

                    # Derive ak, bk
                    ak = np.zeros(len(DS.phalf))
                    bk = np.zeros(len(DS.phalf))
                    ak[-1] = DS.P_TOP[0]  # MarsWRF comes with pressure increasing with N
                    bk[:] = DS.ZNW[0,:]

                    DS['ak'], DS['bk'] = ak, bk
                    DS['ak'].attrs['long_name'] = '(ADDED IN POST PROCESSING)'
                    DS['bk'].attrs['long_name'] = '(ADDED IN POST PROCESSING)'

                # Calculate layer heights (if PH and PHB exist)
                if 'PH' in DS and 'PHB' in DS and 'HGT' in DS and 'G' in DS:
                    zagl_lvl = (DS.PH[:,:-1,:,:] + DS.PHB[0,:-1,:,:]) / DS.G - DS.HGT[0,:,:]
                    zfull3D = 0.5 * (zagl_lvl[:,:-1,:,:] + zagl_lvl[:,1:,:,:])

                # Derive atmospheric temperature
                if all(key in DS for key in ['T', 'T0', 'CP', 'R_D', 'P0', 'PB']):
                    gamma = DS.CP / (DS.CP - DS.R_D)
                    pfull3D = DS.P_TOP + DS.PB[0,:]
                    DS['T'] = (DS.T + DS.T0) * (pfull3D / DS.P0)**((gamma-1.) / gamma)
                    DS['T'].attrs['description'] = '(MODIFIED IN POST PROCESSING) Temperature'
                    DS['T'].attrs['long_name'] = '(MODIFIED IN POST PROCESSING) Temperature'
                    DS['T'].attrs['units'] = 'K'

                # Drop 'Times' variable if present (non-numerical values)
                if 'Times' in DS:
                    DS = DS.drop_vars("Times")

            except Exception as e:
                print(f"{Red}Error processing MarsWRF data: {e}{Nclr}")
                print(f"{Red}Available dimensions: {list(DS.dims)}{Nclr}")
                print(f"{Red}Available variables: {list(DS.variables)}{Nclr}")
                raise
            
        # ==============================================================
        #                  OpenMars Specific Processing
        # ==============================================================
        elif model_type == 'openmars':
            # First save all variable FIELDNAM in attrs longname
            for var_name in DS.data_vars:
                #TODO grab vars and dims here!
                var = DS[var_name]
                if 'FIELDNAM' in var.attrs:
                    var.attrs['long_name'] = var.attrs['FIELDNAM']
                if 'UNITS' in var.attrs:
                    var.attrs['units'] = var.attrs['UNITS']


            # Define Coordinates for New DataFrame
            time        = DS[model.dim_time]         # minutes since simulation start [m]
            lat = DS[model.dim_lat]  #Replace DS.lat
            lon = DS[model.dim_lon]

            DS = DS.assign(pfull=DS[model.dim_pfull]*ref_press)

            DS['pfull'].attrs['long_name']='(MODIFIED IN POST-PROCESSING) ' + DS['lev'].attrs['FIELDNAM']
            DS['pfull'].attrs['long_name']='Pa'

            # add ak,bk as variables
            # add p_half dimensions as vertical grid coordinate

            #Compute sigma values. Swap the sigma array upside down twice  with [::-1] since the layers_mid_point_to_boundary() needs sigma[0]=0, sigma[-1]=1) and then to reorganize the array in the original openMars format with sigma[0]=1, sigma[-1]=0
            bk = layers_mid_point_to_boundary(DS[model.dim_pfull][::-1],1.)[::-1]
            ak = np.zeros(len(DS[model.dim_pfull]) + 1)

            DS[model.phalf]= ak + ref_press*bk
            DS.phalf.attrs['long_name'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
            DS.phalf.attrs['description'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
            DS.phalf.attrs['units'] = 'Pa'

            DS = DS.assign(bk=(model.dim_phalf, np.array(bk)))
            DS = DS.assign(ak=(model.dim_phalf, np.zeros(len(DS[model.dim_pfull]) + 1)))
            
            print('dim phalf=',model.dim_phalf)
            # Update Variable Description & Longname
            DS['ak'].attrs['long_name']='(ADDED IN POST PROCESSING) pressure part of the hybrid coordinate'
            DS['ak'].attrs['units']='Pa'
            DS['bk'].attrs['long_name'] = '(ADDED IN POST PROCESSING) vertical coordinate sigma value'
            DS['bk'].attrs['units']='None'

        # ==============================================================
        #                 Emars Specific Processing
        # ==============================================================
        elif model_type == 'emars':
            # Interpolate U, V, onto Regular Mass Grid (from staggered)
            print(f"{Cyan}Interpolating Staggered Variables to Standard Grid")
            variables_with_latu = [var for var in DS.variables if 'latu' in DS[var].dims]
            variables_with_lonv = [var for var in DS.variables if 'lonv' in DS[var].dims]
            # Loop through, and unstagger. the dims_list process finds the dimensions of the variable and replaces west_east_stag with west_east
            print('     From latu to lat: ' + ', '.join(variables_with_latu ))
            for var_name in variables_with_latu:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt=var.long_name
                units_txt=var.units

                # replace latu in dimensions with lat
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'latu':
                        dims_list[i]='lat'
                        break # Stop the loop once the replacement is made
                new_dims = tuple(dims_list)

                #TODO  Using  'values' as var.isel(latu=slice(None, -1)).values and var.isel(latu=slice(None, -1)) returns array of different size
                #TODO the following reproduce the 'lat' array, but not if the keyword values is ommited
                #latu1_val=ds.latu.isel(latu=slice(None, -1)).values
                #latu2_val=ds.latu.isel(latu=slice(1,  None)).values

                #newlat_val= 0.5 * (latu1 + latu2)
                #newlat_val=np.append(newlat_val,0)
                #newlat_val ==lat
                transformed_var = 0.5 * (var.isel(latu=slice(None, -1)).values + var.isel(latu=slice(1, None)).values)

                # This is equal to lat[0:-1]
                #Add padding at the pole to conserve the same dimension as lat.
                list_pads=[]
                for i, dim in enumerate(new_dims):
                    if dim == 'lat':
                        list_pads.append((0,1)) #(begining,end)=(0,1), meaning 1 padding at the end of the array
                    else:
                        list_pads.append((0,0)) #(begining,end)=(0,0), no pad on that axis
                transformed_var=np.pad(transformed_var,list_pads,mode='constant',constant_values=0)
                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, coords={'lat':DS['lat']})
                DS[var_name].attrs['long_name'] = '(UNSTAGGERED IN POST-PROCESSING) ' + longname_txt
                DS[var_name].attrs['units'] = units_txt
            #==============
            print('     From lonv to lon: ' + ', '.join(variables_with_lonv ))
            for var_name in variables_with_lonv:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt=var.long_name
                units_txt=var.units

                # replace lonv in dimensions with lon
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'lonv':
                        dims_list[i]='lon'
                        break # Stop the loop once the replacement is made
                new_dims = tuple(dims_list)

                transformed_var = 0.5 * (var.isel(lonv=slice(None, -1)).values + var.isel(lonv=slice(1, None)).values)
                # This is equal to lat[0:-1]
                #Add padding at the pole to conserve the same dimension as lat.
                list_pads=[]
                for i, dim in enumerate(new_dims):
                    if dim == 'lon':
                        list_pads.append((0,1)) #(begining,end)=(0,1), meaning 1 padding at the end of the array
                    else:
                        list_pads.append((0,0)) #(begining,end)=(0,0), no pad on that axis
                transformed_var=np.pad(transformed_var,list_pads,mode='wrap') #TODO with this method V[0] =V[-1]: can we add cyclic point before destaggering?

                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, coords={'lon':DS['lon']})
                DS[var_name].attrs['long_name'] = '(UNSTAGGERED IN POST-PROCESSING) ' + longname_txt
                DS[var_name].attrs['units'] = units_txt

        # ==============================================================
        #                 PCM Specific Processing
        # ==============================================================

        elif model_type == 'pcm':

            print(f"{Cyan}Processing pcm file")
            #Adding long_name attibutes
            for var_name in DS.data_vars:
                var = DS[var_name]
                if 'title' in var.attrs:
                    var.attrs['long_name'] = var.attrs['title']

            pfull=DS.aps.values + DS.bps.values*ref_press
            DS['pfull']=(['altitude'],  pfull)
            #Adding a pfull variable
            #DS = DS.assign_coords(pfull=altitude)
            DS['pfull'].attrs['long_name']='(ADDED IN POST-PROCESSING), reference pressure'
            DS['pfull'].attrs['units']='Pa'

            #DS[model.phalf]= ak + ref_press*bk
            #DS.phalf.attrs['long_name'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
            #DS.phalf.attrs['description'] = '(ADDED IN POST PROCESSING) pressure at layer interfaces'
            #DS.phalf.attrs['units'] = 'Pa'

        # ==============================================================
        #                START PROCESSING FOR ALL MODELS
        # ==============================================================
        # check that vertical grid starts at toa with highest level at surface
        if DS[model.dim_pfull][0] != DS[model.dim_pfull].min(): # if toa, lev = 0 is surface then flip
            print(f"{Cyan}Current variables at start of processing for all models:\n{list(DS.variables)}{Nclr}\n")
            DS = DS.isel(**{model.dim_pfull: slice(None, None, -1)})
            DS=DS.isel(**{model.dim_phalf: slice(None, None, -1)}) #Also flip phalf,ak, bk
            print(f"{Red}NOTE: all variables flipped along vertical dimension, so that the top of the atmosphere is now index 0")

        # reorder dimensions
        print(f"{Cyan} Transposing variable dimensions to match order expected in CAP")
        DS = DS.transpose(model.dim_time, model.dim_pfull ,model.dim_lat,model.dim_lon, ...)

        # change longitude from -180-179 to 0-360
        if min(DS[model.dim_lon]) < 0:
                tmp = np.array(DS[model.dim_lon])
                tmp = np.where(tmp<0,tmp+360,tmp)
                DS[model.dim_lon] = tmp
                DS = DS.sortby(model.dim_lon)
                DS[model.lon].attrs['long_name']='(MODIFIED POST-PROCESSING) longitude'
                DS[model.lon].attrs['units']='degrees_E'
                print(f"{Red} NOTE: Longitude changed to 0-360E and all variables appropriately reindexed")

        # add scalar axis to areo [time, scalar_axis])
        inpt_dimlist = DS.dims
        # first check if dimensions are correct and don't need to be modified
        if 'scalar_axis' not in inpt_dimlist:           # first see if scalar axis is a dimension
                scalar_axis = DS.assign_coords(scalar_axis=1)
        if DS[model.areo].dims != (model.time,scalar_axis):
                DS[model.areo] = DS[model.areo].expand_dims('scalar_axis', axis=1)
                DS[model.areo].attrs['long_name'] = '(SCALAR AXIS ADDED IN POST-PROCESSING) ' + DS[model.areo].attrs['long_name']

                print(f"{Red}NOTE: scalar axis added to aerocentric longitude")

        # STANDARDIZED VARIABLES NAMES IF REQUESTED
        if args.retain_names:
            print(f"{Purple}Preserving original names for variable and dimensions")
            ext=ext+'_nat'
        else:
            print(f"{Purple}Using standard FV3 names for variables and dimensions")

            #Model has {'ucomp':'U','temp':'T', 'dim_lon'='XLON'...}
            #Create a reversed dictionaries, e.g. {'U':'ucomp','T':'temp'...}
            #to revert to the original variable names before archiving
            #Note that model.() is constructed only from .amescap_profile and may include names not present in file
            model_dims=dict()
            model_vars=dict()

            #Loop over optential dimensions and variables in model()
            for key_i in model.__dict__.keys():
                val=getattr(model,key_i)
                #Potential variables
                if key_i[0:4]!='dim_':
                    #Check if the key is either a variable (e.g. temp) or coordinate (e.g. lat)
                    if val in list(DS.keys())+list(DS.coords):
                        model_vars[val]=key_i
                #Potential dimensions
                else:
                    if val in list(DS.dims):
                        model_dims[val]=key_i[4:]
            #Sort the dictionaries to remove dupplicate:

            #Remove key/val dupplicates: e.g if {'lat':'lat'}, there is no need to rename that variable
            model_dims = {key: val for key, val in model_dims.items() if key != val}
            model_vars = {key: val for key, val in model_vars.items() if key != val}

            #Update dataset with new variables and names
            print(f"{Cyan}Current variables before renaming:\n{list(DS.variables)}{Nclr}\n")
            
            DS=DS.swap_dims(dims_dict=model_dims)
            DS=DS.rename_vars(name_dict=model_vars)
            
            print(f"{Cyan}Renamed variables:\n{list(DS.variables)}{Nclr}\n")
            #Update CAP's internal variables dictionary
            model=reset_FV3_names(model)

        # Output Processed Data to New **atmos_daily.nc File

        # ==============================================================
        # CREATE ATMOS_DAILY, ATMOS_AVERAGE, & ATMOS_DIURN FILES
        # ==============================================================

        if args.bin_average:
            ext=ext+'_average'
            nday=args.bin_average

            # Output Binned Data to New **atmos_average.nc file
            # Figure out number of timesteps per 5 sol
            dt_in = DS[model.dim_time][1]-DS[model.dim_time][0]

            iperday = int(np.round(1/dt_in)) #At least one per day
            if iperday == 0:
                print(f"{Red}***Error***: Operation not permitted because time sampling in file is less than one time step per day")
                break

            combinedN = int(iperday*nday)
            # Coarsen the 'time' dimension by a factor of 5 and average over each window
            DS_average = DS.coarsen(**{model.dim_time:combinedN}).mean()

            # Update the time coordinate attribute
            DS_average[model.dim_time].attrs['long_name'] = 'time averaged over %s sols'%(nday)

            # Create New File, set time dimension as unlimitted
            fullnameOUT = fullnameIN[:-3]+ext+'.nc'
            DS_average.to_netcdf(fullnameOUT,unlimited_dims=model.dim_time,format='NETCDF4_CLASSIC')


        elif args.bin_diurn:
            ext=ext+'_diurn'

            #Custom number of sol
            if args.bin_average:
                nday=args.bin_average
            else:
                nday=5

            dt_in = DS[model.dim_time][1]-DS[model.dim_time][0]
            iperday = int(np.round(1/dt_in))
            if iperday == 0:
                print(f"{Red}***Error***: Operation not permitted because time sampling in file is less than one time step per day")
                break
            combinedN = int(iperday*nday)
            
            # Output Binned Data to New  **atmos_diurn.nc file
            # create a new time of day dimension
            tod_name = 'time_of_day_%02d' %(iperday)
            days = len(DS[model.dim_time])/iperday

            # initialize the new dataset
            DS_diurn = None

            # loop through coarsened grid, slicing time dimension in 5 sol groups
            for i in range(0, int(days/nday)):

                # slice original dataset in 5 sol periods
                downselect = DS.isel(**{model.dim_time:slice(i*combinedN, i*combinedN+combinedN)})

                # rename the time dimension to the time of day and find the LT equivalent
                downselect = downselect.rename({model.dim_time: tod_name})
                downselect[tod_name] = np.mod(downselect[tod_name]*24, 24).values

                # Group up all instances of the same LT & take the mean
                idx = downselect.groupby(tod_name).mean()

                # add back in the time dimensionn
                idx = idx.expand_dims({model.dim_time: [i]})  # Add 'time' dimension with integer values

                # concatenate into new diurn array with a LT and time dimension (stack along time)
                if DS_diurn is None:
                    DS_diurn = idx
                else:
                    DS_diurn = xr.concat([DS_diurn, idx], dim=model.dim_time)

            # replace the time dimension with the time dimension from DS_average
            DS_time=DS[model.dim_time]
            DS_time_avg= DS_time.coarsen(**{model.dim_time:combinedN}).mean()

            DS_diurn[model.dim_time] = DS_time_avg[model.dim_time]
            # Update the time coordinate attribute
            DS_diurn[model.dim_time].attrs['long_name'] = 'time averaged over %s sols'%(nday)

            # Create New File, set time dimension as unlimitted
            fullnameOUT = fullnameIN[:-3]+ext+'.nc'
            DS_diurn.to_netcdf(fullnameOUT,unlimited_dims=model.dim_time,format='NETCDF4_CLASSIC')

        else:
            ext=ext+'_daily'
            fullnameOUT = fullnameIN[:-3]+ext+'.nc'
            DS.to_netcdf(fullnameOUT,unlimited_dims=model.dim_time,format='NETCDF4_CLASSIC')
        print(f"{Cyan}{fullnameOUT} was created")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
