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
parser=argparse.ArgumentParser(
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
        f"Works on non-MGCM files only.\n"
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
    ext = '' # Initialize empty extension
    if args.gcm_name not in ['marswrf', 'openmars', 'pcm', 'emars']:
        print(f"{Yellow}***Notice***  No operation requested. Use "
              f"'-gcm' and specify openmars, marswrf, pcm, emars")
        exit() # Exit cleanly

    print(f"Running MarsFormat with args: {args}")
    print(f"Current working directory: {os.getcwd()}")
    print(f"Files in input_file: {[f.name for f in args.input_file]}")
    print(f"File exists check: "
          f"{all(os.path.exists(f.name) for f in args.input_file)}")

    path2data = os.getcwd()

    # Load all of the netcdf files
    file_list = [f.name for f in args.input_file]
    model_type = args.gcm_name  # e.g. 'marswrf'
    for filei in file_list:
        # Use os.path.join for platform-independent path handling
        if os.path.isabs(filei):
            fullnameIN = filei
        else:
            fullnameIN = os.path.join(path2data, filei)

        print('Processing...')
        # Load model variables, dimensions
        fNcdf = Dataset(fullnameIN, 'r')
        model = read_variable_dict_amescap_profile(fNcdf)
        fNcdf.close()

        print(f"{Cyan}Reading model attributes from ~.amescap_profile:")
        print(f"{Cyan}{vars(model)}") # Print attributes

        # Open dataset with xarray
        DS = xr.open_dataset(fullnameIN, decode_times=False)

        # --------------------------------------------------------------
        #                    MarsWRF Specific Processing
        # --------------------------------------------------------------
        if model_type == 'marswrf':
            # First, save all variable descriptions in attrs longname
            print(f"{Cyan}Current variables at top of marswrf "
                  f"processing: \n{list(DS.variables)}{Nclr}\n")
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
                    raise KeyError(
                        f"Longitude dimension {model.lon} not found"
                        )

                # Time conversion (minutes to days)
                time = (DS[model.time] / 60 / 24) if 'time' in DS else None

                # Handle latitude and longitude 
                if len(DS[model.lat].shape) > 1:
                    lat = DS[model.lat][0, :, 0]
                    lon = DS[model.lon][0, 0, :]
                else:
                    lat = DS[model.lat]
                    lon = DS[model.lon]
                
                # Convert longitudes to 0-360
                lon360 = (lon + 360)%360

                # Update coordinates 
                if time is not None:
                    DS[model.time] = time
                DS[model.lon] = lon360
                DS[model.lat] = lat

                # Derive pressure levels 
                if (hasattr(DS, 'P_TOP') and 
                    hasattr(DS, 'ZNU') and 
                    hasattr(DS, 'P0')):
                    pfull = DS.P_TOP[0] + DS.ZNU[0,:]*DS.P0
                    phalf = DS.P_TOP[0] + DS.ZNW[0,:]*DS.P0
                    
                    DS = DS.assign_coords(phalf=phalf, pfull=pfull)
                    DS.phalf.attrs['long_name'] = (
                        '(ADDED POST PROCESSING) pressure at layer interfaces'
                        )
                    DS.phalf.attrs['description'] = (
                        '(ADDED POST PROCESSING) pressure at layer interfaces'
                        )
                    DS.phalf.attrs['units'] = ('Pa')

                    # Derive ak, bk
                    ak = np.zeros(len(DS.phalf))
                    bk = np.zeros(len(DS.phalf))
                    ak[-1] = DS.P_TOP[0] # MarsWRF pressure increases w/N
                    bk[:] = DS.ZNW[0,:]

                    DS['ak'], DS['bk'] = ak, bk
                    DS['ak'].attrs['long_name'] = ('(ADDED POST PROCESSING)')
                    DS['bk'].attrs['long_name'] = ('(ADDED POST PROCESSING)')

                # Calculate layer heights (if PH and PHB exist)
                if ('PH' in DS and 
                    'PHB' in DS and 
                    'HGT' in DS and 
                    'G' in DS):
                    zagl_lvl = ((DS.PH[:, :-1, :, :] + DS.PHB[0, :-1, :, :]) 
                                /DS.G - DS.HGT[0, :, :])
                    zfull3D = (
                        0.5*(zagl_lvl[:, :-1, :, :] + zagl_lvl[:, 1:, :, :])
                        )

                # Derive atmospheric temperature
                if all(
                    key in DS for key in ['T', 'T0', 'CP', 'R_D', 'P0', 'PB']
                    ):
                    gamma = DS.CP/(DS.CP - DS.R_D)
                    pfull3D = DS.P_TOP + DS.PB[0,:]
                    DS['T'] = ((DS.T + DS.T0)
                               *(pfull3D / DS.P0)**((gamma-1.) / gamma))
                    DS['T'].attrs['description'] = (
                        '(MODIFIED POST PROCESSING) Temperature'
                        )
                    DS['T'].attrs['long_name'] = (
                        '(MODIFIED POST PROCESSING) Temperature'
                        )
                    DS['T'].attrs['units'] = 'K'

                # Drop 'Times' variable if present (non-numerical values)
                if 'Times' in DS:
                    DS = DS.drop_vars("Times")

            except Exception as e:
                print(f"{Red}Error processing MarsWRF data: {e}{Nclr}")
                print(f"{Red}Available dimensions: {list(DS.dims)}{Nclr}")
                print(f"{Red}Available variables: {list(DS.variables)}{Nclr}")
                raise
            
        # --------------------------------------------------------------
        #                  OpenMars Specific Processing
        # --------------------------------------------------------------
        elif model_type == 'openmars':
            # First save all variable FIELDNAM in attrs longname
            var_list= list(DS.data_vars)+list(DS.coords)
            for var_name in var_list:
                var = DS[var_name]
                if 'FIELDNAM' in var.attrs:
                    var.attrs['long_name'] = var.attrs['FIELDNAM']
                if 'UNITS' in var.attrs:
                    var.attrs['units'] = var.attrs['UNITS']

            # Define Coordinates for New DataFrame
            time = DS[model.dim_time] # min since simulation start [m]
            lat = DS[model.dim_lat] # Replace DS.lat
            lon = DS[model.dim_lon]

            DS = DS.assign(pfull = DS[model.dim_pfull]*ref_press)

            DS['pfull'].attrs['long_name']='(ADDED POST-PROCESSING) reference pressure'
            DS['pfull'].attrs['units']='Pa'


            # add ak,bk as variables
            # add p_half dimensions as vertical grid coordinate

            # Compute sigma values. Swap the sigma array upside down 
            # twice  with [::-1] since layers_mid_point_to_boundary() 
            # needs sigma[0]=0, sigma[-1]=1) and then to reorganize the 
            # array in the original openMars format with sigma[0]=1, 
            # sigma[-1]=0
            bk = layers_mid_point_to_boundary(DS[model.dim_pfull][::-1], 1.)[::-1]
            ak = np.zeros(len(DS[model.dim_pfull]) + 1)

            DS[model.phalf] = (ak + ref_press*bk)
            DS.phalf.attrs['long_name'] = (
                '(ADDED POST PROCESSING) pressure at layer interfaces'
                )
            DS.phalf.attrs['description'] = (
                '(ADDED POST PROCESSING) pressure at layer interfaces'
                )
            DS.phalf.attrs['units'] = ('Pa')

            DS = DS.assign(bk=(model.dim_phalf, np.array(bk)))
            DS = DS.assign(ak=(model.dim_phalf, np.zeros(len(DS[model.dim_pfull]) + 1)))
            
            # Update Variable Description & Longname
            DS['ak'].attrs['long_name'] = ('(ADDED POST PROCESSING) pressure'
                                           'part of the hybrid coordinate')
            DS['ak'].attrs['units'] = ('Pa')
            DS['bk'].attrs['long_name'] = (
                '(ADDED POST PROCESSING) vertical coordinate sigma value'
                )
            DS['bk'].attrs['units'] = ('None')

        # --------------------------------------------------------------
        #                 Emars Specific Processing
        # --------------------------------------------------------------
        elif model_type == 'emars':
            # Interpolate U, V, onto Regular Mass Grid (from staggered)
            print(f"{Cyan}Interpolating Staggered Variables to Standard Grid")

            variables_with_latu = [var for var in DS.variables if 'latu' in DS[var].dims]
            variables_with_lonv = [var for var in DS.variables if 'lonv' in DS[var].dims]
            
            #Replace time array in hours to sols:
            DS[model.time] = DS[model.time].values/24.
            DS[model.time].attrs['long_name'] = 'time'  
            DS[model.time].attrs['units'] = 'days since 0000-00-00 00:00:00'
            
            # Loop through, and unstagger. the dims_list process finds the dimensions of the variable and replaces west_east_stag with west_east
            print('     From latu to lat: ' + ', '.join(variables_with_latu ))
            for var_name in variables_with_latu:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt = var.long_name
                units_txt = var.units

                # Replace latu dims with lat
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'latu':
                        dims_list[i] = 'lat'
                        break # Stop the loop when replacement made
                new_dims = tuple(dims_list)

                #TODO Using 'values' as var.isel(latu=slice(None, -1)).values and var.isel(latu=slice(None, -1)) returns array of different size
                #TODO the following reproduce the 'lat' array, but not if the keyword values is ommited
                #latu1_val=ds.latu.isel(latu=slice(None, -1)).values
                #latu2_val=ds.latu.isel(latu=slice(1,  None)).values

                #newlat_val= 0.5 * (latu1 + latu2)
                #newlat_val=np.append(newlat_val,0)
                #newlat_val ==lat
                transformed_var = (
                    0.5*(var.isel(latu=slice(None, -1)).values 
                         + var.isel(latu=slice(1, None)).values)
                    )

                # This is equal to lat[0:-1]
                # Add padding at the pole to conserve the same dimension as lat
                list_pads = []
                for i, dim in enumerate(new_dims):
                    if dim == 'lat':
                        list_pads.append((0,1))
                        # (begining, end)=(0, 1), meaning 1 padding at 
                        # the end of the array
                    else:
                        list_pads.append((0,0)) 
                        # (begining, end)=(0, 0), no pad on that axis
                
                transformed_var = np.pad(transformed_var, list_pads, 
                                         mode='constant', constant_values=0)
                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, 
                                            coords={'lat':DS['lat']})
                DS[var_name].attrs['long_name'] = (
                    f'(UNSTAGGERED IN POST-PROCESSING) {longname_txt}'
                    )
                DS[var_name].attrs['units'] = units_txt
            
            
            for var_name in variables_with_lonv:
                var = getattr(DS, var_name)
                dims = var.dims
                longname_txt = var.long_name
                units_txt = var.units

                # Replace lonv in dimensions with lon
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == 'lonv':
                        dims_list[i] = 'lon'
                        break # Stop loop once the replacement is made
                new_dims = tuple(dims_list)

                transformed_var = (
                    0.5 * (var.isel(lonv=slice(None, -1)).values 
                           + var.isel(lonv=slice(1, None)).values)
                    )
                # This is equal to lon[0:-1]
                # Add padding 
                list_pads = []
                for i, dim in enumerate(new_dims):
                    if dim == 'lon':
                        list_pads.append((0, 1)) 
                        # (begining, end)=(0, 1), meaning 1 padding at 
                        # the end of the array
                    else:
                        list_pads.append((0, 0)) 
                        # (begining, end)=(0, 0), no pad on that axis
                transformed_var = np.pad(transformed_var, list_pads, 
                                         mode='wrap')
                #TODO with this method V[0] =V[-1]: can we add cyclic point before destaggering?

                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, 
                                            coords={'lon':DS['lon']})
                DS[var_name].attrs['long_name'] = (f'(UNSTAGGERED IN POST-PROCESSING) {longname_txt}')
                DS[var_name].attrs['units'] = units_txt
            DS.drop_vars(['latu','lonv'])
        # --------------------------------------------------------------
        #                 PCM Specific Processing
        # --------------------------------------------------------------
        elif model_type == 'pcm':
            print(f"{Cyan}Processing pcm file")
            # Adding long_name attibutes
            for var_name in DS.data_vars:
                var = DS[var_name]
                if 'title' in var.attrs:
                    var.attrs['long_name'] = var.attrs['title']

            pfull = (DS.aps.values + DS.bps.values*ref_press)
            DS['pfull'] = (['altitude'],  pfull)
            # Adding a pfull variable
            # DS = DS.assign_coords(pfull=altitude)
            DS['pfull'].attrs['long_name'] = (
                '(ADDED POST-PROCESSING), reference pressure'
                )
            DS['pfull'].attrs['units'] = 'Pa'

            # DS[model.phalf] = (ak + ref_press*bk)
            # DS.phalf.attrs['long_name'] = (
                # '(ADDED POST PROCESSING) pressure at layer interfaces'
                # )
            # DS.phalf.attrs['description'] = (
                # '(ADDED POST PROCESSING) pressure at layer interfaces'
                # )
            # DS.phalf.attrs['units'] = 'Pa'

        # --------------------------------------------------------------
        #                START PROCESSING FOR ALL MODELS
        # --------------------------------------------------------------
        # Check that vertical grid starts at TOA w/ largest level at surface
        if DS[model.dim_pfull][0] != DS[model.dim_pfull].min():
            # If TOA, lev = 0 is surface then flip
            print(f"{Cyan}Current variables at start of processing for all "
                  f"models:\n{list(DS.variables)}{Nclr}\n")
            DS = DS.isel(**{model.dim_pfull: slice(None, None, -1)})
            # Flip phalf, ak, bk:
            DS = DS.isel(**{model.dim_phalf: slice(None, None, -1)})
            print(f"{Red}NOTE: all variables flipped along vertical "
                  f"dimension. Top of the atmosphere is now index = 0")

        # Reorder dimensions
        print(f"{Cyan} Transposing variable dimensions to match order "
              f"expected in CAP")
        DS = DS.transpose(model.dim_time, model.dim_pfull, model.dim_lat,
                          model.dim_lon, ...)

        # Change longitude from -180-179 to 0-360
        if min(DS[model.dim_lon]) < 0:
            tmp = np.array(DS[model.dim_lon])
            tmp = np.where(tmp<0, tmp + 360, tmp)
            DS[model.dim_lon] = tmp
            DS = DS.sortby(model.dim_lon)
            DS[model.lon].attrs['long_name'] = (
                '(MODIFIED POST-PROCESSING) longitude'
                )
            DS[model.lon].attrs['units'] = ('degrees_E')
            print(f"{Red} NOTE: Longitude changed to 0-360E and all variables "
                  f"appropriately reindexed")

        # Add scalar axis to areo [time, scalar_axis])
        inpt_dimlist = DS.dims
        # First check if dims are correct - don't need to be modified
        if 'scalar_axis' not in inpt_dimlist:
            # If scalar axis is a dimension
            scalar_axis = DS.assign_coords(scalar_axis=1)
        if DS[model.areo].dims != (model.time,scalar_axis):
            DS[model.areo] = DS[model.areo].expand_dims('scalar_axis', axis=1)
            DS[model.areo].attrs['long_name'] = '(SCALAR AXIS ADDED POST-PROCESSING)'
             
            print(f"{Red}NOTE: scalar axis added to aerocentric longitude")

        # Standardize variables names if requested
        if args.retain_names:
            print(f"{Purple}Preserving original names for variable and "
                  f"dimensions")
            ext = f'{ext}_nat'
        else:
            print(f"{Purple}Using standard FV3 names for variables and "
                  f"dimensions")

            # Model has {'ucomp':'U','temp':'T', 'dim_lon'='XLON'...}
            # Create a reversed dictionary, e.g. {'U':'ucomp','T':'temp'...}
            # to revert to the original variable names before archiving
            # Note that model.() is constructed only from 
            # .amescap_profile and may include names not present in file
            model_dims = dict()
            model_vars = dict()

            # Loop over optential dims and vars in model
            for key_i in model.__dict__.keys():
                val = getattr(model,key_i)
                if key_i[0:4] != 'dim_':
                    # Potential variables
                    if val in (list(DS.keys()) + list(DS.coords)):
                        # Check if the key is a variable (e.g. temp) 
                        # or coordinate (e.g. lat)
                        model_vars[val] = key_i
                else:
                    # Potential dimensions
                    if val in list(DS.dims):
                        model_dims[val] = key_i[4:]
            
            # Sort the dictionaries to remove duplicates
            # Remove key/val duplicates: e.g if {'lat':'lat'}, there is 
            # no need to rename that variable
            model_dims = {key: val for key, val in model_dims.items() if key != val}
            model_vars = {key: val for key, val in model_vars.items() if key != val}
            
            DS = DS.swap_dims(dims_dict = model_dims)
            DS = DS.rename_vars(name_dict = model_vars)
            
            print(f"{Cyan}Renamed variables:\n{list(DS.variables)}{Nclr}\n")
            # Update CAP's internal variables dictionary
            model = reset_FV3_names(model)
            
        # --------------------------------------------------------------
        # CREATE ATMOS_DAILY, ATMOS_AVERAGE, & ATMOS_DIURN FILES
        # --------------------------------------------------------------
        if args.bin_average:
            ext = f'{ext}_average'
            nday = args.bin_average

            # Output Binned Data to New **atmos_average.nc file
            # Figure out number of timesteps per 5 sol
            dt_in = DS[model.dim_time][1] - DS[model.dim_time][0]

            iperday = int(np.round(1/dt_in)) # at least one per day
            if iperday == 0:
                print(f"{Red}***Error***: Operation not permitted because "
                      f"time sampling in file < one time step per day")
                break

            combinedN = int(iperday*nday)
            # Coarsen the 'time' dimension by a factor of 5 and average 
            # over each window
            DS_average = DS.coarsen(**{model.dim_time:combinedN}).mean()

            # Update the time coordinate attribute
            DS_average[model.dim_time].attrs['long_name'] = (
                f'time averaged over {nday} sols'
                )

            # Create New File, set time dimension as unlimitted
            base_name = os.path.splitext(fullnameIN)[0]
            fullnameOUT = f"{base_name}{ext}.nc"
            DS_average.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time, 
                                 format='NETCDF4_CLASSIC')

        elif args.bin_diurn:
            ext = f'{ext}_diurn'

            # Custom number of sols
            if args.bin_average:
                nday = args.bin_average
            else:
                nday = 5

            dt_in = DS[model.dim_time][1] - DS[model.dim_time][0]
            iperday = int(np.round(1/dt_in))
            if iperday == 0:
                print(f"{Red}***Error***: Operation not permitted because "
                      f"time sampling in file < one time step per day")
                break
            
            combinedN = int(iperday * nday)
            
            # Output Binned Data to New  **atmos_diurn.nc file
            # create a new time of day dimension
            tod_name = f"time_of_day_{iperday:02d}"
            days = len(DS[model.dim_time]) / iperday

            # Initialize the new dataset
            DS_diurn = None

            for i in range(0, int(days/nday)):
                # Slice original dataset in 5 sol increments
                downselect = (
                    DS.isel(**{model.dim_time:slice(i*combinedN, 
                                                    i*combinedN + combinedN)})
                    )

                # Rename time dimension to time of day and find the 
                # local time equivalent
                downselect = downselect.rename({model.dim_time: tod_name})
                downselect[tod_name] = (
                    np.mod(downselect[tod_name]*24, 24).values
                    )

                # Average all instances of the same local time
                idx = downselect.groupby(tod_name).mean()

                # Add back in the time dimensionn
                idx = idx.expand_dims({model.dim_time: [i]})

                # Concatenate into new diurn array with a local time 
                # and time dimension (stack along time)
                if DS_diurn is None:
                    DS_diurn = idx
                else:
                    DS_diurn = xr.concat([DS_diurn, idx], dim=model.dim_time)
            
            # ==== Overwrite the ak, bk arrays=== [AK]
            #For some reason I can't track down why the ak('phalf') and bk('phalf') 
            # turn into ak ('time', 'time_of_day_12', 'phalf'), which messes
            # the pressure interpolation. 
            DS_diurn['ak']=DS['ak']
            DS_diurn['bk']=DS['bk']
            
            # replace the time dimension with the time dimension from DS_average
            time_DS=DS[model.dim_time]
            time_avg_DS= time_DS.coarsen(**{model.dim_time:combinedN}).mean()

            DS_diurn[model.dim_time] = time_avg_DS[model.dim_time]

            # Update the time coordinate attribute
            DS_diurn[model.dim_time].attrs['long_name'] = (
                f'time averaged over {nday} sols'
                )

            # Create New File, set time dimension as unlimitted
            fullnameOUT = f'{fullnameIN[:-3]}{ext}.nc'
            DS_diurn.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time, 
                               format='NETCDF4_CLASSIC')

        else:
            ext = f'{ext}_daily'
            fullnameOUT = f'{fullnameIN[:-3]}{ext}.nc'
            DS.to_netcdf(fullnameOUT, unlimited_dims=model.dim_time, 
                         format='NETCDF4_CLASSIC')
        print(f"{Cyan}{fullnameOUT} was created")

if __name__ == "__main__":
    exit_code = main()
    sys.exit(exit_code)
