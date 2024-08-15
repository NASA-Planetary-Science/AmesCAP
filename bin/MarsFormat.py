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
from amescap.Script_utils import (Cyan, Yellow, Nclr, Purple, Red, Green)

# Load generic Python modules
import argparse     # Parse arguments
import os           # Access operating system functions
import numpy as np
import xarray as xr
import pandas as pd
from netCDF4 import Dataset

# Load amesCAP modules
from amescap.Script_utils import (
    read_variable_dict_amescap_profile, filter_vars, get_longname_units
)
from amescap.FV3_utils import (
    daily_to_average, layers_mid_point_to_boundary
)
from amescap.Ncdf_wrapper import Ncdf

xr.set_options(keep_attrs = True)

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

parser.add_argument("input_file", nargs="+",
    help=(f"A netCDF file or list of netCDF files.\n\n"))

parser.add_argument("-t", "--type", type=str, default="legacy",
    help=(
        f"Model that the incoming data comes from. Can be 'openmars', "
        f"'marswrf' or 'legacy' [DEFAULT].\n"
        f"{Green}Usage:\n"
        f"> MarsFormat.py ****.nc\n"
        f"> MarsFormat.py ****.nc -t openmars\n"
        f"{Nclr}\n\n"
    )
)

def main():
    if not (parser.parse_args().type):
        print(f"{Yellow} ***Notice***  No operation requested. Use "
                f"'-type' and specify openmars, marswrf, legacy {Nclr}")
        exit()

    path2data = os.getcwd()

    # Load all of the netcdf files
    file_list = parser.parse_args().input_file
    model_type = parser.parse_args().type
    for filei in file_list:
        # Add path unless full path is provided
        if not ("/" in filei):
            fullnameIN = f"{path2data}/{filei}"
        else:
            fullnameIN = filei
        fullnameOUT = f"{fullnameIN[:-3]}_atmos_daily.nc"

        print("Processing...")
        # Load model variables, dimensions
        fNcdf = Dataset(fullnameIN, "r")
        # model = read_variable_dict_amescap_profile(fNcdf)
        fNcdf.close()
        #print(f"{Cyan}Reading model attributes from ~.amescap_profile:\n"
        #      f"{vars(model)}{Nclr}")
        # dataDIR = (f"{path}{filename}.nc")
        DS = xr.open_dataset(fullnameIN, decode_times=False)

        # ==============================================================
        #               MarsWRF-Specific Processing
        # ==============================================================
        if model_type == "marswrf":
            #TODO long_name is "description" for MarsWRF
            # print("Input File content (description) and (description) "
            #       "attibutes:")
            # print("------")
            # for ivar in DS.keys():
            #     print(ivar, DS[ivar].attrs["description"],
            #           DS[ivar].attrs["units"])
            #     print("------")
            
            # First save all variable descriptions in long_name
            for var_name in DS.data_vars:
                var = DS[var_name]
                if "description" in var.attrs:
                    var.attrs["long_name"] = var.attrs["description"] 
            
            # Reformat Dimension Variables/Coords as Needed
            time = DS["time"]/60/24
            lat = DS["lat"][0, :, 0]
            lon2D = DS["lon"][0, :]
            lon = np.squeeze(lon2D[0, :])
            DS["lon"] = lon
            DS["lat"] = lat
            DS["time"] = time
            
            # Derive half and full reference pressure levels (Pa)
            pfull = DS.P_TOP[0] + DS.ZNU[0,:]*DS.P0 
            phalf = DS.P_TOP[0] + DS.ZNW[0,:]*DS.P0
            DS = DS.assign_coords(phalf = phalf, pfull = pfull)
            DS.phalf.attrs["long_name"] = (
                "(ADDED POST-PROCESSING) pressure at layer interfaces"
                )
            DS.phalf.attrs["description"] = (
                "(ADDED POST-PROCESSING) pressure at layer interfaces"
                )
            DS.phalf.attrs["units"] = "Pa"

            # Update dimensions
            DS = DS.assign_coords(dimensions = "phalf")

            # Update Variable long_name
            DS["lon"].attrs["long_name"] = (
                "(MODIFIED POST-PROCESSING) "
                + DS["lon"].attrs["description"])
            DS["lat"].attrs["long_name"] = (
                "(MODIFIED POST-PROCESSING) "
                + DS["lat"].attrs["description"])
            DS["time"].attrs["long_name"] = (
                "(MODIFIED POST-PROCESSING) days since simulation start "
                "(time/60/24)")
            
            # Update Variable Description & Unit
            DS["lon"].attrs["description"] = (
                "(MODIFIED POST-PROCESSING) "
                + DS["lon"].attrs["description"])
            DS["lat"].attrs["description"] = (
                "(MODIFIED POST-PROCESSING) "
                + DS["lat"].attrs["description"])
            DS["time"].attrs["description"] = (
                "(MODIFIED POST-PROCESSING) days since simulation start "
                "(time/60/24)")
            DS["time"].attrs["units"] = (
                "(MODIFIED POST-PROCESSING) days")

            # ==========================================================
            # Derive ak, bk
            # ==========================================================
            ak = np.zeros(len(DS.phalf))
            bk = np.zeros(len(DS.phalf))
            
            # MarsWRF has pressure increasing with N
            ak[-1] = DS.P_TOP[0]
            bk[:] = DS.ZNW[0, :]

            DS["ak"] = ak
            DS["bk"] = bk

            # Update Variable Description & long_name
            DS["ak"].attrs["long_name"] = "(ADDED POST-PROCESSING)"
            DS["bk"].attrs["long_name"] = "(ADDED POST-PROCESSING)"
            DS["ak"].attrs["description"] = "(ADDED POST-PROCESSING)"
            DS["bk"].attrs["description"] = "(ADDED POST-PROCESSING)"
        
            # ==========================================================
            # Calculate *level* heights above the surface
            # ==========================================================
            zagl_lvl = ((DS.PH[:, :-1, :, :]+DS.PHB[0, :-1, :, :]) 
                        / DS.G 
                        - DS.HGT[0, :, :])

            # ==========================================================
            # Find layer pressures [Pa]
            # ==========================================================
            try:
                # Perturbation pressure + base state pressure 
                # (time-invariant)
                pfull3D = DS.P_TOP + DS.PB[0, :] 
            except NameError:
                pfull3D = DS["ps"][:, :-1, :-1] * DS.ZNU[:, :-1]

            # ==========================================================
            # Derive atmospheric temperature [K]
            # ==========================================================
            gamma = DS.CP / (DS.CP-DS.R_D)
            temp = (DS.T+DS.T0) * (pfull3D/DS.P0) ** ((gamma-1.) / gamma)
            DS = DS.assign(temp = temp)
            DS["temp"].attrs["description"] = (
                "(ADDED POST-PROCESSING) Temperature")
            DS["temp"].attrs["long_name"] = (
                "(ADDED POST-PROCESSING) Temperature")
            DS["temp"].attrs["units"] = "K"

            # ==========================================================
            # Interpolate U, V, W, Zfull onto regular mass grid 
            # (from staggered)
            # ==========================================================
            # For variables staggered in X (lon) [t,z,y,x'] -> regular 
            # mass grid [t,z,y,x]:
            # Step 1: Identify variables with the west_east_stag 
            # dimension where _U is NOT in the variable name (indicates 
            # a staggered grid identifier like the staggered latitude)
            variables_with_west_east_stag = (
                [var for var in DS.variables if "west_east_stag" in 
                 DS[var].dims and "_U" not in var]
                )

            print(f"{Cyan}Interpolating Staggered Variables to Standard Grid"
                  f"{Nclr}")
            # Loop through and unstagger. dims_list finds the dimensions
            # of the variable and replaces west_east_stag with west_east
            print(f"     {Purple}From west_east_stag to west_east: "
                  + ", ".join(variables_with_west_east_stag)
                  + f"{Nclr}")
            
            for var_name in variables_with_west_east_stag:
                var = getattr(DS, var_name)
                dims = var.dims
                # Replace west_east_stag with west_east
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == "west_east_stag":
                        dims_list[i] = "west_east"
                        # Stop the loop once the replacement is made
                        break 
                new_dims = tuple(dims_list)
            
                transformed_var = (
                    0.5 * (var.sel(west_east_stag = slice(None, -1)) 
                           + var.sel(west_east_stag = slice(1, None)))
                    )
                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, 
                                            coords={"XLAT":DS["XLAT"]})
                
                DS[var_name].attrs["description"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])
                DS[var_name].attrs["long_name"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])
                DS[var_name].attrs["stagger"] = "USTAGGERED IN POST-PROCESSING"
        
            # For variables staggered y (lat) [t, z, y", x] 
            # -> regular mass grid [t, z, y, x]: 
            variables_with_south_north_stag = (
                [var for var in DS.variables if "south_north_stag" in 
                 DS[var].dims and "_V" not in var])
            print(f"     {Purple}From south_north_stag to south_north: "
                  + ", ".join(variables_with_south_north_stag)
                  + f"{Nclr}")
            
            for var_name in variables_with_south_north_stag:
                var = getattr(DS, var_name)
                dims = var.dims
                # Replace west_east_stag with west_east
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == "south_north_stag":
                        dims_list[i] = "south_north"
                        # Stop the loop once the replacement is made
                        break 
                new_dims = tuple(dims_list)
                
                transformed_var = (
                    0.5 * (var.sel(south_north_stag = slice(None, -1)) 
                           + var.sel(south_north_stag = slice(1, None)))
                    )
                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims, 
                                            coords={"XLONG":DS["XLONG"]})
                
                DS[var_name].attrs["description"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])

                DS[var_name].attrs["long_name"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])
                DS[var_name].attrs["stagger"] = "USTAGGERED IN POST-PROCESSING"
                
            # For variables staggered p/z (height) [t, z', y, x] 
            # -> regular mass grid [t, z, y, x]:
            variables_with_bottom_top_stag = (
                [var for var in DS.variables if "bottom_top_stag" in 
                 DS[var].dims and "ZNW" not in var and "phalf" not in var])
            print(f"     {Purple}From bottom_top_stag to bottom_top: "
                  + ", ".join(variables_with_bottom_top_stag)
                  + f"{Nclr}")
            
            for var_name in variables_with_bottom_top_stag:
                var = getattr(DS, var_name)
                dims = var.dims
                # Replace bottom_top_stag
                dims_list = list(dims)
                for i, dim in enumerate(dims_list):
                    if dim == "bottom_top_stag":
                        dims_list[i] = "bottom_top"
                        # Stop the loop once the replacement is made
                        break 
                new_dims = tuple(dims_list)
                transformed_var = (
                    0.5 * (var.sel(bottom_top_stag = slice(None, -1)) 
                           + var.sel(bottom_top_stag = slice(1, None)))
                    )
                DS[var_name] = xr.DataArray(transformed_var, dims=new_dims)
                DS[var_name].attrs["description"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])
                DS[var_name].attrs["long_name"] = (
                    "(UNSTAGGERED IN POST-PROCESSING) " 
                    + DS[var_name].attrs["description"])
                DS[var_name].attrs["stagger"] = "USTAGGERED IN POST-PROCESSING"
            
            # ALSO INTERPOLATE TO FIND *LAYER* HEIGHTS ABOVE SURFACE 
            # (i.e., above topography; m)
            zfull3D = (0.5 * (zagl_lvl[:, :-1, :, :] + zagl_lvl[:, 1:, :, :]))

        # ==============================================================
        #                     OpenMars-Specific Processing
        # ==============================================================
        elif model_type == "openmars":
            
            # First save all variable FIELDNAM in attrs long_name
            for var_name in DS.data_vars:
                var = DS[var_name]
                if "FIELDNAM" in var.attrs:
                    var.attrs["long_name"] = var.attrs["FIELDNAM"]

            # Define Coordinates for New DataFrame
            #TODO this is added on to create ak/bk
            ref_press = 720
            time = DS["dim_time"] # Minutes since simulation start
            lat = DS["dim_lat"] # Replace DS.lat
            lon = DS["dim_lon"]

            # Derive half and full reference pressure levels (Pa)
            DS = DS.assign(pfull = DS["dim_pfull"] * ref_press)
            DS["pfull"].attrs["FIELDNAM"] = ("(MODIFIED IN POST-PROCESSING) " 
                                             + DS["pfull"].attrs["FIELDNAM"])
            DS["pfull"].attrs["long_name"] = ("(MODIFIED IN POST-PROCESSING) "
                                              + DS["pfull"].attrs["FIELDNAM"])

            # add ak,bk as variables
            # add p_half dimensions as vertical grid coordinate 

            # Compute sigma values. Swap the sigma array upside down 
            # twice  with [::-1] since the layers_mid_point_to_boundary 
            # needs sigma[0] = 0, sigma[-1] = 1) and then to reorganize 
            # the array in the original openMars format with 
            # sigma[0] = 1, sigma[-1] = 0
            bk = layers_mid_point_to_boundary(
                DS["dim_pfull"][::-1], 1.)[::-1]
            ak = np.zeros(len(DS["dim_pfull"]) + 1)

            DS["phalf"] = ak + ref_press*bk
            DS.phalf.attrs["long_name"] = (
                "(ADDED POST-PROCESSING) pressure at layer interfaces")
            DS.phalf.attrs["description"] = (
                "(ADDED POST-PROCESSING) pressure at layer interfaces")
            DS.phalf.attrs["units"] = "Pa"
            
            DS = DS.assign(bk = ("dim_phalf", np.array(bk)))
            DS = DS.assign(ak = ("dim_phalf", 
                                 np.zeros(len(DS["dim_pfull"]) + 1)))

            # Update Variable Description & long_name
            DS["ak"].attrs["long_name"] = "(ADDED POST-PROCESSING)"
            DS["bk"].attrs["long_name"] = "(ADDED POST-PROCESSING)"
            DS["ak"].attrs["FIELDNAM"] = "(ADDED POST-PROCESSING)"
            DS["bk"].attrs["FIELDNAM"] = "(ADDED POST-PROCESSING)"

        # ==========================================================
        # START PROCESSING FOR ALL MODELS
        # ==========================================================

        # Check that vertical grid starts at TOA with highest level at 
        # surface
        if DS["dim_pfull"][0] != DS["dim_pfull"].min(): 
            # If TOA, lev = 0 = surface, flip dimensions
            DS = DS.isel(**{"dim_pfull": slice(None, None, -1)})
            # Also flip phalf, ak, & bk:
            DS = DS.isel(**{"dim_phalf": slice(None, None, -1)}) 
            print(f"{Red}NOTE: all variables flipped along vertical "
                  f"dimension, so that the top of the atmosphere is now index "
                  f"0{Nclr}")

        # Reorder dimensions
        print(f"{Cyan}Transposing variable dimensions to match order expected "
              f"in CAP{Nclr}") 
        DS = DS.transpose("dim_time", "dim_pfull", 
                          "dim_lat", "dim_lon", ...)

        # Change longitude from -180-179 -> 0-360
        if min(DS["dim_lon"]) < 0:      
            tmp = np.array(DS["dim_lon"])
            tmp = np.where(tmp < 0, tmp + 360, tmp)
            DS["dim_lon"] = tmp
            # DS = DS.assign_coords({"dim_lon":("dim_lon", tmp, DS["dim_lon"].attrs)})
            DS = DS.sortby("dim_lon")
            print(f"{Red}NOTE: Longitude changed to 0-360E and all variables "
                  f"appropriately reindexed{Nclr}")

        # Add scalar axis to areo [time, scalar_axis])
        inpt_dimlist = DS.dims
        # Check if dimensions are correct or need modifying
        if "scalar_axis" not in inpt_dimlist:
            # Is scalar_axis is a dimension?
            scalar_axis = DS.assign_coords(scalar_axis = 1)
        if DS["areo"].dims != ("time", scalar_axis):
            DS["areo"] = DS["areo"].expand_dims("scalar_axis", axis=1)
            DS["areo"].attrs["long_name"] = (
                "(scalar_axis ADDED POST-PROCESSING) " 
                + DS["areo"].attrs["long_name"]
                )
            
            print(f"{Red}NOTE: scalar_axis dimension added to aero{Nclr}")

        # Pipe processed data to new **atmos_daily.nc file
        DS.to_netcdf(fullnameOUT)
        print(f"{Cyan}{fullnameOUT} was created{Nclr}")

        # Create **atmos_average.nc file
        fullnameOUT = f"{fullnameIN[:-3]}_atmos_average.nc"

        # Figure out number of timesteps there are in 5 sols
        dt_in = DS["time"][1] - DS["time"][0]
        iperday = int(np.round(1 / dt_in))
        combinedN = int(iperday * 5)
        time = "dim_time"
        
        # Coarsen the time dimension by a factor of 5 and average 
        # over each window
        DS_average = DS.coarsen(**{"dim_time": combinedN}).mean()

        # Update the time dimension coordinate values to reflect the 
        # new time axis. Save the middle value of the time interval.
        # Start indexing from 2 (3 steps in or center of 5-sol segment) 
        # and jumping in groups of 5.

        # Update the time coordinate attribute
        DS_average["time"].attrs["long_name"] = "time averaged over 5 sols"

        # Create new file
        DS_average.to_netcdf(fullnameOUT)
        print(f"{Cyan}{fullnameOUT} was created{Nclr}")

        # Create **atmos_diurn.nc file
        fullnameOUT = f"{fullnameIN[:-3]}_atmos_diurn.nc"

        # Figure out number of timesteps there are in 5 sols
        dt_in = DS["time"][1] - DS["time"][0]
        iperday = int(np.round(1 / dt_in))
        days = len(DS["time"]) / iperday

        # Create a new time_of_day_XX dimension
        tod_name = f"time_of_day_{iperday:02d}"
        print(f"{Purple}tod_name {tod_name}{Nclr}")
        DS_diurn = DS.copy()

        # Get the time of day in hours
        tod = np.mod(DS["time"][0:iperday] * 24, 24).values
        
        # Specify labels for new dimensions
        ind = pd.MultIndex.from_product((days,tod), names=("time", tod_name))
        DS_reshaped = DS_diurn.assign(time = ind).unstack("time")
        
        print(DS_reshaped)
        # Reshape methood
        DS_reshaped = xr.Dataset(
            {var_name: 
                var.isel(time=slice(0, None, iperday)).stack(time=[days, tod]) 
                for var_name, var in DS_diurn.data_vars.items() 
                if "time" in var.dims}
            )

        print(DS_reshaped["tsurf"][0, :, 0, 0].values)

        # First make a coarser version (1 time step per sol)
        DS_coarse = DS_diurn.coarsen(time = iperday).mean()
        print(DS_diurn["tsurf"][0:11, 0, 0].values)
        print(DS_coarse["tsurf"][0, 0, 0].values)
        
        # Now add a new time_of_day_XX dimension if time is a dimension 
        # for the variable
        DS_coarse.update(
            {var_name: 
                var.expand_dims(**{tod_name: iperday}, axis=1) 
                for var_name, var in DS_coarse.data_vars.items() 
                if "time" in var.dims}
            )
        
        print(DS_coarse["tsurf"][0, :, 0, 0].values)

        for i in range(iperday):
            print(i, DS_diurn["tsurf"][i, 0, 0].values)
            DS_coarse["tsurf"][0, i, :, :].values = (
                DS_diurn["tsurf"][i, :, :].values)
        
        print(DS_coarse["tsurf"][0, :, 0, 0].values)
        
        # Add a new time_of_day_XX dimension if time is a dimension of 
        # the variable
        DS_diurn.update(
            {var_name: 
                var.expand_dims(**{tod_name:iperday}, axis=1) 
                for var_name, var in DS_diurn.data_vars.items() 
                if "time" in var.dims}
            )

        # Get the time of day in hours
        tod = np.mod(DS["time"][0:iperday]*24, 24).values
        print("tod=", tod)

        # Sort the time of day e.g. if TOD = [6, 12, 18, 0] re-arrange 
        # to [0, 6, 12, 18]. Every element in array must be greater 
        # than the one to its left
        if not np.all(tod[1:] >= tod[:-1]):
            # This returns the permutation e.g. if TOD = [6, 12, 18, 0] 
            # then i_sort = [3, 0, 1, 2]
            i_sort = np.argsort(tod)
            print(i_sort)
            # Reorder TOD
            tod = tod[i_sort]
            DS_diurn.update(
                {var_name:
                    var.isel(**{tod_name:i_sort}) for var_name, var in 
                    DS_diurn.data_vars.items() if "time" in var.dims}
                )
            
        print(DS_diurn["tsurf"][0:11, 0, 0, 0])

        # Squeeze time[i:iperday] == i and time_of_day grabs the rest
        DS_coarse = DS_diurn.coarsen(time = iperday).mean()
        print(DS_coarse["tsurf"][0, :, 0, 0])

        # Get the time of day in hours
        tod = np.mod(DS["time"][0:iperday]*24, 24).values
        DS_coarse.coords[tod_name] = (tod_name, tod)

        # Add time_of_day dim if time is a dim of the variable
        DS_coarse.update(
            {var_name: 
                var.expand_dims(**{tod_name:iperday}, axis=1) 
                for var_name, var in DS_coarse.data_vars.items() 
                if "time" in var.dims}
            )

        print(DS_coarse["tsurf"][0:11, 0, 0, 0])
        for i in range(iperday):
            DS_coarse["tsurf"].isel(**{tod_name:i}, time = 0).values = (
                DS["tsurf"].isel(time = i).values)

        print(DS.tsurf[0:11, 0, 0])
        print(DS_coarse.tsurf[0, :, 0, 0])

        # Add a new time_of_day_XX dimension if time is a dimension of 
        # the variable
        DS_diurn.update(
            {var_name: 
                var.expand_dims(**{tod_name:iperday}, axis = 1) 
                for var_name, var in DS_diurn.data_vars.items() 
                if "time" in var.dims})

        print(DS_diurn.tsurf[-1, :, 24, 48])
        # Get the time of day in hours
        tod = np.mod(DS["time"][0:iperday]*24, 24).values
        print("tod=", tod)

        # Sort the time of day e.g. if TOD = [6, 12, 18, 0] re-arrange 
        # to [0, 6, 12, 18]. Every element in array must be greater 
        # than the one to its left
        if not np.all(tod[1:] >= tod[:-1]):
            # This returns the permutation e.g. if TOD = [6, 12, 18, 0] 
            # then i_sort = [3, 0, 1, 2]
            i_sort = np.argsort(tod)
            print(i_sort)
            # Reorder TOD
            tod = tod[i_sort]
            DS_diurn.update(
                {var_name:
                    var.isel(**{tod_name:i_sort}) for var_name, var 
                    in DS_diurn.data_vars.items() if "time" in var.dims})

        # Coarsen to a 5-sol resolution
        # combinedN = number of 5-sol groups in file
        combinedN = int(iperday * 5) 
        time = "dim_time"
        
        # Coarsen the time dimension by a factor of 5 and average 
        # over each window
        DS_diurnave = DS_diurn.coarsen(**{"dim_time": combinedN}).mean()

        # Update the time dimension's coordinate values to reflect the 
        # new time axis. Save the middle value of the time interval. 
        # Start indexing from 2 (3 steps in or center of 5 sol segment) 
        # and jumping in groups of 5

        # Update the time coordinate attribute
        DS_diurnave["time"].attrs["long_name"] = (
            "time averaged over 5 sols")

        print(DS_diurnave)
        # Create new file
        DS_diurnave.to_netcdf(fullnameOUT)
        print(f"{fullnameOUT} was created") 

# Create new file
      # DS_diurn.to_netcdf(fullnameOUT)

# ======================================================================
#               Bin a daily file to an average file 
# ======================================================================
def bin_average(fullnameIN, model):
    nday = 5
    fullnameOUT = f"{fullnameIN[:-3]}_to_average.nc"
    fdaily = Dataset(fullnameIN, "r", format = "NETCDF4_CLASSIC")
    var_list = filter_vars(fdaily)

    # print(fdaily.dimensions) 
    time_in = fdaily.variables["time"][:]
    Nin = len(time_in)
    dt_in = time_in[1] - time_in[0]
    iperday = int(np.round(1 / dt_in))
    combinedN = int(iperday * nday)

    N_even = Nin // combinedN
    N_left = Nin % combinedN

    if N_left != 0:
        print(f"{Yellow}***Warning*** requested {int(nday)}-sol bin period. "
              f"File has {int(iperday)} timesteps/sol and {int(Nin)}/("
              f"{int(nday)} x {int(iperday)}) is not a round number.\n    "
              f"Will use {int(N_even)} bins of ({int(nday)} x {int(iperday)})"
              f"={int(combinedN)} timesteps ({int(N_even*combinedN)}) and "
              f"discard {int(N_left)} timesteps{Nclr}")

    # Define a netcdf object from the netcdf wrapper module
    fnew = Ncdf(fullnameOUT)

    # Copy all dimensions but time from the old file to the new file 
    fnew.copy_all_dims_from_Ncfile(fdaily, exclude_dim = ["time"])

    # Calculate and log the new time array
    fnew.add_dimension("time", None)
    time_out = daily_to_average(time_in[:], dt_in, nday)
    fnew.log_axis1D("time", time_out, "time", longname_txt="sol number",
                    units_txt="days since 0000-00-00 00:00:00", cart_txt="T")

    # Loop over all variables in the file
    for ivar in var_list:
        varNcf = fdaily.variables[ivar]
        # print(ivar, varNcf.dimensions)
        if "time" in varNcf.dimensions:
            print(f"{Cyan}Processing: {ivar}...{Nclr}")
            var_out = daily_to_average(varNcf[:], dt_in, nday)
            longname_txt, units_txt = get_longname_units(fdaily, ivar)
            fnew.log_variable(ivar, var_out, varNcf.dimensions, 
                              longname_txt, units_txt)
        else:
            if ivar in ["pfull", "lat", "lon", "phalf", "pk", "bk", "pstd", 
                        "zstd", "zagl"]:
                print(f"{Cyan}Copying axis: {ivar}...{Nclr}")
                fnew.copy_Ncaxis_with_content(fdaily.variables[ivar])
            else:
                print(f"{Cyan}Copying variable: {ivar}...{Nclr}")
                fnew.copy_Ncvar(fdaily.variables[ivar])
        fnew.close()

# ======================================================================
#                           END OF PROGRAM
# ======================================================================

if __name__ == '__main__':
    main()
