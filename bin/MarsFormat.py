#!/usr/bin/env python3
"""
The MarsFormat executable is a routine that transforms non-MGCM model
output into MGCM-like model output for compatibility with CAP. 

MarsFormat changes variable names, dimension names, dimension order,
and units to the configuration expected by CAP. In some cases, such as
for MarsWRF, variables are derived and regridded onto a standard grid.

The executable requires x arguments:
    * [-openmars --openmars]    convert openMars data to MGCM format
    * [-marswrf --marswrf]      convert MarsWRF data to MGCM format
    

Third-party Requirements:
    * numpy
    * argparse
    * xarray

List of Functions:
    * x
"""

# make print statements appear in color
from amescap.Script_utils import prCyan, Yellow, NoColor, Green

# load generic Python modules
import argparse      # parse arguments
import numpy as np
import xarray as xr
import os            # access operating system functions

# load amesCAP modules
from amescap.FV3_utils import layers_mid_point_to_boundary

xr.set_options(keep_attrs=True)


# ======================================================
#                  ARGUMENT PARSER
# ======================================================

parser = argparse.ArgumentParser(
    description=(
        f"{Yellow}MarsFormat is for converting non-MGCM output "
        f"to MGCM format.{NoColor}\n\n"
    ),
    formatter_class = argparse.RawTextHelpFormatter
)

parser.add_argument(
    'input_file', nargs='+',
    help=(
        f"A netCDF file or list of netCDF files.\n\n"
    )
)

parser.add_argument(
    '-openmars', '--openmars', nargs='+',
    help=(
        f"Produce an MGCM-like daily file from an {Yellow}openMars"
        f"{NoColor} file.\n"
        f"{Green}Usage:\n"
        f"> MarsFormat.py input_file*.nc -openmars daily"
        f"{NoColor}\n\n"
    )
)

parser.add_argument(
    '-marswrf', '--marswrf', nargs='+',
    help=(
        f"Produce an MGCM-like daily file from a "
        f"{Yellow}MarsWRF{NoColor} file.\n"
        f"{Green}Usage:\n"
        f"> MarsFormat.py input_file*.nc -marswrf daily"
        f"{NoColor}\n\n"
    )
)

# parser.add_argument(
#     '-legacy', '--legacy', nargs='+', help=argparse.SUPPRESS
# )

parser.add_argument('--debug', action='store_true',
    help = (f"Debug flag: release the exceptions.\n\n")
)


# ======================================================
#                  DEFINITIONS
# ======================================================

# MarsWRF
def marswrf_to_mgcm(DS):
    #TODO longname is 'description' for MarsWRF

    # Find shape of coordinates. Expecting [t,z,y,x]
    WRF_dims = np.shape(DS.T)
    Xmax = WRF_dims[3]  # x
    Ymax = WRF_dims[2]  # y
    Tmax = WRF_dims[0]  # t
    Zmax = WRF_dims[1]  # z (layer)

    # Define coordinates for new dataFrame
    time = DS.XTIME/60/24  # minutes since simulation start [m]
    lat = DS.XLAT[0, :, 0]
    lon2D = DS.XLONG[0, :]
    lon = np.squeeze(lon2D[0, :])

    # Derive half and full reference pressure levels (phalf and
    # pfull in MGCM; Pa)
    pfull = DS.P_TOP[0] + DS.ZNU[0, :]*DS.P0
    phalf = DS.P_TOP[0] + DS.ZNW[0, :]*DS.P0

    # Calculate level height above the Surface (i.e. above topo)
    zagl_lvl = (
        (DS.PH[:, :Zmax, :, :] + DS.PHB[0, :Zmax, :, :])
        / DS.G
        - DS.HGT[0, :, :]
    )

    # Find layer pressures [Pa]
    try:
        # If perturbation pressure (P_TOP) and base state
        # pressure are present, use them
        pfull3D = DS.P_TOP + DS.PB[0, :]
    except NameError:
        # otherwise, multiply surface pressure by the sigma
        # value (p/psurf) of each layer
        pfull3D = DS.PSFC[:, :Ymax, :Xmax] * DS.ZNU[:, :Zmax]

    # Interpolate variables onto regular mass grid [t,z,y,x]

    # For variables staggered in X (lon)
    #       [t,z,y,x'] -> [t,z,y,x]
    ucomp = 0.5*(DS.U[..., :-1] + DS.U[..., 1:])

    # For variables staggered Y (lat)
    #       [t,z,y',x] -> [t,z,y,x]:
    vcomp = 0.5*(DS.V[:, :, :-1, :] + DS.V[:, :, 1:, :])

    # For variables staggered Z/P (height)
    #       [t,z',y,x] -> [t,z,y,x]:
    w = 0.5*(DS.W[:, :-1, :, :] + DS.W[:, 1:, :, :])

    # Interpolate to find the layer heights above the surface
    # (i.e., above topography; m)
    zfull3D = 0.5*(zagl_lvl[:, :-1, :, :] + zagl_lvl[:, 1:, :, :])

    # Derive atmospheric temperature [K]
    gamma = (DS.CP/(DS.CP - DS.R_D))
    temp = (DS.T + DS.T0)*(pfull3D / DS.P0)**((gamma-1.)/gamma)

    # Derive ak, bk
    ak = np.zeros(len(phalf))
    bk = np.zeros(len(phalf))
    ak[-1] = DS.P_TOP[0]  # in MarsWRF, pressure increases w/N
    bk[:] = DS.ZNW[0, :]

    # Archive variables
    # Each entry has [name, values, dimensions, longname,units]
    var_dict = {
        'ak': [
            ak, ['phalf'],
            'pressure part of the hybrid coordinate', 'Pa'],
        'bk': [
            bk, ['phalf'],
            'vertical coordinate sigma value', 'none'],
        'areo': [
            DS.L_S, ['time'],
            'solar longitude', 'degree'],
        'ps': [
            DS.PSFC, ['time', 'lat', 'lon'],
            'surface pressure', 'Pa'],
        'zsurf': [
            DS.HGT[0, :], ['lat', 'lon'],
            'surface height', 'm'],
        'ucomp': [
            ucomp, ['time', 'pfull', 'lat', 'lon'],
            'zonal winds', 'm/sec'],
        'vcomp': [
            vcomp, ['time', 'pfull', 'lat', 'lon'],
            'meridional winds', 'm/sec'],
        'w': [
            w, ['time', 'pfull', 'lat', 'lon'],
            'vertical winds', 'm/s'],
        'pfull3D': [
            pfull3D, ['time', 'pfull', 'lat', 'lon'],
            'pressure', 'Pa'],
        'temp': [
            temp, ['time', 'pfull', 'lat', 'lon'],
            'temperature', 'K'],
        'h2o_ice_sfc': [
            DS.H2OICE, ['time', 'lat', 'lon'],
            'surface H2O Ice', 'kg/m2'],
        'co2_ice_sfc': [
            DS.CO2ICE, ['time', 'lat', 'lon'],
            'surface CO2 Ice', 'kg/m2'],
        'ts': [
            DS.TSK, ['time', 'lat', 'lon'],
            'surface temperature', 'K'],
    }
    return var_dict, time, lat, lon, phalf, pfull


# OpenMars
def openmars_to_mgcm(DS):
    # Define coordinates for new DataFrame
    ref_press = 720  # TODO this is added on to create ak/bk
    time = DS.time  # minutes since simulation start [m]
    lat = DS.lat
    lon = DS.lon

    # Derive half and full reference pressure levels (phalf and
    # pfull; Pa)
    pfull = DS.lev*ref_press

    # Add p_half dimensions and ak, bk vertical grid coordinates
    # DS.expand_dims({'p_half':len(pfull)+1})

    """
    Compute sigma values.
    Swap the sigma array upside down twice with [::-1].
    The first time is for layers_mid_point_to_boundary() which 
    needs sigma[0]=0, sigma[-1]=1. The second time is to 
    reorganize the array in the original openMars format where 
    sigma[0]=1, sigma[-1]=0
    """
    DS['bk'] = layers_mid_point_to_boundary(DS.lev[::-1], 1.)[::-1]
    # Pure sigma model, set bk=0
    DS['ak'] = np.zeros(len(pfull)+1)
    phalf = np.array(DS['ak']) + ref_press*np.array(DS['bk'])

    # Archive variables
    # Each entry has [name, values, dimensions, longname,units]
    var_dict = {
        'bk': [
            DS.bk, ['phalf'],
            'vertical coordinate sigma value', 'none'],
        'ak': [
            DS.ak, ['phalf'],
            'pressure part of the hybrid coordinate', 'Pa'],
        'areo': [
            DS.Ls, ['time'],
            'solar longitude', 'degree'],
        'ps': [
            DS.ps, ['time', 'lat', 'lon'],
            'surface pressure', 'Pa'],
        'ucomp': [
            DS.u, ['time', 'pfull', 'lat', 'lon'],
            'zonal winds', 'm/sec'],
        'vcomp': [
            DS.v, ['time', 'pfull', 'lat', 'lon'],
            'meridional wind', 'm/sec'],
        'temp': [
            DS.temp, ['time', 'pfull', 'lat', 'lon'],
            'temperature', 'K'],
        'dust_mass_col': [
            DS.dustcol, ['time', 'lat', 'lon'],
            'column integration of dust', 'kg/m2'],
        'co2_ice_sfc': [
            DS.co2ice, ['time', 'lat', 'lon'],
            'surace CO2 ice', 'kg/m2'],
        'ts': [
            DS.tsurf, ['time', 'lat', 'lon'],
            'surface temperature', 'K']
    }
    return var_dict, time, lat, lon, phalf, pfull
    
    
# ======================================================
#                  MAIN PROGRAM
# ======================================================

def main():
    file_list = parser.parse_args().input_file
    data_dir = os.getcwd()
    
    for file_name in file_list:
        # if full path is not already in file_name, add it
        if not ('/' in file_name):
            input_file_name = f"{data_dir}/{file_name}"
        else:
            input_file_name = file_name
        output_file_name = f"{input_file_name[:-3]}_atmos_daily.nc"

        print('Processing...')
        
        input_DS = xr.open_dataset(input_file_name, decode_times=False)

        # If user indicated marsWRF data
        if parser.parse_args().marswrf:
            (archive_vars, time, lat, 
            lon, phalf, pfull) = marswrf_to_mgcm(input_DS)

        # if user indicated openMars data
        elif parser.parse_args().openmars:
            (archive_vars, time, lat, 
            lon, phalf, pfull) = openmars_to_mgcm(input_DS)

        # Create output file (common to all models)
        coord_attributes = {
            'time': ['time', 'days'],
            'pfull': ['ref full pressure level', 'Pa'],
            'lat': ['latitudes', 'degrees_N'],
            'lon': ['longitudes', 'degrees_E'],
            'phalf': ['ref pressure at layer boundaries', 'Pa']
        }

        coord_values = {
            'time': np.array(time),
            'phalf': np.array(phalf),
            'pfull': np.array(pfull),
            'lat': np.array(lat),
            'lon': np.array(lon)
        }
        
        # Empty the dictionary
        var_DataArray = {}
        # Assign description and unit attributes to the dictionary
        for var in archive_vars.keys():
            var_DataArray[var] = xr.DataArray(
                np.array(archive_vars[var][0]), 
                dims = archive_vars[var][1]
            )
            var_DataArray[var].attrs['long_name'] = archive_vars[var][2]
            var_DataArray[var].attrs['units'] =  archive_vars[var][3]

        # Create new DataFrame
        DF = xr.Dataset(var_DataArray, coords=coord_values)

        # Add longname and unit attibutes to the coordiate variables
        for var in coord_attributes.keys():
            DF[var].attrs['long_name'] = coord_attributes[var][0]
            DF[var].attrs['units'] = coord_attributes[var][1]

        # Check that vertical grid starts at top of atmosphere (TOA)
        # and highest pressures are at the surface.
        # If pressure at TOA is not minimized, flip array
        if DF.pfull[0] != DF.pfull.min(): 
            # regrid DS based on pfull
            DF = DF.isel(pfull = slice(None, None, -1))
            # flip phalf, ak, bk
            DF = DF.isel(phalf = slice(None, None, -1))

        # Reorder dimensions
        DF = DF.transpose("phalf", "time", "pfull", "lat", "lon")

        # Convert longitude from -180-179 -> 0-360
        if min(DF.lon)<0:
            tmp_lon = np.array(DF.lon)
            tmp_lon = np.where(tmp_lon<0, tmp_lon+360, tmp_lon)
            DF = DF.assign_coords(
                {
                    'lon': ('lon', tmp_lon, DF.lon.attrs)
                }
            )
            DF = DF.sortby("lon")

        # If areo does not have dimension 'scalar_axis', add it
        input_dimensions = DF.dims
        if 'scalar_axis' not in input_dimensions:
            scalar_axis = DF.assign_coords(scalar_axis=1)
        
        # Ensure areo dimension order is ('time', 'scalar_axis')
        if DF.areo.dims != ('time', scalar_axis):
            DF['areo'] = DF.areo.expand_dims('scalar_axis', axis=1)

        # Pipe processed data to a new netCDF file
        DF.to_netcdf(output_file_name)
        
        prCyan(f"{output_file_name} was created")

# ======================================================
#                  END OF PROGRAM
# ======================================================

if __name__ == "__main__":
    main()
