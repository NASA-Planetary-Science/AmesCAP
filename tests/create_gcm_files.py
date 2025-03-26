#!/usr/bin/env python3
"""
Script to create test NetCDF files for the AMESCAP integration tests.
This script generates emars_test.nc, openmars_test.nc, pcm_test.nc, and marswrf_test.nc
with variables that exactly match the specifications in real files.
"""

import numpy as np
from netCDF4 import Dataset
import os
import datetime

def create_emars_test():
    """Create emars_test.nc with the exact variables and structure as real EMARS files."""
    nc_file = Dataset('emars_test.nc', 'w', format='NETCDF4')
    
    # Define dimensions - using exact dimensions from real EMARS files
    time_dim = nc_file.createDimension('time', 36)
    pfull_dim = nc_file.createDimension('pfull', 30)
    phalf_dim = nc_file.createDimension('phalf', 31)
    lat_dim = nc_file.createDimension('lat', 36)
    latu_dim = nc_file.createDimension('latu', 36)
    lon_dim = nc_file.createDimension('lon', 60)
    lonv_dim = nc_file.createDimension('lonv', 60)
    
    # Create helper function for creating variables
    def create_var(name, dimensions, units, longname, min_val, max_val, data_type=np.float32):
        var = nc_file.createVariable(name, data_type, dimensions)
        var.units = units
        var.long_name = longname
        
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create each variable as found in the real EMARS files
    create_var('Ls', ('time',), 'deg', 'areocentric longitude', 239.9, 269.8)
    create_var('MY', ('time',), 'Martian year', 'Mars Year', 28.0, 28.0)
    create_var('Surface_geopotential', ('lat', 'lon'), 'm^2/s/s', 'surface geopotential height', -24000.0, 26000.0)
    create_var('T', ('time', 'pfull', 'lat', 'lon'), 'K', 'Temperature', 102.4, 291.9)
    create_var('U', ('time', 'pfull', 'latu', 'lon'), 'm/s', 'zonal wind', -257.6, 402.0)
    create_var('V', ('time', 'pfull', 'lat', 'lonv'), 'm/s', 'meridional wind', -278.8, 422.0)
    create_var('ak', ('phalf',), 'pascal', 'pressure part of the hybrid coordinate', 0.0, 2.9)
    create_var('bk', ('phalf',), 'none', 'vertical coordinate sigma value', 0.0, 1.0)
    create_var('earth_day', ('time',), 'Earth day', 'Earth day of the month', 1.0, 31.0)
    create_var('earth_hour', ('time',), 'Earth hour', 'Earth hour of the day', 0.0, 23.0)
    create_var('earth_minute', ('time',), 'Earth minute', 'Earth minute of the hour', 0.0, 59.0)
    create_var('earth_month', ('time',), 'Earth month', 'Earth month of the year', 5.0, 7.0)
    create_var('earth_second', ('time',), 'Earth second', 'Earth second and fractional second of the minute', 0.0, 60.0)
    create_var('earth_year', ('time',), 'Earth year', 'Earth year AD', 2007.0, 2007.0)
    create_var('emars_sol', ('time',), 'Martian sol', 'sols after MY 22 perihelion', 3995.0, 4040.0)
    create_var('lat', ('lat',), 'degree_N', 'latitude', -89.0, 89.0)
    create_var('latu', ('latu',), 'degree_N', 'latitude', -87.4, 87.4)
    create_var('lon', ('lon',), 'degree_E', 'longitude', 3.0, 357.0)
    create_var('lonv', ('lonv',), 'degree_E', 'longitude', 0.0, 354.0)
    create_var('macda_sol', ('time',), 'Martian year', 'sols after the start of MY 24', 3143.0, 3188.0)
    create_var('mars_hour', ('time',), 'Martian hour', 'hour of the Martian day', 0.0, 23.0)
    create_var('mars_soy', ('time',), 'Martian sol', 'sols after the last Martian vernal equinox', 469.0, 514.0)
    create_var('pfull', ('pfull',), 'mb', 'ref full pressure level', 0.0, 7.7)
    create_var('phalf', ('phalf',), 'mb', 'ref half pressure level', 0.0, 7.7)
    create_var('ps', ('time', 'lat', 'lon'), 'pascal', 'surface pressure', 312.1, 1218.1)
    create_var('time', ('time',), 'Martian hour', 'number of hours since start of file', 0.0, 1103.0)
    
    nc_file.close()
    print("Created emars_test.nc")

def create_openmars_test():
    """Create openmars_test.nc with the exact variables and structure as real OpenMARS files."""
    nc_file = Dataset('openmars_test.nc', 'w', format='NETCDF4')
    
    # Define dimensions - using exact dimensions from real OpenMARS files
    time_dim = nc_file.createDimension('time', 24)
    lat_dim = nc_file.createDimension('lat', 36)
    lon_dim = nc_file.createDimension('lon', 72)
    lev_dim = nc_file.createDimension('lev', 40)
    
    # Helper function to create a variable
    def create_var(name, dimensions, units, min_val, max_val, data_type=np.float32):
        var = nc_file.createVariable(name, data_type, dimensions)
        # OpenMARS files appear to have empty units and longname
        if units:
            var.units = units
        
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create each variable as found in real OpenMARS files
    create_var('lon', ('lon',), '', -180.0, 175.0)
    create_var('lat', ('lat',), '', -87.0, 87.0)
    create_var('lev', ('lev',), '', 5.1e-05, 1.0)
    create_var('time', ('time',), '', 3181.1, 3211.0)
    create_var('Ls', ('time',), '', 264.9, 284.1)
    create_var('MY', ('time',), '', 28.0, 28.0)
    create_var('ps', ('time', 'lat', 'lon'), '', 214.5, 1133.5)
    create_var('tsurf', ('time', 'lat', 'lon'), '', 145.5, 309.9)
    create_var('co2ice', ('time', 'lat', 'lon'), '', 0.0, 6860.4)
    create_var('dustcol', ('time', 'lat', 'lon'), '', 6.8e-09, 4.5)
    create_var('u', ('time', 'lev', 'lat', 'lon'), '', -517.1, 384.8)
    create_var('v', ('time', 'lev', 'lat', 'lon'), '', -362.2, 453.3)
    create_var('temp', ('time', 'lev', 'lat', 'lon'), '', 99.3, 299.4)
    
    nc_file.close()
    print("Created openmars_test.nc")

def create_pcm_test():
    """Create pcm_test.nc with the exact variables and structure as real PCM files."""
    nc_file = Dataset('pcm_test.nc', 'w', format='NETCDF4')
    
    # Define dimensions - using exact dimensions from real PCM files
    time_dim = nc_file.createDimension('Time', 6)  # Notice capitalized 'Time'
    altitude_dim = nc_file.createDimension('altitude', 45)
    latitude_dim = nc_file.createDimension('latitude', 36)
    longitude_dim = nc_file.createDimension('longitude', 72)
    interlayer_dim = nc_file.createDimension('interlayer', 46)  # For ap/bp variables
    subsurface_dim = nc_file.createDimension('subsurface_layers', 11)
    index_dim = nc_file.createDimension('index', 10)  # For controle variable
    
    # Helper function to create a variable
    def create_var(name, dimensions, units, longname, min_val, max_val, data_type=np.float32):
        var = nc_file.createVariable(name, data_type, dimensions)
        if units:  # Some variables don't have units
            var.units = units
        if longname:  # Some variables don't have longnames
            var.long_name = longname
        
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create critical PCM variables first (subset - just the key ones)
    create_var('Ls', ('Time',), 'deg', '', 264.5, 280.4)
    create_var('Sols', ('Time',), 'sols', '', 1175.2, 1200.0)
    create_var('Time', ('Time',), 'days since 0000-00-0 00:00:00', 'Time', 488.2, 513.0)
    create_var('altitude', ('altitude',), 'km', 'pseudo-alt', 0.0, 80.5)
    create_var('latitude', ('latitude',), 'degrees_north', 'North latitude', -90.0, 90.0)
    create_var('longitude', ('longitude',), 'degrees_east', 'East longitude', -180.0, 180.0)
    
    # These are the critical variables for vertical grid in PCM
    create_var('ap', ('interlayer',), 'Pa', '', 0.0, 9.1)
    create_var('aps', ('altitude',), 'Pa', '', 0.0, 9.0)
    create_var('bp', ('interlayer',), '', '', 0.0, 1.0)
    create_var('bps', ('altitude',), '', '', 0.0, 1.0)
    
    # Add most important atmospheric variables
    create_var('temp', ('Time', 'altitude', 'latitude', 'longitude'), 'K', '', 108.7, 282.9)
    create_var('ps', ('Time', 'latitude', 'longitude'), 'Pa', '', 166.2, 1081.4)
    create_var('pressure', ('Time', 'altitude', 'latitude', 'longitude'), 'Pa', '', 0.2, 1080.9)
    create_var('u', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -226.1, 309.5)
    create_var('v', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -219.4, 230.9)
    create_var('w', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -3.1, 5.7)
    create_var('h2o_ice', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -4.3e-06, 2.9e-03)
    create_var('h2o_vap', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -8.0e-10, 1.8e-02)
    create_var('tsurf', ('Time', 'latitude', 'longitude'), 'K', '', 145.2, 314.3)
    create_var('co2ice', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', 0.0, 2189.0)
    
    # Add a few additional diagnostic variables
    create_var('tau', ('Time', 'latitude', 'longitude'), 'SI', '', 0.1, 1.7)
    create_var('controle', ('index',), '', '', 0.0, 3.4e+06)
    create_var('soildepth', ('subsurface_layers',), 'm', 'Soil mid-layer depth', 0.0, 18.5)
    
    nc_file.close()
    print("Created pcm_test.nc")

def create_marswrf_test():
    """Create marswrf_test.nc with the exact variables and structure as real MarsWRF files."""
    nc_file = Dataset('marswrf_test.nc', 'w', format='NETCDF4')
    
    # Define dimensions - using exact dimensions from real MarsWRF files
    time_dim = nc_file.createDimension('Time', 12)
    date_str_len_dim = nc_file.createDimension('DateStrLen', 19)  # Length of date strings
    bottom_top_dim = nc_file.createDimension('bottom_top', 30)
    bottom_top_stag_dim = nc_file.createDimension('bottom_top_stag', 31)
    south_north_dim = nc_file.createDimension('south_north', 36)
    south_north_stag_dim = nc_file.createDimension('south_north_stag', 37)
    west_east_dim = nc_file.createDimension('west_east', 60)
    west_east_stag_dim = nc_file.createDimension('west_east_stag', 61)
    soil_layers_stag_dim = nc_file.createDimension('soil_layers_stag', 4)
    
    # Helper function to create a variable
    def create_var(name, dimensions, units, min_val, max_val, data_type=np.float32):
        # Special handling for Times variable
        if name == 'Times':
            var = nc_file.createVariable(name, 'S1', dimensions)
            for t in range(nc_file.dimensions['Time'].size):
                date_str = f'2000-01-{t+1:02d}_00:00:00'
                for c in range(len(date_str)):
                    var[t, c] = date_str[c]
            return var
        
        var = nc_file.createVariable(name, data_type, dimensions)
        if units:  # Some variables don't have units
            var.units = units
        
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create core variables first - this is a subset of all variables in the real file
    # First, handle special string variable 'Times'
    times_var = nc_file.createVariable('Times', 'S1', ('Time', 'DateStrLen'))
    for t in range(nc_file.dimensions['Time'].size):
        date_str = f'2000-01-{t+1:02d}_00:00:00'
        for c in range(len(date_str)):
            times_var[t, c] = date_str[c].encode('utf-8')
    
    # Create the scalar constants
    cp_var = nc_file.createVariable('CP', np.float32, ())
    cp_var.units = 'J kg-1 K-1'
    cp_var[:] = 735.0
    
    g_var = nc_file.createVariable('G', np.float32, ())
    g_var.units = 'm s-2'
    g_var[:] = 3.7
    
    p0_var = nc_file.createVariable('P0', np.float32, ())
    p0_var.units = 'Pa'
    p0_var[:] = 610.0
    
    r_d_var = nc_file.createVariable('R_D', np.float32, ())
    r_d_var.units = 'J kg-1 K-1'
    r_d_var[:] = 192.0
    
    # Create critical grid and coordinate variables
    create_var('ZNU', ('Time', 'bottom_top'), '', 0.0, 1.0)
    create_var('ZNW', ('Time', 'bottom_top_stag'), '', 0.0, 1.0)
    create_var('ZS', ('Time', 'soil_layers_stag'), 'm', 0.0, 19.7)
    create_var('DZS', ('Time', 'soil_layers_stag'), 'm', 0.0, 11.2)
    create_var('U', ('Time', 'bottom_top', 'south_north', 'west_east_stag'), 'm s-1', -514.7, 519.5)
    create_var('V', ('Time', 'bottom_top', 'south_north_stag', 'west_east'), 'm s-1', -205.6, 238.3)
    create_var('W', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm s-1', -40.7, 39.0)
    create_var('PH', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm2 s-2', -73000.0, 16000.0)
    create_var('PHB', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm2 s-2', -28000.0, 270000.0)
    create_var('T', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K', -170.9, 349.1)
    create_var('T_INIT', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K', -101.3, 443.6)
    create_var('MU', ('Time', 'south_north', 'west_east'), 'Pa', -267.5, 87.1)
    create_var('MUB', ('Time', 'south_north', 'west_east'), 'Pa', 128.4, 1244.7)
    create_var('P', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', -267.2, 87.0)
    create_var('PB', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', 1.6, 1243.7)
    
    # Create important Mars-specific variables
    create_var('P_TOP', ('Time',), 'Pa', 5.0, 5.0)  # Fixed value
    create_var('T0', ('Time',), 'K', 170.0, 170.0)  # Fixed value
    create_var('T00', ('Time',), 'K', 0.0, 0.0)     # Fixed value
    create_var('P00', ('Time',), 'Pa', 0.0, 0.0)    # Fixed value
    create_var('PSFC', ('Time', 'south_north', 'west_east'), 'Pa', 86.5, 1331.8)
    create_var('HGT', ('Time', 'south_north', 'west_east'), 'm', -7400.0, 19000.0)
    create_var('TSK', ('Time', 'south_north', 'west_east'), 'K', 142.0, 308.7)
    create_var('TSLB', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), 'K', 141.9, 308.7)
    create_var('XLAT', ('Time', 'south_north', 'west_east'), 'degree_north', -89.0, 89.0)
    create_var('XLONG', ('Time', 'south_north', 'west_east'), 'degree_east', -179.0, 179.0)
    create_var('XLAT_U', ('Time', 'south_north', 'west_east_stag'), 'degree_north', -89.0, 89.0)
    create_var('XLONG_U', ('Time', 'south_north', 'west_east_stag'), 'degree_east', -178.0, 180.0)
    create_var('XLAT_V', ('Time', 'south_north_stag', 'west_east'), 'degree_north', -90.0, 90.0)
    create_var('XLONG_V', ('Time', 'south_north_stag', 'west_east'), 'degree_east', -179.0, 179.0)
    create_var('XTIME', ('Time',), 'minutes', 9.6e5, 2.4e6)
    create_var('JULIAN', ('Time',), 'days', 0.0, 668.0)
    create_var('L_S', ('Time',), 'degrees', 4.6, 360.0)
    create_var('LANDMASK', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0)  # Fixed value
    create_var('ALBEDO', ('Time', 'south_north', 'west_east'), '-', 0.1, 0.8)
    create_var('CO2ICE', ('Time', 'south_north', 'west_east'), 'kg/m^2', 0.0, 1697.5)
    create_var('SWDOWN', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 649.1)
    create_var('GLW', ('Time', 'south_north', 'west_east'), 'W m-2', 1.0, 74.0)
    create_var('HFX', ('Time', 'south_north', 'west_east'), 'W m-2', -30.0, 57.3)
    create_var('UST', ('Time', 'south_north', 'west_east'), 'm s-1', 0.0, 2.6)
    create_var('PBLH', ('Time', 'south_north', 'west_east'), 'm', 6.5, 22000.0)
    create_var('RTHRATSW', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K s-1', 0.0, 0.0)  # Fixed value
    create_var('SHADOWMASK', ('Time', 'south_north', 'west_east'), '-', 0.0, 0.0)  # Fixed value
    create_var('TAU_OD', ('Time', 'bottom_top', 'south_north', 'west_east'), 'unitless', 0.0, 0.9)
    create_var('TAU_OD2D', ('Time', 'south_north', 'west_east'), 'Unitless', 0.0, 2.1)
    
    nc_file.close()
    print("Created marswrf_test.nc")

def main():
    """Main function to create all test files."""
    create_emars_test()
    create_openmars_test()
    create_pcm_test()
    create_marswrf_test()
    
    print("All test NetCDF files created successfully.")

if __name__ == "__main__":
    main()