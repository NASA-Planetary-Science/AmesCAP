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
    
    # Define dimensions - using exact dimensions from real EMARS files (updated)
    time_dim = nc_file.createDimension('time', 1104)
    pfull_dim = nc_file.createDimension('pfull', 28)
    phalf_dim = nc_file.createDimension('phalf', 29)
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
    
    # Define dimensions - using exact dimensions from real OpenMARS files (updated)
    time_dim = nc_file.createDimension('time', 360)
    lat_dim = nc_file.createDimension('lat', 36)
    lon_dim = nc_file.createDimension('lon', 72)
    lev_dim = nc_file.createDimension('lev', 35)
    
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
    
    # Define dimensions - using exact dimensions from real PCM files (updated)
    time_dim = nc_file.createDimension('Time', 100)
    altitude_dim = nc_file.createDimension('altitude', 49)
    latitude_dim = nc_file.createDimension('latitude', 49)
    longitude_dim = nc_file.createDimension('longitude', 65)
    interlayer_dim = nc_file.createDimension('interlayer', 50)
    subsurface_dim = nc_file.createDimension('subsurface_layers', 18)
    index_dim = nc_file.createDimension('index', 100)  # Changed to 100 to match controle variable shape
    
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
    
    # Create all variables from the PCM file
    create_var('Ls', ('Time',), 'deg', '', 264.5, 280.4)
    create_var('Mccntot', ('Time', 'latitude', 'longitude'), 'kg/m2', '', 1.0e-18, 2.2e-03)
    create_var('Nccntot', ('Time', 'latitude', 'longitude'), 'Nbr/m2', '', 8.0e+01, 2.4e+11)
    create_var('Sols', ('Time',), 'sols', '', 1175.2, 1200.0)
    create_var('Time', ('Time',), 'days since 0000-00-0 00:00:00', 'Time', 488.2, 513.0)
    create_var('aire', ('latitude', 'longitude'), '', '', 6.1e+08, 7.4e+10)
    create_var('albedo', ('Time', 'latitude', 'longitude'), '', '', 0.1, 0.9)
    create_var('altitude', ('altitude',), 'km', 'pseudo-alt', 0.0, 80.5)
    create_var('ap', ('interlayer',), 'Pa', '', 0.0, 9.1)
    create_var('aps', ('altitude',), 'Pa', '', 0.0, 9.0)
    create_var('bp', ('interlayer',), '', '', 0.0, 1.0)
    create_var('bps', ('altitude',), '', '', 0.0, 1.0)
    create_var('ccnN', ('Time', 'altitude', 'latitude', 'longitude'), 'part/kg', '', -1.6e+06, 2.0e+09)
    create_var('ccnq', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -4.1e-08, 5.3e-05)
    create_var('co2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', 0.9, 1.0)
    create_var('co2ice', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', 0.0, 2189.0)
    create_var('controle', ('index',), '', '', 0.0, 3.4e+06)
    create_var('dqndust', ('Time', 'latitude', 'longitude'), 'number.m-2.s-1', '', 1.6e+01, 1.8e+05)
    create_var('dqsdust', ('Time', 'latitude', 'longitude'), 'kg.m-2.s-1', '', 1.1e-11, 9.6e-09)
    create_var('dso', ('Time', 'altitude', 'latitude', 'longitude'), 'm2.kg-1', '', 1.1e-16, 4.0e+00)
    create_var('dsodust', ('Time', 'altitude', 'latitude', 'longitude'), 'm2.kg-1', '', 1.1e-16, 4.0e+00)
    create_var('dustN', ('Time', 'altitude', 'latitude', 'longitude'), 'part/kg', '', -2.0e+05, 5.1e+09)
    create_var('dustq', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -9.0e-09, 1.2e-04)
    create_var('fluxsurf_lw', ('Time', 'latitude', 'longitude'), 'W.m-2', '', 8.4, 131.3)
    create_var('fluxsurf_sw', ('Time', 'latitude', 'longitude'), 'W.m-2', '', -0.0, 668.6)
    create_var('fluxtop_lw', ('Time', 'latitude', 'longitude'), 'W.m-2', '', 20.5, 420.7)
    create_var('fluxtop_sw', ('Time', 'latitude', 'longitude'), 'W.m-2', '', -0.0, 299.8)
    create_var('h2o_ice', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -4.3e-06, 2.9e-03)
    create_var('h2o_ice_s', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', -18.2, 8.9)
    create_var('h2o_vap', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -8.0e-10, 1.8e-02)
    create_var('hfmax_th', ('Time', 'latitude', 'longitude'), 'K.m/s', '', 0.0, 4.5)
    create_var('icetot', ('Time', 'latitude', 'longitude'), 'kg/m2', '', -6.7e-22, 1.4e-02)
    create_var('latitude', ('latitude',), 'degrees_north', 'North latitude', -90.0, 90.0)
    create_var('longitude', ('longitude',), 'degrees_east', 'East longitude', -180.0, 180.0)
    create_var('mtot', ('Time', 'latitude', 'longitude'), 'kg/m2', '', 1.9e-05, 9.0e-02)
    create_var('pdqccn2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg.s-1', '', -8.5e-06, 2.3e-05)
    create_var('pdqccnN2', ('Time', 'altitude', 'latitude', 'longitude'), 'nb/kg.s-1', '', -4.5e+08, 7.9e+08)
    create_var('pdqdust2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg.s-1', '', -2.3e-05, 8.5e-06)
    create_var('pdqdustN2', ('Time', 'altitude', 'latitude', 'longitude'), 'nb/kg.s-1', '', -7.9e+08, 4.5e+08)
    create_var('pdqice2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg.s-1', '', -5.1e-07, 2.3e-06)
    create_var('pdqvap2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg.s-1', '', -2.3e-06, 5.1e-07)
    create_var('pdtc_atm', ('Time', 'altitude', 'latitude', 'longitude'), '', '', -3.3e-03, 3.9e-03)
    create_var('phisinit', ('latitude', 'longitude'), '', '', -2.6e+04, 5.2e+04)
    create_var('pressure', ('Time', 'altitude', 'latitude', 'longitude'), 'Pa', '', 0.2, 1080.9)
    create_var('ps', ('Time', 'latitude', 'longitude'), 'Pa', '', 166.2, 1081.4)
    create_var('reffdust', ('Time', 'altitude', 'latitude', 'longitude'), 'm', '', 2.8e-08, 3.3e-05)
    create_var('reffice', ('Time', 'latitude', 'longitude'), 'm', '', 0.0e+00, 1.5e-03)
    create_var('rho', ('Time', 'altitude', 'latitude', 'longitude'), 'kg.m-3', '', 5.9e-06, 3.6e-02)
    create_var('rice', ('Time', 'altitude', 'latitude', 'longitude'), 'm', '', 1.0e-10, 5.0e-04)
    create_var('rmoym', ('Time', 'latitude', 'longitude'), 'm', '', 1.0e-30, 6.4e+01)
    create_var('saturation', ('Time', 'altitude', 'latitude', 'longitude'), 'dimless', '', -4.6e-01, 1.5e+07)
    create_var('soildepth', ('subsurface_layers',), 'm', 'Soil mid-layer depth', 0.0, 18.5)
    create_var('surfccnN', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', 2.3e+09, 2.4e+16)
    create_var('surfccnq', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', 6.8e-05, 8.1e+02)
    create_var('tau', ('Time', 'latitude', 'longitude'), 'SI', '', 0.1, 1.7)
    create_var('tauTES', ('Time', 'latitude', 'longitude'), '', '', 1.0e-19, 2.5e+00)
    create_var('tauref', ('Time', 'latitude', 'longitude'), 'NU', '', 0.0, 1.3)
    create_var('temp', ('Time', 'altitude', 'latitude', 'longitude'), 'K', '', 108.7, 282.9)
    create_var('temp7', ('Time', 'latitude', 'longitude'), 'K', '', 147.6, 273.1)
    create_var('tsurf', ('Time', 'latitude', 'longitude'), 'K', '', 145.2, 314.3)
    create_var('u', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -226.1, 309.5)
    create_var('v', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -219.4, 230.9)
    create_var('vmr_h2oice', ('Time', 'altitude', 'latitude', 'longitude'), 'mol/mol', '', -1.0e-05, 6.9e-03)
    create_var('vmr_h2ovap', ('Time', 'altitude', 'latitude', 'longitude'), 'mol/mol', '', -1.9e-09, 4.5e-02)
    create_var('w', ('Time', 'altitude', 'latitude', 'longitude'), 'm.s-1', '', -3.1, 5.7)
    create_var('wstar', ('Time', 'latitude', 'longitude'), 'm/s', '', 0.0, 8.3)
    create_var('zcondicea', ('Time', 'altitude', 'latitude', 'longitude'), '', '', -4.6e-05, 5.1e-05)
    create_var('zfallice', ('Time', 'latitude', 'longitude'), '', '', 0.0e+00, 1.2e-04)
    create_var('zmax_th', ('Time', 'latitude', 'longitude'), 'm', '', 0.0e+00, 1.2e+04)
    
    nc_file.close()
    print("Created pcm_test.nc")

def create_marswrf_test():
    """Create marswrf_test.nc with the exact variables and structure as real MarsWRF files."""
    nc_file = Dataset('marswrf_test.nc', 'w', format='NETCDF4')
    
    # Define dimensions - using exact dimensions from real MarsWRF files (updated)
    time_dim = nc_file.createDimension('Time', 100)
    date_str_len_dim = nc_file.createDimension('DateStrLen', 19)
    bottom_top_dim = nc_file.createDimension('bottom_top', 43)
    bottom_top_stag_dim = nc_file.createDimension('bottom_top_stag', 44)
    south_north_dim = nc_file.createDimension('south_north', 90)
    south_north_stag_dim = nc_file.createDimension('south_north_stag', 91)
    west_east_dim = nc_file.createDimension('west_east', 180)
    west_east_stag_dim = nc_file.createDimension('west_east_stag', 181)
    soil_layers_stag_dim = nc_file.createDimension('soil_layers_stag', 15)
    
    # Helper function to create a variable
    def create_var(name, dimensions, units, min_val, max_val, data_type=np.float32, is_coordinate=False):
        # Special handling for Times variable
        if name == 'Times':
            var = nc_file.createVariable(name, 'S1', dimensions)
            for t in range(nc_file.dimensions['Time'].size):
                date_str = f'2000-01-{(t % 31) + 1:02d}_00:00:00'
                for c in range(len(date_str)):
                    var[t, c] = date_str[c].encode('utf-8')
            return var
        
        var = nc_file.createVariable(name, data_type, dimensions)
        if units:  # Some variables don't have units
            var.units = units
        
        # Add description attribute to all variables
        var.description = f"{name} variable"
        
        # For coordinate variables, add appropriate attributes
        if is_coordinate:
            if 'LAT' in name:
                var.standard_name = 'latitude'
                var.long_name = 'LATITUDE, SOUTH IS NEGATIVE'
            elif 'LONG' in name:
                var.standard_name = 'longitude'
                var.long_name = 'LONGITUDE, WEST IS NEGATIVE'
            
            # Set the coordinates attribute which helps xarray recognize coordinate variables
            if name == 'XLAT' or name == 'XLONG':
                var.coordinates = 'XLONG XLAT'
                
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create Times variable (special handling)
    times_var = nc_file.createVariable('Times', 'S1', ('Time', 'DateStrLen'))
    for t in range(nc_file.dimensions['Time'].size):
        date_str = f'2000-01-{(t % 31) + 1:02d}_00:00:00'
        for c in range(len(date_str)):
            times_var[t, c] = date_str[c].encode('utf-8')
    
    # Create coordinate variables with is_coordinate=True flag
    create_var('XLAT', ('Time', 'south_north', 'west_east'), 'degree_north', -89.0, 89.0, is_coordinate=True)
    create_var('XLONG', ('Time', 'south_north', 'west_east'), 'degree_east', -179.0, 179.0, is_coordinate=True)
    create_var('XLAT_U', ('Time', 'south_north', 'west_east_stag'), 'degree_north', -89.0, 89.0, is_coordinate=True)
    create_var('XLONG_U', ('Time', 'south_north', 'west_east_stag'), 'degree_east', -178.0, 180.0, is_coordinate=True)
    create_var('XLAT_V', ('Time', 'south_north_stag', 'west_east'), 'degree_north', -90.0, 90.0, is_coordinate=True)
    create_var('XLONG_V', ('Time', 'south_north_stag', 'west_east'), 'degree_east', -179.0, 179.0, is_coordinate=True)
    
    # Create all variables from the MarsWRF file
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
    
    # Add all the other MarsWRF variables
    create_var('FNM', ('Time', 'bottom_top'), '', 0.0, 0.5)
    create_var('FNP', ('Time', 'bottom_top'), '', 0.0, 0.7)
    create_var('RDNW', ('Time', 'bottom_top'), '', -600.0, -40.0)
    create_var('RDN', ('Time', 'bottom_top'), '', -400.0, 0.0)
    create_var('DNW', ('Time', 'bottom_top'), '', -0.025, -0.0017)
    create_var('DN', ('Time', 'bottom_top'), '', -0.025, 0.0)
    create_var('P_HYD', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', 1.1, 1330.7)
    create_var('PSFC', ('Time', 'south_north', 'west_east'), 'Pa', 86.5, 1331.8)
    create_var('T1_5', ('Time', 'south_north', 'west_east'), 'K', 142.1, 281.2)
    create_var('TH1_5', ('Time', 'south_north', 'west_east'), 'K', 128.7, 378.1)
    create_var('U1_5', ('Time', 'south_north', 'west_east'), 'm s-1', -21.6, 21.0)
    create_var('V1_5', ('Time', 'south_north', 'west_east'), 'm s-1', -20.0, 22.6)
    create_var('U1M', ('Time', 'south_north', 'west_east'), 'm s-1', -19.2, 19.0)
    create_var('V1M', ('Time', 'south_north', 'west_east'), 'm s-1', -17.8, 20.2)
    create_var('RHOS', ('Time', 'south_north', 'west_east'), 'kg m-3', 0.00213, 0.0424)
    create_var('RDX', ('Time',), '', 8.5e-06, 8.5e-06)
    create_var('RDY', ('Time',), '', 8.5e-06, 8.5e-06)
    create_var('CF1', ('Time',), '', 1.6, 1.6)
    create_var('CF2', ('Time',), '', -0.7, -0.7)
    create_var('CF3', ('Time',), '', 0.1, 0.1)
    create_var('ITIMESTEP', ('Time',), '', 9.6e+05, 2.4e+06)
    create_var('XTIME', ('Time',), 'minutes', 9.6e+05, 2.4e+06)
    create_var('JULIAN', ('Time',), 'days', 0.0, 668.0)
    create_var('L_S', ('Time',), 'degrees', 4.6, 360.0)
    create_var('LANDMASK', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0)
    create_var('SLOPE', ('Time', 'south_north', 'west_east'), '', 6.4e-06, 1.7e-01)
    create_var('SLP_AZI', ('Time', 'south_north', 'west_east'), 'rad', 0.0, 6.3)
    create_var('TSLB', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), 'K', 141.9, 308.7)
    create_var('SOIL_DEN', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 1500.0, 1500.0)
    create_var('SOIL_CAP', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 837.0, 837.0)
    create_var('SOIL_COND', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 0.0, 0.6)
    create_var('COSZEN', ('Time', 'south_north', 'west_east'), 'dimensionless', 0.0, 1.0)
    create_var('HRANG', ('Time', 'south_north', 'west_east'), 'radians', -3.1, 3.1)
    create_var('DECLIN', ('Time',), 'radians', -0.4, 0.4)
    create_var('SOLCON', ('Time',), 'W m-2', 492.8, 716.7)
    create_var('SUNFRAC', ('Time', 'south_north', 'west_east'), '', 0.0, 1.0)
    create_var('COSZEN_MEAN', ('Time', 'south_north', 'west_east'), '', 0.0, 0.6)
    create_var('SUNBODY', ('Time',), '', 1.4, 1.7)
    create_var('KPBL', ('Time', 'south_north', 'west_east'), '', 4.0, 43.0)
    create_var('MAPFAC_MX', ('Time', 'south_north', 'west_east'), '', 1.0, 57.3)
    create_var('MAPFAC_MY', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0)
    create_var('MAPFAC_UX', ('Time', 'south_north', 'west_east_stag'), '', 1.0, 57.3)
    create_var('MAPFAC_UY', ('Time', 'south_north', 'west_east_stag'), '', 1.0, 1.0)
    create_var('MAPFAC_VX', ('Time', 'south_north_stag', 'west_east'), '', 0.0, 28.7)
    create_var('MF_VX_INV', ('Time', 'south_north_stag', 'west_east'), '', 0.0, 1.0)
    create_var('MAPFAC_VY', ('Time', 'south_north_stag', 'west_east'), '', 1.0, 1.0)
    create_var('F', ('Time', 'south_north', 'west_east'), 's-1', -1.4e-04, 1.4e-04)
    create_var('E', ('Time', 'south_north', 'west_east'), 's-1', 2.5e-06, 1.4e-04)
    create_var('TSK', ('Time', 'south_north', 'west_east'), 'K', 142.0, 308.7)
    create_var('HGT', ('Time', 'south_north', 'west_east'), 'm', -7.4e+03, 1.9e+04)
    create_var('P_FIT_M', ('Time',), 'unitless', 1.1, 1.1)
    create_var('P_TOP', ('Time',), 'Pa', 0.00568, 0.00568)
    create_var('RTHRATEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa K s-1', -14.9, 21.5)
    create_var('RTHRATLW', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K s-1', -2.9e-02, 5.3e-02)
    create_var('RTHRATSW', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K s-1', 0.0, 0.00517)
    create_var('SWDOWN', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 649.1)
    create_var('SWDOWNDIR', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 557.6)
    create_var('GSW', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 554.6)
    create_var('GLW', ('Time', 'south_north', 'west_east'), 'W m-2', 1.0, 74.0)
    create_var('HRVIS', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -0.0e+00, 2.9e-04)
    create_var('HRIR', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -2.7e-02, 3.3e-02)
    create_var('HRAERVIS', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -5.7e-06, 2.4e-03)
    create_var('HRAERIR', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -1.1e-03, 1.1e-03)
    create_var('TOASW', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 716.6)
    create_var('TOALW', ('Time', 'south_north', 'west_east'), 'W m-2', 13.5, 419.4)
    create_var('ALBEDO', ('Time', 'south_north', 'west_east'), '-', 0.1, 0.8)
    create_var('CLAT', ('Time', 'south_north', 'west_east'), 'degree_north', -89.0, 89.0)
    create_var('CLONG', ('Time', 'south_north', 'west_east'), 'degree_east', -179.0, 179.0)
    create_var('ALBBCK', ('Time', 'south_north', 'west_east'), '', 0.1, 0.5)
    create_var('EMBCK', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0)
    create_var('THCBCK', ('Time', 'south_north', 'west_east'), '', 30.0, 877.2)
    create_var('EMISS', ('Time', 'south_north', 'west_east'), '', 0.5, 1.0)
    create_var('RUBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa m s-2', -24.2, 24.4)
    create_var('RVBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa m s-2', -24.6, 25.5)
    create_var('RTHBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa K s-1', -62.3, 115.4)
    create_var('XLAND', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0)
    create_var('ZNT', ('Time', 'south_north', 'west_east'), 'm', 0.0, 0.1)
    create_var('UST', ('Time', 'south_north', 'west_east'), 'm s-1', 0.0, 2.6)
    create_var('PBLH', ('Time', 'south_north', 'west_east'), 'm', 6.5, 22000.0)
    create_var('THC', ('Time', 'south_north', 'west_east'), 'J m-1 K-1 s-0.5', 30.0, 877.2)
    create_var('HFX', ('Time', 'south_north', 'west_east'), 'W m-2', -30.0, 57.3)
    create_var('RNET_2D', ('Time', 'south_north', 'west_east'), 'W m-2', -128.4, 289.9)
    create_var('FLHC', ('Time', 'south_north', 'west_east'), '', 0.0, 2.0)
    create_var('ANGSLOPE', ('Time', 'south_north', 'west_east'), 'radians', 6.4e-06, 1.7e-01)
    create_var('AZMSLOPE', ('Time', 'south_north', 'west_east'), 'radians', 0.0, 6.3)
    create_var('CO2ICE', ('Time', 'south_north', 'west_east'), 'kg/m^2', 0.0, 1697.5)
    create_var('CDOD_SCALE', ('Time', 'south_north', 'west_east'), 'Unitless', 1.0, 1.0)
    create_var('TAU_OD2D', ('Time', 'south_north', 'west_east'), 'Unitless', 0.0, 2.1)
    create_var('FRAC_PERM_CO2', ('Time', 'south_north', 'west_east'), 'fraction', 0.0, 1.0)
    create_var('FRAC_PERM_H2O', ('Time', 'south_north', 'west_east'), 'fraction', 0.0, 1.0)
    create_var('GRD_ICE_PC', ('Time', 'south_north', 'west_east'), 'percent', 0.0, 1.0)
    create_var('GRD_ICE_DP', ('Time', 'south_north', 'west_east'), 'meters', -9999.0, 2.4)
    create_var('TAU_OD', ('Time', 'bottom_top', 'south_north', 'west_east'), 'unitless', 0.0, 0.9)
    
    # Add global scalar constants
    nc_file.setncattr('P0', 610.0)  # Reference pressure in Pa
    nc_file.setncattr('G', 3.72)    # Gravity on Mars in m/s²
    nc_file.setncattr('CP', 770.0)  # Specific heat capacity
    nc_file.setncattr('R_D', 192.0) # Gas constant for Mars atmosphere
    nc_file.setncattr('T0', 300.0)  # Reference temperature in K
    # cp_var = nc_file.createVariable('CP', np.float32, ())
    # cp_var.units = 'J kg-1 K-1'
    # cp_var[:] = 735.0
    
    # g_var = nc_file.createVariable('G', np.float32, ())
    # g_var.units = 'm s-2'
    # g_var[:] = 3.7
    
    # p0_var = nc_file.createVariable('P0', np.float32, ())
    # p0_var.units = 'Pa'
    # p0_var[:] = 610.0
    
    # r_d_var = nc_file.createVariable('R_D', np.float32, ())
    # r_d_var.units = 'J kg-1 K-1'
    # r_d_var[:] = 192.0
    
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