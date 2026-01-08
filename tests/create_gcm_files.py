#!/usr/bin/env python3
"""
Script to create test NetCDF files for the AMESCAP integration tests.
This script generates emars_test.nc, openmars_test.nc, pcm_test.nc, and marswrf_test.nc
with variables that exactly match the specifications in real files.
"""

import numpy as np
from netCDF4 import Dataset

# ----------------------------------------------------------------------
#                      EMARS Dummy File
# ----------------------------------------------------------------------

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
    
    # Create variables with longname and units
    lonv_var = nc_file.createVariable('lonv', 'f4', ('lonv',))
    lonv_var.long_name = "longitude"
    lonv_var.units = "degree_E"
    
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lat_var.long_name = "latitude"
    lat_var.units = "degree_N"
    
    ak_var = nc_file.createVariable('ak', 'f4', ('phalf',))
    ak_var.long_name = "pressure part of the hybrid coordinate"
    ak_var.units = "pascal"
    
    bk_var = nc_file.createVariable('bk', 'f4', ('phalf',))
    bk_var.long_name = "vertical coordinate sigma value"
    bk_var.units = "none"
    
    phalf_var = nc_file.createVariable('phalf', 'f4', ('phalf',))
    phalf_var.long_name = "ref half pressure level"
    phalf_var.units = "mb"
    
    latu_var = nc_file.createVariable('latu', 'f4', ('latu',))
    latu_var.long_name = "latitude"
    latu_var.units = "degree_N"
    
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lon_var.long_name = "longitude"
    lon_var.units = "degree_E"
    
    Ls_var = nc_file.createVariable('Ls', 'f4', ('time',))
    Ls_var.long_name = "areocentric longitude"
    Ls_var.units = "deg"
    
    MY_var = nc_file.createVariable('MY', 'f4', ('time',))
    MY_var.long_name = "Mars Year"
    MY_var.units = "Martian year"
    
    earth_year_var = nc_file.createVariable('earth_year', 'f4', ('time',))
    earth_year_var.long_name = "Earth year AD"
    earth_year_var.units = "Earth year"
    
    earth_month_var = nc_file.createVariable('earth_month', 'f4', ('time',))
    earth_month_var.long_name = "Earth month of the year"
    earth_month_var.units = "Earth month"
    
    earth_day_var = nc_file.createVariable('earth_day', 'f4', ('time',))
    earth_day_var.long_name = "Earth day of the month"
    earth_day_var.units = "Earth day"
    
    earth_hour_var = nc_file.createVariable('earth_hour', 'f4', ('time',))
    earth_hour_var.long_name = "Earth hour of the day"
    earth_hour_var.units = "Earth hour"
    
    earth_minute_var = nc_file.createVariable('earth_minute', 'f4', ('time',))
    earth_minute_var.long_name = "Earth minute of the hour"
    earth_minute_var.units = "Earth minute"
    
    earth_second_var = nc_file.createVariable('earth_second', 'f4', ('time',))
    earth_second_var.long_name = "Earth second and fractional second of the minute"
    earth_second_var.units = "Earth second"
    
    emars_sol_var = nc_file.createVariable('emars_sol', 'f4', ('time',))
    emars_sol_var.long_name = "sols after MY 22 perihelion"
    emars_sol_var.units = "Martian sol"
    
    macda_sol_var = nc_file.createVariable('macda_sol', 'f4', ('time',))
    macda_sol_var.long_name = "sols after the start of MY 24"
    macda_sol_var.units = "Martian year"
    
    mars_hour_var = nc_file.createVariable('mars_hour', 'f4', ('time',))
    mars_hour_var.long_name = "hour of the Martian day"
    mars_hour_var.units = "Martian hour"
    
    mars_soy_var = nc_file.createVariable('mars_soy', 'f4', ('time',))
    mars_soy_var.long_name = "sols after the last Martian vernal equinox"
    mars_soy_var.units = "Martian sol"
    
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    time_var.long_name = "number of hours since start of file"
    time_var.units = "Martian hour"
    
    pfull_var = nc_file.createVariable('pfull', 'f4', ('pfull',))
    pfull_var.long_name = "ref full pressure level"
    pfull_var.units = "mb"

    # --------- Generate realistic values for EMARS-like data ----------
    # earth_month: Values represent month numbers (5=May, 6=June, 7=July)
    def generate_earth_months(length):
        months = []
        
        # Based on the data, first ~339 entries are month 5 (May)
        may_entries = 337  # Exact number from data inspection
        
        # Month 6 (June) entries
        june_entries = 551  # Exact number from data inspection
        
        # Month 7 (July) entries for the rest
        july_entries = length - may_entries - june_entries
        
        # Create the month array
        months.extend([5] * may_entries)
        months.extend([6] * june_entries)
        months.extend([7] * july_entries)
        
        return np.array(months[:length])

    # earth_second: Decreases by precisely 21.0321 seconds each step
    def generate_earth_seconds(length):
        seconds = []
        current = 9.2703  # Starting value
        decrement = 21.0321  # Exact decrement between values
        
        for _ in range(length):
            seconds.append(current)
            current -= decrement
            if current < 0:
                current += 60  # Wrap around when going below 0
        
        return np.array(seconds)

    # earth_day, earth_hour, earth_minute are all unusual patterns.
    # define explicitly
    earth_day = [
        17., 17., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18.,
        18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 19., 19., 19., 19., 19., 
        19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 
        19., 19., 19., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 21., 21., 21., 21., 21., 21., 
        21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 
        21., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 
        22., 22., 22., 22., 22., 22., 22., 22., 22., 23., 23., 23., 23., 23., 23., 23., 
        23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 
        24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 
        24., 24., 24., 24., 24., 24., 24., 25., 25., 25., 25., 25., 25., 25., 25., 25., 
        25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 26., 
        26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 
        26., 26., 26., 26., 26., 26., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 
        27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 28., 28., 28., 
        28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 
        28., 28., 28., 28., 28., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 
        29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 30., 30., 30., 30., 
        30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 
        30., 30., 30., 31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 
        31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 31., 1., 1., 1., 1., 1., 1., 
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 2., 2., 2., 
        2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 
        2., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 
        3., 3., 3., 3., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 4., 
        4., 4., 4., 4., 4., 4., 4., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 
        5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 5., 6., 6., 6., 6., 6., 6., 6., 6., 6., 
        6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 6., 7., 7., 7., 7., 7., 7., 
        7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 7., 8., 8., 8., 
        8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 8., 
        8., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 9., 
        9., 9., 9., 9., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 
        10., 10., 10., 10., 10., 10., 10., 10., 10., 10., 11., 11., 11., 11., 11., 11., 
        11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 
        11., 11., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12., 
        12., 12., 12., 12., 12., 12., 12., 12., 12., 13., 13., 13., 13., 13., 13., 13., 
        13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 13., 
        14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 14., 
        14., 14., 14., 14., 14., 14., 14., 14., 15., 15., 15., 15., 15., 15., 15., 15., 
        15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 15., 16., 
        16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 
        16., 16., 16., 16., 16., 16., 16., 17., 17., 17., 17., 17., 17., 17., 17., 17., 
        17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 18., 18., 
        18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18., 
        18., 18., 18., 18., 18., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 
        19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 19., 20., 20., 20., 
        20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 20., 
        20., 20., 20., 20., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 
        21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 21., 22., 22., 22., 22., 22., 
        22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 
        22., 22., 22., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 
        23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 24., 24., 24., 24., 24., 24., 
        24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 24., 
        24., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 25., 
        25., 25., 25., 25., 25., 25., 25., 25., 25., 26., 26., 26., 26., 26., 26., 26., 
        26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 26., 
        27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 27., 
        27., 27., 27., 27., 27., 27., 27., 28., 28., 28., 28., 28., 28., 28., 28., 28., 
        28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 28., 29., 
        29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 29., 
        29., 29., 29., 29., 29., 29., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 
        30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 30., 1., 1., 
        1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 
        1., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 
        2., 2., 2., 2., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 3., 
        3., 3., 3., 3., 3., 3., 3., 3., 4., 4., 4., 4.
    ]

    earth_hour = [
        22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 13., 14., 15.,
        16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 
        9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 1., 
        2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 
        19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 
        14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 
        7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 
        23., 0., 1., 2., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 
        17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 
        11., 12., 13., 14., 15., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 
        5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 
        22., 23., 0., 1., 2., 3., 4., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 
        16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 
        10., 11., 12., 13., 14., 15., 16., 17., 18., 20., 21., 22., 23., 0., 1., 2., 
        3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 
        20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 9., 10., 11., 12., 13., 
        14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 
        7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 23., 
        0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 
        18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 12., 
        13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 
        6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 
        22., 23., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 
        16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 
        10., 11., 12., 13., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 
        3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 
        20., 21., 22., 23., 0., 1., 2., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 
        14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 
        7., 8., 9., 10., 11., 12., 13., 14., 15., 17., 18., 19., 20., 21., 22., 23., 
        0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 
        18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 7., 8., 9., 10., 11., 
        12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 
        5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 20., 21., 22., 
        23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 
        17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 
        11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 
        3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 
        20., 21., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 
        15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 
        8., 9., 10., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 
        1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 
        18., 19., 20., 21., 22., 23., 0., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 
        12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 
        5., 6., 7., 8., 9., 10., 11., 12., 13., 15., 16., 17., 18., 19., 20., 21., 22., 
        23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 
        17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 4., 5., 6., 7., 8., 9., 10., 
        11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 
        3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 18., 19., 20., 
        21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 
        15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 7., 8., 
        9., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 
        1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 17., 
        18., 19., 21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 
        13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 5., 
        6., 7., 8., 10., 11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 
        23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 
        17., 18., 19., 20., 21., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 
        12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3., 4., 
        5., 6., 7., 8., 9., 10., 11., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 
        23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 15., 16., 
        17., 18., 19., 20., 21., 22., 23., 0., 2., 3., 4., 5., 6., 7., 8., 9., 10., 
        11., 12., 13., 14., 15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 
        3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 15., 16., 17., 18., 19., 20., 
        21., 22., 23., 0., 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14., 
        15., 16., 17., 18., 19., 20., 21., 22., 23., 0., 1., 2., 3.
    ]

    earth_minute = [
        40., 41., 43., 45., 46., 48., 50., 51., 53., 54., 56., 58., 59., 1., 3.,
        4., 6., 8., 9., 11., 13., 14., 16., 18., 19., 21., 23., 24., 26., 27., 29., 
        31., 32., 34., 36., 37., 39., 41., 42., 44., 46., 47., 49., 51., 52., 54., 56., 
        57., 59., 0., 2., 4., 5., 7., 9., 10., 12., 14., 15., 17., 19., 20., 22., 24., 
        25., 27., 29., 30., 32., 33., 35., 37., 38., 40., 42., 43., 45., 47., 48., 50., 
        52., 53., 55., 57., 58., 0., 2., 3., 5., 6., 8., 10., 11., 13., 15., 16., 18., 
        20., 21., 23., 25., 26., 28., 30., 31., 33., 34., 36., 38., 39., 41., 43., 44., 
        46., 48., 49., 51., 53., 54., 56., 58., 59., 1., 3., 4., 6., 7., 9., 11., 12., 
        14., 16., 17., 19., 21., 22., 24., 26., 27., 29., 31., 32., 34., 36., 37., 39., 
        40., 42., 44., 45., 47., 49., 50., 52., 54., 55., 57., 59., 0., 2., 4., 5., 7., 
        9., 10., 12., 13., 15., 17., 18., 20., 22., 23., 25., 27., 28., 30., 32., 33., 
        35., 37., 38., 40., 42., 43., 45., 46., 48., 50., 51., 53., 55., 56., 58., 0., 
        1., 3., 5., 6., 8., 10., 11., 13., 14., 16., 18., 19., 21., 23., 24., 26., 28., 
        29., 31., 33., 34., 36., 38., 39., 41., 43., 44., 46., 47., 49., 51., 52., 54., 
        56., 57., 59., 1., 2., 4., 6., 7., 9., 11., 12., 14., 16., 17., 19., 20., 22., 
        24., 25., 27., 29., 30., 32., 34., 35., 37., 39., 40., 42., 44., 45., 47., 49., 
        50., 52., 53., 55., 57., 58., 0., 2., 3., 5., 7., 8., 10., 12., 13., 15., 17., 
        18., 20., 22., 23., 25., 26., 28., 30., 31., 33., 35., 36., 38., 40., 41., 43., 
        45., 46., 48., 50., 51., 53., 54., 56., 58., 59., 1., 3., 4., 6., 8., 9., 11., 
        13., 14., 16., 18., 19., 21., 23., 24., 26., 27., 29., 31., 32., 34., 36., 37., 
        39., 41., 42., 44., 46., 47., 49., 51., 52., 54., 56., 57., 59., 0., 2., 4., 
        5., 7., 9., 10., 12., 14., 15., 17., 19., 20., 22., 24., 25., 27., 29., 30., 
        32., 33., 35., 37., 38., 40., 42., 43., 45., 47., 48., 50., 52., 53., 55., 57., 
        58., 0., 2., 3., 5., 6., 8., 10., 11., 13., 15., 16., 18., 20., 21., 23., 25., 
        26., 28., 30., 31., 33., 34., 36., 38., 39., 41., 43., 44., 46., 48., 49., 51., 
        53., 54., 56., 58., 59., 1., 3., 4., 6., 7., 9., 11., 12., 14., 16., 17., 19., 
        21., 22., 24., 26., 27., 29., 31., 32., 34., 36., 37., 39., 40., 42., 44., 45., 
        47., 49., 50., 52., 54., 55., 57., 59., 0., 2., 4., 5., 7., 9., 10., 12., 13., 
        15., 17., 18., 20., 22., 23., 25., 27., 28., 30., 32., 33., 35., 37., 38., 40., 
        42., 43., 45., 46., 48., 50., 51., 53., 55., 56., 58., 0., 1., 3., 5., 6., 8., 
        10., 11., 13., 14., 16., 18., 19., 21., 23., 24., 26., 28., 29., 31., 33., 34., 
        36., 38., 39., 41., 43., 44., 46., 47., 49., 51., 52., 54., 56., 57., 59., 1., 
        2., 4., 6., 7., 9., 11., 12., 14., 16., 17., 19., 20., 22., 24., 25., 27., 29., 
        30., 32., 34., 35., 37., 39., 40., 42., 44., 45., 47., 49., 50., 52., 53., 55., 
        57., 58., 0., 2., 3., 5., 7., 8., 10., 12., 13., 15., 17., 18., 20., 22., 23., 
        25., 26., 28., 30., 31., 33., 35., 36., 38., 40., 41., 43., 45., 46., 48., 50., 
        51., 53., 54., 56., 58., 59., 1., 3., 4., 6., 8., 9., 11., 13., 14., 16., 18., 
        19., 21., 23., 24., 26., 27., 29., 31., 32., 34., 36., 37., 39., 41., 42., 44., 
        46., 47., 49., 51., 52., 54., 56., 57., 59., 0., 2., 4., 5., 7., 9., 10., 12., 
        14., 15., 17., 19., 20., 22., 24., 25., 27., 29., 30., 32., 33., 35., 37., 38., 
        40., 42., 43., 45., 47., 48., 50., 52., 53., 55., 57., 58., 0., 2., 3., 5., 6., 
        8., 10., 11., 13., 15., 16., 18., 20., 21., 23., 25., 26., 28., 30., 31., 33., 
        34., 36., 38., 39., 41., 43., 44., 46., 48., 49., 51., 53., 54., 56., 58., 59., 
        1., 3., 4., 6., 7., 9., 11., 12., 14., 16., 17., 19., 21., 22., 24., 26., 27., 
        29., 31., 32., 34., 36., 37., 39., 40., 42., 44., 45., 47., 49., 50., 52., 54., 
        55., 57., 59., 0., 2., 4., 5., 7., 9., 10., 12., 13., 15., 17., 18., 20., 22., 
        23., 25., 27., 28., 30., 32., 33., 35., 37., 38., 40., 42., 43., 45., 46., 48., 
        50., 51., 53., 55., 56., 58., 0., 1., 3., 5., 6., 8., 10., 11., 13., 14., 16., 
        18., 19., 21., 23., 24., 26., 28., 29., 31., 33., 34., 36., 38., 39., 41., 43., 
        44., 46., 47., 49., 51., 52., 54., 56., 57., 59.
    ]
    # emars_sol: Values from 3995 to 4040 with each value repeating 24 times
    def generate_emars_sol(length):
        sol_values = []
        
        # Each Mars sol (day) has 24 hours
        for sol in range(3995, 4041):  # From 3995 to 4040
            sol_values.extend([sol] * 24)  # Repeat each sol 24 times (for each hour)
        
        return np.array(sol_values[:length])

    # macda_sol: Values from 3143 to 3188 with each value repeating 24 times
    def generate_macda_sol(length):
        sol_values = []
        
        # Each Mars sol (day) has 24 hours
        for sol in range(3143, 3189):  # From 3143 to 3188
            sol_values.extend([sol] * 24)  # Repeat each sol 24 times (for each hour)
        
        return np.array(sol_values[:length])

    # mars_hour: Complete 24-hour cycle (0-23) repeated throughout the dataset
    def generate_mars_hours(length):
        # Simple repeating pattern of hours 0-23
        hours = list(range(24))  # Hours 0 through 23
        
        # Repeat the pattern as needed
        full_cycles = length // 24
        remainder = length % 24
        
        mars_hours = []
        for _ in range(full_cycles):
            mars_hours.extend(hours)
        
        # Add any remaining hours to complete the length
        mars_hours.extend(hours[:remainder])
        
        return np.array(mars_hours)

    # mars_soy: Values from 469 to 514 with each value repeating 24 times
    def generate_mars_soy(length):
        soy_values = []
        
        # Each Mars sol of year (SOY) has 24 hours
        for soy in range(469, 515):  # From 469 to 514
            soy_values.extend([soy] * 24)  # Repeat each SOY 24 times (for each hour)
        
        return np.array(soy_values[:length])


    # Generate all the arrays
    earth_month_values = generate_earth_months(1104)
    earth_second_values = generate_earth_seconds(1104)
    emars_sol_values = generate_emars_sol(1104)
    macda_sol_values = generate_macda_sol(1104)
    mars_hour_values = generate_mars_hours(1104)
    mars_soy_values = generate_mars_soy(1104)

    # Linear arrays:
    time_values = np.linspace(0, 1.103e+03, 1104).tolist()
    Ls_values = np.linspace(239.92, 269.82, 1104).tolist()
    MY_values = np.full(1104, 28.0).tolist()
    lat_values = np.linspace(-88.71428571, 88.71428571, 36).tolist()
    latu_values = np.linspace(-87.42857143, 87.42857143, 36).tolist()
    lon_values = np.linspace(3, 357, 60).tolist()
    lonv_values = np.linspace(0, 354, 60).tolist()
    earth_year_values = np.full(1104, 2007.0).tolist()

    # AK: non-linear sequence with 29 values
    ak_values = [2.0000000e-02, 5.7381272e-02, 1.9583981e-01, 5.9229583e-01, 1.5660228e+00,
                2.4454966e+00, 2.7683754e+00, 2.8851693e+00, 2.9172227e+00, 2.9087038e+00,
                2.8598938e+00, 2.7687652e+00, 2.6327055e+00, 2.4509220e+00, 2.2266810e+00,
                1.9684681e+00, 1.6894832e+00, 1.4055812e+00, 1.1324258e+00, 8.8289177e-01,
                6.6548467e-01, 4.8401019e-01, 3.3824119e-01, 2.2510704e-01, 1.3995719e-01,
                7.7611551e-02, 3.3085503e-02, 2.0000001e-03, 0.0000000e+00]

    # BK: non-linear sequence with 29 values
    bk_values = [0.0, 0.0, 0.0, 0.0, 0.0, 0.00193664, 0.00744191, 0.01622727, 0.02707519,
                0.043641, 0.0681068, 0.1028024, 0.14971954, 0.20987134, 0.28270233,
                0.3658161, 0.4552023, 0.545936, 0.6331097, 0.7126763, 0.7819615,
                0.8397753, 0.88620347, 0.9222317, 0.94934535, 0.9691962, 0.98337257,
                0.9932694, 1.0]

    # PFULL: 28 non-linear values
    pfull_values = [3.54665839e-04, 1.12789917e-03, 3.58229615e-03, 1.00147974e-02,
                    2.57178441e-02, 5.92796833e-02, 1.16012250e-01, 1.92695452e-01,
                    2.96839262e-01, 4.52589921e-01, 6.77446304e-01, 9.88319228e-01,
                    1.39717102e+00, 1.90617687e+00, 2.50426698e+00, 3.16685550e+00,
                    3.85940946e+00, 4.54382278e+00, 5.18537078e+00, 5.75801233e+00,
                    6.24681227e+00, 6.64754041e+00, 6.96437863e+00, 7.20689692e+00,
                    7.38721125e+00, 7.51781224e+00, 7.61018405e+00, 7.67406808e+00]

    # PHALF: 29 non-linear values
    phalf_values = [2.00000000e-04, 5.73812730e-04, 1.95839810e-03, 5.92295800e-03,
                    1.56602280e-02, 3.93670884e-02, 8.49864874e-02, 1.53801648e-01,
                    2.37651206e-01, 3.65122739e-01, 5.53021330e-01, 8.19266132e-01,
                    1.17916751e+00, 1.64051846e+00, 2.19907475e+00, 2.83646865e+00,
                    3.52195254e+00, 4.21776294e+00, 4.88626895e+00, 5.49643635e+00,
                    6.02775847e+00, 6.47110991e+00, 6.82714898e+00, 7.10343501e+00,
                    7.31135861e+00, 7.46358670e+00, 7.57229980e+00, 7.64819446e+00,
                    7.70000000e+00]

    earth_month_var[:] = earth_month_values
    earth_second_var[:] = earth_second_values
    emars_sol_var[:] = emars_sol_values
    macda_sol_var[:] = macda_sol_values
    mars_hour_var[:] = mars_hour_values
    mars_soy_var[:] = mars_soy_values
    time_var[:] = time_values
    Ls_var[:] = Ls_values
    MY_var[:] = MY_values
    lat_var[:] = lat_values
    latu_var[:] = latu_values
    lon_var[:] = lon_values
    lonv_var[:] = lonv_values
    earth_year_var[:] = earth_year_values
    ak_var[:] = ak_values
    bk_var[:] = bk_values
    pfull_var[:] = pfull_values
    phalf_var[:] = phalf_values

    # Create helper function for the rest of the variables
    def create_var(name, dimensions, units, longname, min_val, max_val, data_type=np.float32):
        var = nc_file.createVariable(name, data_type, dimensions)
        var.units = units
        var.long_name = longname
        
        shape = tuple(nc_file.dimensions[dim].size for dim in dimensions)
        var[:] = np.random.uniform(min_val, max_val, shape)
        
        return var
    
    # Create each variable as found in the real EMARS files
    create_var('Surface_geopotential', ('lat', 'lon'), 'm^2/s/s', 'surface geopotential height', -24000.0, 26000.0)
    create_var('T', ('time', 'pfull', 'lat', 'lon'), 'K', 'Temperature', 102.4, 291.9)
    create_var('U', ('time', 'pfull', 'latu', 'lon'), 'm/s', 'zonal wind', -257.6, 402.0)
    create_var('V', ('time', 'pfull', 'lat', 'lonv'), 'm/s', 'meridional wind', -278.8, 422.0)
    create_var('ps', ('time', 'lat', 'lon'), 'pascal', 'surface pressure', 312.1, 1218.1)
    
    nc_file.close()
    print("Created emars_test.nc")

# ----------------------------------------------------------------------
#                      OpenMARS Dummy File
# ----------------------------------------------------------------------

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
    
    # Create linear variables
    lon_var = nc_file.createVariable('lon', 'f4', ('lon',))
    lat_var = nc_file.createVariable('lat', 'f4', ('lat',))
    lev_var = nc_file.createVariable('lev', 'f4', ('lev',))
    time_var = nc_file.createVariable('time', 'f4', ('time',))
    Ls_var = nc_file.createVariable('Ls', 'f4', ('time',))
    MY_var = nc_file.createVariable('MY', 'f4', ('time',))
    
    # Use these exact values from the real file
    lon_values = np.linspace(-180., 175., 72).tolist()
    lat_values = np.linspace(87.49999, -87.49999, 36).tolist()
    lev_values = np.linspace(9.9949998e-01, 5.0824954e-05, 35).tolist()
    time_values = np.linspace(3181.0833, 3211., 360).tolist()
    Ls_values = np.linspace(264.93198, 284.14746, 360).tolist()
    MY_values = np.linspace(28.0, 28.0, 360).tolist()
    
    lon_var[:] = lon_values
    lat_var[:] = lat_values
    lev_var[:] = lev_values
    time_var[:] = time_values
    Ls_var[:] = Ls_values
    MY_var[:] = MY_values
    
    # Create each variable as found in real OpenMARS files
    create_var('ps', ('time', 'lat', 'lon'), '', 214.5, 1133.5)
    create_var('tsurf', ('time', 'lat', 'lon'), '', 145.5, 309.9)
    create_var('co2ice', ('time', 'lat', 'lon'), '', 0.0, 6860.4)
    create_var('dustcol', ('time', 'lat', 'lon'), '', 6.8e-09, 4.5)
    create_var('u', ('time', 'lev', 'lat', 'lon'), '', -517.1, 384.8)
    create_var('v', ('time', 'lev', 'lat', 'lon'), '', -362.2, 453.3)
    create_var('temp', ('time', 'lev', 'lat', 'lon'), '', 99.3, 299.4)
    
    nc_file.close()
    print("Created openmars_test.nc")

# ----------------------------------------------------------------------
#                      PCM Dummy File
# ----------------------------------------------------------------------

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
    index_dim = nc_file.createDimension('index', 100)
    
    # Create variables with longname and units
    ap_var = nc_file.createVariable('ap', 'f4', ('interlayer',))
    ap_var.units = "Pa"
    
    bp_var = nc_file.createVariable('bp', 'f4', ('interlayer',))
    bp_var.units = ""
    
    altitude_var = nc_file.createVariable('altitude', 'f4', ('altitude',))
    altitude_var.long_name = "pseudo-alt"
    altitude_var.units = "km"
    
    aps_var = nc_file.createVariable('aps', 'f4', ('altitude',))
    aps_var.units = "Pa"
    
    bps_var = nc_file.createVariable('bps', 'f4', ('altitude',))
    bps_var.units = ""
    
    controle_var = nc_file.createVariable('controle', 'f4', ('index',))
    controle_var.units = ""
    
    longitude_var = nc_file.createVariable('longitude', 'f4', ('longitude',))
    longitude_var.long_name = "East longitude"
    longitude_var.units = "degrees_east"
    
    Ls_var = nc_file.createVariable('Ls', 'f4', ('Time',))
    Ls_var.units = "deg"
    
    Sols_var = nc_file.createVariable('Sols', 'f4', ('Time',))
    Sols_var.units = "sols"
    
    Time_var = nc_file.createVariable('Time', 'f4', ('Time',))
    Time_var.long_name = "Time"
    Time_var.units = "days since 0000-00-0 00:00:00"
    
    soildepth_var = nc_file.createVariable('soildepth', 'f4', ('subsurface_layers',))
    soildepth_var.long_name = "Soil mid-layer depth"
    soildepth_var.units = "m"
    
    latitude_var = nc_file.createVariable('latitude', 'f4', ('latitude',))
    latitude_var.long_name = "North latitude"
    latitude_var.units = "degrees_north"

    # --------- Generate realistic values for PCM-like data ----------
    latitude_values = np.arange(90, -90.1, -3.75).tolist()
    longitude_values = np.arange(-180, 180.1, 5.625).tolist()
    ls_step = (280.42017 - 264.49323) / (100 - 1)
    sols_values = np.linspace(1175.2489, 1199.9989, 100).tolist()
    time_values = np.linspace(488.25, 513.00, 100).tolist()

    # controle: 100 values with first 11 specified and the rest zero
    def generate_controle(length=100):
        # First 11 specific values
        controle_values = [6.40000000e+01, 4.80000000e+01, 4.90000000e+01, 6.87000000e+02,
                        3.39720000e+06, 7.07765139e-05, 3.72000003e+00, 4.34899979e+01,
                        2.56792992e-01, 8.87750000e+04, 9.24739583e+02]
        
        # Pad with zeros to reach desired length
        controle_values.extend([0.0] * (length - len(controle_values)))
        
        return controle_values

    controle_values = generate_controle()

    altitude_values = [4.48063861e-03, 2.35406722e-02, 7.47709209e-02, 1.86963522e-01,
                    3.97702855e-01, 7.38866768e-01, 1.20169999e+00, 1.73769704e+00,
                    2.31003686e+00, 2.91118833e+00, 3.54250751e+00, 4.20557844e+00,
                    4.90199931e+00, 5.63340036e+00, 6.40156335e+00, 7.20835024e+00,
                    8.05566358e+00, 8.94555555e+00, 9.88013681e+00, 1.08616622e+01,
                    1.18925317e+01, 1.29751718e+01, 1.41121552e+01, 1.53062261e+01,
                    1.65602437e+01, 1.78772103e+01, 1.92602494e+01, 2.07126928e+01,
                    2.22379871e+01, 2.38398001e+01, 2.55219722e+01, 2.72884562e+01,
                    2.91434648e+01, 3.10914678e+01, 3.31371307e+01, 3.52852142e+01,
                    3.75408317e+01, 3.99094092e+01, 4.23965466e+01, 4.50081204e+01,
                    4.77503295e+01, 5.06380414e+01, 5.37235185e+01, 5.70999891e+01,
                    6.08520244e+01, 6.50463033e+01, 6.97653980e+01, 7.51124951e+01,
                    8.04595923e+01]

    aps_values = [4.5553492e-03, 2.3924896e-02, 7.5912692e-02, 1.8933944e-01, 4.0066403e-01,
                7.3744106e-01, 1.1827117e+00, 1.6806941e+00, 2.1907256e+00, 2.7015827e+00,
                3.2102795e+00, 3.7139201e+00, 4.2095418e+00, 4.6941628e+00, 5.1648922e+00,
                5.6189032e+00, 6.0534425e+00, 6.4659243e+00, 6.8539081e+00, 7.2151675e+00,
                7.5477109e+00, 7.8497639e+00, 8.1198349e+00, 8.3567305e+00, 8.5595417e+00,
                8.7276497e+00, 8.8607130e+00, 8.9586449e+00, 9.0215597e+00, 9.0497084e+00,
                9.0433722e+00, 9.0027008e+00, 8.9274740e+00, 8.8166723e+00, 8.6677418e+00,
                8.4750948e+00, 8.2265606e+00, 7.8932257e+00, 7.3903952e+00, 6.4670715e+00,
                5.1399899e+00, 3.8560932e+00, 2.8323510e+00, 2.0207324e+00, 1.3885452e+00,
                9.1286123e-01, 5.6945199e-01, 3.3360735e-01, 1.9544031e-01]

    bps_values = [9.9954456e-01, 9.9760950e-01, 9.9242634e-01, 9.8116696e-01, 9.6035337e-01,
                9.2756802e-01, 8.8483083e-01, 8.3773518e-01, 7.9014516e-01, 7.4299800e-01,
                6.9643623e-01, 6.5059197e-01, 6.0560304e-01, 5.6160903e-01, 5.1874298e-01,
                4.7713464e-01, 4.3691111e-01, 3.9818937e-01, 3.6107957e-01, 3.2567981e-01,
                2.9207525e-01, 2.6034081e-01, 2.3053549e-01, 2.0270133e-01, 1.7686437e-01,
                1.5303348e-01, 1.3120057e-01, 1.1133941e-01, 9.3407877e-02, 7.7347368e-02,
                6.3085094e-02, 5.0536096e-02, 3.9604262e-02, 3.0185465e-02, 2.2171425e-02,
                1.5454679e-02, 9.9357506e-03, 5.5426015e-03, 2.2971667e-03, 4.9822254e-04,
                1.1593266e-05, 1.8491624e-09, 1.1037076e-16, 2.7358622e-31, 0.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00]

    ap_values = [0.0, 0.0091107, 0.03873909, 0.11308629, 0.26559258, 0.5357355,
                0.93914664, 1.4262768, 1.9351114, 2.4463396, 2.956826, 3.4637327,
                3.9641075, 4.454976, 4.9333496, 5.3964353, 5.841371, 6.2655144,
                6.6663346, 7.0414815, 7.388854, 7.706568, 7.99296, 8.246711,
                8.466751, 8.652332, 8.802968, 8.918459, 8.998832, 9.044289,
                9.055129, 9.031614, 8.973787, 8.88116, 8.752184, 8.583301,
                8.36689, 8.08623, 7.700221, 7.080569, 5.8535743, 4.4264054,
                3.285781, 2.3789213, 1.6625438, 1.1145465, 0.71117604, 0.4277279,
                0.23948681, 0.0]

    bp_values = [1.0000000e+00, 9.9908912e-01, 9.9612981e-01, 9.8872286e-01, 9.7361106e-01,
                9.4709563e-01, 9.0804040e-01, 8.6162120e-01, 8.1384915e-01, 7.6644123e-01,
                7.1955484e-01, 6.7331761e-01, 6.2786639e-01, 5.8333969e-01, 5.3987837e-01,
                4.9760753e-01, 4.5666179e-01, 4.1716042e-01, 3.7921831e-01, 3.4294087e-01,
                3.0841875e-01, 2.7573177e-01, 2.4494988e-01, 2.1612112e-01, 1.8928154e-01,
                1.6444720e-01, 1.4161976e-01, 1.2078137e-01, 1.0189746e-01, 8.4918290e-02,
                6.9776453e-02, 5.6393731e-02, 4.4678457e-02, 3.4530070e-02, 2.5840862e-02,
                1.8501990e-02, 1.2407369e-02, 7.4641323e-03, 3.6210711e-03, 9.7326230e-04,
                2.3182834e-05, 3.6983245e-09, 2.2074152e-16, 5.4717245e-31, 0.0000000e+00,
                0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00, 0.0000000e+00]

    soildepth_values = [1.41421348e-04, 2.82842695e-04, 5.65685390e-04, 1.13137078e-03,
                        2.26274156e-03, 4.52548312e-03, 9.05096624e-03, 1.81019325e-02,
                        3.62038650e-02, 7.24077299e-02, 1.44815460e-01, 2.89630920e-01,
                        5.79261839e-01, 1.15852368e+00, 2.31704736e+00, 4.63409472e+00,
                        9.26818943e+00, 1.85363789e+01]

    ls_values = [264.49323, 264.65536, 264.81747, 264.97955, 265.14163, 265.30368, 265.46573,
            265.62775, 265.78973, 265.95172, 266.11365, 266.2756, 266.4375, 266.5994,
            266.76126, 266.9231, 267.08493, 267.24673, 267.4085, 267.57028, 267.732,
            267.8937, 268.05542, 268.21707, 268.37872, 268.54034, 268.70193, 268.8635,
            269.02502, 269.18655, 269.34805, 269.50952, 269.67096, 269.8324, 269.99377,
            270.15515, 270.3165, 270.4778, 270.6391, 270.80035, 270.9616, 271.1228,
            271.284, 271.44516, 271.60626, 271.76736, 271.92844, 272.08948, 272.25052,
            272.4115, 272.57245, 272.73337, 272.8943, 273.05515, 273.21597, 273.3768,
            273.53757, 273.69833, 273.85904, 274.01974, 274.1804, 274.34103, 274.50162,
            274.6622, 274.82272, 274.98322, 275.1437, 275.30414, 275.46454, 275.6249,
            275.78525, 275.94556, 276.10583, 276.26608, 276.4263, 276.58646, 276.7466,
            276.9067, 277.0668, 277.22684, 277.38684, 277.5468, 277.70676, 277.86664,
            278.02652, 278.18634, 278.34613, 278.5059, 278.66562, 278.82532, 278.98495,
            279.14456, 279.30414, 279.46368, 279.6232, 279.78265, 279.94208, 280.10147,
            280.26083, 280.42017]

    latitude_var[:] = latitude_values
    longitude_var[:] = longitude_values
    Ls_var[:] = ls_values
    Sols_var[:] = sols_values
    Time_var[:] = time_values
    controle_var[:] = controle_values
    altitude_var[:] = altitude_values
    aps_var[:] = aps_values
    bps_var[:] = bps_values
    ap_var[:] = ap_values
    bp_var[:] = bp_values
    soildepth_var[:] = soildepth_values
    
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
    create_var('Mccntot', ('Time', 'latitude', 'longitude'), 'kg/m2', '', 1.0e-18, 2.2e-03)
    create_var('Nccntot', ('Time', 'latitude', 'longitude'), 'Nbr/m2', '', 8.0e+01, 2.4e+11)
    create_var('aire', ('latitude', 'longitude'), '', '', 6.1e+08, 7.4e+10)
    create_var('albedo', ('Time', 'latitude', 'longitude'), '', '', 0.1, 0.9)
    create_var('ccnN', ('Time', 'altitude', 'latitude', 'longitude'), 'part/kg', '', -1.6e+06, 2.0e+09)
    create_var('ccnq', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', -4.1e-08, 5.3e-05)
    create_var('co2', ('Time', 'altitude', 'latitude', 'longitude'), 'kg/kg', '', 0.9, 1.0)
    create_var('co2ice', ('Time', 'latitude', 'longitude'), 'kg.m-2', '', 0.0, 2189.0)
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

# ----------------------------------------------------------------------
#                      MarsWRF Dummy File
# ----------------------------------------------------------------------

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
    
    # Create and populate ZNU variable (eta values on half/mass levels)
    ZNU_var = nc_file.createVariable('ZNU', 'f4', ('Time', 'bottom_top'))
    znu_values = np.array([0.99916667, 0.99666667, 0.9916667, 0.9816667, 0.9625, 0.9375,
                          0.9125, 0.8875, 0.8625, 0.8375, 0.8125, 0.7875,
                          0.7625, 0.7375, 0.7125, 0.6875, 0.6625, 0.6375,
                          0.6125, 0.5875, 0.5625, 0.5375, 0.5125, 0.4875,
                          0.4625, 0.4375, 0.4125, 0.3875, 0.3625, 0.3375,
                          0.3125, 0.2875, 0.2625, 0.23750001, 0.2125, 0.1875,
                          0.1625, 0.13749999, 0.11250001, 0.08750001, 0.0625, 0.03749999,
                          0.01249999])
    # Repeat the same values for all 100 timesteps
    ZNU_var[:] = np.tile(znu_values, (100, 1))
    ZNU_var.description = "eta values on half (mass) levels"

    # Create and populate FNM variable (upper weight for vertical stretching)
    FNM_var = nc_file.createVariable('FNM', 'f4', ('Time', 'bottom_top'))
    fnm_values = np.array([0., 0.33333334, 0.33333334, 0.33333334, 0.34782556, 0.5000006,
                          0.4999994, 0.5000006, 0.5, 0.4999994, 0.5000006, 0.4999994,
                          0.5000006, 0.5, 0.4999994, 0.5000006, 0.4999994, 0.5000006,
                          0.5, 0.4999994, 0.5000006, 0.4999994, 0.5000006, 0.5,
                          0.4999994, 0.5000006, 0.4999994, 0.5000006, 0.5, 0.4999994,
                          0.5000006, 0.4999994, 0.5000006, 0.5, 0.4999994, 0.5000006,
                          0.4999994, 0.5000006, 0.5, 0.4999994, 0.5000006, 0.4999994,
                          0.5000006])
    FNM_var[:] = np.tile(fnm_values, (100, 1))
    FNM_var.description = "upper weight for vertical stretching"

    # Create and populate FNP variable (lower weight for vertical stretching)
    FNP_var = nc_file.createVariable('FNP', 'f4', ('Time', 'bottom_top'))
    fnp_values = np.array([0., 0.6666667, 0.6666667, 0.6666667, 0.6521745, 0.4999994, 0.5000006,
                          0.4999994, 0.5, 0.5000006, 0.4999994, 0.5000006, 0.4999994, 0.5,
                          0.5000006, 0.4999994, 0.5000006, 0.4999994, 0.5, 0.5000006, 0.4999994,
                          0.5000006, 0.4999994, 0.5, 0.5000006, 0.4999994, 0.5000006, 0.4999994,
                          0.5, 0.5000006, 0.4999994, 0.5000006, 0.4999994, 0.5, 0.5000006,
                          0.4999994, 0.5000006, 0.4999994, 0.5, 0.5000006, 0.4999994, 0.5000006,
                          0.4999994])
    FNP_var[:] = np.tile(fnp_values, (100, 1))
    FNP_var.description = "lower weight for vertical stretching"

    # Create and populate RDNW variable (inverse d(eta) values between full/w levels)
    RDNW_var = nc_file.createVariable('RDNW', 'f4', ('Time', 'bottom_top'))
    rdnw_values = np.array([-600.00055, -300.00027, -150.00014, -75.00007, -39.999943, -40.00004,
                            -39.999943, -40.00004, -40.00004, -39.999943, -40.00004, -39.999943,
                            -40.00004, -40.00004, -39.999943, -40.00004, -39.999943, -40.00004,
                            -40.00004, -39.999943, -40.00004, -39.999943, -40.00004, -40.00004,
                            -39.999943, -40.00004, -39.999943, -40.00004, -40.00004, -39.999943,
                            -40.00004, -39.999943, -40.00004, -40.00004, -39.999943, -40.00004,
                            -39.999943, -40.00004, -40.00004, -39.999943, -40.00004, -39.999943,
                            -40.00004])
    RDNW_var[:] = np.tile(rdnw_values, (100, 1))
    RDNW_var.description = "inverse d(eta) values between full (w) levels"

    # Create and populate RDN variable (inverse d(eta) values between half/mass levels)
    RDN_var = nc_file.createVariable('RDN', 'f4', ('Time', 'bottom_top'))
    rdn_values = np.array([0., -400.0004, -200.0002, -100.0001, -52.17388, -39.999992,
                          -39.999992, -39.999992, -40.00004, -39.999992, -39.999992, -39.999992,
                          -39.999992, -40.00004, -39.999992, -39.999992, -39.999992, -39.999992,
                          -40.00004, -39.999992, -39.999992, -39.999992, -39.999992, -40.00004,
                          -39.999992, -39.999992, -39.999992, -39.999992, -40.00004, -39.999992,
                          -39.999992, -39.999992, -39.999992, -40.00004, -39.999992, -39.999992,
                          -39.999992, -39.999992, -40.00004, -39.999992, -39.999992, -39.999992,
                          -39.999992])
    RDN_var[:] = np.tile(rdn_values, (100, 1))
    RDN_var.description = "inverse d(eta) values between half (mass) levels"

    # Create and populate DNW variable (d(eta) values between full/w levels)
    DNW_var = nc_file.createVariable('DNW', 'f4', ('Time', 'bottom_top'))
    dnw_values = np.array([-0.00166667, -0.00333333, -0.00666666, -0.01333332, -0.02500004, -0.02499998,
                          -0.02500004, -0.02499998, -0.02499998, -0.02500004, -0.02499998, -0.02500004,
                          -0.02499998, -0.02499998, -0.02500004, -0.02499998, -0.02500004, -0.02499998,
                          -0.02499998, -0.02500004, -0.02499998, -0.02500004, -0.02499998, -0.02499998,
                          -0.02500004, -0.02499998, -0.02500004, -0.02499998, -0.02499998, -0.02500004,
                          -0.02499998, -0.02500004, -0.02499998, -0.02499998, -0.02500004, -0.02499998,
                          -0.02500004, -0.02499998, -0.02499998, -0.02500004, -0.02499998, -0.02500004,
                          -0.02499998])
    DNW_var[:] = np.tile(dnw_values, (100, 1))
    DNW_var.description = "d(eta) values between full (w) levels"

    # Create and populate DN variable (d(eta) values between half/mass levels)
    DN_var = nc_file.createVariable('DN', 'f4', ('Time', 'bottom_top'))
    dn_values = np.array([0., -0.0025, -0.005, -0.00999999, -0.01916668, -0.02500001,
                         -0.02500001, -0.02500001, -0.02499998, -0.02500001, -0.02500001, -0.02500001,
                         -0.02500001, -0.02499998, -0.02500001, -0.02500001, -0.02500001, -0.02500001,
                         -0.02499998, -0.02500001, -0.02500001, -0.02500001, -0.02500001, -0.02499998,
                         -0.02500001, -0.02500001, -0.02500001, -0.02500001, -0.02499998, -0.02500001,
                         -0.02500001, -0.02500001, -0.02500001, -0.02499998, -0.02500001, -0.02500001,
                         -0.02500001, -0.02500001, -0.02499998, -0.02500001, -0.02500001, -0.02500001,
                         -0.02500001])
    DN_var[:] = np.tile(dn_values, (100, 1))
    DN_var.description = "d(eta) values between half (mass) levels"

    # Create and populate ZNW variable (eta values on full/w levels)
    ZNW_var = nc_file.createVariable('ZNW', 'f4', ('Time', 'bottom_top_stag'))
    znw_values = np.array([1., 0.99833333, 0.995, 0.98833334, 0.975, 0.95,
                          0.925, 0.9, 0.875, 0.85, 0.825, 0.8,
                          0.775, 0.75, 0.725, 0.7, 0.675, 0.65,
                          0.625, 0.6, 0.575, 0.55, 0.525, 0.5,
                          0.47500002, 0.45, 0.425, 0.39999998, 0.375, 0.35000002,
                          0.325, 0.3, 0.27499998, 0.25, 0.22500002, 0.19999999,
                          0.17500001, 0.14999998, 0.125, 0.10000002, 0.07499999, 0.05000001,
                          0.02499998, 0.])
    ZNW_var[:] = np.tile(znw_values, (100, 1))
    ZNW_var.description = "eta values on full (w) levels"
    
    # Generate data for the Time dimension variables
    
    # RDX and RDY: 100 constant values
    RDX_var = nc_file.createVariable('RDX', 'f4', ('Time',))
    RDX_var[:] = np.full(100, 8.450905e-06)
    RDX_var.description = "INVERSE X GRID LENGTH"
    
    RDY_var = nc_file.createVariable('RDY', 'f4', ('Time',))
    RDY_var[:] = np.full(100, 8.450905e-06)
    RDY_var.description = "INVERSE Y GRID LENGTH"
    
    # DTS, DTSEPS, RESM, ZETATOP, T00, P00, TLP, TISO: 100 zeros
    DTS_var = nc_file.createVariable('DTS', 'f4', ('Time',))
    DTS_var[:] = np.zeros(100)
    DTS_var.description = "SMALL TIMESTEP"
    
    DTSEPS_var = nc_file.createVariable('DTSEPS', 'f4', ('Time',))
    DTSEPS_var[:] = np.zeros(100)
    DTSEPS_var.description = "TIME WEIGHT CONSTANT FOR SMALL STEPS"
    
    RESM_var = nc_file.createVariable('RESM', 'f4', ('Time',))
    RESM_var[:] = np.zeros(100)
    RESM_var.description = "TIME WEIGHT CONSTANT FOR SMALL STEPS"
    
    ZETATOP_var = nc_file.createVariable('ZETATOP', 'f4', ('Time',))
    ZETATOP_var[:] = np.zeros(100)
    ZETATOP_var.description = "ZETA AT MODEL TOP"
    
    T00_var = nc_file.createVariable('T00', 'f4', ('Time',))
    T00_var[:] = np.zeros(100)
    T00_var.description = "BASE STATE TEMPERATURE"
    T00_var.units = "K"
    
    P00_var = nc_file.createVariable('P00', 'f4', ('Time',))
    P00_var[:] = np.zeros(100)
    P00_var.description = "BASE STATE PRESURE"
    P00_var.units = "Pa"
    
    TLP_var = nc_file.createVariable('TLP', 'f4', ('Time',))
    TLP_var[:] = np.zeros(100)
    TLP_var.description = "BASE STATE LAPSE RATE"
    
    TISO_var = nc_file.createVariable('TISO', 'f4', ('Time',))
    TISO_var[:] = np.zeros(100)
    TISO_var.description = "TEMP AT WHICH THE BASE T TURNS CONST"
    TISO_var.units = "K"
    
    # CF1, CF2, CF3: 100 constant values
    CF1_var = nc_file.createVariable('CF1', 'f4', ('Time',))
    CF1_var[:] = np.full(100, 1.5555556)
    CF1_var.description = "2nd order extrapolation constant"
    
    CF2_var = nc_file.createVariable('CF2', 'f4', ('Time',))
    CF2_var[:] = np.full(100, -0.6666667)
    CF2_var.description = "2nd order extrapolation constant"
    
    CF3_var = nc_file.createVariable('CF3', 'f4', ('Time',))
    CF3_var[:] = np.full(100, 0.11111111)
    CF3_var.description = "2nd order extrapolation constant"
    
    # ITIMESTEP and XTIME: 100 values incremented by 14400
    ITIMESTEP_var = nc_file.createVariable('ITIMESTEP', 'f4', ('Time',))
    ITIMESTEP_values = np.arange(961920, 961920 + 14400 * 100, 14400)
    ITIMESTEP_values = [int(x) for x in ITIMESTEP_values]  # Convert to int
    ITIMESTEP_var[:] = ITIMESTEP_values
    
    XTIME_var = nc_file.createVariable('XTIME', 'f4', ('Time',))
    XTIME_var[:] = ITIMESTEP_values  # Same as ITIMESTEP
    XTIME_var.units = "minutes"
    XTIME_var.description = "minutes since simulation start"
    
    # JULIAN: 100 values with specific pattern
    JULIAN_var = nc_file.createVariable('JULIAN', 'f4', ('Time',))
    julian_values = []
    current_value = 668
    
    for i in range(100):
        julian_values.append(current_value)
        
        if current_value == 659:
            current_value = 0
        elif current_value == 668:
            current_value = 9
        else:
            current_value += 10
    
    JULIAN_var[:] = np.array(julian_values)
    JULIAN_var.units = "days"
    JULIAN_var.description = "day of year, 0.0 at 0Z on 1 Jan."
    
    # L_S: 100 values representing solar longitude
    L_S_var = nc_file.createVariable('L_S', 'f4', ('Time',))
    L_S_var[:] = np.array([359.48444, 4.6024184, 9.639707, 14.601253, 19.492262, 24.318121,
                          29.084349, 33.796543, 38.460365, 43.0815, 47.665646, 52.2185,
                          56.745754, 61.25309, 65.74617, 70.23066, 74.71221, 79.19647,
                          83.689095, 88.19573, 92.72206, 97.27376, 101.856514, 106.47602,
                          111.137985, 115.84809, 120.612, 125.43531, 130.32355, 135.2821,
                          140.3162, 145.4308, 150.63055, 155.91972, 161.30203, 166.78061,
                          172.35777, 178.035, 183.81273, 189.69017, 195.66527, 201.73451,
                          207.89287, 214.13376, 220.44897, 226.82884, 233.26224, 239.73683,
                          246.23932, 252.75572, 259.2717, 265.77292, 272.24545, 278.676,
                          285.05234, 291.36343, 297.59964, 303.7529, 309.8167, 315.78613,
                          321.6577, 327.4295, 333.1008, 338.6721, 344.1449, 349.5216,
                          354.8054, 360.0, 5.1096992, 10.139186, 15.09344, 19.97769,
                          24.797337, 29.557905, 34.26501, 38.924305, 43.541485, 48.122246,
                          52.672283, 57.197292, 61.702946, 66.194916, 70.678856, 75.16042,
                          79.64526, 84.13903, 88.647385, 93.175995, 97.730545, 102.31672,
                          106.940216, 111.606735, 116.32197, 121.09157, 125.92112, 130.81615,
                          135.78203, 140.82396, 145.94687, 151.15538])
    L_S_var.units = "degrees"
    L_S_var.description = "Planetocentric solar Longitude"
    
    # DECLIN: 100 values with declination angles
    DECLIN_var = nc_file.createVariable('DECLIN', 'f4', ('Time',))
    DECLIN_var[:] = np.array([-3.86990025e-03, 3.41197848e-02, 7.12934360e-02, 1.07464984e-01,
                              1.42467186e-01, 1.76147699e-01, 2.08365589e-01, 2.38988355e-01,
                              2.67889649e-01, 2.94947445e-01, 3.20042819e-01, 3.43059152e-01,
                              3.63882095e-01, 3.82399738e-01, 3.98503423e-01, 4.12089020e-01,
                              4.23058301e-01, 4.31320816e-01, 4.36795801e-01, 4.39414054e-01,
                              4.39119637e-01, 4.35871571e-01, 4.29645032e-01, 4.20432001e-01,
                              4.08241928e-01, 3.93101573e-01, 3.75054955e-01, 3.54163110e-01,
                              3.30503553e-01, 3.04170370e-01, 2.75274187e-01, 2.43943155e-01,
                              2.10323900e-01, 1.74583584e-01, 1.36912495e-01, 9.75274146e-02,
                              5.66755012e-02, 1.46388495e-02, -2.82607116e-02, -7.16568232e-02,
                              -1.15134262e-01, -1.58225715e-01, -2.00410619e-01, -2.41116881e-01,
                              -2.79727042e-01, -3.15589964e-01, -3.48039240e-01, -3.76418710e-01,
                              -4.00114596e-01, -4.18592125e-01, -4.31432575e-01, -4.38365251e-01,
                              -4.39289480e-01, -4.34281319e-01, -4.23584312e-01, -4.07585740e-01,
                              -3.86782587e-01, -3.61743599e-01, -3.33072543e-01, -3.01376730e-01,
                              -2.67242700e-01, -2.31219843e-01, -1.93810493e-01, -1.55465990e-01,
                              -1.16586491e-01, -7.75237307e-02, -3.85853238e-02, -3.99982528e-05,
                              3.78771052e-02, 7.49585778e-02, 1.11020446e-01, 1.45897120e-01,
                              1.79437563e-01, 2.11501762e-01, 2.41957992e-01, 2.70680368e-01,
                              2.97547251e-01, 3.22439909e-01, 3.45242023e-01, 3.65839422e-01,
                              3.84120494e-01, 3.99976969e-01, 4.13305223e-01, 4.24007773e-01,
                              4.31995004e-01, 4.37187225e-01, 4.39516455e-01, 4.38928276e-01,
                              4.35383230e-01, 4.28858101e-01, 4.19346660e-01, 4.06860083e-01,
                              3.91426831e-01, 3.73092681e-01, 3.51920277e-01, 3.27988863e-01,
                              3.01394105e-01, 2.72248387e-01, 2.40681618e-01, 2.06842363e-01])
    DECLIN_var.units = "radians"
    DECLIN_var.description = "SOLAR DECLINATION"
    
    # SOLCON: 100 values representing solar constant
    SOLCON_var = nc_file.createVariable('SOLCON', 'f4', ('Time',))
    SOLCON_var[:] = np.array([567.3604, 558.22437, 549.61774, 541.57135, 534.1084, 527.2457, 520.99493,
                              515.36365, 510.3564, 505.9753, 502.22095, 499.0928, 496.58984, 494.7108,
                              493.4544, 492.81976, 492.80643, 493.41434, 494.64398, 496.49628, 498.97238,
                              502.0736, 505.80106, 510.1552, 515.13556, 520.7401, 526.96436, 533.80096,
                              541.23846, 549.2602, 557.8433, 566.95734, 576.5625, 586.60864, 597.0333,
                              607.76074, 618.7003, 629.74603, 640.7758, 651.65247, 662.2247, 672.3293,
                              681.79535, 690.4486, 698.11743, 704.6395, 709.8694, 713.6849, 715.994,
                              716.74, 715.90405, 713.5071, 709.6082, 704.30115, 697.70953, 689.9801,
                              681.2758, 671.7687, 661.6329, 651.03906, 640.1497, 629.1154, 618.0727,
                              607.1425, 596.43005, 586.025, 576.00244, 566.424, 557.3393, 548.7875,
                              540.79846, 533.3949, 526.59296, 520.4038, 514.8348, 509.8901, 505.5717,
                              501.87997, 498.81442, 496.3739, 494.55713, 493.36298, 492.7905, 492.83926,
                              493.5093, 494.80118, 496.7158, 499.25436, 502.41818, 506.20825, 510.6251,
                              515.66797, 521.3347, 527.62067, 534.51794, 542.01465, 550.0937, 558.7314,
                              567.8965, 577.5482])
    SOLCON_var.units = "W m-2"
    SOLCON_var.description = "SOLAR CONSTANT"
    
    # SUNBODY: 100 values representing Sun-body distance
    SUNBODY_var = nc_file.createVariable('SUNBODY', 'f4', ('Time',))
    SUNBODY_var[:] = np.array([1.5525658, 1.5652192, 1.5774267, 1.5891018, 1.6001652, 1.6105456, 1.6201782,
                               1.6290059, 1.6369777, 1.6440494, 1.6501831, 1.6553464, 1.6595129, 1.6626616,
                               1.6647768, 1.6658484, 1.665871, 1.6648444, 1.6627738, 1.6596693, 1.6555461,
                               1.6504252, 1.6443326, 1.6373005, 1.6293664, 1.6205746, 1.6109754, 1.600626,
                               1.5895904, 1.57794, 1.5657536, 1.5531176, 1.5401263, 1.5268815, 1.5134925,
                               1.5000758, 1.4867549, 1.4736584, 1.4609202, 1.4486768, 1.4370666, 1.4262266,
                               1.4162911, 1.4073881, 1.3996366, 1.3931441, 1.3880028, 1.3842875, 1.3820535,
                               1.3813341, 1.3821403, 1.38446, 1.3882581, 1.3934788, 1.4000458, 1.4078658,
                               1.416831, 1.4268216, 1.4377091, 1.4493592, 1.4616344, 1.4743968, 1.4875096,
                               1.5008395, 1.5142577, 1.5276415, 1.5408748, 1.5538486, 1.5664614, 1.5786195,
                               1.5902369, 1.6012352, 1.6115435, 1.6210982, 1.6298423, 1.6377261, 1.6447057,
                               1.6507436, 1.6558082, 1.6598738, 1.6629198, 1.664931, 1.665898, 1.6658155,
                               1.6646842, 1.6625097, 1.6593025, 1.6550785, 1.6498592, 1.6436712, 1.636547,
                               1.6285251, 1.6196501, 1.6099732, 1.5995522, 1.5884517, 1.5767441, 1.5645088,
                               1.5518329, 1.5388116])
    SUNBODY_var.description = "Sun-planet distance in AU"
        
    # P_FIT_M: 100 constant values
    P_FIT_M_var = nc_file.createVariable('P_FIT_M', 'f4', ('Time',))
    P_FIT_M_var[:] = np.full(100, 1.1038098)
    P_FIT_M_var.units = "unitless"
    P_FIT_M_var.description = "SCALING OF P FOR CORRECT MASS"
    
    # P_TOP: 100 constant values
    P_TOP_var = nc_file.createVariable('P_TOP', 'f4', ('Time',))
    P_TOP_var[:] = np.full(100, 0.00567928)
    P_TOP_var.units = "Pa"
    P_TOP_var.description = "PRESSURE TOP OF THE MODEL"
    
    # MAX_MSTFX, MAX_MSTFY: 100 zeros
    MAX_MSTFX_var = nc_file.createVariable('MAX_MSTFX', 'f4', ('Time',))
    MAX_MSTFX_var[:] = np.zeros(100)
    MAX_MSTFX_var.description = "Max map factor in domain"
    
    MAX_MSTFY_var = nc_file.createVariable('MAX_MSTFY', 'f4', ('Time',))
    MAX_MSTFY_var[:] = np.zeros(100)
    MAX_MSTFY_var.description = "Max map factor in domain"
    
    # SEED1, SEED2, SAVE_TOPO_FROM_REAL: 100 zeros (integer type)
    SEED1_var = nc_file.createVariable('SEED1', 'i4', ('Time',))
    SEED1_var[:] = np.zeros(100, dtype=np.int32)
    SEED1_var.description = "RANDOM SEED NUMBER 1"
    
    SEED2_var = nc_file.createVariable('SEED2', 'i4', ('Time',))
    SEED2_var[:] = np.zeros(100, dtype=np.int32)
    SEED2_var.description = "RANDOM SEED NUMBER 2"
    
    SAVE_TOPO_FROM_REAL_var = nc_file.createVariable('SAVE_TOPO_FROM_REAL', 'i4', ('Time',))
    SAVE_TOPO_FROM_REAL_var[:] = np.zeros(100, dtype=np.int32)
    SAVE_TOPO_FROM_REAL_var.description = "1=original topo from real/0=topo modified by WRF"
    SAVE_TOPO_FROM_REAL_var.units = "flag"

    # Helper function to create a variable
    def create_var(name, dimensions, units, min_val, max_val, description="", data_type=np.float32, is_coordinate=False):
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
        var.description = description
        
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
        
        # Add stagger information for staggered variables
        if 'stag' in ''.join(dimensions):
            var.stagger = 'X' if 'west_east_stag' in dimensions else ('Y' if 'south_north_stag' in dimensions else 'Z')
        
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
    create_var('XLAT', ('Time', 'south_north', 'west_east'), 'degree_north', -89.0, 89.0, "LATITUDE, SOUTH IS NEGATIVE", is_coordinate=True)
    create_var('XLONG', ('Time', 'south_north', 'west_east'), 'degree_east', -179.0, 179.0, "LONGITUDE, WEST IS NEGATIVE", is_coordinate=True)
    create_var('XLAT_U', ('Time', 'south_north', 'west_east_stag'), 'degree_north', -89.0, 89.0, "LATITUDE, SOUTH IS NEGATIVE", is_coordinate=True)
    create_var('XLONG_U', ('Time', 'south_north', 'west_east_stag'), 'degree_east', -178.0, 180.0, "LONGITUDE, WEST IS NEGATIVE", is_coordinate=True)
    create_var('XLAT_V', ('Time', 'south_north_stag', 'west_east'), 'degree_north', -90.0, 90.0, "LATITUDE, SOUTH IS NEGATIVE", is_coordinate=True)
    create_var('XLONG_V', ('Time', 'south_north_stag', 'west_east'), 'degree_east', -179.0, 179.0, "LONGITUDE, WEST IS NEGATIVE", is_coordinate=True)
    
    # Create all variables from the MarsWRF file with proper descriptions
    create_var('ZS', ('Time', 'soil_layers_stag'), 'm', 0.0, 19.7, "DEPTHS OF CENTERS OF SOIL LAYERS")
    create_var('DZS', ('Time', 'soil_layers_stag'), 'm', 0.0, 11.2, "THICKNESSES OF SOIL LAYERS")
    create_var('U', ('Time', 'bottom_top', 'south_north', 'west_east_stag'), 'm s-1', -514.7, 519.5, "x-wind component")
    create_var('V', ('Time', 'bottom_top', 'south_north_stag', 'west_east'), 'm s-1', -205.6, 238.3, "y-wind component")
    create_var('W', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm s-1', -40.7, 39.0, "z-wind component")
    create_var('PH', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm2 s-2', -73000.0, 16000.0, "perturbation geopotential")
    create_var('PHB', ('Time', 'bottom_top_stag', 'south_north', 'west_east'), 'm2 s-2', -28000.0, 270000.0, "base-state geopotential")
    create_var('T', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K', -170.9, 349.1, "perturbation potential temperature (theta-t0)")
    create_var('T_INIT', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K', -101.3, 443.6, "initial potential temperature")
    create_var('MU', ('Time', 'south_north', 'west_east'), 'Pa', -267.5, 87.1, "perturbation dry air mass in column")
    create_var('MUB', ('Time', 'south_north', 'west_east'), 'Pa', 128.4, 1244.7, "base state dry air mass in column")
    create_var('P', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', -267.2, 87.0, "perturbation pressure")
    create_var('PB', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', 1.6, 1243.7, "BASE STATE PRESSURE")
    
    # Add all the other MarsWRF variables with proper descriptions
    create_var('P_HYD', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa', 1.1, 1330.7, "hydrostatic pressure")
    create_var('PSFC', ('Time', 'south_north', 'west_east'), 'Pa', 86.5, 1331.8, "SFC PRESSURE")
    create_var('T1_5', ('Time', 'south_north', 'west_east'), 'K', 142.1, 281.2, "TEMP at 1.5 M")
    create_var('TH1_5', ('Time', 'south_north', 'west_east'), 'K', 128.7, 378.1, "POT TEMP at 1.5 M")
    create_var('U1_5', ('Time', 'south_north', 'west_east'), 'm s-1', -21.6, 21.0, "U at 1.5 M")
    create_var('V1_5', ('Time', 'south_north', 'west_east'), 'm s-1', -20.0, 22.6, "V at 1.5 M")
    create_var('U1M', ('Time', 'south_north', 'west_east'), 'm s-1', -19.2, 19.0, "U at 1 M")
    create_var('V1M', ('Time', 'south_north', 'west_east'), 'm s-1', -17.8, 20.2, "V at 1 M")
    create_var('RHOS', ('Time', 'south_north', 'west_east'), 'kg m-3', 0.00213, 0.0424, "Surface air density")
    create_var('LANDMASK', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0, "LAND MASK (1 FOR LAND, 0 FOR WATER)")
    create_var('SLOPE', ('Time', 'south_north', 'west_east'), '', 6.4e-06, 1.7e-01, "ELEVATION SLOPE")
    create_var('SLP_AZI', ('Time', 'south_north', 'west_east'), 'rad', 0.0, 6.3, "ELEVATION SLOPE AZIMUTH")
    create_var('TSLB', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), 'K', 141.9, 308.7, "SOIL TEMPERATURE")
    create_var('SOIL_DEN', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 1500.0, 1500.0, "BULK DENSITY OF SOIL")
    create_var('SOIL_CAP', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 837.0, 837.0, "HEAT CAPACITY OF SOIL")
    create_var('SOIL_COND', ('Time', 'soil_layers_stag', 'south_north', 'west_east'), '', 0.0, 0.6, "CONDUCTIVITY OF SOIL")
    create_var('COSZEN', ('Time', 'south_north', 'west_east'), 'dimensionless', 0.0, 1.0, "COS of SOLAR ZENITH ANGLE")
    create_var('HRANG', ('Time', 'south_north', 'west_east'), 'radians', -3.1, 3.1, "SOLAR HOUR ANGLE")
    create_var('SUNFRAC', ('Time', 'south_north', 'west_east'), '', 0.0, 1.0, "Illuminated Fraction")
    create_var('COSZEN_MEAN', ('Time', 'south_north', 'west_east'), '', 0.0, 0.6, "Diurnally Averaged Cos Zenith angle")
    create_var('KPBL', ('Time', 'south_north', 'west_east'), '', 4.0, 43.0, "LEVEL OF PBL TOP")
    create_var('MAPFAC_MX', ('Time', 'south_north', 'west_east'), '', 1.0, 57.3, "Map scale factor on mass grid, x direction")
    create_var('MAPFAC_MY', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0, "Map scale factor on mass grid, y direction")
    create_var('MAPFAC_UX', ('Time', 'south_north', 'west_east_stag'), '', 1.0, 57.3, "Map scale factor on u-grid, x direction")
    create_var('MAPFAC_UY', ('Time', 'south_north', 'west_east_stag'), '', 1.0, 1.0, "Map scale factor on u-grid, y direction")
    create_var('MAPFAC_VX', ('Time', 'south_north_stag', 'west_east'), '', 0.0, 28.7, "Map scale factor on v-grid, x direction")
    create_var('MF_VX_INV', ('Time', 'south_north_stag', 'west_east'), '', 0.0, 1.0, "Inverse map scale factor on v-grid, x direction")
    create_var('MAPFAC_VY', ('Time', 'south_north_stag', 'west_east'), '', 1.0, 1.0, "Map scale factor on v-grid, y direction")
    create_var('F', ('Time', 'south_north', 'west_east'), 's-1', -1.4e-04, 1.4e-04, "Coriolis sine latitude term")
    create_var('E', ('Time', 'south_north', 'west_east'), 's-1', 2.5e-06, 1.4e-04, "Coriolis cosine latitude term")
    create_var('TSK', ('Time', 'south_north', 'west_east'), 'K', 142.0, 308.7, "SURFACE SKIN TEMPERATURE")
    create_var('HGT', ('Time', 'south_north', 'west_east'), 'm', -7.4e+03, 1.9e+04, "Terrain Height")
    create_var('RTHRATEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa K s-1', -14.9, 21.5, "COUPLED THETA TENDENCY DUE TO RADIATION")
    create_var('RTHRATLW', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K s-1', -2.9e-02, 5.3e-02, "UNCOUPLED THETA TENDENCY DUE TO LONG WAVE RADIATION")
    create_var('RTHRATSW', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K s-1', 0.0, 0.00517, "UNCOUPLED THETA TENDENCY DUE TO SHORT WAVE RADIATION")
    create_var('SWDOWN', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 649.1, "DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE")
    create_var('SWDOWNDIR', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 557.6, "DIRECT DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE")
    create_var('GSW', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 554.6, "NET SHORT WAVE FLUX AT GROUND SURFACE")
    create_var('GLW', ('Time', 'south_north', 'west_east'), 'W m-2', 1.0, 74.0, "DOWNWARD LONG WAVE FLUX AT GROUND SURFACE")
    create_var('HRVIS', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -0.0e+00, 2.9e-04, "HEATING RATE IN THE VISIBLE")
    create_var('HRIR', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -2.7e-02, 3.3e-02, "HEATING RATE IN THE INFRARED")
    create_var('HRAERVIS', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -5.7e-06, 2.4e-03, "HEATING RATE DUE TO AEROSOLS IN THE VISIBLE")
    create_var('HRAERIR', ('Time', 'bottom_top', 'south_north', 'west_east'), 'K/s', -1.1e-03, 1.1e-03, "HEATING RATE DUE TO AEROSOLS IN THE INFRARED")
    create_var('TOASW', ('Time', 'south_north', 'west_east'), 'W m-2', 0.0, 716.6, "DOWNWARD SHORT WAVE FLUX AT TOP OF ATMOSPHERE")
    create_var('TOALW', ('Time', 'south_north', 'west_east'), 'W m-2', 13.5, 419.4, "UPWARD LONG WAVE FLUX AT TOP OF ATMOSPHERE")
    create_var('ALBEDO', ('Time', 'south_north', 'west_east'), '-', 0.1, 0.8, "ALBEDO")
    create_var('CLAT', ('Time', 'south_north', 'west_east'), 'degree_north', -89.0, 89.0, "COMPUTATIONAL GRID LATITUDE, SOUTH IS NEGATIVE")
    create_var('CLONG', ('Time', 'south_north', 'west_east'), 'degree_east', -179.0, 179.0, "COMPUTATIONAL GRID LONGITUDE, WEST IS NEGATIVE")
    create_var('ALBBCK', ('Time', 'south_north', 'west_east'), '', 0.1, 0.5, "BACKGROUND ALBEDO")
    create_var('EMBCK', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0, "BACKGROUND EMISSIVITY")
    create_var('THCBCK', ('Time', 'south_north', 'west_east'), '', 30.0, 877.2, "BACKGROUND THERMAL INERTIA")
    create_var('EMISS', ('Time', 'south_north', 'west_east'), '', 0.5, 1.0, "SURFACE EMISSIVITY")
    create_var('RUBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa m s-2', -24.2, 24.4, "COUPLED X WIND TENDENCY DUE TO PBL PARAMETERIZATION")
    create_var('RVBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa m s-2', -24.6, 25.5, "COUPLED Y WIND TENDENCY DUE TO PBL PARAMETERIZATION")
    create_var('RTHBLTEN', ('Time', 'bottom_top', 'south_north', 'west_east'), 'Pa K s-1', -62.3, 115.4, "COUPLED THETA TENDENCY DUE TO PBL PARAMETERIZATION")
    create_var('XLAND', ('Time', 'south_north', 'west_east'), '', 1.0, 1.0, "LAND MASK (1 FOR LAND, 2 FOR WATER)")
    create_var('ZNT', ('Time', 'south_north', 'west_east'), 'm', 0.0, 0.1, "TIME-VARYING ROUGHNESS LENGTH")
    create_var('UST', ('Time', 'south_north', 'west_east'), 'm s-1', 0.0, 2.6, "U* IN SIMILARITY THEORY")
    create_var('PBLH', ('Time', 'south_north', 'west_east'), 'm', 6.5, 22000.0, "PBL HEIGHT")
    create_var('THC', ('Time', 'south_north', 'west_east'), 'J m-1 K-1 s-0.5', 30.0, 877.2, "THERMAL INERTIA")
    create_var('HFX', ('Time', 'south_north', 'west_east'), 'W m-2', -30.0, 57.3, "UPWARD HEAT FLUX AT THE SURFACE")
    create_var('RNET_2D', ('Time', 'south_north', 'west_east'), 'W m-2', -128.4, 289.9, "UPWARD RADIATIVE FLUX AT THE BOTTOM OF THE ATMOSPHERE")
    create_var('FLHC', ('Time', 'south_north', 'west_east'), '', 0.0, 2.0, "SURFACE EXCHANGE COEFFICIENT FOR HEAT")
    create_var('ANGSLOPE', ('Time', 'south_north', 'west_east'), 'radians', 6.4e-06, 1.7e-01, "Slope angle (magnitude)")
    create_var('AZMSLOPE', ('Time', 'south_north', 'west_east'), 'radians', 0.0, 6.3, "Slope azimuth")
    create_var('CO2ICE', ('Time', 'south_north', 'west_east'), 'kg/m^2', 0.0, 1697.5, "Surface CO2 ice")
    create_var('CDOD_SCALE', ('Time', 'south_north', 'west_east'), 'Unitless', 1.0, 1.0, "Column Dust optical Depth scale")
    create_var('TAU_OD2D', ('Time', 'south_north', 'west_east'), 'Unitless', 0.0, 2.1, "Dust Optical Depth normalized to 7mb")
    create_var('FRAC_PERM_CO2', ('Time', 'south_north', 'west_east'), 'fraction', 0.0, 1.0, "fraction of grid point covered in perm co2 ice")
    create_var('FRAC_PERM_H2O', ('Time', 'south_north', 'west_east'), 'fraction', 0.0, 1.0, "fraction of grid point covered in perm h2o ice")
    create_var('GRD_ICE_PC', ('Time', 'south_north', 'west_east'), 'percent', 0.0, 1.0, "% of soil volume occupied by h2o ice")
    create_var('GRD_ICE_DP', ('Time', 'south_north', 'west_east'), 'meters', -9999.0, 2.4, "depth to top of soil h2o ice layer")
    create_var('TAU_OD', ('Time', 'bottom_top', 'south_north', 'west_east'), 'unitless', 0.0, 0.9, "Optictal depth on full eta levels")

    # Add global scalar constants
    nc_file.setncattr('P0', 610.0)  # Reference pressure in Pa
    nc_file.setncattr('G', 3.72)    # Gravity on Mars in m/s²
    nc_file.setncattr('CP', 770.0)  # Specific heat capacity
    nc_file.setncattr('R_D', 192.0) # Gas constant for Mars atmosphere
    nc_file.setncattr('T0', 300.0)  # Reference temperature in K

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